#include <cuda_runtime.h>
#include <stdio.h>
#include "gpu_bridge.h"
#include "bitpacker.h" 

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            printf("CUDA Error: %s at line %d\n", cudaGetErrorString(err), __LINE__); \
            return -1; \
        } \
    } while (0)

PackedGenome *d_pg_plus = NULL;
PackedGenome *d_pg_minus = NULL;
packed_unit_t *d_genome_data_plus = NULL;
packed_unit_t *d_genome_data_minus = NULL;
KmerTask *d_tasks = NULL;
size_t g_num_tasks = 0;

__device__ inline uint8_t get_genome_base(const packed_unit_t *data, size_t index) {
    size_t unit_idx = index / 16;
    size_t bit_offset = (index % 16) * 2;
    return (data[unit_idx] >> bit_offset) & 3;
}

// Функція вставки зі збереженням сортування (від найбільшого до найменшого)
// ДОДАНО: параметр uint8_t is_minus
__device__ void insert_sorted(KmerTask *task, uint32_t position, float score, const uint64_t *packed_seq, uint8_t is_minus) {
    const int MAX_STORE = MAX_TOP_RESULTS; 
    int count = task->worst_score_idx; 
    
    // 1. Пошук позиції вставки
    int insert_pos = -1;
    
    if (count < MAX_STORE) {
         insert_pos = count; 
    }

    for (int i = 0; i < count; i++) {
        if (score > task->top_results[i].similarity_score) {
            insert_pos = i;
            break;
        }
    }

    if (insert_pos == -1 && count >= MAX_STORE) return;

    // 2. Зсув елементів вправо
    int end_idx = (count < MAX_STORE) ? count : MAX_STORE - 1;
    
    for (int j = end_idx; j > insert_pos; j--) {
        task->top_results[j] = task->top_results[j-1];
    }

    // 3. Вставка
    task->top_results[insert_pos].similarity_score = score;
    task->top_results[insert_pos].position_in_genome = position;
    task->top_results[insert_pos].packed_seq[0] = packed_seq[0];
    task->top_results[insert_pos].packed_seq[1] = packed_seq[1];
    task->top_results[insert_pos].is_minus_strand = is_minus; // <--- ДОДАНО: Зберігаємо ланцюг

    // 4. Оновлення лічильника
    if (count < MAX_STORE) {
        task->worst_score_idx = count + 1;
    }
}

__global__ void kmer_scan_kernel(
    const packed_unit_t *genome, 
    size_t genome_len, 
    size_t start_offset, 
    size_t batch_len, 
    KmerTask *tasks, 
    size_t num_tasks, 
    int k,
    float filter_threshold,
    int is_minus_strand
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_tasks) return;

    KmerTask *my_task = &tasks[idx];
    
    int min_match_len = (int)(k * (filter_threshold / 100.0f));
    if (min_match_len < 5) min_match_len = 5;

    size_t end_offset = start_offset + batch_len;
    if (end_offset > genome_len - k) end_offset = genome_len - k;

    uint64_t pat0 = my_task->packed_pattern[0];
    uint64_t pat1 = my_task->packed_pattern[1];

    float cached_cutoff = -1.0f;
    
    if (my_task->worst_score_idx >= MAX_TOP_RESULTS) {
        cached_cutoff = my_task->top_results[MAX_TOP_RESULTS - 1].similarity_score;
    }

    for (size_t i = start_offset; i < end_offset; i++) {
        int match_count = 0;
        
        #pragma unroll
        for (int j = 0; j < k; j++) {
            size_t g_idx = i + j;
            size_t unit_idx = g_idx / 16;
            size_t bit_offset = (g_idx % 16) * 2;
            uint8_t g_base = (genome[unit_idx] >> bit_offset) & 3;
            
            uint8_t t_base;
            if (j < 32) t_base = (pat0 >> (j * 2)) & 3;
            else t_base = (pat1 >> ((j - 32) * 2)) & 3;

            if (g_base == t_base) match_count++;
        }

        if (match_count <= min_match_len) continue;

        float score = (float)match_count / k * 100.0f;

        if (cached_cutoff >= 0.0f && score <= cached_cutoff) continue;

        uint64_t found_packed[2] = {0, 0};
        for (int j = 0; j < k; j++) {
            size_t g_idx = i + j;
            size_t unit_idx = g_idx / 16;
            size_t bit_offset = (g_idx % 16) * 2;
            uint8_t g_base = (genome[unit_idx] >> bit_offset) & 3;

            if (j < 32) found_packed[0] |= ((uint64_t)g_base << (j * 2));
            else found_packed[1] |= ((uint64_t)g_base << ((j - 32) * 2));
        }

        // --- МАТЕМАТИКА МІНУС-ЛАНЦЮГА ---
        // Якщо це мінус-ланцюг, перераховуємо локальний індекс (i) у глобальний старт на прямому ланцюзі
        uint32_t forward_pos = is_minus_strand ? (uint32_t)(genome_len - i - k) : (uint32_t)i;

        // Вставляємо правильну глобальну позицію та ЛАНЦЮГ
        insert_sorted(my_task, forward_pos, score, found_packed, is_minus_strand); // <--- ОНОВЛЕНО

        if (my_task->worst_score_idx >= MAX_TOP_RESULTS) {
             cached_cutoff = my_task->top_results[MAX_TOP_RESULTS - 1].similarity_score;
        }
    }
}

extern "C" int init_gpu_memory(PackedGenome *pg_plus, PackedGenome *pg_minus, KmerTask *tasks, size_t num_tasks) {
    g_num_tasks = num_tasks;
    size_t genome_bytes_plus = pg_plus->data_len * sizeof(packed_unit_t);
    CUDA_CHECK(cudaMalloc((void**)&d_genome_data_plus, genome_bytes_plus));
    CUDA_CHECK(cudaMemcpy(d_genome_data_plus, pg_plus->data, genome_bytes_plus, cudaMemcpyHostToDevice));

    if (pg_minus != NULL) {
        size_t genome_bytes_minus = pg_minus->data_len * sizeof(packed_unit_t);
        CUDA_CHECK(cudaMalloc((void**)&d_genome_data_minus, genome_bytes_minus));
        CUDA_CHECK(cudaMemcpy(d_genome_data_minus, pg_minus->data, genome_bytes_minus, cudaMemcpyHostToDevice));
    } else {
        d_genome_data_minus = NULL;
    }

    size_t tasks_bytes = num_tasks * sizeof(KmerTask);
    CUDA_CHECK(cudaMalloc((void**)&d_tasks, tasks_bytes));
    CUDA_CHECK(cudaMemcpy(d_tasks, tasks, tasks_bytes, cudaMemcpyHostToDevice));
    return 0;
}

extern "C" int run_gpu_batch(size_t genome_offset, size_t batch_len, int k, size_t total_genome_len, float filter_threshold) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (g_num_tasks + threadsPerBlock - 1) / threadsPerBlock;

    // Запуск для ПРЯМОГО ланцюга (is_minus_strand = 0)
    kmer_scan_kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_genome_data_plus, total_genome_len, genome_offset, batch_len, d_tasks, g_num_tasks, k, filter_threshold, 0
    );

    if (d_genome_data_minus != NULL) {
        // Запуск для ЗВОРОТНОГО ланцюга (is_minus_strand = 1)
        kmer_scan_kernel<<<blocksPerGrid, threadsPerBlock>>>(
            d_genome_data_minus, total_genome_len, genome_offset, batch_len, d_tasks, g_num_tasks, k, filter_threshold, 1
        );
    }

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) return -1;
    cudaDeviceSynchronize();
    return 0;
}

extern "C" int retrieve_results_from_gpu(KmerTask *host_tasks, size_t num_tasks) {
    size_t tasks_bytes = num_tasks * sizeof(KmerTask);
    CUDA_CHECK(cudaMemcpy(host_tasks, d_tasks, tasks_bytes, cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" void cleanup_gpu() {
    if (d_genome_data_plus) cudaFree(d_genome_data_plus);
    if (d_genome_data_minus) cudaFree(d_genome_data_minus);
    if (d_tasks) cudaFree(d_tasks);
}
