#ifndef GPU_BRIDGE_H
#define GPU_BRIDGE_H

#include <stddef.h>
#include "bitpacker.h"

#ifdef __cplusplus
extern "C" {
#endif

// Ініціалізація пам'яті GPU
int init_gpu_memory(PackedGenome *pg_plus, PackedGenome *pg_minus, KmerTask *tasks, size_t num_tasks);

// Запуск обробки батчу. Додано параметр filter_threshold.
int run_gpu_batch(size_t genome_offset, size_t batch_len, int k, size_t total_genome_len, float filter_threshold);

// Отримання результатів
int retrieve_results_from_gpu(KmerTask *host_tasks, size_t num_tasks);

// Очищення
void cleanup_gpu();

#ifdef __cplusplus
}
#endif

#endif
