#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <getopt.h> 
#include "genome.h"
#include "bitpacker.h"
#include "gpu_bridge.h"
#include "analysis.h"
#include "checkpoint.h"

#define TOOL_NAME       "DNA K-mer Global Uniqueness Analyzer (dkmer_analyzer)"
#define CREDIT_CONCEPT  "Scientific Concept & Design: Dr. Zaremba (Ukraine), creator of HASDI technology"
#define CREDIT_CODE     "Code Implementation:         Gemini AI"

// --- UI MACROS ---
#define TERM_ALT_SCREEN_ON  "\033[?1049h"
#define TERM_ALT_SCREEN_OFF "\033[?1049l"
#define TERM_HIDE_CURSOR    "\033[?25l"
#define TERM_SHOW_CURSOR    "\033[?25h"
#define TERM_CLEAR_SCREEN   "\033[2J"
#define TERM_HOME_CURSOR    "\033[H"

typedef struct {
    char *target_file;
    char *genome_file;
    char *file_prefix;
    int kmer_len;
    int threads;
    RankingMetric metric;   
    double ui_update_freq;  
    int save_freq_cycles;   
    int enable_snapshot;    
    int scan_minus_strand; 
    float filter_threshold; 
    int test_mode; 
    int no_cache;
} AppConfig;

AppConfig config;
volatile sig_atomic_t stop_requested = 0;

void handle_sigint(int sig) { (void)sig; stop_requested = 1; }

double get_time_diff(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

void format_time_smart(double seconds, char *buffer, size_t size) {
    if (seconds < 60.0) {
        snprintf(buffer, size, "%.1fs", seconds);
    } else if (seconds < 3600.0) {
        int sec_int = (int)seconds;
        snprintf(buffer, size, "%dm %ds", sec_int / 60, sec_int % 60);
    } else {
        int sec_int = (int)seconds;
        int h = sec_int / 3600;
        int m = (sec_int % 3600) / 60;
        snprintf(buffer, size, "%dh %dm", h, m);
    }
}

void print_usage() {
    printf("\nUsage: ./dkmer_analyzer [OPTIONS]\n");
    printf("Required:\n");
    printf("  -t, --target <file>   Target FASTA file\n");
    printf("  -g, --genome <file>   Genome FASTA file\n");
    printf("  -k, --kmer <int>      K-mer length\n");
    printf("\nOptions:\n");
    printf("  -p, --prefix <str>    Output prefix (default: 'dkmer')\n");
    printf("      --no-cache        Disable 2-bit genome caching (Force rebuild)\n");
    printf("  -c, --cpu <int>       CPU threads for packing\n");
    printf("  -f, --filter <float>  Filter threshold %% (default: 25.0)\n");
    printf("  -r, --rank <mode>     0=Worst, 1=Mean, 2=RMS, 3=GC_Asc, 4=GC_Desc\n");
    printf("  -u, --update <sec>    UI/RAM update frequency in seconds (default: 0.5)\n");
    printf("  -s, --save <cycles>   Checkpoint save frequency in batch cycles (default: 60)\n");
    printf("      --no-snapshot     Disable intermediate CSV snapshot generation\n");
    printf("      --plus-only       Scan ONLY the PLUS strand\n");
    printf("      --test            Enable extended Rank #1 details view\n");
    printf("\nCredits:\n  %s\n  %s\n", CREDIT_CONCEPT, CREDIT_CODE);
}

void init_config() {
    config.target_file = NULL;
    config.genome_file = NULL;
    config.file_prefix = "dkmer";
    config.kmer_len = 0;
    config.threads = 1;
    config.metric = SORT_WORST_CASE;
    config.ui_update_freq = 0.5;
    config.save_freq_cycles = 60; 
    config.enable_snapshot = 1;
    config.scan_minus_strand = 1;
    config.filter_threshold = 25.0; 
    config.test_mode = 0;
    config.no_cache = 0;
}

void draw_header_to_buffer(char *buffer, size_t max_len) {
    size_t off = strlen(buffer);
    snprintf(buffer + off, max_len - off, 
        "==========================================================================================\n"
        "                %s               \n"
        "------------------------------------------------------------------------------------------\n"
        " %s   \n"
        " %s                                                   \n"
        "==========================================================================================\n",
        TOOL_NAME, CREDIT_CONCEPT, CREDIT_CODE);
}

void restore_terminal() {
    printf("%s%s", TERM_ALT_SCREEN_OFF, TERM_SHOW_CURSOR);
    fflush(stdout);
}

int main(int argc, char *argv[]) {
    init_config();
    signal(SIGINT, handle_sigint);

    static struct option long_options[] = {
        {"target", required_argument, 0, 't'},
        {"genome", required_argument, 0, 'g'},
        {"kmer",   required_argument, 0, 'k'},
        {"prefix", required_argument, 0, 'p'},
        {"cpu",    required_argument, 0, 'c'},
        {"filter", required_argument, 0, 'f'},
        {"rank",   required_argument, 0, 'r'},
        {"update", required_argument, 0, 'u'}, 
        {"save",   required_argument, 0, 's'}, 
        {"no-cache", no_argument, &config.no_cache, 1},
        {"no-snapshot", no_argument, &config.enable_snapshot, 0}, 
        {"plus-only", no_argument, &config.scan_minus_strand, 0},
        {"test",      no_argument, &config.test_mode, 1},         
        {"help",        no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt; int idx;
    while ((opt = getopt_long(argc, argv, "t:g:k:p:c:f:r:u:s:h", long_options, &idx)) != -1) {
        switch (opt) {
            case 0: break; 
            case 't': config.target_file = optarg; break;
            case 'g': config.genome_file = optarg; break;
            case 'k': config.kmer_len = atoi(optarg); break;
            case 'p': config.file_prefix = optarg; break;
            case 'c': config.threads = atoi(optarg); break;
            case 'f': config.filter_threshold = atof(optarg); break;
            case 'r': config.metric = (RankingMetric)atoi(optarg); break;
            case 'u': config.ui_update_freq = atof(optarg); break; 
            case 's': config.save_freq_cycles = atoi(optarg); break; 
            case 'h': print_usage(); return 0;
            default: print_usage(); return 1;
        }
    }

    if (!config.target_file || !config.genome_file || config.kmer_len <= 0) {
        print_usage(); return 1;
    }

    printf("\n[Configuration]\n");
    printf("Target: %s | Genome: %s | K: %d\n", config.target_file, config.genome_file, config.kmer_len);
    printf("Prefix: %s | Mode: %s | Matches: Top-%d (Hard Limit)\n", 
           config.file_prefix, 
           config.test_mode ? "TEST" : "NORMAL",
           MAX_TOP_RESULTS);
    fflush(stdout);

    // --- LOGIC START ---
    
    char cache_file[512];
    snprintf(cache_file, sizeof(cache_file), "%s_genome_packed.bin", config.file_prefix);
    char map_file[512];
    snprintf(map_file, sizeof(map_file), "%s_genome_map.csv", config.file_prefix);

    PackedGenome *pg_plus = NULL;
    PackedGenome *pg_minus = NULL;
    CleanGenome *genome_plus = NULL; 
    CleanGenome *genome_minus = NULL;
    GenomeMap *map = NULL;
    TargetHitInfo target_info;
    memset(&target_info, 0, sizeof(TargetHitInfo));

    char *target_seq = load_clean_target(config.target_file);
    if (!target_seq) return 1;

    // 1. Try Load Cache
    int cache_loaded = 0;
    if (!config.no_cache) {
        int saved_k = 0;
        if (load_packed_cache(cache_file, &pg_plus, &pg_minus, &map, &saved_k, &target_info)) {
            if (saved_k == config.kmer_len) {
                printf("\n[+] Loaded 2-bit packed genome cache.\n");
                printf("    Target previously found: %s\n", target_info.found ? "YES" : "NO");
                cache_loaded = 1;
            } else {
                printf("\n[!] Cache K-mer mismatch (%d vs %d). Rebuilding.\n", saved_k, config.kmer_len);
                if (pg_plus) free_packed_genome(pg_plus);
                if (pg_minus) free_packed_genome(pg_minus);
                if (map) {
                     if (map->segments) free(map->segments);
                     if (map->chroms) free(map->chroms);
                     free(map);
                     map = NULL;
                }
            }
        }
    }

    // 2. Build if not cached
    if (!cache_loaded) {
        printf("\n[Step 1] Loading genome FASTA...\n");
        genome_plus = load_and_clean_genome(config.genome_file, config.kmer_len);
        if (!genome_plus) return 1;
        map = genome_plus->map; 

        if (config.scan_minus_strand) {
            printf("[Step 2] Generating MINUS strand...\n");
            genome_minus = generate_minus_strand(genome_plus);
        }

        printf("[Step 3] Identifying & Masking Target...\n");
        target_info = exclude_target(genome_plus, genome_minus, target_seq);

        printf("[Step 4] Packing Genome (2-bit)...\n");
        pg_plus = pack_genome_sequence(genome_plus->data, genome_plus->length);
        if (genome_minus) {
            pg_minus = pack_genome_sequence(genome_minus->data, genome_minus->length);
        }

        if (!config.no_cache) {
            printf("[Cache] Saving packed genome to disk...\n");
            save_packed_cache(cache_file, pg_plus, pg_minus, map, config.kmer_len, target_info);
            export_genome_map_csv(map_file, map, target_info);
        }
    }
    
    // 3. Prepare Tasks
    size_t target_len = strlen(target_seq);
    size_t num_kmers = (target_len >= (size_t)config.kmer_len) ? (target_len - config.kmer_len + 1) : 0;
    printf("\n[Info] Target Length: %zu bp | Unique K-mers: %zu\n", target_len, num_kmers);

    KmerTask *tasks = (KmerTask*)malloc(sizeof(KmerTask) * num_kmers);
    for (size_t i = 0; i < num_kmers; i++) {
        tasks[i].id = (uint32_t)i;
        tasks[i].worst_score_idx = 0;
        pack_sequence_chunk(&target_seq[i], tasks[i].packed_pattern, config.kmer_len);
        for (int j = 0; j < MAX_TOP_RESULTS; j++) tasks[i].top_results[j].similarity_score = -1.0f; 
    }

    // 4. Checkpoint Resume Logic
    char checkpoint_file[512];
    snprintf(checkpoint_file, sizeof(checkpoint_file), "%s_checkpoint.bin", config.file_prefix);
    size_t start_offset = 0;
    
    if (has_checkpoint(checkpoint_file)) {
        printf("\n[!] Resume analysis from %s? [y/n]: ", checkpoint_file);
        char answer = getchar();
        int c; while ((c = getchar()) != '\n' && c != EOF); 
        if (answer == 'y' || answer == 'Y') {
            size_t saved_processed;
            load_checkpoint(checkpoint_file, &start_offset, &saved_processed, config.kmer_len, num_kmers, tasks);
        } else {
            unlink(checkpoint_file);
        }
    }

    // 5. GPU Init
    printf("[Step 5] Initializing GPU...\n");
    if (init_gpu_memory(pg_plus, pg_minus, tasks, num_kmers) != 0) return 1;

    // --- MAIN LOOP ---
    size_t batch_size = 256 * 1024; 
    size_t total_len = pg_plus->total_bases; 
    struct timespec start_time, current_time, last_ui_update;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    last_ui_update = start_time;
    size_t batch_counter = 0;
    char frame_buffer[65536];

    printf("%s%s", TERM_ALT_SCREEN_ON, TERM_HIDE_CURSOR);
    fflush(stdout);

    frame_buffer[0] = '\0';
    strcat(frame_buffer, TERM_CLEAR_SCREEN TERM_HOME_CURSOR);
    draw_header_to_buffer(frame_buffer, sizeof(frame_buffer));
    size_t off = strlen(frame_buffer);
    snprintf(frame_buffer + off, sizeof(frame_buffer)-off, "\n[Status] Starting analysis...\n");
    printf("%s", frame_buffer);
    fflush(stdout);

    size_t current_offset = start_offset;
    for (; current_offset < total_len; current_offset += batch_size) {
        batch_counter++;
        if (stop_requested) break;

        size_t current_batch = (current_offset + batch_size > total_len) ? (total_len - current_offset) : batch_size;
        if (run_gpu_batch(current_offset, current_batch, config.kmer_len, total_len, config.filter_threshold) != 0) break;

        clock_gettime(CLOCK_MONOTONIC, &current_time);
        if (get_time_diff(last_ui_update, current_time) >= config.ui_update_freq) {
            retrieve_results_from_gpu(tasks, num_kmers);
            
            frame_buffer[0] = '\0';
            strcat(frame_buffer, TERM_CLEAR_SCREEN TERM_HOME_CURSOR);
            draw_header_to_buffer(frame_buffer, sizeof(frame_buffer));
            
            double pct = (double)(current_offset + current_batch) / total_len * 100.0;
            double elapsed = get_time_diff(start_time, current_time);
            double work_done_pct = pct - ((double)start_offset/total_len*100.0);
            double eta = (work_done_pct > 0.001) ? (elapsed / work_done_pct) * (100 - pct) : 0.0;
            
            char t_str[32], e_str[32];
            format_time_smart(elapsed, t_str, sizeof(t_str));
            format_time_smart(eta, e_str, sizeof(e_str));

            size_t off = strlen(frame_buffer);
            snprintf(frame_buffer + off, sizeof(frame_buffer)-off, 
                     "\n[Running] Offset: %zu / %zu (%.1f%%)\n"
                     "[Timing]  Elapsed: %s | ETA: %s\n", 
                     current_offset + current_batch, total_len, pct, t_str, e_str);
            
            write_interactive_report(frame_buffer, sizeof(frame_buffer), tasks, num_kmers, config.kmer_len, config.metric, config.test_mode, 0, map);
            printf("%s", frame_buffer);
            fflush(stdout);
            last_ui_update = current_time;
        }

        if (batch_counter >= (size_t)config.save_freq_cycles) {
             retrieve_results_from_gpu(tasks, num_kmers);
             save_checkpoint(checkpoint_file, current_offset + current_batch, current_offset + current_batch, config.kmer_len, num_kmers, tasks);
             if (config.enable_snapshot) {
                 char snap_file[512];
                 snprintf(snap_file, sizeof(snap_file), "%s_snapshot_k%d.csv", config.file_prefix, config.kmer_len);
                 // Передаємо target_info для експорту координат
                 save_snapshot_csv(tasks, num_kmers, config.kmer_len, snap_file, config.metric, map, target_info);
             }
             batch_counter = 0; 
        }
    }

    restore_terminal();

    if (stop_requested) {
        printf("\n[!] Interrupted. Saving...\n");
        retrieve_results_from_gpu(tasks, num_kmers);
        save_checkpoint(checkpoint_file, current_offset, current_offset, config.kmer_len, num_kmers, tasks);
    } else {
        printf("\n[Step 6] Analysis Complete. Exporting results...\n");
        
        // 1. Стандартний CSV
        char out_name[512]; 
        snprintf(out_name, sizeof(out_name), "%s_summary_k%d.csv", config.file_prefix, config.kmer_len);
        analyze_and_rank(tasks, num_kmers, config.kmer_len, out_name, config.metric, map, target_info);
        printf("   -> Summary CSV: %s\n", out_name);

        // 2. Повний текстовий звіт
        char full_txt_name[512];
        snprintf(full_txt_name, sizeof(full_txt_name), "%s_full_report_k%d.txt", config.file_prefix, config.kmer_len);
        save_full_text_report(tasks, num_kmers, config.kmer_len, full_txt_name, config.metric, map, target_info);

        // 3. Бінарна база даних
        char bin_name[512];
        snprintf(bin_name, sizeof(bin_name), "%s_full_data_k%d.bin", config.file_prefix, config.kmer_len);
        save_binary_export(tasks, num_kmers, config.kmer_len, bin_name, config.metric, map, target_info);

        unlink(checkpoint_file);
        
        // 4. Очищення снапшота
        char snap_name[512];
        snprintf(snap_name, sizeof(snap_name), "%s_snapshot_k%d.csv", config.file_prefix, config.kmer_len);
        unlink(snap_name);
        printf("   -> Cleanup: Snapshot removed.\n");
    }

    cleanup_gpu();
    free(tasks);
    free_packed_genome(pg_plus);
    if (pg_minus) free_packed_genome(pg_minus);
    free(target_seq);
    
    if (genome_plus) free_genome(genome_plus);
    if (genome_minus) { genome_minus->map = NULL; free_genome(genome_minus); }
    
    if (cache_loaded && map) {
        if (map->segments) free(map->segments);
        if (map->chroms) free(map->chroms);
        free(map);
    }

    return 0;
}
