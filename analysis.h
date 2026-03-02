#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "bitpacker.h"
#include "genome.h" 
#include <stddef.h> 

typedef enum {
    SORT_WORST_CASE = 0, 
    SORT_MEAN_UNIQ,      
    SORT_RMS,            
    SORT_GC_ASC,         
    SORT_GC_DESC         
} RankingMetric;

typedef struct {
    uint32_t original_id;
    char kmer_seq[64];
    char visualized_seq[64];
    
    double gc_content;         
    double mean_uniqueness;    
    double worst_case_score;   
    double rms_score;          
    double min_uniqueness;     

    uint32_t worst_match_genome_pos; 
    char worst_match_seq[64];        
} KmerAnalysisResult;

// --- EXPORT FUNCTIONS (Updated with TargetHitInfo) ---

// 1. Стандартний CSV
void analyze_and_rank(KmerTask *tasks, size_t num_tasks, int k, const char *output_csv, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info);

// 2. Інтерактивний вивід (без змін, це UI)
void write_interactive_report(char *buffer, size_t max_len, KmerTask *tasks, size_t num_tasks, int k, RankingMetric metric, int show_details, int is_loading, GenomeMap *map);

// 3. Снапшот
void save_snapshot_csv(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info);

// 4. Повний текстовий звіт
void save_full_text_report(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info);

// 5. Бінарний експорт
void save_binary_export(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info);

#endif
