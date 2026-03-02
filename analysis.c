#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "analysis.h"
#include "bitpacker.h"
#include "genome.h"

static RankingMetric current_sort_metric = SORT_WORST_CASE;

// --- Comparison functions ---
int compare_results(const void *a, const void *b) {
    const KmerAnalysisResult *ra = (const KmerAnalysisResult *)a;
    const KmerAnalysisResult *rb = (const KmerAnalysisResult *)b;
    switch (current_sort_metric) {
        case SORT_WORST_CASE: if (ra->worst_case_score > rb->worst_case_score) return 1; if (ra->worst_case_score < rb->worst_case_score) return -1; break;
        case SORT_MEAN_UNIQ: if (ra->mean_uniqueness < rb->mean_uniqueness) return 1; if (ra->mean_uniqueness > rb->mean_uniqueness) return -1; break;
        case SORT_RMS: if (ra->rms_score < rb->rms_score) return 1; if (ra->rms_score > rb->rms_score) return -1; break;
        case SORT_GC_ASC: if (ra->gc_content > rb->gc_content) return 1; if (ra->gc_content < rb->gc_content) return -1; break;
        case SORT_GC_DESC: if (ra->gc_content < rb->gc_content) return 1; if (ra->gc_content > rb->gc_content) return -1; break;
    }
    return 0;
}

int compare_overlaps_desc(const void *a, const void *b) {
    const OverlapResult *oa = (const OverlapResult *)a;
    const OverlapResult *ob = (const OverlapResult *)b;
    if (oa->similarity_score < ob->similarity_score) return 1;
    if (oa->similarity_score > ob->similarity_score) return -1;
    return 0;
}

// --- process_kmers ---
static void process_kmers(KmerTask *tasks, size_t num_tasks, int k, int limit_top_n, KmerAnalysisResult *results, GenomeMap *map) {
    char kmer_text[64];
    char match_text[64];
    if (limit_top_n > MAX_TOP_RESULTS) limit_top_n = MAX_TOP_RESULTS;

    for (size_t i = 0; i < num_tasks; i++) {
        KmerTask *t = &tasks[i];
        unpack_sequence(t->packed_pattern, kmer_text, k);
        strcpy(results[i].kmer_seq, kmer_text);
        results[i].original_id = t->id; 
        
        int gc_count = 0;
        for (int p = 0; p < k; p++) {
            char base = toupper(kmer_text[p]);
            if (base == 'G' || base == 'C') gc_count++;
        }
        results[i].gc_content = (double)gc_count / k * 100.0;

        int position_matches[64] = {0}; 
        int valid_entries = 0;
        double sum_uniqueness = 0.0;
        double max_similarity = -1.0;
        int max_sim_idx = -1;
        int count_limit = t->worst_score_idx;
        if (count_limit > MAX_TOP_RESULTS) count_limit = MAX_TOP_RESULTS;

        for (int j = 0; j < count_limit; j++) {
            if (t->top_results[j].similarity_score < 0) continue;
            valid_entries++;
            float score = t->top_results[j].similarity_score;
            double overlap_uniqueness = 100.0 - score;
            sum_uniqueness += overlap_uniqueness;
            if (score > max_similarity) { max_similarity = score; max_sim_idx = j; }

            unpack_sequence(t->top_results[j].packed_seq, match_text, k);
            for (int p = 0; p < k; p++) { if (kmer_text[p] == match_text[p]) position_matches[p]++; }
        }

        for (int p = 0; p < k; p++) {
            if (position_matches[p] == 0) results[i].visualized_seq[p] = toupper(kmer_text[p]);
            else results[i].visualized_seq[p] = tolower(kmer_text[p]);
        }
        results[i].visualized_seq[k] = '\0';

        if (valid_entries > 0) {
            results[i].mean_uniqueness = sum_uniqueness / valid_entries;
            results[i].worst_case_score = max_similarity;
            if (max_sim_idx != -1) {
                uint32_t clean_pos = t->top_results[max_sim_idx].position_in_genome;
                if (map) {
                    GenomicLocation loc = get_genomic_location(map, clean_pos);
                    results[i].worst_match_genome_pos = (uint32_t)loc.global_idx; 
                } else { results[i].worst_match_genome_pos = clean_pos; }
                unpack_sequence(t->top_results[max_sim_idx].packed_seq, results[i].worst_match_seq, k);
            } else { results[i].worst_match_genome_pos = 0; strcpy(results[i].worst_match_seq, "N/A"); }

            double sum_sq_uniqueness = 0.0;
            double min_uniq = 100.0;
            for (int p = 0; p < k; p++) {
                double match_freq = (double)position_matches[p] / valid_entries * 100.0;
                double uniqueness = 100.0 - match_freq;
                sum_sq_uniqueness += (uniqueness * uniqueness);
                if (uniqueness < min_uniq) min_uniq = uniqueness;
            }
            results[i].rms_score = sqrt(sum_sq_uniqueness / k);
            results[i].min_uniqueness = min_uniq;
        } else {
            results[i].mean_uniqueness = 100.0; results[i].worst_case_score = 0.0; results[i].rms_score = 100.0; results[i].min_uniqueness = 100.0;
            results[i].worst_match_genome_pos = 0; strcpy(results[i].worst_match_seq, "NONE");
        }
    }
    qsort(results, num_tasks, sizeof(KmerAnalysisResult), compare_results);
}

// --- Output details for Rank #1 ---
static void append_rank1_details(char *buffer, size_t max_len, KmerTask *tasks, size_t num_tasks, int k, uint32_t rank1_id, int safe_mode, GenomeMap *map) {
    size_t offset = strlen(buffer);
    if (offset >= max_len) return;
    snprintf(buffer + offset, max_len - offset, "--- DETAILS FOR RANK #1: Top Genome Matches ---\n"); offset = strlen(buffer);
    snprintf(buffer + offset, max_len - offset, "%-5s %-24s %-12s %-20s\n", "#", "Genome Sequence", "Matches", "Location"); offset = strlen(buffer);
    snprintf(buffer + offset, max_len - offset, "----------------------------------------------------------------\n"); offset = strlen(buffer);
    if (safe_mode) { for(int j=0; j<10; j++) { snprintf(buffer + offset, max_len - offset, "\n"); offset = strlen(buffer); } return; }

    KmerTask *target_task = NULL;
    for (size_t i = 0; i < num_tasks; i++) { if (tasks[i].id == rank1_id) { target_task = &tasks[i]; break; } }
    
    int count = 0;
    OverlapResult sorted_overlaps[MAX_TOP_RESULTS];
    if (target_task) {
        count = target_task->worst_score_idx;
        if (count > MAX_TOP_RESULTS) count = MAX_TOP_RESULTS; 
        memcpy(sorted_overlaps, target_task->top_results, sizeof(OverlapResult) * count);
        qsort(sorted_overlaps, count, sizeof(OverlapResult), compare_overlaps_desc);
    }

    for (int i = 0; i < 10; i++) {
        if (i < count) {
            char seq_buf[64];
            unpack_sequence(sorted_overlaps[i].packed_seq, seq_buf, k);
            int match_cnt = (int)(sorted_overlaps[i].similarity_score / 100.0 * k + 0.1);
            GenomicLocation loc;
            if (map) loc = get_genomic_location(map, sorted_overlaps[i].position_in_genome);
            else { loc.global_idx = sorted_overlaps[i].position_in_genome; strcpy(loc.chrom_name, "N/A"); loc.local_idx = sorted_overlaps[i].position_in_genome; }
            
            char m_strand = sorted_overlaps[i].is_minus_strand ? '-' : '+';
            char loc_str[128]; snprintf(loc_str, sizeof(loc_str), "%s:%zu(%c)", loc.chrom_name, loc.local_idx, m_strand);
            snprintf(buffer + offset, max_len - offset, "%-5d %-24s %d/%d        %-20s\n", i+1, seq_buf, match_cnt, k, loc_str);
        } else {
             if (count == 0 && i == 0) snprintf(buffer + offset, max_len - offset, "(No significant overlaps found yet)\n");
             else snprintf(buffer + offset, max_len - offset, "\n");
        }
        offset = strlen(buffer);
    }
}

// --- Interactive Report (UI) ---
void write_interactive_report(char *buffer, size_t max_len, KmerTask *tasks, size_t num_tasks, int k, RankingMetric metric, int show_details, int is_loading, GenomeMap *map) {
    KmerAnalysisResult *results = (KmerAnalysisResult*)malloc(sizeof(KmerAnalysisResult) * num_tasks);
    int safe_mode = (results == NULL);
    if (!safe_mode) { current_sort_metric = metric; process_kmers(tasks, num_tasks, k, MAX_TOP_RESULTS, results, map); }

    size_t offset = strlen(buffer);
    snprintf(buffer + offset, max_len - offset, "--- TOP 10 CANDIDATES (Safest First | Metric: %d) ---\n", metric); offset = strlen(buffer);
    snprintf(buffer + offset, max_len - offset, "%-5s %-24s ", "Rank", "Sequence"); offset = strlen(buffer);
    switch (metric) {
        case SORT_WORST_CASE: snprintf(buffer + offset, max_len - offset, "%-12s %-12s %-10s %-6s\n", "Diffs", "MeanUniq", "RMS", "GC%"); break;
        case SORT_RMS: snprintf(buffer + offset, max_len - offset, "%-10s %-12s %-12s %-6s\n", "RMS", "Diffs", "MeanUniq", "GC%"); break;
        case SORT_MEAN_UNIQ: snprintf(buffer + offset, max_len - offset, "%-12s %-12s %-10s %-6s\n", "MeanUniq", "Diffs", "RMS", "GC%"); break;
        default: snprintf(buffer + offset, max_len - offset, "%-6s %-12s %-12s %-10s\n", "GC%", "Diffs", "MeanUniq", "RMS"); break;
    }
    offset = strlen(buffer);
    snprintf(buffer + offset, max_len - offset, "-----------------------------------------------------------------------\n"); offset = strlen(buffer);
    
    for (int i = 0; i < 10; i++) {
        if (!safe_mode && i < (int)num_tasks) {
            if (is_loading) { snprintf(buffer + offset, max_len - offset, "%-5d %-24s %-12s %-12s %-10s %-6s\n", i+1, results[i].kmer_seq, "...", "...", "...", "..."); } else {
                int match_cnt = (int)(results[i].worst_case_score / 100.0 * k + 0.1);
                int diff_cnt = k - match_cnt;
                char diff_display[32]; snprintf(diff_display, sizeof(diff_display), "%d/%d", diff_cnt, k);
                char s_mean[16], s_rms[16], s_gc[16];
                snprintf(s_mean, sizeof(s_mean), "%.2f%%", results[i].mean_uniqueness);
                snprintf(s_rms, sizeof(s_rms), "%.2f%%", results[i].rms_score);
                snprintf(s_gc, sizeof(s_gc), "%.1f%%", results[i].gc_content);
                snprintf(buffer + offset, max_len - offset, "%-5d %-24s ", i+1, results[i].visualized_seq); offset = strlen(buffer);
                switch (metric) {
                    case SORT_WORST_CASE: snprintf(buffer + offset, max_len - offset, "%-12s %-12s %-10s %-6s\n", diff_display, s_mean, s_rms, s_gc); break;
                    case SORT_RMS: snprintf(buffer + offset, max_len - offset, "%-10s %-12s %-12s %-6s\n", s_rms, diff_display, s_mean, s_gc); break;
                    case SORT_MEAN_UNIQ: snprintf(buffer + offset, max_len - offset, "%-12s %-12s %-10s %-6s\n", s_mean, diff_display, s_rms, s_gc); break;
                    default: snprintf(buffer + offset, max_len - offset, "%-6s %-12s %-12s %-10s\n", s_gc, diff_display, s_mean, s_rms); break;
                }
            }
        } else { snprintf(buffer + offset, max_len - offset, "\n"); }
        offset = strlen(buffer);
    }
    if (show_details) {
        snprintf(buffer + offset, max_len - offset, "\n"); offset = strlen(buffer);
        if (is_loading) snprintf(buffer + offset, max_len - offset, "--- DETAILS (Waiting for analysis...) ---\n\n\n"); 
        else { if (!safe_mode && num_tasks > 0) append_rank1_details(buffer, max_len, tasks, num_tasks, k, results[0].original_id, 0, map); }
    }
    if (!safe_mode) free(results);
}

// =========================================================================================
// EXPORT FUNCTIONS 
// =========================================================================================

// --- 1. STANDARD CSV EXPORT ---
void analyze_and_rank(KmerTask *tasks, size_t num_tasks, int k, const char *output_csv, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info) {
    printf("\n[Step 6] Finalizing Analysis (Metric Code: %d)...\n", metric);
    KmerAnalysisResult *results = (KmerAnalysisResult*)malloc(sizeof(KmerAnalysisResult) * num_tasks);
    if (!results) return;

    current_sort_metric = metric;
    process_kmers(tasks, num_tasks, k, MAX_TOP_RESULTS, results, map);

    FILE *fp = fopen(output_csv, "w");
    if (!fp) { perror("Failed to open output file"); free(results); return; }

    fprintf(fp, "# Sorted by Metric: %d\n", metric);
    fprintf(fp, "Rank,Target_Kmer_Index,Target_Chrom,Target_Abs_Start,Target_Abs_End,Target_Strand,Sequence_Raw,Sequence_Visual,GC_Content,Differences,Mean_Uniqueness,RMS_Score,Min_Pos_Uniqueness,Worst_Match_Sequence,Worst_Match_Genome_Pos,Worst_Match_Chromosome,Worst_Match_Local_Pos\n");
    
    for (size_t i = 0; i < num_tasks; i++) {
        int match_cnt = (int)(results[i].worst_case_score / 100.0 * k + 0.1);
        int diff_cnt = k - match_cnt;
        
        char t_chrom[64] = "N/A";
        char t_strand[8] = "N/A";
        size_t t_start = 0;
        size_t t_end = 0;

        if (target_info.found) {
            strncpy(t_chrom, target_info.chrom_name, 63); t_chrom[63]='\0';
            strncpy(t_strand, target_info.strand, 7); t_strand[7]='\0';
            
            if (strcmp(target_info.strand, "MINUS") == 0) {
                t_start = target_info.local_end - results[i].original_id - k + 1;
                t_end   = target_info.local_end - results[i].original_id;
            } else {
                t_start = target_info.local_start + results[i].original_id;
                t_end   = t_start + k - 1;
            }
        }

        GenomicLocation loc;
        char worst_strand = '+';
        if (map) {
             KmerTask *t = &tasks[results[i].original_id];
             int max_sim_idx = -1; float max_score = -1.0;
             int count = t->worst_score_idx; 
             if (count > MAX_TOP_RESULTS) count = MAX_TOP_RESULTS;
             
             for(int j=0; j<count; j++) { if(t->top_results[j].similarity_score > max_score) { max_score = t->top_results[j].similarity_score; max_sim_idx = j; } }
             if(max_sim_idx != -1) {
                 loc = get_genomic_location(map, t->top_results[max_sim_idx].position_in_genome);
                 worst_strand = t->top_results[max_sim_idx].is_minus_strand ? '-' : '+';
             } else { loc.global_idx = 0; strcpy(loc.chrom_name, "N/A"); loc.local_idx = 0; }
        } else {
             loc.global_idx = results[i].worst_match_genome_pos; strcpy(loc.chrom_name, "N/A"); loc.local_idx = 0;
        }

        fprintf(fp, "%zu,%u,%s,%zu,%zu,%s,%s,%s,%.2f,%d/%d,%.4f,%.4f,%.4f,%s,%zu,%s,%zu(%c)\n", 
                i+1, results[i].original_id,
                t_chrom, t_start, t_end, t_strand,
                results[i].kmer_seq, results[i].visualized_seq,
                results[i].gc_content, diff_cnt, k, results[i].mean_uniqueness, 
                results[i].rms_score, results[i].min_uniqueness,
                results[i].worst_match_seq, 
                loc.global_idx, loc.chrom_name, loc.local_idx, worst_strand);
    }
    fclose(fp);
    free(results);
}

void save_snapshot_csv(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info) {
    analyze_and_rank(tasks, num_tasks, k, filename, metric, map, target_info);
}

// --- 2. FULL TEXT REPORT (HUMAN READABLE) ---
void save_full_text_report(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info) {
    printf("[Export] Generating Full Text Report '%s'...\n", filename);
    KmerAnalysisResult *results = (KmerAnalysisResult*)malloc(sizeof(KmerAnalysisResult) * num_tasks);
    if (!results) return;

    current_sort_metric = metric;
    process_kmers(tasks, num_tasks, k, MAX_TOP_RESULTS, results, map);

    FILE *fp = fopen(filename, "w");
    if (!fp) { perror("Failed to open text report file"); free(results); return; }

    fprintf(fp, "================================================================================\n");
    fprintf(fp, "  FULL GENOME UNIQUENESS REPORT\n");
    fprintf(fp, "  Sort Metric: %d | K-mer: %d | Total K-mers: %zu\n", metric, k, num_tasks);
    
    if (target_info.found) {
        fprintf(fp, "  TARGET: Found in %s (%s Strand)\n", target_info.chrom_name, target_info.strand);
        fprintf(fp, "  REGION: %zu - %zu (Local)\n", target_info.local_start, target_info.local_end);
    } else {
        fprintf(fp, "  TARGET: Not found or masked manually.\n");
    }
    fprintf(fp, "================================================================================\n\n");

    char seq_buf[64];

    for (size_t i = 0; i < num_tasks; i++) {
        KmerTask *t = &tasks[results[i].original_id];
        
        char t_loc_str[128] = "N/A";
        if (target_info.found) {
            size_t t_start = 0;
            size_t t_end = 0;
            char t_strand_char = (strcmp(target_info.strand, "MINUS") == 0) ? '-' : '+';
            if (t_strand_char == '-') {
                t_start = target_info.local_end - t->id - k + 1;
                t_end   = target_info.local_end - t->id;
            } else {
                t_start = target_info.local_start + t->id;
                t_end   = t_start + k - 1;
            }
            snprintf(t_loc_str, sizeof(t_loc_str), "%s:%zu-%zu(%c)", target_info.chrom_name, t_start, t_end, t_strand_char);
        }

        fprintf(fp, "№%zu | ID: %u | Seq: %s | Target: %s\n", i+1, t->id, results[i].kmer_seq, t_loc_str);
        fprintf(fp, "Stats: MeanUniq: %.2f%% | RMS: %.2f%% | MinUniq: %.2f%% | GC: %.1f%%\n", 
                results[i].mean_uniqueness, results[i].rms_score, results[i].min_uniqueness, results[i].gc_content);
        
        fprintf(fp, "Matches (Top-%d):\n", MAX_TOP_RESULTS);
        fprintf(fp, "   %-5s %-24s %-16s %-18s %-20s %-25s\n", 
                "Idx", "Sequence", "Matches", "Clean_Range", "Global_Orig_Range", "Chrom_Location");
        fprintf(fp, "   ----- ------------------------ ---------------- ------------------ -------------------- -------------------------\n");

        int count = t->worst_score_idx;
        if (count > MAX_TOP_RESULTS) count = MAX_TOP_RESULTS;

        OverlapResult sorted_overlaps[MAX_TOP_RESULTS];
        memcpy(sorted_overlaps, t->top_results, sizeof(OverlapResult) * count);
        qsort(sorted_overlaps, count, sizeof(OverlapResult), compare_overlaps_desc);

        for (int j = 0; j < count; j++) {
            unpack_sequence(sorted_overlaps[j].packed_seq, seq_buf, k);
            
            GenomicLocation loc;
            if (map) {
                loc = get_genomic_location(map, sorted_overlaps[j].position_in_genome);
            } else {
                loc.global_idx = 0; strcpy(loc.chrom_name, "N/A"); loc.local_idx = 0;
            }

            uint32_t clean_start = sorted_overlaps[j].position_in_genome;
            uint32_t clean_end   = clean_start + k - 1;
            size_t global_start  = loc.global_idx;
            size_t global_end    = loc.global_idx + k - 1;
            size_t local_start   = loc.local_idx;
            size_t local_end     = loc.local_idx + k - 1;

            char range_clean[64], range_global[64], range_local[128], score_str[64];
            
            int match_cnt = (int)(sorted_overlaps[j].similarity_score / 100.0 * k + 0.1);
            snprintf(score_str, sizeof(score_str), "%d/%d (%.1f%%)", match_cnt, k, sorted_overlaps[j].similarity_score);
            snprintf(range_clean, sizeof(range_clean), "%u-%u", clean_start, clean_end);
            snprintf(range_global, sizeof(range_global), "%zu-%zu", global_start, global_end);
            
            char m_strand = sorted_overlaps[j].is_minus_strand ? '-' : '+';
            snprintf(range_local, sizeof(range_local), "%s:%zu-%zu(%c)", loc.chrom_name, local_start, local_end, m_strand);

            fprintf(fp, "   %-5d %-24s %-16s %-18s %-20s %s\n",
                    j+1, seq_buf, score_str, range_clean, range_global, range_local);
        }
        fprintf(fp, "\n--------------------------------------------------------------------------------\n\n");
    }

    fclose(fp);
    free(results);
    printf("   -> Text report saved.\n");
}

// --- 3. BINARY EXPORT ---

typedef struct {
    char magic[6];      
    uint32_t k;
    uint64_t num_tasks;
    uint32_t max_matches; 
    
    int target_found;
    char target_chrom[64];
    char target_strand[8];
    uint64_t target_start;
    uint64_t target_end;
} BinaryHeader;

typedef struct {
    uint32_t rank;
    uint32_t id;
    char sequence[64];
    double mean_uniq;
    double rms;
    double min_uniq;
    double gc;
    uint32_t match_count; 
    
    uint64_t kmer_target_start;
    uint64_t kmer_target_end;
} BinaryKmerHeader;

typedef struct {
    char sequence[64];
    float similarity;
    uint32_t clean_start;
    uint32_t clean_end;
    uint64_t global_orig_start;
    uint64_t global_orig_end;
    char chrom_name[64];
    uint64_t local_start;
    uint64_t local_end;
    char strand; // '+' або '-'
} BinaryMatchRecord;

void save_binary_export(KmerTask *tasks, size_t num_tasks, int k, const char *filename, RankingMetric metric, GenomeMap *map, TargetHitInfo target_info) {
    printf("[Export] Generating Binary Database '%s'...\n", filename);
    KmerAnalysisResult *results = (KmerAnalysisResult*)malloc(sizeof(KmerAnalysisResult) * num_tasks);
    if (!results) return;

    current_sort_metric = metric;
    process_kmers(tasks, num_tasks, k, MAX_TOP_RESULTS, results, map);

    FILE *fp = fopen(filename, "wb");
    if (!fp) { perror("Failed to open binary file"); free(results); return; }

    BinaryHeader header;
    memset(&header, 0, sizeof(BinaryHeader));
    memcpy(header.magic, "DKBIN1", sizeof(header.magic));
    header.k = (uint32_t)k;
    header.num_tasks = (uint64_t)num_tasks;
    header.max_matches = MAX_TOP_RESULTS; 
    
    header.target_found = target_info.found;
    if (target_info.found) {
        strncpy(header.target_chrom, target_info.chrom_name, 63);
        strncpy(header.target_strand, target_info.strand, 7);
        header.target_start = target_info.local_start;
        header.target_end = target_info.local_end;
    }

    fwrite(&header, sizeof(BinaryHeader), 1, fp);

    char seq_buf[64];

    for (size_t i = 0; i < num_tasks; i++) {
        KmerTask *t = &tasks[results[i].original_id];
        
        int count = t->worst_score_idx;
        if (count > MAX_TOP_RESULTS) count = MAX_TOP_RESULTS;

        BinaryKmerHeader k_header;
        memset(&k_header, 0, sizeof(BinaryKmerHeader));
        k_header.rank = (uint32_t)(i + 1);
        k_header.id = t->id;
        strncpy(k_header.sequence, results[i].kmer_seq, 63);
        k_header.mean_uniq = results[i].mean_uniqueness;
        k_header.rms = results[i].rms_score;
        k_header.min_uniq = results[i].min_uniqueness;
        k_header.gc = results[i].gc_content;
        k_header.match_count = (uint32_t)count;
        
        if (target_info.found) {
            if (strcmp(target_info.strand, "MINUS") == 0) {
                k_header.kmer_target_start = target_info.local_end - t->id - k + 1;
                k_header.kmer_target_end   = target_info.local_end - t->id;
            } else {
                k_header.kmer_target_start = target_info.local_start + t->id;
                k_header.kmer_target_end   = k_header.kmer_target_start + k - 1;
            }
        }
        
        fwrite(&k_header, sizeof(BinaryKmerHeader), 1, fp);

        OverlapResult sorted_overlaps[MAX_TOP_RESULTS];
        memcpy(sorted_overlaps, t->top_results, sizeof(OverlapResult) * count);
        qsort(sorted_overlaps, count, sizeof(OverlapResult), compare_overlaps_desc);

        for (int j = 0; j < count; j++) {
            BinaryMatchRecord match_rec;
            memset(&match_rec, 0, sizeof(BinaryMatchRecord));
            unpack_sequence(sorted_overlaps[j].packed_seq, seq_buf, k);
            strncpy(match_rec.sequence, seq_buf, 63);
            match_rec.similarity = sorted_overlaps[j].similarity_score;
            match_rec.clean_start = sorted_overlaps[j].position_in_genome;
            match_rec.clean_end   = match_rec.clean_start + k - 1;
            match_rec.strand = sorted_overlaps[j].is_minus_strand ? '-' : '+'; // <--- ДОДАНО
            
            GenomicLocation loc;
            if (map) {
                loc = get_genomic_location(map, sorted_overlaps[j].position_in_genome);
                match_rec.global_orig_start = loc.global_idx;
                match_rec.global_orig_end   = loc.global_idx + k - 1;
                strncpy(match_rec.chrom_name, loc.chrom_name, 63);
                match_rec.local_start = loc.local_idx;
                match_rec.local_end   = loc.local_idx + k - 1;
            }
            fwrite(&match_rec, sizeof(BinaryMatchRecord), 1, fp);
        }
    }
    fclose(fp);
    free(results);
    printf("   -> Binary database saved.\n");
}
