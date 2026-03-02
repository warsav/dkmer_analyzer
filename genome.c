#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "genome.h"

#define BUFFER_SIZE 65536 

// Магічне число для перевірки формату кешу (ASCII "DK2B")
static const uint32_t CACHE_MAGIC = 0x42324B44; 

// --- Допоміжні функції ---

int is_major_chromosome(const char *name) {
    if (strncmp(name, "CM", 2) == 0) return 1;
    if (strncmp(name, "J0", 2) == 0) return 1; 
    return 0;
}

void flush_scaffold_buffer(char *buffer, int *count) {
    if (*count > 0) {
        printf("   [+ Segments] %d non-chromosomal sequences: (%s)\n", *count, buffer);
        buffer[0] = '\0'; 
        *count = 0;
    }
}

void add_segment(GenomeMap *map, size_t c_off, size_t o_off, size_t len) {
    if (len == 0) return;
    if (map->count >= map->capacity) {
        map->capacity = (map->capacity == 0) ? 1024 : map->capacity * 2;
        GenomeSegment *new_seg = realloc(map->segments, sizeof(GenomeSegment) * map->capacity);
        if (!new_seg) { perror("Realloc failed"); exit(1); }
        map->segments = new_seg;
    }
    map->segments[map->count].clean_offset = c_off;
    map->segments[map->count].original_offset = o_off;
    map->segments[map->count].length = len;
    map->count++;
}

void add_chromosome(GenomeMap *map, const char *name, size_t start_idx) {
    if (map->chrom_count >= map->chrom_capacity) {
        map->chrom_capacity = (map->chrom_capacity == 0) ? 64 : map->chrom_capacity * 2;
        ChromosomeEntry *new_chroms = realloc(map->chroms, sizeof(ChromosomeEntry) * map->chrom_capacity);
        if (!new_chroms) { perror("Realloc failed"); exit(1); }
        map->chroms = new_chroms;
    }
    ChromosomeEntry *entry = &map->chroms[map->chrom_count];
    
    strncpy(entry->name, name, 63);
    entry->name[63] = '\0';
    
    entry->global_start_idx = start_idx;
    entry->length = 0; 
    
    if (map->chrom_count > 0) {
        ChromosomeEntry *prev = &map->chroms[map->chrom_count - 1];
        prev->length = start_idx - prev->global_start_idx;
    }
    
    map->chrom_count++;
}

// --- ЗАВАНТАЖЕННЯ ГЕНОМУ (FASTA) ---

CleanGenome* load_and_clean_genome(const char *filename, int k) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("File open error"); return NULL; }

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);

    CleanGenome *cg = (CleanGenome*)malloc(sizeof(CleanGenome));
    cg->data = (char*)malloc(fsize + 1); 
    cg->kmer_len = k;
    cg->map = (GenomeMap*)calloc(1, sizeof(GenomeMap));
    
    size_t write_idx = 0;
    
    // original_idx тепер рахує БАЗИ, а не байти файлу.
    // Це виправляє зсув координат.
    size_t original_base_count = 0; 
    
    char line[BUFFER_SIZE];
    char current_chrom[256] = "";
    int in_gap = 1; 
    
    size_t segment_start_clean = 0;
    size_t segment_start_orig = 0;

    char scaffold_names[1024] = "";
    int scaffold_count = 0;

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') line[--len] = '\0';

        if (line[0] == '>') {
            // Закриваємо попередній сегмент
            if (!in_gap) {
                add_segment(cg->map, segment_start_clean, segment_start_orig, write_idx - segment_start_clean);
                in_gap = 1;
            }

            char *name_start = line + 1;
            char *space = strchr(name_start, ' ');
            if (space) *space = '\0';
            
            strncpy(current_chrom, name_start, 255);
            current_chrom[255] = '\0';
            
            // Тут ми використовуємо поточний лічильник баз як старт хромосоми
            add_chromosome(cg->map, current_chrom, original_base_count);

            if (is_major_chromosome(current_chrom)) {
                flush_scaffold_buffer(scaffold_names, &scaffold_count);
                printf("   [Chr Detected] %s at genomic index %zu\n", current_chrom, original_base_count);
            } else {
                if (scaffold_count < 10) {
                    size_t current_len = strlen(scaffold_names);
                    size_t remaining = sizeof(scaffold_names) - current_len - 1;
                    if (remaining > 0) {
                        if (scaffold_count > 0) snprintf(scaffold_names + current_len, remaining, ", %.20s", current_chrom);
                        else snprintf(scaffold_names + current_len, remaining, "%.20s", current_chrom);
                    }
                } else if (scaffold_count == 10) {
                    strncat(scaffold_names, ", ...", 5);
                }
                scaffold_count++;
            }
            
            // УВАГА: Ми НЕ збільшуємо original_base_count на довжину заголовка!
            // Координати не враховують хедери.
            continue;
        }

        for (size_t i = 0; i < len; i++) {
            char c = toupper(line[i]);
            
            // Ми інкрементуємо лічильник баз (оригінальних координат)
            // тільки якщо це нуклеотид або N.
            if ((c >= 'A' && c <= 'Z')) {
                original_base_count++;
            }

            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                if (in_gap) {
                    in_gap = 0;
                    segment_start_clean = write_idx;
                    // start_orig має вказувати на індекс першої літери цього сегменту
                    segment_start_orig = original_base_count - 1; 
                }
                cg->data[write_idx++] = c;
            } else {
                // N або інші символи (Gap)
                if (!in_gap) {
                    add_segment(cg->map, segment_start_clean, segment_start_orig, write_idx - segment_start_clean);
                    in_gap = 1;
                }
            }
        }
        // УВАГА: Ми НЕ збільшуємо original_base_count на символ \n!
    }
    
    if (!in_gap) {
        add_segment(cg->map, segment_start_clean, segment_start_orig, write_idx - segment_start_clean);
    }
    
    if (cg->map->chrom_count > 0) {
        ChromosomeEntry *last = &cg->map->chroms[cg->map->chrom_count - 1];
        last->length = original_base_count - last->global_start_idx;
    }
    
    flush_scaffold_buffer(scaffold_names, &scaffold_count);

    cg->length = write_idx;
    cg->data[write_idx] = '\0'; 
    
    char *final_data = realloc(cg->data, write_idx + 1);
    if (final_data) cg->data = final_data;
    
    cg->map->total_original_length = original_base_count;

    fclose(fp);
    printf("   [Map] Built %zu segments across %zu chromosomes/scaffolds.\n", cg->map->count, cg->map->chrom_count);
    return cg;
}

// --- КЕШУВАННЯ (2-BIT) ---

int save_packed_cache(const char *filename, PackedGenome *pg_plus, PackedGenome *pg_minus, GenomeMap *map, int k, TargetHitInfo target_info) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) { perror("Save cache failed"); return 0; }

    fwrite(&CACHE_MAGIC, sizeof(uint32_t), 1, fp);
    fwrite(&k, sizeof(int), 1, fp);
    
    fwrite(&map->count, sizeof(size_t), 1, fp);
    fwrite(&map->total_original_length, sizeof(size_t), 1, fp);
    if (map->count > 0) fwrite(map->segments, sizeof(GenomeSegment), map->count, fp);

    fwrite(&map->chrom_count, sizeof(size_t), 1, fp);
    if (map->chrom_count > 0) fwrite(map->chroms, sizeof(ChromosomeEntry), map->chrom_count, fp);

    fwrite(&target_info, sizeof(TargetHitInfo), 1, fp);

    fwrite(&pg_plus->total_bases, sizeof(size_t), 1, fp);
    fwrite(&pg_plus->data_len, sizeof(size_t), 1, fp);
    fwrite(pg_plus->data, sizeof(packed_unit_t), pg_plus->data_len, fp);

    int has_minus = (pg_minus != NULL);
    fwrite(&has_minus, sizeof(int), 1, fp);
    if (has_minus) {
        fwrite(&pg_minus->total_bases, sizeof(size_t), 1, fp);
        fwrite(&pg_minus->data_len, sizeof(size_t), 1, fp);
        fwrite(pg_minus->data, sizeof(packed_unit_t), pg_minus->data_len, fp);
    }

    fclose(fp);
    
    double size_mb = (double)((pg_plus->data_len * 4) + (has_minus ? pg_minus->data_len * 4 : 0)) / (1024*1024);
    printf("   [Cache] Saved 2-bit packed genome to '%s' (%.2f MB)\n", filename, size_mb);
    return 1;
}

int load_packed_cache(const char *filename, PackedGenome **pg_plus_out, PackedGenome **pg_minus_out, GenomeMap **map_out, int *k_out, TargetHitInfo *target_info_out) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;

    uint32_t magic;
    if (fread(&magic, sizeof(uint32_t), 1, fp) != 1 || magic != CACHE_MAGIC) {
        fclose(fp); return 0; 
    }

    if (fread(k_out, sizeof(int), 1, fp) != 1) { fclose(fp); return 0; }

    GenomeMap *map = (GenomeMap*)calloc(1, sizeof(GenomeMap));
    if (fread(&map->count, sizeof(size_t), 1, fp) != 1) goto load_err;
    if (fread(&map->total_original_length, sizeof(size_t), 1, fp) != 1) goto load_err;
    
    map->capacity = map->count;
    if (map->count > 0) {
        map->segments = malloc(sizeof(GenomeSegment) * map->count);
        if (fread(map->segments, sizeof(GenomeSegment), map->count, fp) != map->count) goto load_err;
    }
    
    if (fread(&map->chrom_count, sizeof(size_t), 1, fp) != 1) goto load_err;
    map->chrom_capacity = map->chrom_count;
    if (map->chrom_count > 0) {
        map->chroms = malloc(sizeof(ChromosomeEntry) * map->chrom_count);
        if (fread(map->chroms, sizeof(ChromosomeEntry), map->chrom_count, fp) != map->chrom_count) goto load_err;
    }
    *map_out = map;

    if (fread(target_info_out, sizeof(TargetHitInfo), 1, fp) != 1) goto load_err;

    PackedGenome *pg_plus = (PackedGenome*)malloc(sizeof(PackedGenome));
    if (fread(&pg_plus->total_bases, sizeof(size_t), 1, fp) != 1) { free(pg_plus); goto load_err; }
    if (fread(&pg_plus->data_len, sizeof(size_t), 1, fp) != 1) { free(pg_plus); goto load_err; }
    
    pg_plus->data = (packed_unit_t*)malloc(sizeof(packed_unit_t) * pg_plus->data_len);
    if (fread(pg_plus->data, sizeof(packed_unit_t), pg_plus->data_len, fp) != pg_plus->data_len) {
        free(pg_plus->data); free(pg_plus); goto load_err;
    }
    *pg_plus_out = pg_plus;

    int has_minus;
    if (fread(&has_minus, sizeof(int), 1, fp) != 1) goto load_err;
    
    if (has_minus) {
        PackedGenome *pg_minus = (PackedGenome*)malloc(sizeof(PackedGenome));
        if (fread(&pg_minus->total_bases, sizeof(size_t), 1, fp) != 1) { free(pg_minus); goto load_err; }
        if (fread(&pg_minus->data_len, sizeof(size_t), 1, fp) != 1) { free(pg_minus); goto load_err; }
        
        pg_minus->data = (packed_unit_t*)malloc(sizeof(packed_unit_t) * pg_minus->data_len);
        if (fread(pg_minus->data, sizeof(packed_unit_t), pg_minus->data_len, fp) != pg_minus->data_len) {
            free(pg_minus->data); free(pg_minus); goto load_err;
        }
        *pg_minus_out = pg_minus;
    } else {
        *pg_minus_out = NULL;
    }

    fclose(fp);
    return 1;

load_err:
    fclose(fp);
    if (map->segments) free(map->segments);
    if (map->chroms) free(map->chroms);
    free(map);
    return 0;
}

int export_genome_map_csv(const char *filename, GenomeMap *map, TargetHitInfo target_info) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return 0;

    fprintf(fp, "# Genome Map Export\n");
    fprintf(fp, "# Target Status: %s\n", target_info.found ? "FOUND & MASKED" : "NOT FOUND");
    
    if (target_info.found) {
        fprintf(fp, "# Target Location: %s (%s Strand)\n", target_info.chrom_name, target_info.strand);
        fprintf(fp, "# Global Range: %zu - %zu\n", target_info.global_start, target_info.global_end);
        fprintf(fp, "# Local Range:  %zu - %zu\n\n", target_info.local_start, target_info.local_end);
    } else {
        fprintf(fp, "\n");
    }

    fprintf(fp, "Type,Name,Length_bp,Global_Start_Idx,Global_End_Idx\n");
    
    if (target_info.found) {
        fprintf(fp, "TARGET,%s_TARGET_REGION,%zu,%zu,%zu\n", 
                target_info.chrom_name,
                target_info.global_end - target_info.global_start + 1,
                target_info.global_start,
                target_info.global_end);
    }

    for (size_t i = 0; i < map->chrom_count; i++) {
        ChromosomeEntry *chr = &map->chroms[i];
        size_t end_idx = chr->global_start_idx + chr->length;
        if (i < map->chrom_count - 1) end_idx = map->chroms[i+1].global_start_idx - 1;

        fprintf(fp, "CHROM,%s,%zu,%zu,%zu\n", 
                chr->name, chr->length, chr->global_start_idx, end_idx);
    }

    fclose(fp);
    printf("   [Map] Enhanced genome map exported to '%s'\n", filename);
    return 1;
}

// --- ІНШІ ФУНКЦІЇ ---

void free_genome(CleanGenome *genome) {
    if (!genome) return;
    if (genome->data) free(genome->data);
    if (genome->map) {
        if (genome->map->segments) free(genome->map->segments);
        if (genome->map->chroms) free(genome->map->chroms);
        free(genome->map);
    }
    free(genome);
}

GenomicLocation get_genomic_location(GenomeMap *map, size_t clean_idx) {
    GenomicLocation loc;
    loc.global_idx = 0;
    strcpy(loc.chrom_name, "Unknown");
    loc.local_idx = 0;

    int left = 0, right = map->count - 1;
    int seg_idx = -1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        GenomeSegment *s = &map->segments[mid];
        if (clean_idx >= s->clean_offset && clean_idx < s->clean_offset + s->length) {
            seg_idx = mid; break;
        }
        if (clean_idx < s->clean_offset) right = mid - 1;
        else left = mid + 1;
    }

    if (seg_idx != -1) {
        GenomeSegment *s = &map->segments[seg_idx];
        size_t offset_in_segment = clean_idx - s->clean_offset;
        loc.global_idx = s->original_offset + offset_in_segment;
        
        left = 0; right = map->chrom_count - 1;
        int chr_idx = -1;
        while (left <= right) {
            int mid = left + (right - left) / 2;
            ChromosomeEntry *c = &map->chroms[mid];
            int is_start_ok = (loc.global_idx >= c->global_start_idx);
            int is_end_ok = 1;
            if (mid < (int)map->chrom_count - 1) {
                if (loc.global_idx >= map->chroms[mid+1].global_start_idx) is_end_ok = 0;
            }
            if (is_start_ok && is_end_ok) { chr_idx = mid; break; }
            if (!is_start_ok) right = mid - 1;
            else left = mid + 1;
        }

        if (chr_idx != -1) {
            ChromosomeEntry *c = &map->chroms[chr_idx];
            strncpy(loc.chrom_name, c->name, 63);
            loc.chrom_name[63] = '\0';
            loc.local_idx = (loc.global_idx - c->global_start_idx); 
        }
    }
    return loc;
}

CleanGenome* generate_minus_strand(CleanGenome *plus) {
    CleanGenome *minus = (CleanGenome*)malloc(sizeof(CleanGenome));
    minus->length = plus->length;
    minus->kmer_len = plus->kmer_len;
    minus->map = plus->map; 
    minus->data = (char*)malloc(plus->length + 1);
    
    for (size_t i = 0; i < plus->length; i++) {
        char base = plus->data[plus->length - 1 - i];
        char comp = 'N';
        switch (base) {
            case 'A': comp = 'T'; break;
            case 'T': comp = 'A'; break;
            case 'C': comp = 'G'; break;
            case 'G': comp = 'C'; break;
            default:  comp = 'N'; break;
        }
        minus->data[i] = comp;
    }
    minus->data[plus->length] = '\0';
    return minus;
}

// ВИПРАВЛЕНА: Пропускає заголовки FASTA
char* load_clean_target(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("Error opening target file"); return NULL; }
    
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);
    
    char *clean = (char*)malloc(fsize + 1);
    size_t idx = 0;
    
    char line[BUFFER_SIZE];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') continue; // SKIP HEADER
        
        for (int i = 0; line[i]; i++) {
            char c = toupper(line[i]);
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                clean[idx++] = c;
            }
        }
    }
    clean[idx] = '\0';
    fclose(fp);
    
    if (idx == 0) {
        fprintf(stderr, "Error: Target sequence is empty or invalid.\n");
        free(clean);
        return NULL;
    }
    
    char *final_clean = realloc(clean, idx + 1);
    if (final_clean) clean = final_clean;
    
    return clean;
}

TargetHitInfo exclude_target(CleanGenome *plus, CleanGenome *minus, const char *target_seq) {
    TargetHitInfo info;
    memset(&info, 0, sizeof(TargetHitInfo));
    info.found = 0;

    size_t t_len = strlen(target_seq);
    char *found = strstr(plus->data, target_seq);
    char *found_strand_name = "PLUS";
    size_t found_clean_idx = 0;
    
    // 1. Search PLUS
    if (found) {
        found_clean_idx = found - plus->data;
    } else if (minus) {
        // 2. Search MINUS
        found = strstr(minus->data, target_seq);
        if (found) {
            size_t found_minus_idx = found - minus->data;
            found_clean_idx = plus->length - (found_minus_idx + t_len);
            found_strand_name = "MINUS";
        }
    }

    if (found) {
        info.found = 1;
        strcpy(info.strand, found_strand_name);
        
        GenomicLocation loc = get_genomic_location(plus->map, found_clean_idx);
        strcpy(info.chrom_name, loc.chrom_name);
        info.global_start = loc.global_idx;
        info.global_end = loc.global_idx + t_len - 1;
        info.local_start = loc.local_idx;
        info.local_end = loc.local_idx + t_len - 1;

        printf(" -> Target found on %s strand.\n", info.strand);
        printf("    Position: %s:%zu-%zu\n", info.chrom_name, info.local_start, info.local_end);

        size_t mask_start = (found_clean_idx > t_len) ? (found_clean_idx - t_len) : 0;
        size_t mask_end = found_clean_idx + t_len + t_len;
        if (mask_end > plus->length) mask_end = plus->length;
        
        for (size_t i = mask_start; i < mask_end; i++) plus->data[i] = 'N';
        
        if (minus) {
            size_t m_mask_start = plus->length - mask_end;
            size_t m_mask_end = plus->length - mask_start;
            if (m_mask_end > minus->length) m_mask_end = minus->length;
            
            for (size_t i = m_mask_start; i < m_mask_end; i++) minus->data[i] = 'N';
        }
        printf(" -> Region masked (+/- %zu bp).\n", t_len);
    } else {
        printf(" -> Target NOT found in genome.\n");
    }

    return info;
}
