#ifndef GENOME_H
#define GENOME_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "bitpacker.h" // Потрібно для PackedGenome

// ШАР 1: Відрізки чистого ДНК
typedef struct {
    size_t clean_offset;    
    size_t original_offset; 
    size_t length;          
} GenomeSegment;

// ШАР 2: Хромосоми
typedef struct {
    char name[64];           
    size_t global_start_idx; 
    size_t length;           
} ChromosomeEntry;

// Головна карта
typedef struct {
    GenomeSegment *segments;
    size_t count;
    size_t capacity;
    size_t total_original_length; 

    ChromosomeEntry *chroms;
    size_t chrom_count;
    size_t chrom_capacity;
} GenomeMap;

typedef struct {
    char *data;         
    size_t length;      
    int kmer_len;
    GenomeMap *map;     
} CleanGenome;

typedef struct {
    size_t global_idx;
    char chrom_name[64];
    size_t local_idx; 
} GenomicLocation;

// Структура для збереження інфо про знайдену ціль (для CSV та візуалізації)
typedef struct {
    int found;              // 1 якщо знайдено
    char strand[8];         // "PLUS" або "MINUS"
    size_t global_start;    // Глобальний індекс початку
    size_t global_end;      // Глобальний індекс кінця
    char chrom_name[64];    
    size_t local_start;
    size_t local_end;
} TargetHitInfo;

// Основні функції
CleanGenome* load_and_clean_genome(const char *filename, int k);
void free_genome(CleanGenome *genome);
CleanGenome* generate_minus_strand(CleanGenome *plus_strand);
char* load_clean_target(const char *filename);
GenomicLocation get_genomic_location(GenomeMap *map, size_t clean_idx);

// Модифікована функція exclude_target тепер повертає інфо про знайдену ціль
TargetHitInfo exclude_target(CleanGenome *plus, CleanGenome *minus, const char *target_seq);

// --- НОВІ ФУНКЦІЇ (2-BIT CACHING) ---

// Зберігає вже ЗАПАКОВАНІ (2-bit) геноми + мапу + інфо про ціль
int save_packed_cache(const char *filename, PackedGenome *pg_plus, PackedGenome *pg_minus, GenomeMap *map, int k, TargetHitInfo target_info);

// Завантажує все назад. Повертає 1 якщо успішно.
int load_packed_cache(const char *filename, PackedGenome **pg_plus, PackedGenome **pg_minus, GenomeMap **map, int *k, TargetHitInfo *target_info);

// Експорт карти у CSV (теперь включає info про ціль)
int export_genome_map_csv(const char *filename, GenomeMap *map, TargetHitInfo target_info);

#endif
