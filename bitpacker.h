#ifndef BITPACKER_H
#define BITPACKER_H

#include <stdint.h>
#include <stddef.h>

// Максимально можлива кількість (Hard Limit).
// Дефолт в main.c = 100. Юзер може підняти до 1000 через -n.
#define MAX_TOP_RESULTS 1000 

typedef uint32_t packed_unit_t;

typedef struct {
    packed_unit_t *data;
    size_t total_bases;
    size_t data_len;
} PackedGenome;

typedef struct {
    uint64_t sequence_head; 
    uint32_t position_in_genome; 
    float similarity_score;      
    uint64_t packed_seq[2]; 
} OverlapResult;

typedef struct {
    uint64_t packed_pattern[2]; 
    uint32_t id;                
    OverlapResult top_results[MAX_TOP_RESULTS]; 
    uint32_t worst_score_idx; // Використовується як лічильник кількості (count)
} KmerTask;

PackedGenome* pack_genome_sequence(const char *raw_seq, size_t length);
void pack_sequence_chunk(const char *seq, uint64_t *out_buffer, int k);
void free_packed_genome(PackedGenome *pg);
void unpack_sequence(const uint64_t *packed, char *out_str, int k);

#endif
