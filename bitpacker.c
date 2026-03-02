#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "bitpacker.h"

// A=00, C=01, G=10, T=11
// N (та інші) = 00 (A) - безпечно, бо ми їх фільтрували
uint8_t encode_base(char c) {
    switch (toupper(c)) {
        case 'A': return 0; // 00
        case 'C': return 1; // 01
        case 'G': return 2; // 10
        case 'T': return 3; // 11
        default:  return 0; // Treat N as A
    }
}

char decode_base(uint8_t b) {
    switch (b & 3) { // Беремо тільки останні 2 біти
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return '?';
    }
}

PackedGenome* pack_genome_sequence(const char *raw_seq, size_t length) {
    PackedGenome *pg = (PackedGenome*)malloc(sizeof(PackedGenome));
    pg->total_bases = length;
    
    // 1 unit (32 bits) = 16 bases
    // Size = ceil(length / 16)
    pg->data_len = (length + 15) / 16;
    pg->data = (packed_unit_t*)calloc(pg->data_len, sizeof(packed_unit_t));
    
    if (!pg->data) return NULL;

    for (size_t i = 0; i < length; i++) {
        uint8_t code = encode_base(raw_seq[i]);
        size_t unit_idx = i / 16;
        size_t bit_offset = (i % 16) * 2; // 2 біти на базу
        
        // Записуємо 2 біти у відповідну позицію
        // Використовуємо OR
        pg->data[unit_idx] |= ((packed_unit_t)code << bit_offset);
    }
    
    return pg;
}

// Пакує короткий шматок (k-mer) у масив uint64
void pack_sequence_chunk(const char *seq, uint64_t *out_buffer, int k) {
    // Очищуємо буфер (припускаємо 2 x uint64, що вистачить на 64 нуклеотиди)
    out_buffer[0] = 0;
    out_buffer[1] = 0;
    
    for (int i = 0; i < k; i++) {
        uint8_t code = encode_base(seq[i]);
        
        if (i < 32) {
            // Перші 32 бази -> out_buffer[0]
            out_buffer[0] |= ((uint64_t)code << (i * 2));
        } else if (i < 64) {
            // Наступні 32 бази -> out_buffer[1]
            out_buffer[1] |= ((uint64_t)code << ((i - 32) * 2));
        }
    }
}

void unpack_sequence(const uint64_t *packed, char *out_str, int k) {
    for (int i = 0; i < k; i++) {
        uint8_t code;
        if (i < 32) {
            code = (packed[0] >> (i * 2)) & 3;
        } else {
            code = (packed[1] >> ((i - 32) * 2)) & 3;
        }
        out_str[i] = decode_base(code);
    }
    out_str[k] = '\0';
}

void free_packed_genome(PackedGenome *pg) {
    if (pg) {
        free(pg->data);
        free(pg);
    }
}
