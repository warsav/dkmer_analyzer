#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "bitpacker.h"

// Структура заголовка чекпоінту
typedef struct {
    size_t offset;          // Де зупинились
    size_t total_processed; // Скільки пройшли (для UI)
    int k;                  // Параметр K (щоб не переплутати запуски)
    size_t num_tasks;       // Кількість k-мерів
    uint64_t timestamp;     // Час збереження
} CheckpointHeader;

// Збереження стану
int save_checkpoint(const char *filename, size_t offset, size_t total_processed, int k, size_t num_tasks, KmerTask *tasks);

// Завантаження стану
int load_checkpoint(const char *filename, size_t *offset, size_t *total_processed, int k, size_t num_tasks, KmerTask *tasks);

// Перевірка наявності валідного чекпоінту
int has_checkpoint(const char *filename);

#endif
