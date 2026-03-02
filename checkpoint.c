#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "checkpoint.h"

int save_checkpoint(const char *filename, size_t offset, size_t total_processed, int k, size_t num_tasks, KmerTask *tasks) {
    // Спочатку пишемо у тимчасовий файл, щоб не пошкодити основний при збої
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", filename);

    FILE *fp = fopen(temp_filename, "wb");
    if (!fp) {
        perror("Failed to create checkpoint file");
        return 0;
    }

    CheckpointHeader header;
    header.offset = offset;
    header.total_processed = total_processed;
    header.k = k;
    header.num_tasks = num_tasks;
    header.timestamp = (uint64_t)time(NULL);

    // 1. Header
    fwrite(&header, sizeof(CheckpointHeader), 1, fp);

    // 2. Tasks data (Top-100 lists mostly)
    fwrite(tasks, sizeof(KmerTask), num_tasks, fp);

    fclose(fp);

    // Атомарна заміна файлу
    if (rename(temp_filename, filename) != 0) {
        perror("Failed to commit checkpoint");
        return 0;
    }

    return 1;
}

int has_checkpoint(const char *filename) {
    return (access(filename, F_OK) != -1);
}

int load_checkpoint(const char *filename, size_t *offset, size_t *total_processed, int k, size_t num_tasks, KmerTask *tasks) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;

    CheckpointHeader header;
    if (fread(&header, sizeof(CheckpointHeader), 1, fp) != 1) {
        fclose(fp);
        return 0;
    }

    // Валідація: чи підходить цей сейв до поточного запуску?
    if (header.k != k || header.num_tasks != num_tasks) {
        fprintf(stderr, "Checkpoint mismatch! Saved K=%d, Tasks=%zu. Current K=%d, Tasks=%zu.\n", 
                header.k, header.num_tasks, k, num_tasks);
        fclose(fp);
        return -1; // -1 означає "файл є, але не підходить"
    }

    *offset = header.offset;
    *total_processed = header.total_processed;

    // Читаємо масив завдань
    if (fread(tasks, sizeof(KmerTask), num_tasks, fp) != num_tasks) {
        fprintf(stderr, "Checkpoint file truncated or corrupted.\n");
        fclose(fp);
        return 0;
    }

    fclose(fp);
    return 1;
}
