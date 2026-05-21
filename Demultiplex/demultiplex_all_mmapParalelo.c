#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define MAX_SEQ_LENGTH 2048
#define PATTERN_LENGTH 200

// =====================
// ESTRUCTURA PARA HILOS
// =====================
typedef struct {
    const char *data;
    size_t size;
    char pattern1[PATTERN_LENGTH];
    char pattern2[PATTERN_LENGTH];
    char output_file[256];
    int index;
} ThreadData;

// =====================
// BUSCAR PATRÓN
// =====================
int find_pattern_in_region(const char *sequence, const char *pattern) {
    const char *pos = strstr(sequence, pattern);
    return (pos != NULL) ? (pos - sequence) : -1;
}

// =====================
// MMAP
// =====================
char *map_fastq(const char *filename, size_t *size) {
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        perror("open");
        exit(1);
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        perror("fstat");
        exit(1);
    }

    *size = sb.st_size;

    char *mapped = mmap(NULL, *size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        perror("mmap");
        exit(1);
    }

    close(fd);
    return mapped;
}

// =====================
// FUNCIÓN DEL HILO
// =====================
void *process_index(void *arg) {
    ThreadData *td = (ThreadData *)arg;

    FILE *output = fopen(td->output_file, "w");
    if (!output) {
        perror("fopen");
        pthread_exit(NULL);
    }

    const char *ptr = td->data;
    const char *end = td->data + td->size;

    while (ptr < end) {

        const char *line1 = ptr;
        const char *nl1 = memchr(ptr, '\n', end - ptr);
        if (!nl1) break;

        const char *line2 = nl1 + 1;
        const char *nl2 = memchr(line2, '\n', end - line2);
        if (!nl2) break;

        const char *line3 = nl2 + 1;
        const char *nl3 = memchr(line3, '\n', end - line3);
        if (!nl3) break;

        const char *line4 = nl3 + 1;
        const char *nl4 = memchr(line4, '\n', end - line4);
        if (!nl4) break;

        int seq_len = nl2 - line2;
        if (seq_len <= 0 || seq_len >= MAX_SEQ_LENGTH) {
            ptr = nl4 + 1;
            continue;
        }

        char sequence[MAX_SEQ_LENGTH];
        memcpy(sequence, line2, seq_len);
        sequence[seq_len] = '\0';

        char first_bases[PATTERN_LENGTH + 1] = {0};
        char last_bases[PATTERN_LENGTH + 1] = {0};

        if (seq_len > PATTERN_LENGTH) {
            strncpy(first_bases, sequence, PATTERN_LENGTH);
            strncpy(last_bases, sequence + seq_len - PATTERN_LENGTH, PATTERN_LENGTH);
        } else {
            strncpy(first_bases, sequence, seq_len);
            strncpy(last_bases, sequence, seq_len);
        }

        int p1 = find_pattern_in_region(first_bases, td->pattern1);
        int p2 = find_pattern_in_region(last_bases, td->pattern2);

        if (p1 != -1 || p2 != -1) {
            fwrite(line1, 1, nl1 - line1 + 1, output);
            fwrite(line2, 1, nl2 - line2 + 1, output);
            fwrite(line3, 1, nl3 - line3 + 1, output);
            fwrite(line4, 1, nl4 - line4 + 1, output);
        }

        ptr = nl4 + 1;
    }

    fclose(output);
    pthread_exit(NULL);
}

// =====================
// MAIN PROCESS
// =====================
void process_patient_patterns_parallel(const char *pattern_file, const char *fastq_file) {

    size_t size;
    char *data = map_fastq(fastq_file, &size);

    FILE *pf = fopen(pattern_file, "r");
    if (!pf) {
        fprintf(stderr, "Error abriendo archivo de patrones\n");
        exit(1);
    }

    pthread_t threads[200];   // suficiente para 96
    ThreadData tdata[200];

    int i = 0;

    while (fscanf(pf, "%[^,],%s\n", tdata[i].pattern1, tdata[i].pattern2) != EOF) {

        tdata[i].data = data;
        tdata[i].size = size;
        tdata[i].index = i + 1;

        sprintf(tdata[i].output_file,
                "index_%d_filtered_sequences.fastq", i + 1);

        printf("Creando hilo para index %d\n", i + 1);

        pthread_create(&threads[i], NULL, process_index, &tdata[i]);

        i++;
    }

    // Esperar todos los hilos
    for (int j = 0; j < i; j++) {
        pthread_join(threads[j], NULL);
    }

    munmap(data, size);
    fclose(pf);
}

// =====================
// MAIN
// =====================
int main() {

    const char *fastq_file =
    "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/complete_filtradas.fastq";

    const char *pattern_file =
    "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/lista_barcodes.txt";

    process_patient_patterns_parallel(pattern_file, fastq_file);

    return 0;
}