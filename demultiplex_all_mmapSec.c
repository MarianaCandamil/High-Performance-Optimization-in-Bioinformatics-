#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define MAX_SEQ_LENGTH 2048
#define PATTERN_LENGTH 200

// Buscar patrón en región
int find_pattern_in_region(const char *sequence, const char *pattern) {
    const char *pos = strstr(sequence, pattern);
    if (pos != NULL) {
        return pos - sequence;
    }
    return -1;
}

// Mapear archivo FASTQ en memoria
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

// Filtrar usando mmap (misma lógica que tu versión)
void filter_reads_mmap(const char *data, size_t size,
                      const char *pattern1, const char *pattern2,
                      const char *output_file) {

    FILE *output = fopen(output_file, "w");
    if (!output) {
        perror("fopen output");
        exit(1);
    }

    const char *ptr = data;
    const char *end = data + size;

    while (ptr < end) {

        // Línea 1 (header)
        const char *line1 = ptr;
        const char *nl1 = memchr(ptr, '\n', end - ptr);
        if (!nl1) break;

        // Línea 2 (sequence)
        const char *line2 = nl1 + 1;
        const char *nl2 = memchr(line2, '\n', end - line2);
        if (!nl2) break;

        // Línea 3 (+)
        const char *line3 = nl2 + 1;
        const char *nl3 = memchr(line3, '\n', end - line3);
        if (!nl3) break;

        // Línea 4 (quality)
        const char *line4 = nl3 + 1;
        const char *nl4 = memchr(line4, '\n', end - line4);
        if (!nl4) break;

        // Copiar secuencia
        int seq_len = nl2 - line2;
        if (seq_len <= 0 || seq_len >= MAX_SEQ_LENGTH) {
            ptr = nl4 + 1;
            continue;
        }

        char sequence[MAX_SEQ_LENGTH];
        memcpy(sequence, line2, seq_len);
        sequence[seq_len] = '\0';

        // Regiones
        char first_bases[PATTERN_LENGTH + 1] = {0};
        char last_bases[PATTERN_LENGTH + 1] = {0};

        if (seq_len > PATTERN_LENGTH) {
            strncpy(first_bases, sequence, PATTERN_LENGTH);
            strncpy(last_bases,
                    sequence + seq_len - PATTERN_LENGTH,
                    PATTERN_LENGTH);
        } else {
            strncpy(first_bases, sequence, seq_len);
            strncpy(last_bases, sequence, seq_len);
        }

        int p1 = find_pattern_in_region(first_bases, pattern1);
        int p2 = find_pattern_in_region(last_bases, pattern2);

        if (p1 != -1 || p2 != -1) {
            fwrite(line1, 1, nl1 - line1 + 1, output);
            fwrite(line2, 1, nl2 - line2 + 1, output);
            fwrite(line3, 1, nl3 - line3 + 1, output);
            fwrite(line4, 1, nl4 - line4 + 1, output);
        }

        ptr = nl4 + 1;
    }

    fclose(output);
}

// Procesar patrones (mantiene tus 96 archivos)
void process_patient_patterns(const char *pattern_file, const char *fastq_file) {

    size_t size;
    char *data = map_fastq(fastq_file, &size);

    FILE *pf = fopen(pattern_file, "r");
    if (!pf) {
        fprintf(stderr, "Error abriendo archivo de patrones\n");
        exit(1);
    }

    char pattern1[PATTERN_LENGTH];
    char pattern2[PATTERN_LENGTH];
    int i = 1;

    while (fscanf(pf, "%[^,],%s\n", pattern1, pattern2) != EOF) {

        char output_file[256];
        sprintf(output_file, "index_%d_filtered_sequences.fastq", i);

        fprintf(stdout, "Procesando index %d\n", i);

        filter_reads_mmap(data, size, pattern1, pattern2, output_file);

        i++;
    }

    fclose(pf);
    munmap(data, size);
}

int main() {

    const char *fastq_file = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/complete_filtradas.fastq";

    const char *pattern_file = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/lista_barcodes.txt";

    process_patient_patterns(pattern_file, fastq_file);

    return 0;
}