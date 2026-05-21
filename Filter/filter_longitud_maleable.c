#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

#define LINES_PER_CHUNK 400000 // 100,000 reads por bloque
#define WINDOW_SIZE 128        // Macro-lote independiente del número de hilos

void filtrar_secuencias_fastq_maleable(const char* input_fastq, const char* output_fastq, size_t min_len, size_t max_len) {
    int fd_in = open(input_fastq, O_RDONLY);
    if (fd_in < 0) { perror("Error abriendo input"); return; }

    struct stat st;
    if (fstat(fd_in, &st) < 0) { perror("Error stat"); close(fd_in); return; }
    size_t file_size = st.st_size;

    char *data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd_in, 0);
    if (data == MAP_FAILED) { perror("Error mmap"); close(fd_in); return; }

    FILE *out_fp = fopen(output_fastq, "w");
    if (!out_fp) { perror("Error abriendo output"); munmap(data, file_size); close(fd_in); return; }

    // --- FASE 1: INDEXACIÓN ESTRUCTURAL ---
    size_t capacity = 1024;
    char **chunk_starts = malloc(capacity * sizeof(char*));
    char **chunk_ends   = malloc(capacity * sizeof(char*));
    int chunk_count = 0;
    char *ptr = data;
    char *end = data + file_size;
    chunk_starts[0] = ptr;
    size_t line_count = 0;

    // Escaneo secuencial ultrarrápido (simil a fase Map)
    while (ptr < end) {
        ptr = memchr(ptr, '\n', end - ptr);
        if (ptr) {
            ptr++;
            line_count++;
            if (line_count == LINES_PER_CHUNK) {
                chunk_ends[chunk_count] = ptr;
                chunk_count++;
                if (chunk_count >= capacity) {
                    capacity *= 2;
                    chunk_starts = realloc(chunk_starts, capacity * sizeof(char*));
                    chunk_ends   = realloc(chunk_ends, capacity * sizeof(char*));
                }
                chunk_starts[chunk_count] = ptr;
                line_count = 0;
            }
        } else {
            break;
        }
    }
    if (line_count > 0 || chunk_starts[chunk_count] < end) {
        chunk_ends[chunk_count] = end;
        chunk_count++;
    }

    // --- FASE 2: PROCESAMIENTO MALEABLE ---
    // Activamos la maleabilidad a nivel de sistema (Runtime de OpenMP)
    omp_set_dynamic(1); 

    // Los buffers dependen de la cantidad de trabajo (WINDOW_SIZE), NO de los hilos.
    // Esto permite que el SO reduzca o aumente los hilos en tiempo real sin fallos.
    char **out_buffers = malloc(WINDOW_SIZE * sizeof(char*));
    size_t *out_sizes  = malloc(WINDOW_SIZE * sizeof(size_t));

    for (int b = 0; b < chunk_count; b += WINDOW_SIZE) {
        int current_window = (b + WINDOW_SIZE < chunk_count) ? WINDOW_SIZE : (chunk_count - b);

        // schedule(dynamic, 1) asigna los chunks 1 a 1 al hilo que esté libre.
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < current_window; i++) {
            int global_i = b + i;
            char *c_start = chunk_starts[global_i];
            char *c_end = chunk_ends[global_i];
            size_t max_out = c_end - c_start;

            // Cada iteración maneja su propia memoria. Cero colisiones.
            char *local_buf = malloc(max_out); 
            char *out_ptr = local_buf;
            char *p = c_start;

            while (p < c_end) {
                char *rec_start = p;
                
                char *nl1 = memchr(p, '\n', c_end - p); if (!nl1) break; p = nl1 + 1;
                char *nl2 = memchr(p, '\n', c_end - p); if (!nl2) break;
                size_t seq_len = nl2 - p; p = nl2 + 1;
                char *nl3 = memchr(p, '\n', c_end - p); if (!nl3) break; p = nl3 + 1;
                char *nl4 = memchr(p, '\n', c_end - p); 
                
                if (!nl4) p = c_end; else p = nl4 + 1;

                if (seq_len >= min_len && seq_len <= max_len) {
                    size_t rec_len = p - rec_start;
                    memcpy(out_ptr, rec_start, rec_len);
                    out_ptr += rec_len;
                }
            }
            // Guardamos el puntero local en el array global de la ventana
            out_buffers[i] = local_buf;
            out_sizes[i] = out_ptr - local_buf;
        }

        // --- FASE 3: ESCRITURA ESTRICTA (Reduce) ---
        // Volcado secuencial para garantizar el orden byte a byte.
        for (int i = 0; i < current_window; i++) {
            if (out_sizes[i] > 0) {
                fwrite(out_buffers[i], 1, out_sizes[i], out_fp);
            }
            free(out_buffers[i]); // Liberación inmediata de RAM
        }
    }

    free(out_buffers);
    free(out_sizes);
    free(chunk_starts);
    free(chunk_ends);
    fclose(out_fp);
    munmap(data, file_size);
    close(fd_in);
}

int main() {
    const char *input_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete.fastq";
    const char *output_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtradas_maleable.fastq";
    
    size_t min_len = 300;
    size_t max_len = 1200;

    printf("Iniciando filtrado FASTQ con Arquitectura Maleable...\n");
    
    double start_time = omp_get_wtime();
    filtrar_secuencias_fastq_maleable(input_fastq, output_fastq, min_len, max_len);
    double end_time = omp_get_wtime();

    printf("Filtrado completado. Tiempo: %f segundos.\n", end_time - start_time);

    return 0;
}