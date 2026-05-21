#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

#define LINES_PER_CHUNK 400000 // 100,000 reads por bloque (aprox 20MB - 50MB)

void filtrar_secuencias_fastq_omp(const char* input_fastq, const char* output_fastq, size_t min_len, size_t max_len) {
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

    // Escaneo ultrarrápido SIMD para demarcar bloques exactos sin corromper reads
    while (ptr < end) {
        ptr = memchr(ptr, '\n', end - ptr);
        if (ptr) {
            ptr++; // Avanzar el salto de línea
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
    // Manejar el último bloque residual
    if (line_count > 0 || chunk_starts[chunk_count] < end) {
        chunk_ends[chunk_count] = end;
        chunk_count++;
    }

    // --- FASE 2 & 3: PROCESAMIENTO PARALELO EN LOTES ---
    // Usamos lotes para limitar el uso de RAM del nodo SLURM
    int num_threads = omp_get_max_threads();
    int batch_size = num_threads * 4; // Lote dinámico basado en cores disponibles
    
    char **out_buffers = malloc(batch_size * sizeof(char*));
    size_t *out_sizes  = malloc(batch_size * sizeof(size_t));

    for (int b = 0; b < chunk_count; b += batch_size) {
        int current_batch = (b + batch_size < chunk_count) ? batch_size : (chunk_count - b);

        // Procesamiento asíncrono de bloques
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < current_batch; i++) {
            int global_i = b + i;
            char *c_start = chunk_starts[global_i];
            char *c_end = chunk_ends[global_i];
            size_t max_out = c_end - c_start;

            out_buffers[i] = malloc(max_out); // Buffer local para no cruzar hilos
            char *out_ptr = out_buffers[i];
            char *p = c_start;

            while (p < c_end) {
                char *rec_start = p;
                
                // Extraer 4 líneas
                char *nl1 = memchr(p, '\n', c_end - p); if (!nl1) break; p = nl1 + 1;
                char *nl2 = memchr(p, '\n', c_end - p); if (!nl2) break;
                size_t seq_len = nl2 - p; p = nl2 + 1;
                char *nl3 = memchr(p, '\n', c_end - p); if (!nl3) break; p = nl3 + 1;
                char *nl4 = memchr(p, '\n', c_end - p); 
                
                if (!nl4) p = c_end; else p = nl4 + 1;

                // Filtrar y empaquetar en memoria temporal
                if (seq_len >= min_len && seq_len <= max_len) {
                    size_t rec_len = p - rec_start;
                    memcpy(out_ptr, rec_start, rec_len);
                    out_ptr += rec_len;
                }
            }
            out_sizes[i] = out_ptr - out_buffers[i];
        }

        // Escritura Secuencial Estricta (Garantiza el orden original)
        for (int i = 0; i < current_batch; i++) {
            if (out_sizes[i] > 0) {
                fwrite(out_buffers[i], 1, out_sizes[i], out_fp);
            }
            free(out_buffers[i]); // Liberar RAM del nodo
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
    const char *output_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtradas_omp.fastq";
    
    size_t min_len = 300;
    size_t max_len = 1200;

    int threads = omp_get_max_threads();
    printf("Iniciando filtrado FASTQ...\n");
    printf("- Hilos OpenMP detectados: %d\n", threads);
    
    double start_time = omp_get_wtime();
    filtrar_secuencias_fastq_omp(input_fastq, output_fastq, min_len, max_len);
    double end_time = omp_get_wtime();

    printf("Filtrado completado con éxito en %f segundos.\n", end_time - start_time);

    return 0;
}