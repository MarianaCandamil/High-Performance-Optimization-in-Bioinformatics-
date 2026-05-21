#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

void filtrar_secuencias_fastq_mmap(const char* input_fastq, const char* output_fastq, size_t min_len, size_t max_len) {
    // 1. Abrir archivo de entrada directamente con llamadas al sistema (POSIX)
    int fd_in = open(input_fastq, O_RDONLY);
    if (fd_in < 0) {
        perror("Error fatal: No se pudo abrir el archivo FASTQ de entrada");
        return;
    }

    // Obtener el tamaño exacto del archivo
    struct stat st;
    if (fstat(fd_in, &st) < 0) {
        perror("Error fatal: No se pudo obtener el tamaño del archivo");
        close(fd_in);
        return;
    }
    size_t file_size = st.st_size;

    // 2. Mapear el archivo entero en memoria (Zero-copy read)
    // PROT_READ: Solo lectura. MAP_PRIVATE: Los cambios no se escriben en el disco.
    char *data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd_in, 0);
    if (data == MAP_FAILED) {
        perror("Error fatal: Fallo al ejecutar mmap");
        close(fd_in);
        return;
    }

    // 3. Abrir archivo de salida usando buffers de alto rendimiento de C
    FILE *out_fp = fopen(output_fastq, "w");
    if (!out_fp) {
        perror("Error fatal: No se pudo abrir el archivo de salida");
        munmap(data, file_size);
        close(fd_in);
        return;
    }

    // Punteros para recorrer la memoria mapeada
    char *ptr = data;
    char *end = data + file_size;

    // 4. Bucle principal de procesamiento de memoria
    while (ptr < end) {
        char *start_record = ptr; // Guardamos el inicio de las 4 líneas

        // Línea 1: Header (@...)
        // memchr usa instrucciones vectorizadas (AVX/SSE) internamente en glibc
        char *nl1 = memchr(ptr, '\n', end - ptr);
        if (!nl1) break; // Archivo truncado o fin de archivo

        // Línea 2: Secuencia (La que nos importa para la longitud)
        ptr = nl1 + 1;
        char *nl2 = memchr(ptr, '\n', end - ptr);
        if (!nl2) break;
        
        // Calculamos la longitud usando aritmética de punteros
        size_t seq_len = nl2 - ptr;

        // Línea 3: Separador (+)
        ptr = nl2 + 1;
        char *nl3 = memchr(ptr, '\n', end - ptr);
        if (!nl3) break;

        // Línea 4: Calidad
        ptr = nl3 + 1;
        char *nl4 = memchr(ptr, '\n', end - ptr);
        
        if (!nl4) {
            // Manejo de caso borde: Última línea del archivo sin salto de línea al final
            ptr = end; 
        } else {
            ptr = nl4 + 1; // Avanzamos al inicio del siguiente record
        }

        // 5. Lógica de filtrado y escritura en bloque
        if (seq_len >= min_len && seq_len <= max_len) {
            // ¡La magia! En lugar de escribir 4 líneas por separado, escribimos
            // el bloque entero de memoria de una sola vez.
            size_t block_size = ptr - start_record;
            fwrite(start_record, 1, block_size, out_fp);
        }
    }

    // 6. Limpieza y liberación de recursos
    fclose(out_fp);
    munmap(data, file_size);
    close(fd_in);
}

int main() {
    // Tus rutas originales
    const char *input_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete.fastq";
    const char *output_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtradas_secuencial_time.fastq";
    
    // Parámetros
    size_t min_len = 300;
    size_t max_len = 1200;

    printf("Iniciando filtrado FASTQ (Nivel Avanzado: mmap + Vectorización)...\n");
    filtrar_secuencias_fastq_mmap(input_fastq, output_fastq, min_len, max_len);
    printf("Filtrado completado con éxito.\n");

    return 0;
}