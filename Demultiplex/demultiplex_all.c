#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LENGTH 2048
#define PATTERN_LENGTH 200

// Función para verificar si un patrón está en una secuencia y devolver la posición
int find_pattern_in_region(const char *sequence, const char *pattern) {
    const char *pos = strstr(sequence, pattern);
    if (pos != NULL) {
        return pos - sequence;  // Retornar el índice donde comienza el patrón
    }
    return -1;  // Patrón no encontrado
}

// Filtra lecturas por dos patrones
void filter_reads_by_two_patterns(const char *fastq_file, const char *pattern1, const char *pattern2, const char *output_file) {
    FILE *input = fopen(fastq_file, "r");
    FILE *output = fopen(output_file, "w");

    if (!input || !output) {
        fprintf(stderr, "Error abriendo archivo\n");
        exit(1);
    }

    char line[MAX_SEQ_LENGTH];
    char sequence[MAX_SEQ_LENGTH];

    while (fgets(line, sizeof(line), input)) {
        // Leer la secuencia (asume que está en la línea correcta de un archivo FASTQ)
        fgets(sequence, sizeof(sequence), input);

        // Comprobar patrones en las primeras y últimas 100 bases
        char first_bases[PATTERN_LENGTH + 1] = {0};
        char last_bases[PATTERN_LENGTH + 1] = {0};
        strncpy(first_bases, sequence, PATTERN_LENGTH);
        strncpy(last_bases, sequence + strlen(sequence) - PATTERN_LENGTH - 1, PATTERN_LENGTH);

        int pattern1_found = find_pattern_in_region(first_bases, pattern1);
        int pattern2_found = find_pattern_in_region(last_bases, pattern2);

        // Si se encuentra uno de los patrones, escribir en el archivo de salida
        if (pattern1_found!=-1 || pattern2_found!=-1) {
            fprintf(stdout, "%s",line);
            fprintf(stdout, "%d,%d\n",pattern1_found,pattern2_found);
            fputs(line, output);  // Escribir la línea de encabezado
            fputs(sequence, output);  // Escribir la secuencia
            fgets(line, sizeof(line), input);  // Saltar la línea +
            fputs(line, output);
            fgets(line, sizeof(line), input);  // Saltar la línea de calidad
            fputs(line, output);
        }
    }

    fclose(input);
    fclose(output);
}

// Procesa patrones para cada paciente
void process_patient_patterns(const char *pattern_file, const char *fastq_file) {
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
	fprintf(stdout, "index %d\n",i);
        sprintf(output_file, "index_%d_filtered_sequences.fastq", i);
        filter_reads_by_two_patterns(fastq_file, pattern1, pattern2, output_file);
        i++;
    }

    fclose(pf);
}

int main() {
    const char *fastq_file = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/complete_filtradas.fastq";  // Ruta al archivo FASTQ
    const char *pattern_file = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/2.demultiplexacion/forward/lista_barcodes.txt";  // Ruta al archivo de patrones

    // Llamar a la función para procesar patrones para todos los pacientes
    process_patient_patterns(pattern_file, fastq_file);

    return 0;
}
