from Bio import SeqIO

def filtrar_secuencias_fastq(input_fastq, output_fastq, min_len=300, max_len=1200):
    """
    Filtra las secuencias FASTQ basadas en la longitud y las guarda en un nuevo archivo.
    
    Parámetros:
        input_fastq (str): Ruta del archivo FASTQ de entrada.
        output_fastq (str): Ruta del archivo FASTQ de salida con las secuencias filtradas.
        min_len (int): Longitud mínima de las secuencias a filtrar.
        max_len (int): Longitud máxima de las secuencias a filtrar.
    """
    with open(input_fastq, "r") as input_handle, open(output_fastq, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            seq_len = len(record.seq)
            if min_len <= seq_len <= max_len:
                SeqIO.write(record, output_handle, "fastq")

# Ejemplo de uso
input_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete.fastq"  # Archivo FASTQ de entrada
output_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtradas_secuencial_time.fastq"  # Archivo FASTQ de salida
filtrar_secuencias_fastq(input_fastq, output_fastq)
