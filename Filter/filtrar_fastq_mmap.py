import mmap


def filtrar_fastq_mmap(input_fastq, output_fastq, min_len=300, max_len=1200):
    """
    Versión optimizada usando mmap + parsing manual FASTQ
    Mucho más rápida que Biopython para filtrado por longitud.
    """

    with open(input_fastq, "rb") as f_in, open(output_fastq, "wb") as f_out:
        mm = mmap.mmap(f_in.fileno(), 0, access=mmap.ACCESS_READ)

        write_buffer = []
        buffer_limit = 10000  # número de reads antes de escribir

        i = 0
        size = mm.size()

        while i < size:
            # Leer 4 líneas FASTQ
            l1 = mm.readline()
            if not l1:
                break
            l2 = mm.readline()
            l3 = mm.readline()
            l4 = mm.readline()

            # Longitud de la secuencia (sin \n)
            seq_len = len(l2.rstrip(b"\n"))

            if min_len <= seq_len <= max_len:
                write_buffer.append(l1)
                write_buffer.append(l2)
                write_buffer.append(l3)
                write_buffer.append(l4)

            # Flush buffer (reduce syscalls)
            if len(write_buffer) >= buffer_limit * 4:
                f_out.writelines(write_buffer)
                write_buffer.clear()

            i = mm.tell()

        # Escribir lo restante
        if write_buffer:
            f_out.writelines(write_buffer)

        mm.close()


# Ejemplo de uso
input_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete.fastq"
output_fastq = "/Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtradas_mmap.fastq"

filtrar_fastq_mmap(input_fastq, output_fastq)
