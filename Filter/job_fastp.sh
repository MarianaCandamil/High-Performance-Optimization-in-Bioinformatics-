#!/bin/bash
#SBATCH -J Fastp
#SBATCH -D .
#SBATCH -e QC_fastp_%j.err
#SBATCH -o QC_fastp_%j.out
#SBATCH -n 27
##SBATCH --partition=fat
##SBATCH --nodelist=ibcu12


source /Tayra-Share/miniconda/bin/activate
conda activate preprocessing


# Process with fastp
fastp -i /Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete.fastq --thread 27 --length_required 300 --length_limit 1200 -o /Tayra-Share/home/Biotec/Projects/Proyectos_Ejecucion/202202PE_ORIGENCALDAS/LIB_C_Complete/1.filter/complete_filtered_fastp2.fastq


echo "Finished"
