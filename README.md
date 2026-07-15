# README

## Project Description

This project aimed to optimize the preprocessing of sequencing data generated using Oxford Nanopore Technologies (ONT). Specifically, two critical stages of the preprocessing pipeline were improved:

1. **Filtering sequencing reads by length**
2. **Demultiplexing reads for each patient/sample**

The proposed implementations focus on improving execution performance through optimized sequential and parallel approaches using Python and C, including the use of `mmap`-based memory access strategies for efficient FASTQ file processing.

---

## Data Availability

The datasets used in this study cannot be publicly shared because they correspond to sensitive patient information and are therefore subject to privacy and ethical restrictions.

---

# Compilation Instructions

## Parallel Version

```bash
gcc -O3 demultiplex_all_mmapSec.c -lpthread -o demux
```

## Sequential Version

```bash
gcc -O3 demultiplex_all_mmapSec.c -o demux
```

## CPU Profiling

To generate a CPU execution profile, compile the C implementation with profiling support enabled:

```bash
gcc -pg -fopenmp -O2 -o demultiplex_all_profile demultiplex_all.c
```

For the Python implementation, execute the script with cProfile:
```python
python -m cProfile -o perfil_cpu.prof filter_longitud.py
```
---

# Repository Structure

## Filter Folder

This folder contains the implementations related to FASTQ read-length filtering.

| File | Description |
|---|---|
| `filter_longitud.py` | Sequential Baseline Biopython Filter |
| `filtrar_fastq_mmap.py` | Sequential Python `mmap` Filter |
| `filter_longitud_sec1.c` | Sequential C `mmap` Filter |
| `filter_longitud_par1.c` | Parallel C `mmap` Filter |
| `filter_longitud_maleable.c` | Parallel C `mmap` Malleable Filter |
| `job.fastp.sh` | Slurm job using fastp as reference parallel filter program |

---

## Demultiplex Folder

This folder contains the implementations related to read demultiplexing.

| File | Description |
|---|---|
| `demultiplex_all.c` | C Sequential Baseline Demultiplex |
| `demultiplex_all_mmapSec.c` | C Sequential Demultiplex using `mmap` |
| `demultiplex_all_mmapParalelo.c` | C Parallel Demultiplex using `mmap` |

---

# Objective

The main objective of this repository is to provide optimized preprocessing strategies for ONT sequencing data, reducing computational overhead and improving scalability when handling large FASTQ datasets. The implementations explore both sequential and parallel paradigms, as well as low-level memory optimization techniques for high-performance bioinformatics workflows.
