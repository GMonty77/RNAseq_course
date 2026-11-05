---Readme---
# RNA-Seq Analysis: Effect of *Toxoplasma gondii* Infection on Mouse Blood

## Project Overview
This repository contains the scripts and documentation for the RNA-seq data analysis workflow conducted on the IBU HPC cluster at the University of Bern. The goal of this study is to investigate differential gene expression and enriched biological processes in mouse (*Mus musculus*) blood samples following infection with *Toxoplasma gondii* (Infected vs. Wildtype Control).

## Directory Structure
The project on the IBU cluster (`/data/users/$USER/rnaseq_course`) is organized as follows:

```text
rnaseq_course/
├── config/          # Configuration files
├── logs/            # SLURM standard output (.out) and error (.err) logs (ignored by Git)
├── raw_data/        # Symlink to shared FASTQ files (/data/courses/rnaseq_course/toxoplasma_de/reads_Blood)
├── resources/       # Reference genome (FASTA), annotation (GTF), and HISAT2 indices (ignored by Git)
├── results/         # Output directories for each analysis step (ignored by Git)
│   ├── fastqc/      # Individual FastQC HTML reports and ZIP archives
│   └── multiqc/     # Aggregated MultiQC summary report
├── scripts/         # Numbered SLURM and R scripts for reproducible execution
├── .gitignore       # Excludes large data files, logs, and results from version control
└── README.md        # Project documentation
