#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Stop the script if something is wrong
set -euo pipefail
# ---- Variables and Settings ----
WORKDIR=/data/users/$USER/rnaseq_course
READS_DIR=$WORKDIR/raw_data
OUT_DIR=$WORKDIR/results/fastqc
CONTAINER=/containers/apptainer/fastqc-0.12.1.sif
THREADS=$SLURM_CPUS_PER_TASK
# ---- Commands ----
mkdir -p "$OUT_DIR"
apptainer exec --bind /data "$CONTAINER" fastqc \
    -t "$THREADS" \
    -o "$OUT_DIR" \
    "$READS_DIR"/*.fastq.gz
