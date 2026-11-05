#!/usr/bin/env bash
#SBATCH --job-name=multiqc
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Stop script immediately if any command fails
set +u
set -euo pipefail

# ---- Settings & Variables ----
WORKDIR=/data/users/$USER/rnaseq_course
INPUT_DIR=$WORKDIR/results
OUT_DIR=$WORKDIR/results/multiqc_fastp
CONTAINER=/containers/apptainer/multiqc-1.33.sif

# Ensure log and output directories exist
mkdir -p logs "$OUT_DIR"

# ---- Run MultiQC ----
echo "Starting MultiQC aggregation..."
apptainer exec --bind /data "$CONTAINER" multiqc \
    "$INPUT_DIR" \
    -o "$OUT_DIR"
echo "MultiQC report successfully generated in $OUT_DIR"
