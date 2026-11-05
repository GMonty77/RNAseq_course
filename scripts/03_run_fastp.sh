#!/usr/bin/env bash
#SBATCH --job-name=fastp
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Allow unbound variables briefly if /etc/bashrc is triggered, then enforce strict mode
set +u
set -euo pipefail

# ---- Settings & Variables ----
WORKDIR=/data/users/$USER/rnaseq_course
READS_DIR=$WORKDIR/raw_data
OUT_DIR=$WORKDIR/results/fastp
CONTAINER=/containers/apptainer/fastp_0.24.1.sif
THREADS=$SLURM_CPUS_PER_TASK

# Ensure log and output directories exist
mkdir -p logs "$OUT_DIR"

# ---- Run fastp on Paired-End Reads ----
echo "Starting fastp trimming and quality filtering..."

for R1 in "$READS_DIR"/*_1.fastq.gz; do
    # Define corresponding Mate 2 file and sample base name
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    echo "Processing sample: $SAMPLE"

    apptainer exec --bind /data "$CONTAINER" fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT_DIR/${SAMPLE}_1.trimmed.fastq.gz" \
        -O "$OUT_DIR/${SAMPLE}_2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "$OUT_DIR/${SAMPLE}_fastp.html" \
        --json "$OUT_DIR/${SAMPLE}_fastp.json"
done

echo "fastp processing complete! Results saved in $OUT_DIR"
