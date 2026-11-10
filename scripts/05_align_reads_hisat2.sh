#!/usr/bin/env bash
#SBATCH --job-name=hisat2_align
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Allow unbound variables briefly if /etc/bashrc is triggered, then enforce strict mode
set +u
set -euo pipefail

# ---- Settings & Variables ----
WORKDIR=/data/users/$USER/rnaseq_course
INDEX_PREFIX=$WORKDIR/resources/hisat2_index/GRCm39
OUT_DIR=$WORKDIR/results/hisat2
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
THREADS=$SLURM_CPUS_PER_TASK

# Choose input directory:
# Option A: Use fastp trimmed reads (if you ran the fastp step)
READS_DIR=$WORKDIR/results/fastp
R1_PATTERN="*_1.trimmed.fastq.gz"

# Option B: Use raw reads directly (uncomment the two lines below if NOT using fastp)
# READS_DIR=$WORKDIR/raw_data
# R1_PATTERN="*_1.fastq.gz"

# Ensure log and output directories exist
mkdir -p logs "$OUT_DIR"

# ---- Run HISAT2 Alignment & SAMtools Sorting ----
echo "Starting HISAT2 alignment and BAM sorting with $THREADS threads..."

for R1 in "$READS_DIR"/$R1_PATTERN; do
    # Define Mate 2 and clean sample name
    if [[ "$R1" == *".trimmed.fastq.gz" ]]; then
        R2="${R1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
        SAMPLE=$(basename "$R1" _1.trimmed.fastq.gz)
    else
        R2="${R1/_1.fastq.gz/_2.fastq.gz}"
        SAMPLE=$(basename "$R1" _1.fastq.gz)
    fi

    BAM_OUT="$OUT_DIR/${SAMPLE}.sorted.bam"
    SUMMARY_LOG="$OUT_DIR/${SAMPLE}_hisat2_summary.txt"

    echo "=========================================="
    echo "Mapping and sorting sample: $SAMPLE"
    echo "=========================================="

    # Map paired-end reads with HISAT2 (strandness RF) and pipe directly to samtools sort
    apptainer exec --bind /data,/scratch "$CONTAINER" bash -c "
        hisat2 -p $THREADS \
            --rna-strandness RF \
            -x $INDEX_PREFIX \
            -1 $R1 \
            -2 $R2 \
            --summary-file $SUMMARY_LOG \
            --new-summary | \
        samtools sort -@ $THREADS -T /scratch/${SAMPLE}_temp -o $BAM_OUT -
    "

    # Index the sorted BAM file for downstream visualization (e.g., IGV) and QC
    echo "Indexing $BAM_OUT..."
    apptainer exec --bind /data "$CONTAINER" samtools index -@ "$THREADS" "$BAM_OUT"

    echo "Finished $SAMPLE"
done

echo "All samples successfully mapped, sorted, and indexed in $OUT_DIR!"
