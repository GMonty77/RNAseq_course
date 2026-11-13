#!/usr/bin/env bash
#SBATCH --job-name=featurecounts
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Temporarily allow unbound variables if /etc/bashrc is triggered, then enforce strict mode
set +u
set -euo pipefail

# ---- Settings & Variables ----
WORKDIR=/data/users/$USER/rnaseq_course
RES_DIR=$WORKDIR/resources
BAM_DIR=$WORKDIR/results/hisat2
OUT_DIR=$WORKDIR/results/featurecounts
OUT_FILE=$OUT_DIR/gene_counts.txt
CLEAN_OUT=$OUT_DIR/gene_counts_matrix.txt
THREADS=$SLURM_CPUS_PER_TASK

# Locate the Ensembl GTF annotation file in resources/
ENSEMBL_REL=113
GTF=$RES_DIR/Mus_musculus.GRCm39.${ENSEMBL_REL}.gtf
if [ ! -f "$GTF" ]; then
    GTF=$(ls "$RES_DIR"/*.gtf 2>/dev/null | head -n 1)
fi

# Locate the subread (featureCounts) Apptainer image in /containers/apptainer/
CONTAINER=/containers/apptainer/subread_2.0.1--hed695b0_0.sif
if [ ! -f "$CONTAINER" ]; then
    CONTAINER=$(ls /containers/apptainer/*subread*.sif 2>/dev/null | head -n 1)
fi

# Ensure log and output directories exist
mkdir -p logs "$OUT_DIR"

echo "Using container: $CONTAINER"
echo "Using annotation: $GTF"
echo "Counting mapped reads/fragments per gene across all BAM files..."

# ---- 1. Run featureCounts ----
apptainer exec --bind /data "$CONTAINER" featureCounts \
    -T "$THREADS" \
    -p \
    -s 2 \
    -Q 10 \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$OUT_FILE" \
    "$BAM_DIR"/*.sorted.bam

# ---- 2. Format a clean count matrix for DESeq2 in R ----
# Removes the first '#' comment line and extracts column 1 (Geneid) + columns 7 onwards (sample counts)
tail -n +2 "$OUT_FILE" | cut -f 1,7- > "$CLEAN_OUT"

echo "Feature quantification complete!"
echo "Full output table: $OUT_FILE"
echo "Summary stats:     ${OUT_FILE}.summary"
echo "Clean R matrix:    $CLEAN_OUT"
