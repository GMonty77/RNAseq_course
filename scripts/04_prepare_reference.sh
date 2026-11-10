#!/usr/bin/env bash
#SBATCH --job-name=prepare_reference
#SBATCH --partition=pibu_el8
#SBATCH --account=class-467713-ws-2026
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Allow unbound variables briefly if /etc/bashrc is triggered, then enforce strict mode
set +u
set -euo pipefail

# ---- Settings & Variables ----
WORKDIR=/data/users/$USER/rnaseq_course
RES_DIR=$WORKDIR/resources
INDEX_DIR=$RES_DIR/hisat2_index
THREADS=$SLURM_CPUS_PER_TASK

# Adjust the container filename if a different version is in /containers/apptainer/
# Run: ls /containers/apptainer/*hisat* to verify the exact .sif name
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

# Ensembl Release (GRCm39) URLs for Mouse Genome (FASTA) and Annotation (GTF)
ENSEMBL_REL=113
FASTA_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_REL}.gtf.gz"

FASTA_GZ="$RES_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
FASTA="$RES_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_GZ="$RES_DIR/Mus_musculus.GRCm39.${ENSEMBL_REL}.gtf.gz"
GTF="$RES_DIR/Mus_musculus.GRCm39.${ENSEMBL_REL}.gtf"
INDEX_PREFIX="$INDEX_DIR/GRCm39"

# Ensure directories exist
mkdir -p logs "$RES_DIR" "$INDEX_DIR"

# ---- 1. Download Reference Genome & Annotation ----
echo "Downloading Ensembl release ${ENSEMBL_REL} FASTA and GTF to $RES_DIR..."

if [ ! -f "$FASTA" ]; then
    wget -O "$FASTA_GZ" "$FASTA_URL"
    echo "Decompressing FASTA..."
    gunzip -f "$FASTA_GZ"
else
    echo "FASTA file already exists, skipping download."
fi

if [ ! -f "$GTF" ]; then
    wget -O "$GTF_GZ" "$GTF_URL"
    echo "Decompressing GTF..."
    gunzip -f "$GTF_GZ"
else
    echo "GTF file already exists, skipping download."
fi

# ---- 2. Build HISAT2 Index ----
echo "Building HISAT2 index with $THREADS threads..."

apptainer exec --bind /data "$CONTAINER" hisat2-build \
    -p "$THREADS" \
    "$FASTA" \
    "$INDEX_PREFIX"

echo "HISAT2 index successfully built at ${INDEX_PREFIX}.*.ht2"
