#!/bin/bash
# Run selenoprofiles on a genome, extract transcript sequences, and compare
# against the reference annotation.
#
# Usage:
#   ./03_run_selenoprofiles.sh <genome_fasta> <annotation_gff3> <output_dir>
#
# Arguments:
#   genome_fasta    - Path to the genome FASTA file (decompressed)
#   annotation_gff3 - Path to the reference annotation GFF3 file (decompressed,
#                     aliased so chromosome IDs match the genome)
#   output_dir      - Directory where results will be written; created if absent
#
# Outputs written to <output_dir>:
#   Selenoprofiles.gtf              - GTF produced by selenoprofiles
#   Selenoprofiles.fasta            - Transcript sequences extracted by gffread
#   Selenoprofiles_annotation_result.csv - Comparison vs reference annotation

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------
if [ "$#" -lt 3 ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 <genome_fasta> <annotation_gff3> <output_dir>"
    exit 1
fi

GENOME=$(realpath "$1")
ANNOTATION=$(realpath "$2")
OUTDIR="$3"

if [ ! -f "$GENOME" ]; then
    echo "Error: Genome file '$GENOME' not found"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo "Error: Annotation file '$ANNOTATION' not found"
    exit 1
fi

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")

GTF_OUT="$OUTDIR/all_predictions.gtf"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running selenoprofiles on genome: $GENOME"

# ---------------------------------------------------------------------------
# Run selenoprofiles
# The container exposes the 'selenoprofiles' command with profiles pre-installed.
# -o  output folder
# -t  target genome FASTA
# -s  species name (tells selenoprofiles which pre-built blast db to use)
# -p  comma-separated profile list to search (eukarya covers common selenoproteins)
# ---------------------------------------------------------------------------
if command -v selenoprofiles &> /dev/null; then
    # Direct execution (selenoprofiles installed in PATH via conda or pip)
    selenoprofiles \
        -o "$OUTDIR/selenoprofiles_run" \
        -t "$GENOME" \
        -s eukarya \
        -p eukarya \
        -output_gtf_file "$GTF_OUT"
else
    # Run via Docker container
    # All inputs are mounted from a shared staging directory so the container
    # can reach both the genome and write output atomically.
    STAGING=$(mktemp -d)
    trap 'rm -rf "$STAGING"' EXIT
    cp "$GENOME"     "$STAGING/genome.fa"
    cp "$ANNOTATION" "$STAGING/annotation.gff3"
    mkdir -p "$STAGING/sp_out"

    docker run --rm \
        -v "$STAGING:/data" \
        maxtico/selenoprofiles_container:latest \
        selenoprofiles \
            -o /data/sp_out \
            -t /data/genome.fa \
            -s eukarya \
            -p eukarya \
            -output_gtf_file /data/sp_out/all_predictions.gtf

    cp "$STAGING/sp_out/all_predictions.gtf" "$GTF_OUT" 2>/dev/null || true
    cp -r "$STAGING/sp_out/." "$OUTDIR/selenoprofiles_run/" 2>/dev/null || true
fi

# ---------------------------------------------------------------------------
# Copy / rename the GTF to the canonical output name
# ---------------------------------------------------------------------------
if [ ! -f "$GTF_OUT" ]; then
    echo "Error: selenoprofiles did not produce $GTF_OUT"
    exit 1
fi
cp "$GTF_OUT" "$OUTDIR/Selenoprofiles.gtf"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Selenoprofiles GTF written to $OUTDIR/Selenoprofiles.gtf"

# ---------------------------------------------------------------------------
# Extract transcript sequences with gffread
# ---------------------------------------------------------------------------
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Extracting selenoprotein transcript sequences with gffread"
gffread "$GTF_OUT" -g "$GENOME" -w "$OUTDIR/Selenoprofiles.fasta"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Transcript FASTA written to $OUTDIR/Selenoprofiles.fasta"

# ---------------------------------------------------------------------------
# Compare selenoprofiles predictions vs the reference annotation
# selenoprofiles assess  -s <predictions.gtf>  -e <reference.gff3>
#                        -f <genome.fa>         -o <output.csv>
# ---------------------------------------------------------------------------
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running selenoprofiles assess vs reference annotation"
if command -v selenoprofiles &> /dev/null; then
    selenoprofiles assess \
        -s "$GTF_OUT" \
        -e "$ANNOTATION" \
        -f "$GENOME" \
        -o "$OUTDIR/Selenoprofiles_annotation_result.csv"
else
    docker run --rm \
        -v "$OUTDIR:/sp_out" \
        -v "$(dirname "$ANNOTATION"):/annot_dir:ro" \
        -v "$(dirname "$GENOME"):/genome_dir:ro" \
        maxtico/selenoprofiles_container:latest \
        selenoprofiles assess \
            -s /sp_out/Selenoprofiles.gtf \
            -e /annot_dir/$(basename "$ANNOTATION") \
            -f /genome_dir/$(basename "$GENOME") \
            -o /sp_out/Selenoprofiles_annotation_result.csv
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Assessment written to $OUTDIR/Selenoprofiles_annotation_result.csv"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done. Output directory: $OUTDIR"