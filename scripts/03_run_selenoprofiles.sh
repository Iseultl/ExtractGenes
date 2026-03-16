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
#   all_predictions.gtf              - Raw GTF produced by selenoprofiles
#   Selenoprofiles.fasta             - Transcript sequences extracted by gffread
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
    # Run via Docker container.
    # Mount genome and annotation as read-only and OUTDIR as the output volume.
    # This avoids copying large files to a staging directory and eliminates the
    # root-ownership problem that arises when the container writes to a host-owned
    # temp directory (which then cannot be cleaned up by the non-root runner).
    docker run --rm \
        -v "$GENOME":/input/genome.fa:ro \
        -v "$ANNOTATION":/input/annotation.gff3:ro \
        -v "$OUTDIR":/output \
        maxtico/selenoprofiles_container:latest \
        selenoprofiles \
            -o /output/selenoprofiles_run \
            -t /input/genome.fa \
            -s eukarya \
            -p eukarya \
            -output_gtf_file /output/all_predictions.gtf

    # The container runs as root so output files are root-owned; fix permissions
    # so the host runner can read and copy them.
    docker run --rm \
        -v "$OUTDIR":/output \
        busybox \
        chmod -R a+rw /output
fi

# ---------------------------------------------------------------------------
# Verify the GTF was produced
# ---------------------------------------------------------------------------
if [ ! -f "$GTF_OUT" ]; then
    echo "Error: selenoprofiles did not produce $GTF_OUT"
    exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Selenoprofiles GTF: $GTF_OUT"

# ---------------------------------------------------------------------------
# Expand the region of the selenoprotein predictions by 300bp upstream and downstream
# This is needed to capture the full transcript sequence for gffread extraction
# ---------------------------------------------------------------------------
if [ -x "$(command -v python)" ]; then
    python "$(dirname "$0")/expand_gtf_regions.py" \
        --input_gtf "$GTF_OUT" \
        --output_gtf "$GTF_OUT" \
        --genome_fasta "$GENOME" \
        --expand_upstream 300 \
        --expand_downstream 300
else
    echo "Error: Python is not available in PATH"
    exit 1
fi

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
        -v "$OUTDIR":/sp_out \
        -v "$(dirname "$ANNOTATION")":/annot_dir:ro \
        -v "$(dirname "$GENOME")":/genome_dir:ro \
        maxtico/selenoprofiles_container:latest \
        selenoprofiles assess \
            -s /sp_out/all_predictions.gtf \
            -e /annot_dir/$(basename "$ANNOTATION") \
            -f /genome_dir/$(basename "$GENOME") \
            -o /sp_out/Selenoprofiles_annotation_result.csv

    docker run --rm \
        -v "$OUTDIR":/sp_out \
        busybox \
        chmod -R a+rw /sp_out
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Assessment written to $OUTDIR/Selenoprofiles_annotation_result.csv"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done. Output directory: $OUTDIR"