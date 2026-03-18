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
#   all_predictions_expanded.gtf     - Expanded-coordinate GTF used for FASTA extraction
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

GENOME="$1"
ANNOTATION="$2"
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

GTF_OUT="$OUTDIR/all_predictions.gtf"
GTF_EDITED="$OUTDIR/all_predictions_expanded.gtf"

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
        -o "$OUTDIR" \
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
            -o /output \
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
        --output_gtf "$GTF_EDITED" \
        --genome_fasta "$GENOME" \
        --expand_upstream 300 \
        --expand_downstream 300
else
    echo "Error: Python is not available in PATH"
    exit 1
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Expanded GTF regions written to $GTF_EDITED"
# ---------------------------------------------------------------------------
# Extract transcript sequences with gffread
# ---------------------------------------------------------------------------
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Extracting selenoprotein transcript sequences with gffread"
gffread "$GTF_EDITED" -g "$GENOME" -w "$OUTDIR/Selenoprofiles.fasta"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Transcript FASTA written to $OUTDIR/Selenoprofiles.fasta"

# ---------------------------------------------------------------------------
# Compare selenoprofiles predictions vs the reference annotation.
# Skip assess (and write an empty CSV) when the GTF has no predictions —
# this is valid for genomes with no selenoproteins and should not be treated
# as a pipeline failure.
# selenoprofiles assess  -s <predictions.gtf>  -e <reference.gff3>
#                        -f <genome.fa>         -o <output.csv>
# ---------------------------------------------------------------------------
CSV_OUT="$OUTDIR/Selenoprofiles_annotation_result.csv"
ASSESS_STDERR="$OUTDIR/Selenoprofiles_assess.stderr.log"
ASSESS_GTF="$OUTDIR/all_predictions.assess.gtf"
ASSESS_ANNOT="$OUTDIR/annotation.assess.gff3"
GENOME_CONTIGS="$OUTDIR/genome.contigs.txt"

# Only run assess when the GTF contains actual selenocysteine predictions.
# Selenoprofiles may still write a GTF for other residue types, but the
# assess step is only meaningful when "selenocysteine" is present.
if ! grep -qi 'selenocysteine' "$GTF_OUT"; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] No selenocysteine predictions in GTF; skipping assess step"
    printf 'transcript_id\ttranscript_id_ens\tType_annotation\n' > "$CSV_OUT"
else
    # Build a set of contig IDs present in the FASTA and filter both
    # prediction/reference files to those contigs. This prevents pyfaidx
    # crashes (e.g., "MT not in assembly.fna") from assess.
    awk '/^>/{h=$1; sub(/^>/, "", h); print h}' "$GENOME" | sort -u > "$GENOME_CONTIGS"

    awk -v contigs="$GENOME_CONTIGS" '
        BEGIN {
            while ((getline line < contigs) > 0) {
                keep[line] = 1
            }
            close(contigs)
        }
        /^#/ { print; next }
        ($1 in keep) { print }
    ' "$GTF_OUT" > "$ASSESS_GTF"

    awk -v contigs="$GENOME_CONTIGS" '
        BEGIN {
            while ((getline line < contigs) > 0) {
                keep[line] = 1
            }
            close(contigs)
        }
        /^#/ { print; next }
        ($1 in keep) { print }
    ' "$ANNOTATION" > "$ASSESS_ANNOT"

    if ! awk 'BEGIN{ok=0} !/^#/ && NF>0 {ok=1} END{exit(ok?0:1)}' "$ASSESS_GTF"; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] No assessable predictions after contig filtering; writing empty assessment"
        printf 'transcript_id\ttranscript_id_ens\tType_annotation\n' > "$CSV_OUT"
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] Assessment written to $CSV_OUT"
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done. Output directory: $OUTDIR"
        exit 0
    fi

    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Running selenoprofiles assess vs reference annotation"
    assess_ok=1
    if command -v selenoprofiles &> /dev/null; then
        if ! selenoprofiles assess \
            -s "$ASSESS_GTF" \
            -e "$ASSESS_ANNOT" \
            -f "$GENOME" \
            -o "$CSV_OUT" \
            2> "$ASSESS_STDERR"; then
            assess_ok=0
        fi
    else
        if ! docker run --rm \
            -v "$OUTDIR":/sp_out \
            -v "$(dirname "$GENOME")":/genome_dir:ro \
            maxtico/selenoprofiles_container:latest \
            selenoprofiles assess \
                -s /sp_out/$(basename "$ASSESS_GTF") \
                -e /sp_out/$(basename "$ASSESS_ANNOT") \
                -f /genome_dir/$(basename "$GENOME") \
                -o /sp_out/Selenoprofiles_annotation_result.csv \
                2> "$ASSESS_STDERR"; then
            assess_ok=0
        fi

        docker run --rm \
            -v "$OUTDIR":/sp_out \
            busybox \
            chmod -R a+rw /sp_out
    fi

    if [ "$assess_ok" -ne 1 ]; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] Warning: selenoprofiles assess failed; writing empty assessment and continuing"
        if grep -qi "unexpected keyword argument 'int64'" "$ASSESS_STDERR"; then
            echo "[$(date +'%Y-%m-%d %H:%M:%S')] Detected PyRanges API incompatibility during assess"
        fi
        if grep -qi "not in .*assembly" "$ASSESS_STDERR"; then
            echo "[$(date +'%Y-%m-%d %H:%M:%S')] Detected chromosome mismatch during assess"
        fi
        printf 'transcript_id\ttranscript_id_ens\tType_annotation\n' > "$CSV_OUT"
    fi
fi
echo "[$(date +'%Y-%m-%d %H:%M:%S')] Assessment written to $CSV_OUT"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] Done. Output directory: $OUTDIR"