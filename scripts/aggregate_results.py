#!/usr/bin/env python3
"""
Aggregate per-job result fragments into the shared TSV files and copy
per-annotation output files (BUSCO.gff, BUSCO.fasta, Selenoprofiles.*) into
the repository output tree.

Usage:
    python aggregate_results.py <artifacts_dir> <busco_tsv> <retry_tsv> [<outputs_dir>]

Each run-analysis job uploads up to two fragment files plus an optional
annotation output directory:
    result_<annotation_id>.tsv        -- one BUSCO row (header + data) on success
    log_<annotation_id>.tsv           -- one retry row (header + data) on failure
    <annotation_id>/                  -- per-annotation output files (BUSCO.gff,
                                         BUSCO.fasta, Selenoprofiles.*)

This script:
  - Scans <artifacts_dir> recursively for fragment TSVs and appends new rows,
    skipping any already present.
  - Copies per-annotation output directories to <outputs_dir>/<annotation_id>/
    (default: outputs/ relative to CWD).
  - BUSCO.tsv    dedup key: annotation_id       (one success row per annotation)
  - .retry.log   dedup key: (annotation_id, run_at)  (full history of failures)
"""
import sys
import csv
import logging
import shutil
from pathlib import Path

from utils import HEADER, RETRY_HEADER

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_existing_ids(tsv_path):
    """Return set of annotation_ids already in a TSV file."""
    p = Path(tsv_path)
    if not p.exists():
        return set()
    ids = set()
    with open(p, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('annotation_id'):
                ids.add(row['annotation_id'])
    return ids


def load_existing_retry_entries(tsv_path):
    """Return set of (annotation_id, run_at) tuples already in .retry.log."""
    p = Path(tsv_path)
    if not p.exists():
        return set()
    entries = set()
    with open(p, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('annotation_id') and row.get('run_at'):
                entries.add((row['annotation_id'], row['run_at']))
    return entries


def ensure_header(tsv_path, header):
    """Ensure TSV exists and matches expected header schema."""
    p = Path(tsv_path)
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, 'w', newline='') as f:
            csv.writer(f, delimiter='\t', lineterminator='\n').writerow(header)
        return

    with open(p, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)

    if not rows:
        with open(p, 'w', newline='') as f:
            csv.writer(f, delimiter='\t', lineterminator='\n').writerow(header)
        return

    current_header = rows[0]
    if current_header == header:
        return

    # Re-write file to expected schema, carrying forward overlapping columns.
    data_rows = rows[1:]
    migrated_rows = []
    if data_rows:
        for row in data_rows:
            row_map = {current_header[i]: row[i] if i < len(row) else ''
                       for i in range(len(current_header))}
            migrated_rows.append([row_map.get(col, '') for col in header])

    with open(p, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(migrated_rows)

    logger.info(f"Updated header/schema for {p}")


def append_rows(tsv_path, rows):
    """Append a list of dicts to a TSV file (no header written)."""
    if not rows:
        return
    with open(tsv_path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writerows(rows)


def read_fragment(fragment_path, expected_header):
    """Read a fragment TSV. Returns list of row dicts matching expected_header."""
    rows = []
    with open(fragment_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if set(expected_header).issubset(row.keys()):
                rows.append({k: row[k] for k in expected_header})
    return rows


# Filenames produced per annotation that identify an annotation output directory.
_ANNOTATION_OUTPUT_FILES = frozenset({
    "BUSCO.fasta", "BUSCO.gtf",
    "All_predictions.gtf", "Selenoprofiles.fasta",
    "Selenoprofiles_annotation_result.csv",
})


def copy_annotation_outputs(artifacts_dir, outputs_dir):
    """
    Locate per-annotation output directories anywhere inside artifacts_dir by
    searching for directories that contain known output file names, then copy
    their contents to <outputs_dir>/<annotation_id>/.

    This is intentionally content-based rather than depth-based so it works
    regardless of how actions/download-artifact@v4 nests the artifact directory
    structure (e.g. artifacts/analysis-batch-0/batch_output/<annotation_id>/).
    """
    outputs_dir = Path(outputs_dir)
    copied_count = 0

    # Collect annotation output dirs by finding any of the known output files.
    seen_dirs: set[Path] = set()
    for fname in sorted(_ANNOTATION_OUTPUT_FILES):
        for found_file in sorted(artifacts_dir.rglob(fname)):
            ann_dir = found_file.parent
            if ann_dir in seen_dirs:
                continue
            seen_dirs.add(ann_dir)

            annotation_id = ann_dir.name
            dest = outputs_dir / annotation_id
            dest.mkdir(parents=True, exist_ok=True)

            for src_file in sorted(ann_dir.iterdir()):
                if not src_file.is_file():
                    continue
                dst_file = dest / src_file.name
                # Skip if identical file already present (same size heuristic)
                if dst_file.exists() and dst_file.stat().st_size == src_file.stat().st_size:
                    logger.debug(f"  ~ skip (same size): {dst_file}")
                    continue
                shutil.copy2(src_file, dst_file)
                logger.info(f"  + {annotation_id}/{src_file.name}")
                copied_count += 1

    logger.info(f"Copied {copied_count} annotation output file(s) to {outputs_dir}")
    return copied_count


def main():
    if len(sys.argv) not in (4, 5):
        print("Usage: python aggregate_results.py "
              "<artifacts_dir> <busco_tsv> <retry_tsv> [<outputs_dir>]")
        sys.exit(1)

    artifacts_dir = Path(sys.argv[1])
    busco_tsv     = sys.argv[2]
    retry_tsv     = sys.argv[3]
    outputs_dir   = sys.argv[4] if len(sys.argv) == 5 else "outputs"
    print(f"Artifacts dir: {artifacts_dir}")
    print(f"BUSCO.tsv    : {busco_tsv}")
    print(f".retry.log   : {retry_tsv}")
    print(f"Outputs dir   : {outputs_dir}") 
    
    if not artifacts_dir.is_dir():
        logger.error(f"Artifacts directory not found: {artifacts_dir}")
        sys.exit(1)

    existing_busco_ids      = load_existing_ids(busco_tsv)
    existing_retry_entries  = load_existing_retry_entries(retry_tsv)
    logger.info(f"Existing BUSCO rows   : {len(existing_busco_ids)}")
    logger.info(f"Existing retry rows   : {len(existing_retry_entries)}")

    ensure_header(busco_tsv, HEADER)
    ensure_header(retry_tsv, RETRY_HEADER)

    busco_new = []
    retry_new = []

    result_fragments = sorted(artifacts_dir.rglob("result_*.tsv"))
    log_fragments    = sorted(artifacts_dir.rglob("log_*.tsv"))

    for frag in result_fragments:
        rows = read_fragment(frag, HEADER)
        print(rows)
        for row in rows:
            if row['annotation_id'] not in existing_busco_ids:
                busco_new.append(row)
                existing_busco_ids.add(row['annotation_id'])
                logger.info(f"  + BUSCO: {row['annotation_id']}")
            else:
                logger.info(f"  ~ skip (already exists): {row['annotation_id']}")

    for frag in log_fragments:
        rows = read_fragment(frag, RETRY_HEADER)
        for row in rows:
            key = (row['annotation_id'], row['run_at'])
            if key not in existing_retry_entries:
                retry_new.append(row)
                existing_retry_entries.add(key)
                logger.info(f"  + retry: {row['annotation_id']} @ {row['run_at']}")

    append_rows(busco_tsv, busco_new)
    append_rows(retry_tsv, retry_new)
    logger.info(f"Appended {len(busco_new)} BUSCO rows and {len(retry_new)} retry rows.")

    # Copy per-annotation output files (BUSCO.gff, BUSCO.fasta, Selenoprofiles.*)
    copy_annotation_outputs(artifacts_dir, outputs_dir)


if __name__ == "__main__":
    main()
