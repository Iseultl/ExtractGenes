#!/usr/bin/env python3
"""
Run a strided slice of pending BUSCO analyses for one matrix job.

Usage:
    python run_busco_batch.py <annotations_tsv> <busco_tsv> <retry_tsv> \
        <chunk_index> <chunk_count> <output_dir> [max_per_job]

Slicing: pending_sorted[chunk_index::chunk_count]
  e.g.  chunk 0 of 4 on [a,b,c,d,e,f,g,h] → [a,e]
        chunk 1 of 4                        → [b,f]
        chunk 2 of 4                        → [c,g]
        chunk 3 of 4                        → [d,h]

Per-annotation failures are recorded in the retry TSV fragment and
execution continues — the batch script always exits 0.
"""
import sys
import csv
import argparse
import subprocess
import logging
from datetime import datetime
from pathlib import Path

from utils import load_ids, compute_pending_ids, RETRY_HEADER

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_annotations(tsv_path):
    """Return dict of annotation_id → {annotation_url, assembly_url}."""
    annotations = {}
    with open(tsv_path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or not row[0].strip():
                continue
            # Skip header row if present
            if row[0].strip() == 'annotation_id':
                continue
            annotations[row[0].strip()] = {
                'annotation_url': row[1].strip(),
                'assembly_url':   row[2].strip(),
            }
    return annotations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('annotations_tsv', help='Path to annotations.tsv')
    parser.add_argument('busco_tsv',       help='Path to BUSCO.tsv (successful results)')
    parser.add_argument('retry_tsv',       help='Path to .retry.log (failed runs)')
    parser.add_argument('chunk_index',     type=int, help='Index of this chunk (0-based)')
    parser.add_argument('chunk_count',     type=int, help='Total number of chunks')
    parser.add_argument('output_dir',      help='Directory to write result/log fragments')
    parser.add_argument('max_per_job',     type=int, nargs='?', default=None,
                        help='Cap on annotations processed by this chunk')
    parser.add_argument('--giveup-tsv',    default=None,
                        help='Path to .giveup.log (given-up annotations to exclude)')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotations = load_annotations(args.annotations_tsv)
    success_ids = load_ids(args.busco_tsv)
    retry_ids   = load_ids(args.retry_tsv)
    giveup_ids  = load_ids(args.giveup_tsv) if args.giveup_tsv else set()

    pending_ids = compute_pending_ids(set(annotations.keys()), success_ids, retry_ids, giveup_ids)

    my_slice = pending_ids[args.chunk_index::args.chunk_count]
    if args.max_per_job is not None:
        my_slice = my_slice[:args.max_per_job]

    logger.info(f"Chunk {args.chunk_index}/{args.chunk_count}: "
                f"{len(my_slice)} annotations to process"
                + (f" (capped at {args.max_per_job})" if args.max_per_job else ""))

    script = Path(__file__).parent / 'run_busco_analysis.py'
    succeeded  = 0
    fail_count = 0

    for i, annotation_id in enumerate(my_slice, 1):
        ann = annotations[annotation_id]
        result_tsv   = str(output_dir / f"result_{annotation_id}.tsv")
        log_fragment = str(output_dir / f"log_{annotation_id}.tsv")

        logger.info(f"[{i}/{len(my_slice)}] Processing {annotation_id}")

        ann_output_dir = str(output_dir / annotation_id)

        try:
            ret = subprocess.run(
                [sys.executable, str(script),
                 ann['annotation_url'],
                 ann['assembly_url'],
                 annotation_id,
                 result_tsv,
                 log_fragment,
                 ann_output_dir],
                check=False   # do not raise on non-zero exit
            )
            if ret.returncode == 0:
                succeeded += 1
                logger.info(f"  ✓ {annotation_id}")
            else:
                fail_count += 1
                logger.warning(f"  ✗ {annotation_id} (exit {ret.returncode})")
        except Exception as e:
            fail_count += 1
            logger.error(f"  ✗ {annotation_id} — unexpected error: {e}")
            # Write a log fragment so the aggregator records this failure
            # and the annotation is not silently rescheduled
            try:
                with open(log_fragment, 'w', newline='') as lf:
                    w = csv.writer(lf, delimiter='\t', lineterminator='\n')
                    w.writerow(RETRY_HEADER)
                    w.writerow([annotation_id,
                                datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                'unexpected_error'])
            except Exception as write_err:
                logger.error(f"  Could not write log fragment: {write_err}")

    logger.info(f"Chunk {args.chunk_index} complete: "
                f"{succeeded} succeeded, {fail_count} failed out of {len(my_slice)}")
    # Always exit 0 — individual failures are recorded in log fragments
    sys.exit(0)


if __name__ == '__main__':
    main()
