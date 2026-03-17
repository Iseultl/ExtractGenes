#!/usr/bin/env python3
"""Shared utilities for BUSCO-tracker scripts."""
import csv
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

HEADER = ['annotation_id', 'assembly_accession', 'species', 'lineage', 'busco_count', 'complete',
                'single', 'duplicated', 'fragmented', 'missing', 'Well_annotated', 'Upstream', 'Downstream', 'Skipped', 'Out_of_frame', 'Spliced']
RETRY_HEADER = ['annotation_id', 'run_at', 'step']
GIVEUP_HEADER = RETRY_HEADER


def compute_pending_ids(all_ids, success_ids, error_ids, giveup_ids=None):
    """Return pending IDs in priority order: never-run first, then failed retries."""
    if giveup_ids is None:
        giveup_ids = set()
    never_run = all_ids - success_ids - error_ids - giveup_ids
    failed    = (error_ids - success_ids) & all_ids
    return sorted(never_run) + sorted(failed)


def load_ids(tsv_path, column='annotation_id'):
    """Return set of values from a TSV column. Returns empty set if file missing."""
    p = Path(tsv_path)
    if not p.exists():
        logger.info(f"{tsv_path} not found — treating as empty")
        return set()
    ids = set()
    with open(p, newline='') as f:
        first_line = f.readline()
        if not first_line:
            return ids
        f.seek(0)
        has_header = first_line.split('\t')[0].strip() == column
        if has_header:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get(column):
                    ids.add(row[column])
        else:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row and row[0].strip():
                    ids.add(row[0].strip())
    return ids

