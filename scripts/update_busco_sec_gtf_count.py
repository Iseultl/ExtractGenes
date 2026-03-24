#!/usr/bin/env python3
"""Recompute GTF-derived selenocysteine columns in BUSCO.tsv.

This script updates:
  - Selenocysteine_GTF_Count
  - Selenocysteine_GTF_Count_Genes

from feature rows with type "Selenocysteine" in:

    <dataset_root>/<annotation_id>/All_predictions.gtf
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TARGET_COLUMN_INSTANCES = "Selenocysteine_GTF_Count"
TARGET_COLUMN_GENES = "Selenocysteine_GTF_Count_Genes"


def _get_gtf_attr(attr_string: str, key: str) -> str | None:
    """Return a GTF attribute value for key from either key "value" or key=value."""
    m = re.search(rf'(?:^|;\s*){re.escape(key)}\s+"([^"]+)"', attr_string)
    if m:
        return m.group(1).strip()
    m = re.search(rf'(?:^|;\s*){re.escape(key)}=([^;]+)', attr_string)
    if m:
        return m.group(1).strip().strip('"')
    return None


def count_selenocysteine_sites_and_genes(gtf_path: Path) -> tuple[int, int]:
    """Return (selenocysteine_instances, unique_genes_with_selenocysteine)."""
    instances = 0
    gene_ids: set[str] = set()
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "Selenocysteine":
                instances += 1
                gene_id = _get_gtf_attr(parts[8], "gene_id")
                if gene_id:
                    gene_id = re.sub(r"^Sec\d+:", "", gene_id)
                    gene_ids.add(gene_id)
    return instances, len(gene_ids)


def update_busco_tsv(
    busco_tsv: Path,
    dataset_root: Path,
    backup: bool,
    set_missing_to_zero: bool,
) -> tuple[int, int, int]:
    """Update BUSCO.tsv in place.

    Returns:
        (updated_rows, missing_gtf_rows, unchanged_rows)
    """
    with busco_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if not fieldnames:
        raise ValueError(f"No header found in {busco_tsv}")

    if TARGET_COLUMN_INSTANCES not in fieldnames:
        fieldnames.append(TARGET_COLUMN_INSTANCES)
    if TARGET_COLUMN_GENES not in fieldnames:
        fieldnames.append(TARGET_COLUMN_GENES)

    updated = 0
    missing = 0
    unchanged = 0

    for row in rows:
        annotation_id = (row.get("annotation_id") or "").strip()
        if not annotation_id:
            unchanged += 1
            continue

        gtf_path = dataset_root / annotation_id / "All_predictions.gtf"
        if gtf_path.exists():
            instances, genes = count_selenocysteine_sites_and_genes(gtf_path)
            row[TARGET_COLUMN_INSTANCES] = str(instances)
            row[TARGET_COLUMN_GENES] = str(genes)
            updated += 1
        else:
            missing += 1
            if set_missing_to_zero:
                row[TARGET_COLUMN_INSTANCES] = "0"
                row[TARGET_COLUMN_GENES] = "0"

    if backup:
        backup_path = busco_tsv.with_suffix(busco_tsv.suffix + ".bak")
        backup_path.write_text(busco_tsv.read_text())

    tmp_path = busco_tsv.with_suffix(busco_tsv.suffix + ".tmp")
    with tmp_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    tmp_path.replace(busco_tsv)
    return updated, missing, unchanged


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Update Selenocysteine_GTF_Count and "
            "Selenocysteine_GTF_Count_Genes in BUSCO.tsv from All_predictions.gtf "
            "per annotation_id."
        )
    )
    parser.add_argument(
        "--busco-tsv",
        type=Path,
        default=Path("BUSCO/eukaryota_odb12/BUSCO.tsv"),
        help="Path to BUSCO.tsv to update",
    )
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path("/Users/iseult/Desktop/ML_Dataset"),
        help="Root folder that contains annotation_id folders",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip writing BUSCO.tsv.bak before updating",
    )
    parser.add_argument(
        "--set-missing-to-zero",
        action="store_true",
        help=(
            "Set Selenocysteine_GTF_Count and Selenocysteine_GTF_Count_Genes "
            "to 0 if All_predictions.gtf is missing"
        ),
    )
    args = parser.parse_args()

    updated, missing, unchanged = update_busco_tsv(
        busco_tsv=args.busco_tsv,
        dataset_root=args.dataset_root,
        backup=not args.no_backup,
        set_missing_to_zero=args.set_missing_to_zero,
    )

    print(f"Updated rows: {updated}")
    print(f"Missing GTF rows: {missing}")
    print(f"Unchanged rows: {unchanged}")
    print(f"Updated file: {args.busco_tsv}")


if __name__ == "__main__":
    main()
