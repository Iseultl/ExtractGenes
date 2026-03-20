#!/usr/bin/env python3
"""Backfill Stop_codon counts into BUSCO.tsv from ML_Dataset outputs."""

import argparse
import csv
from pathlib import Path


def count_stop_codons(result_csv: Path) -> int:
    if not result_csv.exists():
        return 0
    stop_count = 0
    with result_csv.open(newline="") as sf:
        reader = csv.DictReader(sf, delimiter="\t")
        for row in reader:
            annotation = (row.get("Type_annotation") or "").strip().replace(" ", "_")
            if annotation == "Stop_codon":
                stop_count += 1
    return stop_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("busco_tsv", help="Path to BUSCO.tsv")
    parser.add_argument("ml_dataset_root", help="Path to ML_Dataset root")
    args = parser.parse_args()

    busco_path = Path(args.busco_tsv)
    ml_root = Path(args.ml_dataset_root)

    header = [
        "annotation_id",
        "assembly_accession",
        "species",
        "lineage",
        "busco_count",
        "complete",
        "single",
        "duplicated",
        "fragmented",
        "missing",
        "Well_annotated",
        "Upstream",
        "Downstream",
        "Skipped",
        "Out_of_frame",
        "Spliced",
        "Stop_codon",
        "Selenocysteine_Gene_Count",
        "Selenocysteine_GTF_Count",
    ]

    with busco_path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    for row in rows:
        annotation_id = (row.get("annotation_id") or "").strip()
        result_csv = ml_root / annotation_id / "Selenoprofiles_annotation_result.csv"
        row["Stop_codon"] = str(count_stop_codons(result_csv))

    with busco_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "NA") for k in header})

    print(f"Updated {len(rows)} rows in {busco_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
