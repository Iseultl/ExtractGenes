#!/usr/bin/env python3
"""Expand GTF intervals upstream/downstream while clamping to chromosome bounds."""

import argparse
from pathlib import Path


def load_fasta_lengths(fasta_path):
    """Return {sequence_id: sequence_length} from a FASTA file."""
    lengths = {}
    current_id = None
    current_len = 0

    with open(fasta_path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)

    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def main():
    parser = argparse.ArgumentParser(
        description="Expand GTF regions upstream and downstream."
    )
    parser.add_argument("--input_gtf", required=True, help="Input GTF file")
    parser.add_argument("--output_gtf", required=True, help="Output GTF file")
    parser.add_argument("--genome_fasta", required=True, help="Genome FASTA file")
    parser.add_argument("--expand_upstream", type=int, default=300)
    parser.add_argument("--expand_downstream", type=int, default=300)
    args = parser.parse_args()

    seq_lengths = load_fasta_lengths(args.genome_fasta)

    in_path = Path(args.input_gtf)
    out_path = Path(args.output_gtf)

    with in_path.open() as src, out_path.open("w") as dst:
        for raw in src:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                dst.write(raw)
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                dst.write(raw)
                continue

            seqname = cols[0]
            strand = cols[6]

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                dst.write(raw)
                continue

            seq_len = seq_lengths.get(seqname)
            if seq_len is None:
                dst.write(raw)
                continue

            if strand == "+":
                new_start = max(1, start - args.expand_upstream)
                new_end = min(seq_len, end + args.expand_downstream)
            else:
                new_start = max(1, start - args.expand_downstream)
                new_end = min(seq_len, end + args.expand_upstream)

            cols[3] = str(new_start)
            cols[4] = str(new_end)
            dst.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    main()