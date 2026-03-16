#!/usr/bin/env python3
"""
Per-annotation analysis pipeline: BUSCO + Selenoprofiles

For each annotation this script:
  1. Downloads GFF + assembly FASTA from the AnnoTrEive API URLs.
  2. Aliases sequence IDs with annocli so GFF chromosome names match the FASTA.
  3. Extracts the longest protein isoform per gene (01_extract_proteins.sh).
  4. Runs BUSCO in protein mode (02_run_BUSCO.sh).
  5. Extracts transcript sequences + GFF records for BUSCO-matched genes.
  6. Runs selenoprofiles on the genome to identify selenoproteins
     (03_run_selenoprofiles.sh).
  7. Writes all output files to <output_dir>/<annotation_id>/.
  8. Writes a result fragment TSV on success or a log fragment TSV on failure.

Usage:
    python run_busco_analysis.py \\
        <annotation_url> <assembly_url> <annotation_id> \\
        <result_tsv> <log_tsv> <output_dir>
"""
import csv
import gzip
import logging
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path

from utils import HEADER, RETRY_HEADER

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers: downloading
# ---------------------------------------------------------------------------

def download_file(url, dest_path):
    """Download *url* to *dest_path*. Returns (ok: bool, error_message: str)."""
    logger.info(f"Downloading {url} -> {dest_path}")
    try:
        urllib.request.urlretrieve(url, dest_path)
        logger.info(f"Download complete: {dest_path}")
        return True, "NA"
    except urllib.error.HTTPError as e:
        msg = f"HTTP {e.code} downloading {url}: {e.reason}"
        logger.error(msg)
        return False, msg
    except urllib.error.URLError as e:
        msg = f"URL error downloading {url}: {e.reason}"
        logger.error(msg)
        return False, msg
    except Exception as e:
        msg = f"Unexpected error downloading {url}: {e}"
        logger.error(msg)
        return False, msg


# ---------------------------------------------------------------------------
# Helpers: running shell scripts
# ---------------------------------------------------------------------------

def run_shell_script(script_path, args, step_name):
    """Run a shell script; return (ok: bool, stdout: str, stderr: str)."""
    # Run scripts via bash so they do not depend on executable file mode.
    cmd = ["bash", str(script_path)] + args
    logger.info(f"Running {step_name}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"{step_name} completed successfully")
        return True, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} failed (exit {e.returncode})")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        return False, e.stdout, e.stderr
    except FileNotFoundError as e:
        logger.error(f"Script not found: {script_path}")
        return False, "", str(e)


# ---------------------------------------------------------------------------
# Helpers: BUSCO result parsing
# ---------------------------------------------------------------------------

def parse_busco_results(busco_output_dir):
    """Parse BUSCO short-summary; return dict with lineage/counts."""
    logger.info(f"Parsing BUSCO results from {busco_output_dir}")
    summary_files = list(Path(busco_output_dir).glob("short_summary.*.txt"))
    if not summary_files:
        raise ValueError(f"BUSCO summary file not found in {busco_output_dir}")

    content = summary_files[0].read_text()
    results: dict = {
        "lineage": "", "busco_count": None,
        "complete": None, "single": None,
        "duplicated": None, "fragmented": None, "missing": None,
    }
    for pattern, key, cast in [
        (r"lineage dataset is: (\S+)", "lineage", str),
        (r"C:(\d+(?:\.\d+)?)%", "complete", float),
        (r"S:(\d+(?:\.\d+)?)%", "single", float),
        (r"D:(\d+(?:\.\d+)?)%", "duplicated", float),
        (r"F:(\d+(?:\.\d+)?)%", "fragmented", float),
        (r"M:(\d+(?:\.\d+)?)%", "missing", float),
        (r"(\d+)\s+total BUSCO", "busco_count", int),
    ]:
        m = re.search(pattern, content, re.IGNORECASE)
        if m:
            results[key] = cast(m.group(1))

    logger.info(f"BUSCO results: {results}")
    return results


def get_busco_matched_ids(busco_output_dir):
    """
    Return the set of protein/transcript IDs that BUSCO matched
    (status Complete, Duplicated, or Fragmented) from full_table.tsv.
    """
    tables = list(Path(busco_output_dir).glob("run_*/full_table.tsv"))
    if not tables:
        logger.warning("full_table.tsv not found; BUSCO sequence extraction skipped")
        return set()
    matched = set()
    with open(tables[0]) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            status, seq_id = cols[1], cols[2]
            if status in ("Complete", "Duplicated", "Fragmented"):
                matched.add(seq_id)
    logger.info(f"BUSCO matched IDs: {len(matched)}")
    return matched

# ----------------------------------------------------------------------------
# Selenoprofiles result parsing
# ----------------------------------------------------------------------------
def parse_selenoprofiles_results(seleno_result_csv):
    """
    Parse Selenoprofiles_annotation_result.csv and return a dict with counts
    for each annotation category in the third column.
    Returns a dict with keys: Downstream, Well_annotated, Upstream, Out_of_frame, Skipped.
    If the file does not exist, returns all counts as 0.
    """
    categories = ["Downstream", "Well_annotated", "Upstream", "Out_of_frame", "Skipped"]
    counts = {cat: 0 for cat in categories}
    if not Path(seleno_result_csv).exists():
        return counts

    with open(seleno_result_csv, newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) < 3:
                continue
            cat = row[2].strip()
            if cat in counts:
                counts[cat] += 1
    return counts


# ---------------------------------------------------------------------------
# Helpers: GFF3 filtering for BUSCO gene extraction
# ---------------------------------------------------------------------------

def _get_attr(attr_string, key):
    """Return the value of a GFF3 attribute, or None."""
    m = re.search(rf'(?:^|;){re.escape(key)}=([^;]+)', attr_string)
    return m.group(1).strip() if m else None


def _open_gff(path):
    """Open a plain or gzip-compressed GFF file for reading."""
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def extract_busco_sequences(
    busco_output_dir, alias_gff_file, fasta_file, annotation_id, out_dir
):
    """
    Produce BUSCO.gff (filtered reference annotation) and BUSCO.fasta
    (transcript sequences via gffread) for the BUSCO-matched genes.

    Parameters
    ----------
    busco_output_dir : path-like  BUSCO run output directory.
    alias_gff_file   : path-like  Aliased GFF3 (may be .gz).
    fasta_file       : path-like  Decompressed genome FASTA.
    annotation_id    : str        Used only for logging.
    out_dir          : Path       Destination for BUSCO.gff and BUSCO.fasta.
    """
    matched_ids = get_busco_matched_ids(busco_output_dir)
    if not matched_ids:
        logger.warning(f"{annotation_id}: no BUSCO-matched IDs — skipping extraction")
        return

    # -----------------------------------------------------------------------
    # Two-pass GFF filter
    # Pass 1: find parent gene IDs for all matched mRNA IDs
    # -----------------------------------------------------------------------
    gene_ids: set[str] = set()
    mrna_ids: set[str] = set(matched_ids)

    with _open_gff(alias_gff_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            ftype, attrs = parts[2], parts[8]
            if ftype in ("mRNA", "transcript"):
                feat_id = _get_attr(attrs, "ID")
                if feat_id and feat_id in mrna_ids:
                    parent = _get_attr(attrs, "Parent")
                    if parent:
                        gene_ids.add(parent)

    logger.info(
        f"{annotation_id}: filtering GFF — "
        f"{len(mrna_ids)} mRNAs across {len(gene_ids)} genes"
    )

    # -----------------------------------------------------------------------
    # Pass 2: write filtered GFF3 to a temporary file, then convert to GTF
    # and extract transcript sequences.  The temp GFF3 is not kept in the
    # output directory — only BUSCO.gtf and BUSCO.fasta are retained.
    # -----------------------------------------------------------------------
    busco_gff3_tmp = out_dir / "_busco_filtered.gff3"
    with _open_gff(alias_gff_file) as fh, open(busco_gff3_tmp, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            ftype, attrs = parts[2], parts[8]
            feat_id = _get_attr(attrs, "ID")
            parent   = _get_attr(attrs, "Parent")

            keep = False
            if ftype in ("gene", "pseudogene") and feat_id in gene_ids:
                keep = True
            elif ftype in ("mRNA", "transcript") and feat_id in mrna_ids:
                keep = True
            elif ftype not in ("gene", "pseudogene", "mRNA", "transcript") \
                    and parent in mrna_ids:
                keep = True

            if keep:
                out.write(line)

    # Convert filtered GFF3 → GTF format
    busco_gtf = out_dir / "BUSCO.gtf"
    cmd_gtf = [
        "gffread", str(busco_gff3_tmp),
        "-T",
        "-o", str(busco_gtf),
    ]
    logger.info(f"Converting filtered GFF3 to GTF: {' '.join(cmd_gtf)}")
    try:
        subprocess.run(cmd_gtf, check=True, capture_output=True, text=True)
        logger.info(f"{annotation_id}: BUSCO.gtf written to {busco_gtf}")
    except subprocess.CalledProcessError as e:
        logger.error(f"gffread -T failed: {e.stderr}")
        raise

    # -----------------------------------------------------------------------
    # Extract transcript sequences with gffread (w = exon-stitched RNA seqs)
    # Use the filtered GFF3 (not GTF) as gffread reads GFF3 reliably for
    # sequence extraction.
    # -----------------------------------------------------------------------
    busco_fasta = out_dir / "BUSCO.fasta"
    cmd = [
        "gffread", str(busco_gff3_tmp),
        "-g", str(fasta_file),
        "-w", str(busco_fasta),
    ]
    logger.info(f"Extracting BUSCO transcript sequences: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"{annotation_id}: BUSCO.fasta written to {busco_fasta}")
    except subprocess.CalledProcessError as e:
        logger.error(f"gffread failed: {e.stderr}")
        raise
    finally:
        # Remove the intermediate GFF3 — only GTF and FASTA are kept
        busco_gff3_tmp.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Helpers: TSV fragment writers
# ---------------------------------------------------------------------------

def write_result_tsv(result_tsv, annotation_id, busco_results):
    """Write a single-row BUSCO result fragment TSV."""
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerow([
            annotation_id,
            busco_results["lineage"],
            busco_results["busco_count"] if busco_results["busco_count"] is not None else "NA",
            busco_results["complete"]    if busco_results["complete"]    is not None else "NA",
            busco_results["single"]      if busco_results["single"]      is not None else "NA",
            busco_results["duplicated"]  if busco_results["duplicated"]  is not None else "NA",
            busco_results["fragmented"]  if busco_results["fragmented"]  is not None else "NA",
            busco_results["missing"]     if busco_results["missing"]     is not None else "NA",
            busco_results.get("Well_annotated", 0),
            busco_results.get("Upstream", 0),
            busco_results.get("Downstream", 0),
            busco_results.get("Skipped", 0),
            busco_results.get("Out_of_frame", 0),
        ])


def write_log_tsv(log_tsv, annotation_id, step):
    """Write a single-row failure log fragment TSV."""
    run_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    step = " ".join(step.splitlines()).strip()
    with open(log_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(RETRY_HEADER)
        writer.writerow([annotation_id, run_at, step])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 7:
        print(
            "Usage: python run_busco_analysis.py "
            "<annotation_url> <assembly_url> <annotation_id> "
            "<result_tsv> <log_tsv> <output_dir>"
        )
        sys.exit(1)

    annotation_url = sys.argv[1]
    assembly_url   = sys.argv[2]
    annotation_id  = sys.argv[3]
    result_tsv     = sys.argv[4]
    log_tsv        = sys.argv[5]
    output_dir     = Path(sys.argv[6])
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting analysis for {annotation_id}")

    # Locate shell scripts relative to this file
    script_dir        = Path(__file__).parent
    extract_script    = script_dir / "01_extract_proteins.sh"
    busco_script      = script_dir / "02_run_BUSCO.sh"
    seleno_script     = script_dir / "03_run_selenoprofiles.sh"

    for script in (extract_script, busco_script, seleno_script):
        if not script.exists():
            logger.error(f"Required script not found: {script}")
            write_log_tsv(log_tsv, annotation_id, f"script_missing:{script.name}")
            return 1

    work_dir = Path(tempfile.mkdtemp(prefix=f"analysis_{annotation_id}_"))
    logger.info(f"Working directory: {work_dir}")

    try:
        # ------------------------------------------------------------------
        # Step 1: Download annotation GFF and assembly FASTA
        # ------------------------------------------------------------------
        logger.info("STEP 1: Download files")
        gff_file   = work_dir / "annotation.gff.gz"
        fasta_file = work_dir / "assembly.fna.gz"

        ok, err = download_file(annotation_url, gff_file)
        if not ok:
            write_log_tsv(log_tsv, annotation_id, err)
            return 1

        ok, err = download_file(assembly_url, fasta_file)
        if not ok:
            write_log_tsv(log_tsv, annotation_id, err)
            return 1

        # ------------------------------------------------------------------
        # Step 2: Alias sequence IDs with annocli
        # ------------------------------------------------------------------
        logger.info("STEP 2: Alias sequence IDs (annocli alias)")
        alias_gff_file = work_dir / "annotation.aliasMatch.gff3.gz"
        cmd = [
            "annocli", "alias",
            str(gff_file), str(fasta_file),
            "--output", str(alias_gff_file),
        ]
        try:
            alias_result = subprocess.run(
                cmd, capture_output=True, text=True, check=True
            )
            alias_stderr = alias_result.stderr
            logger.info("annocli alias completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"annocli alias failed (exit {e.returncode})")
            logger.error(f"stderr: {e.stderr}")
            write_log_tsv(
                log_tsv, annotation_id,
                e.stderr if e.stderr else "alias_ids_failed"
            )
            return 1
        except FileNotFoundError:
            logger.error("annocli not found in PATH")
            write_log_tsv(log_tsv, annotation_id, "annocli_not_found")
            return 1

        if not alias_gff_file.exists():
            logger.error(f"Expected aliasMatch file not found: {alias_gff_file}")
            write_log_tsv(
                log_tsv, annotation_id,
                alias_stderr if alias_stderr else "alias_output_missing"
            )
            return 1

        # ------------------------------------------------------------------
        # Step 3: Decompress files for downstream steps
        # ------------------------------------------------------------------
        logger.info("STEP 3: Decompress GFF and genome FASTA")
        alias_gff_decompressed = work_dir / "annotation.aliasMatch.gff3"
        fasta_decompressed     = work_dir / "assembly.fna"

        with gzip.open(alias_gff_file, "rb") as src, \
                open(alias_gff_decompressed, "wb") as dst:
            shutil.copyfileobj(src, dst)

        with gzip.open(fasta_file, "rb") as src, \
                open(fasta_decompressed, "wb") as dst:
            shutil.copyfileobj(src, dst)

        # ------------------------------------------------------------------
        # Step 4: Extract proteins (longest isoform per gene)
        # ------------------------------------------------------------------
        logger.info("STEP 4: Extract proteins")
        ok, _, stderr = run_shell_script(
            extract_script,
            [str(alias_gff_file), str(fasta_file)],
            "extract_proteins",
        )
        if not ok:
            write_log_tsv(
                log_tsv, annotation_id,
                stderr if stderr else "extract_proteins_failed"
            )
            return 1

        # 01_extract_proteins.sh writes <gff_basename>_proteins.faa beside the GFF
        protein_file = work_dir / "annotation.aliasMatch_proteins.faa"
        if not protein_file.exists():
            logger.error(f"Expected protein file not found: {protein_file}")
            write_log_tsv(log_tsv, annotation_id, "protein_file_missing")
            return 1

        # ------------------------------------------------------------------
        # Step 5: Run BUSCO (protein mode, offline)
        # ------------------------------------------------------------------
        logger.info("STEP 5: Run BUSCO")
        lineage_candidates = [
            Path("assets/busco_downloads/lineages/eukaryota_odb12"),
            Path("busco_downloads/lineages/eukaryota_odb12"),
            Path("eukaryota_odb12"),
        ]
        lineage_path = next((p for p in lineage_candidates if p.exists()), None)
        if lineage_path is None:
            logger.error(
                "Lineage folder not found. Tried: "
                + ", ".join(str(p) for p in lineage_candidates)
            )
            write_log_tsv(log_tsv, annotation_id, "lineage_missing")
            return 1

        busco_output = str(work_dir / f"busco_{annotation_id}")
        ok, _, stderr = run_shell_script(
            busco_script,
            [str(protein_file), str(lineage_path), busco_output],
            "run_busco",
        )
        if not ok:
            write_log_tsv(
                log_tsv, annotation_id,
                stderr if stderr else "busco_failed"
            )
            return 1

        # ------------------------------------------------------------------
        # Step 6: Extract BUSCO transcript sequences + filter GFF
        # ------------------------------------------------------------------
        logger.info("STEP 6: Extract BUSCO gene sequences and filter reference GFF")
        try:
            extract_busco_sequences(
                busco_output,
                alias_gff_decompressed,
                fasta_decompressed,
                annotation_id,
                output_dir,
            )
        except Exception as e:
            # Non-fatal: log but continue to selenoprofiles
            logger.warning(f"BUSCO sequence extraction error (non-fatal): {e}")

        # ------------------------------------------------------------------
        # Step 7: Run selenoprofiles on the genome
        # ------------------------------------------------------------------
        logger.info("STEP 7: Run selenoprofiles")
        seleno_outdir = work_dir / f"selenoprofiles_{annotation_id}"
        seleno_outdir.mkdir(parents=True, exist_ok=True)

        ok, _, stderr = run_shell_script(
            seleno_script,
            [
                str(fasta_decompressed),
                str(alias_gff_decompressed),
                str(seleno_outdir),
            ],
            "run_selenoprofiles",
        )
        if not ok:
            # Non-fatal: record in log but still write BUSCO result
            logger.warning(
                f"Selenoprofiles failed for {annotation_id} (non-fatal): {stderr}"
            )
        else:
            # Copy selenoprofiles outputs to per-annotation output dir.
            # all_predictions.gtf is the raw GTF produced by selenoprofiles;
            # we store it as All_predictions.gtf to match the expected output name.
            for src_name, dst_name in (
                ("all_predictions.gtf",               "All_predictions.gtf"),
                ("Selenoprofiles.fasta",               "Selenoprofiles.fasta"),
                ("Selenoprofiles_annotation_result.csv", "Selenoprofiles_annotation_result.csv"),
            ):
                src = seleno_outdir / src_name
                if src.exists():
                    shutil.copy2(src, output_dir / dst_name)
                    logger.info(f"Copied {src_name} -> {dst_name} in {output_dir}")
                else:
                    logger.warning(f"Expected selenoprofiles output missing: {src_name}")

        # ------------------------------------------------------------------
        # Step 8: Parse BUSCO results and write result TSV fragment
        # ------------------------------------------------------------------
        logger.info("STEP 8: Parse BUSCO results and write result fragment")
        busco_results = parse_busco_results(busco_output)
        seleno_results = parse_selenoprofiles_results(seleno_outdir / "Selenoprofiles_annotation_result.csv")
        busco_results.update(seleno_results)
        write_result_tsv(result_tsv, annotation_id, busco_results)

        logger.info(f"Analysis complete for {annotation_id}")
        return 0

    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        write_log_tsv(log_tsv, annotation_id, f"unexpected_error: {e}")
        return 1

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
        logger.info(f"Cleaned up working directory: {work_dir}")


if __name__ == "__main__":
    sys.exit(main())

