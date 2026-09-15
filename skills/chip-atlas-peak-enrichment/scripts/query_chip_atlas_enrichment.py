"""
Submit, poll, and retrieve ChIP-Atlas Enrichment Analysis jobs.

Submits a gene list to the WABI endpoint, polls for completion, and
parses the headerless 11-column TSV result into a tidy DataFrame with
the columns expected by the rest of this skill.

Inputs
------
* `gene_list`   : list of gene symbols (>=3)
* `genome`      : one of hg38, hg19, mm10, mm9, rn6, dm6, dm3, ce11, ce10, sacCer3
* `antigen_class`, `cell_class`, `threshold`, `distance_up`, `distance_down`

Outputs
-------
* `request_id`              from `submit_enrichment_job`
* raw API DataFrame         from `retrieve_results`
* contract-compliant frame  from `parse_api_results`

Example
-------
    python query_chip_atlas_enrichment.py --self-test
"""

from __future__ import annotations

import argparse
import io
import re
import sys
import time
from typing import List, Optional, Tuple

try:
    import numpy as np
    import pandas as pd
    import requests
except ImportError:
    sys.exit("Install with: pip install pandas numpy requests")


WABI_URL = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"

SUPPORTED_GENOMES = (
    "hg38", "hg19", "mm10", "mm9", "rn6",
    "dm6", "dm3", "ce11", "ce10", "sacCer3",
)

VALID_ANTIGEN_CLASSES = (
    "TFs and others", "Histone", "RNA polymerase", "Input control",
    "ATAC-Seq", "DNase-seq", "Bisulfite-Seq", "Annotation tracks",
)

VALID_THRESHOLDS = {50, 100, 200, 500}

RAW_TSV_COLUMNS = [
    "experiment_id", "antigen_class", "antigen",
    "cell_type_class", "cell_type",
    "num_peaks", "overlap_a", "overlap_b",
    "log10_pvalue", "log10_qvalue", "fold_enrichment",
]


def _sanitize_symbols(symbols: List[str]) -> List[str]:
    sanitized: List[str] = []
    for symbol in symbols:
        cleaned = "".join(c if (c.isalnum() or c == "_") else "_" for c in str(symbol).strip())
        if cleaned:
            sanitized.append(cleaned)
    return sanitized


def submit_enrichment_job(
    gene_list: List[str],
    genome: str = "hg38",
    antigen_class: str = "TFs and others",
    cell_class: str = "All cell types",
    threshold: int = 50,
    distance_up: int = 5000,
    distance_down: int = 5000,
    title: str = "chip_atlas_enrichment",
    description_a: str = "input_genes",
    description_b: str = "RefSeq_background",
    request_timeout: int = 60,
) -> str:
    if len(gene_list) < 3:
        raise ValueError(f"Need at least 3 gene symbols (got {len(gene_list)})")
    if genome not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'. Supported: {SUPPORTED_GENOMES}")
    if antigen_class not in VALID_ANTIGEN_CLASSES:
        raise ValueError(
            f"antigen_class '{antigen_class}' is not in the accepted vocabulary "
            f"{VALID_ANTIGEN_CLASSES}"
        )
    if threshold not in VALID_THRESHOLDS:
        raise ValueError(f"threshold must be one of {sorted(VALID_THRESHOLDS)}")
    if not cell_class or not str(cell_class).strip():
        raise ValueError("cell_class must be non-empty (use 'All cell types' to skip filtering)")

    sanitized = _sanitize_symbols(gene_list)

    payload = {
        "address": "",
        "format": "json",
        "result": "www",
        "genome": genome,
        "antigenClass": antigen_class,
        "cellClass": cell_class,
        "threshold": str(threshold),
        "typeA": "gene",
        "bedAFile": "\n".join(sanitized),
        "typeB": "refseq",
        "bedBFile": "empty",
        "permTime": "1",
        "title": title,
        "descriptionA": description_a,
        "descriptionB": description_b,
        "distanceUp": str(distance_up),
        "distanceDown": str(distance_down),
    }

    try:
        response = requests.post(WABI_URL, data=payload, timeout=request_timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(
            f"Could not reach the ChIP-Atlas API ({exc}). "
            "Check connectivity to https://dtn1.ddbj.nig.ac.jp."
        ) from exc

    try:
        body = response.json()
    except ValueError as exc:
        raise RuntimeError(
            f"API did not return JSON. First 500 chars:\n{response.text.strip()[:500]}"
        ) from exc

    request_id = body.get("requestId")
    if not request_id:
        error_msg = body.get("error-message") or body.get("Message") or "Unknown error"
        raise RuntimeError(f"API rejected the submission: {error_msg}")

    print(f"✓ Enrichment job submitted: request_id={request_id}")
    return request_id


def poll_job_status(
    request_id: str,
    timeout: int = 600,
    poll_interval: int = 15,
) -> str:
    deadline = time.time() + timeout
    url = f"{WABI_URL}{request_id}?info=status"

    while time.time() < deadline:
        elapsed_minutes = (timeout - (deadline - time.time())) / 60
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"   ! poll error ({exc}); retrying in {poll_interval}s")
            time.sleep(poll_interval)
            continue

        status = _extract_status(response.text)
        if status == "finished":
            print(f"✓ Job {request_id} finished after {elapsed_minutes:.1f} min")
            return status
        if status in {"error", "failed"}:
            raise RuntimeError(
                f"ChIP-Atlas job {request_id} failed. Response:\n{response.text}"
            )

        print(f"   ... still running ({elapsed_minutes:.1f} min)")
        time.sleep(poll_interval)

    raise RuntimeError(
        f"Job {request_id} did not finish within {timeout} seconds. "
        "Increase timeout, narrow cell_class, or rerun later."
    )


def _extract_status(body: str) -> Optional[str]:
    for line in body.splitlines():
        line = line.strip().lower()
        if line.startswith("status"):
            parts = re.split(r"[:\t]\s*", line, maxsplit=1)
            if len(parts) == 2:
                return parts[1].strip()
    return None


def retrieve_results(request_id: str, request_timeout: int = 120) -> pd.DataFrame:
    url = f"{WABI_URL}{request_id}?info=result&format=tsv"
    try:
        response = requests.get(url, timeout=request_timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"Failed to download results for {request_id}: {exc}") from exc

    if not response.text.strip():
        raise RuntimeError(f"Empty TSV for {request_id}. The job may need a resubmit.")

    frame = pd.read_csv(io.StringIO(response.text), sep="\t", header=None, names=RAW_TSV_COLUMNS)
    print(f"✓ Retrieved {len(frame)} rows from the API")
    return frame


def parse_api_results(api_df: pd.DataFrame) -> pd.DataFrame:
    """Convert the raw 11-column TSV into the skill's standard contract."""
    out = pd.DataFrame()
    out["experiment_id"] = api_df["experiment_id"]
    out["antigen"] = api_df["antigen"]
    out["cell_type"] = api_df["cell_type"]
    out["antigen_class"] = api_df["antigen_class"]
    out["cell_type_class"] = api_df["cell_type_class"]
    out["fold_enrichment"] = pd.to_numeric(api_df["fold_enrichment"], errors="coerce").fillna(0)

    log_p = pd.to_numeric(api_df["log10_pvalue"], errors="coerce").fillna(0)
    log_q = pd.to_numeric(api_df["log10_qvalue"], errors="coerce").fillna(0)
    out["p_value"] = np.where(log_p < 0, np.power(10.0, log_p), 1.0)
    out["q_value"] = np.where(log_q < 0, np.power(10.0, log_q), 1.0)
    out["significant"] = out["q_value"] < 0.05

    overlap_a = api_df["overlap_a"].apply(_parse_fraction)
    out["regions_with_overlaps"] = overlap_a.apply(lambda pair: pair[0])
    out["total_regions"] = overlap_a.apply(lambda pair: pair[1])
    out["overlap_rate"] = out["regions_with_overlaps"] / out["total_regions"].replace(0, 1)

    overlap_b = api_df["overlap_b"].apply(_parse_fraction)
    out["background_overlaps"] = overlap_b.apply(lambda pair: pair[0])
    out["background_size"] = overlap_b.apply(lambda pair: pair[1])

    out["total_peaks"] = pd.to_numeric(api_df["num_peaks"], errors="coerce").fillna(0).astype(int)

    out = out.sort_values("q_value", ascending=True).reset_index(drop=True)
    return out


def _parse_fraction(value) -> Tuple[int, int]:
    try:
        numerator, denominator = str(value).split("/")
        return int(numerator), int(denominator)
    except (ValueError, AttributeError):
        return 0, 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--genes", nargs="+")
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--antigen-class", default="TFs and others", choices=list(VALID_ANTIGEN_CLASSES))
    parser.add_argument("--cell-class", default="All cell types")
    parser.add_argument("--threshold", type=int, default=50, choices=sorted(VALID_THRESHOLDS))
    parser.add_argument("--distance-up", type=int, default=5000)
    parser.add_argument("--distance-down", type=int, default=5000)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args(argv)

    if args.self_test:
        sample = ["TP53", "CDKN1A", "BAX", "MDM2", "BBC3"]
        sanitized = _sanitize_symbols(sample)
        assert sanitized == sample, "Sanitizer should pass plain alphanumerics"
        print("✓ Self-test passed: sanitizer round-trips standard symbols")
        return 0

    if not args.genes:
        parser.error("Provide --genes or --self-test")
    request_id = submit_enrichment_job(
        gene_list=args.genes,
        genome=args.genome,
        antigen_class=args.antigen_class,
        cell_class=args.cell_class,
        threshold=args.threshold,
        distance_up=args.distance_up,
        distance_down=args.distance_down,
    )
    poll_job_status(request_id, timeout=args.timeout)
    raw = retrieve_results(request_id)
    parsed = parse_api_results(raw)
    print(f"✓ Parsed {len(parsed)} experiments; significant: {int(parsed['significant'].sum())}")
    print(parsed[["antigen", "cell_type", "fold_enrichment", "q_value"]].head())
    return 0


if __name__ == "__main__":
    sys.exit(main())
