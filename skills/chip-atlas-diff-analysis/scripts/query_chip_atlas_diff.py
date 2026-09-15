"""
Submit, poll, and retrieve ChIP-Atlas Diff Analysis jobs.

The ChIP-Atlas WABI endpoint accepts two newline-separated lists of
experiment accessions and runs either edgeR (DPR mode, for ChIP/ATAC/DNase-
seq) or metilene (DMR mode, for Bisulfite-seq). The result is a ZIP
archive containing an 8-column BED, an IGV-friendly BED9 companion, an
IGV session XML, and a log.

Inputs
------
* `experiments_a`, `experiments_b` : lists of SRX/ERX/DRX/GSM accessions
  (>=2 per group)
* `genome`         : one of hg38, hg19, mm10, mm9, rn6, dm6, dm3, ce11,
                     ce10, sacCer3
* `analysis_type`  : 'diffbind' (default) or 'dmr'

Outputs
-------
* A `request_id` string returned by `submit_diff_job`
* A dict of extracted file paths returned by `retrieve_zip_results`

Example
-------
    python query_chip_atlas_diff.py --self-test
"""

from __future__ import annotations

import argparse
import io
import os
import re
import sys
import time
import zipfile
from typing import Dict, Iterable, List, Optional

try:
    import requests
except ImportError:
    sys.exit("Install with: pip install requests")


WABI_URL = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"

SUPPORTED_GENOMES = (
    "hg38", "hg19", "mm10", "mm9", "rn6",
    "dm6", "dm3", "ce11", "ce10", "sacCer3",
)

ANALYSIS_TYPES = {
    "diffbind": "Differential peak regions via edgeR (ChIP/ATAC/DNase-seq)",
    "dmr": "Differentially methylated regions via metilene (Bisulfite-seq)",
}

ACCESSION_RE = re.compile(r"^(SRX|ERX|DRX|GSM)\d+$")


def _validate_accessions(accessions: Iterable[str], group_label: str) -> List[str]:
    cleaned: List[str] = []
    for acc in accessions:
        acc = str(acc).strip()
        if not acc:
            continue
        if not ACCESSION_RE.match(acc):
            print(f"   ! {group_label}: '{acc}' does not match SRX/ERX/DRX/GSM — keeping anyway")
        cleaned.append(acc)
    if len(cleaned) < 2:
        raise ValueError(f"Group {group_label} needs at least 2 experiment IDs (got {len(cleaned)})")
    return cleaned


def submit_diff_job(
    experiments_a: List[str],
    experiments_b: List[str],
    genome: str = "hg38",
    analysis_type: str = "diffbind",
    title: str = "chip_atlas_diff_run",
    description_a: str = "group_A",
    description_b: str = "group_B",
    request_timeout: int = 60,
) -> str:
    """Submit a diff job and return its request ID."""
    if genome not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'. Supported: {SUPPORTED_GENOMES}")
    if analysis_type not in ANALYSIS_TYPES:
        raise ValueError(
            f"Unknown analysis_type '{analysis_type}'. Expected one of {list(ANALYSIS_TYPES)}"
        )

    group_a = _validate_accessions(experiments_a, "A")
    group_b = _validate_accessions(experiments_b, "B")

    payload = {
        "address": "",
        "format": "text",
        "result": "www",
        "genome": genome,
        "antigenClass": analysis_type,
        "cellClass": "empty",
        "threshold": "50",
        "typeA": "srx",
        "bedAFile": "\n".join(group_a),
        "typeB": "srx",
        "bedBFile": "\n".join(group_b),
        "permTime": "1",
        "title": title,
        "descriptionA": description_a,
        "descriptionB": description_b,
    }

    try:
        response = requests.post(WABI_URL, data=payload, timeout=request_timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(
            f"Could not reach the ChIP-Atlas API ({exc}). "
            "Check connectivity to https://dtn1.ddbj.nig.ac.jp before retrying."
        ) from exc

    request_id = _extract_request_id(response.text)
    if not request_id:
        snippet = response.text.strip()[:500]
        raise RuntimeError(
            f"Could not parse a request ID from the API response. First 500 chars:\n{snippet}"
        )

    print(f"✓ Diff job submitted: request_id={request_id}")
    return request_id


def _extract_request_id(body: str) -> Optional[str]:
    for line in body.splitlines():
        line = line.strip()
        if not line.lower().startswith("requestid"):
            continue
        parts = re.split(r"[:\t]\s*", line, maxsplit=1)
        if len(parts) == 2 and parts[1].strip():
            return parts[1].strip()
    return None


def poll_job_status(
    request_id: str,
    timeout: int = 900,
    poll_interval: int = 15,
) -> str:
    """Block until the job finishes, fails, or `timeout` seconds elapse."""
    deadline = time.time() + timeout
    url = f"{WABI_URL}{request_id}?info=status"

    while time.time() < deadline:
        elapsed = (timeout - (deadline - time.time())) / 60
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"   ! Poll error ({exc}); retrying in {poll_interval}s")
            time.sleep(poll_interval)
            continue

        status = _extract_status(response.text)
        if status == "finished":
            print(f"✓ Job {request_id} finished after {elapsed:.1f} min")
            return status
        if status in {"error", "failed"}:
            raise RuntimeError(
                f"ChIP-Atlas job {request_id} failed. Server response:\n{response.text}"
            )

        print(f"   ... still running ({elapsed:.1f} min)")
        time.sleep(poll_interval)

    raise RuntimeError(
        f"Job {request_id} did not finish within {timeout} seconds. "
        "Increase timeout or reduce the number of experiments per group."
    )


def _extract_status(body: str) -> Optional[str]:
    for line in body.splitlines():
        line = line.strip().lower()
        if line.startswith("status"):
            parts = re.split(r"[:\t]\s*", line, maxsplit=1)
            if len(parts) == 2:
                return parts[1].strip()
    return None


def retrieve_zip_results(
    request_id: str,
    output_dir: str,
    fetch_timeout: int = 300,
) -> Dict[str, str]:
    """Download the result ZIP and extract it into `output_dir/raw_results/`."""
    url = f"{WABI_URL}{request_id}?info=result&format=zip"

    try:
        response = requests.get(url, timeout=fetch_timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"Failed to download results for {request_id}: {exc}") from exc

    if not response.content:
        raise RuntimeError(
            f"The API returned an empty body for {request_id}. "
            "Wait a few seconds and try again — the result may not be fully staged."
        )

    raw_dir = os.path.join(output_dir, "raw_results")
    os.makedirs(raw_dir, exist_ok=True)

    try:
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            archive.extractall(raw_dir)
    except zipfile.BadZipFile as exc:
        raise RuntimeError(
            f"Invalid ZIP returned for {request_id}. The job may still be finalizing."
        ) from exc

    paths = _index_extracted_files(raw_dir)
    if "bed" not in paths:
        listing = []
        for root, _, files in os.walk(raw_dir):
            listing.extend(files)
        raise RuntimeError(
            f"No main BED file found in the result archive. Extracted files: {listing}"
        )
    print(f"✓ Retrieved {len(paths)} result files into {raw_dir}")
    return paths


def _index_extracted_files(root: str) -> Dict[str, str]:
    index: Dict[str, str] = {}
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            path = os.path.join(dirpath, name)
            if name.endswith(".igv.bed"):
                index["igv_bed"] = path
            elif name.endswith(".bed") and "igv" not in name:
                index["bed"] = path
            elif name.endswith(".log"):
                index["log"] = path
            elif name.endswith(".xml"):
                index["igv_xml"] = path
            elif name.lower() == "readme.md":
                index["readme"] = path
    return index


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--self-test", action="store_true", help="Validate inputs without contacting the API")
    parser.add_argument("--experiments-a", nargs="+", help="Group A accessions")
    parser.add_argument("--experiments-b", nargs="+", help="Group B accessions")
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--analysis-type", default="diffbind", choices=list(ANALYSIS_TYPES))
    parser.add_argument("--output-dir", default="diff_results")
    parser.add_argument("--description-a", default="group_A")
    parser.add_argument("--description-b", default="group_B")
    parser.add_argument("--title", default="chip_atlas_diff_run")
    parser.add_argument("--timeout", type=int, default=900, help="Poll timeout in seconds")
    args = parser.parse_args(argv)

    if args.self_test:
        ids_a = ["SRX5865959", "SRX5865960"]
        ids_b = ["SRX5865967", "SRX5865968"]
        _validate_accessions(ids_a, "A")
        _validate_accessions(ids_b, "B")
        print("✓ Self-test passed: accession validators accept SRX-style IDs")
        return 0

    if not args.experiments_a or not args.experiments_b:
        parser.error("Provide --experiments-a and --experiments-b, or use --self-test")

    request_id = submit_diff_job(
        experiments_a=args.experiments_a,
        experiments_b=args.experiments_b,
        genome=args.genome,
        analysis_type=args.analysis_type,
        title=args.title,
        description_a=args.description_a,
        description_b=args.description_b,
    )
    poll_job_status(request_id, timeout=args.timeout)
    paths = retrieve_zip_results(request_id, args.output_dir)
    print(f"✓ Files: {paths}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
