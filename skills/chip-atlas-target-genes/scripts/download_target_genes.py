"""
Download and parse a ChIP-Atlas target-genes TSV.

Layout
------
Files live at
    https://chip-atlas.dbcls.jp/data/{genome}/target/{Protein}.{Distance}.tsv

with one row per gene and columns:
    Target_genes
    {Protein}|Average
    {SRX_ID}|{cell_type}  (repeats)
    STRING

Inputs : protein name, genome, distance window
Outputs: (summary DataFrame, wide-format experiment DataFrame)

Example
-------
    python download_target_genes.py --check TP53 --genome hg38 --distance 5
"""

from __future__ import annotations

import argparse
import io
import re
import sys
from typing import Tuple

try:
    import pandas as pd
    import requests
except ImportError:
    sys.exit("Install with: pip install pandas requests")


BASE_URL = "https://chip-atlas.dbcls.jp/data"

SUPPORTED_GENOMES = (
    "hg38", "hg19", "mm10", "mm9", "rn6",
    "dm6", "dm3", "ce11", "ce10", "sacCer3",
)
VALID_DISTANCES = (1, 5, 10)

EXPERIMENT_COL_RE = re.compile(r"^[SED]RX\d+\|")


def _build_url(protein: str, genome: str, distance: int) -> str:
    return f"{BASE_URL}/{genome}/target/{protein}.{distance}.tsv"


def check_antigen_available(
    protein: str,
    genome: str = "hg38",
    distance: int = 5,
    request_timeout: int = 15,
) -> bool:
    """Return True if the TSV exists on the data server."""
    if genome not in SUPPORTED_GENOMES:
        print(f"  ! invalid genome '{genome}'. Supported: {SUPPORTED_GENOMES}")
        return False
    if distance not in VALID_DISTANCES:
        print(f"  ! invalid distance {distance}. Supported: {VALID_DISTANCES}")
        return False

    url = _build_url(protein, genome, distance)
    try:
        response = requests.head(url, timeout=request_timeout, allow_redirects=True)
    except requests.RequestException as exc:
        print(f"  ! network error checking {url}: {exc}")
        return False

    if response.status_code == 200:
        return True
    if response.status_code == 404:
        print(f"  ! no target-gene data for '{protein}' ({genome}, ±{distance}kb)")
        print("    - protein symbols are case-sensitive (use 'TP53', not 'tp53')")
        print("    - histone marks (H3K4me3, H3K27ac, …) are NOT served by this endpoint")
        print(f"    - check https://chip-atlas.org/target_genes for available antigens")
        return False
    print(f"  ! unexpected HTTP {response.status_code} for {url}")
    return False


def download_target_genes(
    protein: str,
    genome: str = "hg38",
    distance: int = 5,
    request_timeout: int = 120,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch the TSV and return (summary, experiment_wide) DataFrames."""
    url = _build_url(protein, genome, distance)
    print(f"  downloading {url}")
    try:
        response = requests.get(url, timeout=request_timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc

    body = response.content.decode("utf-8")
    print(f"  received {len(body) / 1_048_576:.1f} MB")

    frame = pd.read_csv(io.StringIO(body), sep="\t")
    if frame.empty:
        raise ValueError(f"{url} returned an empty table")

    summary, experiments = parse_target_genes_tsv(frame, protein)
    print(f"  ✓ parsed {len(summary)} genes across {experiments.shape[1] - 1} experiments")
    return summary, experiments


def parse_target_genes_tsv(frame: pd.DataFrame, protein: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    gene_column = "Target_genes"
    if gene_column not in frame.columns:
        raise ValueError(
            f"Missing expected '{gene_column}' column. First columns: {list(frame.columns[:5])}"
        )

    avg_column = f"{protein}|Average"
    if avg_column not in frame.columns:
        candidates = [col for col in frame.columns if col.lower().endswith("|average")]
        if not candidates:
            raise ValueError(
                f"Could not find an Average column. Header sample: {list(frame.columns[:5])}"
            )
        avg_column = candidates[0]
        print(f"  using fallback Average column: {avg_column}")

    experiment_cols = [
        col for col in frame.columns if EXPERIMENT_COL_RE.match(col) and col != avg_column
    ]

    has_string = "STRING" in frame.columns

    experiment_df = frame[[gene_column] + experiment_cols].rename(columns={gene_column: "gene"}).copy()

    experiment_values = frame[experiment_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    summary = pd.DataFrame({
        "gene": frame[gene_column].astype(str),
        "avg_score": pd.to_numeric(frame[avg_column], errors="coerce").fillna(0),
        "string_score": pd.to_numeric(frame["STRING"], errors="coerce").fillna(0) if has_string else 0,
        "num_experiments": len(experiment_cols),
    })
    summary["num_bound"] = (experiment_values > 0).sum(axis=1).values
    summary["binding_rate"] = summary["num_bound"] / max(len(experiment_cols), 1)
    summary["max_score"] = experiment_values.max(axis=1).values

    summary = summary.sort_values("avg_score", ascending=False).reset_index(drop=True)

    cell_types = {col.split("|", 1)[1] for col in experiment_cols if "|" in col}
    print(f"  unique cell types: {len(cell_types)}")
    return summary, experiment_df


def get_cell_types(experiment_df: pd.DataFrame) -> list:
    """Return sorted unique cell types parsed from the experiment column headers."""
    types = set()
    for column in experiment_df.columns:
        if column == "gene":
            continue
        parts = column.split("|", 1)
        if len(parts) == 2:
            types.add(parts[1])
    return sorted(types)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--genome", default="hg38", choices=SUPPORTED_GENOMES)
    parser.add_argument("--distance", type=int, default=5, choices=VALID_DISTANCES)
    parser.add_argument("--check", action="store_true", help="HEAD-check availability without downloading")
    parser.add_argument("--preview", action="store_true", help="Print the top 5 rows after download")
    args = parser.parse_args(argv)

    if args.check:
        available = check_antigen_available(args.protein, args.genome, args.distance)
        print("available" if available else "not available")
        return 0 if available else 1

    summary, _ = download_target_genes(args.protein, args.genome, args.distance)
    if args.preview:
        print(summary.head().to_string(index=False))
    print(f"✓ {len(summary)} target genes loaded for {args.protein} ({args.genome}, ±{args.distance}kb)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
