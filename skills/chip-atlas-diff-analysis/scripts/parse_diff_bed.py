"""
Parse the 8-column BED produced by the ChIP-Atlas Diff Analysis job.

Each row of `wabi_result.bed` is:
    chrom  chromStart  chromEnd  counts_a  counts_b  logFC  pvalue  qvalue

This module:
* validates the column count (tolerating the BED9 IGV companion file as a
  fallback),
* coerces numeric columns,
* derives `region_size`, `direction`, `significant`, per-group mean counts,
  and a BED-compatible `score`,
* sorts ascending by qvalue with ties broken by `|logFC|` then position.

Inputs : a path to `wabi_result.bed`
Outputs: a pandas DataFrame with the original 8 columns plus the derived ones.

Example
-------
    python parse_diff_bed.py --bed raw_results/wabi_result.bed --preview
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas numpy")


PRIMARY_COLUMNS = [
    "chrom", "chromStart", "chromEnd",
    "counts_a", "counts_b",
    "logFC", "pvalue", "qvalue",
]


def parse_diff_bed(bed_path: str) -> pd.DataFrame:
    """Read and decorate the diff-analysis BED at `bed_path`."""
    try:
        raw = pd.read_csv(bed_path, sep="\t", header=None, comment="#", dtype=str)
    except pd.errors.ParserError as exc:
        raise RuntimeError(f"Could not parse {bed_path}: {exc}") from exc

    mask = ~raw.iloc[:, 0].astype(str).str.lower().str.startswith(("browser", "track"))
    raw = raw[mask].reset_index(drop=True)

    if raw.empty:
        raise RuntimeError(f"BED file {bed_path} is empty")

    column_count = raw.shape[1]
    if column_count == len(PRIMARY_COLUMNS):
        raw.columns = PRIMARY_COLUMNS
        frame = _coerce_primary(raw)
    elif column_count == 9:
        raw.columns = ["chrom", "chromStart", "chromEnd", "name", "score",
                       "strand", "thickStart", "thickEnd", "itemRgb"]
        frame = _parse_igv_bed(raw)
    else:
        raise RuntimeError(
            f"Unexpected column count {column_count} in {bed_path}. "
            f"Expected 8 (primary BED) or 9 (IGV BED)."
        )

    frame = frame.sort_values(
        ["qvalue", "_abs_logFC", "chrom", "chromStart"],
        ascending=[True, False, True, True],
    ).drop(columns=["_abs_logFC"]).reset_index(drop=True)
    return frame


def summarize_regions(frame: pd.DataFrame) -> None:
    """Print a short text summary of the parsed regions."""
    n_total = len(frame)
    print(f"   Regions parsed: {n_total}")
    if "direction" in frame.columns:
        a = int((frame["direction"] == "A_enriched").sum())
        b = int((frame["direction"] == "B_enriched").sum())
        print(f"   A-enriched: {a}; B-enriched: {b}")
    if "significant" in frame.columns:
        print(f"   Significant (FDR<0.05): {int(frame['significant'].sum())}")
    if "region_size" in frame.columns and n_total > 0:
        print(f"   Median region size: {int(frame['region_size'].median())} bp")
    if "qvalue" in frame.columns and n_total > 0:
        print(f"   Min Q-value: {frame['qvalue'].min():.2e}")


def _coerce_primary(frame: pd.DataFrame) -> pd.DataFrame:
    frame["chromStart"] = pd.to_numeric(frame["chromStart"], errors="coerce")
    frame["chromEnd"] = pd.to_numeric(frame["chromEnd"], errors="coerce")
    frame["logFC"] = pd.to_numeric(frame["logFC"], errors="coerce")
    frame["pvalue"] = pd.to_numeric(frame["pvalue"], errors="coerce")
    frame["qvalue"] = pd.to_numeric(frame["qvalue"], errors="coerce")

    frame["region_size"] = (frame["chromEnd"] - frame["chromStart"]).astype("Int64")
    frame["direction"] = np.where(frame["logFC"].fillna(0) > 0, "A_enriched", "B_enriched")
    frame["significant"] = frame["qvalue"].fillna(1.0) < 0.05
    frame["mean_count_a"] = frame["counts_a"].apply(_mean_of_csv)
    frame["mean_count_b"] = frame["counts_b"].apply(_mean_of_csv)

    safe_q = frame["qvalue"].clip(lower=1e-300)
    frame["score"] = np.minimum(-np.log10(safe_q) * 100, 1000).fillna(0).astype(int)
    frame["_abs_logFC"] = frame["logFC"].abs()
    return frame


def _parse_igv_bed(frame: pd.DataFrame) -> pd.DataFrame:
    from urllib.parse import unquote

    frame["chromStart"] = pd.to_numeric(frame["chromStart"], errors="coerce")
    frame["chromEnd"] = pd.to_numeric(frame["chromEnd"], errors="coerce")
    frame["region_size"] = (frame["chromEnd"] - frame["chromStart"]).astype("Int64")
    frame["counts_a"] = ""
    frame["counts_b"] = ""

    def _extract(field: str, key: str) -> float:
        try:
            decoded = unquote(str(field))
            for piece in decoded.split(";"):
                if piece.startswith(f"{key}="):
                    return float(piece.split("=", 1)[1])
        except (ValueError, IndexError):
            pass
        return float("nan")

    frame["logFC"] = frame["name"].apply(lambda x: _extract(x, "LogFC"))
    frame["pvalue"] = frame["name"].apply(lambda x: _extract(x, "P-value"))
    frame["qvalue"] = frame["name"].apply(lambda x: _extract(x, "Q-value"))

    frame["direction"] = frame["itemRgb"].apply(_direction_from_rgb)
    frame["significant"] = frame["qvalue"].fillna(1.0) < 0.05
    frame["mean_count_a"] = float("nan")
    frame["mean_count_b"] = float("nan")
    frame["score"] = pd.to_numeric(frame.get("score", 0), errors="coerce").fillna(0).astype(int)
    frame["_abs_logFC"] = frame["logFC"].abs()
    return frame


def _mean_of_csv(value) -> float:
    try:
        parts = [float(p) for p in str(value).split(",") if p]
        return sum(parts) / len(parts) if parts else 0.0
    except ValueError:
        return 0.0


def _direction_from_rgb(rgb: Optional[str]) -> str:
    try:
        components = [int(x) for x in str(rgb).split(",")[:3]]
        if len(components) < 3:
            return "unknown"
        r, _, b = components
        if r > b:
            return "A_enriched"
        if b > r:
            return "B_enriched"
        return "unchanged"
    except (ValueError, TypeError):
        return "unknown"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bed", required=True, help="Path to wabi_result.bed")
    parser.add_argument("--preview", action="store_true", help="Print the first 5 rows")
    args = parser.parse_args(argv)
    frame = parse_diff_bed(args.bed)
    summarize_regions(frame)
    if args.preview:
        print(frame.head().to_string(index=False))
    print(f"✓ Parsed {len(frame)} regions from {args.bed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
