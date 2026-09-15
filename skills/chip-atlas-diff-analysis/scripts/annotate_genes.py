"""
Annotate the top differential regions with UCSC RefSeq gene symbols.

The annotator queries the UCSC REST endpoint
(`https://api.genome.ucsc.edu/getData/track`) for each region. If the
region overlaps a gene, the annotation reads "SYMBOL"; otherwise it
reports the nearest gene within 5 kb as "SYMBOL (+1.3kb)".

Failures are absorbed gracefully — after three consecutive HTTP errors the
annotator stops and returns whatever it gathered. Statistical results are
unaffected.

Inputs : DataFrame with `chrom`, `chromStart`, `chromEnd`
Outputs: pandas Series of annotations, aligned to the input index
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

try:
    import pandas as pd
    import requests
except ImportError:
    sys.exit("Install with: pip install pandas requests")


UCSC_URL = "https://api.genome.ucsc.edu/getData/track"
FLANK_BP = 5000


def annotate_nearest_genes(
    frame: pd.DataFrame,
    genome: str = "hg38",
    max_rows: int = 50,
    request_timeout: int = 10,
) -> pd.Series:
    """Return gene-symbol annotations for the first `max_rows` of `frame`."""
    annotations = pd.Series("", index=frame.index, dtype=object)
    subset = frame.head(max_rows)
    if subset.empty:
        return annotations

    print(f"   Annotating top {len(subset)} regions via UCSC REST ({genome})…")
    failures = 0
    for index, row in subset.iterrows():
        if failures >= 3:
            print("   ! UCSC API failing repeatedly — skipping further annotations")
            break

        try:
            symbol = _gene_for_region(
                row["chrom"], int(row["chromStart"]), int(row["chromEnd"]),
                genome=genome, request_timeout=request_timeout,
            )
            if symbol:
                annotations.at[index] = symbol
        except Exception:
            failures += 1
    annotated = int((annotations != "").sum())
    print(f"   annotated {annotated}/{len(subset)} regions with nearest gene")
    return annotations


def _gene_for_region(
    chrom: str,
    start: int,
    end: int,
    genome: str,
    request_timeout: int,
) -> str:
    overlap = _gene_symbols(chrom, start, end, genome=genome, request_timeout=request_timeout)
    if overlap:
        return ", ".join(sorted(overlap))

    flank_start = max(0, start - FLANK_BP)
    flank_end = end + FLANK_BP
    records = _gene_records(chrom, flank_start, flank_end, genome=genome, request_timeout=request_timeout)
    if not records:
        return ""

    mid = (start + end) / 2
    closest = None
    closest_distance = float("inf")
    for record in records:
        gene_mid = (record.get("txStart", 0) + record.get("txEnd", 0)) / 2
        distance = abs(mid - gene_mid)
        if record.get("name2") and distance < closest_distance:
            closest = record["name2"]
            closest_distance = distance

    if closest is None:
        return ""
    if closest_distance == 0:
        return closest
    return f"{closest} ({closest_distance / 1000:+.1f}kb)"


def _gene_symbols(chrom: str, start: int, end: int, genome: str, request_timeout: int) -> set:
    return {
        record["name2"]
        for record in _gene_records(chrom, start, end, genome=genome, request_timeout=request_timeout)
        if record.get("name2")
    }


def _gene_records(chrom: str, start: int, end: int, genome: str, request_timeout: int) -> list:
    params = {
        "genome": genome,
        "track": "ncbiRefSeq",
        "chrom": chrom,
        "start": int(start),
        "end": int(end),
    }
    response = requests.get(UCSC_URL, params=params, timeout=request_timeout)
    response.raise_for_status()
    return response.json().get("ncbiRefSeq", []) or []


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="Parsed CSV with chrom/chromStart/chromEnd")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--out", default="annotated_top.csv")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.csv)
    annotations = annotate_nearest_genes(frame, genome=args.genome, max_rows=args.top)
    out_frame = frame.head(args.top).copy()
    out_frame["nearest_gene"] = annotations.reindex(out_frame.index)
    out_frame.to_csv(args.out, index=False)
    print(f"✓ Wrote {len(out_frame)} annotated rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
