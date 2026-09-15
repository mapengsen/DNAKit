"""
Genome-aware filtering helpers for diff-analysis result tables.

Two layers of filtering are exposed:

1. `filter_standard_chromosomes` / `filter_min_region_size`
   These implement the default QC trims (drop non-standard contigs,
   drop regions smaller than 10 bp). They are called automatically by
   `run_diff_workflow.py`.

2. `filter_regions` is a user-facing filter that accepts qvalue, |logFC|,
   region-size, direction, score, and chromosome predicates.

Inputs : a parsed DataFrame (see `parse_diff_bed.py`)
Outputs: a filtered DataFrame.

Example
-------
    python filter_regions.py --csv diff_regions_all.csv --genome hg38 --max-qvalue 0.05
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


def _autosomes(count: int) -> set:
    return {f"chr{i}" for i in range(1, count + 1)}


STANDARD_CHROMS = {
    "hg38": _autosomes(22) | {"chrX", "chrY", "chrM"},
    "hg19": _autosomes(22) | {"chrX", "chrY", "chrM"},
    "mm10": _autosomes(19) | {"chrX", "chrY", "chrM"},
    "mm9":  _autosomes(19) | {"chrX", "chrY", "chrM"},
    "rn6":  _autosomes(20) | {"chrX", "chrY", "chrM"},
    "dm6": {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"},
    "dm3": {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"},
    "ce11": {"chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"},
    "ce10": {"chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"},
    "sacCer3": {"chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
                "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI",
                "chrM"},
}


def get_standard_chroms(genome: Optional[str]) -> set:
    """Return the set of canonical chromosome names for `genome`."""
    if genome and genome in STANDARD_CHROMS:
        return STANDARD_CHROMS[genome]
    return STANDARD_CHROMS["hg38"]


def filter_standard_chromosomes(frame: pd.DataFrame, genome: Optional[str] = None) -> pd.DataFrame:
    """Drop rows whose `chrom` is not in the standard set for `genome`."""
    standard = get_standard_chroms(genome)
    before = len(frame)
    kept = frame[frame["chrom"].isin(standard)].copy()
    removed = before - len(kept)
    if removed:
        print(f"   - removed {removed} regions on non-standard contigs")
    return kept.reset_index(drop=True)


def filter_min_region_size(frame: pd.DataFrame, min_size: int = 10) -> pd.DataFrame:
    """Drop rows with `region_size < min_size`."""
    if "region_size" not in frame.columns:
        return frame
    before = len(frame)
    kept = frame[frame["region_size"] >= min_size].copy()
    removed = before - len(kept)
    if removed:
        print(f"   - removed {removed} regions smaller than {min_size} bp")
    return kept.reset_index(drop=True)


def filter_regions(
    frame: pd.DataFrame,
    max_qvalue: Optional[float] = None,
    min_abs_logfc: Optional[float] = None,
    min_score: int = 0,
    min_region_size: int = 0,
    max_region_size: Optional[int] = None,
    direction: Optional[str] = None,
    chromosomes: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Apply a configurable cascade of filters."""
    current = frame.copy()
    starting = len(current)

    if max_qvalue is not None and "qvalue" in current.columns:
        current = current[current["qvalue"] < max_qvalue]
        print(f"   - qvalue < {max_qvalue}: {len(current)} rows")
    if min_abs_logfc is not None and "logFC" in current.columns:
        current = current[current["logFC"].abs() >= min_abs_logfc]
        print(f"   - |logFC| >= {min_abs_logfc}: {len(current)} rows")
    if min_score > 0 and "score" in current.columns:
        current = current[current["score"] >= min_score]
        print(f"   - score >= {min_score}: {len(current)} rows")
    if min_region_size > 0 and "region_size" in current.columns:
        current = current[current["region_size"] >= min_region_size]
        print(f"   - region_size >= {min_region_size}: {len(current)} rows")
    if max_region_size is not None and "region_size" in current.columns:
        current = current[current["region_size"] <= max_region_size]
        print(f"   - region_size <= {max_region_size}: {len(current)} rows")
    if direction is not None and "direction" in current.columns:
        current = current[current["direction"] == direction]
        print(f"   - direction = {direction}: {len(current)} rows")
    if chromosomes:
        chroms = list(chromosomes)
        current = current[current["chrom"].isin(chroms)]
        print(f"   - chrom in {chroms}: {len(current)} rows")

    if len(current) != starting:
        print(f"   filter: {starting} -> {len(current)} rows")
    return current.reset_index(drop=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="Parsed CSV (diff_regions_*.csv)")
    parser.add_argument("--genome", default=None)
    parser.add_argument("--max-qvalue", type=float, default=None)
    parser.add_argument("--min-abs-logfc", type=float, default=None)
    parser.add_argument("--min-region-size", type=int, default=0)
    parser.add_argument("--max-region-size", type=int, default=None)
    parser.add_argument("--direction", choices=["A_enriched", "B_enriched"], default=None)
    parser.add_argument("--out", default="filtered_regions.csv")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.csv)
    if args.genome:
        frame = filter_standard_chromosomes(frame, genome=args.genome)
        frame = filter_min_region_size(frame)
    frame = filter_regions(
        frame,
        max_qvalue=args.max_qvalue,
        min_abs_logfc=args.min_abs_logfc,
        min_region_size=args.min_region_size,
        max_region_size=args.max_region_size,
        direction=args.direction,
    )
    frame.to_csv(args.out, index=False)
    print(f"✓ Wrote {len(frame)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
