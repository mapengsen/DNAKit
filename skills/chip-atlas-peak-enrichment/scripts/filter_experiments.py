"""
Post-processing filters for ChIP-Atlas enrichment results.

Useful when the broad `antigenClass`/`cellClass` categories return too
many rows or you want to focus on a particular factor / cell type / FE
range without resubmitting the API job.

Inputs : a DataFrame from `parse_api_results`
Outputs: a filtered DataFrame

Example
-------
    python filter_experiments.py --csv enrichment_results_all.csv \
        --antigen TP53 MYC --min-fold-enrichment 5 --max-qvalue 0.05
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


def filter_experiments(
    frame: pd.DataFrame,
    antigen_filter: Optional[Iterable[str]] = None,
    cell_type_filter: Optional[Iterable[str]] = None,
    cell_class_filter: Optional[Iterable[str]] = None,
    min_fold_enrichment: float = 1.0,
    max_qvalue: Optional[float] = None,
    max_pvalue: Optional[float] = None,
    min_overlap: int = 0,
) -> pd.DataFrame:
    current = frame.copy()
    starting = len(current)

    if antigen_filter:
        wanted = {value.lower() for value in antigen_filter}
        current = current[current["antigen"].str.lower().isin(wanted)]
        print(f"   - antigen ∈ {sorted(wanted)}: {len(current)} rows")

    if cell_type_filter:
        wanted = {value.lower() for value in cell_type_filter}
        current = current[current["cell_type"].str.lower().isin(wanted)]
        print(f"   - cell_type ∈ {sorted(wanted)}: {len(current)} rows")

    if cell_class_filter and "cell_type_class" in current.columns:
        wanted = {value.lower() for value in cell_class_filter}
        current = current[current["cell_type_class"].str.lower().isin(wanted)]
        print(f"   - cell_type_class ∈ {sorted(wanted)}: {len(current)} rows")

    if min_fold_enrichment > 1.0:
        current = current[current["fold_enrichment"] >= min_fold_enrichment]
        print(f"   - fold_enrichment >= {min_fold_enrichment}: {len(current)} rows")

    if max_qvalue is not None and "q_value" in current.columns:
        current = current[current["q_value"] <= max_qvalue]
        print(f"   - q_value <= {max_qvalue}: {len(current)} rows")

    if max_pvalue is not None and "p_value" in current.columns:
        current = current[current["p_value"] <= max_pvalue]
        print(f"   - p_value <= {max_pvalue}: {len(current)} rows")

    if min_overlap > 0:
        current = current[current["regions_with_overlaps"] >= min_overlap]
        print(f"   - regions_with_overlaps >= {min_overlap}: {len(current)} rows")

    if len(current) != starting:
        print(f"   filter: {starting} -> {len(current)} rows")
    return current.reset_index(drop=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--antigen", nargs="*", default=None)
    parser.add_argument("--cell-type", nargs="*", default=None)
    parser.add_argument("--cell-class", nargs="*", default=None)
    parser.add_argument("--min-fold-enrichment", type=float, default=1.0)
    parser.add_argument("--max-qvalue", type=float, default=None)
    parser.add_argument("--max-pvalue", type=float, default=None)
    parser.add_argument("--min-overlap", type=int, default=0)
    parser.add_argument("--out", default="filtered_enrichment.csv")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.csv)
    filtered = filter_experiments(
        frame,
        antigen_filter=args.antigen,
        cell_type_filter=args.cell_type,
        cell_class_filter=args.cell_class,
        min_fold_enrichment=args.min_fold_enrichment,
        max_qvalue=args.max_qvalue,
        max_pvalue=args.max_pvalue,
        min_overlap=args.min_overlap,
    )
    filtered.to_csv(args.out, index=False)
    print(f"✓ Wrote {len(filtered)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
