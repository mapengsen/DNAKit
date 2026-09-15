"""
Post-download filters for ChIP-Atlas target-genes results.

Two functions are exported:

* `filter_targets`         — score / cell-type / STRING / binding-rate /
                             top-N filtering. When `cell_types` is set,
                             `avg_score`, `max_score`, `binding_rate`,
                             and `num_bound` are **recomputed** from
                             only the matching experiment columns.

* `detect_colocated_genes` — flags genes that share identical
                             per-experiment score vectors (i.e. genes
                             that sit at the same locus within the TSS
                             window). Adds a `colocated_group` column.

Inputs : summary DataFrame, optional wide experiment DataFrame
Outputs: filtered DataFrames (+ info dict for the colocation pass)

Example
-------
    python filter_targets.py --summary target_genes_all.csv \
        --experiments experiment_scores_top50.csv --top-n 100
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


def filter_targets(
    summary_df: pd.DataFrame,
    experiment_df: Optional[pd.DataFrame] = None,
    min_avg_score: float = 0.0,
    cell_types: Optional[Iterable[str]] = None,
    min_string_score: float = 0.0,
    top_n: Optional[int] = None,
    min_binding_rate: float = 0.0,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    starting = len(summary_df)
    current = summary_df.copy()
    experiments = experiment_df.copy() if experiment_df is not None else None

    if cell_types and experiments is not None:
        current, experiments = _filter_by_cell_type(current, experiments, list(cell_types))

    if min_avg_score > 0:
        current = current[current["avg_score"] >= min_avg_score]
        print(f"   - avg_score >= {min_avg_score}: {len(current)}")

    if min_string_score > 0:
        current = current[current["string_score"] >= min_string_score]
        print(f"   - string_score >= {min_string_score}: {len(current)}")

    if min_binding_rate > 0:
        current = current[current["binding_rate"] >= min_binding_rate]
        print(f"   - binding_rate >= {min_binding_rate}: {len(current)}")

    if top_n is not None and top_n > 0:
        current = current.head(top_n)
        print(f"   - top_n = {top_n}: {len(current)}")

    if experiments is not None:
        experiments = experiments[experiments["gene"].isin(current["gene"])]
        experiments = experiments.reset_index(drop=True)

    current = current.reset_index(drop=True)
    if len(current) != starting:
        print(f"   filter: {starting} -> {len(current)} target genes")
    return current, experiments


def _filter_by_cell_type(
    summary_df: pd.DataFrame,
    experiment_df: pd.DataFrame,
    cell_types: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    desired = {value.lower() for value in cell_types}
    keep_columns = ["gene"]
    for column in experiment_df.columns:
        if column == "gene":
            continue
        parts = column.split("|", 1)
        if len(parts) == 2 and parts[1].lower() in desired:
            keep_columns.append(column)

    matched = len(keep_columns) - 1
    if matched == 0:
        print(f"   ! no experiments match cell types {cell_types}; keeping originals")
        sample = _sample_cell_types(experiment_df)
        print(f"   available cell types (sample): {sample}")
        return summary_df, experiment_df

    print(f"   cell_type filter: {matched} experiments retained")
    experiments = experiment_df[keep_columns].copy()
    exp_cols = [c for c in keep_columns if c != "gene"]
    values = experiments[exp_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    summary_df = summary_df.copy()
    summary_df["avg_score"] = values.mean(axis=1).values
    summary_df["num_experiments"] = matched
    summary_df["num_bound"] = (values > 0).sum(axis=1).values
    summary_df["binding_rate"] = summary_df["num_bound"] / max(matched, 1)
    summary_df["max_score"] = values.max(axis=1).values
    summary_df = summary_df.sort_values("avg_score", ascending=False).reset_index(drop=True)
    return summary_df, experiments


def detect_colocated_genes(
    summary_df: pd.DataFrame,
    experiment_df: Optional[pd.DataFrame],
) -> Tuple[pd.DataFrame, Dict]:
    """Tag genes that share identical per-experiment score vectors."""
    annotated = summary_df.copy()
    info: Dict = {"n_groups": 0, "n_genes_affected": 0, "n_independent_loci": len(annotated), "groups": []}

    if experiment_df is None or len(experiment_df) == 0:
        annotated["colocated_group"] = 0
        return annotated, info

    exp_cols = [c for c in experiment_df.columns if c != "gene"]
    if not exp_cols:
        annotated["colocated_group"] = 0
        return annotated, info

    score_matrix = experiment_df.set_index("gene")[exp_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    common = annotated["gene"][annotated["gene"].isin(score_matrix.index)].tolist()
    score_matrix = score_matrix.loc[common]

    buckets: Dict[Tuple, List[str]] = {}
    for gene in score_matrix.index:
        key = tuple(score_matrix.loc[gene].values)
        buckets.setdefault(key, []).append(gene)

    gene_to_group: Dict[str, int] = {}
    next_group = 1
    groups: List[Dict] = []
    for genes in buckets.values():
        if len(genes) < 2:
            continue
        for gene in genes:
            gene_to_group[gene] = next_group
        groups.append({"group_id": next_group, "genes": genes, "size": len(genes)})
        next_group += 1

    annotated["colocated_group"] = annotated["gene"].map(gene_to_group).fillna(0).astype(int)
    info["n_groups"] = len(groups)
    info["n_genes_affected"] = sum(group["size"] for group in groups)
    info["n_independent_loci"] = len(annotated) - info["n_genes_affected"] + info["n_groups"]
    info["groups"] = groups

    if groups:
        print(f"   colocated genes: {info['n_genes_affected']} genes in {len(groups)} groups "
              f"({info['n_independent_loci']} independent loci)")
    else:
        print("   colocated genes: none detected")
    return annotated, info


def _sample_cell_types(experiment_df: pd.DataFrame, limit: int = 10) -> List[str]:
    seen = set()
    for column in experiment_df.columns:
        if column == "gene":
            continue
        parts = column.split("|", 1)
        if len(parts) == 2:
            seen.add(parts[1])
    return sorted(seen)[:limit]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary", required=True, help="target_genes_all.csv")
    parser.add_argument("--experiments", help="experiment_scores_top50.csv (wide format)")
    parser.add_argument("--min-avg-score", type=float, default=0.0)
    parser.add_argument("--min-string-score", type=float, default=0.0)
    parser.add_argument("--min-binding-rate", type=float, default=0.0)
    parser.add_argument("--cell-types", nargs="*", default=None)
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--out", default="filtered_targets.csv")
    args = parser.parse_args(argv)

    summary = pd.read_csv(args.summary)
    experiments = pd.read_csv(args.experiments) if args.experiments else None
    filtered_summary, _ = filter_targets(
        summary,
        experiment_df=experiments,
        min_avg_score=args.min_avg_score,
        cell_types=args.cell_types,
        min_string_score=args.min_string_score,
        top_n=args.top_n,
        min_binding_rate=args.min_binding_rate,
    )
    filtered_summary.to_csv(args.out, index=False)
    print(f"✓ Wrote {len(filtered_summary)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
