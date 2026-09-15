"""
Write the full target-genes output set to disk.

Files
-----
| Path                              | Notes |
|-----------------------------------|-------|
| `analysis_object.pkl`             | Reload with `pickle.load(...)` |
| `target_genes_all.csv`            | All retained target genes |
| `target_genes_top50.csv`          | Top 50 by `avg_score` |
| `target_genes_with_string.csv`    | Genes with `string_score > 0` (only when present) |
| `experiment_scores_top50.csv`     | Wide-format per-experiment scores for the top 50 |
| `summary_report.md`               | Human-readable analysis report |

Inputs : the dict returned by `run_target_genes_workflow`
Outputs: filesystem only
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from datetime import datetime
from typing import Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas")


def export_all(results: dict, output_dir: str = "target_genes_results") -> None:
    os.makedirs(output_dir, exist_ok=True)
    target_genes: pd.DataFrame = results["target_genes"]
    experiment_data: pd.DataFrame = results["experiment_data"]
    protein: str = results["protein"]
    parameters: dict = results["parameters"]
    metadata: dict = results["metadata"]

    print(f"\nExporting target-genes results to {output_dir}/")

    pickle_path = os.path.join(output_dir, "analysis_object.pkl")
    payload = {
        "target_genes": target_genes,
        "experiment_data": experiment_data,
        "cell_types": results.get("cell_types", []),
        "protein": protein,
        "parameters": parameters,
        "metadata": metadata,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    with open(pickle_path, "wb") as handle:
        pickle.dump(payload, handle)
    print(f"   wrote {pickle_path}")

    all_path = os.path.join(output_dir, "target_genes_all.csv")
    target_genes.to_csv(all_path, index=False)
    print(f"   wrote {all_path} ({len(target_genes)} rows)")

    top50 = target_genes.head(50)
    top50_path = os.path.join(output_dir, "target_genes_top50.csv")
    top50.to_csv(top50_path, index=False)
    print(f"   wrote {top50_path}")

    with_string = target_genes[target_genes["string_score"] > 0]
    if len(with_string):
        string_path = os.path.join(output_dir, "target_genes_with_string.csv")
        with_string.to_csv(string_path, index=False)
        print(f"   wrote {string_path} ({len(with_string)} genes)")

    if experiment_data is not None and len(experiment_data):
        top50_genes = top50["gene"].tolist()
        wide = experiment_data[experiment_data["gene"].isin(top50_genes)]
        wide_path = os.path.join(output_dir, "experiment_scores_top50.csv")
        wide.to_csv(wide_path, index=False)
        print(f"   wrote {wide_path} ({len(wide)} rows)")

    report_path = os.path.join(output_dir, "summary_report.md")
    _write_markdown(report_path, results)
    print(f"   wrote {report_path}")
    print("✓ Export complete")


def _write_markdown(report_path: str, results: dict) -> None:
    target_genes: pd.DataFrame = results["target_genes"]
    protein: str = results["protein"]
    metadata: dict = results["metadata"]
    parameters: dict = results["parameters"]
    experiment_data: Optional[pd.DataFrame] = results.get("experiment_data")
    top10 = target_genes.head(10)
    avg_scores = target_genes["avg_score"]

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(f"# {protein} Target Genes Report\n\n")
        handle.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        handle.write(
            "**Data source:** ChIP-Atlas — Zou Z. *et al.* (2024). "
            "ChIP-Atlas 3.0. *Nucleic Acids Research* 52(D1):D45–D53. "
            "[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)\n\n"
        )

        handle.write("## Analysis Parameters\n\n")
        handle.write("| Parameter | Value |\n")
        handle.write("|-----------|-------|\n")
        handle.write(f"| Protein | {protein} |\n")
        handle.write(f"| Genome | {metadata['genome']} |\n")
        handle.write(f"| Distance from TSS | ±{metadata['distance_kb']} kb |\n")
        handle.write(f"| Total genes (before filter) | {metadata['total_genes_before_filter']} |\n")
        handle.write(f"| Total genes (after filter) | {metadata['total_genes_after_filter']} |\n")
        colocated = metadata.get("colocated_info", {})
        if colocated.get("n_groups", 0):
            handle.write(
                f"| Independent loci | {colocated['n_independent_loci']} "
                f"({colocated['n_genes_affected']} genes in "
                f"{colocated['n_groups']} co-located groups) |\n"
            )
        handle.write(f"| Total experiments | {metadata['total_experiments']} |\n")
        handle.write(f"| Unique cell types | {metadata['unique_cell_types']} |\n")
        cell_filter = parameters.get("cell_types")
        if cell_filter:
            handle.write(f"| Cell-type filter | {cell_filter} |\n")

        handle.write("\n## Summary Statistics\n\n")
        handle.write(f"- Median binding score: {avg_scores.median():.1f}\n")
        handle.write(f"- Mean binding score: {avg_scores.mean():.1f}\n")
        handle.write(f"- Max binding score: {avg_scores.max():.1f}\n")
        n_500 = int((avg_scores >= 500).sum())
        n_200 = int((avg_scores >= 200).sum())
        n_100 = int((avg_scores >= 100).sum())
        n_50 = int((avg_scores >= 50).sum())
        handle.write(f"- Genes with score ≥500 (per-experiment Q ≤ 1e-50): {n_500}\n")
        handle.write(f"- Genes with score ≥200 (per-experiment Q ≤ 1e-20): {n_200}\n")
        handle.write(f"- Genes with score ≥100 (per-experiment Q ≤ 1e-10): {n_100}\n")
        handle.write(f"- Genes with score ≥50 (per-experiment Q ≤ 1e-5): {n_50}\n")
        n_below_50 = len(avg_scores) - n_50
        if n_below_50:
            handle.write(
                f"- Genes with score <50 (weak evidence): {n_below_50} "
                f"({100 * n_below_50 / len(avg_scores):.0f}% of returned genes)\n"
            )
        top10_rate = top10["binding_rate"].median() if len(top10) else 0.0
        handle.write(
            f"\n> Average scores are diluted by zero-binding experiments. The top 10 genes are "
            f"bound in a median of {top10_rate:.1%} of experiments; the remaining "
            f"~{100 - top10_rate * 100:.1f}% contribute zeros that pull the average down. "
            "Use `binding_rate` and `max_score` for complementary views.\n\n"
        )

        n_with_string = int((target_genes["string_score"] > 0).sum())
        if n_with_string:
            pct_string = 100 * n_with_string / len(target_genes)
            min_nonzero = float(target_genes.loc[target_genes["string_score"] > 0, "string_score"].min())
            handle.write(
                f"- Genes with STRING interactions: {n_with_string}/{len(target_genes)} "
                f"({pct_string:.0f}%) — STRING scores are pre-embedded in the ChIP-Atlas data "
                "(not queried separately). Coverage is sparse "
                f"(lowest non-zero score: {min_nonzero:.0f}). A STRING score of 0 does NOT mean "
                "'not a target.'\n"
            )
        else:
            handle.write(
                f"- Genes with STRING interactions: 0/{len(target_genes)} — no STRING data "
                f"available for {protein} in ChIP-Atlas.\n"
            )

        before = metadata.get("total_genes_before_filter", len(target_genes))
        after = metadata.get("total_genes_after_filter", len(target_genes))
        if after < before and len(avg_scores):
            handle.write(
                f"- Filtering cutoff: kept top {after} genes "
                f"(lowest included avg_score: {avg_scores.iloc[-1]:.1f})\n"
            )

        handle.write("\n## Top 10 Target Genes\n\n")
        handle.write(
            "> Rankings are data-derived. Use these exact genes, ranks, and scores when summarizing "
            "results — do not substitute well-known targets that the data places lower.\n\n"
        )
        total_experiments = metadata.get("total_experiments", 0)
        handle.write("| Rank | Gene | Avg score | Max score | Binding rate (n/N) | STRING |\n")
        handle.write("|------|------|-----------|-----------|--------------------|--------|\n")
        for rank, (_, row) in enumerate(top10.iterrows(), 1):
            gene = row["gene"]
            if row.get("colocated_group", 0) > 0:
                gene = f"{gene} *"
            handle.write(
                f"| {rank} | {gene} | {row['avg_score']:.1f} | {row['max_score']:.0f} | "
                f"{row['binding_rate']:.1%} ({int(row['num_bound'])}/{total_experiments}) | "
                f"{row['string_score']:.0f} |\n"
            )

        top20 = target_genes.head(20)
        if len(top20) > 10:
            handle.write("\n### Ranks 11–20\n\n")
            handle.write("| Rank | Gene | Avg score | Max score | Binding rate (n/N) | STRING |\n")
            handle.write("|------|------|-----------|-----------|--------------------|--------|\n")
            for rank, (_, row) in enumerate(top20.iloc[10:].iterrows(), 11):
                gene = row["gene"]
                if row.get("colocated_group", 0) > 0:
                    gene = f"{gene} *"
                handle.write(
                    f"| {rank} | {gene} | {row['avg_score']:.1f} | {row['max_score']:.0f} | "
                    f"{row['binding_rate']:.1%} ({int(row['num_bound'])}/{total_experiments}) | "
                    f"{row['string_score']:.0f} |\n"
                )
            if "colocated_group" in top20.columns and (top20["colocated_group"] > 0).any():
                handle.write("\n\\* Gene is part of a co-located group (see Caveats).\n")

        if experiment_data is not None:
            exp_cols = [c for c in experiment_data.columns if c != "gene"]
            if exp_cols:
                cell_type_counts: dict = {}
                for column in exp_cols:
                    parts = column.split("|", 1)
                    name = parts[1] if len(parts) == 2 else column
                    cell_type_counts[name] = cell_type_counts.get(name, 0) + 1
                ordered = sorted(cell_type_counts.items(), key=lambda item: -item[1])
                total_exp = len(exp_cols)
                handle.write("\n## Experiment Composition\n\n")
                handle.write(f"Total experiments: {total_exp}\n\n")
                handle.write("Top cell types by experiment count:\n\n")
                handle.write("| Cell type | Experiments | Share |\n")
                handle.write("|-----------|------------|-------|\n")
                top5_sum = 0
                for name, count in ordered[:5]:
                    pct = 100 * count / total_exp
                    top5_sum += count
                    handle.write(f"| {name} | {count} | {pct:.1f}% |\n")
                top5_pct = 100 * top5_sum / total_exp
                n_rare = sum(1 for _, count in ordered if 100 * count / total_exp < 1)
                n_dominant = sum(1 for _, count in ordered if 100 * count / total_exp >= 5)
                handle.write(
                    f"\nUneven distribution: the top 5 cell types account for {top5_pct:.0f}% of "
                    f"experiments; {n_rare} of {len(ordered)} cell types contribute <1% each.\n"
                )
                handle.write(f"- {n_dominant} cell types contribute ≥5% each.\n")
                handle.write(
                    "- To remove this bias, re-run with `cell_types=[...]` matching your "
                    "system of interest.\n"
                )

        handle.write("\n## Interpretation Reference\n\n")
        handle.write("Binding score (MACS2): −10 × log10(Q-value). Higher = stronger evidence.\n")
        handle.write("- Score ≥500: very strong (per-experiment Q ≤ 1e-50)\n")
        handle.write("- Score 100–500: strong\n")
        handle.write("- Score 50–100: moderate\n")
        handle.write("- Score <50: weak\n\n")
        handle.write(
            "Note that those thresholds apply to **individual experiment** scores. The average "
            "in the ranking includes zeros from non-binding experiments, so an average of 500 "
            "means the *consensus* reaches that level, not that every experiment hit Q ≤ 1e-50.\n\n"
        )
        handle.write(
            "Binding rate: fraction of experiments with any binding. High rate = consistent "
            "across cell types; low rate + high `max_score` = strong but cell-type-specific.\n\n"
        )
        handle.write(
            "STRING score: 0–1000. Genes with both high binding and high STRING are the "
            "highest-confidence targets, but STRING=0 does NOT mean a gene is not a target — "
            "even canonical targets (BBC3/PUMA for TP53) can have STRING score 0.\n"
        )

        handle.write("\n## Caveats\n\n")
        handle.write(
            "- **Averaging includes zeros.** Genes bound in fewer experiments get diluted "
            "averages. Use `max_score` and `binding_rate` columns for complementary views.\n"
        )
        handle.write(
            "- **Cell-type bias.** Experiment counts are uneven (see the Experiment Composition "
            "table). Rankings skew toward well-studied lines; re-run with `cell_types=[...]` "
            "for a system-specific view.\n"
        )
        handle.write(
            "- **External annotations.** ChIP-Atlas only provides binding data. Any biological "
            "function description in the user-facing summary comes from background knowledge, "
            "not from this analysis.\n"
        )
        if colocated.get("n_groups", 0):
            n_affected = colocated["n_genes_affected"]
            n_groups = colocated["n_groups"]
            n_independent = colocated["n_independent_loci"]
            total = metadata["total_genes_after_filter"]
            pct = 100 * n_affected / total if total else 0.0
            handle.write(
                f"- **Co-located genes.** {n_affected} of {total} genes ({pct:.1f}%) share "
                f"identical binding scores with at least one other gene because they sit at "
                f"the same locus within the ±{metadata['distance_kb']}kb TSS window. The "
                f"effective number of independent loci is {n_independent} "
                f"(not {total}). For pathway enrichment, consider collapsing co-located "
                f"groups (the `colocated_group` column marks them) to avoid double-counting "
                f"loci.\n"
            )
            top_groups = sorted(colocated["groups"], key=lambda g: g["size"], reverse=True)[:3]
            if top_groups:
                examples = []
                for group in top_groups:
                    names = ", ".join(group["genes"][:4])
                    if group["size"] > 4:
                        names += f", … ({group['size']} total)"
                    examples.append(names)
                handle.write(f"  - Largest groups: {'; '.join(examples)}\n")

        handle.write("\n## Files Generated\n\n")
        handle.write("| File | Description |\n")
        handle.write("|------|-------------|\n")
        handle.write("| `analysis_object.pkl` | Full results object |\n")
        handle.write("| `target_genes_all.csv` | All retained target genes |\n")
        handle.write("| `target_genes_top50.csv` | Top 50 by `avg_score` |\n")
        if (target_genes["string_score"] > 0).any():
            handle.write("| `target_genes_with_string.csv` | Genes with STRING score > 0 |\n")
        handle.write("| `experiment_scores_top50.csv` | Per-experiment scores for the top 50 |\n")
        handle.write("| `target_genes_top_targets.{png,svg}` | Top targets bar plot |\n")
        handle.write("| `target_genes_score_distribution.{png,svg}` | Score distribution |\n")
        handle.write("| `target_genes_heatmap.{png,svg}` | Clustered heatmap |\n")
        handle.write("| `target_genes_string_vs_binding.{png,svg}` | STRING vs binding scatter |\n")
        handle.write("| `summary_report.md` | This report |\n")


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pickle", required=True)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args(argv)
    with open(args.pickle, "rb") as handle:
        results = pickle.load(handle)
    output_dir = args.output_dir or os.path.dirname(os.path.abspath(args.pickle))
    export_all(results, output_dir=output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
