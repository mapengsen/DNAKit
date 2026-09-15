"""
Write the full output set for ChIP-Atlas peak enrichment.

Files
-----
| Path                                | Notes |
|-------------------------------------|-------|
| `analysis_object.pkl`               | Reload with `pickle.load(...)` |
| `enrichment_results_all.csv`        | Every API row |
| `enrichment_results_significant.csv`| qvalue<0.05 subset |
| `enrichment_results_top20.csv`      | Top 20 by qvalue, overlap >= 2 |
| `failed_genes.txt`                  | Ensembl-failed symbols (only when present) |
| `summary_report.md`                 | Human-readable analysis report |

Inputs : the dict returned by `run_enrichment_workflow`
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


THRESHOLD_PVAL = {50: "1e-5", 100: "1e-10", 200: "1e-20", 500: "1e-50"}


def export_all(results: dict, output_dir: str = "enrichment_results") -> None:
    os.makedirs(output_dir, exist_ok=True)
    print("\nExporting enrichment results…")

    frame: pd.DataFrame = results["enrichment_results"]
    params: dict = results.get("parameters", {})
    input_genes = list(results.get("input_genes", []))
    input_regions = list(results.get("input_regions", []))
    failed_genes = list(results.get("failed_genes", []))
    api_total_regions = int(results.get("api_total_regions", len(input_genes)))

    significant = frame[frame["q_value"] < 0.05] if len(frame) else frame
    rankable = frame[frame["regions_with_overlaps"] >= 2] if len(frame) else frame

    pickle_path = os.path.join(output_dir, "analysis_object.pkl")
    payload = {
        "enrichment_results": frame,
        "input_genes": input_genes,
        "input_regions": input_regions,
        "failed_genes": failed_genes,
        "metadata": results.get("metadata", frame),
        "parameters": params,
        "api_total_regions": api_total_regions,
        "n_genes": len(input_genes),
        "n_regions": len(input_regions),
        "n_experiments": len(frame),
        "n_significant": len(significant),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    with open(pickle_path, "wb") as handle:
        pickle.dump(payload, handle)
    print(f"   wrote {pickle_path}")

    all_path = os.path.join(output_dir, "enrichment_results_all.csv")
    frame.to_csv(all_path, index=False)
    print(f"   wrote {all_path} ({len(frame)} rows)")

    if len(significant):
        sig_path = os.path.join(output_dir, "enrichment_results_significant.csv")
        significant.to_csv(sig_path, index=False)
        print(f"   wrote {sig_path} ({len(significant)} significant rows)")
    else:
        print("   no rows with q<0.05")

    top20_source = rankable if len(rankable) else frame
    top20 = top20_source.sort_values("q_value").head(20)
    top20_path = os.path.join(output_dir, "enrichment_results_top20.csv")
    top20.to_csv(top20_path, index=False)
    print(f"   wrote {top20_path} ({len(top20)} rows)")

    if failed_genes:
        failed_path = os.path.join(output_dir, "failed_genes.txt")
        with open(failed_path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(failed_genes))
        print(f"   wrote {failed_path} ({len(failed_genes)} symbols)")

    report_path = os.path.join(output_dir, "summary_report.md")
    _write_markdown(
        report_path=report_path,
        frame=frame,
        significant=significant,
        rankable=rankable,
        params=params,
        input_genes=input_genes,
        input_regions=input_regions,
        failed_genes=failed_genes,
        api_total_regions=api_total_regions,
    )
    print(f"   wrote {report_path}")
    print("✓ Export complete")


def _write_markdown(
    report_path: str,
    frame: pd.DataFrame,
    significant: pd.DataFrame,
    rankable: pd.DataFrame,
    params: dict,
    input_genes: list,
    input_regions: list,
    failed_genes: list,
    api_total_regions: int,
) -> None:
    threshold = int(params.get("threshold", 50))
    threshold_pval = THRESHOLD_PVAL.get(threshold, "?")
    n_genes = len(input_genes)
    n_ensembl = len(input_regions)

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("# ChIP-Atlas Peak Enrichment Report\n\n")
        handle.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        handle.write("## Analysis Parameters\n\n")
        handle.write(f"- Genome: `{params.get('genome', 'hg38')}`\n")
        handle.write(f"- Antigen class: `{params.get('antigen_class', 'TFs and others')}`\n")
        handle.write(f"- Cell class: `{params.get('cell_class', 'All cell types')}`\n")
        handle.write(
            f"- Threshold: `{threshold}` (MACS2 -10·log10(p) ≥ {threshold}; ≈ p < {threshold_pval}). "
            "Default 50 balances sensitivity and specificity.\n"
        )
        handle.write(
            f"- TSS window: {params.get('distance_up', 5000)} bp upstream, "
            f"{params.get('distance_down', 5000)} bp downstream\n"
        )
        if params.get("antigen_filter"):
            handle.write(f"- Antigen filter: {params['antigen_filter']}\n")
        if params.get("cell_type_filter"):
            handle.write(f"- Cell-type filter: {params['cell_type_filter']}\n")

        if n_genes <= 10:
            handle.write(
                f"\n> ⚠ **Exploratory analysis (N={n_genes} genes).** Each gene contributes "
                f"~{100.0 / n_genes:.0f}% of the overlap signal. Lead the user-facing narrative with "
                "this caveat — every moderate fold-enrichment finding could be eliminated by a single "
                "gene's inclusion or exclusion.\n\n"
            )

        handle.write("\n## Summary Statistics\n\n")
        handle.write(f"- Input genes submitted: {n_genes}\n")
        if api_total_regions != n_genes:
            handle.write(
                f"- API regions analysed: {api_total_regions} "
                f"({n_genes - api_total_regions} gene(s) dropped — possible causes: gene not in RefSeq, "
                "unrecognised symbol, or duplicate promoter region). The API does not report which "
                "specific gene was dropped.\n"
            )
        handle.write(f"- Ensembl regions resolved: {n_ensembl}\n")
        handle.write(f"- Failed Ensembl lookups: {len(failed_genes)}\n")
        handle.write(f"- Experiments returned: {len(frame)}\n")
        handle.write(f"- Significant (q<0.05): {len(significant)}\n")

        if n_ensembl == 0 and n_genes > 0:
            handle.write(
                "\n⚠ Ensembl coordinate verification returned 0 genes. The ChIP-Atlas API uses RefSeq "
                "internally and the enrichment results are valid, but they lack independent coordinate "
                "verification. State this explicitly in the user summary; do not dismiss it as "
                "\"optional.\"\n"
            )
        elif 0 < n_ensembl < api_total_regions:
            handle.write(
                f"\n⚠ Ensembl mapped {n_ensembl}/{n_genes} genes while the API analysed "
                f"{api_total_regions}. The two systems use different gene databases (Ensembl vs "
                "RefSeq) — surface both counts in the summary and note that enrichment results rely "
                "on the API's RefSeq mapping.\n"
            )

        handle.write("\n## Top 10 Experiments (per-experiment ranking)\n\n")
        handle.write("*Individual ChIP-seq experiments. The same factor may appear multiple times.*\n\n")
        handle.write("| Rank | Factor | Cell type | Q-value | Fold enrichment | Overlap | Regions |\n")
        handle.write("|------|--------|-----------|---------|-----------------|---------|---------|\n")
        top10 = rankable.head(10) if len(rankable) else frame.head(10)
        for rank, (_, row) in enumerate(top10.iterrows(), 1):
            fe_display = f"{row['fold_enrichment']:.1f}" if row['fold_enrichment'] < 100_000 else ">100,000"
            handle.write(
                f"| {rank} | {row['antigen']} | {row['cell_type']} | "
                f"{row['q_value']:.2e} | {fe_display}× | {row['overlap_rate']:.0%} | "
                f"{int(row['regions_with_overlaps'])}/{int(row['total_regions'])} |\n"
            )
        handle.write("\n*Minimum 2 gene overlaps required for ranking.*\n")

        handle.write("\n## Top Factors (aggregated by unique factor)\n\n")
        agg_source = rankable if len(rankable) else frame
        if len(agg_source):
            sig_source = agg_source[agg_source["q_value"] < 0.05]
            n_unique_sig = int(sig_source["antigen"].nunique()) if len(sig_source) else 0
            agg = agg_source.groupby("antigen").agg(
                n_experiments=("experiment_id", "count"),
                best_qvalue=("q_value", "min"),
                median_overlap=("overlap_rate", "median"),
                max_overlap=("overlap_rate", "max"),
                median_fe=("fold_enrichment", "median"),
            ).reset_index()

            if len(sig_source):
                sig_agg = sig_source.groupby("antigen").agg(
                    n_significant=("experiment_id", "count"),
                    median_fe_sig=("fold_enrichment", "median"),
                ).reset_index()
                agg = agg.merge(sig_agg, on="antigen", how="left")
                agg["n_significant"] = agg["n_significant"].fillna(0).astype(int)
                agg["median_fe_sig"] = agg["median_fe_sig"].fillna(0)
            else:
                agg["n_significant"] = 0
                agg["median_fe_sig"] = 0.0

            agg = agg.sort_values("best_qvalue").head(10)

            handle.write(
                f"*Top 10 of {n_unique_sig} significantly enriched factor(s). "
                "`Experiments` counts experiments with ≥2 overlaps; `Sig` counts those reaching q<0.05. "
                "Quote `Median FE (Sig)` when discussing enrichment strength — `Median FE (All)` "
                "averages in non-significant experiments and dilutes the signal.*\n\n"
            )
            handle.write(
                "| Rank | Factor | Experiments | Sig | Best Q-value | Median FE (Sig) | Median FE (All) | Max overlap |\n"
            )
            handle.write(
                "|------|--------|-------------|-----|--------------|-----------------|-----------------|-------------|\n"
            )
            for rank, (_, row) in enumerate(agg.iterrows(), 1):
                fe_all = f"{row['median_fe']:.1f}" if row['median_fe'] < 100_000 else ">100,000"
                if row["n_significant"] > 0:
                    fe_sig = f"{row['median_fe_sig']:.1f}" if row['median_fe_sig'] < 100_000 else ">100,000"
                    fe_sig_cell = f"{fe_sig}×"
                else:
                    fe_sig_cell = "—"
                handle.write(
                    f"| {rank} | {row['antigen']} | {int(row['n_experiments'])} | "
                    f"{int(row['n_significant'])} | {row['best_qvalue']:.2e} | {fe_sig_cell} | "
                    f"{fe_all}× | {row['max_overlap']:.0%} |\n"
                )

            handle.write(
                "\n⚠ Per-factor multiple testing: `Best Q-value` is the BH-corrected q from each "
                "factor's most significant experiment. Interpreting all top-N factors as independent "
                "discoveries overstates confidence; q-values are not corrected across factors.\n"
            )
            if n_unique_sig > 10:
                handle.write(
                    f"\n⚠ {n_unique_sig - 10} additional significantly enriched factor(s) sit "
                    f"below the top 10 shown here. Report the total ({n_unique_sig}) when "
                    "summarizing, not just the visible rows.\n"
                )
            if n_genes <= 10:
                handle.write(
                    f"\n⚠ Small gene set ({n_genes}): every moderate fold-enrichment factor in this "
                    "table needs an explicit caveat that removing a single gene could eliminate the "
                    "signal.\n"
                )

        if len(significant):
            handle.write("\n## Overlap Rate Summary\n\n")
            handle.write(
                f"Across {len(significant)} significant experiments (q<0.05): median "
                f"{significant['overlap_rate'].median():.0%}, range "
                f"{significant['overlap_rate'].min():.0%}–{significant['overlap_rate'].max():.0%}, "
                f"mean {significant['overlap_rate'].mean():.0%}. "
                "Use these numbers verbatim in the user summary.\n"
            )
            if len(rankable):
                top10_data = rankable.head(10)
                handle.write(
                    f"\nAcross the top 10 experiments (overlap ≥ 2): median "
                    f"{top10_data['overlap_rate'].median():.0%}, range "
                    f"{top10_data['overlap_rate'].min():.0%}–{top10_data['overlap_rate'].max():.0%}.\n"
                )

        handle.write("\n## Interpretation Reference\n\n")
        handle.write("- q < 0.001 — highly significant after correction\n")
        handle.write("- q < 0.01  — significant\n")
        handle.write("- q < 0.05  — genome-wide significant\n")
        handle.write("- q ≥ 0.05  — not significant\n\n")
        handle.write("- Fold enrichment > 10× — very strong; 5–10× strong; 2–5× moderate; <2× weak\n")
        handle.write("- Fold enrichment > 100,000× — sentinel (zero background overlap); ignore\n\n")
        handle.write(
            f"- Threshold {threshold} corresponds to MACS2 peak calling at ~p < {threshold_pval}. "
            "Higher thresholds yield fewer but more confident peaks.\n"
        )

        handle.write("\n## Important Caveats\n\n")
        handle.write(
            "1. **Data availability bias.** Well-studied factors (TP53, RELA, CTCF) have hundreds "
            "to thousands of experiments and naturally show up in any broad enrichment. Compare "
            "`Experiments` vs `Sig` columns when reporting a factor.\n"
        )
        handle.write(
            "2. **Redundant experiments.** The per-experiment table contains the same factor in "
            "multiple cell types and labs. Use the aggregated factor table for deduplicated counts.\n"
        )
        handle.write(
            "3. **Validation required.** Validate top hits with orthogonal methods (expression "
            "correlation, motif analysis, perturbation experiments).\n"
        )
        handle.write(
            "4. **Threshold sensitivity.** Rerun at `threshold=100` to confirm robustness; if the "
            "top factors survive, the call is stable.\n"
        )
        if n_genes <= 10:
            pct_per_gene = 100.0 / max(api_total_regions, 1)
            handle.write(
                f"5. **Small gene set ({api_total_regions} regions).** Each gene contributes "
                f"~{pct_per_gene:.1f}% to overlap rates. Lead with exploratory framing and validate "
                "any finding with a larger gene set (20–100 genes recommended).\n"
            )

        handle.write("\n## Files Generated\n\n")
        handle.write("- `analysis_object.pkl` — full results object\n")
        handle.write("- `enrichment_results_all.csv` — every API row\n")
        if len(significant):
            handle.write("- `enrichment_results_significant.csv` — q<0.05 subset\n")
        handle.write("- `enrichment_results_top20.csv` — top 20 by q-value, overlap ≥ 2\n")
        if failed_genes:
            handle.write(f"- `failed_genes.txt` — {len(failed_genes)} unresolved Ensembl symbols\n")
        handle.write("- `chipatlas_enrichment.{png,svg}` — 4-panel summary figure\n")

        handle.write("\n## Citations\n\n")
        handle.write(
            "- Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore chromosome "
            "architecture together with large-scale regulome data. *Nucleic Acids Research* "
            "52(D1):D45–D53. [doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)\n"
        )
        handle.write(
            "- Oki S. *et al.* (2018). ChIP-Atlas: a data-mining suite powered by full integration "
            "of public ChIP-seq data. *EMBO Reports* 19(12):e46255. "
            "[doi:10.15252/embr.201846255](https://doi.org/10.15252/embr.201846255)\n"
        )


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
