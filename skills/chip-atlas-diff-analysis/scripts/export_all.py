"""
Write the full set of ChIP-Atlas diff-analysis artefacts to disk.

Files written into `output_dir`:

| File                             | Notes |
|----------------------------------|-------|
| `analysis_object.pkl`            | Reload everything with `pickle.load(...)` |
| `diff_regions_all.csv`           | Default QC-filtered table |
| `diff_regions_significant.csv`   | qvalue < 0.05 subset |
| `diff_regions_top50.csv`         | Top 50 by Q; carries `nearest_gene` when UCSC annotation succeeds |
| `diff_regions_unfiltered.csv`    | Pre-QC table (only when QC removed rows) |
| `summary_report.md`              | Human-readable analysis report |

Inputs : the dict returned by `run_diff_workflow`
Outputs: filesystem-only; the function returns nothing useful.

Example
-------
    python export_all.py --pickle diff_results/analysis_object.pkl
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


_HAS_UCSC_ANNOTATOR = True
try:
    try:
        from scripts.annotate_genes import annotate_nearest_genes
    except ImportError:
        from annotate_genes import annotate_nearest_genes
except ImportError:
    _HAS_UCSC_ANNOTATOR = False


def export_all(
    results: dict,
    output_dir: str = "diff_results",
    qvalue_threshold: float = 0.05,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    print("\nExporting results…")

    frame: pd.DataFrame = results["diff_regions"]
    unfiltered: pd.DataFrame = results.get("diff_regions_unfiltered", frame)
    params: dict = results.get("parameters", {})
    warnings = results.get("qc_warnings", [])

    pickle_path = _write_pickle(results, frame, unfiltered, output_dir)
    print(f"   wrote {pickle_path}")

    all_path = os.path.join(output_dir, "diff_regions_all.csv")
    frame.to_csv(all_path, index=False)
    print(f"   wrote {all_path} ({len(frame)} rows)")

    if "qvalue" in frame.columns:
        sig = frame[frame["qvalue"] < qvalue_threshold]
        if len(sig):
            sig_path = os.path.join(output_dir, "diff_regions_significant.csv")
            sig.to_csv(sig_path, index=False)
            print(f"   wrote {sig_path} ({len(sig)} rows, qvalue<{qvalue_threshold})")
        else:
            print(f"   no rows pass qvalue<{qvalue_threshold}")

    top50 = _select_top(frame, n=50)
    annotations = _annotate(top50, params)
    if annotations is not None:
        top50 = top50.copy()
        top50["nearest_gene"] = annotations.reindex(top50.index).fillna("")
    top50_path = os.path.join(output_dir, "diff_regions_top50.csv")
    top50.to_csv(top50_path, index=False)
    print(f"   wrote {top50_path}")

    if len(unfiltered) != len(frame):
        unfiltered_path = os.path.join(output_dir, "diff_regions_unfiltered.csv")
        unfiltered.to_csv(unfiltered_path, index=False)
        print(f"   wrote {unfiltered_path} (pre-QC, {len(unfiltered)} rows)")

    report_path = os.path.join(output_dir, "summary_report.md")
    _write_markdown(
        report_path=report_path,
        frame=frame,
        unfiltered=unfiltered,
        params=params,
        warnings=warnings,
        experiments_a=results.get("experiments_a", []),
        experiments_b=results.get("experiments_b", []),
        top10_annotations=annotations,
    )
    print(f"   wrote {report_path}")
    print("✓ Export complete")


def _write_pickle(results, frame, unfiltered, output_dir):
    payload = {
        "diff_regions": frame,
        "diff_regions_unfiltered": unfiltered,
        "qc_warnings": results.get("qc_warnings", []),
        "experiments_a": results.get("experiments_a", []),
        "experiments_b": results.get("experiments_b", []),
        "raw_files": results.get("raw_files", {}),
        "log_content": results.get("log_content", ""),
        "parameters": results.get("parameters", {}),
        "n_regions": len(frame),
        "n_regions_unfiltered": len(unfiltered),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    path = os.path.join(output_dir, "analysis_object.pkl")
    with open(path, "wb") as handle:
        pickle.dump(payload, handle)
    return path


def _select_top(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    if frame.empty:
        return frame
    if "qvalue" in frame.columns and "logFC" in frame.columns:
        ranked = frame.copy()
        ranked["_abs_logFC"] = ranked["logFC"].abs()
        return (
            ranked.sort_values(
                ["qvalue", "_abs_logFC", "chrom", "chromStart"],
                ascending=[True, False, True, True],
            )
            .head(n)
            .drop(columns=["_abs_logFC"])
        )
    if "score" in frame.columns:
        return frame.sort_values("score", ascending=False).head(n)
    return frame.head(n)


def _annotate(top: pd.DataFrame, params: dict) -> Optional[pd.Series]:
    if not _HAS_UCSC_ANNOTATOR or top.empty:
        return None
    genome = params.get("genome", "hg38")
    try:
        return annotate_nearest_genes(top, genome=genome, max_rows=len(top))
    except Exception as exc:
        print(f"   ! gene annotation failed ({exc}); continuing without it")
        return None


def _write_markdown(
    report_path: str,
    frame: pd.DataFrame,
    unfiltered: pd.DataFrame,
    params: dict,
    warnings: list,
    experiments_a: list,
    experiments_b: list,
    top10_annotations: Optional[pd.Series],
) -> None:
    desc_a = params.get("description_a", "Group A")
    desc_b = params.get("description_b", "Group B")
    analysis_type = params.get("analysis_type", "diffbind")
    type_label = "Differential Peak Regions" if analysis_type == "diffbind" else "Differentially Methylated Regions"

    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(f"# ChIP-Atlas Diff Analysis Report\n\n")
        handle.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        handle.write("## Analysis Parameters\n\n")
        handle.write(f"- **Analysis type:** {type_label}\n")
        handle.write(f"- **Genome:** {params.get('genome', 'N/A')}\n")
        handle.write(f"- **Group A ({desc_a}):** {', '.join(experiments_a)} (n={len(experiments_a)})\n")
        handle.write(f"- **Group B ({desc_b}):** {', '.join(experiments_b)} (n={len(experiments_b)})\n")

        design_caveats = params.get("design_caveats") or []
        if design_caveats:
            handle.write("\n## Experimental Design Caveats\n\n")
            for caveat in design_caveats:
                handle.write(f"- {caveat}\n")

        if warnings:
            handle.write("\n## QC Warnings\n\n")
            for warning in warnings:
                tag = "**MAJOR**" if warning["severity"] == "major" else "Minor"
                handle.write(f"### {tag}: {warning['issue']}\n\n")
                handle.write(f"{warning['details']}\n\n")
                handle.write(f"**Recommendation:** {warning['recommendation']}\n\n")

        if len(unfiltered) != len(frame):
            handle.write("\n## QC Filtering Applied\n\n")
            handle.write(f"- Pre-QC rows: {len(unfiltered)}\n")
            handle.write(f"- Post-QC rows: {len(frame)}\n")
            handle.write(f"- Removed: {len(unfiltered) - len(frame)} (non-standard contigs and <10 bp regions)\n")
            handle.write(f"- Pre-QC data preserved in `diff_regions_unfiltered.csv`\n")

        n_a_exp = len(experiments_a)
        n_b_exp = len(experiments_b)
        handle.write("\n## Statistical Caveats\n\n")
        handle.write(f"- Sample size: n={n_a_exp} (A) vs n={n_b_exp} (B)\n")
        if min(n_a_exp, n_b_exp) < 5:
            handle.write("- With n<5 per group, edgeR FDR estimates may be unstable for moderate effects.\n")
        handle.write("- Validate key findings with orthogonal methods or additional replicates.\n")
        handle.write("- Results reflect the quality of the public ChIP-Atlas tracks.\n")

        n_a = int((frame.get("direction", pd.Series()) == "A_enriched").sum()) if len(frame) else 0
        n_b = int((frame.get("direction", pd.Series()) == "B_enriched").sum()) if len(frame) else 0
        handle.write("\n## Summary Statistics\n\n")
        handle.write(f"- Total regions (post-QC): {len(frame)}\n")
        handle.write(f"- Enriched in {desc_a}: {n_a}\n")
        handle.write(f"- Enriched in {desc_b}: {n_b}\n")
        if "significant" in frame.columns and len(frame):
            n_sig = int(frame["significant"].sum())
            handle.write(f"- Significant (FDR<0.05): {n_sig}\n")
            sig_subset = frame[frame["significant"]]
            n_sig_a = int((sig_subset["direction"] == "A_enriched").sum())
            n_sig_b = int((sig_subset["direction"] == "B_enriched").sum())
            handle.write(f"  - Significant in {desc_a}: {n_sig_a}\n")
            handle.write(f"  - Significant in {desc_b}: {n_sig_b}\n")
        if "logFC" in frame.columns and len(frame):
            handle.write(f"- log2FC range: {frame['logFC'].min():.2f} to {frame['logFC'].max():.2f}\n")
        if "qvalue" in frame.columns and len(frame):
            handle.write(f"- Min Q-value: {frame['qvalue'].min():.2e}\n")
        if "region_size" in frame.columns and len(frame):
            handle.write(f"- Median region size: {int(frame['region_size'].median())} bp\n")

        sex_artifact_chroms = {
            w["issue"].split(" ")[0]
            for w in warnings
            if "sex chromosome confound" in w.get("issue", "")
        }

        top10 = _select_top(frame, n=10)
        if len(top10):
            handle.write("\n## Top 10 Regions (strict Q-value ranking)\n\n")
            include_genes = top10_annotations is not None and (top10_annotations.fillna("") != "").any()
            if include_genes:
                handle.write("| Rank | Location | Size (bp) | logFC | Q-value | Direction | Nearest gene |\n")
                handle.write("|------|----------|-----------|-------|---------|-----------|--------------|\n")
            else:
                handle.write("| Rank | Location | Size (bp) | logFC | Q-value | Direction |\n")
                handle.write("|------|----------|-----------|-------|---------|-----------|\n")

            for rank, (index, row) in enumerate(top10.iterrows(), 1):
                location = f"{row['chrom']}:{int(row['chromStart']):,}-{int(row['chromEnd']):,}"
                size = int(row.get("region_size", row["chromEnd"] - row["chromStart"]))
                logfc = f"{row['logFC']:.2f}" if "logFC" in row else "N/A"
                qval = f"{row['qvalue']:.2e}" if "qvalue" in row else "N/A"
                direction = row.get("direction", "unknown")
                if row["chrom"] in sex_artifact_chroms:
                    direction = f"{direction} *"
                if include_genes:
                    gene = (top10_annotations.get(index, "") if top10_annotations is not None else "") or ""
                    handle.write(f"| {rank} | {location} | {size:,} | {logfc} | {qval} | {direction} | {gene} |\n")
                else:
                    handle.write(f"| {rank} | {location} | {size:,} | {logfc} | {qval} | {direction} |\n")

            handle.write("\nRanked strictly by Q-value. Some regions may lack nearest-gene annotations, which "
                         "does not diminish their statistical significance.\n")
            if sex_artifact_chroms:
                handle.write(f"\n\\* Probable sex-chromosome artifact ({', '.join(sorted(sex_artifact_chroms))}).\n")

        handle.write("\n## Files Generated\n\n")
        handle.write("- `analysis_object.pkl` — full results object\n")
        handle.write("- `diff_regions_all.csv` — QC-filtered table\n")
        handle.write("- `diff_regions_significant.csv` — qvalue<0.05 subset (if any)\n")
        handle.write("- `diff_regions_top50.csv` — top 50 by Q-value\n")
        if len(unfiltered) != len(frame):
            handle.write("- `diff_regions_unfiltered.csv` — pre-QC table\n")
        handle.write("- `volcano_plot.{png,svg}` — logFC vs -log10(Q)\n")
        handle.write("- `chromosome_distribution.{png,svg}` — significant regions per chromosome\n")
        handle.write("- `region_size_distribution.{png,svg}` — log-scaled region-size histogram\n")
        handle.write("- `ma_plot.{png,svg}` — log2 mean count vs logFC\n")
        handle.write("- `raw_results/` — original ChIP-Atlas archive\n")

        handle.write("\n## Suggested Next Steps\n\n")
        handle.write("1. Inspect top regions in IGV using `raw_results/wabi_result.xml`\n")
        handle.write("2. Motif enrichment with HOMER or MEME-ChIP on the significant set\n")
        handle.write("3. Cross-reference with RNA-seq differential expression\n")
        handle.write("4. Region-to-gene linking via GREAT or `chip-atlas-target-genes`\n")


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pickle", required=True, help="Path to analysis_object.pkl")
    parser.add_argument("--output-dir", default=None, help="Defaults to the pickle directory")
    parser.add_argument("--qvalue", type=float, default=0.05)
    args = parser.parse_args(argv)

    with open(args.pickle, "rb") as handle:
        results = pickle.load(handle)
    output_dir = args.output_dir or os.path.dirname(os.path.abspath(args.pickle))
    export_all(results, output_dir=output_dir, qvalue_threshold=args.qvalue)
    return 0


if __name__ == "__main__":
    sys.exit(main())
