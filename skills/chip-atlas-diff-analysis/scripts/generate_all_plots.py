"""
Publication-style plots for ChIP-Atlas diff-analysis results.

Plots generated (PNG + SVG, plotnine + Prism theme):
1. Volcano plot          : logFC vs -log10(qvalue), coloured by direction
2. Chromosome distribution: significant regions per standard chromosome,
                           split by direction. Sex-chromosome confounds
                           flagged with an asterisk.
3. Region size histogram : log-scaled region size by direction
4. MA plot               : log2 mean count vs logFC, coloured by direction

Inputs : the dict returned by `run_diff_workflow`
Outputs: PNGs/SVGs written into `output_dir`

Example
-------
    python generate_all_plots.py --pickle diff_results/analysis_object.pkl
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from typing import Iterable, List, Optional

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas numpy")

try:
    from plotnine import (
        aes, element_text, geom_col, geom_histogram, geom_hline,
        geom_point, geom_vline, ggplot, labs, position_dodge,
        scale_color_manual, scale_fill_manual, scale_x_discrete,
        scale_x_log10, theme,
    )
    from plotnine_prism import theme_prism
except ImportError:
    sys.exit("Install with: pip install plotnine plotnine-prism")

try:
    from scripts.filter_regions import get_standard_chroms
except ImportError:
    from filter_regions import get_standard_chroms


COLOUR_A = "#E74C3C"
COLOUR_B = "#3498DB"
COLOUR_NS = "#95A5A6"


def generate_all_plots(
    results: dict,
    output_dir: str = "diff_results",
    top_n: int = 20,
) -> Optional[str]:
    """Render all four plots; return the path to `volcano_plot.png`."""
    os.makedirs(output_dir, exist_ok=True)
    print("\nGenerating visualizations…")

    frame = results.get("diff_regions")
    params = results.get("parameters", {})
    desc_a = params.get("description_a", "Group A")
    desc_b = params.get("description_b", "Group B")
    genome = params.get("genome", "hg38")
    qc_warnings = results.get("qc_warnings", [])

    if frame is None or len(frame) == 0:
        print("   no regions to plot")
        print("✓ Visualizations skipped (empty result)")
        return None

    _save(_plot_volcano(frame, desc_a, desc_b), output_dir, "volcano_plot")
    _save(
        _plot_chromosome_distribution(frame, desc_a, desc_b, qc_warnings, genome),
        output_dir, "chromosome_distribution", width=10, height=6,
    )
    _save(_plot_region_size(frame, desc_a, desc_b), output_dir, "region_size_distribution")
    _save(_plot_ma(frame, desc_a, desc_b), output_dir, "ma_plot")

    print("✓ All visualizations generated")
    return os.path.join(output_dir, "volcano_plot.png")


def _plot_volcano(frame: pd.DataFrame, desc_a: str, desc_b: str):
    data = frame.copy()
    data["neg_log10_q"] = -np.log10(data["qvalue"].clip(lower=1e-300))
    sig = data["significant"] if "significant" in data.columns else data["qvalue"] < 0.05
    data["group"] = np.where(
        ~sig, "Not significant",
        np.where(data["logFC"] > 0, f"{desc_a} enriched", f"{desc_b} enriched"),
    )

    return (
        ggplot(data, aes(x="logFC", y="neg_log10_q", color="group"))
        + geom_point(size=1.5, alpha=0.6)
        + geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="gray", size=0.5)
        + geom_vline(xintercept=0, linetype="solid", color="gray", size=0.3, alpha=0.5)
        + scale_color_manual(values={
            "Not significant": COLOUR_NS,
            f"{desc_a} enriched": COLOUR_A,
            f"{desc_b} enriched": COLOUR_B,
        })
        + labs(title="Volcano plot", x="log2 fold change",
               y="-log10(Q-value)", color="Direction")
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_chromosome_distribution(
    frame: pd.DataFrame,
    desc_a: str,
    desc_b: str,
    qc_warnings: Iterable[dict],
    genome: str,
):
    sig_frame = frame[frame["significant"]] if "significant" in frame.columns else frame[frame["qvalue"] < 0.05]
    standard = get_standard_chroms(genome)
    standard_sig = sig_frame[sig_frame["chrom"].isin(standard)]
    other_count = len(sig_frame) - len(standard_sig)

    flagged_sex = set()
    for warning in qc_warnings or []:
        if "sex chromosome confound" in warning.get("issue", ""):
            chrom_token = warning["issue"].split(" ")[0]
            flagged_sex.add(chrom_token.replace("chr", ""))

    ordering = _chrom_ordering(standard_sig["chrom"].unique())
    rows: List[dict] = []
    for chrom in ordering:
        label = chrom.replace("chr", "")
        if label in flagged_sex:
            label = f"{label}*"
        a_count = int(((standard_sig["chrom"] == chrom) & (standard_sig["direction"] == "A_enriched")).sum())
        b_count = int(((standard_sig["chrom"] == chrom) & (standard_sig["direction"] == "B_enriched")).sum())
        if a_count or b_count:
            rows.append({"chrom": label, "count": a_count, "group": f"{desc_a} enriched"})
            rows.append({"chrom": label, "count": b_count, "group": f"{desc_b} enriched"})

    if not rows:
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No significant differential regions by chromosome")
            + theme_prism()
        )

    plot_df = pd.DataFrame(rows)
    active_labels: List[str] = []
    for chrom in ordering:
        label = chrom.replace("chr", "")
        if label in flagged_sex:
            label = f"{label}*"
        if label in plot_df["chrom"].values:
            active_labels.append(label)

    title = "Significant regions by chromosome (FDR<0.05)"
    if other_count:
        title += f"\n({other_count} regions on non-standard contigs not shown)"
    if flagged_sex:
        title += "\n(* probable sex-chromosome artifact)"

    return (
        ggplot(plot_df, aes(x="chrom", y="count", fill="group"))
        + geom_col(position=position_dodge(preserve="single"), alpha=0.85)
        + scale_fill_manual(values={
            f"{desc_a} enriched": COLOUR_A,
            f"{desc_b} enriched": COLOUR_B,
        })
        + scale_x_discrete(limits=active_labels)
        + labs(title=title, x="Chromosome", y="Significant regions", fill="Direction")
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_region_size(frame: pd.DataFrame, desc_a: str, desc_b: str):
    if "region_size" not in frame.columns or frame["region_size"].dropna().empty:
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No region-size data available")
            + theme_prism()
        )
    data = frame[["region_size", "direction"]].dropna().copy()
    data = data[data["region_size"] > 0]
    data["group"] = np.where(
        data["direction"] == "A_enriched",
        f"{desc_a} enriched",
        f"{desc_b} enriched",
    )
    return (
        ggplot(data, aes(x="region_size", fill="group"))
        + geom_histogram(bins=30, alpha=0.85, position="dodge")
        + scale_x_log10()
        + scale_fill_manual(values={
            f"{desc_a} enriched": COLOUR_A,
            f"{desc_b} enriched": COLOUR_B,
        })
        + labs(title="Region size distribution", x="Region size (bp, log10)",
               y="Frequency", fill="Direction")
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _plot_ma(frame: pd.DataFrame, desc_a: str, desc_b: str):
    if "mean_count_a" not in frame.columns or "mean_count_b" not in frame.columns:
        return (
            ggplot(pd.DataFrame({"x": [0], "y": [0]}), aes("x", "y"))
            + labs(title="No count data for MA plot")
            + theme_prism()
        )
    data = frame.copy()
    mean_expr = (data["mean_count_a"] + data["mean_count_b"]) / 2
    data["log2_mean_count"] = np.log2(mean_expr.clip(lower=0.1))
    sig = data["significant"] if "significant" in data.columns else data["qvalue"] < 0.05
    data["group"] = np.where(
        ~sig, "Not significant",
        np.where(data["logFC"] > 0, f"{desc_a} enriched", f"{desc_b} enriched"),
    )
    return (
        ggplot(data, aes(x="log2_mean_count", y="logFC", color="group"))
        + geom_point(size=1.2, alpha=0.5)
        + geom_hline(yintercept=0, linetype="solid", color="gray", size=0.3, alpha=0.5)
        + scale_color_manual(values={
            "Not significant": COLOUR_NS,
            f"{desc_a} enriched": COLOUR_A,
            f"{desc_b} enriched": COLOUR_B,
        })
        + labs(title="MA plot", x="log2 mean normalized count",
               y="log2 fold change", color="Direction")
        + theme_prism()
        + theme(legend_position="bottom")
    )


def _chrom_ordering(chromosomes: Iterable[str]) -> List[str]:
    def sort_key(chrom: str):
        name = chrom.replace("chr", "")
        if name.isdigit():
            return (0, int(name), "")
        return (1, 0, name)

    return sorted(chromosomes, key=sort_key)


def _save(plot, output_dir: str, base_name: str, width: int = 8, height: int = 6) -> None:
    png_path = os.path.join(output_dir, f"{base_name}.png")
    plot.save(png_path, dpi=300, width=width, height=height, verbose=False)
    print(f"   saved {png_path}")
    svg_path = os.path.join(output_dir, f"{base_name}.svg")
    try:
        plot.save(svg_path, width=width, height=height, verbose=False)
        print(f"   saved {svg_path}")
    except Exception:
        print(f"   ! SVG export failed for {base_name}; PNG only")


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pickle", required=True, help="Path to analysis_object.pkl")
    parser.add_argument("--output-dir", default="diff_results")
    args = parser.parse_args(argv)

    with open(args.pickle, "rb") as handle:
        results = pickle.load(handle)
    path = generate_all_plots(results, output_dir=args.output_dir)
    print(f"✓ Plots written to {path or args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
