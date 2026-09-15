"""
Four-panel summary plot for ChIP-Atlas enrichment results.

Panels
------
1. Top-N enriched factors by -log10(q-value), coloured by significance
   tier and annotated with fold enrichment + overlap.
2. P-value distribution (-log10 histogram with a dashed p=0.05 marker).
3. Overlap rate vs fold enrichment scatter, coloured by -log10(q).
4. Volcano plot (log2 fold enrichment vs -log10(q)).

Plots are composited into one PNG (and SVG when possible). Each panel is
generated with plotnine + Prism theme, then rendered into a 2×2 matplotlib
grid via the canvas buffer trick so we get a single combined figure.

Inputs : the dict returned by `run_enrichment_workflow`
Outputs: `chipatlas_enrichment.{png,svg}` in `output_dir`

Example
-------
    python generate_all_plots.py --pickle enrichment_results/analysis_object.pkl
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import warnings
from typing import Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas numpy matplotlib")

try:
    from plotnine import (
        aes, annotate, coord_flip, element_text, geom_col, geom_histogram,
        geom_hline, geom_point, geom_text, geom_vline, ggplot, labs,
        scale_color_cmap, scale_color_manual, scale_fill_manual,
        scale_y_continuous, theme,
    )
    from plotnine_prism import theme_prism
except ImportError:
    sys.exit("Install with: pip install plotnine plotnine-prism")

warnings.filterwarnings("ignore")


def generate_all_plots(
    results: dict,
    output_dir: str = "enrichment_results",
    top_n: int = 15,
    base_name: str = "chipatlas_enrichment",
) -> Optional[str]:
    """Render the 4-panel figure; return the PNG path."""
    os.makedirs(output_dir, exist_ok=True)
    print("\nGenerating visualizations…")

    frame: pd.DataFrame = results["enrichment_results"]
    if frame.empty:
        print("   nothing to plot — enrichment_results is empty")
        return None

    panels = [
        _plot_top_factors(frame, top_n),
        _plot_pvalue_histogram(frame),
        _plot_overlap_scatter(frame),
        _plot_volcano(frame),
    ]

    png_path, svg_path = _combine_and_save(panels, results, output_dir, base_name)
    print("✓ All visualizations generated")
    if svg_path is not None:
        print(f"   wrote {png_path} and {svg_path}")
    else:
        print(f"   wrote {png_path} (SVG skipped)")
    return png_path


def _plot_top_factors(frame: pd.DataFrame, top_n: int):
    rankable = frame[frame["regions_with_overlaps"] >= 2].copy()
    if rankable.empty:
        rankable = frame.copy()
    rankable = rankable.head(top_n).reset_index(drop=True)
    rankable["label"] = rankable.apply(
        lambda row: f"{str(row['antigen'])[:20]} ({str(row['cell_type'])[:15]})",
        axis=1,
    )
    rankable["label"] = _unique_labels(rankable["label"].tolist())
    rankable["neg_log10_q"] = -np.log10(rankable["q_value"].clip(lower=1e-300))
    rankable["fold_str"] = rankable["fold_enrichment"].apply(
        lambda value: f"{value:.0f}×" if value < 100_000 else ">100k×"
    )
    rankable["annotation"] = rankable.apply(
        lambda row: f"{row['fold_str']} ({int(row['regions_with_overlaps'])}/{int(row['total_regions'])})",
        axis=1,
    )
    rankable["label"] = pd.Categorical(
        rankable["label"], categories=rankable["label"].iloc[::-1].tolist(),
    )
    rankable["sig_tier"] = rankable["q_value"].apply(_sig_tier)

    return (
        ggplot(rankable, aes(x="label", y="neg_log10_q", fill="sig_tier"))
        + geom_col(color="#2C5F8A", size=0.4)
        + geom_text(
            aes(label="annotation"),
            ha="left",
            nudge_y=max(rankable["neg_log10_q"]) * 0.02 if len(rankable) else 0,
            size=7,
        )
        + coord_flip()
        + scale_fill_manual(values={
            "q < 0.001": "#1B4F72",
            "q < 0.01":  "#2E86C1",
            "q < 0.05":  "#85C1E9",
            "n.s.":      "#D5D8DC",
        })
        + geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="#E74C3C", alpha=0.5)
        + scale_y_continuous(expand=(0.05, 0, 0.20, 0))
        + labs(title=f"Top {len(rankable)} enriched factors", x="", y="-log10(Q-value)", fill="")
        + theme_prism()
        + theme(figure_size=(8, 6), legend_position="bottom")
    )


def _plot_pvalue_histogram(frame: pd.DataFrame):
    data = pd.DataFrame({
        "neg_log10_p": -np.log10(frame["p_value"].clip(lower=1e-300)),
    })
    return (
        ggplot(data, aes(x="neg_log10_p"))
        + geom_histogram(bins=min(30, len(frame)), fill="#E8747C", color="#8B1A1A", alpha=0.75)
        + geom_vline(xintercept=-np.log10(0.05), linetype="dashed", color="#E74C3C")
        + annotate("text", x=-np.log10(0.05) + 0.3, y=float("inf"),
                   label="p = 0.05", ha="left", va="top", size=8, color="#E74C3C")
        + labs(title="P-value distribution", x="-log10(P-value)", y="Frequency")
        + theme_prism()
        + theme(figure_size=(8, 6))
    )


def _plot_overlap_scatter(frame: pd.DataFrame):
    data = frame.copy()
    data["neg_log10_q"] = -np.log10(data["q_value"].clip(lower=1e-300))
    data["fold_display"] = data["fold_enrichment"].clip(upper=1000)
    return (
        ggplot(data, aes(x="overlap_rate", y="fold_display", color="neg_log10_q"))
        + geom_point(size=2, alpha=0.7)
        + geom_hline(yintercept=2, linetype="dashed", color="#7F8C8D", alpha=0.5)
        + scale_color_cmap("viridis")
        + labs(
            title="Overlap rate vs fold enrichment",
            x="Overlap rate (genes with peaks / submitted)",
            y="Fold enrichment (capped at 1000×)",
            color="-log10(Q)",
        )
        + theme_prism()
        + theme(figure_size=(8, 6))
    )


def _plot_volcano(frame: pd.DataFrame):
    data = frame.copy()
    data["fold_capped"] = data["fold_enrichment"].clip(lower=1e-3, upper=1000)
    data["log2_fold"] = np.log2(data["fold_capped"])
    data["neg_log10_q"] = -np.log10(data["q_value"].clip(lower=1e-300))
    data["sig_label"] = np.where(data["q_value"] < 0.05, "Significant (q<0.05)", "Not significant")
    return (
        ggplot(data, aes(x="log2_fold", y="neg_log10_q", color="sig_label"))
        + geom_point(size=1.5, alpha=0.55)
        + geom_vline(xintercept=np.log2(2), linetype="dashed", color="#7F8C8D", alpha=0.5)
        + geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="#E74C3C", alpha=0.5)
        + scale_color_manual(values={
            "Significant (q<0.05)": "#E74C3C",
            "Not significant": "#BDC3C7",
        })
        + labs(title="Volcano plot",
               x="log2(fold enrichment, capped at 1000×)",
               y="-log10(Q-value)", color="")
        + theme_prism()
        + theme(figure_size=(8, 6), legend_position="bottom")
    )


def _sig_tier(q_value: float) -> str:
    if q_value < 0.001:
        return "q < 0.001"
    if q_value < 0.01:
        return "q < 0.01"
    if q_value < 0.05:
        return "q < 0.05"
    return "n.s."


def _unique_labels(labels):
    seen = {}
    unique = []
    for label in labels:
        if label in seen:
            seen[label] += 1
            unique.append(f"{label} #{seen[label]}")
        else:
            seen[label] = 1
            unique.append(label)
    return unique


def _combine_and_save(panels, results, output_dir, base_name):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    for axis, panel in zip(axes.flat, panels):
        rendered = panel.draw()
        rendered.set_size_inches(8, 6)
        rendered.set_dpi(300)
        rendered.canvas.draw()
        width, height = rendered.canvas.get_width_height()
        buffer = np.frombuffer(rendered.canvas.buffer_rgba(), dtype=np.uint8)
        axis.imshow(buffer.reshape(height, width, 4))
        axis.axis("off")
        plt.close(rendered)

    n_genes = len(results["input_genes"])
    n_exp = len(results["enrichment_results"])
    genome = results["parameters"].get("genome", "hg38")
    fig.suptitle(
        f"ChIP-Atlas peak enrichment — {n_genes} genes, {n_exp} experiments, genome {genome}",
        fontsize=15, fontweight="bold", y=0.995,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    png_path = os.path.join(output_dir, f"{base_name}.png")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")

    svg_path = os.path.join(output_dir, f"{base_name}.svg")
    try:
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
    except Exception:
        svg_path = None
    plt.close(fig)
    return png_path, svg_path


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pickle", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--top-n", type=int, default=15)
    args = parser.parse_args(argv)
    with open(args.pickle, "rb") as handle:
        results = pickle.load(handle)
    output_dir = args.output_dir or os.path.dirname(os.path.abspath(args.pickle))
    path = generate_all_plots(results, output_dir=output_dir, top_n=args.top_n)
    print(f"✓ PNG: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
