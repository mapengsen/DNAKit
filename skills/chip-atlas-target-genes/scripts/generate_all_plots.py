"""
Plots for ChIP-Atlas target-genes results.

Four files are produced (PNG + SVG):
* `target_genes_top_targets`            — bar plot of top-N targets
* `target_genes_score_distribution`     — log-scaled distribution of avg
* `target_genes_heatmap`                — clustered heatmap (genes × exp)
* `target_genes_string_vs_binding`      — STRING vs binding scatter
                                          (binding-rate fallback when
                                          STRING is unavailable)

Inputs : the dict returned by `run_target_genes_workflow`
Outputs: PNG / SVG files in `output_dir`
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from typing import Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
except ImportError:
    sys.exit("Install with: pip install pandas numpy matplotlib seaborn")

try:
    from plotnine import (
        aes, element_text, geom_col, geom_histogram, geom_label,
        geom_point, geom_vline, ggplot, coord_flip, labs,
        scale_color_manual, scale_fill_gradient, theme,
    )
    from plotnine_prism import theme_prism
except ImportError:
    sys.exit("Install with: pip install plotnine plotnine-prism")


def generate_all_plots(
    results: dict,
    output_dir: str = "target_genes_results",
    top_n: int = 25,
) -> Optional[str]:
    os.makedirs(output_dir, exist_ok=True)
    target_genes: pd.DataFrame = results["target_genes"]
    experiment_data: pd.DataFrame = results["experiment_data"]
    protein: str = results["protein"]
    metadata: dict = results["metadata"]

    if target_genes.empty:
        print("nothing to plot — target_genes is empty")
        return None

    base = os.path.join(output_dir, "target_genes")

    print("\nRendering plots…")
    _plot_top_targets(target_genes, top_n, protein, metadata, base + "_top_targets")
    _plot_score_distribution(target_genes, protein, base + "_score_distribution")
    _plot_heatmap(target_genes, experiment_data, protein, base + "_heatmap")
    _plot_string_vs_binding(target_genes, protein, base + "_string_vs_binding")
    print("✓ All visualizations generated")
    return base + "_top_targets.png"


def _plot_top_targets(target_genes, top_n, protein, metadata, base_path):
    top = target_genes.head(min(top_n, len(target_genes))).copy()
    top["gene"] = pd.Categorical(top["gene"], categories=top["gene"][::-1], ordered=True)
    plot = (
        ggplot(top, aes(x="gene", y="avg_score", fill="binding_rate"))
        + geom_col(alpha=0.85)
        + coord_flip()
        + scale_fill_gradient(low="#FEE08B", high="#D73027")
        + labs(
            title=f"{protein} — top {len(top)} target genes",
            subtitle=(
                f"{metadata['total_experiments']} experiments, "
                f"{metadata['genome']} ±{metadata['distance_kb']}kb"
            ),
            x="",
            y="Average MACS2 binding score",
            fill="Binding\nrate",
        )
        + theme_prism()
        + theme(
            axis_text_y=element_text(size=7 if len(top) > 15 else 8),
            plot_title=element_text(size=12, weight="bold"),
            plot_subtitle=element_text(size=9),
        )
    )
    _save_plot(plot, base_path)


def _plot_score_distribution(target_genes, protein, base_path):
    data = target_genes[target_genes["avg_score"] > 0].copy()
    if data.empty:
        return
    data["log_score"] = np.log10(data["avg_score"] + 1)
    n_500 = int((target_genes["avg_score"] >= 500).sum())
    n_100 = int((target_genes["avg_score"] >= 100).sum())
    n_50 = int((target_genes["avg_score"] >= 50).sum())

    thresholds = pd.DataFrame({
        "x": [np.log10(51), np.log10(101), np.log10(501)],
    })
    thresholds = thresholds[thresholds["x"] <= data["log_score"].max()]

    subtitle = (
        f"{len(data)} genes with binding (median: {target_genes['avg_score'].median():.0f}) | "
        f"≥500: {n_500}, ≥100: {n_100}, ≥50: {n_50}"
    )

    plot = (
        ggplot(data, aes(x="log_score"))
        + geom_histogram(bins=50, fill="#4C72B0", alpha=0.85, color="white")
        + geom_vline(
            data=thresholds, mapping=aes(xintercept="x"),
            linetype="dashed", color="#7F8C8D", alpha=0.7,
        )
        + labs(
            title=f"{protein} — binding score distribution",
            subtitle=subtitle,
            x="log10(avg_score + 1)",
            y="Genes",
        )
        + theme_prism()
        + theme(
            plot_title=element_text(size=12, weight="bold"),
            plot_subtitle=element_text(size=9),
        )
    )
    _save_plot(plot, base_path)


def _plot_heatmap(target_genes, experiment_data, protein, base_path, n_genes=20, n_experiments=20):
    top_genes = target_genes.head(min(n_genes, len(target_genes)))["gene"].values
    if len(top_genes) == 0:
        return
    exp_cols = [c for c in experiment_data.columns if c != "gene"]
    if not exp_cols:
        return

    subset = experiment_data[experiment_data["gene"].isin(top_genes)].set_index("gene")
    values = subset[exp_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    if values.empty:
        return

    column_means = values.mean(axis=0)
    selected = column_means.nlargest(min(n_experiments, len(exp_cols))).index.tolist()
    values = values[selected].reindex(top_genes)

    labels = []
    for column in selected:
        parts = column.split("|", 1)
        labels.append(parts[1] if len(parts) == 2 else column[:10])
    values.columns = labels

    total_cell_types = {column.split("|", 1)[1] if "|" in column else column for column in exp_cols}
    selected_cell_types = set(labels)

    grid = sns.clustermap(
        values,
        cmap="viridis",
        figsize=(12, max(8, len(values) * 0.4)),
        dendrogram_ratio=0.1,
        cbar_kws={"label": "MACS2 binding score"},
        row_cluster=True,
        col_cluster=True,
        xticklabels=True,
        yticklabels=True,
    )
    plt.setp(grid.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(grid.ax_heatmap.get_yticklabels(), fontsize=7)
    grid.ax_heatmap.set_xlabel("Experiment (cell type)", fontsize=10)
    grid.ax_heatmap.set_ylabel("Gene", fontsize=10)
    grid.fig.suptitle(
        f"{protein} — binding heatmap (top genes × experiments)",
        fontsize=12, fontweight="bold", y=1.02,
    )
    grid.fig.text(
        0.5, 0.99,
        f"Top {len(values)} genes × {len(selected)} highest-scoring experiments "
        f"({len(selected_cell_types)} of {len(total_cell_types)} cell types)",
        ha="center", fontsize=8, color="#555555",
    )
    _save_plot(grid.fig, base_path)


def _plot_string_vs_binding(target_genes, protein, base_path):
    data = target_genes.copy()
    has_string = (data["string_score"] > 0).any()

    if not has_string:
        plot = (
            ggplot(data, aes(x="avg_score", y="binding_rate"))
            + geom_point(alpha=0.45, size=1.5, color="#4C72B0")
            + labs(
                title=f"{protein} — binding score vs binding rate",
                subtitle="No STRING interaction data available for this protein",
                x="Average MACS2 binding score",
                y="Binding rate (fraction of experiments)",
            )
            + theme_prism()
            + theme(
                plot_title=element_text(size=12, weight="bold"),
                plot_subtitle=element_text(size=9),
            )
        )
        _save_plot(plot, base_path)
        return

    data["has_string"] = np.where(data["string_score"] > 0, "STRING interaction", "No STRING data")

    with_string = data[data["string_score"] > 0].copy()
    if len(with_string):
        with_string["combined_rank"] = (
            with_string["avg_score"].rank(ascending=False)
            + with_string["string_score"].rank(ascending=False)
        )
        top_labels = with_string.nsmallest(min(8, len(with_string)), "combined_rank")
    else:
        top_labels = pd.DataFrame()

    plot = (
        ggplot(data, aes(x="avg_score", y="string_score", color="has_string"))
        + geom_point(alpha=0.5, size=1.5)
        + scale_color_manual(
            values={"No STRING data": "#CCCCCC", "STRING interaction": "#E74C3C"},
            name="",
        )
        + labs(
            title=f"{protein} — binding vs STRING",
            subtitle=f"{len(with_string)} gene(s) with STRING interaction data",
            x="Average MACS2 binding score",
            y="STRING interaction score",
        )
        + theme_prism()
        + theme(
            plot_title=element_text(size=12, weight="bold"),
            plot_subtitle=element_text(size=9),
            legend_position="bottom",
        )
    )
    if not top_labels.empty:
        plot = plot + geom_label(
            data=top_labels,
            mapping=aes(label="gene"),
            size=7, alpha=0.8, nudge_y=20,
            show_legend=False,
        )
    _save_plot(plot, base_path)


def _save_plot(plot, base_path, width=10, height=8, dpi=300):
    png_path = base_path + ".png"
    svg_path = base_path + ".svg"
    if hasattr(plot, "save"):
        plot.save(png_path, dpi=dpi, width=width, height=height, verbose=False)
        print(f"   wrote {png_path}")
        try:
            plot.save(svg_path, width=width, height=height, verbose=False)
            print(f"   wrote {svg_path}")
        except Exception:
            print(f"   ! SVG export failed; PNG only for {os.path.basename(base_path)}")
    else:
        plot.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"   wrote {png_path}")
        try:
            plot.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
            print(f"   wrote {svg_path}")
        except Exception:
            print(f"   ! SVG export failed; PNG only for {os.path.basename(base_path)}")
        plt.close(plot)


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pickle", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args(argv)
    with open(args.pickle, "rb") as handle:
        results = pickle.load(handle)
    output_dir = args.output_dir or os.path.dirname(os.path.abspath(args.pickle))
    generate_all_plots(results, output_dir=output_dir, top_n=args.top_n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
