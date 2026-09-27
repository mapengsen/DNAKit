"""Optional PNG output for collection evaluation results."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendUnavailableError

from .results import EvaluationReport


def _pyplot() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise BackendUnavailableError(
            "Plot output requires matplotlib.",
            hint='Install it with: python -m pip install "dnakit[bio-plots]"',
        ) from exc
    return plt


def _save(figure: Any, plt: Any, filename: str) -> Path:
    target = Path.cwd() / filename
    try:
        figure.tight_layout()
        figure.savefig(target, dpi=190)
    finally:
        plt.close(figure)
    return target


def _finish_axes(ax: Any) -> None:
    ax.grid(axis="y", color="#cbd5e1", linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_uniqueness_chart(
    score: float,
    unique_count: int,
    record_count: int,
    *,
    equivalence: str,
) -> Path:
    """Save the configured uniqueness ratio as one bar chart in the working directory."""
    plt = _pyplot()
    from matplotlib.ticker import PercentFormatter

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    if record_count:
        ax.bar([0], [1], width=0.58, color="#e2e8f0")
        ax.bar([0], [score], width=0.58, color="#0f766e")
        ax.text(
            0,
            score / 2 if score else 0.05,
            f"{score:.0%}",
            ha="center",
            va="center",
            color="white" if score >= 0.2 else "#0f172a",
            fontsize=18,
            weight="bold",
        )
        ax.text(
            0,
            1.035,
            f"{unique_count} distinct / {record_count} records",
            ha="center",
            va="bottom",
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No records (score 1 by convention)",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(0, 1.13)
    ax.set_xticks([0], ["Evaluated set"])
    ax.set_ylabel("Distinct groups / all records")
    ax.set_title("Exact uniqueness" if equivalence == "exact" else f"Uniqueness ({equivalence})")
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    _finish_axes(ax)
    return _save(fig, plt, "uniqueness_ratio.png")


def save_diversity_charts(
    mean_distance: float | None, pair_distances: list[float]
) -> tuple[Path, Path]:
    """Save mean distance and all-pairs distance distribution in the working directory."""
    plt = _pyplot()
    max_distance = max(pair_distances, default=0.0)

    mean_fig, mean_ax = plt.subplots(figsize=(5.2, 4.2))
    if mean_distance is None:
        mean_ax.text(
            0.5,
            0.5,
            "Undefined: fewer than two records",
            ha="center",
            va="center",
            transform=mean_ax.transAxes,
        )
    else:
        mean_ax.bar([0], [mean_distance], width=0.58, color="#0f766e")
        mean_ax.text(
            0,
            mean_distance + max(max_distance, 1.0) * 0.025,
            f"{mean_distance:.2f} edits",
            ha="center",
            va="bottom",
            fontsize=14,
            weight="bold",
        )
    mean_ax.set_xlim(-0.8, 0.8)
    mean_ax.set_ylim(0, max(1.0, max_distance * 1.12))
    mean_ax.set_xticks([0], ["Evaluated set"])
    mean_ax.set_ylabel("Mean distance (edit operations)")
    mean_ax.set_title("Mean pairwise Levenshtein distance")
    _finish_axes(mean_ax)
    mean_path = _save(mean_fig, plt, "diversity_mean_bar.png")

    hist_fig, hist_ax = plt.subplots(figsize=(7.2, 4.2))
    if pair_distances:
        minimum = min(pair_distances)
        maximum = max_distance
        integer_distances = all(distance.is_integer() for distance in pair_distances)
        if integer_distances and maximum - minimum <= 50:
            bins: int | list[float] = [
                distance - 0.5 for distance in range(int(minimum), int(maximum) + 2)
            ]
        else:
            bins = min(50, max(10, math.ceil(math.sqrt(len(pair_distances)))))
        hist_ax.hist(pair_distances, bins=bins, color="#3b82f6", edgecolor="white")
        if mean_distance is not None:
            hist_ax.axvline(
                mean_distance,
                color="#dc2626",
                linestyle="--",
                linewidth=2,
                label=f"Mean = {mean_distance:.2f}",
            )
            hist_ax.legend(loc="upper right")
    else:
        hist_ax.text(
            0.5, 0.5, "No sequence pairs", ha="center", va="center", transform=hist_ax.transAxes
        )
    hist_ax.set_xlabel("Pairwise Levenshtein distance (edit operations)")
    hist_ax.set_ylabel("Number of sequence pairs")
    hist_ax.set_title("All pairwise distances")
    _finish_axes(hist_ax)
    hist_path = _save(hist_fig, plt, "diversity_pairwise_distance_histogram.png")
    return mean_path, hist_path


def _novelty_histogram_bins(values: list[float]) -> int | list[float]:
    minimum = min(values)
    maximum = max(values)
    if all(value.is_integer() for value in values) and maximum - minimum <= 50:
        return [value - 0.5 for value in range(int(minimum), int(maximum) + 2)]
    return min(50, max(10, math.ceil(math.sqrt(len(values)))))


def save_novelty_charts(
    result: EvaluationReport,
    nearest_distances: list[float],
    *,
    batch_labels: tuple[str, ...],
) -> tuple[Path, Path, Path]:
    """Save summary, per-query novelty, and nearest-reference edit-distance PNGs."""
    plt = _pyplot()
    from matplotlib.ticker import PercentFormatter

    entries = result.entries
    if len(batch_labels) != len(entries):
        raise ValueError("batch_labels must match the number of novelty entries.")
    calculation = result.parameters["novelty_calculation"]
    novelty_values = [float(entry.metrics["novelty"]) for entry in entries]
    exact_matches = [entry.metrics["exact_match"] is True for entry in entries]
    grouped: dict[str, list[int]] = {}
    for index, label in enumerate(batch_labels):
        grouped.setdefault(label, []).append(index)
    if not grouped:
        grouped["Evaluated set"] = []
    labels = list(grouped)
    group_means = [
        math.fsum(novelty_values[index] for index in indices) / len(indices) if indices else None
        for indices in grouped.values()
    ]
    group_exact_rates = [
        sum(exact_matches[index] for index in indices) / len(indices) if indices else None
        for indices in grouped.values()
    ]

    summary_fig, summary_axes = plt.subplots(1, 2, figsize=(max(10.5, 1.2 * len(labels) + 6), 4.3))
    for ax, values, color, is_rate in (
        (summary_axes[0], group_means, "#0f766e", False),
        (summary_axes[1], group_exact_rates, "#dc6b4a", True),
    ):
        for index, value in enumerate(values):
            if value is None:
                continue
            ax.bar(index, value, width=0.65, color=color)
            label = f"{value:.1%}" if is_rate or calculation == "similarity" else f"{value:.2f}"
            offset = 0.025 if is_rate or calculation == "similarity" else max(value * 0.025, 0.05)
            ax.text(index, value + offset, label, ha="center", va="bottom")
        if all(value is None for value in values):
            ax.text(0.5, 0.5, "No queries", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks(range(len(labels)), labels, rotation=30 if len(labels) > 3 else 0)
        ax.set_xlim(-0.7, len(labels) - 0.3)
        if is_rate or calculation == "similarity":
            ax.set_ylim(0, 1.08)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
        else:
            upper = max((value for value in values if value is not None), default=0.0)
            ax.set_ylim(0, max(1.0, upper * 1.18))
        _finish_axes(ax)
    summary_axes[0].set_title(
        "Mean novelty"
        if calculation == "similarity"
        else "Mean nearest-reference Levenshtein distance"
    )
    summary_axes[0].set_ylabel("Mean novelty" if calculation == "similarity" else "Edit operations")
    summary_axes[1].set_title("Exact match rate")
    summary_axes[1].set_ylabel("Fraction of queries")
    summary_path = _save(summary_fig, plt, "novelty_batch_summary.png")

    novelty_fig, novelty_ax = plt.subplots(figsize=(8.3, 4.6))
    if novelty_values:
        novelty_bins: int | list[float] = (
            [index / 20 for index in range(21)]
            if calculation == "similarity"
            else _novelty_histogram_bins(novelty_values)
        )
        novelty_ax.hist(novelty_values, bins=novelty_bins, color="#3b82f6", edgecolor="white")
        if result.score is not None:
            novelty_ax.axvline(
                result.score,
                color="#111827",
                linestyle="--",
                linewidth=1.7,
                label=f"Mean = {result.score:.2f}",
            )
            novelty_ax.legend(loc="upper right")
        if calculation == "similarity":
            novelty_ax.set_xlim(0, 1)
    else:
        novelty_ax.text(
            0.5, 0.5, "No queries", ha="center", va="center", transform=novelty_ax.transAxes
        )
    novelty_ax.set_xlabel(
        "Per-query novelty = 1 - nearest similarity"
        if calculation == "similarity"
        else "Per-query novelty (edit operations)"
    )
    novelty_ax.set_ylabel("Number of queries")
    novelty_ax.set_title("Per-query novelty distribution")
    _finish_axes(novelty_ax)
    novelty_path = _save(novelty_fig, plt, "novelty_distribution.png")

    distance_fig, distance_ax = plt.subplots(figsize=(8.3, 4.6))
    if nearest_distances:
        distance_ax.hist(
            nearest_distances,
            bins=_novelty_histogram_bins(nearest_distances),
            color="#3b82f6",
            edgecolor="white",
        )
        average = math.fsum(nearest_distances) / len(nearest_distances)
        distance_ax.axvline(
            average,
            color="#111827",
            linestyle="--",
            linewidth=1.7,
            label=f"Mean = {average:.2f}",
        )
        distance_ax.legend(loc="upper right")
    else:
        distance_ax.text(
            0.5,
            0.5,
            "No reference distances available",
            ha="center",
            va="center",
            transform=distance_ax.transAxes,
        )
    distance_ax.set_xlabel("Nearest-reference Levenshtein distance (edit operations)")
    distance_ax.set_ylabel("Number of queries")
    distance_ax.set_title("Nearest-reference edit distance")
    _finish_axes(distance_ax)
    distance_path = _save(distance_fig, plt, "novelty_nearest_levenshtein.png")
    return summary_path, novelty_path, distance_path
