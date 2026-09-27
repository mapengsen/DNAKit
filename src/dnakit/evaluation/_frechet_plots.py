"""Optional PNG figures for the Fréchet DNA distance evaluation."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization.figure import save_figure


def require_plot_dependencies(*, pca: bool) -> None:
    """Check optional plotting dependencies before model inference starts."""
    try:
        importlib.import_module("matplotlib.pyplot")
    except ImportError as exc:
        raise BackendUnavailableError(
            "Fréchet image output requires Matplotlib.",
            code="MISSING_FRECHET_PLOT_DEPENDENCY",
            hint='Install the optional dependency with: python -m pip install "dnakit[bio-plots]"',
        ) from exc
    if pca:
        try:
            importlib.import_module("sklearn.decomposition")
        except ImportError as exc:
            raise BackendUnavailableError(
                "Fréchet PCA image output requires scikit-learn.",
                code="MISSING_FRECHET_PCA_DEPENDENCY",
                hint='Install the optional dependency with: python -m pip install "dnakit[neural]"',
            ) from exc


def save_distance_bar(
    distance: float,
    target: Path,
    *,
    encoder: str,
    overwrite: bool,
) -> Path:
    """Save one scalar Fréchet distance as a labeled bar chart."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.bar(["Left vs right"], [distance], width=0.48, color="#287d8e")
    upper = max(1.0, distance * 1.2)
    ax.set_ylim(0.0, upper)
    ax.text(0, distance + upper * 0.025, f"{distance:.6g}", ha="center", va="bottom")
    ax.set_ylabel("Fréchet distance")
    ax.set_title(f"Fréchet DNA distance · {encoder}")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    return save_figure(fig, target, overwrite=overwrite, dpi=180)


def save_pca_scatter(
    matrix: Any,
    left_count: int,
    target: Path,
    *,
    encoder: str,
    overwrite: bool,
) -> Path:
    """Fit one PCA to both groups in the same space used by the distance."""
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.decomposition import PCA

    vectors = np.asarray(matrix, dtype=np.float64)
    component_count = min(2, vectors.shape[0] - 1, vectors.shape[1])
    if np.all(vectors == vectors[0]):
        coordinates = np.zeros((vectors.shape[0], 2), dtype=np.float64)
        explained = np.zeros(2, dtype=np.float64)
    else:
        fitted = PCA(n_components=component_count, svd_solver="auto", random_state=0)
        projected = fitted.fit_transform(vectors)
        coordinates = np.zeros((vectors.shape[0], 2), dtype=np.float64)
        coordinates[:, :component_count] = projected
        explained = np.zeros(2, dtype=np.float64)
        explained[:component_count] = np.nan_to_num(fitted.explained_variance_ratio_)

    fig, ax = plt.subplots(figsize=(6.5, 5.1))
    for group, points, color in (
        ("Left", coordinates[:left_count], "#287d8e"),
        ("Right", coordinates[left_count:], "#dd7045"),
    ):
        ax.scatter(
            points[:, 0],
            points[:, 1],
            s=42,
            color=color,
            alpha=0.8,
            edgecolors="white",
            linewidths=0.5,
            label=f"{group} (n={len(points)})",
        )
    ax.set_xlabel(f"PC1 ({explained[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1%} variance)")
    ax.set_title(f"DNA representation PCA · {encoder}")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    return save_figure(fig, target, overwrite=overwrite, dpi=180)


__all__ = ["require_plot_dependencies", "save_distance_bar", "save_pca_scatter"]
