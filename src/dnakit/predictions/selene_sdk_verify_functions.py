"""Selene plots of regulatory model scores that DNAKit draws without the calls that fail.

selene_sdk 0.6.0 cannot draw two of its interpret plots with the installed matplotlib 3.11 and plotly 6.9:

- ``sequence_logo`` builds the letter glyphs and then attaches a path effect whose ``draw_path`` calls
  ``Affine2D.identity``, which is gone; rendering the figure fails. Here the glyphs are placed with the
  same bar-to-glyph geometry and drawn as ordinary patches in data coordinates. The glyph outlines are
  selene's own ``_SVG_PATHS`` table, read at call time.
- ``variant_diffs_scatter_plot`` builds its layout with ``titlefont``, which plotly 6 rejects. The same
  layout is written with ``title`` objects. Upstream also indexes the unsorted scores when a feature
  filter is given, which pairs labels with the wrong rows; here the rows and the feature columns are
  both taken in label order.

Both plots write a new file (PNG for the logo, HTML for the scatter) and return the data drawn.
selene_sdk, matplotlib and plotly are imported only when a function runs.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import math
from collections.abc import Sequence
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.dadi_functions import _new_file

_SELENE_EXTRA = "bio-selene"
_DNA_BASES = ("A", "C", "G", "T")
# seaborn's "Set1" palette, first four colours (selene's default when no colour scheme is given)
_DEFAULT_COLOURS = ("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
_ORDERS = ("alpha", "value")
_PLOT_DPI = 100
_MAX_POSITIONS = 5000


def _selene_module(dotted: str) -> Any:
    return require(dotted, _SELENE_EXTRA)


def _version() -> str:
    try:
        return importlib.metadata.version("selene-sdk")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - selene is installed here
        return "unknown"


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must hold finite numbers.")
    return float(value)


def _matrix(value: object, name: str, columns: int) -> list[list[float]]:
    """A non-empty rectangular nested list of finite numbers with exactly ``columns`` columns."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of rows.")
    rows: list[list[float]] = []
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise ConfigurationError(f"{name} must be a list of rows of numbers.")
        if len(row) != columns:
            raise ConfigurationError(f"{name} must have {columns} columns, found {len(row)}.")
        rows.append([_number(item, f"{name} entry") for item in row])
    if len(rows) > _MAX_POSITIONS:
        raise ConfigurationError(f"{name} has more than {_MAX_POSITIONS} positions.")
    return rows


def _colours(value: object) -> list[str]:
    if value is None:
        return list(_DEFAULT_COLOURS)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("color_scheme must be a list of colour names or hex codes.")
    colours = [str(item) for item in value]
    if len(colours) < len(_DNA_BASES) or any(not item for item in colours):
        raise ConfigurationError("color_scheme must give one colour for each of A, C, G and T.")
    return colours


def _png_target(image_path: object) -> Path:
    target = _new_file(image_path, "image_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("image_path must end with .png.")
    return target


def _html_target(output_path: object) -> Path:
    target = _new_file(output_path, "output_path")
    if target.suffix.lower() != ".html":
        raise ConfigurationError("output_path must end with .html.")
    return target


def _stack(scores: Any, np: Any, bases: Any, order: str) -> None:
    """Sort the bases of each position in place, as selene's sequence_logo does (order 'value')."""
    base_names = np.array(_DNA_BASES)
    if order == "alpha":
        for i in range(scores.shape[0]):
            bases[i, :] = base_names[i]
        return
    if np.sum(scores < 0) != 0:
        sorted_scores = np.zeros_like(scores)
        for j in range(scores.shape[1]):
            div = np.sum(scores[:, j] < 0.0)
            negative_idx = np.argwhere(scores[:, j] < 0.0).flatten()
            negative_sort_idx = np.argsort(scores[negative_idx, j], axis=None)
            sorted_scores[:div, j] = scores[negative_idx[negative_sort_idx], j]
            bases[:div, j] = base_names[negative_idx[negative_sort_idx]].flatten()
            positive_idx = np.argwhere(scores[:, j] >= 0.0).flatten()
            positive_sort_idx = np.argsort(scores[positive_idx, j], axis=None)
            sorted_scores[div:, j] = scores[positive_idx[positive_sort_idx], j]
            bases[div:, j] = base_names[positive_idx[positive_sort_idx]].flatten()
        scores[:, :] = sorted_scores
    else:
        for j in range(scores.shape[1]):
            sort_idx = np.argsort(scores[:, j], axis=None)[::-1]
            bases[:, j] = base_names[sort_idx]
            scores[:, j] = scores[sort_idx, j]


def _offsets(scores: Any, np: Any) -> tuple[Any, Any]:
    positive = np.zeros_like(scores)
    negative = np.zeros_like(scores)
    for i in range(scores.shape[0] - 1):
        y = scores[i, :]
        if i > 0:
            negative[i + 1, :] = negative[i, :]
            positive[i + 1, :] = positive[i, :]
        neg_idx = np.argwhere(y < 0.0)
        pos_idx = np.argwhere(y >= 0.0)
        negative[i + 1, neg_idx] += y[neg_idx]
        positive[i + 1, pos_idx] += y[pos_idx]
    return positive, negative


def plot_sequence_logo_scores(
    score_matrix: Sequence[Sequence[float]],
    image_path: str | Path,
    *,
    order: str = "value",
    width: float = 1.0,
    color_scheme: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a DNA sequence logo from an L x 4 score matrix (for example attribution scores) to PNG.

    ``score_matrix`` has one row per position and one column per base (A, C, G, T). ``order='value'``
    stacks the bases of each position by value, the largest at the bottom, and negative scores below
    zero; ``order='alpha'`` keeps the alphabet order. Each stacked bar is filled with its base's
    glyph outline, which is selene's own. ``color_scheme`` defaults to the first four colours of
    seaborn's Set1 palette, as selene does.
    """
    if order not in _ORDERS:
        raise ConfigurationError(f"order must be one of {list(_ORDERS)}.")
    width_value = _number(width, "width")
    if width_value < 0:
        raise ConfigurationError("width must be non-negative.")
    target = _png_target(image_path)
    colours = _colours(color_scheme)
    np = require("numpy", _SELENE_EXTRA)
    plt = require("matplotlib.pyplot", _SELENE_EXTRA)
    matplotlib_path = require("matplotlib.path", _SELENE_EXTRA)
    patches = require("matplotlib.patches", _SELENE_EXTRA)
    transforms = require("matplotlib.transforms", _SELENE_EXTRA)
    vis = _selene_module("selene_sdk.interpret.vis")
    glyphs = vis._SVG_PATHS
    scores = np.array(_matrix(score_matrix, "score_matrix", len(_DNA_BASES)), dtype=float).T.copy()
    n_bases, length = scores.shape
    bases = np.empty(scores.shape, dtype=object)
    bases[:, :] = "?"
    _stack(scores, np, bases, order)
    positive_offsets, negative_offsets = _offsets(scores, np)
    figure, ax = plt.subplots(figsize=(max(length * 0.4, 3.0), 2.5))
    try:
        for i in range(n_bases):
            x_coords = np.arange(length) + 0.5
            y = scores[i, :]
            offsets = np.zeros(length)
            negative_idx = np.argwhere(y < 0.0)
            positive_idx = np.argwhere(y >= 0.0)
            offsets[negative_idx] = negative_offsets[i, negative_idx]
            offsets[positive_idx] = positive_offsets[i, positive_idx]
            bars = ax.bar(x_coords, y, color="black", width=width_value, bottom=offsets)
            for j, bar in enumerate(bars):
                base = bases[i, j]
                bar.set_color(colours[_DNA_BASES.index(base)])
                bar.set_edgecolor(None)
        new_patches = []
        for index, bar in enumerate(list(ax.patches)):
            base = bases[index // length, index % length]
            vertices, codes = glyphs[base][0], glyphs[base][1]
            glyph = matplotlib_path.Path(vertices, codes)
            b_x, b_y, b_w, b_h = bar.get_bbox().bounds
            t_x, t_y, t_w, t_h = glyph.get_extents().bounds
            scale = (b_w / t_w, b_h / t_h)
            patch = patches.PathPatch(glyph, facecolor=bar.get_facecolor(), lw=0.0)
            bar.set_facecolor("none")
            # glyph units -> data units: move the glyph's lower-left corner to the origin, scale it to
            # the bar box, then move it to the bar's lower-left corner (a negative bar flips the glyph)
            placed = (
                transforms.Affine2D().translate(-t_x, -t_y).scale(*scale).translate(b_x, b_y)
                + ax.transData
            )
            patch.set_transform(placed)
            new_patches.append(patch)
        for patch in new_patches:
            ax.add_artist(patch)
        ax.set_xlim(0, length)
        ax.set_xticks(np.arange(length) + 0.5)
        ax.set_xticklabels([str(item) for item in range(length)])
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except (ValueError, TypeError, AttributeError, ZeroDivisionError, KeyError) as exc:
        raise ConfigurationError(f"matplotlib could not draw the sequence logo: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_sequence_logo_scores",
        "image_path": str(target),
        "image_format": "png",
        "order": order,
        "bases": [[str(item) for item in row] for row in bases],
        "scores": [[float(item) for item in row] for row in scores],
        "provenance": {
            "package": f"selene-sdk {_version()}",
            "convention": "sequence logo drawn as in selene interpret sequence_logo; glyph path effect replaced "
            "by an ordinary patch transform; matplotlib (Agg)",
        },
    }


def _ordered_variants(labels: Sequence[Sequence[Any]]) -> tuple[list[tuple[str, int, str, str]], list[int]]:
    """Variant labels sorted by chromosome then position, as selene's ordered_variants_and_indices."""
    vis = _selene_module("selene_sdk.interpret.vis")
    groups: dict[str, list[tuple[int, int, str, str]]] = {}
    for index, label in enumerate(labels):
        if len(label) != 5:
            raise ConfigurationError("each label must be (chrom, pos, name, ref, alt).")
        chrom, pos, _name, ref, alt = label
        groups.setdefault(str(chrom), []).append((index, int(pos), str(ref), str(alt)))
    ordered_labels: list[tuple[str, int, str, str]] = []
    ordered_indices: list[int] = []
    for chrom in sorted(groups, key=vis.sort_standard_chrs):
        for index, pos, ref, alt in sorted(groups[chrom], key=lambda item: item[1:]):
            ordered_indices.append(index)
            ordered_labels.append((chrom, pos, ref, alt))
    return ordered_labels, ordered_indices


def plot_variant_diffs_scatter(
    scores: Sequence[Sequence[float]],
    labels: Sequence[Sequence[Any]],
    features: Sequence[str],
    output_path: str | Path,
    *,
    feature_indices: Sequence[int] | None = None,
    nth_percentile: float | None = None,
    threshold_line: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scatter each variant's largest absolute difference score as a point, written to a new HTML file.

    ``scores`` is a variants by features matrix of absolute difference scores; ``labels`` gives each
    variant as (chrom, pos, name, ref, alt) and ``features`` names the model's features.
    ``feature_indices`` keeps only those feature columns before taking the maximum. ``nth_percentile``
    keeps the variants whose maximum is at or above that percentile. ``threshold_line`` draws a
    horizontal line. Points are ordered by chromosome and position, as selene does by default.
    """
    if isinstance(scores, (str, bytes)) or not isinstance(scores, Sequence) or not scores:
        raise ConfigurationError("scores must be a non-empty variants by features matrix.")
    if isinstance(features, (str, bytes)) or not isinstance(features, Sequence):
        raise ConfigurationError("features must be a list of feature names.")
    if len(labels) != len(scores):
        raise ConfigurationError("labels must have one entry per row of scores.")
    matrix = _matrix(scores, "scores", len(features))
    target = _html_target(output_path)
    np = require("numpy", _SELENE_EXTRA)
    go = require("plotly.graph_objs", _SELENE_EXTRA)
    plotly_offline = require("plotly.offline", _SELENE_EXTRA)
    ordered_labels, indices = _ordered_variants(labels)
    ordered = np.array(matrix, dtype=float)[indices, :]
    if feature_indices is not None:
        if isinstance(feature_indices, (str, bytes)) or not isinstance(feature_indices, Sequence):
            raise ConfigurationError("feature_indices must be a list of column numbers.")
        chosen = [int(item) for item in feature_indices]
        if not chosen or any(item < 0 or item >= len(features) for item in chosen):
            raise ConfigurationError("feature_indices must name existing feature columns.")
        ordered = ordered[:, chosen]
    maxima = np.amax(ordered, axis=1)
    kept = list(range(len(ordered_labels)))
    if nth_percentile is not None:
        percentile = _number(nth_percentile, "nth_percentile")
        if not 0 <= percentile <= 100:
            raise ConfigurationError("nth_percentile must be between 0 and 100.")
        cut = np.percentile(maxima, percentile)
        kept = [i for i in kept if maxima[i] >= cut]
        if not kept:
            raise ConfigurationError("no variant is at or above that percentile.")
    display = [ordered_labels[i] for i in kept]
    values = [float(maxima[i]) for i in kept]
    text_labels = [
        f"max diff score: {value}<br />{chrom} {pos}, {ref}/{alt}"
        for (chrom, pos, ref, alt), value in zip(display, values, strict=True)
    ]
    label_x = [f"{chrom} {pos}" for chrom, pos, _ref, _alt in display]
    layout_options: dict[str, Any] = {
        "title": {"text": "Max probability difference scores"},
        "hovermode": "closest",
        "hoverlabel": {"font": {"size": 16}},
        "xaxis": {
            "title": {"text": "Genome coordinates", "font": {"family": "Arial, sans-serif", "size": 16}},
            "showticklabels": True,
            "tickangle": 35,
            "nticks": 25,
            "tickmode": "auto",
            "automargin": True,
        },
        "yaxis": {
            "title": {"text": "Absolute difference", "font": {"family": "Arial, sans-serif", "size": 16}},
        },
    }
    if threshold_line is not None:
        threshold = _number(threshold_line, "threshold_line")
        layout_options["shapes"] = [
            {
                "type": "line",
                "x0": label_x[0],
                "y0": threshold,
                "x1": label_x[-1],
                "y1": threshold,
                "line": {"color": "rgba(255, 99, 71, 1)", "width": 2},
            }
        ]
    trace = go.Scatter(
        x=label_x,
        y=values,
        mode="markers",
        marker={"color": "#39CCCC", "line": {"width": 1}},
        text=text_labels,
        hoverinfo="text",
    )
    figure = go.Figure(data=[trace], layout=go.Layout(**layout_options))
    try:
        plotly_offline.plot(figure, filename=str(target), auto_open=False)
    except (ValueError, TypeError, OSError) as exc:
        raise ConfigurationError(f"plotly could not write the scatter plot: {exc}") from exc
    return {
        "status": "completed",
        "operation": "plot_variant_diffs_scatter",
        "output_path": str(target),
        "output_format": "html",
        "features": [str(features[i]) for i in (feature_indices or range(len(features)))],
        "variants": [
            {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "max_abs_diff": value}
            for (chrom, pos, ref, alt), value in zip(display, values, strict=True)
        ],
        "provenance": {
            "package": f"selene-sdk {_version()}",
            "convention": "variant max absolute difference scatter as in selene variant_diffs_scatter_plot, "
            "titles written with plotly 6 title objects",
        },
    }


__all__ = ["plot_sequence_logo_scores", "plot_variant_diffs_scatter"]

_register_public_verbose(globals())
