"""Category and histogram counts normalized by the number of valid samples."""

from __future__ import annotations

import colorsys
import csv
import math
import xml.etree.ElementTree as ET
from bisect import bisect_right
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import replace
from itertools import pairwise
from pathlib import Path
from typing import Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core.issues import Issue
from dnakit.exceptions import ConfigurationError, InputFormatError

from ._svg import SVGBuilder, number, validate_xml_text
from .config import DistributionPlotConfig
from .results import SVGArtifact

DistributionTask = Literal["classification", "regression", "probability"]
DistributionValue = str | int | float | None


def _bar_colors(count: int, config: DistributionPlotConfig) -> list[str]:
    if config.bar_colors is not None:
        return [config.bar_colors[index % len(config.bar_colors)] for index in range(count)]
    colors = [config.theme.accent]
    candidates = ("#f28e2b", "#59a14f", "#e15759", "#b07aa1", "#76b7b2", "#edc949")
    used = {config.theme.accent.lower()}
    index = 0
    while len(colors) < count:
        if index < len(candidates):
            color = candidates[index]
        else:
            # Extend the palette deterministically without repeating RGB colors.
            rgb = colorsys.hls_to_rgb((index * 0.61803398875) % 1, 0.5 + index % 3 * 0.06, 0.65)
            color = "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)
        if color.lower() not in used:
            colors.append(color)
            used.add(color.lower())
        index += 1
    return colors


def _validate_task(task: str, config: DistributionPlotConfig) -> None:
    if task not in ("classification", "regression", "probability"):
        raise ConfigurationError("task must be 'classification', 'regression' or 'probability'.")
    if task == "classification" and config.bin_edges is not None:
        raise ConfigurationError("bin_edges only applies to regression or probability.")
    if task != "classification" and config.categories is not None:
        raise ConfigurationError("categories only applies to classification.")


def _clean_values(
    values: Iterable[DistributionValue], task: DistributionTask, limit: int
) -> tuple[list[str], list[float], int, int]:
    labels: list[str] = []
    numbers: list[float] = []
    total = missing = 0
    if isinstance(values, (str, bytes)):
        raise TypeError("values must be an iterable of samples, not a single string.")
    for total, value in enumerate(values, start=1):
        if total > limit:
            raise ConfigurationError("Distribution exceeds max_values; no samples were truncated.")
        if value is None or (isinstance(value, str) and not value.strip()):
            missing += 1
            continue
        if isinstance(value, bool) or not isinstance(value, (str, int, float)):
            raise InputFormatError(f"Sample {total} must be a label or a finite number.")
        if task == "classification":
            label = str(value).strip()
            if isinstance(value, float) and not math.isfinite(value):
                raise InputFormatError(f"Sample {total} is not finite.")
            if len(label) > 1024:
                raise InputFormatError(f"Sample {total} label exceeds 1024 characters.")
            validate_xml_text(label, name="category")
            labels.append(label)
        else:
            try:
                numeric = float(value)
            except (ValueError, OverflowError) as exc:
                raise InputFormatError(
                    f"Sample {total} is not a finite numeric value: {value!r}."
                ) from exc
            if not math.isfinite(numeric):
                raise InputFormatError(f"Sample {total} is not finite: {value!r}.")
            if task == "probability" and not 0 <= numeric <= 1:
                raise InputFormatError(
                    f"Sample {total} probability must be in [0, 1]; got {numeric}.",
                    code="INVALID_PROBABILITY",
                    context={"sample_index": total - 1, "value": numeric},
                    hint="Select a probability column, or use task='regression' for logits/scores.",
                )
            numbers.append(numeric)
    if total == missing:
        raise InputFormatError("No valid samples remain after skipping blank values.")
    return labels, numbers, total, missing


def _histogram(
    values: list[float], config: DistributionPlotConfig
) -> tuple[list[float], list[int]]:
    low, high = min(values), max(values)
    if config.bin_edges is not None:
        edges = list(config.bin_edges)
        if low < edges[0] or high > edges[-1]:
            raise InputFormatError("bin_edges must cover every valid sample; none may be dropped.")
    else:
        if low == high:
            # A constant column is one occupied interval, regardless of requested bin count.
            padding = max(0.5, abs(low) * 0.01)
            edges = [low - padding, high + padding]
        else:
            # Convex interpolation avoids overflow in high - low for very large values.
            edges = [
                low * (1 - index / config.bins) + high * (index / config.bins)
                for index in range(config.bins + 1)
            ]
            edges[0], edges[-1] = low, high
        if any(not math.isfinite(edge) for edge in edges) or any(
            left >= right for left, right in pairwise(edges)
        ):
            raise ConfigurationError(
                "Cannot represent distinct finite bin edges; use fewer bins or explicit edges."
            )
    counts = [0] * (len(edges) - 1)
    for value in values:
        index = min(bisect_right(edges, value) - 1, len(counts) - 1)
        counts[index] += 1
    return edges, counts


def _render_distribution(
    labels: list[str],
    counts: list[int],
    edges: list[float],
    *,
    task: DistributionTask,
    total: int,
    missing: int,
    config: DistributionPlotConfig,
) -> SVGArtifact:
    valid = total - missing
    proportions = [count / valid for count in counts]
    colors = _bar_colors(len(counts), config)
    metadata: dict[str, object] = {
        "task": task,
        "normalization": "count / valid_count",
        "total_count": total,
        "valid_count": valid,
        "missing_count": missing,
        "labels": labels,
        "counts": counts,
        "proportions": proportions,
        "bar_colors": colors,
        "bin_edges": edges,
        "bin_closure": "[left, right), final bin [left, right]" if edges else None,
    }
    theme = config.theme
    builder = SVGBuilder(
        config.size,
        config.size,
        title=config.title or "Sample distribution",
        description=f"{task}: {valid} valid samples, {missing} blank samples skipped.",
        kind="sample_distribution",
        theme=theme,
        metadata=metadata,
    )
    left = max(80.0, config.tick_font_size * 4.5 + config.font_size * 2)
    right = config.size - max(40, config.tick_font_size * 3)
    top = max(36, config.value_font_size * 2)
    bottom_margin = max(80, config.font_size * 3 + config.tick_font_size * 2)
    if config.label_rotation:
        bottom_margin += config.tick_font_size * 7
    bottom = config.size - bottom_margin
    if right - left < 100 or bottom - top < 100:
        raise ConfigurationError(
            "Canvas is too small for these font sizes; increase size or reduce fonts."
        )
    width, height = right - left, bottom - top
    target_max = min(1.0, max(proportions) * 1.1)
    magnitude = 10 ** math.floor(math.log10(target_max / 6))
    tick_step = next(
        step * magnitude for step in (1, 2, 5, 10) if step * magnitude >= target_max / 6
    )
    ymax = min(1.0, math.ceil(target_max / tick_step) * tick_step)
    tick_count = round(ymax / tick_step)
    text_style = {"font_family": theme.font_family, "fill": theme.foreground}
    for index in range(tick_count + 1):
        fraction = index / tick_count
        y = bottom - height * fraction
        if config.show_grid:
            builder.line(left, y, right, y, stroke=theme.grid, class_="grid-line")
        builder.text(
            left - 10,
            y + config.tick_font_size * 0.35,
            f"{ymax * fraction * 100:.4g}%",
            text_anchor="end",
            font_size=config.tick_font_size,
            **text_style,
        )
    builder.line(left, top, left, bottom, stroke=theme.foreground)
    builder.line(left, bottom, right, bottom, stroke=theme.foreground)

    def position(value: float) -> float:
        # Scaling first also supports ranges whose unscaled difference overflows.
        scale = max(abs(edges[0]), abs(edges[-1]), 1.0)
        fraction = (value / scale - edges[0] / scale) / (edges[-1] / scale - edges[0] / scale)
        return left + width * fraction

    stride = max(1, math.ceil(len(counts) / config.max_tick_labels))
    for index, (label, count, proportion) in enumerate(
        zip(labels, counts, proportions, strict=True)
    ):
        start = position(edges[index]) if edges else left + width * index / len(counts)
        end = position(edges[index + 1]) if edges else left + width * (index + 1) / len(counts)
        center = (start + end) / 2
        bar_width = (end - start) * config.bar_width
        bar_height = height * proportion / ymax
        node = builder.rect(
            center - bar_width / 2,
            bottom - bar_height,
            bar_width,
            bar_height,
            fill=colors[index],
            class_="distribution-bar",
            data_label=label,
            data_count=count,
            data_proportion=proportion,
        )
        # Full labels and exact counts remain accessible even if tick text is shortened.
        ET.SubElement(
            node, "{http://www.w3.org/2000/svg}title"
        ).text = f"{label}: {count}/{valid} ({proportion:.2%})"
        if config.show_values:
            builder.text(
                center,
                bottom - bar_height - 6,
                f"{proportion:.1%}",
                text_anchor="middle",
                font_size=config.value_font_size,
                class_="bar-value",
                **text_style,
            )
        if not edges and index % stride == 0:
            tick_y = bottom + config.tick_font_size + 8
            display = label if len(label) <= 18 else label[:15] + "..."
            builder.text(
                center,
                tick_y,
                display,
                font_size=config.tick_font_size,
                text_anchor="end" if config.label_rotation else "middle",
                transform=f"rotate(-{config.label_rotation} {number(center)} {number(tick_y)})",
                class_="category-tick",
                **text_style,
            )
    if edges:
        tick_stride = max(1, math.ceil((len(edges) - 1) / config.max_tick_labels))
        tick_indices = sorted(set(range(0, len(edges), tick_stride)) | {len(edges) - 1})
        # Reserve the rightmost boundary, then omit labels that would overlap
        # at the requested font size. Every bin remains drawn and counted.
        visible_ticks = []
        next_left = math.inf
        angle = math.radians(config.label_rotation)
        for index in reversed(tick_indices):
            text_width = len(f"{edges[index]:.5g}") * config.tick_font_size * 0.65
            half_width = (
                text_width * math.cos(angle) + config.tick_font_size * math.sin(angle) + 12
            ) / 2
            x = position(edges[index])
            if x + half_width <= next_left:
                visible_ticks.append(index)
                next_left = x - half_width
        for index in reversed(visible_ticks):
            x, y = position(edges[index]), bottom + config.tick_font_size + 8
            builder.text(
                x,
                y,
                f"{edges[index]:.5g}",
                font_size=config.tick_font_size,
                text_anchor="end" if config.label_rotation else "middle",
                transform=f"rotate(-{config.label_rotation} {number(x)} {number(y)})",
                class_="bin-tick",
                **text_style,
            )
    builder.text(
        (left + right) / 2,
        config.size - config.font_size,
        config.x_label or ("Category" if task == "classification" else "Value"),
        text_anchor="middle",
        font_size=config.font_size,
        **text_style,
    )
    axis_x, axis_y = config.font_size * 1.5, (top + bottom) / 2
    builder.text(
        axis_x,
        axis_y,
        config.y_label,
        text_anchor="middle",
        font_size=config.font_size,
        transform=f"rotate(-90 {number(axis_x)} {number(axis_y)})",
        **text_style,
    )
    issues = []
    if missing:
        issues.append(
            Issue(
                "VIZ_MISSING_VALUES",
                "warning",
                "Blank samples were excluded from the denominator.",
                details={"missing_count": missing, "valid_count": valid},
            )
        )
    return builder.artifact(issues=issues)


def plot_distribution(
    values: Iterable[DistributionValue],
    *,
    task: DistributionTask,
    config: DistributionPlotConfig | None = None,
    verbose: bool = True,
) -> SVGArtifact:
    """Plot sample proportions; blanks are skipped and other invalid values raise.

    Classification uses first-seen label order unless ``categories`` is supplied.
    Bars use distinct colors, narrower widths and larger text by default;
    ``DistributionPlotConfig`` controls the palette, width and font sizes.
    Regression uses equal-width bins unless ``bin_edges`` is supplied. Heights
    are proportions, including for unequal bin widths, never probability density.
    Probability uses the same histogram with explicit [0, 1] value validation.
    """
    resolved = DistributionPlotConfig() if config is None else config
    if not isinstance(resolved, DistributionPlotConfig):
        raise TypeError("config must be DistributionPlotConfig.")
    _validate_task(task, resolved)
    labels, numbers, total, missing = _clean_values(values, task, resolved.max_values)
    edges: list[float] = []
    if task == "classification":
        counter = Counter(labels)
        labels = list(resolved.categories) if resolved.categories is not None else list(counter)
        if len(labels) > 1000:
            raise ConfigurationError("At most 1000 categories may be plotted.")
        if set(counter) - set(labels):
            raise InputFormatError("categories must include every observed label.")
        counts = [counter[label] for label in labels]
    else:
        edges, counts = _histogram(numbers, resolved)
        labels = [
            f"[{a:.8g}, {b:.8g}{']' if i == len(counts) - 1 else ')'}"
            for i, (a, b) in enumerate(pairwise(edges))
        ]
    return _render_distribution(
        labels, counts, edges, task=task, total=total, missing=missing, config=resolved
    )


def plot_csv_distributions(
    source: str | Path,
    *,
    columns: Sequence[str],
    task: DistributionTask,
    config: DistributionPlotConfig | None = None,
    encoding: str = "utf-8-sig",
    delimiter: str = ",",
    verbose: bool = True,
) -> dict[str, SVGArtifact]:
    """Read named CSV columns once, returning a separately normalized plot per column.

    CSV headers must be unique and rows must match the header width. UTF-8 BOM
    is accepted by default. Only empty/whitespace cells count as missing values.
    The caller chooses an output filename for each artifact with ``save_image``.
    """
    resolved = DistributionPlotConfig() if config is None else config
    if not isinstance(resolved, DistributionPlotConfig):
        raise TypeError("config must be DistributionPlotConfig.")
    _validate_task(task, resolved)
    if (
        isinstance(columns, str)
        or not columns
        or any(not isinstance(column, str) or not column for column in columns)
        or len(set(columns)) != len(columns)
    ):
        raise ConfigurationError("columns must be a non-empty sequence of unique column names.")
    if not isinstance(delimiter, str) or len(delimiter) != 1 or delimiter in "\r\n\0":
        raise ConfigurationError("delimiter must be one character, excluding newline and NUL.")
    selected: dict[str, list[DistributionValue]] = {column: [] for column in columns}
    with Path(source).open("r", encoding=encoding, newline="") as stream:
        reader = csv.reader(stream, delimiter=delimiter, strict=True)
        try:
            header = next(reader, None)
            if not header or len(set(header)) != len(header):
                raise InputFormatError("CSV must have a non-empty header with unique column names.")
            absent = [column for column in columns if column not in header]
            if absent:
                raise InputFormatError(
                    f"CSV columns not found: {absent}; available columns: {header}."
                )
            indices = {column: header.index(column) for column in columns}
            for row_index, row in enumerate(reader, start=1):
                if row_index > resolved.max_values:
                    raise ConfigurationError("CSV exceeds max_values; no rows were truncated.")
                if len(row) != len(header):
                    raise InputFormatError(
                        f"CSV record {row_index + 1} has {len(row)} fields; expected {len(header)}."
                    )
                for column, index in indices.items():
                    selected[column].append(row[index])
        except csv.Error as exc:
            raise InputFormatError(f"Malformed CSV at line {reader.line_num}: {exc}.") from exc
    plots: dict[str, SVGArtifact] = {}
    for column, values in selected.items():
        try:
            plots[column] = plot_distribution(
                values, task=task, config=replace(resolved, x_label=resolved.x_label or column)
            )
        except InputFormatError as exc:
            raise InputFormatError(f"Column {column!r}: {exc}") from exc
    return plots


__all__ = ["DistributionTask", "plot_csv_distributions", "plot_distribution"]

_register_public_verbose(globals())
