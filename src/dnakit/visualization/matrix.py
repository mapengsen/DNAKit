"""Deterministic SVG heatmaps for precomputed similarity matrices."""

from __future__ import annotations

import math

from dnakit.core.enums import IssueSeverity
from dnakit.core.issues import Issue
from dnakit.exceptions import ConfigurationError
from dnakit.similarity.results import SimilarityMatrixResult

from ._svg import SVGBuilder, interpolate_color, number, validate_xml_text
from .config import HeatmapConfig, LimitPolicy
from .results import SVGArtifact

_MAX_LABEL_LENGTH = 1_024
_DISPLAY_LABEL_LENGTH = 32


def _validate_labels(labels: tuple[str, ...]) -> None:
    for index, label in enumerate(labels):
        if len(label) > _MAX_LABEL_LENGTH:
            raise ConfigurationError(
                "A matrix label exceeds the visualization safety limit.",
                code="VISUALIZATION_SIZE_LIMIT",
                context={"label_index": index, "max_label_length": _MAX_LABEL_LENGTH},
            )
        try:
            validate_xml_text(label, name="matrix label")
        except ValueError as exc:
            raise ConfigurationError(
                "A matrix label contains a character forbidden by XML.",
                code="INVALID_VISUALIZATION_LABEL",
                context={"label_index": index},
            ) from exc


def _short_label(value: str) -> str:
    if len(value) <= _DISPLAY_LABEL_LENGTH:
        return value
    return value[: _DISPLAY_LABEL_LENGTH - 1] + "…"


def _stride_indices(length: int, maximum: int) -> tuple[int, ...]:
    if length <= maximum:
        return tuple(range(length))
    step = math.ceil(length / maximum)
    return tuple(range(0, length, step))


def _scale(values: tuple[tuple[float, ...], ...], config: HeatmapConfig) -> tuple[float, float]:
    flat = tuple(value for row in values for value in row)
    data_min = min(flat, default=0.0)
    data_max = max(flat, default=1.0)
    lower = data_min if config.value_min is None else float(config.value_min)
    upper = data_max if config.value_max is None else float(config.value_max)
    if lower >= upper:
        if config.value_min is not None and config.value_max is None:
            upper = lower + 1.0
        elif config.value_min is None and config.value_max is not None:
            lower = upper - 1.0
        elif lower == 0:
            upper = 1.0
        else:
            padding = max(abs(lower) * 0.05, 0.5)
            lower -= padding
            upper += padding
    return lower, upper


def plot_similarity_matrix(
    result: SimilarityMatrixResult,
    *,
    config: HeatmapConfig | None = None,
) -> SVGArtifact:
    """Render an existing matrix in its original order without mirroring it."""

    if not isinstance(result, SimilarityMatrixResult):
        raise TypeError("result must be SimilarityMatrixResult.")
    resolved = HeatmapConfig() if config is None else config
    if not isinstance(resolved, HeatmapConfig):
        raise TypeError("config must be HeatmapConfig or None.")
    _validate_labels(result.labels)

    size = result.item_count
    display_limit = min(resolved.max_rows, resolved.max_columns)
    if size > display_limit and resolved.limit_policy is LimitPolicy.ERROR:
        raise ConfigurationError(
            "Similarity heatmap exceeds its configured dimensions.",
            code="VISUALIZATION_SIZE_LIMIT",
            context={
                "item_count": size,
                "max_columns": resolved.max_columns,
                "max_rows": resolved.max_rows,
            },
            hint="Raise the limits or explicitly choose limit_policy='stride'.",
        )
    indices = _stride_indices(size, display_limit)
    strided = len(indices) < size
    lower, upper = _scale(result.values, resolved)
    issues: list[Issue] = []
    if size == 0:
        issues.append(
            Issue(
                "VIZ_EMPTY_MATRIX",
                IssueSeverity.INFO,
                "The similarity result contains an empty matrix.",
            )
        )
    if strided:
        issues.append(
            Issue(
                "VIZ_MATRIX_STRIDED",
                IssueSeverity.WARNING,
                "Rows and columns were explicitly display-sampled with the same stride.",
                details={"input_items": size, "displayed_items": len(indices)},
            )
        )
    clipped_count = sum(value < lower or value > upper for row in result.values for value in row)
    if clipped_count:
        issues.append(
            Issue(
                "VIZ_MATRIX_VALUES_CLIPPED",
                IssueSeverity.WARNING,
                "Values outside the configured color range were clipped.",
                details={"count": clipped_count},
            )
        )

    labels = tuple(result.labels[index] for index in indices)
    title = resolved.title or f"{result.name} ({result.method})"
    label_gutter = 170 if resolved.show_labels and labels else 24
    top_gutter = 190 if resolved.show_labels and labels else 60
    legend_width = 80
    matrix_size = max(1, len(indices)) * resolved.cell_size
    width = label_gutter + matrix_size + legend_width + 28
    height = top_gutter + matrix_size + 42
    metadata = {
        "display_indices": indices,
        "displayed_items": len(indices),
        "input_items": size,
        "matrix_method": result.method,
        "strided": strided,
        "symmetric": result.symmetric,
        "value_kind": result.value_kind,
        "value_max": upper,
        "value_min": lower,
    }
    builder = SVGBuilder(
        width,
        height,
        title=title,
        description="Heatmap of a precomputed DNA similarity or distance matrix.",
        kind="similarity-matrix",
        theme=resolved.theme,
        metadata=metadata,
    )
    builder.text(
        width / 2,
        28,
        title,
        fill=resolved.theme.foreground,
        font_family=resolved.theme.font_family,
        font_size=16,
        font_weight="bold",
        text_anchor="middle",
        class_="plot-title",
    )

    if not indices:
        builder.text(
            width / 2,
            height / 2,
            "Empty matrix",
            fill=resolved.theme.muted,
            font_family=resolved.theme.font_family,
            font_size=14,
            text_anchor="middle",
            class_="empty-state",
        )
        return builder.artifact(issues=issues)

    for display_row, source_row in enumerate(indices):
        y = top_gutter + display_row * resolved.cell_size
        if resolved.show_labels:
            builder.text(
                label_gutter - 8,
                y + resolved.cell_size * 0.72,
                _short_label(result.labels[source_row]),
                fill=resolved.theme.foreground,
                font_family=resolved.theme.font_family,
                font_size=10,
                text_anchor="end",
                class_="row-label",
                data_label=result.labels[source_row],
                data_source_index=source_row,
            )
        for display_column, source_column in enumerate(indices):
            value = result.values[source_row][source_column]
            fraction = (min(upper, max(lower, value)) - lower) / (upper - lower)
            cell_x = label_gutter + display_column * resolved.cell_size
            builder.rect(
                cell_x,
                y,
                resolved.cell_size,
                resolved.cell_size,
                fill=interpolate_color(resolved.low_color, resolved.high_color, fraction),
                stroke=resolved.theme.background,
                stroke_width=0.5,
                class_="heatmap-cell",
                data_column=source_column,
                data_row=source_row,
                data_value=number(value),
            )
            if resolved.show_values and resolved.cell_size >= 18:
                builder.text(
                    cell_x + resolved.cell_size / 2,
                    y + resolved.cell_size * 0.7,
                    number(value),
                    fill=resolved.theme.foreground,
                    font_family=resolved.theme.font_family,
                    font_size=max(7, min(10, resolved.cell_size - 9)),
                    text_anchor="middle",
                    class_="cell-value",
                )
    if resolved.show_labels:
        for display_column, source_column in enumerate(indices):
            label_x = label_gutter + (display_column + 0.5) * resolved.cell_size
            label_y = top_gutter - 8
            builder.text(
                label_x,
                label_y,
                _short_label(result.labels[source_column]),
                fill=resolved.theme.foreground,
                font_family=resolved.theme.font_family,
                font_size=10,
                text_anchor="start",
                transform=f"rotate(-45 {number(label_x)} {number(label_y)})",
                class_="column-label",
                data_label=result.labels[source_column],
                data_source_index=source_column,
            )

    legend_x = label_gutter + matrix_size + 24
    legend_steps = 20
    legend_height = min(200, matrix_size)
    for step in range(legend_steps):
        fraction = step / (legend_steps - 1)
        legend_y = top_gutter + (1 - fraction) * legend_height
        builder.rect(
            legend_x,
            legend_y,
            14,
            legend_height / legend_steps + 1,
            fill=interpolate_color(resolved.low_color, resolved.high_color, fraction),
            class_="legend-step",
        )
    builder.text(
        legend_x + 20,
        top_gutter + 8,
        number(upper),
        fill=resolved.theme.muted,
        font_family=resolved.theme.font_family,
        font_size=9,
        class_="legend-maximum",
    )
    builder.text(
        legend_x + 20,
        top_gutter + legend_height,
        number(lower),
        fill=resolved.theme.muted,
        font_family=resolved.theme.font_family,
        font_size=9,
        class_="legend-minimum",
    )
    return builder.artifact(issues=issues)


__all__ = ["plot_similarity_matrix"]
