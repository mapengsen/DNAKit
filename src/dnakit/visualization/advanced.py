"""Dependency-free SVG renderers for advanced DNAKit result relationships."""

from __future__ import annotations

import math

from dnakit.alignment import AlignmentResult
from dnakit.core import (
    CompoundLocation,
    DNARecord,
    DNASequence,
    Interval,
    Topology,
    UnresolvedLocation,
)
from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.exceptions import ConfigurationError

from ._svg import SVGBuilder, number
from .config import SVGTheme
from .results import SVGArtifact

_MAX_CANVAS_DIMENSION = 1_000_000
_MAX_MAP_PARTS = 100_000
_MAX_ALIGNMENT_COLUMNS = 100_000


def _bounded_integer(
    value: object,
    name: str,
    *,
    minimum: int = 1,
    maximum: int,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _validate_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > 10_000:
        raise ConfigurationError(
            f"{name} exceeds the visualization text limit.",
            code="VISUALIZATION_TEXT_LIMIT",
            context={"name": name, "max_characters": 10_000},
        )
    return value


def _resolve_theme(theme: SVGTheme | None) -> SVGTheme:
    if theme is None:
        return SVGTheme()
    if not isinstance(theme, SVGTheme):
        raise ConfigurationError("theme must be SVGTheme or None.")
    return theme


def _sequence_and_features(
    value: DNA | DNASequence | DNARecord,
) -> tuple[DNASequence, tuple[object, ...]]:
    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return sequence, () if record is None else record.features
    raise ConfigurationError("value must be DNASequence or DNARecord.")


def _parts(location: object) -> tuple[Interval, ...]:
    if isinstance(location, Interval):
        return (location,)
    if isinstance(location, CompoundLocation):
        return location.parts
    if isinstance(location, UnresolvedLocation):
        return ()
    return ()


def _feature_parts(
    features: tuple[object, ...], max_features: int
) -> tuple[tuple[int, Interval], ...]:
    if len(features) > max_features:
        raise ConfigurationError("Feature count exceeds max_features.")
    output: list[tuple[int, Interval]] = []
    for index, feature in enumerate(features):
        for part in _parts(getattr(feature, "location", None)):
            if len(output) >= max_features:
                raise ConfigurationError(
                    "Rendered feature-part count exceeds max_features.",
                    code="VISUALIZATION_FEATURE_PART_LIMIT",
                )
            output.append((index, part))
    return tuple(output)


def plot_linear_map(
    value: DNA | DNASequence | DNARecord,
    *,
    width: int = 900,
    height: int = 220,
    max_features: int = 10_000,
    title: str = "Linear DNA map",
    theme: SVGTheme | None = None,
) -> SVGArtifact:
    """Draw a resolved linear coordinate axis and feature intervals."""

    _bounded_integer(width, "width", minimum=120, maximum=_MAX_CANVAS_DIMENSION)
    _bounded_integer(height, "height", minimum=100, maximum=_MAX_CANVAS_DIMENSION)
    _bounded_integer(max_features, "max_features", maximum=_MAX_MAP_PARTS)
    _validate_text(title, "title")
    sequence, features = _sequence_and_features(value)
    if sequence.topology is not Topology.LINEAR or sequence.coordinate_span is None:
        raise ConfigurationError("Linear map requires a resolved linear sequence.")
    parts = _feature_parts(features, max_features)
    resolved_theme = _resolve_theme(theme)
    span = sequence.coordinate_span
    builder = SVGBuilder(
        width,
        height,
        title=title,
        description="DNAKit linear DNA coordinate and feature map",
        kind="linear-map",
        theme=resolved_theme,
        metadata={
            "span": span,
            "feature_count": len(features),
            "rendered_part_count": len(parts),
            "unresolved_feature_count": sum(
                isinstance(getattr(feature, "location", None), UnresolvedLocation)
                for feature in features
            ),
        },
    )
    left, right, axis_y = 60, width - 40, height // 2
    builder.line(left, axis_y, right, axis_y, stroke=resolved_theme.foreground, stroke_width=3)
    scale = (right - left) / max(1, span)
    builder.text(left, axis_y + 28, 0, fill=resolved_theme.muted)
    builder.text(right, axis_y + 28, span, fill=resolved_theme.muted, text_anchor="end")
    for index, part in parts:
        builder.rect(
            left + part.start * scale,
            axis_y - 18 - 18 * (index % 3),
            max(1.0, len(part) * scale),
            12,
            fill=resolved_theme.accent,
            class_="feature",
            data_feature_index=index,
            data_start=part.start,
            data_end=part.end,
        )
    return builder.artifact()


def plot_circular_map(
    value: DNA | DNASequence | DNARecord,
    *,
    size: int = 520,
    max_features: int = 10_000,
    title: str = "Circular DNA map",
    theme: SVGTheme | None = None,
) -> SVGArtifact:
    """Draw a resolved circular molecule and angular feature arcs."""

    _bounded_integer(size, "size", minimum=100, maximum=_MAX_CANVAS_DIMENSION)
    _bounded_integer(max_features, "max_features", maximum=_MAX_MAP_PARTS)
    _validate_text(title, "title")
    sequence, features = _sequence_and_features(value)
    if sequence.topology is not Topology.CIRCULAR or sequence.coordinate_span is None:
        raise ConfigurationError("Circular map requires a resolved circular sequence.")
    parts = _feature_parts(features, max_features)
    resolved_theme = _resolve_theme(theme)
    span = sequence.coordinate_span
    center, radius = size / 2, size * 0.34
    builder = SVGBuilder(
        size,
        size,
        title=title,
        description="DNAKit circular DNA feature map",
        kind="circular-map",
        theme=resolved_theme,
        metadata={
            "span": span,
            "feature_count": len(features),
            "rendered_part_count": len(parts),
            "unresolved_feature_count": sum(
                isinstance(getattr(feature, "location", None), UnresolvedLocation)
                for feature in features
            ),
        },
    )
    builder.element(
        "circle",
        cx=center,
        cy=center,
        r=radius,
        fill="none",
        stroke=resolved_theme.foreground,
        stroke_width=4,
        class_="backbone",
    )
    for index, part in parts:
        start_angle = 2 * math.pi * part.start / span - math.pi / 2
        x1 = center + radius * math.cos(start_angle)
        y1 = center + radius * math.sin(start_angle)
        if len(part) == span:
            builder.element(
                "circle",
                cx=center,
                cy=center,
                r=radius,
                fill="none",
                stroke=resolved_theme.accent,
                stroke_width=12,
                class_="feature",
                data_feature_index=index,
                data_start=part.start,
                data_end=part.end,
            )
            continue
        if len(part) == 0:
            builder.element(
                "circle",
                cx=x1,
                cy=y1,
                r=6,
                fill=resolved_theme.accent,
                class_="feature",
                data_feature_index=index,
                data_start=part.start,
                data_end=part.end,
            )
            continue
        end_angle = 2 * math.pi * part.end / span - math.pi / 2
        x2, y2 = center + radius * math.cos(end_angle), center + radius * math.sin(end_angle)
        large = int(len(part) > span / 2)
        path = (
            f"M {number(x1)} {number(y1)} A {number(radius)} {number(radius)} "
            f"0 {large} 1 {number(x2)} {number(y2)}"
        )
        builder.element(
            "path",
            d=path,
            fill="none",
            stroke=resolved_theme.accent,
            stroke_width=12,
            class_="feature",
            data_feature_index=index,
            data_start=part.start,
            data_end=part.end,
        )
    builder.text(center, center, f"{span} bp", text_anchor="middle", fill=resolved_theme.foreground)
    return builder.artifact()


def plot_alignment(
    result: AlignmentResult,
    *,
    columns_per_line: int = 80,
    max_columns: int = 20_000,
    theme: SVGTheme | None = None,
) -> SVGArtifact:
    """Render aligned strings in fixed-width blocks with a match line."""

    if not isinstance(result, AlignmentResult):
        raise ConfigurationError("result must be AlignmentResult.")
    _bounded_integer(
        max_columns,
        "max_columns",
        maximum=_MAX_ALIGNMENT_COLUMNS,
    )
    _bounded_integer(
        columns_per_line,
        "columns_per_line",
        maximum=max_columns,
    )
    if result.alignment_length > max_columns:
        raise ConfigurationError("Alignment exceeds max_columns.")
    resolved_theme = _resolve_theme(theme)
    blocks = max(1, math.ceil(result.alignment_length / columns_per_line))
    width, height = max(640, columns_per_line * 10 + 100), 70 + blocks * 72
    builder = SVGBuilder(
        width,
        height,
        title="DNA alignment",
        description="DNAKit pairwise alignment",
        kind="alignment",
        theme=resolved_theme,
        metadata={
            "mode": result.method,
            "score": result.score,
            "columns": result.alignment_length,
        },
    )
    for block in range(blocks):
        start = block * columns_per_line
        end = min(result.alignment_length, start + columns_per_line)
        query = result.aligned_query[start:end]
        target = result.aligned_target[start:end]
        matches = "".join(
            "|" if left == right and left != "-" else " "
            for left, right in zip(query, target, strict=True)
        )
        y = 50 + block * 72
        builder.text(24, y, query, class_="aligned-query", fill=resolved_theme.foreground)
        builder.text(24, y + 20, matches, class_="match-line", fill=resolved_theme.accent)
        builder.text(24, y + 40, target, class_="aligned-target", fill=resolved_theme.foreground)
    return builder.artifact()


__all__ = [
    "plot_alignment",
    "plot_circular_map",
    "plot_linear_map",
]
