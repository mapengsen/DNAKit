"""Minimal dependency-free SVG scatter plots and heatmaps."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path

from dnakit.exceptions import ConfigurationError

# Okabe-Ito colour-blind-safe qualitative palette, then extra distinct colours.
PALETTE = (
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000",
    "#999999", "#882255", "#44AA99", "#117733", "#332288", "#AA4499", "#DDCC77", "#88CCEE",
)
WIDTH, HEIGHT, MARGIN = 640, 480, 40


def _scale(values: Sequence[float], size: float) -> list[float]:
    low, high = min(values), max(values)
    span = high - low
    if span == 0:
        return [size / 2.0] * len(values)
    return [(v - low) / span * size for v in values]


def _sequential_colour(fraction: float) -> str:
    """White-to-dark-blue ramp for ``fraction`` in [0, 1]."""

    fraction = min(1.0, max(0.0, fraction))
    red = round(247 + (8 - 247) * fraction)
    green = round(251 + (48 - 251) * fraction)
    blue = round(255 + (107 - 255) * fraction)
    return f"#{red:02x}{green:02x}{blue:02x}"


def scatter_svg(
    xs: Sequence[float],
    ys: Sequence[float],
    *,
    names: Sequence[str],
    labels: Sequence[str] | None = None,
    values: Sequence[float] | None = None,
    title: str = "",
    invert_y: bool = False,
    radius: float = 3.0,
) -> str:
    """Scatter plot with categorical colours (``labels``) or a continuous ramp (``values``)."""

    plot_w, plot_h = WIDTH - 2 * MARGIN, HEIGHT - 2 * MARGIN
    px = _scale(xs, plot_w)
    py = _scale(ys, plot_h)
    legend: list[str] = []
    colours: list[str]
    if labels is not None:
        order = sorted(set(labels))
        lookup = {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(order)}
        colours = [lookup[label] for label in labels]
        legend = [
            f'<rect x="{WIDTH - MARGIN + 6}" y="{MARGIN + 16 * i}" width="10" height="10" fill="{lookup[label]}"/>'
            f'<text x="{WIDTH - MARGIN + 20}" y="{MARGIN + 16 * i + 9}" font-size="10">{escape(label)}</text>'
            for i, label in enumerate(order[:20])
        ]
    elif values is not None:
        low, high = min(values), max(values)
        span = (high - low) or 1.0
        colours = [_sequential_colour((v - low) / span) for v in values]
    else:
        colours = [PALETTE[0]] * len(xs)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH + 80}" height="{HEIGHT}" viewBox="0 0 {WIDTH + 80} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{MARGIN}" y="22" font-size="14">{escape(title)}</text>',
        f'<rect x="{MARGIN}" y="{MARGIN}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#888888"/>',
    ]
    for name, x, y, colour in zip(names, px, py, colours):
        cy = MARGIN + (y if invert_y else plot_h - y)
        parts.append(
            f'<circle cx="{MARGIN + x:.2f}" cy="{cy:.2f}" r="{radius}" fill="{colour}"><title>{escape(name)}</title></circle>'
        )
    parts.extend(legend)
    parts.append("</svg>")
    return "".join(parts)


def heatmap_svg(
    matrix: Sequence[Sequence[float]],
    row_names: Sequence[str],
    col_names: Sequence[str],
    *,
    title: str = "",
) -> str:
    """Heatmap with one rectangle per value, colour scaled between the matrix minimum and maximum."""

    flat = [v for row in matrix for v in row]
    low, high = min(flat), max(flat)
    span = (high - low) or 1.0
    label_w = 90
    cell_w = max(2.0, min(20.0, (WIDTH - label_w) / max(1, len(col_names))))
    cell_h = max(6.0, min(20.0, (HEIGHT - 60) / max(1, len(row_names))))
    width = label_w + cell_w * len(col_names) + 20
    height = 40 + cell_h * len(row_names) + 20
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="8" y="20" font-size="14">{escape(title)}</text>',
    ]
    for i, (name, row) in enumerate(zip(row_names, matrix)):
        y = 40 + i * cell_h
        parts.append(f'<text x="{label_w - 4}" y="{y + cell_h * 0.75:.1f}" font-size="{min(10, cell_h):.0f}" text-anchor="end">{escape(name)}</text>')
        for j, value in enumerate(row):
            parts.append(
                f'<rect x="{label_w + j * cell_w:.2f}" y="{y:.2f}" width="{cell_w:.2f}" height="{cell_h:.2f}" '
                f'fill="{_sequential_colour((value - low) / span)}"><title>{escape(name)} / {escape(col_names[j])}: {value:.4g}</title></rect>'
            )
    parts.append("</svg>")
    return "".join(parts)


def write_svg(svg: str, options: Mapping[str, object]) -> str | None:
    """Write the SVG only when ``options['output_path']`` is set; never overwrite an existing file."""

    target = options.get("output_path")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or not str(target).strip():
        raise ConfigurationError("options['output_path'] must be a path.", code="INVALID_OPTION", hint="Pass a file path ending in .svg.")
    path = Path(target).expanduser()
    if path.suffix.lower() != ".svg":
        raise ConfigurationError("output_path must end with .svg.", code="INVALID_OPTION", hint="Use a .svg file name.")
    if path.exists():
        raise ConfigurationError(f"Output file already exists: {path}", code="OUTPUT_EXISTS", hint="Choose a new output_path.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}", code="OUTPUT_DIR", hint="Create the directory first.")
    path.write_text(svg, encoding="utf-8")
    return str(path)


def finite(values: Sequence[float]) -> bool:
    return all(math.isfinite(v) for v in values)
