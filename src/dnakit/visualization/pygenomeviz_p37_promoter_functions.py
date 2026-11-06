"""Promoter glyph rendering with pyGenomeViz 0.4 and 1.x (pack P37)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from importlib.metadata import version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def draw_genome_promoters(
    promoters: Sequence[Mapping[str, Any]],
    output_path: str,
    *,
    segment_size: int,
    segment_name: str = "genome",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw promoter arrows using FeatureSegment.add_promoter and return image path plus glyph data.

    Each promoter has ``position`` (0-based), optional ``strand`` (+1/-1),
    ``length`` (relative to segment, default 0.02), and ``color`` (default black).
    Uses FeatureSegment.add_promoter on pyGenomeViz 1.x. On 0.4 the same
    right-angle promoter geometry is drawn on the feature track's matplotlib
    axes, preserving compatibility with existing DNAKit genome plots.
    """
    pgv = require("pygenomeviz", "bio-pygenomeviz")
    modern = int(version("pygenomeviz").split(".")[0]) >= 1
    if isinstance(segment_size, bool) or not isinstance(segment_size, int) or segment_size < 1:
        raise ConfigurationError("segment_size must be a positive integer.")
    if isinstance(promoters, (str, bytes)) or not isinstance(promoters, Sequence) or not promoters:
        raise ConfigurationError("promoters must be a non-empty list.")
    target = Path(output_path).expanduser().resolve()
    if (
        target.suffix.lower() not in (".png", ".svg", ".pdf")
        or target.exists()
        or not target.parent.is_dir()
    ):
        raise ConfigurationError(
            "output_path must be a new PNG, SVG or PDF file in an existing directory."
        )
    graph = pgv.GenomeViz()
    track = graph.add_feature_track(segment_name, segment_size)
    segment = track.get_segment() if modern else None
    rows = []
    for item in promoters:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Each promoter must be a mapping.")
        position, strand, length = (
            item.get("position"),
            item.get("strand", 1),
            item.get("length", 0.02),
        )
        if (
            isinstance(position, bool)
            or not isinstance(position, (int, float))
            or not 0 <= position <= segment_size
        ):
            raise ConfigurationError("position must lie inside the segment.")
        if isinstance(strand, bool) or strand not in (-1, 1):
            raise ConfigurationError("strand must be +1 or -1.")
        if isinstance(length, bool) or not isinstance(length, (int, float)) or not 0 < length < 1:
            raise ConfigurationError("length must be a fraction in (0, 1).")
        endpoint = float(position) + strand * float(length) * segment_size
        if not 0 <= endpoint <= segment_size:
            raise ConfigurationError("The promoter arrow extends outside the segment.")
        color = str(item.get("color", "black"))
        if modern:
            segment.add_promoter(float(position), float(length), strand=strand, color=color)
        rows.append(
            {
                "position": position,
                "strand": strand,
                "length": length,
                "color": color,
                "arrow_start": [float(position), 0.0],
                "arrow_end": [endpoint, 1.0],
            }
        )
    pyplot = require("matplotlib.pyplot", "bio-pygenomeviz")
    try:
        figure = graph.plotfig()
    except ValueError as exc:
        raise ConfigurationError(f"A promoter glyph extends outside the segment: {exc}") from exc
    try:
        if not modern:
            for row in rows:
                track.ax.annotate(
                    "",
                    xy=row["arrow_end"],
                    xytext=row["arrow_start"],
                    arrowprops={
                        "arrowstyle": "-|>,head_length=1,head_width=0.5",
                        "connectionstyle": "angle",
                        "shrinkA": 0,
                        "shrinkB": 0,
                        "color": row["color"],
                    },
                    annotation_clip=False,
                )
        figure.savefig(target, bbox_inches="tight")
    finally:
        pyplot.close(figure)
    return {
        "status": "completed",
        "image_path": str(target),
        "promoters": rows,
        "segment_size": segment_size,
        "provenance": {
            "pygenomeviz": version("pygenomeviz"),
            "method": "FeatureSegment.add_promoter"
            if modern
            else "feature-track matplotlib promoter geometry",
        },
    }


__all__ = ["draw_genome_promoters"]
_register_public_verbose(globals())
