"""Validated configuration for dependency-free SVG visualizations."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import cast

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

_COLOR_PATTERN = re.compile(r"^(?:#[0-9A-Fa-f]{3,8}|[A-Za-z][A-Za-z0-9_-]*)$")
_FONT_PATTERN = re.compile(r"^[A-Za-z0-9 _,-]+$")
_DNA_DISPLAY_SYMBOLS = frozenset("ACGTRYSWKMBDHVN")


class LimitPolicy(str, Enum):
    """Explicit behavior when a plot exceeds its configured display limit."""

    ERROR = "error"
    TRUNCATE = "truncate"
    STRIDE = "stride"


def _positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code="INVALID_VISUALIZATION_CONFIG",
            context={name: value},
        )


def _non_negative_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(
            f"{name} must be a non-negative integer.",
            code="INVALID_VISUALIZATION_CONFIG",
            context={name: value},
        )


def _finite_optional(value: float | None, name: str) -> None:
    if value is not None and (
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
    ):
        raise ConfigurationError(
            f"{name} must be a finite number or None.",
            code="INVALID_VISUALIZATION_CONFIG",
            context={name: value},
        )


def _color(value: str, name: str) -> None:
    if not isinstance(value, str) or not _COLOR_PATTERN.fullmatch(value):
        raise ConfigurationError(
            f"{name} must be a simple CSS color name or hexadecimal color.",
            code="INVALID_VISUALIZATION_COLOR",
            context={name: value},
        )


def _title(value: str | None, name: str = "title") -> None:
    if value is not None and (not isinstance(value, str) or not value.strip()):
        raise ConfigurationError(
            f"{name} must be a non-empty string or None.",
            code="INVALID_VISUALIZATION_CONFIG",
        )


def _is_xml_character(value: str) -> bool:
    codepoint = ord(value)
    return (
        0x20 <= codepoint <= 0xD7FF
        or 0xE000 <= codepoint <= 0xFFFD
        or 0x10000 <= codepoint <= 0x10FFFF
    )


def _symbol_map(value: Mapping[str, str]) -> FrozenDict:
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            "symbol_map must be a mapping.",
            code="INVALID_VISUALIZATION_SYMBOL_MAP",
            context={"value_type": type(value).__name__},
        )
    resolved: dict[str, str] = {}
    for source_symbol, display_symbol in value.items():
        if (
            not isinstance(source_symbol, str)
            or len(source_symbol) != 1
            or source_symbol not in _DNA_DISPLAY_SYMBOLS
        ):
            raise ConfigurationError(
                "symbol_map keys must be uppercase DNA or IUPAC symbols.",
                code="INVALID_VISUALIZATION_SYMBOL_MAP",
                context={"source_symbol": repr(source_symbol)},
            )
        if (
            not isinstance(display_symbol, str)
            or len(display_symbol) != 1
            or display_symbol.isspace()
            or not display_symbol.isprintable()
            or not _is_xml_character(display_symbol)
        ):
            raise ConfigurationError(
                "symbol_map values must each be one visible XML character.",
                code="INVALID_VISUALIZATION_SYMBOL_MAP",
                context={"source_symbol": source_symbol},
            )
        resolved[source_symbol] = display_symbol
    return FrozenDict({symbol: resolved[symbol] for symbol in sorted(resolved)})


def _coerce_policy(value: LimitPolicy | str, *, allowed: set[LimitPolicy]) -> LimitPolicy:
    try:
        policy = value if isinstance(value, LimitPolicy) else LimitPolicy(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Unknown visualization limit policy.",
            code="INVALID_VISUALIZATION_LIMIT_POLICY",
            context={"policy": value},
        ) from exc
    if policy not in allowed:
        raise ConfigurationError(
            "The selected limit policy is not valid for this plot type.",
            code="INVALID_VISUALIZATION_LIMIT_POLICY",
            context={"policy": policy.value, "allowed": sorted(item.value for item in allowed)},
        )
    return policy


@dataclass(frozen=True, slots=True)
class SVGTheme:
    """Shared colors and typography embedded directly in an SVG artifact."""

    background: str = "#ffffff"
    foreground: str = "#1f2937"
    muted: str = "#6b7280"
    grid: str = "#d1d5db"
    accent: str = "#2563eb"
    missing: str = "#e5e7eb"
    font_family: str = "monospace"

    def __post_init__(self) -> None:
        for name in ("background", "foreground", "muted", "grid", "accent", "missing"):
            _color(getattr(self, name), name)
        if not isinstance(self.font_family, str) or not _FONT_PATTERN.fullmatch(self.font_family):
            raise ConfigurationError(
                "font_family contains unsupported characters.",
                code="INVALID_VISUALIZATION_FONT",
                context={"font_family": self.font_family},
            )


@dataclass(frozen=True, slots=True)
class SequencePlotConfig:
    """Layout and safety limits for the gap-aware sequence text plot.

    ``column_spacing`` adds horizontal pixels between adjacent text columns,
    while ``line_spacing`` adds vertical pixels between adjacent sequence rows.
    ``symbol_map`` replaces rendered DNA/IUPAC symbols without changing sequence
    coordinates, highlights, complements, or base colors. ``title`` supplies
    SVG accessibility text and is not drawn inside the square plot canvas.
    """

    bases_per_line: int = 60
    max_symbols: int = 10_000
    max_gaps: int = 1_000
    limit_policy: LimitPolicy = LimitPolicy.ERROR
    start_coordinate: int = 0
    show_coordinates: bool = True
    show_complement: bool = False
    font_size: int = 14
    cell_width: int = 12
    line_height: int = 62
    column_spacing: int = 0
    line_spacing: int = 0
    margin: int = 24
    title: str | None = None
    theme: SVGTheme = field(default_factory=SVGTheme)
    symbol_map: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in (
            "bases_per_line",
            "max_symbols",
            "max_gaps",
            "font_size",
            "cell_width",
            "line_height",
            "margin",
        ):
            _positive_int(getattr(self, name), name)
        for name in ("column_spacing", "line_spacing"):
            _non_negative_int(getattr(self, name), name)
        _non_negative_int(self.start_coordinate, "start_coordinate")
        for name in ("show_coordinates", "show_complement"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(
                    f"{name} must be a boolean.", code="INVALID_VISUALIZATION_CONFIG"
                )
        _title(self.title)
        if not isinstance(self.theme, SVGTheme):
            raise ConfigurationError("theme must be SVGTheme.")
        object.__setattr__(
            self,
            "symbol_map",
            cast(Mapping[str, str], _symbol_map(self.symbol_map)),
        )
        object.__setattr__(
            self,
            "limit_policy",
            _coerce_policy(self.limit_policy, allowed={LimitPolicy.ERROR, LimitPolicy.TRUNCATE}),
        )


@dataclass(frozen=True, slots=True)
class HeatmapConfig:
    """Limits and color scale for square matrix and fingerprint heatmaps.

    ``title`` supplies SVG accessibility text and is not drawn inside the plot.
    """

    cell_size: int = 18
    max_rows: int = 200
    max_columns: int = 256
    limit_policy: LimitPolicy = LimitPolicy.ERROR
    value_min: float | None = None
    value_max: float | None = None
    low_color: str = "#f7fbff"
    high_color: str = "#08306b"
    show_values: bool = False
    show_labels: bool = True
    title: str | None = None
    theme: SVGTheme = field(default_factory=SVGTheme)

    def __post_init__(self) -> None:
        for name in ("cell_size", "max_rows", "max_columns"):
            _positive_int(getattr(self, name), name)
        _finite_optional(self.value_min, "value_min")
        _finite_optional(self.value_max, "value_max")
        if (
            self.value_min is not None
            and self.value_max is not None
            and self.value_min >= self.value_max
        ):
            raise ConfigurationError(
                "value_min must be smaller than value_max.",
                code="INVALID_VISUALIZATION_CONFIG",
            )
        for name in ("low_color", "high_color"):
            _color(getattr(self, name), name)
            if not re.fullmatch(r"#(?:[0-9A-Fa-f]{3}|[0-9A-Fa-f]{6})", getattr(self, name)):
                raise ConfigurationError(
                    f"{name} must use #RGB or #RRGGBB for interpolation.",
                    code="INVALID_VISUALIZATION_COLOR",
                )
        for name in ("show_values", "show_labels"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"{name} must be a boolean.")
        _title(self.title)
        if not isinstance(self.theme, SVGTheme):
            raise ConfigurationError("theme must be SVGTheme.")
        object.__setattr__(
            self,
            "limit_policy",
            _coerce_policy(self.limit_policy, allowed={LimitPolicy.ERROR, LimitPolicy.STRIDE}),
        )


@dataclass(frozen=True, slots=True)
class SaveConfig:
    """Filesystem policy for deterministic SVG persistence."""

    overwrite: bool = False
    create_parents: bool = False

    def __post_init__(self) -> None:
        for name in ("overwrite", "create_parents"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(
                    f"{name} must be a boolean.", code="INVALID_VISUALIZATION_CONFIG"
                )


@dataclass(frozen=True, slots=True)
class ImageExportConfig:
    """Safe raster/PDF export settings for the optional graphics backend.

    When both ``width`` and ``height`` are supplied, they must be equal so
    exports from DNAKit's square plot artifacts remain square.
    """

    dpi: int = 600
    width: int | None = None
    height: int | None = None
    transparent: bool = False
    max_output_pixels: int = 100_000_000
    overwrite: bool = False
    create_parents: bool = False

    def __post_init__(self) -> None:
        _positive_int(self.dpi, "dpi")
        _positive_int(self.max_output_pixels, "max_output_pixels")
        for name in ("width", "height"):
            value = getattr(self, name)
            if value is not None:
                _positive_int(value, name)
        if self.width is not None and self.height is not None and self.width != self.height:
            raise ConfigurationError(
                "width and height must be equal for square image output.",
                code="INVALID_VISUALIZATION_CONFIG",
            )
        for name in ("transparent", "overwrite", "create_parents"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(
                    f"{name} must be a boolean.", code="INVALID_VISUALIZATION_CONFIG"
                )


__all__ = [
    "HeatmapConfig",
    "ImageExportConfig",
    "LimitPolicy",
    "SVGTheme",
    "SaveConfig",
    "SequencePlotConfig",
]
