"""Dependency-free gap-aware DNA sequence text visualization."""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

from dnakit.core.enums import IssueSeverity, Topology
from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.core.gap import Gap
from dnakit.core.issues import Issue
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError

from ._svg import SVGBuilder, validate_xml_text
from .config import LimitPolicy, SequencePlotConfig, _color
from .results import SVGArtifact

_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


@dataclass(frozen=True, slots=True)
class Highlight:
    """A zero-based, half-open highlight in nucleotide-symbol coordinates."""

    start: int
    end: int
    label: str | None = None
    color: str = "#fde68a"
    foreground: str | None = None
    priority: int = 0
    opacity: float = 0.7

    def __post_init__(self) -> None:
        for name, value in (("start", self.start), ("end", self.end)):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(
                    f"Highlight {name} must be a non-negative integer.",
                    code="INVALID_HIGHLIGHT",
                )
        if self.end <= self.start:
            raise ConfigurationError(
                "Highlight end must be greater than start.", code="INVALID_HIGHLIGHT"
            )
        if self.label is not None:
            if not isinstance(self.label, str) or not self.label.strip():
                raise ConfigurationError(
                    "Highlight label must be non-empty or None.", code="INVALID_HIGHLIGHT"
                )
            try:
                validate_xml_text(self.label, name="highlight label")
            except ValueError as exc:
                raise ConfigurationError(
                    "Highlight label contains a character forbidden by XML.",
                    code="INVALID_HIGHLIGHT",
                ) from exc
        _color(self.color, "highlight color")
        if self.foreground is not None:
            _color(self.foreground, "highlight foreground")
        if isinstance(self.priority, bool) or not isinstance(self.priority, int):
            raise ConfigurationError(
                "Highlight priority must be an integer.", code="INVALID_HIGHLIGHT"
            )
        if (
            isinstance(self.opacity, bool)
            or not isinstance(self.opacity, (int, float))
            or not math.isfinite(self.opacity)
            or not 0.0 <= self.opacity <= 1.0
        ):
            raise ConfigurationError(
                "Highlight opacity must be finite and between zero and one.",
                code="INVALID_HIGHLIGHT",
            )


@dataclass(frozen=True, slots=True)
class _Cell:
    text: str
    columns: int
    symbol_index: int | None
    coordinate_start: int | None
    coordinate_end: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class _Row:
    cells: tuple[_Cell, ...]
    columns: int


def _coerce_input(value: DNA | DNASequence | DNARecord) -> tuple[DNASequence, str | None, int]:
    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return (
            sequence,
            None if record is None else record.id,
            0 if record is None else len(record.features),
        )
    raise TypeError("plot_sequence() accepts DNA, DNASequence, or DNARecord.")


def _gap_label(gap: Gap) -> str:
    return "[… bp]" if gap.length is None else f"[{gap.length} bp]"


def _cells(sequence: DNASequence, display_limit: int) -> tuple[_Cell, ...]:
    cells: list[_Cell] = []
    symbol_index = 0
    coordinate: int | None = 0
    for part in sequence.parts:
        if isinstance(part, Gap):
            if symbol_index >= display_limit and sequence.symbol_length > display_limit:
                break
            label = _gap_label(part)
            coordinate_end = (
                coordinate + part.length
                if coordinate is not None and part.length is not None
                else None
            )
            cells.append(
                _Cell(
                    label,
                    max(4, len(label) + 1),
                    None,
                    coordinate,
                    coordinate_end,
                    True,
                )
            )
            coordinate = coordinate_end
            continue
        for symbol in part:
            if symbol_index >= display_limit:
                return tuple(cells)
            coordinate_end = None if coordinate is None else coordinate + 1
            cells.append(_Cell(symbol, 1, symbol_index, coordinate, coordinate_end, False))
            symbol_index += 1
            coordinate = coordinate_end
    return tuple(cells)


def _rows(cells: tuple[_Cell, ...], columns_per_line: int) -> tuple[_Row, ...]:
    rows: list[_Row] = []
    current: list[_Cell] = []
    columns = 0
    for cell in cells:
        if current and columns + cell.columns > columns_per_line:
            rows.append(_Row(tuple(current), columns))
            current = []
            columns = 0
        current.append(cell)
        columns += cell.columns
    if current:
        rows.append(_Row(tuple(current), columns))
    return tuple(rows)


def _active_highlight(
    symbol_index: int, highlights: tuple[Highlight, ...]
) -> tuple[int, Highlight] | None:
    candidates = tuple(
        (index, highlight)
        for index, highlight in enumerate(highlights)
        if highlight.start <= symbol_index < highlight.end
    )
    if not candidates:
        return None
    return max(candidates, key=lambda item: (item[1].priority, item[0]))


def _coordinate_label(value: int | None, offset: int) -> str:
    return "?" if value is None else str(value + offset)


def plot_sequence(
    value: DNA | DNASequence | DNARecord,
    *,
    highlights: Iterable[Highlight] = (),
    config: SequencePlotConfig | None = None,
) -> SVGArtifact:
    """Render symbols, explicit gaps, coordinates, and optional highlights.

    Highlight coordinates always refer to ``DNASequence.symbols`` and therefore
    exclude :class:`~dnakit.core.Gap` spans.  Record features are deliberately
    not rendered implicitly; linear feature maps remain an advanced module. The
    result uses a square canvas without a visible decorative title.
    """

    sequence, sequence_id, feature_count = _coerce_input(value)
    resolved = SequencePlotConfig() if config is None else config
    if not isinstance(resolved, SequencePlotConfig):
        raise TypeError("config must be SequencePlotConfig or None.")
    highlight_tuple = tuple(highlights)
    if any(not isinstance(item, Highlight) for item in highlight_tuple):
        raise TypeError("highlights must contain only Highlight objects.")
    if any(item.end > sequence.symbol_length for item in highlight_tuple):
        raise ConfigurationError(
            "A highlight exceeds sequence symbol_length.",
            code="HIGHLIGHT_OUT_OF_RANGE",
            context={"symbol_length": sequence.symbol_length},
        )

    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    if len(gaps) > resolved.max_gaps:
        raise ConfigurationError(
            "Sequence plot exceeds max_gaps.",
            code="VISUALIZATION_SIZE_LIMIT",
            context={"gap_count": len(gaps), "max_gaps": resolved.max_gaps},
        )
    truncated = sequence.symbol_length > resolved.max_symbols
    if truncated and resolved.limit_policy is LimitPolicy.ERROR:
        raise ConfigurationError(
            "Sequence plot exceeds max_symbols.",
            code="VISUALIZATION_SIZE_LIMIT",
            context={
                "symbol_length": sequence.symbol_length,
                "max_symbols": resolved.max_symbols,
            },
            hint="Raise max_symbols or explicitly choose limit_policy='truncate'.",
        )
    display_limit = min(sequence.symbol_length, resolved.max_symbols)
    cells = _cells(sequence, display_limit)
    rows = _rows(cells, resolved.bases_per_line)

    issues: list[Issue] = []
    if truncated:
        issues.append(
            Issue(
                "VIZ_SEQUENCE_TRUNCATED",
                IssueSeverity.WARNING,
                "The sequence text plot was explicitly truncated.",
                details={
                    "symbol_length": sequence.symbol_length,
                    "displayed_symbols": display_limit,
                },
            )
        )
    clipped_highlights = sum(item.end > display_limit for item in highlight_tuple)
    if clipped_highlights:
        issues.append(
            Issue(
                "VIZ_HIGHLIGHT_CLIPPED",
                IssueSeverity.WARNING,
                "At least one highlight extends beyond the displayed sequence prefix.",
                details={"count": clipped_highlights},
            )
        )
    if sequence.topology is Topology.CIRCULAR:
        issues.append(
            Issue(
                "VIZ_CIRCULAR_LINEARIZED",
                IssueSeverity.INFO,
                "The circular sequence is shown as a linear text view.",
            )
        )

    title = resolved.title or (
        f"DNA sequence: {sequence_id}" if sequence_id is not None else "DNA sequence"
    )
    coordinate_gutter = 72 if resolved.show_coordinates else 0
    max_columns = max((row.columns for row in rows), default=20)
    content_width = (
        max_columns * resolved.cell_width + max(0, max_columns - 1) * resolved.column_spacing
    )
    natural_width = max(
        420,
        resolved.margin * 2
        + coordinate_gutter * (2 if resolved.show_coordinates else 0)
        + content_width,
    )
    row_height = resolved.line_height + (resolved.font_size + 8 if resolved.show_complement else 0)
    row_count = max(1, len(rows))
    content_height = row_count * row_height + max(0, len(rows) - 1) * resolved.line_spacing
    natural_height = max(
        120,
        resolved.margin * 2 + content_height,
    )
    canvas_size = max(natural_width, natural_height)
    layout_left = (canvas_size - natural_width) / 2
    content_top = (canvas_size - content_height) / 2
    metadata = {
        "coordinate_space": "symbol-highlights; sequence-span-labels",
        "column_spacing": resolved.column_spacing,
        "display_start_coordinate": resolved.start_coordinate,
        "displayed_symbols": display_limit,
        "feature_count_not_rendered": feature_count,
        "gap_count": len(gaps),
        "highlight_count": len(highlight_tuple),
        "line_spacing": resolved.line_spacing,
        "sequence_id": sequence_id,
        "symbol_map": dict(resolved.symbol_map),
        "symbol_length": sequence.symbol_length,
        "topology": sequence.topology.value,
        "truncated": truncated,
        "unknown_gap_count": sum(gap.length is None for gap in gaps),
    }
    builder = SVGBuilder(
        canvas_size,
        canvas_size,
        title=title,
        description="Gap-aware DNA symbol text plot in input order.",
        kind="sequence",
        theme=resolved.theme,
        metadata=metadata,
    )

    if not rows:
        builder.text(
            canvas_size / 2,
            canvas_size / 2,
            "Empty sequence",
            fill=resolved.theme.muted,
            font_family=resolved.theme.font_family,
            font_size=resolved.font_size,
            text_anchor="middle",
            class_="empty-state",
        )
        return builder.artifact(issues=issues)

    base_colors = {
        "A": "#16a34a",
        "C": "#2563eb",
        "G": "#d97706",
        "T": "#dc2626",
    }
    for row_index, row in enumerate(rows):
        y = content_top + row_index * (row_height + resolved.line_spacing) + resolved.font_size
        x = layout_left + resolved.margin + coordinate_gutter
        if resolved.show_coordinates:
            builder.text(
                layout_left + resolved.margin,
                y,
                _coordinate_label(row.cells[0].coordinate_start, resolved.start_coordinate),
                fill=resolved.theme.muted,
                font_family=resolved.theme.font_family,
                font_size=max(9, resolved.font_size - 2),
                class_="coordinate coordinate-start",
            )
        for cell_index, cell in enumerate(row.cells):
            cell_width = (
                cell.columns * resolved.cell_width
                + max(0, cell.columns - 1) * resolved.column_spacing
            )
            if cell.is_gap:
                builder.rect(
                    x,
                    y - resolved.font_size,
                    cell_width - 2,
                    resolved.font_size + 6,
                    fill=resolved.theme.missing,
                    stroke=resolved.theme.grid,
                    stroke_dasharray="3 2",
                    rx=3,
                    class_="gap",
                    data_coordinate_start=_coordinate_label(
                        cell.coordinate_start, resolved.start_coordinate
                    ),
                    data_coordinate_end=_coordinate_label(
                        cell.coordinate_end, resolved.start_coordinate
                    ),
                )
                builder.text(
                    x + cell_width / 2,
                    y,
                    cell.text,
                    fill=resolved.theme.muted,
                    font_family=resolved.theme.font_family,
                    font_size=max(9, resolved.font_size - 2),
                    text_anchor="middle",
                    class_="gap-label",
                )
            else:
                assert cell.symbol_index is not None
                active = _active_highlight(cell.symbol_index, highlight_tuple)
                foreground = base_colors.get(cell.text, resolved.theme.muted)
                display_symbol = resolved.symbol_map.get(cell.text, cell.text)
                if active is not None:
                    highlight_index, highlight = active
                    builder.rect(
                        x,
                        y - resolved.font_size,
                        cell_width,
                        resolved.font_size + 6,
                        fill=highlight.color,
                        fill_opacity=highlight.opacity,
                        class_="highlight",
                        data_highlight_index=highlight_index,
                        data_priority=highlight.priority,
                    )
                    if highlight.foreground is not None:
                        foreground = highlight.foreground
                    if highlight.label is not None and cell.symbol_index == highlight.start:
                        builder.text(
                            x,
                            y - resolved.font_size - 3,
                            highlight.label,
                            fill=resolved.theme.foreground,
                            font_family=resolved.theme.font_family,
                            font_size=max(8, resolved.font_size - 4),
                            class_="highlight-label",
                        )
                builder.text(
                    x + cell_width / 2,
                    y,
                    display_symbol,
                    fill=foreground,
                    font_family=resolved.theme.font_family,
                    font_size=resolved.font_size,
                    font_weight="bold",
                    text_anchor="middle",
                    class_="base",
                    data_source_symbol=cell.text,
                    data_symbol_index=cell.symbol_index,
                    data_coordinate=_coordinate_label(
                        cell.coordinate_start, resolved.start_coordinate
                    ),
                )
                if resolved.show_complement:
                    complement_symbol = cell.text.translate(_COMPLEMENT)
                    builder.text(
                        x + cell_width / 2,
                        y + resolved.font_size + 8,
                        resolved.symbol_map.get(complement_symbol, complement_symbol),
                        fill=foreground,
                        font_family=resolved.theme.font_family,
                        font_size=resolved.font_size,
                        text_anchor="middle",
                        class_="complement",
                        data_source_symbol=complement_symbol,
                        data_symbol_index=cell.symbol_index,
                    )
            x += cell_width
            if cell_index + 1 < len(row.cells):
                x += resolved.column_spacing
        if resolved.show_coordinates:
            builder.text(
                x + 8,
                y,
                _coordinate_label(row.cells[-1].coordinate_end, resolved.start_coordinate),
                fill=resolved.theme.muted,
                font_family=resolved.theme.font_family,
                font_size=max(9, resolved.font_size - 2),
                class_="coordinate coordinate-end",
            )
    return builder.artifact(issues=issues)


__all__ = ["Highlight", "plot_sequence"]
