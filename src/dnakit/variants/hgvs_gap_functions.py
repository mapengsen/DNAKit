"""hgvs operations not covered by hgvs_functions: sequence-context text, variant equivalents and alignment context.

The text helpers need only the hgvs parser. ``hgvs_full_house`` and ``hgvs_variant_context_with_alignment``
take an ``AlignmentMapper`` that the user builds on a UTA database and a sequence store (hgvs does not ship
either); the alignment context also needs the ``uta_align`` package from PyPI. Nothing here opens a
connection on its own.
"""

from __future__ import annotations

import importlib
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-hgvs"
_UTA_HINT = (
    "Install the uta_align package from PyPI (python -m pip install uta-align) and build the AlignmentMapper "
    "on a UTA database connection, as described in the hgvs documentation."
)
_MAX_GROUP_SIZE = 88  # hgvs lays out 100 columns: 9 for the location, 3 for the separator


def _parser() -> Any:
    return require("hgvs.parser", _EXTRA).Parser()


def _parse(variant: str, name: str = "variant") -> Any:
    if not isinstance(variant, str) or not variant.strip():
        raise ConfigurationError(
            f"{name} must be an HGVS string such as NM_000399.3:c.2689_2690insA."
        )
    try:
        return _parser().parse_hgvs_variant(variant.strip())
    except Exception as exc:  # hgvs raises several parser-specific exceptions
        raise ConfigurationError(
            f"{name} is not a valid HGVS variant: {variant!r} ({exc})"
        ) from exc


def _span(span: Sequence[int], name: str = "span") -> tuple[int, int]:
    if len(span) != 2 or any(
        isinstance(value, bool) or not isinstance(value, int) for value in span
    ):
        raise ConfigurationError(f"{name} must be a pair of integers (start, end).")
    start, end = int(span[0]), int(span[1])
    if start < 0 or end < start:
        raise ConfigurationError(f"{name} must satisfy 0 <= start <= end.")
    return start, end


def _mapper(mapper: Any) -> Any:
    if mapper is None:
        raise ConfigurationError(
            "mapper must be an hgvs AlignmentMapper built on a UTA database connection."
        )
    return mapper


def hgvs_sequence_line(
    variant: str,
    span: Sequence[int],
    content: str,
    *,
    direction: str = "",
    verbose: bool = True,
) -> str:
    """Format one sequence line of a variant-context display (hgvs.utils.context.seq_line_fmt).

    ``variant`` is an HGVS string that supplies the accession, the variant type and the trailing comment.
    ``span`` is the (start, end) coordinate pair printed at the two ends of the line; ``direction`` (">" or "<")
    marks the strand. The result is one text line.
    """
    context = require("hgvs.utils.context", _EXTRA)
    parsed = _parse(variant)
    start, end = _span(span)
    if direction not in ("", ">", "<"):
        raise ConfigurationError("direction must be '', '>' or '<'.")
    return context.seq_line_fmt(var=parsed, span=(start, end), content=content, dir=direction)


def hgvs_pointer_line(variant: str, span: Sequence[int], *, verbose: bool = True) -> str:
    """Draw the pointer line under a variant inside a displayed span (hgvs.utils.context.pointer_line).

    The variant must be a genomic, coding or non-coding HGVS variant with a position. Insertions are marked with
    ``><``; other edits with one ``*`` per base. The result is one text line.
    """
    context = require("hgvs.utils.context", _EXTRA)
    parsed = _parse(variant)
    if getattr(parsed, "posedit", None) is None or parsed.type not in ("g", "c", "n"):
        raise ConfigurationError(
            "pointer lines need a genomic, coding or non-coding variant (g., c. or n.)."
        )
    start, end = _span(span)
    try:
        return context.pointer_line(parsed, (start, end))
    except AttributeError as exc:
        raise ConfigurationError(
            f"{variant!r} has no simple position to point at ({exc})."
        ) from exc


def format_sequence_blocks(
    sequence: str,
    *,
    start: int | None = None,
    end: int | None = None,
    group_size: int = 3,
    verbose: bool = True,
) -> str:
    """Format a DNA sequence in 3-base groups with coordinates, 1-based, from [start, end) (hgvs format_sequence).

    The output has a ruler line above each row of bases and the first base position at the left of each row.
    ``group_size`` sets the number of bases per group (hgvs default 3). Returns the full text.
    """
    context = require("hgvs.utils.context", _EXTRA)
    if not isinstance(sequence, str) or not sequence:
        raise ConfigurationError("sequence must be a non-empty string.")
    if (
        isinstance(group_size, bool)
        or not isinstance(group_size, int)
        or not 1 <= group_size <= _MAX_GROUP_SIZE
    ):
        raise ConfigurationError(f"group_size must be an integer from 1 to {_MAX_GROUP_SIZE}.")
    first = 0 if start is None else start
    last = len(sequence) if end is None else end
    if not 0 <= first < last <= len(sequence):
        raise ConfigurationError("start and end must satisfy 0 <= start < end <= len(sequence).")
    return "".join(context.format_sequence(sequence, first, last, group_size=group_size))


def hgvs_full_house(
    mapper: Any,
    variant: str,
    *,
    transcript: str | None = None,
    verbose: bool = True,
) -> dict[str, str]:
    """Express one variant as g., c., n. and p. HGVS strings through an AlignmentMapper (hgvs full_house).

    A genomic variant needs ``transcript`` when more than one transcript overlaps it (the mapper's
    ``relevant_transcripts`` decides). The mapper must be built on a UTA database connection and a sequence store.
    Returns the four HGVS strings keyed by ``g``, ``c``, ``n`` and ``p``.
    """
    context = require("hgvs.utils.context", _EXTRA)
    parsed = _parse(variant)
    try:
        result = context.full_house(_mapper(mapper), parsed, tx_ac=transcript)
    except RuntimeError as exc:
        raise ConfigurationError(str(exc)) from exc
    return {key: str(value) for key, value in result.items()}


def hgvs_variant_context_with_alignment(
    mapper: Any,
    variant: str,
    *,
    margin: int = 20,
    transcript: str | None = None,
    verbose: bool = True,
) -> str:
    """Show a coding variant against the transcript alignment with genomic and cDNA context (hgvs variant_context_w_alignment).

    The genomic and transcript sequences are aligned with uta_align (global alignment), and the result is a
    multi-line text figure. ``margin`` adds that many bases on each side. Needs the ``uta_align`` package and a
    UTA-backed mapper; without the package a BackendUnavailableError with an install hint is raised.
    """
    try:
        importlib.import_module("uta_align.align.algorithms")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"uta_align could not be imported: {exc}", hint=_UTA_HINT
        ) from exc
    context = require("hgvs.utils.context", _EXTRA)
    parsed = _parse(variant)
    if margin < 0:
        raise ConfigurationError("margin must not be negative.")
    try:
        return context.variant_context_w_alignment(
            _mapper(mapper), parsed, margin=margin, tx_ac=transcript
        )
    except RuntimeError as exc:
        raise ConfigurationError(str(exc)) from exc


__all__ = [
    "format_sequence_blocks",
    "hgvs_full_house",
    "hgvs_pointer_line",
    "hgvs_sequence_line",
    "hgvs_variant_context_with_alignment",
]
_register_public_verbose(globals())
