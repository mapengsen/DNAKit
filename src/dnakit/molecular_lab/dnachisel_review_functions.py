"""Review wrappers for DnaChisel interval and codon-table helpers (G18n).

Each function imports DnaChisel (and Biopython) only when it runs. A missing package
raises ``BackendUnavailableError`` with the ``dnakit[bio-dnachisel]`` install hint.
Windows and segments are ``(start, end)`` pairs with 0-based ``start`` inclusive and
``end`` exclusive. Every function returns a JSON-safe dictionary with ``status``,
``operation``, the result fields and ``provenance``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnachisel"
_MAX_COORDINATE = 10**12
_MAX_SPAN = 10**12
_MAX_SEGMENTS = 100_000
_MAX_GAP = 10**12


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _optional_integer(value: object, name: str) -> int | None:
    if value is None:
        return None
    return _integer(value, name, minimum=1, maximum=_MAX_GAP)


def _window(value: object, name: str) -> tuple[int, int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError(f"{name} must be a (start, end) pair of integers.")
    start = _integer(value[0], f"{name}[0]", minimum=0, maximum=_MAX_COORDINATE)
    end = _integer(value[1], f"{name}[1]", minimum=0, maximum=_MAX_COORDINATE)
    if start > end:
        raise ConfigurationError(f"{name} must satisfy start <= end.")
    return start, end


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnachisel")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnachisel", "version": installed, **notes}


def dna_window_overlap(
    window1: Sequence[int],
    window2: Sequence[int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the overlap of two genomic windows given as ``(start, end)`` pairs.

    ``overlap`` is ``[start, end]`` for the shared span, or ``None`` when the windows
    do not meet. Two windows that only touch at one coordinate give a zero-length
    overlap, as DnaChisel reports it; ``overlap_length`` is then ``0`` and
    ``overlaps`` is ``False``.
    """

    first = _window(window1, "window1")
    second = _window(window2, "window2")
    result = require("dnachisel.biotools", _EXTRA).windows_overlap(first, second)
    if result is None:
        overlap = None
        length = 0
    else:
        overlap = [int(result[0]), int(result[1])]
        length = overlap[1] - overlap[0]
    return {
        "status": "completed",
        "operation": "dna_window_overlap",
        "window1": list(first),
        "window2": list(second),
        "overlap": overlap,
        "overlap_length": length,
        "overlaps": length > 0,
        "provenance": _provenance(method="dnachisel.biotools.windows_overlap"),
    }


def subdivide_dna_window(
    window: Sequence[int],
    *,
    max_span: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a ``(start, end)`` window into consecutive sub-windows of at most ``max_span`` bases.

    The sub-windows are contiguous and cover the whole window, for example
    ``(0, 10)`` with ``max_span=4`` gives ``[0, 4]``, ``[4, 8]`` and ``[8, 10]``.
    An empty window gives no sub-windows.
    """

    span_value = _integer(max_span, "max_span", minimum=1, maximum=_MAX_SPAN)
    start, end = _window(window, "window")
    module = require("dnachisel.biotools", _EXTRA)
    parts = [[int(a), int(b)] for a, b in module.subdivide_window((start, end), span_value)]
    return {
        "status": "completed",
        "operation": "subdivide_dna_window",
        "window": [start, end],
        "max_span": span_value,
        "window_count": len(parts),
        "windows": parts,
        "provenance": _provenance(method="dnachisel.biotools.subdivide_window"),
    }


def group_nearby_dna_segments(
    segments: Sequence[Sequence[int]],
    *,
    max_start_gap: int | None = None,
    max_start_spread: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group ``(start, end)`` segments whose start positions lie close together.

    Segments are sorted by start and then end. A segment joins the current group when
    its start is less than ``max_start_gap`` after the previous segment's start and
    less than ``max_start_spread`` after the group's first start. Leave a limit as
    ``None`` to skip that test. Each group is a list of ``[start, end]`` pairs.
    """

    if isinstance(segments, (str, bytes)) or not isinstance(segments, Sequence):
        raise ConfigurationError("segments must be a list of (start, end) pairs.")
    if len(segments) > _MAX_SEGMENTS:
        raise ConfigurationError(f"segments may contain at most {_MAX_SEGMENTS} items.")
    pairs = [_window(item, f"segments[{index}]") for index, item in enumerate(segments)]
    gap = _optional_integer(max_start_gap, "max_start_gap")
    spread = _optional_integer(max_start_spread, "max_start_spread")
    module = require("dnachisel.biotools", _EXTRA)
    groups = module.group_nearby_segments(pairs, max_start_gap=gap, max_start_spread=spread)
    rows = [[[int(a), int(b)] for a, b in group] for group in groups]
    return {
        "status": "completed",
        "operation": "group_nearby_dna_segments",
        "segment_count": len(pairs),
        "max_start_gap": gap,
        "max_start_spread": spread,
        "group_count": len(rows),
        "groups": rows,
        "provenance": _provenance(method="dnachisel.biotools.group_nearby_segments"),
    }


def codon_backtranslation_table(
    table_name: str = "Standard", *, verbose: bool = True
) -> dict[str, Any]:
    """Return the back-translation table of a Biopython codon table as codon lists.

    ``table_name`` is a Biopython unambiguous DNA codon table, such as ``"Standard"``
    or ``"Vertebrate Mitochondrial"``. Each key is an amino-acid letter, ``"*"`` for
    stops, ``"START"`` for start codons or an ambiguity letter (``"X"``, ``"B"``,
    ``"J"``, ``"Z"``); each value is the sorted list of codons for that key.
    """

    codon_table = require("Bio.Data.CodonTable", _EXTRA).unambiguous_dna_by_name
    if not isinstance(table_name, str) or table_name not in codon_table:
        options = ", ".join(sorted(codon_table))
        raise ConfigurationError(f"table_name must be one of: {options}.")
    module = require("dnachisel.biotools", _EXTRA)
    raw = module.get_backtranslation_table(table_name)
    table = {str(key): sorted(str(codon) for codon in raw[key]) for key in sorted(raw)}
    return {
        "status": "completed",
        "operation": "codon_backtranslation_table",
        "table_name": table_name,
        "entry_count": len(table),
        "table": table,
        "provenance": _provenance(method="dnachisel.biotools.get_backtranslation_table"),
    }


__all__ = [
    "codon_backtranslation_table",
    "dna_window_overlap",
    "group_nearby_dna_segments",
    "subdivide_dna_window",
]

_register_public_verbose(globals(), tuple(__all__))
