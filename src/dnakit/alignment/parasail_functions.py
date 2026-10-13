"""Parasail DNA alignments beyond the global, local and semi-global adapters in backends.py.

Each wrapper stands for one parasail family and selects its member with keyword arguments:

- ``parasail_sg_<ends>_align()`` covers ten semi-global families. ``<ends>`` names the sequence
  ends that may stay unaligned: ``q`` query, ``d`` database, ``b`` begin, ``e`` end. The member
  used is ``sg_<ends>_trace_<method>_<width>``; ``method`` and ``width`` choose the SIMD variant.
- ``parasail_nw_banded_score()`` returns a global score restricted to a diagonal band.
- ``parasail_ssw_align()`` returns the best local alignment from Striped Smith-Waterman (SSW).

Plain global, local and semi-global alignments are already exposed by ``align_dna_parasail()``
and ``score_dna_parasail()``. parasail is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment.backends import (
    _DEFAULT_MAX_CELLS,
    _DNA_LETTERS,
    _MAX_CELLS,
    _PARASAIL_ALPHABET,
    _PARASAIL_EXTRA,
    _SEQUENCE_LIMIT,
    _aligned_region,
    _choose,
    _column_counts,
    _integer,
    _letters,
    _package,
)
from dnakit.bioinformatics._common import require
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_METHODS = ("diag", "scan", "striped")
_WIDTHS = (8, 16, 32, 64)
_BAND_LIMIT = 100_000
_SSW_OPERATIONS = "MIDNSHP=X"
_SEMIGLOBAL_ENDS: dict[str, tuple[str, ...]] = {
    "qb": ("query_begin",),
    "qe": ("query_end",),
    "qx": ("query_begin", "query_end"),
    "db": ("database_begin",),
    "de": ("database_end",),
    "dx": ("database_begin", "database_end"),
    "qb_db": ("query_begin", "database_begin"),
    "qb_de": ("query_begin", "database_end"),
    "qe_db": ("query_end", "database_begin"),
    "qe_de": ("query_end", "database_end"),
}


def _scoring(
    match: object, mismatch: object, gap_open: object, gap_extend: object
) -> dict[str, int]:
    return {
        "match": _integer(match, "match", minimum=-1_000, maximum=1_000),
        "mismatch": _integer(mismatch, "mismatch", minimum=-1_000, maximum=1_000),
        "gap_open": _integer(gap_open, "gap_open", minimum=0, maximum=1_000),
        "gap_extend": _integer(gap_extend, "gap_extend", minimum=0, maximum=1_000),
    }


def _width(value: object) -> int | str:
    if isinstance(value, str) and value == "sat":
        return "sat"
    if isinstance(value, bool) or not isinstance(value, int) or value not in _WIDTHS:
        raise ConfigurationError("width must be 8, 16, 32, 64 or 'sat'.")
    return value


def _sequences(query: object, target: object) -> tuple[str, str]:
    query_text = _letters(query, "query", allowed=_DNA_LETTERS, maximum=_SEQUENCE_LIMIT)
    target_text = _letters(target, "target", allowed=_DNA_LETTERS, maximum=_SEQUENCE_LIMIT)
    return query_text, target_text


def _cell_limit(value: object) -> int:
    return _integer(value, "max_cells", minimum=1, maximum=_MAX_CELLS)


def _run(parasail: Any, name: str, *args: Any) -> Any:
    function = getattr(parasail, name, None)
    if not callable(function):
        raise ConfigurationError(f"The installed parasail lacks {name}().")
    try:
        result = function(*args)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"parasail could not align the sequences: {exc}") from exc
    if result is None:
        raise BackendExecutionError(
            f"parasail {name}() returned no result; for ssw this means the score is near the "
            "16-bit limit (about 32767)."
        )
    return result


def _matrix(parasail: Any, scoring: Mapping[str, int]) -> Any:
    return parasail.matrix_create(_PARASAIL_ALPHABET, scoring["match"], scoring["mismatch"])


def _not_found_fields() -> dict[str, Any]:
    return {
        "found": False,
        "query_start": None,
        "query_end": None,
        "target_start": None,
        "target_end": None,
        "alignment_columns": 0,
        "matches": 0,
        "mismatches": 0,
        "insertions": 0,
        "deletions": 0,
        "identity": None,
    }


def _semiglobal_align(
    ends: str,
    query: object,
    target: object,
    *,
    match: object,
    mismatch: object,
    gap_open: object,
    gap_extend: object,
    method: object,
    width: object,
    max_cells: object,
) -> dict[str, Any]:
    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    method_name = _choose(method, _METHODS, "method")
    width_value = _width(width)
    width_label = width_value if isinstance(width_value, str) else str(width_value)
    query_text, target_text = _sequences(query, target)
    if len(query_text) * len(target_text) > _cell_limit(max_cells):
        raise ConfigurationError(f"alignment needs more than max_cells={max_cells} DP cells.")
    parasail = require("parasail", _PARASAIL_EXTRA)
    result = _run(
        parasail,
        f"sg_{ends}_trace_{method_name}_{width_label}",
        query_text,
        target_text,
        scoring["gap_open"],
        scoring["gap_extend"],
        _matrix(parasail, scoring),
    )
    if bool(result.saturated):
        raise BackendExecutionError(
            f"parasail width {width_label} saturated for these inputs; use a wider width or 'sat'."
        )
    trace = result.get_traceback()
    cigar = result.get_cigar()
    query_aligned = str(trace.query)
    target_aligned = str(trace.ref)
    offset_query, offset_target = int(cigar.beg_query), int(cigar.beg_ref)
    q_begin, q_finish, t_begin, t_finish, first, last = _aligned_region(
        query_aligned, target_aligned
    )
    columns = last - first + 1 if last >= first else 0
    found = columns > 0
    counts = (
        _column_counts(query_aligned, target_aligned, first, last)
        if found
        else {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0}
    )
    free_ends = list(_SEMIGLOBAL_ENDS[ends])
    span = (
        {
            "query_start": q_begin + offset_query,
            "query_end": q_finish + offset_query,
            "target_start": t_begin + offset_target,
            "target_end": t_finish + offset_target,
        }
        if found
        else {"query_start": None, "query_end": None, "target_start": None, "target_end": None}
    )
    return {
        "status": "completed",
        "operation": f"parasail_sg_{ends}_align",
        "free_end_gaps": free_ends,
        "method": method_name,
        "width": width_value,
        "score": int(result.score),
        "query_length": len(query_text),
        "target_length": len(target_text),
        "query_aligned": query_aligned,
        "target_aligned": target_aligned,
        "cigar": cigar.decode.decode("ascii"),
        "found": found,
        **span,
        "matches": counts["matches"],
        "mismatches": counts["mismatches"],
        "insertions": counts["insertions"],
        "deletions": counts["deletions"],
        "alignment_columns": columns,
        "identity": counts["matches"] / columns if found else None,
        "provenance": {
            "package": _package("parasail"),
            "algorithm": f"parasail sg_{ends}_trace (semi-global; free end gaps: "
            f"{', '.join(free_ends)})",
            "variant": f"method={method_name}, width={width_label}",
            "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
            "cigar": "'=' match, 'X' mismatch, 'I' query base absent from target, "
            "'D' target base absent from query; the CIGAR covers the whole traceback, "
            "including free end gaps",
            "coordinate_system": "zero-based half-open spans; spans exclude free end gaps",
            "alphabet": "A, C, G, T and N only",
        },
    }


def parasail_sg_qb_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qb: the query begin is free, so a query prefix may stay unaligned."""

    return _semiglobal_align(
        "qb",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qe_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qe: the query end is free, so a query suffix may stay unaligned."""

    return _semiglobal_align(
        "qe",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qx_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qx: both query ends are free, so the query may overhang the target."""

    return _semiglobal_align(
        "qx",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_db_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_db: the database begin is free, so a target prefix may be skipped."""

    return _semiglobal_align(
        "db",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_de_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_de: the database end is free, so a target suffix may be skipped."""

    return _semiglobal_align(
        "de",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_dx_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_dx: both target ends are free, so the query fits inside the target."""

    return _semiglobal_align(
        "dx",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qb_db_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qb_db: query begin and database begin are free (suffix-suffix)."""

    return _semiglobal_align(
        "qb_db",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qb_de_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qb_de: query begin and database end are free (query suffix overlaps
    the target prefix)."""

    return _semiglobal_align(
        "qb_de",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qe_db_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qe_db: query end and database begin are free (query prefix overlaps
    the target suffix)."""

    return _semiglobal_align(
        "qe_db",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_sg_qe_de_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    method: str = "scan",
    width: int | str = "sat",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align with parasail sg_qe_de: query end and database end are free (prefix-prefix)."""

    return _semiglobal_align(
        "qe_de",
        query,
        target,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
        method=method,
        width=width,
        max_cells=max_cells,
    )


def parasail_nw_banded_score(
    query: str | DNA,
    target: str | DNA,
    *,
    band: int = 10,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score a global DNA alignment whose path stays within a diagonal band of half-width band."""

    half_width = _integer(band, "band", minimum=1, maximum=_BAND_LIMIT)
    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    query_text, target_text = _sequences(query, target)
    if len(query_text) * (2 * half_width + 1) > _cell_limit(max_cells):
        raise ConfigurationError(f"banded alignment needs more than max_cells={max_cells} cells.")
    parasail = require("parasail", _PARASAIL_EXTRA)
    result = _run(
        parasail,
        "nw_banded",
        query_text,
        target_text,
        scoring["gap_open"],
        scoring["gap_extend"],
        half_width,
        _matrix(parasail, scoring),
    )
    if bool(result.saturated):
        raise BackendExecutionError("parasail reported a saturated score for these inputs.")
    return {
        "status": "completed",
        "operation": "parasail_nw_banded_score",
        "score": int(result.score),
        "band": half_width,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "provenance": {
            "package": _package("parasail"),
            "algorithm": "parasail nw_banded (global, affine gaps, score only)",
            "band": "cells with |i - j| <= band around the main diagonal; parasail widens the "
            "band when the lengths differ so that the final cell stays reachable",
            "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
            "alphabet": "A, C, G, T and N only",
        },
    }


def _ssw_columns(
    values: list[int], query_text: str, target_text: str, query_index: int, target_index: int
) -> tuple[str, str, str, int, int]:
    query_parts: list[str] = []
    target_parts: list[str] = []
    cigar_parts: list[str] = []
    for value in values:
        length, operation = value >> 4, value & 0xF
        if operation >= len(_SSW_OPERATIONS):
            raise BackendExecutionError(f"parasail ssw returned an unknown CIGAR code {value}.")
        symbol = _SSW_OPERATIONS[operation]
        cigar_parts.append(f"{length}{symbol}")
        if symbol in "M=X":
            query_next, target_next = query_index + length, target_index + length
            if query_next > len(query_text) or target_next > len(target_text):
                raise BackendExecutionError("parasail ssw returned a CIGAR longer than its input.")
            query_parts.append(query_text[query_index:query_next])
            target_parts.append(target_text[target_index:target_next])
            query_index, target_index = query_next, target_next
        elif symbol == "I":
            query_parts.append(query_text[query_index : query_index + length])
            target_parts.append("-" * length)
            query_index += length
        elif symbol == "D":
            query_parts.append("-" * length)
            target_parts.append(target_text[target_index : target_index + length])
            target_index += length
        else:
            raise BackendExecutionError(f"parasail ssw returned unsupported CIGAR op {symbol!r}.")
    return (
        "".join(query_parts),
        "".join(target_parts),
        "".join(cigar_parts),
        query_index,
        target_index,
    )


def parasail_ssw_align(
    query: str | DNA,
    target: str | DNA,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the best local DNA alignment with parasail's Striped Smith-Waterman (SSW) kernel."""

    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    query_text, target_text = _sequences(query, target)
    if len(query_text) * len(target_text) > _cell_limit(max_cells):
        raise ConfigurationError(f"alignment needs more than max_cells={max_cells} DP cells.")
    parasail = require("parasail", _PARASAIL_EXTRA)
    result = _run(
        parasail,
        "ssw",
        query_text,
        target_text,
        scoring["gap_open"],
        scoring["gap_extend"],
        _matrix(parasail, scoring),
    )
    score = int(result.score1)
    count = int(result.cigarLen)
    base = {
        "status": "completed",
        "operation": "parasail_ssw_align",
        "score": score,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "provenance": {
            "package": _package("parasail"),
            "algorithm": "parasail ssw (Striped Smith-Waterman, local, best hit only)",
            "score_limit": "SSW keeps 16-bit scores; a score near 32767 returns no result",
            "cigar": "SAM-style operations: '=' match, 'X' mismatch, 'I' query base absent from "
            "target, 'D' target base absent from query",
            "coordinate_system": "zero-based half-open spans",
            "alphabet": "A, C, G, T and N only",
        },
    }
    if score <= 0 or count == 0:
        return {
            **base,
            **_not_found_fields(),
            "query_aligned": "",
            "target_aligned": "",
            "cigar": "",
        }
    values = [int(result.cigar[index]) for index in range(count)]
    read_begin = int(result.read_begin1)
    ref_begin = int(result.ref_begin1)
    query_aligned, target_aligned, cigar_text, query_index, target_index = _ssw_columns(
        values, query_text, target_text, read_begin, ref_begin
    )
    if query_index != int(result.read_end1) + 1 or target_index != int(result.ref_end1) + 1:
        raise BackendExecutionError("parasail ssw returned a CIGAR that does not match its spans.")
    columns = len(query_aligned)
    counts = _column_counts(query_aligned, target_aligned, 0, columns - 1)
    return {
        **base,
        "found": True,
        "query_start": read_begin,
        "query_end": query_index,
        "target_start": ref_begin,
        "target_end": target_index,
        "query_aligned": query_aligned,
        "target_aligned": target_aligned,
        "cigar": cigar_text,
        "matches": counts["matches"],
        "mismatches": counts["mismatches"],
        "insertions": counts["insertions"],
        "deletions": counts["deletions"],
        "alignment_columns": columns,
        "identity": counts["matches"] / columns,
    }


__all__ = [
    "parasail_nw_banded_score",
    "parasail_sg_db_align",
    "parasail_sg_de_align",
    "parasail_sg_dx_align",
    "parasail_sg_qb_align",
    "parasail_sg_qb_db_align",
    "parasail_sg_qb_de_align",
    "parasail_sg_qe_align",
    "parasail_sg_qe_db_align",
    "parasail_sg_qe_de_align",
    "parasail_sg_qx_align",
    "parasail_ssw_align",
]

_register_public_verbose(globals(), tuple(__all__))
