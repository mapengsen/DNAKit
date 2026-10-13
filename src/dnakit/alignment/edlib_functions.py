"""Edlib NICE alignment view for DNA: aligned query, match line and aligned target.

The wrapper calls ``edlib.align(task="path")`` and then ``edlib.getNiceAlignment()``. It reports the
same edit distance and extended CIGAR as ``align_dna_edlib()``, and adds the three text rows that
show how the two sequences line up. edlib is imported only when the function runs.
"""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment.backends import (
    _EDLIB_EXTRA,
    _IUPAC_EQUALITIES,
    _IUPAC_LETTERS,
    _SEQUENCE_LIMIT,
    _choose,
    _cigar_counts,
    _edlib_error,
    _flag,
    _integer,
    _letters,
    _package,
)
from dnakit.bioinformatics._common import require
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError


def _gap_symbol(value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 1
        or not value.isprintable()
        or value.isalnum()
        or value == " "
    ):
        raise ConfigurationError(
            "gap_symbol must be one printable non-alphanumeric, non-space character."
        )
    return value


def edlib_nice_alignment(
    query: str | DNA,
    target: str | DNA,
    *,
    mode: str = "NW",
    k: int = -1,
    iupac: bool = False,
    gap_symbol: str = "-",
    max_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two DNA sequences with Edlib and return the NICE rows: query, match line, target."""

    mode_name = _choose(mode, ("HW", "NW"), "mode")
    max_distance = _integer(k, "k", minimum=-1, maximum=2_147_483_647)
    ambiguous = _flag(iupac, "iupac")
    gap = _gap_symbol(gap_symbol)
    length_limit = _integer(max_length, "max_length", minimum=1, maximum=_SEQUENCE_LIMIT)
    query_text = _letters(query, "query", allowed=_IUPAC_LETTERS, maximum=length_limit)
    target_text = _letters(target, "target", allowed=_IUPAC_LETTERS, maximum=length_limit)
    edlib = require("edlib", _EDLIB_EXTRA)
    equalities = _IUPAC_EQUALITIES if ambiguous else None
    try:
        result = edlib.align(
            query_text,
            target_text,
            mode=mode_name,
            task="path",
            k=max_distance,
            additionalEqualities=equalities,
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise _edlib_error(exc, "align the sequences") from exc
    distance = int(result["editDistance"])
    base = {
        "status": "completed",
        "operation": "edlib_nice_alignment",
        "mode": mode_name,
        "k": max_distance,
        "iupac": ambiguous,
        "gap_symbol": gap,
        "query_length": len(query_text),
        "target_length": len(target_text),
        "provenance": {
            "package": _package("edlib"),
            "format": "NICE rows from edlib.getNiceAlignment: '|' equal bases, '.' mismatch, "
            "' ' gap column; gaps in query_aligned and target_aligned use gap_symbol",
            "cigar": "extended: '=' match, 'X' mismatch, 'I' query base absent from target, "
            "'D' target base absent from query",
            "coordinate_system": "zero-based half-open spans",
            "ambiguity": "IUPAC codes equal their bases when iupac=True; literal otherwise",
        },
    }
    if distance < 0:
        return {
            **base,
            "found": False,
            "edit_distance": None,
            "query_aligned": None,
            "matched_aligned": None,
            "target_aligned": None,
            "cigar": None,
            "target_start": None,
            "target_end": None,
            "identity": None,
        }
    try:
        rows = edlib.getNiceAlignment(result, query_text, target_text, gapSymbol=gap)
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as exc:
        raise BackendExecutionError(f"Edlib returned no NICE rows: {exc}") from exc
    cigar = str(result["cigar"])
    counts = _cigar_counts(cigar)
    columns = sum(counts.values())
    locations = result.get("locations") or []
    if not locations:
        raise BackendExecutionError("Edlib returned an alignment without a target location.")
    start, end = (int(value) for value in locations[0])
    return {
        **base,
        "found": True,
        "edit_distance": distance,
        "query_aligned": str(rows["query_aligned"]),
        "matched_aligned": str(rows["matched_aligned"]),
        "target_aligned": str(rows["target_aligned"]),
        "cigar": cigar,
        "target_start": start,
        "target_end": end + 1,
        "identity": counts["="] / columns if columns else None,
    }


__all__ = ["edlib_nice_alignment"]

_register_public_verbose(globals(), tuple(__all__))
