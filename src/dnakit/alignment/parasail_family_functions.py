"""parasail alignment kernels selected by name parts, and SSW with a reusable query profile.

``parasail`` ships about 2,400 alignment functions whose names are built from the same parts:
``<algorithm>[_<output>][_<method>][_profile][_<width>]``. ``parasail_dp_variant()`` builds that
name from keyword arguments, so every such member is reachable without a wrapper of its own.
``parasail_ssw_profile_align()`` builds the SSW query profile once and aligns it to many targets.
parasail is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment.backends import (
    _DEFAULT_MAX_CELLS,
    _MAX_CELLS,
    _PARASAIL_ALPHABET,
    _PARASAIL_EXTRA,
    _choose,
    _column_counts,
    _integer,
    _package,
)
from dnakit.alignment.parasail_functions import (
    _not_found_fields,
    _scoring,
    _sequences,
    _ssw_columns,
)
from dnakit.bioinformatics._common import require
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_ALGORITHMS = (
    "nw",
    "sw",
    "sg",
    "sg_qb",
    "sg_qe",
    "sg_qx",
    "sg_db",
    "sg_de",
    "sg_dx",
    "sg_qb_db",
    "sg_qb_de",
    "sg_qe_db",
    "sg_qe_de",
)
_OUTPUTS: dict[str, str] = {
    "score": "",
    "stats": "_stats",
    "table": "_table",
    "rowcol": "_rowcol",
    "stats_table": "_stats_table",
    "stats_rowcol": "_stats_rowcol",
    "trace": "_trace",
}
_METHODS = ("reference", "scan", "striped", "diag")
_WIDTHS = ("8", "16", "32", "64", "sat")
_TABLE_CELLS = 250_000
_SSW_SCORE_SIZES = (0, 1, 2)
_MAX_TARGETS = 1_000


def _width_label(value: object) -> str:
    if value is None:
        return ""
    label = str(value) if isinstance(value, int) and not isinstance(value, bool) else value
    if not isinstance(label, str) or label not in _WIDTHS:
        raise ConfigurationError("width must be 8, 16, 32, 64 or 'sat'.")
    return label


def _function_name(algorithm: str, output: str, method: str, profile: bool, width: str) -> str:
    """Compose the parasail function name from its parts (the single source of the mapping)."""

    name = algorithm + _OUTPUTS[output]
    if method != "reference":
        name += f"_{method}"
        if profile:
            name += "_profile"
        if width:
            name += f"_{width}"
    return name


def _ints(values: Any) -> list[int]:
    return [int(value) for value in values]


def _int_table(table: Any, rows: int) -> list[list[int]]:
    return [_ints(table[index]) for index in range(rows)]


def parasail_dp_variant(
    query: str | DNA,
    target: str | DNA,
    *,
    algorithm: str = "nw",
    output: str = "score",
    method: str = "striped",
    width: int | str | None = "sat",
    profile: bool = False,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one parasail dynamic-programming kernel chosen by name parts.

    ``algorithm`` is ``nw`` (global), ``sw`` (local), ``sg`` (semi-global) or ``sg_<ends>`` with the
    free ends ``qb``, ``qe``, ``qx``, ``db``, ``de``, ``dx`` and the pairs ``qb_db``, ``qb_de``,
    ``qe_db``, ``qe_de``. ``output`` selects what the kernel returns: ``score``, ``stats``
    (matches, similar columns and alignment length at the best cell), ``table`` (the whole score
    matrix), ``rowcol`` (last row and column), ``stats_table``, ``stats_rowcol`` or ``trace``
    (the aligned strings and CIGAR). ``method`` is ``reference`` (the plain kernel; ``width`` and
    ``profile`` do not apply), ``scan``, ``striped`` or ``diag``. ``profile=True`` builds a query
    profile first and is available for ``scan`` and ``striped``.
    ``width`` is 8, 16, 32, 64, ``sat`` (retry wider on saturation) or ``None`` for the kernels
    without a width suffix such as ``nw_scan``.
    A saturated narrow kernel raises ``BackendExecutionError``.
    """

    algorithm_name = _choose(algorithm, _ALGORITHMS, "algorithm")
    output_name = _choose(output, tuple(_OUTPUTS), "output")
    method_name = _choose(method, _METHODS, "method")
    if not isinstance(profile, bool):
        raise ConfigurationError("profile must be True or False.")
    if profile and method_name not in ("scan", "striped"):
        raise ConfigurationError("profile=True needs method='scan' or method='striped'.")
    if profile and "stats" in output_name:
        raise ConfigurationError(
            "profile=True is refused for stats outputs: the parasail stats profile kernels "
            "crash the interpreter. Use profile=False."
        )
    width_text = _width_label(width)
    if profile and not width_text:
        raise ConfigurationError("profile=True needs a width.")
    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    query_text, target_text = _sequences(query, target)
    cells = len(query_text) * len(target_text)
    limit = _integer(max_cells, "max_cells", minimum=1, maximum=_MAX_CELLS)
    if cells > limit:
        raise ConfigurationError(f"alignment needs more than max_cells={max_cells} DP cells.")
    if "table" in output_name and cells > _TABLE_CELLS:
        raise ConfigurationError(
            f"table outputs hold every DP cell; {cells} cells exceed the {_TABLE_CELLS} limit."
        )
    name = _function_name(algorithm_name, output_name, method_name, profile, width_text)
    parasail = require("parasail", _PARASAIL_EXTRA)
    function = getattr(parasail, name, None)
    if not callable(function):
        raise ConfigurationError(f"The installed parasail has no {name}().")
    matrix = parasail.matrix_create(_PARASAIL_ALPHABET, scoring["match"], scoring["mismatch"])
    try:
        if profile:
            builder = getattr(parasail, f"profile_create_{width_text}", None)
            if not callable(builder):
                raise ConfigurationError(
                    f"The installed parasail has no profile_create_{width_text}()."
                )
            first = builder(query_text, matrix)
            result = function(first, target_text, scoring["gap_open"], scoring["gap_extend"])
        else:
            result = function(
                query_text, target_text, scoring["gap_open"], scoring["gap_extend"], matrix
            )
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"parasail could not align the sequences: {exc}") from exc
    if result is None:
        raise BackendExecutionError(f"parasail {name}() returned no result.")
    if bool(result.saturated):
        raise BackendExecutionError(
            f"parasail {name}() saturated for these inputs; use a wider width or 'sat'."
        )
    payload: dict[str, Any] = {
        "status": "completed",
        "operation": "parasail_dp_variant",
        "parasail_function": name,
        "algorithm": algorithm_name,
        "output": output_name,
        "method": method_name,
        "width": width_text or None if method_name != "reference" else None,
        "profile": profile,
        "score": int(result.score),
        "end_query": int(result.end_query),
        "end_target": int(result.end_ref),
        "query_length": len(query_text),
        "target_length": len(target_text),
    }
    if "stats" in output_name:
        payload["matches"] = int(result.matches)
        payload["similar"] = int(result.similar)
        payload["length"] = int(result.length)
    if "rowcol" in output_name:
        payload["score_row"] = _ints(result.score_row)
        payload["score_col"] = _ints(result.score_col)
        if "stats" in output_name:
            for field in ("matches", "similar", "length"):
                payload[f"{field}_row"] = _ints(getattr(result, f"{field}_row"))
                payload[f"{field}_col"] = _ints(getattr(result, f"{field}_col"))
    if "table" in output_name:
        rows = len(query_text)
        payload["score_table"] = _int_table(result.score_table, rows)
        if "stats" in output_name:
            for field in ("matches", "similar", "length"):
                payload[f"{field}_table"] = _int_table(getattr(result, f"{field}_table"), rows)
    if output_name == "trace":
        trace = result.get_traceback()
        cigar = result.get_cigar()
        payload["query_aligned"] = str(trace.query)
        payload["target_aligned"] = str(trace.ref)
        payload["cigar"] = cigar.decode.decode("ascii")
        payload["query_offset"] = int(cigar.beg_query)
        payload["target_offset"] = int(cigar.beg_ref)
    payload["provenance"] = {
        "package": _package("parasail"),
        "algorithm": f"parasail {name}",
        "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
        "coordinate_system": "end_query and end_target are zero-based indices of the best cell",
        "alphabet": "A, C, G, T and N only",
    }
    return payload


def parasail_ssw_profile_align(
    query: str | DNA,
    targets: Sequence[str | DNA],
    *,
    score_size: int = 2,
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align one query to many targets with SSW, building the query profile once.

    ``score_size`` is the width parasail uses to build the profile: 0 (8-bit), 1 (16-bit) or 2
    (both). Other values are refused because parasail does not check them. Results keep the order
    of ``targets``; each hit has the same fields as ``parasail_ssw_align()``.
    """

    if isinstance(targets, (str, bytes, DNA)) or not isinstance(targets, Sequence):
        raise ConfigurationError("targets must be a list of DNA sequences.")
    if not 1 <= len(targets) <= _MAX_TARGETS:
        raise ConfigurationError(f"targets must hold 1 to {_MAX_TARGETS} sequences.")
    size = _integer(score_size, "score_size", minimum=0, maximum=2)
    if size not in _SSW_SCORE_SIZES:
        raise ConfigurationError("score_size must be 0, 1 or 2.")
    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    limit = _integer(max_cells, "max_cells", minimum=1, maximum=_MAX_CELLS)
    pairs = [_sequences(query, target) for target in targets]
    query_text = pairs[0][0]
    if any(len(query_text) * len(target_text) > limit for _, target_text in pairs):
        raise ConfigurationError(f"an alignment needs more than max_cells={max_cells} DP cells.")
    parasail = require("parasail", _PARASAIL_EXTRA)
    matrix = parasail.matrix_create(_PARASAIL_ALPHABET, scoring["match"], scoring["mismatch"])
    try:
        built = parasail.ssw_init(query_text, matrix, size)
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"parasail could not build the SSW profile: {exc}") from exc
    hits: list[dict[str, Any]] = []
    for index, (_, target_text) in enumerate(pairs):
        try:
            result = parasail.ssw_profile(
                built, target_text, scoring["gap_open"], scoring["gap_extend"]
            )
        except (TypeError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(f"parasail could not align target {index}: {exc}") from exc
        if result is None:
            raise BackendExecutionError(
                f"parasail ssw returned no result for target {index}; the score is near the "
                "16-bit limit (about 32767)."
            )
        score = int(result.score1)
        count = int(result.cigarLen)
        hit: dict[str, Any] = {
            "target_index": index,
            "score": score,
            "target_length": len(target_text),
        }
        if score <= 0 or count == 0:
            hits.append({**hit, **_not_found_fields(), "cigar": ""})
            continue
        values = [int(result.cigar[position]) for position in range(count)]
        read_begin, ref_begin = int(result.read_begin1), int(result.ref_begin1)
        query_aligned, target_aligned, cigar_text, query_index, target_index = _ssw_columns(
            values, query_text, target_text, read_begin, ref_begin
        )
        if query_index != int(result.read_end1) + 1 or target_index != int(result.ref_end1) + 1:
            raise BackendExecutionError(
                "parasail ssw returned a CIGAR that does not match its spans."
            )
        columns = len(query_aligned)
        counts = _column_counts(query_aligned, target_aligned, 0, columns - 1)
        hits.append(
            {
                **hit,
                "found": True,
                "query_start": read_begin,
                "query_end": query_index,
                "target_start": ref_begin,
                "target_end": target_index,
                "cigar": cigar_text,
                "matches": counts["matches"],
                "mismatches": counts["mismatches"],
                "insertions": counts["insertions"],
                "deletions": counts["deletions"],
                "alignment_columns": columns,
                "identity": counts["matches"] / columns,
            }
        )
    return {
        "status": "completed",
        "operation": "parasail_ssw_profile_align",
        "score_size": size,
        "query_length": len(query_text),
        "target_count": len(hits),
        "hits": hits,
        "provenance": {
            "package": _package("parasail"),
            "algorithm": "parasail ssw_init + ssw_profile (Striped Smith-Waterman, best hit only)",
            "score_limit": "SSW keeps 16-bit scores; a score near 32767 returns no result",
            "coordinate_system": "zero-based half-open spans",
            "alphabet": "A, C, G, T and N only",
        },
    }


__all__ = ["parasail_dp_variant", "parasail_ssw_profile_align"]

_register_public_verbose(globals(), tuple(__all__))
