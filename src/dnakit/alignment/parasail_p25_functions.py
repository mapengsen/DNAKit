"""parasail traceback rendering: the three-line text view and the SAM CIGAR of one DNA alignment.

``parasail_dp_variant(output="trace")`` already runs every algorithm x output kernel and returns the
two aligned strings and the CIGAR. This module adds the part it does not return: parasail's own
three-line text view (query line, comparison line, database line) produced by
``Result.get_traceback()`` and the SAM-style CIGAR from ``Result.get_cigar()``. parasail is imported
only when the function runs.
"""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.alignment.backends import (
    _DEFAULT_MAX_CELLS,
    _MAX_CELLS,
    _PARASAIL_ALPHABET,
    _PARASAIL_EXTRA,
    _choose,
    _integer,
    _package,
)
from dnakit.alignment.parasail_family_functions import _ALGORITHMS
from dnakit.alignment.parasail_functions import _scoring, _sequences
from dnakit.bioinformatics._common import require
from dnakit.core import DNA
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_VIEWS = ("text", "cigar", "both")

__all__ = ["parasail_render_alignment"]


def parasail_render_alignment(
    query: str | DNA,
    target: str | DNA,
    *,
    algorithm: str = "sw",
    view: str = "both",
    match: int = 5,
    mismatch: int = -4,
    gap_open: int = 10,
    gap_extend: int = 1,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two DNA sequences with parasail and render the traceback as text and/or a SAM CIGAR.

    ``algorithm`` is ``nw``, ``sw`` or any ``sg`` variant (``sg``, ``sg_qb``, ``sg_qe``, ``sg_qx``,
    ``sg_db``, ``sg_de``, ``sg_dx``, ``sg_qb_db``, ``sg_qb_de``, ``sg_qe_db``, ``sg_qe_de``).
    ``view`` selects ``text`` (parasail ``Result.get_traceback()``: query line, comparison line
    where ``|`` marks an identity, and database line), ``cigar`` (``Result.get_cigar()``: SAM
    operations ``=``, ``X``, ``I``, ``D`` plus the start offsets) or ``both``.
    Covers the work-list rows "traceback alignment rendered as SAM CIGAR string" and "as three-line
    text view". Uses the striped kernel with saturation retry, like the parasail Python examples.
    """

    algorithm_name = _choose(algorithm, _ALGORITHMS, "algorithm")
    view_name = _choose(view, _VIEWS, "view")
    scoring = _scoring(match, mismatch, gap_open, gap_extend)
    query_text, target_text = _sequences(query, target)
    limit = _integer(max_cells, "max_cells", minimum=1, maximum=_MAX_CELLS)
    if len(query_text) * len(target_text) > limit:
        raise ConfigurationError(f"alignment needs more than max_cells={max_cells} DP cells.")
    parasail = require("parasail", _PARASAIL_EXTRA)
    name = f"{algorithm_name}_trace_striped_sat"
    function = getattr(parasail, name, None)
    if not callable(function):
        raise ConfigurationError(f"The installed parasail has no {name}().")
    matrix = parasail.matrix_create(_PARASAIL_ALPHABET, scoring["match"], scoring["mismatch"])
    try:
        result = function(
            query_text, target_text, scoring["gap_open"], scoring["gap_extend"], matrix
        )
    except (TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"parasail could not align the sequences: {exc}") from exc
    if result is None:
        raise BackendExecutionError(f"parasail {name}() returned no result.")
    payload: dict[str, Any] = {
        "status": "completed",
        "operation": "parasail_render_alignment",
        "parasail_function": name,
        "algorithm": algorithm_name,
        "view": view_name,
        "score": int(result.score),
        "end_query": int(result.end_query),
        "end_target": int(result.end_ref),
        "query_length": len(query_text),
        "target_length": len(target_text),
    }
    if view_name in ("text", "both"):
        trace = result.get_traceback()
        payload["query_line"] = str(trace.query)
        payload["comparison_line"] = str(trace.comp)
        payload["target_line"] = str(trace.ref)
        payload["text"] = "\n".join((str(trace.query), str(trace.comp), str(trace.ref)))
    if view_name in ("cigar", "both"):
        cigar = result.get_cigar()
        payload["cigar"] = cigar.decode.decode("ascii")
        payload["query_offset"] = int(cigar.beg_query)
        payload["target_offset"] = int(cigar.beg_ref)
    payload["provenance"] = {
        "package": _package("parasail"),
        "algorithm": f"parasail {name}",
        "gap_penalty": "a gap of length k costs gap_open + (k - 1) * gap_extend",
        "alphabet": "A, C, G, T and N only",
    }
    return payload


_register_public_verbose(globals(), tuple(__all__))
