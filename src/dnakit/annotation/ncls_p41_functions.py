"""Limited-overlap NCLS queries on genomic integer intervals (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def ncls_limited_overlaps(
    starts: list[int],
    ends: list[int],
    query_starts: list[int],
    query_ends: list[int],
    *,
    selection: str = "first",
    k: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """First, last or up to k stored overlaps per query, exactly as NCLS64 returns them.

    Inputs are half-open integer genomic intervals on a single chromosome. IDs
    are input row positions. selection is first, last or k; upstream order is kept.
    """
    h.choice(selection, ("first", "last", "k"), "selection")
    if len(starts) != len(ends) or len(query_starts) != len(query_ends):
        raise ConfigurationError("Start and end arrays must have equal lengths.")
    if any(type(x) is not int or x < 0 for x in starts + ends + query_starts + query_ends) or any(
        a >= b for a, b in zip(starts + query_starts, ends + query_ends, strict=False)
    ):
        raise ConfigurationError("Require integer intervals with 0 <= start < end.")
    np, module = require("numpy", "bio-analysis"), require("ncls", "bio-ncls")
    index = module.NCLS(
        np.array(starts, dtype="int64"),
        np.array(ends, dtype="int64"),
        np.arange(len(starts), dtype="int64"),
    )
    args = (
        np.array(query_starts, dtype="int64"),
        np.array(query_ends, dtype="int64"),
        np.arange(len(query_starts), dtype="int64"),
    )
    if selection == "k":
        result = index.k_overlaps_both(*args, h.integer(k, "k", minimum=1))
    else:
        result = (index.first_overlap_both if selection == "first" else index.last_overlap_both)(
            *args
        )
    return {
        "status": "completed",
        "backend": "ncls",
        "selection": selection,
        "query_indices": result[0].tolist(),
        "stored_indices": result[1].tolist(),
    }


__all__ = ["ncls_limited_overlaps"]

_register_public_verbose(globals(), tuple(__all__))
