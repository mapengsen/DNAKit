"""Edlib SHW prefix alignment tasks, preserving exact native locations (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require


def edlib_prefix_alignment(
    query: str, target: str, *, task: str = "distance", k: int = -1, verbose: bool = True
) -> dict[str, Any]:
    """SHW prefix edit distance, locations or extended CIGAR path (task=distance/locations/path).

    Native locations are zero-based inclusive; unlike the older DNAKit alignment
    backend this function explicitly supports SHW and does not normalize case.
    """
    h.choice(task, ("distance", "locations", "path"), "task")
    module = require("edlib", "bio-edlib")
    result = module.align(query, target, mode="SHW", task=task, k=h.integer(k, "k", minimum=-1))
    return {
        "status": "completed",
        "backend": "edlib",
        "mode": "SHW",
        "task": task,
        "edit_distance": result["editDistance"],
        "locations": [list(x) for x in result.get("locations", [])],
        "cigar": result.get("cigar"),
    }


__all__ = ["edlib_prefix_alignment"]

_register_public_verbose(globals(), tuple(__all__))
