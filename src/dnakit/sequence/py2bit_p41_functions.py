"""py2bit regional composition and exact masking coordinates (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require


def twobit_region_composition(
    path: str,
    chrom: str,
    *,
    start: int = 0,
    end: int | None = None,
    fraction: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return py2bit A/C/G/T counts or fractions over a 0-based half-open region."""
    module = require("py2bit", "bio-py2bit")
    handle = module.open(str(h.existing(path, "path")))
    try:
        last = handle.chroms(chrom) if end is None else end
        values = handle.bases(chrom, start, last, fraction=fraction)
    finally:
        handle.close()
    return {"status": "completed", "backend": "py2bit", "bases": values, "fraction": fraction}


def twobit_masked_blocks(
    path: str, chrom: str, *, kind: str = "hard", verbose: bool = True
) -> dict[str, Any]:
    """List exact hard-masked N blocks or soft-masked lowercase blocks in 0-based half-open coordinates."""
    h.choice(kind, ("hard", "soft"), "kind")
    module = require("py2bit", "bio-py2bit")
    handle = module.open(str(h.existing(path, "path")), storeMasked=True)
    try:
        blocks = (
            handle.hardMaskedBlocks(chrom) if kind == "hard" else handle.softMaskedBlocks(chrom)
        )
    finally:
        handle.close()
    return {
        "status": "completed",
        "backend": "py2bit",
        "kind": kind,
        "blocks": [list(b) for b in blocks],
    }


__all__ = ["twobit_masked_blocks", "twobit_region_composition"]

_register_public_verbose(globals(), tuple(__all__))
