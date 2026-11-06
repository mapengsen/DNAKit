"""bigWig/bigBed type and schema inspection (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require


def inspect_bigwig_file(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Identify bigWig/bigBed and return a bigBed autoSQL schema, if present."""
    module = require("pyBigWig", "bio-pybigwig")
    handle = module.open(str(h.existing(path, "path")))
    try:
        bigwig, bigbed = bool(handle.isBigWig()), bool(handle.isBigBed())
        schema = handle.SQL() if bigbed else None
        if isinstance(schema, bytes):
            schema = schema.decode("utf-8")
        return {
            "status": "completed",
            "backend": "pyBigWig",
            "is_bigwig": bigwig,
            "is_bigbed": bigbed,
            "autosql": schema,
        }
    finally:
        handle.close()


__all__ = ["inspect_bigwig_file"]
_register_public_verbose(globals(), tuple(__all__))
