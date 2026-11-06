"""Bounded sequential SLOW5/BLOW5 read and ID enumeration (P41)."""

from __future__ import annotations

from itertools import islice
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def slow5_read_records(
    path: str,
    *,
    include_signal: bool = True,
    max_reads: int = 1000,
    max_signal_samples: int = 1000000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate read IDs or raw read records without constructing a side-effecting index.

    include_signal=False returns IDs; True iterates raw integer signal records.
    At most max_reads records and max_signal_samples total samples are returned.
    """
    source = h.existing(path, "path")
    limit = h.integer(max_reads, "max_reads", minimum=1, maximum=100000)
    sample_limit = h.integer(max_signal_samples, "max_signal_samples", minimum=1, maximum=10000000)
    if not isinstance(include_signal, bool):
        raise ConfigurationError("include_signal must be boolean.")
    handle = require("pyslow5", "bio-pyslow5").Open(str(source), "r")
    records = []
    samples = 0
    truncated = False
    try:
        for read in islice(handle.seq_reads(pA=False), limit + 1):
            if len(records) == limit:
                truncated = True
                break
            if include_signal:
                samples += len(read["signal"])
                if samples > sample_limit:
                    raise ConfigurationError("Raw signals exceed max_signal_samples.")
                records.append(
                    {str(k): (v.tolist() if hasattr(v, "tolist") else v) for k, v in read.items()}
                )
            else:
                records.append(str(read["read_id"]))
    finally:
        handle.close()
    return {
        "status": "completed",
        "backend": "pyslow5",
        "records": records,
        "include_signal": include_signal,
        "truncated": truncated,
    }


__all__ = ["slow5_read_records"]
_register_public_verbose(globals(), tuple(__all__))
