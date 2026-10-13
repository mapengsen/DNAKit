"""cutadapt functions added by review group G18o (sequencing_qc area).

The one wrapper here is the untrimmed interval that cutadapt reports for adapter matches
(``cutadapt.adapters.remainder``). The report functions of cutadapt need the statistics object of
a finished command-line run, and the other helpers only format those statistics, so they are not
wrapped. cutadapt is imported only when the function runs.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require, sequence_text
from dnakit.exceptions import ConfigurationError

_CUTADAPT_EXTRA = "bio-cutadapt"
_PACKAGE = "cutadapt"
_POSITIONS = ("3prime", "5prime", "anywhere")
_MAX_ADAPTERS = 10
_MAX_ERROR_RATE = 0.5
_MAX_OVERLAP = 1_000


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _adapter_text(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("each adapter must be a non-empty DNA string.")
    return sequence_text(value)


def _adapter_object(
    module: Any, sequence: str, position: str, max_errors: float, min_overlap: int
) -> Any:
    options = {"max_errors": max_errors, "min_overlap": min_overlap, "name": "adapter"}
    try:
        if position == "3prime":
            return module.BackAdapter(sequence, **options)
        if position == "5prime":
            return module.FrontAdapter(sequence, **options)
        return module.AnywhereAdapter(sequence, **options)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"cutadapt rejected the adapter: {exc}") from exc


def adapter_untrimmed_interval(
    sequence: str,
    adapters: Sequence[str],
    *,
    position: str = "3prime",
    max_errors: float = 0.1,
    min_overlap: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the part of a read that remains after trimming one or more adapters.

    The adapters are applied in the given order. Each one is matched to the part of the read left
    after the previous adapter, as cutadapt does for chained adapters. The untrimmed interval is
    cutadapt's ``remainder`` of those matches, reported as 0-based half-open coordinates of the
    original read. An adapter that does not match is skipped. With no match, the whole read
    remains.
    """

    text = sequence_text(sequence)
    if not isinstance(adapters, Sequence) or isinstance(adapters, str):
        raise ConfigurationError("adapters must be a sequence of DNA strings.")
    if not 1 <= len(adapters) <= _MAX_ADAPTERS:
        raise ConfigurationError(f"adapters must contain 1 to {_MAX_ADAPTERS} sequences.")
    specs = [_adapter_text(item) for item in adapters]
    if position not in _POSITIONS:
        raise ConfigurationError(f"position must be one of {list(_POSITIONS)}.")
    if isinstance(max_errors, bool) or not isinstance(max_errors, (int, float)):
        raise ConfigurationError("max_errors must be a number.")
    if not 0 <= float(max_errors) <= _MAX_ERROR_RATE:
        raise ConfigurationError(f"max_errors must be in [0, {_MAX_ERROR_RATE}].")
    if (
        isinstance(min_overlap, bool)
        or not isinstance(min_overlap, int)
        or not 1 <= min_overlap <= _MAX_OVERLAP
    ):
        raise ConfigurationError(f"min_overlap must be an integer in [1, {_MAX_OVERLAP}].")
    adapters_module = require("cutadapt.adapters", _CUTADAPT_EXTRA)
    matches: list[Any] = []
    details: list[dict[str, Any]] = []
    start, stop = 0, len(text)
    for spec in specs:
        current = text[start:stop]
        matcher = _adapter_object(adapters_module, spec, position, float(max_errors), min_overlap)
        match = matcher.match_to(current)
        if match is None:
            continue
        matches.append(match)
        details.append(
            {
                "adapter": spec,
                "start": start + int(match.rstart),
                "stop": start + int(match.rstop),
                "errors": int(match.errors),
            }
        )
        relative_start, relative_stop = match.remainder_interval()
        start, stop = start + int(relative_start), start + int(relative_stop)
    if matches:
        untrimmed_start, untrimmed_stop = adapters_module.remainder(matches)
        if (untrimmed_start, untrimmed_stop) != (start, stop):
            raise ConfigurationError("cutadapt remainder disagrees with the chained matches.")
    else:
        untrimmed_start, untrimmed_stop = 0, len(text)
    return {
        "status": "completed",
        "operation": "adapter_untrimmed_interval",
        "sequence": text,
        "position": position,
        "matched_count": len(matches),
        "start": int(untrimmed_start),
        "stop": int(untrimmed_stop),
        "length": int(untrimmed_stop - untrimmed_start),
        "untrimmed_sequence": text[untrimmed_start:untrimmed_stop],
        "adapter_matches": details,
        "provenance": _provenance(method="cutadapt.adapters.remainder over chained matches"),
    }


__all__ = ["adapter_untrimmed_interval"]

_register_public_verbose(globals(), tuple(__all__))
