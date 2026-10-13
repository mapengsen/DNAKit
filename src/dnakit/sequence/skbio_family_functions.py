"""scikit-bio function added by family group H3 (sequence area).

``nucleotide_distance(model="jc69")`` computes the Jukes-Cantor distance from two sequences.
:func:`correct_jc69_distance` applies the same correction to p-distances (observed fractions of
differing sites) that were measured elsewhere. scikit-bio is imported only when the function runs.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-skibio"
_PACKAGE = "scikit-bio"
_MAX_VALUES = 1_000_000
_MAX_ALPHABET = 1_000


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def correct_jc69_distance(
    p_distances: Sequence[float],
    *,
    alphabet_size: int = 4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply the Jukes-Cantor (JC69) correction to observed p-distances.

    Calls ``skbio.sequence.distance.jc69_correct``: ``D = -a ln(1 - p / a)`` with
    ``a = (n - 1) / n`` for an alphabet of ``n`` symbols (4 for DNA). Each ``p`` must lie in
    [0, 1]. A value at or above ``a`` has no finite distance (saturation); it is reported as
    ``None`` in ``corrected`` and the result status is ``"partial"``.
    """

    if isinstance(p_distances, (str, bytes)) or not isinstance(p_distances, Sequence):
        raise ConfigurationError("p_distances must be a non-empty list of numbers.")
    if not 1 <= len(p_distances) <= _MAX_VALUES:
        raise ConfigurationError(f"p_distances must contain 1 to {_MAX_VALUES} numbers.")
    values: list[float] = []
    for item in p_distances:
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
            raise ConfigurationError("each p-distance must be a finite number.")
        if not 0.0 <= float(item) <= 1.0:
            raise ConfigurationError("each p-distance must lie in [0, 1].")
        values.append(float(item))
    if (
        isinstance(alphabet_size, bool)
        or not isinstance(alphabet_size, int)
        or not 2 <= alphabet_size <= _MAX_ALPHABET
    ):
        raise ConfigurationError(f"alphabet_size must be an integer in [2, {_MAX_ALPHABET}].")
    distance = require("skbio.sequence.distance", _EXTRA)
    numpy = require("numpy", _EXTRA)
    result = distance.jc69_correct(numpy.array(values, dtype="float64"), chars=alphabet_size)
    corrected = [None if math.isnan(float(item)) else float(item) for item in result]
    saturated = sum(item is None for item in corrected)
    return {
        "status": "partial" if saturated else "completed",
        "operation": "correct_jc69_distance",
        "alphabet_size": alphabet_size,
        "saturation_limit": (alphabet_size - 1) / alphabet_size,
        "count": len(corrected),
        "saturated_count": saturated,
        "p_distances": values,
        "corrected": corrected,
        "provenance": _provenance(
            parser="skbio.sequence.distance.jc69_correct",
            undefined="p at or above (n - 1) / n is reported as null",
        ),
    }


__all__ = ["correct_jc69_distance"]

_register_public_verbose(globals(), tuple(__all__))
