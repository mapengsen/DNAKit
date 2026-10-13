"""DnaChisel difference functions that the family check found to have no existing wrapper (H2).

Each function imports DnaChisel (and Biopython) only when it runs. A missing package raises
``BackendUnavailableError`` with the ``dnakit[bio-dnachisel]`` install hint. Positions are
0-based; segments are half-open ``[start, end)``. Every function returns a JSON-safe
dictionary with ``status``, ``operation``, the result fields and ``provenance``.
"""

from __future__ import annotations

import importlib.metadata
import re
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnachisel"
_MAX_DNA = 200_000
_MAX_FEATURES = 10_000
_FEATURE_TYPE = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,49}")


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnachisel")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnachisel", "version": installed, **notes}


def _label(value: object) -> str:
    if isinstance(value, (list, tuple)):
        return str(value[0]) if value else ""
    return str(value)


def _pair(sequence_a: object, sequence_b: object) -> tuple[str, str]:
    first = _dna(sequence_a, "sequence")
    second = _dna(sequence_b, "reference")
    if len(first) != len(second):
        raise ConfigurationError("sequence and reference must have the same length.")
    return first, second


def find_sequence_difference_positions(
    sequence: str,
    reference: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a per-position difference vector for two equal-length DNA sequences.

    ``difference_array`` has one entry for each position: ``1`` where ``sequence`` differs
    from ``reference`` and ``0`` where they agree. ``positions`` lists the differing
    0-based positions. Unlike ``find_sequence_difference_segments``, the result keeps
    one value per position instead of merged segments. No alignment is performed.
    """

    module = require("dnachisel.biotools.sequences_differences", _EXTRA)
    first, second = _pair(sequence, reference)
    flags = [int(item) for item in module.sequences_differences_array(first, second)]
    positions = [index for index, flag in enumerate(flags) if flag]
    return {
        "status": "completed",
        "operation": "find_sequence_difference_positions",
        "sequence_length": len(first),
        "difference_count": len(positions),
        "positions": positions,
        "difference_array": flags,
        "provenance": _provenance(method="dnachisel.biotools.sequences_differences_array"),
    }


def annotate_sequence_differences(
    sequence: str,
    reference: str,
    *,
    feature_type: str = "misc_feature",
    prefix: str = "#",
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate where a DNA sequence differs from a reference as labelled features.

    Each run of adjacent differing positions becomes one feature with ``start`` inclusive
    and ``end`` exclusive. Its label is ``prefix`` followed by the bases of the
    ``reference`` at that location (for example ``#A`` or ``#TT``), which is what
    ``find_sequence_difference_segments`` does not report. The two sequences must have
    the same length. Identical sequences give zero features.
    """

    biotools = require("dnachisel.biotools", _EXTRA)
    seq_record = require("Bio.SeqRecord", _EXTRA).SeqRecord
    seq_class = require("Bio.Seq", _EXTRA).Seq
    first, second = _pair(sequence, reference)
    if not isinstance(feature_type, str) or not _FEATURE_TYPE.fullmatch(feature_type):
        raise ConfigurationError("feature_type must be a name such as misc_feature.")
    if not isinstance(prefix, str) or len(prefix) > 20 or not prefix.isprintable():
        raise ConfigurationError("prefix must be printable text of at most 20 characters.")
    rows: list[dict[str, Any]] = []
    if first != second:
        try:
            annotated = biotools.annotate_differences(
                seq_record(seq_class(first), id="sequence"),
                seq_record(seq_class(second), id="reference"),
                feature_type=feature_type,
                prefix=prefix,
            )
        except (IndexError, RuntimeError, TypeError, ValueError) as exc:
            raise ConfigurationError(
                f"DnaChisel could not annotate the differences: {exc}"
            ) from exc
        rows = sorted(
            (
                {
                    "start": int(feature.location.start),
                    "end": int(feature.location.end),
                    "type": str(feature.type),
                    "label": _label(feature.qualifiers["label"]),
                }
                for feature in annotated.features
            ),
            key=lambda row: (row["start"], row["end"]),
        )
    if len(rows) > _MAX_FEATURES:
        raise ConfigurationError(f"The sequences differ in more than {_MAX_FEATURES} segments.")
    return {
        "status": "completed",
        "operation": "annotate_sequence_differences",
        "sequence_length": len(first),
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(method="dnachisel.biotools.annotate_differences"),
    }


__all__ = [
    "annotate_sequence_differences",
    "find_sequence_difference_positions",
]

_register_public_verbose(globals(), tuple(__all__))
