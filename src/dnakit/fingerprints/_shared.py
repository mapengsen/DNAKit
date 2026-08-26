"""Shared public policies and validation for native fingerprints."""

from __future__ import annotations

from enum import Enum
from typing import TypeAlias, TypeVar

from dnakit.core.facade import DNA, resolve_single_dna
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError

SequenceInput: TypeAlias = DNA | DNASequence | DNARecord


class GapEncodingPolicy(str, Enum):
    """Treatment of explicit :class:`~dnakit.core.Gap` parts."""

    ERROR = "error"
    OMIT = "omit"
    EXPAND = "expand"


class IntegerAmbiguityPolicy(str, Enum):
    """Treatment of non-canonical IUPAC symbols in integer encoding."""

    ERROR = "error"
    IUPAC = "iupac"
    SENTINEL = "sentinel"


class OneHotAmbiguityPolicy(str, Enum):
    """Treatment of non-canonical IUPAC symbols in one-hot encoding."""

    ERROR = "error"
    FRACTIONAL = "fractional"
    ZERO = "zero"


class FingerprintAmbiguityPolicy(str, Enum):
    """Treatment of IUPAC symbols in exact k-mer fingerprints."""

    ERROR = "error"
    IGNORE = "ignore"


class KmerFingerprintMode(str, Enum):
    """Aggregation applied to exact k-mer observations."""

    COUNT = "count"
    FREQUENCY = "frequency"
    BINARY = "binary"


class FingerprintRepresentation(str, Enum):
    """Materialized fingerprint storage representation."""

    DENSE = "dense"
    SPARSE = "sparse"


EnumType = TypeVar("EnumType", bound=Enum)


def sequence_and_id(value: SequenceInput) -> tuple[DNASequence, str | None]:
    """Return a sequence and the optional identifier of its source record."""

    if isinstance(value, (DNA, DNARecord, DNASequence)):
        sequence, record = resolve_single_dna(value)
        return sequence, None if record is None else record.id
    raise ConfigurationError(
        "Fingerprint input must be a DNASequence or DNARecord.",
        context={"input_type": type(value).__name__},
    )


def coerce_enum(value: EnumType | str, enum_type: type[EnumType], name: str) -> EnumType:
    """Coerce a string enum value with a stable configuration diagnostic."""

    try:
        return value if isinstance(value, enum_type) else enum_type(value)
    except (TypeError, ValueError) as exc:
        choices = ", ".join(repr(item.value) for item in enum_type)
        raise ConfigurationError(
            f"Unknown {name}.",
            context={name: value},
            hint=f"Choose one of: {choices}.",
        ) from exc


def validate_bool(value: bool, name: str) -> None:
    """Reject truthy integers where a real boolean is required."""

    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.", context={name: value})


def validate_positive_int(value: int, name: str) -> None:
    """Validate a positive integer while rejecting booleans."""

    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            context={name: value},
        )


__all__ = [
    "FingerprintAmbiguityPolicy",
    "FingerprintRepresentation",
    "GapEncodingPolicy",
    "IntegerAmbiguityPolicy",
    "KmerFingerprintMode",
    "OneHotAmbiguityPolicy",
    "SequenceInput",
]
