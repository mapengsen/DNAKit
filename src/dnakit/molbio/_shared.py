"""Shared validation and audit helpers for molecular-biology simulations."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import TypeVar

from dnakit.core import (
    DNA,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    Issue,
    OriginClass,
    Provenance,
    ReferenceInfo,
    Topology,
)
from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.facade import resolve_single_dna
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError, UnsupportedGapOperationError

_IUPAC = frozenset("ACGTRYSWKMBDHVN")
_CANONICAL = frozenset("ACGT")
T = TypeVar("T")
_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)

IUPAC_BASES: Mapping[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def validate_positive_int(value: object, name: str, *, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive integer.", context={name: value})
    if maximum is not None and value > maximum:
        raise ConfigurationError(
            f"{name} exceeds its supported maximum.",
            context={name: value, "maximum": maximum},
        )
    return value


def validate_non_negative_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer.", context={name: value})
    return value


def finite_fraction(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not 0.0 <= value <= 1.0
    ):
        raise ConfigurationError(f"{name} must be finite and between 0 and 1.")
    return float(value)


def validate_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def materialize_bounded(
    values: Iterable[T],
    *,
    max_items: int,
    name: str,
    reject_text: bool = True,
) -> tuple[T, ...]:
    """Materialize a finite iterable under an explicit item limit."""

    validate_positive_int(max_items, f"max_{name}", maximum=100_000_000)
    if reject_text and isinstance(values, (str, bytes)):
        raise ConfigurationError(f"{name} must be a non-text iterable.")
    resolved: list[T] = []
    for item in values:
        if len(resolved) >= max_items:
            raise ConfigurationError(f"{name} exceeds its configured item limit.")
        resolved.append(item)
    return tuple(resolved)


def validate_iupac_text(value: object, name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a DNA string.")
    resolved = value.upper()
    if not resolved and not allow_empty:
        raise ConfigurationError(f"{name} must be non-empty.")
    invalid = sorted(set(resolved) - _IUPAC)
    if invalid:
        raise InvalidAlphabetError(
            f"{name} contains symbols outside the DNA IUPAC alphabet.",
            context={"invalid_symbols": invalid},
        )
    return resolved


def require_sequence(
    sequence: DNA | DNASequence,
    *,
    operation: str,
    max_length: int,
    canonical: bool = False,
    allow_circular: bool = True,
    allow_empty: bool = False,
) -> str:
    if isinstance(sequence, DNA):
        sequence, _ = resolve_single_dna(sequence)
    if not isinstance(sequence, DNASequence):
        raise ConfigurationError(
            f"{operation} requires a DNASequence.",
            context={"input_type": type(sequence).__name__},
        )
    validate_positive_int(max_length, "max_length", maximum=100_000_000)
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            f"{operation} cannot silently omit explicit Gap objects.",
            code="UNSUPPORTED_GAPPED_MOLBIO_OPERATION",
        )
    if not allow_circular and sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            f"{operation} requires linear topology.",
            code="CIRCULAR_MOLBIO_OPERATION_UNSUPPORTED",
        )
    symbols = sequence.symbols
    if not symbols and not allow_empty:
        raise ConfigurationError(f"{operation} requires a non-empty sequence.")
    if len(symbols) > max_length:
        raise ConfigurationError(
            f"{operation} input exceeds max_length.",
            code="MOLBIO_SEQUENCE_LIMIT_EXCEEDED",
            context={"sequence_length": len(symbols), "max_length": max_length},
        )
    if canonical and set(symbols) - _CANONICAL:
        raise InvalidAlphabetError(
            f"{operation} requires canonical A/C/G/T symbols.",
            code="AMBIGUOUS_MOLBIO_SEQUENCE_UNSUPPORTED",
        )
    return symbols


def reverse_complement_text(value: str) -> str:
    return value.translate(_COMPLEMENT)[::-1]


def iupac_compatible(left: str, right: str) -> bool:
    return bool(IUPAC_BASES[left] & IUPAC_BASES[right])


def circular_slice(symbols: str, start: int, length: int) -> str:
    if length < 0 or length > len(symbols):
        raise ConfigurationError("Circular slice length must be within one template traversal.")
    if not symbols:
        return ""
    start %= len(symbols)
    return (symbols + symbols)[start : start + length]


def native_provenance(
    *,
    reimplementation: bool = False,
    reference_name: str | None = None,
    reference_version: str | None = None,
    reference_checksum: str | None = None,
) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=(
                ImplementationLabel.REIMPLEMENTATION
                if reimplementation
                else ImplementationLabel.NATIVE
            ),
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=(
                OriginClass.PUBLISHED_ALGORITHM if reimplementation else OriginClass.DNAKIT
            ),
        ),
        reference=(
            None
            if reference_name is None
            else ReferenceInfo(
                reference_name,
                version=reference_version,
                checksum=reference_checksum,
            )
        ),
    )


def adapter_provenance(*, reference_name: str, reference_version: str | None = None) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.EXTERNAL,
            origin_class=OriginClass.INTEGRATION,
        ),
        reference=ReferenceInfo(reference_name, version=reference_version),
    )


def freeze_parameters(values: Mapping[str, object]) -> FrozenDict:
    return freeze_mapping(values)


def validate_issues(issues: Iterable[Issue]) -> tuple[Issue, ...]:
    resolved = tuple(issues)
    if any(not isinstance(issue, Issue) for issue in resolved):
        raise ConfigurationError("issues must contain only Issue objects.")
    return resolved


__all__ = [
    "IUPAC_BASES",
    "adapter_provenance",
    "circular_slice",
    "finite_fraction",
    "freeze_parameters",
    "iupac_compatible",
    "materialize_bounded",
    "native_provenance",
    "require_sequence",
    "reverse_complement_text",
    "validate_issues",
    "validate_iupac_text",
    "validate_non_negative_int",
    "validate_positive_int",
    "validate_text",
]
