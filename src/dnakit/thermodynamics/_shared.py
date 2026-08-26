"""Shared validation and provenance helpers for native thermodynamics."""

from __future__ import annotations

from dnakit.core import (
    DNA,
    Citation,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    Topology,
)
from dnakit.core.facade import resolve_single_dna
from dnakit.core.gap import Gap
from dnakit.exceptions import (
    ConfigurationError,
    InvalidAlphabetError,
    UnsupportedGapOperationError,
)


def canonical_linear_symbols(
    sequence: DNA | DNASequence,
    *,
    operation: str,
    min_length: int,
    max_length: int,
) -> str:
    if isinstance(sequence, DNA):
        sequence, _ = resolve_single_dna(sequence)
    if not isinstance(sequence, DNASequence):
        raise ConfigurationError(
            f"{operation} requires a DNASequence.",
            code="INVALID_THERMODYNAMIC_SEQUENCE",
        )
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            f"{operation} does not define terminal corrections for circular DNA.",
            code="CIRCULAR_THERMODYNAMICS_UNSUPPORTED",
        )
    if any(isinstance(part, Gap) for part in sequence.parts):
        raise UnsupportedGapOperationError(
            f"{operation} cannot silently omit or bridge explicit Gap objects.",
            code="THERMODYNAMIC_GAP_UNSUPPORTED",
        )
    symbols = sequence.symbols
    if any(symbol not in "ACGT" for symbol in symbols):
        raise InvalidAlphabetError(
            f"{operation} requires unambiguous A/C/G/T symbols.",
            code="AMBIGUOUS_THERMODYNAMICS_UNSUPPORTED",
        )
    if len(symbols) < min_length:
        raise ConfigurationError(
            f"{operation} requires at least {min_length} nucleotide(s).",
            code="THERMODYNAMIC_SEQUENCE_TOO_SHORT",
            context={"sequence_length": len(symbols), "minimum": min_length},
        )
    if len(symbols) > max_length:
        raise ConfigurationError(
            f"{operation} input exceeds its explicit applicability or resource limit.",
            code="THERMODYNAMIC_SEQUENCE_TOO_LONG",
            context={"sequence_length": len(symbols), "maximum": max_length},
        )
    return symbols


def native_provenance(*, citation: bool = True) -> Provenance:
    citations = (
        (
            Citation(
                "santalucia1998",
                title="A unified view of DNA nearest-neighbor thermodynamics",
                doi="10.1073/pnas.95.4.1460",
            ),
        )
        if citation
        else ()
    )
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.REIMPLEMENTATION,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.PUBLISHED_ALGORITHM,
            citations=citations,
        )
    )


__all__ = ["canonical_linear_symbols", "native_provenance"]
