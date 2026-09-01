"""Codon-aware evolutionary selection analysis."""

from __future__ import annotations

from typing import TypeAlias

from dnakit.backends.scientific import _run_scientific_function
from dnakit.core import DNA, DNARecord, DNASequence, ProviderResult
from dnakit.exceptions import ConfigurationError
from dnakit.molbio._shared import require_sequence

DNDSInput: TypeAlias = DNA | DNARecord | DNASequence


def _coding_sequence(value: DNDSInput, name: str) -> str:
    sequence = value.sequence if isinstance(value, DNARecord) else value
    return require_sequence(
        sequence,
        operation=f"calculate_dn_ds {name}",
        max_length=3_000_000,
        canonical=True,
        allow_circular=False,
    )


def calculate_dn_ds(sequence_a: DNDSInput, sequence_b: DNDSInput) -> ProviderResult:
    """Estimate Nei-Gojobori dN/dS for two codon-aligned coding sequences."""

    first = _coding_sequence(sequence_a, "sequence_a")
    second = _coding_sequence(sequence_b, "sequence_b")
    if len(first) != len(second):
        raise ConfigurationError(
            "dN/dS inputs must have equal lengths.",
            code="DNDS_ALIGNMENT_LENGTH_MISMATCH",
            context={"sequence_a_length": len(first), "sequence_b_length": len(second)},
        )
    if len(first) < 3 or len(first) % 3:
        raise ConfigurationError(
            "dN/dS inputs must contain complete aligned codons.",
            code="INVALID_DNDS_CODON_ALIGNMENT",
            context={"sequence_length": len(first)},
        )
    return _run_scientific_function(
        "calculate_dn_ds",
        {"seq1": first, "seq2": second},
        parameters={"sequence_length": len(first), "codon_count": len(first) // 3},
    )


__all__ = ["DNDSInput", "calculate_dn_ds"]
