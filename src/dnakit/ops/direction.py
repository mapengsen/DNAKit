"""Direction-changing operations for immutable DNA sequences."""

from __future__ import annotations

from dnakit.core import DNASequence
from dnakit.ops._common import require_sequence


def reverse(sequence: DNASequence) -> DNASequence:
    """Return a new sequence with part order and nucleotide order reversed."""

    return require_sequence(sequence).reverse()


def complement(sequence: DNASequence) -> DNASequence:
    """Return a new sequence using the complete DNA IUPAC complement table."""

    return require_sequence(sequence).complement()


def reverse_complement(sequence: DNASequence) -> DNASequence:
    """Return a new reverse-complemented sequence, preserving explicit gaps."""

    return require_sequence(sequence).reverse_complement()


__all__ = ["complement", "reverse", "reverse_complement"]
