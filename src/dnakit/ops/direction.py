"""Direction-changing operations for immutable DNA sequences."""

from __future__ import annotations

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core import DNASequence
from dnakit.ops._common import require_sequence


def reverse(sequence: DNASequence, *, verbose: bool = True) -> DNASequence:
    """Return a new sequence with part order and nucleotide order reversed."""

    return require_sequence(sequence).reverse()


def complement(sequence: DNASequence, *, verbose: bool = True) -> DNASequence:
    """Return a new sequence using the complete DNA IUPAC complement table."""

    return require_sequence(sequence).complement()


def reverse_complement(sequence: DNASequence, *, verbose: bool = True) -> DNASequence:
    """Return a new reverse-complemented sequence, preserving explicit gaps."""

    return require_sequence(sequence).reverse_complement()


__all__ = ["complement", "reverse", "reverse_complement"]

_register_public_verbose(globals())
