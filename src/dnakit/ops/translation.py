"""Deterministic transcription and genetic-code translation."""

from __future__ import annotations

from typing import Literal, TypeAlias

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError, UnsupportedGapOperationError
from dnakit.ops._common import require_sequence

StrandPolicy: TypeAlias = Literal["forward", "reverse"]
StopPolicy: TypeAlias = Literal["include", "truncate", "error"]
UnknownCodonPolicy: TypeAlias = Literal["x", "error"]
IncompleteCodonPolicy: TypeAlias = Literal["ignore", "error"]

_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_IUPAC_RNA = frozenset("ACGURYSWKMBDHVN")
_DNA_COMPLEMENT = str.maketrans(
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

# NCBI genetic code table 1, represented with DNA codons.
_TABLE_1 = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def _coerce_strand(strand: StrandPolicy) -> StrandPolicy:
    if strand not in ("forward", "reverse"):
        raise ConfigurationError(
            "strand must be 'forward' or 'reverse'.",
            code="INVALID_STRAND_POLICY",
            context={"strand": strand},
        )
    return strand


def _ungapped_symbols(sequence: DNASequence, *, operation: str) -> str:
    resolved = require_sequence(sequence)
    if resolved.is_gapped:
        raise UnsupportedGapOperationError(
            f"{operation} does not silently omit Gap objects.",
            code="UNSUPPORTED_GAPPED_OPERATION",
            context={"operation": operation},
            hint="Translate or transcribe each nucleotide fragment separately.",
        )
    return resolved.symbols


def transcribe(sequence: DNASequence, *, strand: StrandPolicy = "forward") -> str:
    """Transcribe DNA to RNA on the selected strand.

    The return value is uppercase RNA text.  Explicit gaps are rejected rather
    than silently removed.  ``strand='reverse'`` transcribes the reverse
    complement of the supplied DNA representation.
    """

    symbols = _ungapped_symbols(sequence, operation="transcribe")
    if _coerce_strand(strand) == "reverse":
        symbols = symbols.translate(_DNA_COMPLEMENT)[::-1]
    return symbols.replace("T", "U")


def _coerce_nucleic_text(sequence: DNASequence | str) -> str:
    if isinstance(sequence, DNASequence):
        return _ungapped_symbols(sequence, operation="translate")
    if not isinstance(sequence, str):
        raise InvalidAlphabetError(
            "Translation input must be DNASequence or uppercase DNA/RNA text.",
            context={"type": type(sequence).__name__},
        )
    symbols = set(sequence)
    if "T" in symbols and "U" in symbols:
        raise InvalidAlphabetError(
            "Translation input cannot mix thymine and uracil.",
            code="MIXED_NUCLEIC_ACID_ALPHABET",
        )
    allowed = _IUPAC_RNA if "U" in symbols else _IUPAC_DNA
    invalid = sorted(symbols - allowed)
    if invalid:
        raise InvalidAlphabetError(
            "Translation input contains invalid or non-normalized symbols.",
            context={"invalid_symbols": invalid},
            hint="Uppercase and clean raw input before translation.",
        )
    return sequence.replace("U", "T")


def translate(
    sequence: DNASequence | str,
    *,
    frame: int = 0,
    table: int = 1,
    strand: StrandPolicy = "forward",
    stop_policy: StopPolicy = "include",
    unknown_policy: UnknownCodonPolicy = "x",
    incomplete_policy: IncompleteCodonPolicy = "ignore",
) -> str:
    """Translate DNA or RNA using NCBI genetic code table 1.

    Coordinates use a zero-based frame (0, 1, or 2).  Ambiguous codons become
    ``X`` or raise according to ``unknown_policy``.  Stop codons are included as
    ``*``, truncate translation, or raise according to ``stop_policy``.  A
    trailing partial codon is ignored or rejected according to
    ``incomplete_policy``.
    """

    if isinstance(frame, bool) or not isinstance(frame, int) or frame not in (0, 1, 2):
        raise ConfigurationError(
            "frame must be one of 0, 1, or 2.",
            code="INVALID_READING_FRAME",
            context={"frame": frame},
        )
    if isinstance(table, bool) or not isinstance(table, int) or table != 1:
        raise ConfigurationError(
            "Only NCBI genetic code table 1 is available in the MVP.",
            code="UNSUPPORTED_GENETIC_CODE",
            context={"table": table},
        )
    _coerce_strand(strand)
    if stop_policy not in ("include", "truncate", "error"):
        raise ConfigurationError(
            "stop_policy must be 'include', 'truncate', or 'error'.",
            code="INVALID_STOP_POLICY",
            context={"stop_policy": stop_policy},
        )
    if unknown_policy not in ("x", "error"):
        raise ConfigurationError(
            "unknown_policy must be 'x' or 'error'.",
            code="INVALID_UNKNOWN_CODON_POLICY",
            context={"unknown_policy": unknown_policy},
        )
    if incomplete_policy not in ("ignore", "error"):
        raise ConfigurationError(
            "incomplete_policy must be 'ignore' or 'error'.",
            code="INVALID_INCOMPLETE_CODON_POLICY",
            context={"incomplete_policy": incomplete_policy},
        )

    dna = _coerce_nucleic_text(sequence)
    if strand == "reverse":
        dna = dna.translate(_DNA_COMPLEMENT)[::-1]
    coding = dna[frame:]
    remainder = len(coding) % 3
    if remainder and incomplete_policy == "error":
        raise InvalidAlphabetError(
            "Translation input has a trailing incomplete codon.",
            code="INCOMPLETE_CODON",
            context={"frame": frame, "trailing_symbols": remainder},
        )

    protein: list[str] = []
    coding_end = len(coding) - remainder
    for offset in range(0, coding_end, 3):
        codon = coding[offset : offset + 3]
        amino_acid = _TABLE_1.get(codon)
        if amino_acid is None:
            if unknown_policy == "error":
                raise InvalidAlphabetError(
                    "Ambiguous codon cannot be translated under unknown_policy='error'.",
                    code="AMBIGUOUS_CODON",
                    context={"codon": codon, "codon_index": offset // 3},
                )
            amino_acid = "X"
        if amino_acid == "*":
            if stop_policy == "truncate":
                break
            if stop_policy == "error":
                raise InvalidAlphabetError(
                    "Stop codon encountered under stop_policy='error'.",
                    code="STOP_CODON",
                    context={"codon": codon, "codon_index": offset // 3},
                )
        protein.append(amino_acid)
    return "".join(protein)


__all__ = ["transcribe", "translate"]
