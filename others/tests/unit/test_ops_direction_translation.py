"""Human-verifiable tests for direction, transcription, and translation."""

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError, UnsupportedGapOperationError
from dnakit.ops import complement, reverse, reverse_complement, transcribe, translate


def test_direction_functions_delegate_without_mutating_input() -> None:
    gap = Gap(2)
    sequence = DNASequence(["ARY", gap, "GC"], alphabet=DNAAlphabet.IUPAC)

    assert reverse(sequence).parts == ("CG", gap, "YRA")
    assert complement(sequence).parts == ("TYR", gap, "CG")
    assert reverse_complement(sequence).parts == ("GC", gap, "RYT")
    assert sequence.parts == ("ARY", gap, "GC")


def test_transcription_is_strand_explicit_and_rejects_gap_omission() -> None:
    sequence = DNASequence("ATGCC")

    assert transcribe(sequence) == "AUGCC"
    assert transcribe(sequence, strand="reverse") == "GGCAU"

    with pytest.raises(ConfigurationError) as error:
        transcribe(sequence, strand="both")  # type: ignore[arg-type]
    assert error.value.code == "INVALID_STRAND_POLICY"
    with pytest.raises(UnsupportedGapOperationError):
        transcribe(DNASequence(["ATG", Gap(2), "CCC"]))


def test_table_one_translation_stop_and_frame_policies() -> None:
    sequence = DNASequence("AATGGGCTAAC")

    assert translate(sequence, frame=0) == "NGL"
    assert translate(sequence, frame=1) == "MG*"
    assert translate(sequence, frame=1, stop_policy="truncate") == "MG"
    assert translate("AUGGCCUAA") == "MA*"

    with pytest.raises(InvalidAlphabetError) as error:
        translate(sequence, frame=1, stop_policy="error")
    assert error.value.code == "STOP_CODON"


def test_reverse_translation_applies_reverse_complement_before_frame() -> None:
    # RC(AAACAT) == ATGTTT; frame 0 therefore translates to MF.
    assert translate(DNASequence("AAACAT"), strand="reverse", frame=0) == "MF"
    # RC(AAACATC) == GATGTTT; frame 1 selects ATGTTT.
    assert translate(DNASequence("AAACATC"), strand="reverse", frame=1) == "MF"


def test_translation_ambiguity_partial_codon_and_configuration_are_explicit() -> None:
    ambiguous = DNASequence("ATGNNN", alphabet=DNAAlphabet.IUPAC)

    assert translate(ambiguous) == "MX"
    with pytest.raises(InvalidAlphabetError) as ambiguous_error:
        translate(ambiguous, unknown_policy="error")
    assert ambiguous_error.value.code == "AMBIGUOUS_CODON"
    with pytest.raises(InvalidAlphabetError) as partial_error:
        translate("ATGA", incomplete_policy="error")
    assert partial_error.value.code == "INCOMPLETE_CODON"
    with pytest.raises(ConfigurationError) as table_error:
        translate("ATG", table=2)
    assert table_error.value.code == "UNSUPPORTED_GENETIC_CODE"
    with pytest.raises(ConfigurationError):
        translate("ATG", stop_policy="bad")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        translate("ATG", unknown_policy="bad")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        translate("ATG", incomplete_policy="bad")  # type: ignore[arg-type]
    with pytest.raises(InvalidAlphabetError):
        translate("atg")
    with pytest.raises(UnsupportedGapOperationError):
        translate(DNASequence(["ATG", Gap(1), "TAA"]))
