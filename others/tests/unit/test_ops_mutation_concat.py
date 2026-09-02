"""Tests for reproducible mutation and sequence concatenation."""

import random

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.ops import concat, concat_overlap, mutate


def test_specified_mutation_is_one_based_neither_implicitly_nor_combinatorial() -> None:
    source = DNASequence("ACGT")

    result = mutate(source, position=1, replacement="T")
    assert result.sequence.symbols == "ATGT"
    assert result.mode == "specified"
    assert result.edit.start == 1
    assert result.edit.removed_symbols == "C"
    assert result.seed is None


def test_random_mutation_requires_explicit_source_and_is_reproducible() -> None:
    source = DNASequence("AAAA")

    first = mutate(source, seed=19)
    second = mutate(source, seed=19)
    assert first == second
    assert first.sequence != source
    assert first.seed == 19

    generator = random.Random(7)
    from_rng = mutate(source, rng=generator, allowed_bases="CT")
    assert from_rng.edit.replacement_symbols in "CT"
    assert from_rng.random_source == "rng"
    assert from_rng.rng_state_before is not None
    assert from_rng.rng_state_after is not None
    replay = random.Random()
    replay.setstate(from_rng.rng_state_before)
    replayed = mutate(source, rng=replay, allowed_bases="CT")
    assert replayed.sequence == from_rng.sequence
    assert replayed.rng_state_after == from_rng.rng_state_after
    with pytest.raises(ConfigurationError) as missing:
        mutate(source)
    assert missing.value.code == "RANDOM_SOURCE_REQUIRED"
    with pytest.raises(ConfigurationError):
        mutate(source, seed=1, rng=random.Random(1))
    with pytest.raises(ConfigurationError):
        mutate(source, position=1)
    with pytest.raises(ConfigurationError):
        mutate(source, position=1, replacement="TT")
    with pytest.raises(ConfigurationError) as no_op:
        mutate(source, position=1, replacement="A")
    assert no_op.value.code == "NO_OP_MUTATION"


def test_random_mutation_can_resolve_iupac_symbol_to_allowed_canonical_base() -> None:
    source = DNASequence("N", alphabet=DNAAlphabet.IUPAC)

    result = mutate(source, seed=3, allowed_bases="AC")
    assert result.sequence.symbols in {"A", "C"}
    assert result.edit.removed_symbols == "N"


def test_concat_with_linker_or_gap_preserves_input_and_promotes_alphabet() -> None:
    first = DNASequence("AA")
    second = DNASequence("CN", alphabet=DNAAlphabet.IUPAC)

    linked = concat([first, second], linker="T")
    assert linked.symbols == "AATCN"
    assert linked.alphabet is DNAAlphabet.IUPAC
    gap = Gap(None, crossable=False)
    gapped = concat([first, DNASequence("CC"), DNASequence("GG")], gap=gap)
    assert gapped.parts == ("AA", gap, "CC", gap, "GG")
    assert first.symbols == "AA"


def test_concat_rejects_ambiguous_separator_and_circular_or_mismatched_inputs() -> None:
    with pytest.raises(ConfigurationError):
        concat([DNASequence("AA")])
    with pytest.raises(ConfigurationError):
        concat([DNASequence("AA"), DNASequence("CC")], linker="T", gap=Gap(1))
    with pytest.raises(ConfigurationError) as circular:
        concat([DNASequence("AA", topology=Topology.CIRCULAR), DNASequence("CC")])
    assert circular.value.code == "CIRCULAR_CONCAT_NOT_SUPPORTED"
    with pytest.raises(ConfigurationError) as later_circular:
        concat([DNASequence("AA"), DNASequence("CC", topology=Topology.CIRCULAR)])
    assert later_circular.value.code == "CIRCULAR_FRAGMENT_NOT_SUPPORTED"
    with pytest.raises(ConfigurationError) as circular_linker:
        concat(
            [DNASequence("AA"), DNASequence("CC")],
            linker=DNASequence("T", topology=Topology.CIRCULAR),
        )
    assert circular_linker.value.code == "CIRCULAR_FRAGMENT_NOT_SUPPORTED"


def test_concat_overlap_keeps_the_longest_exact_junction_once() -> None:
    first = DNASequence("AAACCC")
    second = DNASequence("CCCGN", alphabet=DNAAlphabet.IUPAC)

    joined = concat_overlap([first, second])

    assert joined.symbols == "AAACCCGN"
    assert joined.alphabet is DNAAlphabet.IUPAC
    assert first.symbols == "AAACCC"
    assert second.symbols == "CCCGN"


def test_concat_overlap_respects_overlap_bounds_and_requires_two_fragments() -> None:
    assert concat_overlap(["AAACCC", "CCCGG"], max_overlap=2).symbols == "AAACCCCGG"
    with pytest.raises(ConfigurationError) as no_overlap:
        concat_overlap(["AACC", "CCGG"], min_overlap=3, max_overlap=3)
    assert no_overlap.value.code == "OVERLAP_NOT_FOUND"

    with pytest.raises(ConfigurationError):
        concat_overlap(["AA"])
    with pytest.raises(ConfigurationError):
        concat_overlap(["AA", "AC", "CG"])
    with pytest.raises(ConfigurationError):
        concat_overlap(["AA", "AC"], min_overlap=3, max_overlap=2)


def test_concat_overlap_rejects_gapped_or_circular_inputs() -> None:
    with pytest.raises(UnsupportedGapOperationError) as gapped:
        concat_overlap([DNASequence(["AA", Gap(1)]), DNASequence("AACC")])
    assert gapped.value.code == "OVERLAP_CONCAT_GAPPED_INPUT"

    with pytest.raises(ConfigurationError) as circular:
        concat_overlap([DNASequence("AA", topology=Topology.CIRCULAR), DNASequence("AACC")])
    assert circular.value.code == "CIRCULAR_CONCAT_NOT_SUPPORTED"
