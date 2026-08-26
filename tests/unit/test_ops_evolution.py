"""Tests for reproducible EvoAug-inspired sequence generation."""

from collections import Counter

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.ops import (
    crossover,
    evolution_generate,
    indel_generate,
    kmer_shuffle,
    rearrange_generate,
)


def test_evolution_generation_is_reproducible_and_uses_evoaug_priority() -> None:
    source = DNASequence("ACGT" * 10)
    first = evolution_generate(
        source,
        augmentations=("mutation", "translocation", "inversion"),
        max_augmentations=3,
        seed=19,
        mut_frac=0.1,
        shift_min=1,
        shift_max=5,
        invert_min=2,
        invert_max=5,
    )
    second = evolution_generate(
        source,
        augmentations=("mutation", "translocation", "inversion"),
        max_augmentations=3,
        seed=19,
        mut_frac=0.1,
        shift_min=1,
        shift_max=5,
        invert_min=2,
        invert_max=5,
    )

    assert first == second
    assert first.augmentations == ("inversion", "translocation", "mutation")
    assert first.sequence.symbols != source.symbols
    assert len(first.sequence) == len(source)
    assert source.symbols == "ACGT" * 10
    assert first.algorithm_version == "dnakit-evoaug-v1"
    assert first.rng_state_before != first.rng_state_after


def test_evolution_indels_support_variable_length_and_evoaug_padding() -> None:
    source = DNASequence("ACGT" * 10)

    deleted = evolution_generate(
        source,
        augmentations=("deletion",),
        seed=7,
        delete_min=2,
        delete_max=2,
        pad_indels=False,
    )
    inserted = evolution_generate(
        source,
        augmentations=("insertion",),
        seed=7,
        insert_min=2,
        insert_max=5,
        pad_indels=False,
    )
    padded_delete = evolution_generate(
        source,
        augmentations=("deletion",),
        seed=7,
        delete_min=2,
        delete_max=5,
        pad_indels=True,
    )
    padded_insert = evolution_generate(
        source,
        augmentations=("insertion",),
        seed=7,
        insert_min=2,
        insert_max=5,
        pad_indels=True,
    )

    deleted_length = deleted.steps[0].length
    inserted_length = inserted.steps[0].length
    assert deleted_length is not None
    assert inserted_length is not None
    assert len(deleted.sequence) == len(source) - deleted_length
    assert len(inserted.sequence) == len(source) + inserted_length
    assert len(padded_delete.sequence) == len(source)
    assert len(padded_insert.sequence) == len(source) + 5


def test_evolution_generation_supports_soft_operation_count_and_reverse_complement() -> None:
    source = DNASequence("AACG")
    result = evolution_generate(
        source,
        augmentations=("reverse_complement", "translocation"),
        max_augmentations=2,
        hard_aug=False,
        seed=3,
        rc_prob=1.0,
        shift_min=1,
        shift_max=1,
    )

    assert 1 <= len(result.steps) <= 2
    assert set(result.augmentations) <= {"reverse_complement", "translocation"}
    if result.augmentations == ("reverse_complement",):
        assert result.sequence.symbols == "CGTT"


def test_evolution_generation_rejects_unsupported_sequence_states_and_noise() -> None:
    with pytest.raises(ConfigurationError) as missing_source:
        evolution_generate(DNASequence("ACGT"))
    assert missing_source.value.code == "EVOLUTION_RANDOM_SOURCE_REQUIRED"

    with pytest.raises(ConfigurationError) as noise:
        evolution_generate(
            DNASequence("ACGT"),
            augmentations=("noise",),  # type: ignore[arg-type]
            seed=1,
        )
    assert noise.value.code == "INVALID_EVOLUTION_AUGMENTATION"

    with pytest.raises(ConfigurationError) as ambiguous:
        evolution_generate(
            DNASequence("ACGN", alphabet=DNAAlphabet.IUPAC),
            augmentations=("mutation",),
            seed=1,
        )
    assert ambiguous.value.code == "EVOLUTION_CANONICAL_DNA_REQUIRED"

    with pytest.raises(UnsupportedGapOperationError) as gapped:
        evolution_generate(
            DNASequence.from_fragments(("AC", "GT"), (Gap(2),)),
            augmentations=("mutation",),
            seed=1,
        )
    assert gapped.value.code == "EVOLUTION_GAPPED_SEQUENCE_NOT_SUPPORTED"


def test_indel_generate_exposes_natural_variable_length_edits() -> None:
    source = DNASequence("ACGT" * 10)
    inserted = indel_generate(
        source,
        operation="insertion",
        min_length=3,
        max_length=3,
        seed=11,
    )
    deleted = indel_generate(
        source,
        operation="deletion",
        min_length=2,
        max_length=2,
        seed=11,
    )

    assert len(inserted.sequence) == len(source) + 3
    assert len(deleted.sequence) == len(source) - 2
    assert inserted == indel_generate(
        source,
        operation="insertion",
        min_length=3,
        max_length=3,
        seed=11,
    )
    assert source.symbols == "ACGT" * 10


def test_rearrange_generate_supports_exchange_inversion_and_duplication() -> None:
    source = DNASequence("AAGTCC")
    exchanged = rearrange_generate(source, operation="exchange", segment_count=3, seed=7)
    inverted = rearrange_generate(source, operation="inversion", segment_count=1, seed=7)
    duplicated = rearrange_generate(source, operation="duplication", segment_count=1, seed=7)

    assert Counter(exchanged.sequence.symbols) == Counter(source.symbols)
    assert exchanged.permutation is not None
    assert inverted.sequence.symbols == source.reverse_complement().symbols
    assert duplicated.sequence.symbols == source.symbols * 2
    assert exchanged == rearrange_generate(source, operation="exchange", segment_count=3, seed=7)


def test_kmer_shuffle_preserves_exact_overlapping_counts() -> None:
    source = DNASequence("AAGATCGATCGGATC")

    first = kmer_shuffle(source, k=2, seed=23)
    second = kmer_shuffle(source, k=2, seed=23)
    for k in (1, 2, 3):
        result = kmer_shuffle(source, k=k, seed=23)
        expected = Counter(
            source.symbols[index : index + k] for index in range(len(source) - k + 1)
        )
        observed = Counter(
            result.sequence.symbols[index : index + k]
            for index in range(len(result.sequence) - k + 1)
        )
        assert observed == expected
        assert result.sequence.symbols != source.symbols
    assert first == second


def test_kmer_shuffle_reports_unique_reconstruction() -> None:
    with pytest.raises(ConfigurationError) as error:
        kmer_shuffle(DNASequence("AAAA"), k=2, seed=1, max_attempts=2)
    assert error.value.code == "KMER_SHUFFLE_NO_ALTERNATIVE"

    unchanged = kmer_shuffle(
        DNASequence("AAAA"),
        k=2,
        seed=1,
        ensure_different=False,
    )
    assert unchanged.sequence.symbols == "AAAA"


def test_crossover_combines_two_equal_length_parents() -> None:
    first = DNASequence("AAAACCCC")
    second = DNASequence("GGGGTTTT")
    specified = crossover(first, second, position=4)
    randomized = crossover(first, second, seed=31)

    assert specified.sequence.symbols == "AAAATTTT"
    assert specified.position == 4
    assert specified.random_source == "none"
    assert randomized == crossover(first, second, seed=31)
    assert (
        randomized.sequence.symbols[: randomized.position] == first.symbols[: randomized.position]
    )
    assert (
        randomized.sequence.symbols[randomized.position :] == second.symbols[randomized.position :]
    )


def test_crossover_rejects_unequal_parent_lengths() -> None:
    with pytest.raises(ConfigurationError) as error:
        crossover(DNASequence("AAAA"), DNASequence("CCCCCC"), position=2)
    assert error.value.code == "CROSSOVER_LENGTH_MISMATCH"
