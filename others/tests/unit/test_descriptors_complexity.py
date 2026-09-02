"""Tests for linguistic complexity and exact-repeat coverage."""

import json

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap, Topology
from dnakit.descriptors import exact_repeat_fraction, linguistic_complexity
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError


def test_linguistic_complexity_distinguishes_repetitive_and_diverse_sequences() -> None:
    repetitive = linguistic_complexity(DNASequence("AAAAAAAA"), max_word_size=3)
    diverse = linguistic_complexity(DNASequence("ACGTAGCT"), max_word_size=3)

    assert 0 <= repetitive.score < diverse.score <= 1
    assert repetitive.by_k["1"] == pytest.approx(0.25)
    assert repetitive.formula.startswith("product_k")
    json.dumps(diverse.to_dict(), sort_keys=True)


def test_complexity_ambiguity_gap_and_work_limits_are_explicit() -> None:
    value = DNASequence("ANAA", alphabet=DNAAlphabet.IUPAC)
    with pytest.raises(InvalidAlphabetError):
        linguistic_complexity(value)
    ignored = linguistic_complexity(value, ambiguity_policy="ignore")
    assert ignored.score >= 0

    separated = linguistic_complexity(DNASequence(["AA", Gap(2), "AA"]), max_word_size=2)
    crossed = linguistic_complexity(
        DNASequence(["AA", Gap(2), "AA"]), max_word_size=2, cross_gaps=True
    )
    assert crossed.observation_count > separated.observation_count
    with pytest.raises(ConfigurationError) as limit:
        linguistic_complexity(DNASequence("ACGT"), max_word_size=2, max_observations=1)
    assert limit.value.code == "COMPLEXITY_OBSERVATION_LIMIT"


def test_exact_repeat_fraction_uses_union_coverage() -> None:
    repeated = exact_repeat_fraction(DNASequence("ATATATGC"), min_unit_length=2)
    none = exact_repeat_fraction(DNASequence("ACGT"), min_unit_length=2)

    assert repeated.repeat_fraction == pytest.approx(6 / 8)
    assert repeated.repeated_base_count == 6
    assert repeated.runs[0].unit == "AT"
    assert repeated.runs[0].repeat_count == 3
    assert none.repeat_fraction == 0


def test_exact_repeat_descriptor_bounds_parameters() -> None:
    with pytest.raises(ConfigurationError):
        exact_repeat_fraction(DNASequence("AAAA"), min_repeats=1)
    with pytest.raises(ConfigurationError) as limit:
        exact_repeat_fraction(DNASequence("AAAAAA"), max_comparisons=1)
    assert limit.value.code == "REPEAT_COMPARISON_LIMIT"


def test_complexity_descriptors_reject_implicit_circular_linearization() -> None:
    circular = DNASequence("ACGT", topology=Topology.CIRCULAR)
    with pytest.raises(ConfigurationError) as complexity:
        linguistic_complexity(circular)
    assert complexity.value.code == "COMPLEXITY_CIRCULAR_UNSUPPORTED"
    with pytest.raises(ConfigurationError) as repeats:
        exact_repeat_fraction(circular)
    assert repeats.value.code == "REPEAT_CIRCULAR_UNSUPPORTED"
