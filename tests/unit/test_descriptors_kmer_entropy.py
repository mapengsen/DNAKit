"""Human-verifiable tests for DESC-006 and DESC-007."""

from __future__ import annotations

import math

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap
from dnakit.descriptors import canonical_kmer, kmer_statistics, shannon_entropy
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError


def test_overlapping_count_frequency_presence_and_large_k() -> None:
    result = kmer_statistics(DNASequence("ACGT"), 2)
    too_large = kmer_statistics(DNASequence("ACGT"), 5)

    assert result.counts == {"AC": 1, "CG": 1, "GT": 1}
    assert result.frequencies == {"AC": 1 / 3, "CG": 1 / 3, "GT": 1 / 3}
    assert result.presence == ("AC", "CG", "GT")
    assert result.denominator == 3
    assert too_large.counts == {}
    assert too_large.frequencies == {}
    assert too_large.presence == ()
    assert too_large.denominator == 0


def test_non_overlapping_and_canonical_kmers_are_explicit() -> None:
    non_overlapping = kmer_statistics(DNASequence("ACGTAC"), 2, overlapping=False)
    canonical = kmer_statistics(DNASequence("ACGT"), 2, canonical=True)

    assert non_overlapping.counts == {"AC": 2, "GT": 1}
    assert canonical.counts == {"AC": 2, "CG": 1}
    assert canonical_kmer("GT") == "AC"
    with pytest.raises(ConfigurationError):
        canonical_kmer("AN")
    with pytest.raises(ConfigurationError):
        canonical_kmer("")


def test_kmers_do_not_bridge_gap_or_ambiguity_by_default() -> None:
    gapped = DNASequence(["AC", Gap(5), "GT"])
    ambiguous = DNASequence("ACNGT", alphabet=DNAAlphabet.IUPAC)

    assert kmer_statistics(gapped, 2).counts == {"AC": 1, "GT": 1}
    assert kmer_statistics(gapped, 2, cross_gaps=True).counts == {
        "AC": 1,
        "CG": 1,
        "GT": 1,
    }
    ignored = kmer_statistics(ambiguous, 2, ambiguity_policy="ignore")
    assert ignored.counts == {"AC": 1, "GT": 1}
    assert ignored.ignored_ambiguity_count == 1
    with pytest.raises(InvalidAlphabetError):
        kmer_statistics(ambiguous, 2)


def test_kmers_never_bridge_a_non_crossable_gap() -> None:
    sequence = DNASequence(["AC", Gap(1, crossable=False), "GT"])

    result = kmer_statistics(sequence, 2, cross_gaps=True)

    assert result.counts == {"AC": 1, "GT": 1}


def test_shannon_base_and_kmer_entropy_are_hand_checkable() -> None:
    uniform = shannon_entropy(DNASequence("ACGT"))
    constant = shannon_entropy(DNASequence("AAAA"))
    kmer = shannon_entropy(DNASequence("ACAC"), unit="kmer", k=2)
    empty = shannon_entropy(DNASequence(""))

    assert uniform.entropy == 2.0
    assert constant.entropy == 0.0
    assert math.isclose(kmer.entropy, 0.9182958340544896)
    assert kmer.observation_count == 3
    assert kmer.category_count == 2
    assert empty.entropy == 0.0
    assert empty.observation_count == 0


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"unit": "word"}, "unit"),
        ({"log_base": 1.0}, "log_base"),
        ({"k": 0, "unit": "kmer"}, "positive integer"),
        ({"cross_gaps": 1}, "cross_gaps"),
    ],
)
def test_entropy_rejects_invalid_configuration(
    kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ConfigurationError, match=message):
        shannon_entropy(DNASequence("ACGT"), **kwargs)  # type: ignore[arg-type]
