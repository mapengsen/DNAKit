"""Human-verifiable tests for fixed-schema exact k-mer fingerprints."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.fingerprints import (
    FingerprintRepresentation,
    KmerFingerprintMode,
    kmer,
    kmer_fingerprint,
)


def test_dense_count_fingerprint_has_complete_lexicographic_schema() -> None:
    result = kmer(DNARecord(DNASequence("ACGT"), "record-1"), k=2)

    assert result.dimension == 16
    assert result.feature_names[:5] == ("AA", "AC", "AG", "AT", "CA")
    assert result.feature_names[-1] == "TT"
    assert result.dense_values()[result.feature_names.index("AC")] == 1
    assert result.dense_values()[result.feature_names.index("CG")] == 1
    assert result.dense_values()[result.feature_names.index("GT")] == 1
    assert sum(result.dense_values()) == 3
    assert result.observation_count == 3
    assert result.sequence_id == "record-1"
    assert result.schema_version == "dnakit.kmer.acgt.v1"


def test_sparse_binary_and_presence_alias_preserve_the_same_schema() -> None:
    binary = kmer(DNASequence("AAAA"), k=2, mode="binary", representation="sparse")
    presence = kmer_fingerprint(DNASequence("AAAA"), k=2, mode="presence")

    assert binary.mode is KmerFingerprintMode.BINARY
    assert binary.representation is FingerprintRepresentation.SPARSE
    assert binary.values == {"AA": 1}
    assert binary.sparse_values() == {"AA": 1}
    assert binary.dense_values()[0] == 1
    assert binary.dimension == presence.dimension == 16
    assert presence.mode is KmerFingerprintMode.BINARY


def test_frequency_uses_valid_observation_count_as_denominator() -> None:
    result = kmer(DNASequence("AAAA"), k=2, mode="frequency")

    assert result.observation_count == 3
    assert result.dense_values()[0] == 1.0
    assert sum(result.dense_values()) == 1.0


def test_canonical_schema_and_counts_are_reverse_complement_collapsed() -> None:
    result = kmer(DNASequence("ACGT"), k=2, canonical=True, representation="sparse")

    assert result.dimension == 10
    assert result.feature_names == tuple(sorted(result.feature_names))
    assert result.values == {"AC": 2, "CG": 1}
    assert "GT" not in result.feature_names


def test_empty_and_k_larger_than_sequence_return_fixed_all_zero_vectors() -> None:
    empty = kmer(DNASequence(""), k=2)
    too_large = kmer(DNASequence("A"), k=2, representation="sparse")

    assert empty.dimension == 16
    assert empty.observation_count == 0
    assert empty.dense_values() == (0,) * 16
    assert too_large.values == {}
    assert too_large.dense_values() == (0,) * 16


def test_ambiguity_and_gap_traversal_are_explicit() -> None:
    ambiguous = DNASequence("ACNGT", alphabet=DNAAlphabet.IUPAC)
    gapped = DNASequence(["AC", Gap(5), "GT"])

    ignored = kmer(ambiguous, k=2, ambiguity_policy="ignore", representation="sparse")
    split = kmer(gapped, k=2, representation="sparse")
    crossed = kmer(gapped, k=2, cross_gaps=True, representation="sparse")

    assert ignored.values == {"AC": 1, "GT": 1}
    assert ignored.ignored_ambiguity_count == 1
    assert split.values == {"AC": 1, "GT": 1}
    assert crossed.values == {"AC": 1, "CG": 1, "GT": 1}
    with pytest.raises(InvalidAlphabetError):
        kmer(ambiguous, k=2)


def test_fingerprint_never_bridges_a_non_crossable_gap() -> None:
    sequence = DNASequence(["AC", Gap(1, crossable=False), "GT"])

    result = kmer(sequence, k=2, cross_gaps=True, representation="sparse")

    assert result.values == {"AC": 1, "GT": 1}


@pytest.mark.parametrize(
    "kwargs",
    [
        {"k": 0},
        {"k": True},
        {"k": 2, "canonical": 1},
        {"k": 2, "overlapping": 1},
        {"k": 2, "cross_gaps": 0},
        {"k": 2, "mode": "weighted"},
        {"k": 2, "representation": "array"},
        {"k": 2, "ambiguity_policy": "fractional"},
        {"k": 2, "max_dimension": 0},
        {"k": 3, "max_dimension": 63},
    ],
)
def test_kmer_fingerprint_rejects_invalid_or_unsafe_configuration(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ConfigurationError):
        kmer(DNASequence("ACGT"), **kwargs)  # type: ignore[arg-type]


def test_kmer_result_is_json_serializable_with_order_and_parameters() -> None:
    result = kmer(DNASequence("ACGT"), k=2, canonical=True, representation="sparse")
    payload = json.loads(json.dumps(result.to_dict()))

    assert payload["feature_names"] == list(result.feature_names)
    assert payload["values"] == {"AC": 2, "CG": 1}
    assert payload["canonical"] is True
    assert payload["representation"] == "sparse"
