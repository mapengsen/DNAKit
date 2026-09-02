"""Tests for fixed-length hashed k-mer and named-panel bit fingerprints."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.fingerprints import hashed_kmer_fingerprint, panel_fingerprint
from dnakit.similarity import fingerprint_similarity


def test_hashed_kmer_fingerprint_is_fixed_length_and_deterministic() -> None:
    value = DNARecord(DNASequence("ACGTAC"), "record-1")
    result = hashed_kmer_fingerprint(
        value,
        k=3,
        n_bits=16,
        seed=7,
        representation="sparse",
    )
    repeated = hashed_kmer_fingerprint(
        value,
        k=3,
        n_bits=16,
        seed=7,
        representation="sparse",
    )

    assert result.dimension == 16
    assert result.values == {"bit:1": 1, "bit:10": 1}
    assert result.observation_count == 4
    assert result.set_bit_count == 2
    assert result.sequence_id == "record-1"
    assert result == repeated
    json.dumps(result.to_dict(), sort_keys=True)


def test_hashed_kmer_canonical_mode_is_reverse_complement_invariant() -> None:
    forward = hashed_kmer_fingerprint(DNASequence("ACGTAC"), k=3, n_bits=64)
    reverse = hashed_kmer_fingerprint(DNASequence("GTACGT"), k=3, n_bits=64)

    assert forward.dense_values() == reverse.dense_values()
    assert fingerprint_similarity(forward, reverse, metric="tanimoto").value == 1.0


def test_hashed_kmer_ambiguity_policy_is_explicit() -> None:
    value = DNASequence("ACNGT", alphabet=DNAAlphabet.IUPAC)

    ignored = hashed_kmer_fingerprint(value, k=2, n_bits=32, ambiguity_policy="ignore")

    assert ignored.observation_count == 2
    assert ignored.ignored_ambiguity_count == 1
    with pytest.raises(InvalidAlphabetError):
        hashed_kmer_fingerprint(value, k=2, n_bits=32)


def test_panel_fingerprint_has_one_interpretable_bit_per_named_pattern() -> None:
    result = panel_fingerprint(
        DNARecord(DNASequence("ATGGAATTC"), "record-1"),
        {"start": "ATG", "eco": "GAATTC", "tata": "TATA"},
    )

    assert result.feature_names == ("panel:eco", "panel:start", "panel:tata")
    assert result.dense_values() == (1, 1, 0)
    assert result.observation_count == 2
    assert result.set_bit_count == 2
    assert result.parameters["panel"] == {
        "eco": "GAATTC",
        "start": "ATG",
        "tata": "TATA",
    }
    json.dumps(result.to_dict(), sort_keys=True)


def test_panel_fingerprint_scans_both_strands_and_validates_schema_for_similarity() -> None:
    reverse_strand = panel_fingerprint(DNASequence("CAT"), {"start": "ATG"})
    different_panel = panel_fingerprint(DNASequence("CAT"), {"start": "CAT"})

    assert reverse_strand.dense_values() == (1,)
    with pytest.raises(ConfigurationError, match="schemas must match"):
        fingerprint_similarity(reverse_strand, different_panel)


@pytest.mark.parametrize(
    ("function", "kwargs"),
    [
        (hashed_kmer_fingerprint, {"k": 0}),
        (hashed_kmer_fingerprint, {"k": 2, "n_bits": 0}),
        (hashed_kmer_fingerprint, {"k": 2, "n_bits": 1_000_001}),
        (hashed_kmer_fingerprint, {"k": 2, "seed": -1}),
        (hashed_kmer_fingerprint, {"k": 2, "canonical": 1}),
        (panel_fingerprint, {"panel": {}}),
        (panel_fingerprint, {"panel": {"x": "AC"}, "mode": "regex"}),
        (panel_fingerprint, {"panel": {"x": "AC"}, "max_panel_size": 10_001}),
    ],
)
def test_bit_fingerprints_reject_invalid_configuration(
    function: object,
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ConfigurationError):
        function(DNASequence("ACGT"), **kwargs)  # type: ignore[operator]
