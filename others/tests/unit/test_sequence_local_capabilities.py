"""Regression tests for the dependency-free formal sequence APIs."""

import pytest

from dnakit.sequence import (
    calculate_sequence_metrics,
    count_kmers,
    deduplicate_sequences,
    reverse_complement_dna,
    sequence_length_gc_summary,
)


def test_kmer_counts_are_overlapping_and_sorted() -> None:
    assert count_kmers("ACACA", k=2) == {"AC": 2, "CA": 2}


def test_metrics_include_ambiguity_without_inflating_gc() -> None:
    metrics = calculate_sequence_metrics(" acgtn ")
    assert metrics["length"] == 5
    assert metrics["gc_percent"] == pytest.approx(40.0)
    assert metrics["n_count"] == 1
    assert metrics["ambiguous_count"] == 1


def test_deduplication_preserves_first_normalized_sequence() -> None:
    assert deduplicate_sequences(["acgt", " ACGT ", "TGCA"]) == ["ACGT", "TGCA"]
    assert deduplicate_sequences(["ACG", "CGT"], canonical=True) == ["ACG"]


def test_reverse_complement_supports_iupac_symbols() -> None:
    assert reverse_complement_dna("ACGTRN") == "NYACGT"


def test_length_gc_summary_handles_empty_and_ambiguous_sequences() -> None:
    assert sequence_length_gc_summary([])["count"] == 0
    summary = sequence_length_gc_summary(["ACGT", "GGNN"])
    assert summary["lengths"] == [4, 4]
    assert summary["mean_gc_fraction"] == pytest.approx(0.5)


def test_length_gc_summary_rejects_non_dna() -> None:
    with pytest.raises(ValueError, match="invalid DNA symbols"):
        sequence_length_gc_summary(["ACX"])
