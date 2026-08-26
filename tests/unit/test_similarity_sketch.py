"""Tests for sketch similarity."""

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints import fracminhash, minhash
from dnakit.similarity import sketch_similarity


def test_sketch_jaccard_and_directional_containment() -> None:
    left = minhash(DNARecord(DNASequence("ACGTAC"), "left"), k=2, num_hashes=100)
    right = minhash(DNARecord(DNASequence("ACGTTC"), "right"), k=2, num_hashes=100)

    jaccard = sketch_similarity(left, right, min_shared_hashes=1)
    containment = sketch_similarity(left, right, metric="containment")

    assert jaccard.left_id == "left"
    assert jaccard.right_id == "right"
    assert jaccard.value == pytest.approx(jaccard.shared_hash_count / jaccard.union_hash_count)
    assert containment.value == pytest.approx(
        containment.shared_hash_count / containment.left_hash_count
    )
    assert jaccard.passed_min_shared_hashes


def test_empty_sketch_semantics_are_explicit() -> None:
    empty = minhash(DNASequence("A"), k=2)
    nonempty = minhash(DNASequence("AC"), k=2)

    assert sketch_similarity(empty, empty).value == 1
    assert sketch_similarity(empty, nonempty).value == 0
    assert sketch_similarity(empty, nonempty, metric="containment").value == 0


def test_incompatible_sketches_and_invalid_parameters_are_rejected() -> None:
    left = minhash(DNASequence("ACGT"), k=2, num_hashes=10)
    with pytest.raises(ConfigurationError) as schema:
        sketch_similarity(left, minhash(DNASequence("ACGT"), k=3, num_hashes=10))
    assert schema.value.code == "SKETCH_SCHEMA_MISMATCH"
    with pytest.raises(ConfigurationError):
        sketch_similarity(left, fracminhash(DNASequence("ACGT"), k=2, scaled=1))
    with pytest.raises(ConfigurationError):
        sketch_similarity(left, left, min_shared_hashes=True)
