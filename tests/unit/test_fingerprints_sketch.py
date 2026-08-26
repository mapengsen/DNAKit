"""Tests for deterministic native sketches."""

import json

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.fingerprints import fracminhash, minhash


def test_minhash_is_deterministic_bounded_and_auditable() -> None:
    value = DNARecord(DNASequence("ACGTACGT"), "record")

    first = minhash(value, k=3, num_hashes=3, seed=7)
    second = minhash(value, k=3, num_hashes=3, seed=7)

    assert first == second
    assert first.sequence_id == "record"
    assert len(first.hashes) == 2
    assert first.hashes == tuple(sorted(first.hashes))
    assert first.observation_count == 6
    assert first.unique_hash_count == 2
    json.dumps(first.to_dict(), sort_keys=True)


def test_canonical_minhash_matches_reverse_complement() -> None:
    forward = DNASequence("AACG")
    reverse = forward.reverse_complement()

    assert (
        minhash(forward, k=2, num_hashes=100).hashes == minhash(reverse, k=2, num_hashes=100).hashes
    )
    assert (
        minhash(forward, k=2, num_hashes=100, canonical=False).hashes
        != minhash(reverse, k=2, num_hashes=100, canonical=False).hashes
    )


def test_fracminhash_threshold_selection_and_schema() -> None:
    all_hashes = fracminhash(DNASequence("ACGTACGT"), k=2, scaled=1)
    sparse = fracminhash(DNASequence("ACGTACGT"), k=2, scaled=2)

    assert all_hashes.hashes
    assert set(sparse.hashes) <= set(all_hashes.hashes)
    assert sparse.selection == "scaled"
    assert sparse.threshold == 2**63
    assert sparse.schema_version == "dnakit.sketch.fracminhash.v1"


def test_sketch_rejects_ambiguity_gaps_and_invalid_limits() -> None:
    with pytest.raises(ConfigurationError) as ambiguity:
        minhash(DNASequence("AN", alphabet=DNAAlphabet.IUPAC), k=1)
    assert ambiguity.value.code == "SKETCH_AMBIGUITY_NOT_ALLOWED"
    with pytest.raises(UnsupportedGapOperationError) as gap:
        minhash(DNASequence(["A", Gap(2), "C"]), k=1)
    assert gap.value.code == "SKETCH_GAP_NOT_ALLOWED"
    with pytest.raises(ConfigurationError):
        minhash(DNASequence("AC"), k=1, seed=True)
    with pytest.raises(ConfigurationError) as limit:
        minhash(DNASequence("AC"), k=1, num_hashes=2, max_hashes=1)
    assert limit.value.code == "SKETCH_HASH_LIMIT_EXCEEDED"
    with pytest.raises(ConfigurationError) as unique_limit:
        minhash(
            DNASequence("ACGT"),
            k=1,
            num_hashes=2,
            canonical=False,
            max_unique_hashes=2,
        )
    assert unique_limit.value.code == "SKETCH_UNIQUE_HASH_LIMIT_EXCEEDED"
    with pytest.raises(ConfigurationError) as circular:
        minhash(DNASequence("ACGT", topology=Topology.CIRCULAR), k=2)
    assert circular.value.code == "SKETCH_CIRCULAR_UNSUPPORTED"
