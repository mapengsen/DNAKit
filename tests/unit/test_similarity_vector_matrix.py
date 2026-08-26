"""Tests for k-mer, fingerprint, unified comparison, and pairwise matrices."""

import json
import math
from collections.abc import Iterator

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.fingerprints import kmer as kmer_fingerprint
from dnakit.similarity import (
    compare,
    exact_similarity,
    fingerprint_similarity,
    kmer_similarity,
    kmer_vector_similarity,
    similarity_matrix,
)


def test_kmer_jaccard_uses_literal_iupac_symbols() -> None:
    left = DNASequence("AA")
    right = DNASequence("AN", alphabet=DNAAlphabet.IUPAC)

    result = kmer_similarity(left, right, k=1)

    assert result.value == pytest.approx(0.5)
    assert result.iupac_matching == "literal"
    assert result.parameters["mode"] == "set"
    assert result.components["shared_weight"] == 1


def test_kmer_modes_metrics_and_canonical_reverse_complement() -> None:
    set_result = kmer_similarity(DNASequence("AAA"), DNASequence("AAC"), k=1)
    count_result = kmer_similarity(
        DNASequence("AAA"),
        DNASequence("AAC"),
        k=1,
        mode="count",
    )
    containment = kmer_similarity(
        DNASequence("A"),
        DNASequence("AC"),
        k=1,
        metric="containment",
    )
    canonical = kmer_similarity(
        DNASequence("ARY", alphabet=DNAAlphabet.IUPAC),
        DNASequence("RYT", alphabet=DNAAlphabet.IUPAC),
        k=2,
        canonical=True,
    )

    assert set_result.value == pytest.approx(0.5)
    assert count_result.value == pytest.approx(0.5)
    assert containment.value == 1
    assert canonical.value == 1


def test_kmer_empty_boundaries_and_gap_rejection_are_explicit() -> None:
    both_empty = kmer_similarity(DNASequence("A"), DNASequence("T"), k=2)
    one_empty = kmer_similarity(DNASequence("A"), DNASequence("TT"), k=2)

    assert both_empty.value == 1
    assert one_empty.value == 0
    with pytest.raises(ConfigurationError) as zero_error:
        kmer_similarity(
            DNASequence(""),
            DNASequence(""),
            k=1,
            zero_vector_policy="error",
        )
    assert zero_error.value.code == "ZERO_VECTOR_UNDEFINED"
    with pytest.raises(UnsupportedGapOperationError):
        kmer_similarity(DNASequence(["A", Gap(1), "C"]), DNASequence("AC"), k=1)


def test_named_kmer_vectors_support_directional_count_containment() -> None:
    direct = kmer_vector_similarity(
        {"AA": 2, "AC": 1},
        {"AA": 1, "AC": 1, "TT": 3},
        metric="containment",
        mode="count",
    )
    dispatched = compare(
        {"AA": 2, "AC": 1},
        {"AA": 1, "AC": 1, "TT": 3},
        method="kmer_containment",
        mode="count",
    )

    assert direct.value == pytest.approx(2 / 3)
    assert dispatched.value == pytest.approx(2 / 3)  # type: ignore[union-attr]


def test_kmer_float_accumulation_is_stable_across_mapping_order() -> None:
    feature_names = [f"K{index:04d}" for index in range(1_000)]
    left_items = [("HUGE", 1e16), *((name, 1.0) for name in feature_names)]
    right_items = [("HUGE", 1e16), *((name, 1.0) for name in feature_names[:500])]

    forward = kmer_vector_similarity(
        dict(left_items),
        dict(right_items),
        metric="jaccard",
        mode="count",
    )
    reversed_order = kmer_vector_similarity(
        dict(reversed(left_items)),
        dict(reversed(right_items)),
        metric="jaccard",
        mode="count",
    )

    assert forward.to_dict() == reversed_order.to_dict()
    assert forward.components["left_weight"] == math.fsum((1e16, *([1.0] * 1_000)))
    assert forward.components["right_weight"] == math.fsum((1e16, *([1.0] * 500)))


def test_numeric_fingerprint_metrics_have_declared_formulas() -> None:
    left = (1.0, 1.0)
    right = (1.0, 0.0)

    assert fingerprint_similarity(left, right, metric="tanimoto").value == pytest.approx(0.5)
    assert fingerprint_similarity(left, right, metric="jaccard").value == pytest.approx(0.5)
    assert fingerprint_similarity(left, right, metric="cosine").value == pytest.approx(
        1 / math.sqrt(2)
    )
    assert fingerprint_similarity(left, right, metric="euclidean").value == 1
    assert fingerprint_similarity(left, right, metric="manhattan").value == 1
    assert fingerprint_similarity((1.0, -1.0), (-1.0, 1.0), metric="cosine").value == -1


def test_zero_vectors_negative_tanimoto_and_weighting_are_explicit() -> None:
    assert fingerprint_similarity((0, 0), (0, 0), metric="cosine").value == 1
    assert fingerprint_similarity((0, 0), (1, 0), metric="cosine").value == 0
    assert fingerprint_similarity({}, {"A": 1}, metric="tanimoto").value == 0
    weighted = fingerprint_similarity(
        {"A": 1, "B": 1},
        {"A": 1, "B": 0},
        metric="jaccard",
        weights={"A": 2, "B": 1},
    )
    assert weighted.value == pytest.approx(2 / 3)
    with pytest.raises(ConfigurationError) as zero_error:
        fingerprint_similarity((0, 0), (0, 0), metric="cosine", zero_vector_policy="error")
    assert zero_error.value.code == "ZERO_VECTOR_UNDEFINED"
    with pytest.raises(ConfigurationError) as negative_error:
        fingerprint_similarity((-1, 0), (1, 0), metric="tanimoto")
    assert negative_error.value.code == "NEGATIVE_VECTOR_NOT_ALLOWED"


def test_scalar_weight_audit_records_resolved_named_mapping_json_stably() -> None:
    result = fingerprint_similarity(
        {"B": 1, "A": 1},
        {"A": 1, "B": 0},
        metric="jaccard",
        weights={"B": 2.25},
    )

    assert result.parameters["weighted"] is True
    assert result.parameters["weights"] == {"A": 1.0, "B": 2.25}
    payload = result.to_dict()
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    decoded = json.loads(encoded)

    assert decoded["parameters"]["weights"] == {"A": 1.0, "B": 2.25}
    assert json.dumps(decoded, sort_keys=True, separators=(",", ":")) == encoded

    positional = fingerprint_similarity(
        (1, 0),
        (1, 1),
        metric="cosine",
        weights=(0.5, 4),
    ).to_dict()
    assert positional["parameters"]["weights"] == [0.5, 4.0]


def test_versioned_fingerprint_schema_and_gap_provenance_are_preserved() -> None:
    first = kmer_fingerprint(DNARecord(DNASequence("AAC"), "first"), k=1)
    second = kmer_fingerprint(DNARecord(DNASequence("ACC"), "second"), k=1)
    compared = fingerprint_similarity(first, second)

    assert compared.left_id == "first"
    assert compared.right_id == "second"
    with pytest.raises(ConfigurationError) as schema_error:
        fingerprint_similarity(first, kmer_fingerprint(DNASequence("AAC"), k=2))
    assert schema_error.value.code == "FINGERPRINT_SCHEMA_MISMATCH"

    gapped = kmer_fingerprint(DNASequence(["A", Gap(2), "C"]), k=1)
    gapped_comparison = fingerprint_similarity(gapped, gapped)
    assert gapped.gap_count == 1
    assert gapped_comparison.value == 1


def test_compare_dispatches_exact_distance_kmer_and_vector_methods() -> None:
    exact = exact_similarity(DNASequence("ATG"), DNASequence("CAT"), reverse_complement=True)
    hamming = compare(DNASequence("AC"), DNASequence("AT"), method="hamming")
    kmer = compare(DNASequence("AA"), DNASequence("AC"), method="kmer_jaccard", k=1)
    vector = compare((1, 0), (1, 1), method="cosine")

    assert exact.value == 1
    assert hamming.distance == 1  # type: ignore[union-attr]
    assert kmer.value == pytest.approx(0.5)  # type: ignore[union-attr]
    assert vector.value == pytest.approx(1 / math.sqrt(2))  # type: ignore[union-attr]
    with pytest.raises(ConfigurationError):
        compare(DNASequence("A"), DNASequence("A"), method="kmer_jaccard")
    with pytest.raises(ConfigurationError):
        compare(DNASequence("A"), DNASequence("A"), method="unknown")


def test_pairwise_exact_and_distance_matrices_preserve_labels_and_order() -> None:
    records = [
        DNARecord(DNASequence("AA"), "a"),
        DNARecord(DNASequence("AA"), "b"),
        DNARecord(DNASequence("AT"), "c"),
    ]

    exact = similarity_matrix(records, method="exact")
    hamming = similarity_matrix(records, method="hamming")

    assert exact.labels == ("a", "b", "c")
    assert exact.values == ((1.0, 1.0, 0.0), (1.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    assert hamming.value_kind == "distance"
    assert hamming.values == ((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), (1.0, 1.0, 0.0))
    assert exact.symmetric and hamming.symmetric


def test_matrix_weight_audit_records_positional_values_through_json_roundtrip() -> None:
    result = similarity_matrix(
        [(1, 0), (0, 1)],
        method="cosine",
        weights=(0.125, 3.75),
    )

    assert result.parameters["weighted"] is True
    assert result.parameters["weights"] == (0.125, 3.75)
    payload = result.to_dict()
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    decoded = json.loads(encoded)

    assert decoded["parameters"]["weights"] == [0.125, 3.75]
    assert json.dumps(decoded, sort_keys=True, separators=(",", ":")) == encoded


def test_matrix_weight_audit_merges_resolved_named_feature_defaults() -> None:
    result = similarity_matrix(
        [{"A": 1}, {"A": 1, "B": 1}],
        method="jaccard",
        weights={"A": 2},
    )

    assert result.parameters["weights"] == {"A": 2.0, "B": 1.0}
    encoded = json.dumps(
        result.to_dict(),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    assert json.loads(encoded)["parameters"]["weights"] == {"A": 2.0, "B": 1.0}


def test_weighted_matrix_requires_an_auditable_vector_schema() -> None:
    with pytest.raises(ConfigurationError) as empty_error:
        similarity_matrix([], method="cosine", weights=(1, 2))
    assert empty_error.value.code == "EMPTY_WEIGHTED_MATRIX"

    with pytest.raises(ConfigurationError) as ignored_error:
        similarity_matrix([DNASequence("A")], method="exact", weights=(1,))
    assert ignored_error.value.code == "UNSUPPORTED_MATRIX_WEIGHTS"

    with pytest.raises(ConfigurationError) as dimension_error:
        similarity_matrix([(1, 0)], method="cosine", weights=(1,))
    assert dimension_error.value.code == "WEIGHT_DIMENSION_MISMATCH"

    with pytest.raises(ConfigurationError) as key_error:
        fingerprint_similarity(
            {"A": 1},
            {"A": 1},
            weights={1: 2, "A": 1},  # type: ignore[dict-item]
        )
    assert key_error.value.code == "INVALID_VECTOR_WEIGHTS"


def test_kmer_containment_matrix_is_directional_and_matrix_size_is_bounded() -> None:
    result = similarity_matrix(
        [DNASequence("A"), DNASequence("AC")],
        method="kmer_containment",
        k=1,
    )

    assert not result.symmetric
    assert result.values == ((1.0, 1.0), (0.5, 1.0))
    assert similarity_matrix([], method="exact").values == ()
    with pytest.raises(ConfigurationError) as size_error:
        similarity_matrix(
            [DNASequence("A"), DNASequence("C")],
            method="exact",
            max_items=1,
        )
    assert size_error.value.code == "SIMILARITY_MATRIX_SIZE_LIMIT"
    with pytest.raises(ConfigurationError):
        similarity_matrix([], method="unknown")


def test_matrix_size_limit_consumes_only_one_item_beyond_limit() -> None:
    consumed = 0

    def unbounded_sequences() -> Iterator[DNASequence]:
        nonlocal consumed
        while True:
            consumed += 1
            yield DNASequence("A")

    with pytest.raises(ConfigurationError) as size_error:
        similarity_matrix(unbounded_sequences(), method="exact", max_items=3)

    assert size_error.value.code == "SIMILARITY_MATRIX_SIZE_LIMIT"
    assert size_error.value.context == {
        "item_count": 4,
        "item_count_is_lower_bound": True,
        "max_items": 3,
    }
    assert consumed == 4
