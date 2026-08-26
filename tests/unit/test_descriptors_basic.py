"""Human-verifiable tests for DESC-001..005."""

from __future__ import annotations

import math
from dataclasses import FrozenInstanceError

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.descriptors import (
    base_composition,
    base_skew,
    cpg_features,
    gc_at_content,
    length_features,
)
from dnakit.exceptions import InvalidAlphabetError


def test_length_features_keep_symbol_gap_and_unknown_span_distinct() -> None:
    known = length_features(DNASequence(["ACN", Gap(4), "GT"], alphabet=DNAAlphabet.IUPAC))
    unknown = length_features(DNASequence(["AC", Gap(None), "GT"]))

    assert known.symbol_length == 5
    assert known.coordinate_span == 9
    assert known.canonical_base_count == 4
    assert known.ambiguity_length == 1
    assert known.known_gap_length == 4
    assert known.gap_count == 1
    assert known.unknown_gap_count == 0
    assert unknown.coordinate_span is None
    assert unknown.known_gap_length == 0
    assert unknown.unknown_gap_count == 1


def test_record_identifier_and_composition_denominator_are_auditable() -> None:
    record = DNARecord(
        DNASequence("AACGTN", alphabet=DNAAlphabet.IUPAC),
        "example",
    )
    result = base_composition(record, ambiguity_policy="ignore")

    assert result.sequence_id == "example"
    assert result.counts == {"A": 2, "C": 1, "G": 1, "T": 1}
    assert result.denominator == 5
    assert result.fractions["A"] == 0.4
    assert result.ignored_ambiguity_count == 1
    with pytest.raises(InvalidAlphabetError) as error:
        base_composition(record)
    assert error.value.code == "DESCRIPTOR_AMBIGUITY_NOT_ALLOWED"


def test_gc_at_and_skew_use_declared_formulas_and_zero_is_undefined() -> None:
    content = gc_at_content(DNASequence("AACG"))
    skew = base_skew(DNASequence("AACG"))
    empty_skew = base_skew(DNASequence(""))

    assert content.gc_count == 2
    assert content.at_count == 2
    assert content.gc_fraction == 0.5
    assert content.at_fraction == 0.5
    assert skew.gc_skew == 0.0
    assert skew.at_skew == 1.0
    assert empty_skew.gc_skew is None
    assert empty_skew.at_skew is None


def test_cpg_count_density_and_oe_have_explicit_denominators() -> None:
    result = cpg_features(DNASequence("ACGCGT"))

    assert result.cpg_count == 2
    assert result.adjacent_pair_denominator == 5
    assert math.isclose(result.density or 0.0, 0.4)
    assert result.expected_length_denominator == 6
    assert result.observed_expected == 3.0
    assert result.density_formula == "count(CG)/eligible_adjacent_canonical_pairs"


def test_cpg_does_not_bridge_gap_or_ignored_iupac_unless_requested() -> None:
    gapped = DNASequence(["C", Gap(None), "G"])
    ambiguous = DNASequence("CNG", alphabet=DNAAlphabet.IUPAC)

    default = cpg_features(gapped)
    crossed = cpg_features(gapped, cross_gaps=True)
    ignored = cpg_features(ambiguous, ambiguity_policy="ignore")

    assert default.cpg_count == 0
    assert default.density is None
    assert crossed.cpg_count == 1
    assert crossed.density == 1.0
    assert crossed.unknown_gap_count == 1
    assert ignored.cpg_count == 0
    assert ignored.adjacent_pair_denominator == 0


def test_cpg_never_bridges_a_non_crossable_gap() -> None:
    sequence = DNASequence(["C", Gap(1, crossable=False), "G"])

    result = cpg_features(sequence, cross_gaps=True)

    assert result.cpg_count == 0
    assert result.adjacent_pair_denominator == 0


def test_descriptor_results_and_nested_mappings_are_immutable() -> None:
    result = base_composition(DNASequence("ACGT"))

    with pytest.raises(FrozenInstanceError):
        result.denominator = 0  # type: ignore[misc]
    with pytest.raises(TypeError):
        result.counts["A"] = 99  # type: ignore[index]
    assert result.to_dict()["counts"] == {"A": 1, "C": 1, "G": 1, "T": 1}
