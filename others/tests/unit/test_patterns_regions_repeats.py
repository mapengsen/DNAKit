from __future__ import annotations

import pytest

from dnakit import DNASequence, Gap
from dnakit.core import Interval, UnresolvedLocation
from dnakit.exceptions import ConfigurationError, SequenceError
from dnakit.patterns import (
    find_cpg_islands,
    find_inverted_repeats,
    find_low_complexity_regions,
    find_microsatellites,
    find_reverse_complement_palindromes,
    find_tandem_repeats,
)


def test_cpg_island_uses_explicit_formula_thresholds_and_merges_windows() -> None:
    result = find_cpg_islands(
        DNASequence("CG" * 20),
        window_size=10,
        step=2,
        min_gc=0.8,
        min_observed_expected=1.5,
        min_region_length=10,
    )

    assert len(result.hits) == 1
    hit = result.hits[0]
    assert hit.symbol_location == Interval(0, 40)
    assert hit.attributes["gc_fraction"] == 1.0
    assert hit.attributes["observed_expected"] == 2.0
    assert result.parameters["observed_expected_formula"] == "CpG * L / (C * G)"


def test_cpg_island_skips_ambiguous_windows_and_never_crosses_gap() -> None:
    sequence = DNASequence(["CGCG", Gap(None), "CGCGCG"], alphabet="iupac")

    result = find_cpg_islands(
        sequence,
        window_size=4,
        min_region_length=4,
        min_gc=1.0,
        min_observed_expected=1.0,
    )

    assert len(result.hits) == 2
    assert result.hits[0].symbol_location == Interval(0, 4)
    assert isinstance(result.hits[1].coordinate_location, UnresolvedLocation)


def test_low_complexity_entropy_regions_mask_while_preserving_gaps() -> None:
    sequence = DNASequence(["A" * 12, Gap(3), "ACGT" * 3], alphabet="strict")

    result = find_low_complexity_regions(
        sequence,
        window_size=6,
        min_region_length=6,
        max_entropy=0.5,
    )

    assert result.analysis.hits[0].symbol_location == Interval(0, 12)
    assert result.masked_sequence.parts[0] == "N" * 12
    assert isinstance(result.masked_sequence.parts[1], Gap)
    assert result.masked_sequence.parts[2] == "ACGT" * 3
    assert result.analysis.parameters["is_dust"] is False


def test_reverse_complement_palindromes_support_iupac_and_limits() -> None:
    result = find_reverse_complement_palindromes(
        DNASequence("GAATTC", alphabet="iupac"),
        min_length=4,
        max_length=6,
    )
    ambiguous = find_reverse_complement_palindromes(
        DNASequence("ANNT", alphabet="iupac"),
        min_length=4,
        max_length=4,
    )

    assert any(hit.sequence == "GAATTC" for hit in result.hits)
    assert ambiguous.hits[0].sequence == "ANNT"
    with pytest.raises(ConfigurationError, match="max_comparisons"):
        find_reverse_complement_palindromes(
            DNASequence("A" * 20),
            min_length=2,
            max_length=10,
            max_comparisons=1,
        )
    with pytest.raises(ConfigurationError, match="max_comparison_cells"):
        find_reverse_complement_palindromes(
            DNASequence("GAATTC"),
            min_length=4,
            max_length=6,
            max_comparison_cells=3,
        )


def test_inverted_repeat_reports_arms_loop_and_locations() -> None:
    result = find_inverted_repeats(
        DNASequence("ACGTAAACGT"),
        min_arm_length=4,
        max_arm_length=4,
        min_loop_length=2,
        max_loop_length=2,
    )

    assert len(result.hits) == 1
    hit = result.hits[0]
    assert hit.left_arm == "ACGT"
    assert hit.right_arm == "ACGT"
    assert hit.loop_length == 2
    assert hit.symbol_location == Interval(0, 10)
    assert hit.loop_symbol_location == Interval(4, 6)


def test_tandem_repeat_uses_smallest_period_and_does_not_cross_gap() -> None:
    sequence = DNASequence(["ATATAT", Gap(2), "ATAT"])

    result = find_tandem_repeats(
        sequence,
        min_unit_length=1,
        max_unit_length=3,
        min_repeats=2,
    )

    assert [(hit.unit, hit.repeat_count) for hit in result.hits] == [("AT", 3), ("AT", 2)]
    assert result.parameters["mismatches"] == 0
    assert result.parameters["indels"] == 0


def test_microsatellite_thresholds_are_unit_specific_and_audited() -> None:
    result = find_microsatellites(DNASequence("AAAAAACACACA"))

    assert [(hit.unit, hit.repeat_count) for hit in result.hits] == [("A", 6), ("CA", 3)]
    assert result.name == "microsatellite_scan"
    assert result.parameters["interruptions_allowed"] is False
    assert result.parameters["min_repeats_by_unit"] == {
        "1": 6,
        "2": 3,
        "3": 3,
        "4": 3,
        "5": 3,
        "6": 3,
    }


def test_repeat_and_region_circular_boundaries_are_explicit() -> None:
    circular = DNASequence("ATATAT", topology="circular")

    for operation in (
        find_cpg_islands,
        find_low_complexity_regions,
        find_reverse_complement_palindromes,
        find_inverted_repeats,
        find_tandem_repeats,
    ):
        with pytest.raises(SequenceError, match="CIRCULAR_PATTERN_UNSUPPORTED"):
            operation(circular)


def test_region_and_repeat_configuration_validation() -> None:
    with pytest.raises(ConfigurationError, match="max_windows"):
        find_cpg_islands(DNASequence("CG" * 20), window_size=4, max_windows=1)
    with pytest.raises(ConfigurationError, match="at least 2"):
        find_tandem_repeats(DNASequence("AAAA"), min_repeats=1)
    with pytest.raises(ConfigurationError, match="every accepted unit length"):
        find_microsatellites(DNASequence("AAAAAA"), min_repeats_by_unit={1: 6})
    with pytest.raises(ConfigurationError) as unit_limit:
        find_tandem_repeats(DNASequence("AAAA"), max_unit_length=10**9)
    assert unit_limit.value.code == "TANDEM_REPEAT_UNIT_LIMIT"
