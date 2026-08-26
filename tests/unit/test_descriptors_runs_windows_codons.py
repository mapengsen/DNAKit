"""Human-verifiable tests for DESC-009, DESC-011, and DESC-012."""

from __future__ import annotations

import math
from collections.abc import Iterator

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Gap
from dnakit.descriptors import codon_statistics, homopolymer_runs, window_descriptors
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError


def test_homopolymer_coordinates_and_longest_by_base() -> None:
    result = homopolymer_runs(DNASequence("AAACCGTTTT"), min_run_length=2)

    assert result.longest_length == 4
    assert result.longest_by_base == {"A": 3, "C": 2, "G": 1, "T": 4}
    assert [(run.base, run.length, run.symbol_start, run.symbol_end) for run in result.runs] == [
        ("A", 3, 0, 3),
        ("C", 2, 3, 5),
        ("T", 4, 6, 10),
    ]
    assert result.runs[-1].coordinate_start == 6
    assert result.runs[-1].coordinate_end == 10


def test_homopolymer_gap_and_iupac_boundaries_are_not_hidden() -> None:
    gapped = DNASequence(["AA", Gap(3), "AA"])
    unknown = DNASequence(["AA", Gap(None), "AA"])
    ambiguous = DNASequence("AANAA", alphabet=DNAAlphabet.IUPAC)

    assert [run.length for run in homopolymer_runs(gapped).runs] == [2, 2]
    crossed = homopolymer_runs(gapped, cross_gaps=True)
    assert [run.length for run in crossed.runs] == [4]
    assert crossed.runs[0].coordinate_start == 0
    assert crossed.runs[0].coordinate_end == 7
    assert crossed.runs[0].crossed_gap_count == 1
    crossed_unknown = homopolymer_runs(unknown, cross_gaps=True).runs[0]
    assert crossed_unknown.coordinate_end is None
    assert crossed_unknown.crossed_unknown_gap
    ignored = homopolymer_runs(ambiguous, ambiguity_policy="ignore")
    assert [run.length for run in ignored.runs] == [2, 2]
    with pytest.raises(InvalidAlphabetError):
        homopolymer_runs(ambiguous)


def test_window_gc_entropy_and_cpg_values_are_hand_checkable() -> None:
    result = window_descriptors(
        DNASequence("ACGT"),
        ["gc", "entropy", "cpg"],
        window_size=2,
        step=2,
    )

    assert [(row.symbol_start, row.symbol_end) for row in result.windows] == [(0, 2), (2, 4)]
    first = result.windows[0]
    assert first.values["gc_fraction"] == 0.5
    assert first.values["shannon_entropy"] == 1.0
    assert first.values["cpg_count"] == 0
    assert first.coordinate_start == 0
    assert first.coordinate_end == 2


def test_window_result_records_entropy_log_base() -> None:
    result = window_descriptors(
        DNASequence("ACGT"),
        ["entropy"],
        window_size=4,
        entropy_log_base=10,
    )

    assert result.entropy_log_base == 10.0
    assert result.to_dict()["entropy_log_base"] == 10.0
    assert result.windows[0].values["shannon_entropy"] == pytest.approx(math.log10(4))


def test_windows_respect_gap_series_and_partial_boundary_policy() -> None:
    sequence = DNASequence(["ACG", Gap(5), "TTA"])
    separated = window_descriptors(sequence, ["gc"], window_size=2, step=2)
    partial = window_descriptors(
        sequence,
        ["gc"],
        window_size=2,
        step=2,
        include_partial=True,
    )
    crossed = window_descriptors(sequence, ["gc"], window_size=4, cross_gaps=True)

    assert [(row.symbol_start, row.symbol_end) for row in separated.windows] == [(0, 2), (3, 5)]
    assert [(row.symbol_start, row.symbol_end) for row in partial.windows] == [
        (0, 2),
        (2, 3),
        (3, 5),
        (5, 6),
    ]
    assert [row.is_partial for row in partial.windows] == [False, True, False, True]
    assert crossed.windows[0].crossed_gap_count == 1
    assert crossed.windows[0].coordinate_end == 9


def test_window_ignore_policy_keeps_ambiguity_in_position_but_not_denominator() -> None:
    sequence = DNASequence("ANG", alphabet=DNAAlphabet.IUPAC)
    result = window_descriptors(
        sequence,
        ["gc", "entropy", "cpg"],
        window_size=3,
        ambiguity_policy="ignore",
    )
    values = result.windows[0].values

    assert values["canonical_base_denominator"] == 2
    assert values["gc_fraction"] == 0.5
    assert values["shannon_entropy"] == 1.0
    assert values["cpg_count"] == 0
    assert values["cpg_pair_denominator"] == 0


def test_codon_counts_frequencies_and_start_stop_density() -> None:
    result = codon_statistics(DNASequence("ATGAAATAA"))

    assert result.counts == {"AAA": 1, "ATG": 1, "TAA": 1}
    assert result.codon_count == 3
    assert result.start_count == 1
    assert result.stop_count == 1
    assert math.isclose(result.start_density or 0.0, 1 / 3)
    assert math.isclose(result.stop_density or 0.0, 1 / 3)
    assert result.incomplete_base_count == 0


def test_codon_frame_ambiguity_gap_and_trailing_bases_are_explicit() -> None:
    framed = codon_statistics(DNASequence("AATGAAATAAC"), frame=1)
    ambiguous = codon_statistics(
        DNASequence("ATGNNTTAA", alphabet=DNAAlphabet.IUPAC),
        ambiguity_policy="ignore",
    )
    gapped_sequence = DNASequence(["AT", Gap(4), "GAAA"])
    gapped = codon_statistics(gapped_sequence)
    crossed = codon_statistics(gapped_sequence, cross_gaps=True)

    assert framed.counts == {"AAA": 1, "ATG": 1, "TAA": 1}
    assert framed.incomplete_base_count == 1
    assert ambiguous.counts == {"ATG": 1, "TAA": 1}
    assert ambiguous.ignored_ambiguity_codon_count == 1
    assert gapped.counts == {"GAA": 1}
    assert gapped.gap_interrupted_codon_count == 1
    assert gapped.incomplete_base_count == 1
    assert gapped.phase_coordinate_system == "sequence-coordinate"
    assert crossed.counts == {"AAA": 1, "ATG": 1}
    assert crossed.gap_interrupted_codon_count == 0
    assert crossed.phase_coordinate_system == "concatenated-symbol"


def test_symbol_algorithms_never_cross_a_non_crossable_gap() -> None:
    homopolymer = homopolymer_runs(
        DNASequence(["AA", Gap(1, crossable=False), "AA"]),
        cross_gaps=True,
    )
    windows = window_descriptors(
        DNASequence(["AC", Gap(1, crossable=False), "GT"]),
        ["gc"],
        window_size=3,
        cross_gaps=True,
    )
    codons = codon_statistics(
        DNASequence(["AT", Gap(1, crossable=False), "G"]),
        cross_gaps=True,
    )

    assert [run.length for run in homopolymer.runs] == [2, 2]
    assert windows.windows == ()
    assert codons.codon_count == 0
    assert codons.gap_interrupted_codon_count == 1


def test_codon_phase_becomes_unresolved_after_unknown_length_gap() -> None:
    sequence = DNASequence(["AT", Gap(None), "GAAA"])

    result = codon_statistics(sequence)

    assert result.counts == {}
    assert result.gap_interrupted_codon_count == 1
    assert result.phase_unresolved_after_gap
    assert result.unresolved_downstream_base_count == 4


def test_known_gap_length_advances_codon_phase() -> None:
    result = codon_statistics(DNASequence(["A", Gap(1), "TGAAA"]))

    assert result.counts == {"GAA": 1}
    assert result.gap_interrupted_codon_count == 1
    assert result.incomplete_base_count == 1


@pytest.mark.parametrize(
    "kwargs",
    [
        {"frame": 3},
        {"genetic_code": 2},
    ],
)
def test_codon_rejects_unsupported_configuration(kwargs: dict[str, int]) -> None:
    with pytest.raises(ConfigurationError):
        codon_statistics(DNASequence("ATG"), **kwargs)  # type: ignore[arg-type]


def test_window_rejects_unknown_or_duplicate_descriptor() -> None:
    with pytest.raises(ConfigurationError):
        window_descriptors(DNASequence("ACGT"), ["repeat"], window_size=2)
    with pytest.raises(ConfigurationError):
        window_descriptors(DNASequence("ACGT"), ["gc", "gc"], window_size=2)


def test_window_descriptor_iterable_is_bounded() -> None:
    consumed = 0

    def names() -> Iterator[str]:
        nonlocal consumed
        for name in ("gc", "entropy", "cpg", "gc"):
            consumed += 1
            yield name

    with pytest.raises(ConfigurationError):
        window_descriptors(DNASequence("ACGT"), names(), window_size=2)
    assert consumed == 4
