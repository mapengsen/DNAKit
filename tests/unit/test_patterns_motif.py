from __future__ import annotations

from collections.abc import Iterator

import pytest

from dnakit import DNARecord, DNASequence, Gap
from dnakit.core import CompoundLocation, Interval, Strand, UnresolvedLocation
from dnakit.exceptions import ConfigurationError, SequenceError
from dnakit.patterns import (
    PWM,
    scan_motif,
    scan_promoter_motifs,
    scan_pwm,
    scan_tf_pwm,
)


def test_exact_iupac_and_regex_motif_scans_are_overlapping_and_auditable() -> None:
    record = DNARecord(DNASequence("AAAA", alphabet="iupac"), "seq-1")

    exact = scan_motif(record, "AA", strand="forward")
    iupac = scan_motif(DNASequence("AGAT"), "AR", mode="iupac", strand="forward")
    regex = scan_motif(DNASequence("ACGATG"), "A.G", mode="regex", strand="forward")

    assert [hit.symbol_location for hit in exact.hits] == [
        Interval(0, 2),
        Interval(1, 3),
        Interval(2, 4),
    ]
    assert exact.sequence_id == "seq-1"
    assert exact.parameters["target_iupac_rule"] == "literal"
    assert [hit.matched_sequence for hit in iupac.hits] == ["AG"]
    assert [hit.matched_sequence for hit in regex.hits] == ["ACG", "ATG"]
    assert regex.provenance.implementation.label.value == "reimplementation"


def test_reverse_strand_and_merged_palindromic_hits() -> None:
    reverse = scan_motif(DNASequence("CCAT"), "ATG", strand="reverse")
    merged = scan_motif(
        DNASequence("GAATTC"),
        "GAATTC",
        strand="both",
        merge_strands=True,
    )

    assert len(reverse.hits) == 1
    assert reverse.hits[0].strand is Strand.REVERSE
    assert reverse.hits[0].symbol_location == Interval(1, 4)
    assert len(merged.hits) == 1
    assert merged.hits[0].strand is Strand.BOTH


def test_motif_does_not_cross_gaps_and_unknown_coordinate_is_explicit() -> None:
    sequence = DNASequence(
        ["AC", Gap(3), "GT", Gap(None), "ACG"],
        alphabet="iupac",
    )

    crossing = scan_motif(sequence, "CGT", strand="forward")
    downstream = scan_motif(sequence, "ACG", strand="forward")

    assert crossing.hits == ()
    assert len(downstream.hits) == 1
    assert downstream.hits[0].symbol_location == Interval(4, 7)
    assert isinstance(downstream.hits[0].coordinate_location, UnresolvedLocation)
    assert downstream.gap_count == 2
    assert downstream.unknown_gap_count == 1


def test_fixed_motif_wraps_circular_origin_but_regex_boundary_is_rejected() -> None:
    sequence = DNASequence("AACC", topology="circular")

    result = scan_motif(sequence, "CAA", strand="forward")

    assert len(result.hits) == 1
    assert result.hits[0].wraps_origin
    assert result.hits[0].symbol_location == CompoundLocation((Interval(3, 4), Interval(0, 2)))
    with pytest.raises(SequenceError, match="CIRCULAR_REGEX_UNSUPPORTED"):
        scan_motif(sequence, "C.A", mode="regex", strand="forward")


def test_pwm_uses_explicit_log_odds_background_threshold_and_skips_ambiguity() -> None:
    pwm = PWM(
        "AC-box",
        {
            "A": [10, 0],
            "C": [0, 10],
            "G": [0, 0],
            "T": [0, 0],
        },
    )

    result = scan_pwm(
        DNASequence("ACNC", alphabet="iupac"),
        pwm,
        threshold=3.5,
        background={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        strand="forward",
    )

    assert len(result.hits) == 1
    assert result.hits[0].score == pytest.approx(4.0)
    assert result.hits[0].threshold == 3.5
    assert result.parameters["ambiguous_target_policy"] == "skip-window"
    assert result.parameters["background"] == {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25,
    }


def test_promoter_and_tf_scans_preserve_non_prediction_boundary() -> None:
    promoter = scan_promoter_motifs(DNASequence("GGTATAATCC"), strand="forward")
    pwm = PWM(
        "input",
        {"A": [4], "C": [0], "G": [0], "T": [0]},
    )
    tf = scan_tf_pwm(
        DNASequence("AA"),
        "TF-X",
        pwm,
        threshold=1.0,
        strand="forward",
    )

    assert any(hit.motif_name == "bacterial_minus_10_consensus" for hit in promoter.hits)
    assert promoter.parameters["activity_prediction"] is False
    assert tf.name == "tf_motif_scan"
    assert tf.parameters["binding_strength_prediction"] is False
    assert all(hit.motif_name == "TF-X" for hit in tf.hits)


def test_motif_limits_and_invalid_definitions_fail_cleanly() -> None:
    with pytest.raises(ConfigurationError, match="max_scan_length"):
        scan_motif(DNASequence("AAAA"), "A", max_scan_length=3)
    with pytest.raises(ConfigurationError, match="empty matches"):
        scan_motif(DNASequence("AAAA"), "A*", mode="regex")
    with pytest.raises(ConfigurationError, match="safe DNA-regex subset"):
        scan_motif(DNASequence("AAAA"), "(A+)+", mode="regex")
    with pytest.raises(ConfigurationError, match="at most one unbounded quantifier"):
        scan_motif(DNASequence("AAAA"), "A+A+", mode="regex")
    with pytest.raises(ConfigurationError, match="PWM column"):
        PWM("bad", {"A": [0], "C": [0], "G": [0], "T": [0]})


def test_pwm_bounds_row_iterables(monkeypatch: pytest.MonkeyPatch) -> None:
    import dnakit.patterns.motif as motif_module

    monkeypatch.setattr(motif_module, "MAX_PWM_LENGTH", 2)
    consumed = 0

    def values() -> Iterator[float]:
        nonlocal consumed
        while True:
            consumed += 1
            yield 1.0

    with pytest.raises(ConfigurationError) as error:
        PWM("bounded", {base: values() for base in "ACGT"})
    assert error.value.code == "PWM_LENGTH_LIMIT"
    assert consumed == 3


def test_motif_result_truncation_is_explicit() -> None:
    result = scan_motif(DNASequence("AAAA"), "A", strand="forward", max_matches=2)

    assert len(result.hits) == 2
    assert result.truncated
    assert result.max_matches == 2
    assert result.issues[0].code == "PATTERN_MATCH_LIMIT_REACHED"
