from __future__ import annotations

from collections.abc import Iterator

import pytest

from dnakit import DNASequence, Gap
from dnakit.core import CompoundLocation, Interval, Strand, UnresolvedLocation
from dnakit.exceptions import ConfigurationError, SequenceError
from dnakit.patterns import (
    PAMRule,
    RestrictionEnzyme,
    scan_codon_sites,
    scan_orfs,
    scan_pam_candidates,
    scan_restriction_sites,
)


def test_codon_sites_and_complete_orf_use_six_frame_coordinates() -> None:
    sequence = DNASequence("ATGAAATAA")

    sites = scan_codon_sites(sequence, strand="forward")
    orfs = scan_orfs(sequence, strand="forward")

    observed_sites = {(hit.kind, hit.codon, hit.frame) for hit in sites.hits}
    assert ("start", "ATG", 1) in observed_sites
    assert ("stop", "TAA", 1) in observed_sites
    assert len(orfs.hits) == 1
    assert orfs.hits[0].symbol_location == Interval(0, 9)
    assert orfs.hits[0].translation == "MK*"
    assert orfs.hits[0].complete


def test_reverse_orf_nested_starts_and_incomplete_policy() -> None:
    reverse = scan_orfs(DNASequence("TTATTTCAT"), strand="reverse")
    nested = scan_orfs(DNASequence("ATGATGTAA"), strand="forward")
    incomplete = scan_orfs(
        DNASequence("CCCATGAAA"),
        strand="forward",
        require_complete=False,
    )

    assert reverse.hits[0].strand is Strand.REVERSE
    assert reverse.hits[0].symbol_location == Interval(0, 9)
    assert [hit.start_codon for hit in nested.hits] == ["ATG", "ATG"]
    assert [hit.nucleotide_length for hit in nested.hits] == [9, 6]
    assert incomplete.hits[0].stop_codon is None
    assert incomplete.hits[0].translation == "MK"


def test_orfs_split_at_gaps_and_reject_circular_origin_ambiguity() -> None:
    gapped = DNASequence(["ATG", Gap(3), "TAA"])

    assert scan_orfs(gapped, strand="forward").hits == ()
    with pytest.raises(SequenceError, match="CIRCULAR_PATTERN_UNSUPPORTED"):
        scan_orfs(DNASequence("ATGTAA", topology="circular"))


def test_orf_after_unknown_gap_retains_symbol_location_and_unresolved_coordinate() -> None:
    sequence = DNASequence(["AAA", Gap(None), "ATGAAATAA"])

    result = scan_orfs(sequence, strand="forward")

    assert result.hits[0].symbol_location == Interval(3, 12)
    assert isinstance(result.hits[0].coordinate_location, UnresolvedLocation)


def test_code_11_alternative_start_is_supported_and_audited() -> None:
    result = scan_orfs(
        DNASequence("GTGAAATAA"),
        genetic_code=11,
        strand="forward",
    )

    assert result.hits[0].start_codon == "GTG"
    assert result.hits[0].translation == "MK*"
    assert result.parameters["genetic_code"] == 11


def test_restriction_builtin_and_custom_definitions_report_cuts() -> None:
    sequence = DNASequence("AGAATTCCCGGG")
    custom = RestrictionEnzyme("Custom", "CCCGGG", 3, 3)

    result = scan_restriction_sites(sequence, ["EcoRI", custom])

    assert [(hit.enzyme, hit.top_cut, hit.bottom_cut) for hit in result.hits] == [
        ("EcoRI", 2, 6),
        ("Custom", 9, 9),
    ]
    assert result.parameters["complete_rebase_catalog"] is False
    assert result.hits[0].strand is Strand.BOTH


def test_restriction_unknown_gap_and_small_catalog_boundary() -> None:
    sequence = DNASequence(["AAA", Gap(None), "GAATTC"])

    result = scan_restriction_sites(sequence, ["EcoRI"])

    assert isinstance(result.hits[0].coordinate_location, UnresolvedLocation)
    assert result.hits[0].top_cut is None
    with pytest.raises(ConfigurationError, match="small built-in catalog"):
        scan_restriction_sites(DNASequence("AAAA"), ["ImaginaryI"])
    with pytest.raises(ConfigurationError, match="Type IIS"):
        RestrictionEnzyme("Outside", "AAAA", 5, 1)


def test_non_palindromic_user_restriction_definition_keeps_orientation() -> None:
    enzyme = RestrictionEnzyme("Directional", "AAG", 1, 2)

    result = scan_restriction_sites(DNASequence("AAGCTT"), [enzyme])

    assert [(hit.strand, hit.top_cut, hit.bottom_cut) for hit in result.hits] == [
        (Strand.FORWARD, 1, 2),
        (Strand.REVERSE, 4, 5),
    ]


def test_circular_restriction_site_wraps_and_cut_coordinates_are_modular() -> None:
    sequence = DNASequence("AATTCCCG", topology="circular")

    result = scan_restriction_sites(sequence, ["EcoRI"])

    assert len(result.hits) == 1
    assert isinstance(result.hits[0].symbol_location, CompoundLocation)
    assert result.hits[0].wraps_origin
    assert result.hits[0].top_cut == 0
    assert result.hits[0].bottom_cut == 4


def test_spcas9_candidate_has_guide_pam_gc_and_no_prediction_claim() -> None:
    sequence = DNASequence("A" * 20 + "TGG")

    result = scan_pam_candidates(sequence, "SpCas9", strand="forward")

    assert len(result.hits) == 1
    hit = result.hits[0]
    assert hit.guide_sequence == "A" * 20
    assert hit.pam_sequence == "TGG"
    assert hit.guide_symbol_location == Interval(0, 20)
    assert hit.pam_symbol_location == Interval(20, 23)
    assert hit.gc_fraction == 0.0
    assert result.parameters["efficiency_prediction"] is False
    assert result.parameters["off_target_prediction"] is False


def test_custom_five_prime_pam_reverse_and_filters() -> None:
    rule = PAMRule("Tiny", "TT", "5prime", 4)
    forward = scan_pam_candidates(
        DNASequence("TTACGT"),
        rule,
        strand="forward",
        min_gc=0.5,
        exclude_motifs=("AAAA",),
    )
    reverse = scan_pam_candidates(
        DNASequence("ACGTAA"),
        rule,
        strand="reverse",
    )

    assert forward.hits[0].guide_sequence == "ACGT"
    assert reverse.hits[0].strand is Strand.REVERSE


def test_circular_pam_candidate_wraps_and_ambiguous_guides_are_configurable() -> None:
    rule = PAMRule("Tiny", "GG", "3prime", 4)
    circular = scan_pam_candidates(
        DNASequence("AGGAAA", topology="circular"),
        rule,
        strand="forward",
    )
    ambiguous = scan_pam_candidates(
        DNASequence("NNNNGG", alphabet="iupac"),
        rule,
        strand="forward",
        allow_ambiguous_guides=True,
    )

    assert circular.hits[0].wraps_origin
    assert isinstance(circular.hits[0].guide_symbol_location, CompoundLocation)
    assert ambiguous.hits[0].guide_sequence == "NNNN"


def test_pam_scan_does_not_cross_gap_and_marks_unknown_downstream_coordinates() -> None:
    rule = PAMRule("Tiny", "GG", "3prime", 4)
    crossing = scan_pam_candidates(
        DNASequence(["AAAA", Gap(2), "GG"]),
        rule,
        strand="forward",
    )
    downstream = scan_pam_candidates(
        DNASequence(["AAA", Gap(None), "AAAAGG"]),
        rule,
        strand="forward",
    )

    assert crossing.hits == ()
    assert isinstance(downstream.hits[0].guide_coordinate_location, UnresolvedLocation)
    assert isinstance(downstream.hits[0].pam_coordinate_location, UnresolvedLocation)


def test_coding_and_crispr_resource_limits_fail_before_unbounded_work() -> None:
    with pytest.raises(ConfigurationError, match="max_codon_checks"):
        scan_codon_sites(DNASequence("A" * 100), max_codon_checks=1)
    with pytest.raises(ConfigurationError, match="max_scan_length"):
        scan_pam_candidates(DNASequence("A" * 30), "SpCas9", max_scan_length=10)

    consumed = 0

    def codons() -> Iterator[str]:
        nonlocal consumed
        while True:
            consumed += 1
            yield "ATG"

    with pytest.raises(ConfigurationError) as codon_limit:
        scan_codon_sites(DNASequence("ATG"), start_codons=codons())
    assert codon_limit.value.code == "CODON_SET_SIZE_LIMIT"
    assert consumed == 65

    motif_consumed = 0

    def motifs() -> Iterator[str]:
        nonlocal motif_consumed
        while True:
            motif_consumed += 1
            yield "A"

    with pytest.raises(ConfigurationError, match="max_exclude_motifs"):
        scan_pam_candidates(
            DNASequence("A" * 20 + "TGG"),
            "SpCas9",
            exclude_motifs=motifs(),
            max_exclude_motifs=1,
        )
    assert motif_consumed == 2
