"""Human-checkable tests for bounded global and local DNA alignment."""

import json

import pytest

from dnakit.alignment import AlignmentConfig, align_pairwise
from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError


def test_global_alignment_score_identity_coverage_and_coordinates() -> None:
    result = align_pairwise(
        DNARecord(DNASequence("ACGT"), "query"),
        DNARecord(DNASequence("AGT"), "target"),
    )

    assert result.aligned_query == "ACGT"
    assert result.aligned_target == "A-GT"
    assert result.score == 2.0
    assert result.matches == 3
    assert result.deletions == 1
    assert result.identity == 0.75
    assert result.query_coverage == 1.0
    assert result.target_coverage == 1.0
    assert result.query_id == "query"
    assert result.columns[1].query_position == 1
    assert result.columns[1].target_position is None
    json.dumps(result.to_dict(), sort_keys=True)


def test_local_alignment_reports_source_span_and_partial_coverage() -> None:
    result = align_pairwise(
        DNASequence("TTACGTAA"),
        DNASequence("GGACGTCC"),
        config=AlignmentConfig(mode="local", mismatch_score=-2, gap_score=-2),
    )

    assert result.aligned_query == "ACGT"
    assert result.aligned_target == "ACGT"
    assert (result.query_start, result.query_end) == (2, 6)
    assert (result.target_start, result.target_end) == (2, 6)
    assert result.identity == 1.0
    assert result.query_coverage == 0.5
    assert result.target_coverage == 0.5


def test_empty_and_literal_iupac_alignment_are_explicit() -> None:
    empty = align_pairwise(DNASequence(""), DNASequence(""))
    literal = align_pairwise(
        DNASequence("AN", alphabet=DNAAlphabet.IUPAC),
        DNASequence("AA"),
    )

    assert empty.identity is None
    assert empty.query_coverage == empty.target_coverage == 1.0
    assert literal.matches == 1
    assert literal.mismatches == 1
    assert literal.parameters["iupac_matching"] == "literal"


def test_alignment_rejects_gaps_invalid_config_and_cell_overflow() -> None:
    with pytest.raises(UnsupportedGapOperationError):
        align_pairwise(DNASequence(["A", Gap(2), "T"]), DNASequence("AT"))
    with pytest.raises(ConfigurationError):
        align_pairwise(DNASequence("A"), DNASequence("A"), config={})  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError) as cell_error:
        align_pairwise(
            DNASequence("A" * 10),
            DNASequence("T" * 10),
            config=AlignmentConfig(max_cells=100),
        )
    assert cell_error.value.code == "ALIGNMENT_CELL_LIMIT"
    with pytest.raises(ConfigurationError) as circular_error:
        align_pairwise(
            DNASequence("AC", topology=Topology.CIRCULAR),
            DNASequence("AC"),
        )
    assert circular_error.value.code == "ALIGNMENT_CIRCULAR_UNSUPPORTED"


def test_global_tie_break_is_deterministic() -> None:
    first = align_pairwise(
        DNASequence("A"),
        DNASequence("C"),
        config=AlignmentConfig(mismatch_score=-2, gap_score=-1),
    )
    second = align_pairwise(
        DNASequence("A"),
        DNASequence("C"),
        config=AlignmentConfig(mismatch_score=-2, gap_score=-1),
    )

    assert first.to_dict() == second.to_dict()
    assert first.aligned_query == "A"
    assert first.aligned_target == "C"


def test_semi_global_alignment_has_free_terminal_gaps() -> None:
    result = align_pairwise(
        DNASequence("ACGT"),
        DNASequence("TTACGTAA"),
        config=AlignmentConfig(mode="semi_global", mismatch_score=-2, gap_score=-2),
    )

    assert result.aligned_query == "ACGT"
    assert result.aligned_target == "ACGT"
    assert (result.query_start, result.query_end) == (0, 4)
    assert (result.target_start, result.target_end) == (2, 6)
    assert result.score == 4.0
    assert result.parameters["end_gap_policy"] == "free-query-and-target-ends"


def test_affine_gap_scoring_prefers_one_long_gap() -> None:
    result = align_pairwise(
        DNASequence("AAAAAA"),
        DNASequence("AAA"),
        config=AlignmentConfig(gap_open_score=-2, gap_extend_score=-0.25),
    )

    assert result.aligned_query.replace("-", "") == "AAAAAA"
    assert result.aligned_target.replace("-", "") == "AAA"
    assert result.deletions == 3
    assert result.score == pytest.approx(0.5)
    assert result.algorithm_version == "affine-gap-dp-v1"
    assert result.parameters["dp_cells"] == 84
    assert result.provenance.implementation.label.value == "reimplementation"
    assert {citation.key for citation in result.provenance.implementation.citations} == {
        "needleman-wunsch-1970",
        "gotoh-1982",
    }


def test_local_affine_alignment_can_start_after_either_prefix() -> None:
    target_prefix = align_pairwise(
        DNASequence("A"),
        DNASequence("CA"),
        config=AlignmentConfig(mode="local", gap_open_score=-2, gap_extend_score=-0.25),
    )
    query_prefix = align_pairwise(
        DNASequence("CA"),
        DNASequence("A"),
        config=AlignmentConfig(mode="local", gap_open_score=-2, gap_extend_score=-0.25),
    )

    for result in (target_prefix, query_prefix):
        assert result.score == 1.0
        assert result.aligned_query == result.aligned_target == "A"
        assert result.parameters["end_gap_policy"] == "excluded-local-flanks"
        assert result.provenance.implementation.label.value == "reimplementation"
    assert (target_prefix.query_start, target_prefix.target_start) == (0, 1)
    assert (query_prefix.query_start, query_prefix.target_start) == (1, 0)


def test_affine_configuration_and_cell_limit_are_strict() -> None:
    with pytest.raises(ConfigurationError):
        AlignmentConfig(gap_open_score=-2)
    with pytest.raises(ConfigurationError) as error:
        align_pairwise(
            DNASequence("AAA"),
            DNASequence("AAA"),
            config=AlignmentConfig(gap_open_score=-2, gap_extend_score=-1, max_cells=47),
        )
    assert error.value.code == "ALIGNMENT_CELL_LIMIT"
