"""Tests for bounded approximate substring matching."""

from collections.abc import Iterator

import pytest

from dnakit.core import DNARecord, DNASequence, Gap, Strand, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.similarity import approximate_search


def test_approximate_search_reports_exact_and_mismatch_hits() -> None:
    result = approximate_search(
        DNARecord(DNASequence("ACG"), "query"),
        DNARecord(DNASequence("TTACGATG"), "target"),
        max_distance=1,
    )

    assert result.query_id == "query"
    assert result.target_count == 1
    assert any(
        (match.start, match.end, match.distance, match.target_id) == (2, 5, 0, "target")
        for match in result.matches
    )
    assert any((match.start, match.end, match.distance) == (5, 8, 1) for match in result.matches)
    assert result.dp_cells == 4 * 9


def test_approximate_search_supports_indels_and_reverse_complement() -> None:
    indel = approximate_search(DNASequence("ACG"), DNASequence("TACCGT"), max_distance=1)
    reverse = approximate_search(
        DNASequence("ATG"),
        DNASequence("GGCATCC"),
        max_distance=0,
        reverse_complement=True,
    )

    assert any((hit.start, hit.end, hit.distance) == (1, 4, 1) for hit in indel.matches)
    assert [(hit.start, hit.end, hit.strand) for hit in reverse.matches] == [(2, 5, Strand.REVERSE)]


def test_approximate_search_bounds_targets_cells_and_matches() -> None:
    consumed = 0

    def targets() -> Iterator[DNASequence]:
        nonlocal consumed
        while True:
            consumed += 1
            yield DNASequence("AAA")

    with pytest.raises(ConfigurationError) as targets_error:
        approximate_search(DNASequence("A"), targets(), max_distance=0, max_targets=2)
    assert targets_error.value.code == "SEARCH_TARGET_LIMIT_EXCEEDED"
    assert consumed == 3
    with pytest.raises(ConfigurationError) as cells_error:
        approximate_search(DNASequence("AAA"), DNASequence("AAA"), max_distance=0, max_cells=15)
    assert cells_error.value.code == "APPROXIMATE_MATCH_CELL_LIMIT"
    with pytest.raises(ConfigurationError) as matches_error:
        approximate_search(DNASequence("A"), DNASequence("AAA"), max_distance=0, max_matches=2)
    assert matches_error.value.code == "APPROXIMATE_MATCH_LIMIT_EXCEEDED"


def test_approximate_search_rejects_ambiguous_semantics() -> None:
    with pytest.raises(ConfigurationError) as empty:
        approximate_search(DNASequence(""), DNASequence("A"), max_distance=0)
    assert empty.value.code == "EMPTY_APPROXIMATE_QUERY"
    with pytest.raises(ConfigurationError):
        approximate_search(
            DNASequence("A", topology=Topology.CIRCULAR), DNASequence("A"), max_distance=0
        )
    with pytest.raises(ConfigurationError):
        approximate_search(DNASequence("A"), DNASequence("A"), max_distance=float("nan"))


@pytest.mark.parametrize(
    "query,target,role",
    [
        (DNASequence(["A", Gap(2), "C"]), DNASequence("AC"), "query"),
        (DNASequence("AC"), DNASequence(["A", Gap(2), "C"]), "target"),
    ],
)
def test_approximate_search_never_silently_drops_gaps(
    query: DNASequence, target: DNASequence, role: str
) -> None:
    with pytest.raises(UnsupportedGapOperationError) as error:
        approximate_search(query, target, max_distance=0)
    assert error.value.code == "APPROXIMATE_GAP_NOT_ALLOWED"
    assert error.value.context["role"] == role
