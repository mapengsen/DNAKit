"""Human-verifiable exact and reverse-complement search tests."""

from collections.abc import Callable, Iterator

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, DNASet, Gap, Strand
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.similarity import exact_search, reverse_complement_search, subsequence_search


def test_exact_search_distinguishes_whole_sequence_from_subsequence() -> None:
    query = DNARecord(DNASequence("AC"), "query")
    targets = DNASet.from_records(
        [
            DNARecord(DNASequence("AC"), "equal"),
            DNARecord(DNASequence("TACG"), "contains"),
            DNARecord(DNASequence("AG"), "different"),
        ]
    )

    result = exact_search(query, targets)

    assert result.found
    assert result.query_id == "query"
    assert result.target_count == 3
    assert [(hit.target_id, hit.start, hit.end) for hit in result.matches] == [("equal", 0, 2)]
    assert result.full_length
    assert result.coordinate_system == "0-based-half-open"


def test_exact_search_marks_forward_reverse_and_palindromic_hits() -> None:
    query = DNASequence("AN", alphabet=DNAAlphabet.IUPAC)
    targets = DNASet.from_records(
        [
            DNARecord(DNASequence("AN", alphabet=DNAAlphabet.IUPAC), "forward"),
            DNARecord(DNASequence("NT", alphabet=DNAAlphabet.IUPAC), "reverse"),
        ]
    )

    result = exact_search(query, targets, reverse_complement=True)
    assert [(hit.target_id, hit.strand) for hit in result.matches] == [
        ("forward", Strand.FORWARD),
        ("reverse", Strand.REVERSE),
    ]

    palindrome = exact_search(
        DNASequence("AT"),
        DNASequence("AT"),
        reverse_complement=True,
    )
    assert palindrome.matches[0].strand is Strand.BOTH


def test_subsequence_search_overlap_and_literal_iupac_behavior() -> None:
    query = DNASequence("ANA", alphabet=DNAAlphabet.IUPAC)
    target = DNASequence("ANANA", alphabet=DNAAlphabet.IUPAC)

    overlapping = subsequence_search(query, target)
    non_overlapping = subsequence_search(query, target, overlapping=False)

    assert [(hit.start, hit.end) for hit in overlapping.matches] == [(0, 3), (2, 5)]
    assert [(hit.start, hit.end) for hit in non_overlapping.matches] == [(0, 3)]
    assert not subsequence_search(
        DNASequence("N", alphabet=DNAAlphabet.IUPAC),
        DNASequence("A"),
    ).found
    assert overlapping.iupac_matching == "literal"

    across_targets = subsequence_search(
        DNASequence("AC"),
        [
            DNARecord(DNASequence("TAC"), "first"),
            DNARecord(DNASequence("ACAC"), "second"),
        ],
    )
    assert across_targets.target_count == 2
    assert [(hit.target_index, hit.target_id, hit.start) for hit in across_targets.matches] == [
        (0, "first", 1),
        (1, "second", 0),
        (1, "second", 2),
    ]


def test_reverse_complement_search_and_palindrome_merge_are_explicit() -> None:
    reverse_only = reverse_complement_search(DNASequence("ATG"), DNASequence("GGCATCC"))
    assert [(hit.start, hit.end, hit.strand) for hit in reverse_only.matches] == [
        (2, 5, Strand.REVERSE)
    ]

    merged = reverse_complement_search(DNASequence("AT"), DNASequence("ATAT"))
    unmerged = reverse_complement_search(
        DNASequence("AT"),
        DNASequence("ATAT"),
        merge_strands=False,
    )
    assert [(hit.start, hit.strand) for hit in merged.matches] == [
        (0, Strand.BOTH),
        (2, Strand.BOTH),
    ]
    assert len(unmerged.matches) == 4


def test_empty_query_matches_each_subsequence_boundary_but_only_empty_exact_target() -> None:
    subsequence = subsequence_search(DNASequence(""), DNASequence("AC"))
    exact = exact_search(DNASequence(""), [DNASequence(""), DNASequence("A")])

    assert [(hit.start, hit.end) for hit in subsequence.matches] == [(0, 0), (1, 1), (2, 2)]
    assert [hit.target_index for hit in exact.matches] == [0]


def test_exact_search_bounds_infinite_target_iterables_without_silent_truncation() -> None:
    consumed = 0

    def infinite_targets() -> Iterator[DNASequence]:
        nonlocal consumed
        while True:
            consumed += 1
            yield DNASequence("A")

    with pytest.raises(ConfigurationError) as error:
        exact_search(DNASequence("A"), infinite_targets(), max_targets=3)

    assert error.value.code == "SEARCH_TARGET_LIMIT_EXCEEDED"
    assert error.value.context == {
        "target_count": 4,
        "target_count_is_lower_bound": True,
        "max_targets": 3,
    }
    assert consumed == 4


@pytest.mark.parametrize("search", [subsequence_search, reverse_complement_search])
def test_subsequence_searches_bound_infinite_target_iterables(
    search: Callable[..., object],
) -> None:
    consumed = 0

    def infinite_targets() -> Iterator[DNASequence]:
        nonlocal consumed
        while True:
            consumed += 1
            yield DNASequence("A")

    with pytest.raises(ConfigurationError) as error:
        search(DNASequence("A"), infinite_targets(), max_targets=2)

    assert error.value.code == "SEARCH_TARGET_LIMIT_EXCEEDED"
    assert consumed == 3


def test_search_match_limit_rejects_empty_query_and_reverse_hits() -> None:
    with pytest.raises(ConfigurationError) as empty_error:
        subsequence_search(DNASequence(""), DNASequence("ACGT"), max_matches=4)
    assert empty_error.value.code == "SEARCH_MATCH_LIMIT_EXCEEDED"
    assert empty_error.value.context == {
        "match_count": 5,
        "match_count_is_lower_bound": True,
        "max_matches": 4,
    }

    with pytest.raises(ConfigurationError) as reverse_error:
        reverse_complement_search(
            DNASequence("AT"),
            DNASequence("ATAT"),
            merge_strands=False,
            max_matches=3,
        )
    assert reverse_error.value.code == "SEARCH_MATCH_LIMIT_EXCEEDED"


def test_search_limits_are_auditable_and_palindrome_merging_counts_returned_hits() -> None:
    result = reverse_complement_search(
        DNASequence("AT"),
        DNASequence("ATAT"),
        max_targets=7,
        max_matches=2,
    )

    assert len(result.matches) == 2
    assert result.max_targets == 7
    assert result.max_matches == 2
    assert result.to_dict()["max_targets"] == 7
    assert result.to_dict()["max_matches"] == 2


def test_search_limits_must_be_positive_integers() -> None:
    with pytest.raises(ConfigurationError):
        exact_search(DNASequence("A"), [], max_targets=0)
    with pytest.raises(ConfigurationError):
        subsequence_search(DNASequence("A"), DNASequence("A"), max_matches=True)


def test_search_rejects_gap_and_invalid_runtime_policies() -> None:
    gapped = DNASequence(["A", Gap(2), "C"])

    with pytest.raises(UnsupportedGapOperationError) as error:
        exact_search(DNASequence("AC"), gapped)
    assert error.value.code == "SIMILARITY_GAP_NOT_ALLOWED"
    with pytest.raises(UnsupportedGapOperationError):
        subsequence_search(gapped, DNASequence("AC"))
    with pytest.raises(ConfigurationError):
        subsequence_search(DNASequence("A"), DNASequence("A"), strand="unknown")
    with pytest.raises(ConfigurationError):
        exact_search(DNASequence("A"), "A")  # type: ignore[arg-type]


def test_search_rejects_first_invalid_generated_target_without_overconsuming() -> None:
    consumed = 0

    def targets() -> Iterator[object]:
        nonlocal consumed
        consumed += 1
        yield "invalid"
        while True:
            consumed += 1
            yield DNASequence("A")

    with pytest.raises(ConfigurationError) as error:
        exact_search(DNASequence("A"), targets(), max_targets=10)  # type: ignore[arg-type]

    assert error.value.code == "INVALID_SEARCH_TARGET"
    assert error.value.context["target_index"] == 0
    assert consumed == 1
