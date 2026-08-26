"""Literal exact sequence and subsequence search."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from dnakit.core import DNASequence, DNASet, Strand
from dnakit.exceptions import ConfigurationError
from dnakit.similarity._shared import (
    SequenceInput,
    materialize_targets,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from dnakit.similarity.results import SearchResult, SequenceMatch

DEFAULT_MAX_SEARCH_TARGETS = 100_000
DEFAULT_MAX_SEARCH_MATCHES = 1_000_000


def _coerce_strand(value: Strand | str) -> Strand:
    try:
        resolved = value if isinstance(value, Strand) else Strand(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Unknown search strand.",
            code="INVALID_SEARCH_STRAND",
            context={"strand": value},
            hint="Choose 'forward', 'reverse', or 'both'.",
        ) from exc
    if resolved is Strand.UNKNOWN:
        raise ConfigurationError(
            "Search strand cannot be 'unknown'.",
            code="INVALID_SEARCH_STRAND",
        )
    return resolved


def _literal_occurrences(query: str, target: str, *, overlapping: bool) -> Iterator[int]:
    if query == "":
        yield from range(len(target) + 1)
        return
    cursor = 0
    step = 1 if overlapping else len(query)
    while cursor <= len(target) - len(query):
        position = target.find(query, cursor)
        if position < 0:
            break
        yield position
        cursor = position + step


def _orientations(query: DNASequence, strand: Strand) -> tuple[tuple[Strand, str], ...]:
    selected: list[tuple[Strand, str]] = []
    if strand in (Strand.FORWARD, Strand.BOTH):
        selected.append((Strand.FORWARD, query.symbols))
    if strand in (Strand.REVERSE, Strand.BOTH):
        selected.append((Strand.REVERSE, query.reverse_complement().symbols))
    return tuple(selected)


def _ordered_matches(matches: list[SequenceMatch]) -> tuple[SequenceMatch, ...]:
    strand_order = {Strand.FORWARD: 0, Strand.REVERSE: 1, Strand.BOTH: 2}
    return tuple(
        sorted(
            matches,
            key=lambda match: (
                match.target_index,
                match.start,
                match.end,
                strand_order[match.strand],
            ),
        )
    )


class _BoundedMatchCollector:
    """Collect returned matches while enforcing a strict non-truncating limit."""

    def __init__(self, *, max_matches: int, merge_strands: bool) -> None:
        self._max_matches = max_matches
        self._merge_strands = merge_strands
        self._matches: list[SequenceMatch] = []
        self._groups: dict[tuple[int, str | None, int, int], set[Strand]] = {}

    def add(self, match: SequenceMatch) -> None:
        if not self._merge_strands:
            if len(self._matches) >= self._max_matches:
                self._raise_limit()
            self._matches.append(match)
            return
        key = (match.target_index, match.target_id, match.start, match.end)
        strands = self._groups.get(key)
        if strands is None:
            if len(self._groups) >= self._max_matches:
                self._raise_limit()
            self._groups[key] = {match.strand}
        else:
            strands.add(match.strand)

    def materialize(self) -> tuple[SequenceMatch, ...]:
        if not self._merge_strands:
            return _ordered_matches(self._matches)
        matches = [
            SequenceMatch(
                target_index,
                target_id,
                start,
                end,
                Strand.BOTH if len(strands) > 1 else next(iter(strands)),
            )
            for (target_index, target_id, start, end), strands in self._groups.items()
        ]
        return _ordered_matches(matches)

    def _raise_limit(self) -> None:
        raise ConfigurationError(
            "Search result exceeds max_matches.",
            code="SEARCH_MATCH_LIMIT_EXCEEDED",
            context={
                "match_count": self._max_matches + 1,
                "match_count_is_lower_bound": True,
                "max_matches": self._max_matches,
            },
            hint="Narrow the search or increase max_matches explicitly.",
        )


def exact_search(
    query: SequenceInput,
    targets: SequenceInput | DNASet | Iterable[SequenceInput],
    *,
    reverse_complement: bool = False,
    merge_strands: bool = True,
    max_targets: int = DEFAULT_MAX_SEARCH_TARGETS,
    max_matches: int = DEFAULT_MAX_SEARCH_MATCHES,
) -> SearchResult:
    """Find targets whose complete symbol strings equal the query.

    Record metadata, topology, and strandedness do not change literal symbol
    equality.  With ``reverse_complement=True``, a reverse-only match is marked
    ``Strand.REVERSE``; a palindromic hit is merged to ``Strand.BOTH`` by
    default.  Empty query matches only an empty target.
    """

    validate_bool(reverse_complement, "reverse_complement")
    validate_bool(merge_strands, "merge_strands")
    validate_positive_int(max_targets, "max_targets")
    validate_positive_int(max_matches, "max_matches")
    query_sequence, query_id = sequence_and_id(query, role="query")
    target_items = materialize_targets(targets, max_targets=max_targets)
    strand = Strand.BOTH if reverse_complement else Strand.FORWARD
    matches = _BoundedMatchCollector(
        max_matches=max_matches,
        merge_strands=merge_strands,
    )
    for target_index, target_input in enumerate(target_items):
        target, target_id = sequence_and_id(target_input, role=f"target[{target_index}]")
        for orientation, oriented_query in _orientations(query_sequence, strand):
            if oriented_query == target.symbols:
                matches.add(
                    SequenceMatch(
                        target_index,
                        target_id,
                        0,
                        target.symbol_length,
                        orientation,
                    )
                )
    materialized = matches.materialize()
    return SearchResult(
        name="exact_search",
        method="literal_full_sequence_equality",
        query_id=query_id,
        query_length=query_sequence.symbol_length,
        target_count=len(target_items),
        matches=materialized,
        overlapping=False,
        reverse_complement=reverse_complement,
        merge_strands=merge_strands,
        full_length=True,
        max_targets=max_targets,
        max_matches=max_matches,
    )


def subsequence_search(
    query: SequenceInput,
    target: SequenceInput | DNASet | Iterable[SequenceInput],
    *,
    strand: Strand | str = Strand.FORWARD,
    overlapping: bool = True,
    merge_strands: bool = True,
    max_targets: int = DEFAULT_MAX_SEARCH_TARGETS,
    max_matches: int = DEFAULT_MAX_SEARCH_MATCHES,
) -> SearchResult:
    """Find literal query occurrences across one or more targets.

    Coordinates are internal zero-based half-open intervals and never wrap a
    circular origin. IUPAC symbols are ordinary exact symbols: ``N`` matches
    ``N``, not every canonical nucleotide. An empty query matches every target
    boundary ``0..len(target)``.
    """

    resolved_strand = _coerce_strand(strand)
    validate_bool(overlapping, "overlapping")
    validate_bool(merge_strands, "merge_strands")
    validate_positive_int(max_targets, "max_targets")
    validate_positive_int(max_matches, "max_matches")
    query_sequence, query_id = sequence_and_id(query, role="query")
    target_items = materialize_targets(target, max_targets=max_targets)
    matches = _BoundedMatchCollector(
        max_matches=max_matches,
        merge_strands=merge_strands,
    )
    for target_index, target_input in enumerate(target_items):
        target_sequence, target_id = sequence_and_id(
            target_input,
            role=f"target[{target_index}]",
        )
        for orientation, oriented_query in _orientations(query_sequence, resolved_strand):
            for start in _literal_occurrences(
                oriented_query,
                target_sequence.symbols,
                overlapping=overlapping,
            ):
                matches.add(
                    SequenceMatch(
                        target_index,
                        target_id,
                        start,
                        start + len(oriented_query),
                        orientation,
                    )
                )
    materialized = matches.materialize()
    return SearchResult(
        name="subsequence_search",
        method="literal_exact_subsequence",
        query_id=query_id,
        query_length=query_sequence.symbol_length,
        target_count=len(target_items),
        matches=materialized,
        overlapping=overlapping,
        reverse_complement=resolved_strand in (Strand.REVERSE, Strand.BOTH),
        merge_strands=merge_strands,
        full_length=False,
        max_targets=max_targets,
        max_matches=max_matches,
    )


def reverse_complement_search(
    query: SequenceInput,
    target: SequenceInput | DNASet | Iterable[SequenceInput],
    *,
    overlapping: bool = True,
    merge_strands: bool = True,
    max_targets: int = DEFAULT_MAX_SEARCH_TARGETS,
    max_matches: int = DEFAULT_MAX_SEARCH_MATCHES,
) -> SearchResult:
    """Search both the supplied query and its reverse complement."""

    return subsequence_search(
        query,
        target,
        strand=Strand.BOTH,
        overlapping=overlapping,
        merge_strands=merge_strands,
        max_targets=max_targets,
        max_matches=max_matches,
    )


find_subsequence = subsequence_search


__all__ = [
    "DEFAULT_MAX_SEARCH_MATCHES",
    "DEFAULT_MAX_SEARCH_TARGETS",
    "exact_search",
    "find_subsequence",
    "reverse_complement_search",
    "subsequence_search",
]
