"""Concatenation of immutable linear DNA sequences."""

from __future__ import annotations

from collections.abc import Iterable

from dnakit.core import DNAAlphabet, DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.ops._common import coerce_fragment, combine_parts, promote_alphabet


def concat(
    sequences: Iterable[DNASequence | str],
    *,
    linker: DNASequence | str | None = None,
    gap: Gap | None = None,
) -> DNASequence:
    """Concatenate at least two linear sequences with an optional separator.

    ``linker`` and ``gap`` are mutually exclusive.  No rotation, implicit
    reverse-complement, or feature synchronization is performed.
    """

    raw = tuple(sequences)
    if len(raw) < 2:
        raise ConfigurationError(
            "concat requires at least two sequences.",
            code="INSUFFICIENT_CONCAT_INPUTS",
            context={"count": len(raw)},
        )
    if linker is not None and gap is not None:
        raise ConfigurationError(
            "linker and gap are mutually exclusive.",
            code="MULTIPLE_CONCAT_SEPARATORS",
        )
    if gap is not None and not isinstance(gap, Gap):
        raise ConfigurationError(
            "gap must be a Gap object.",
            code="INVALID_CONCAT_GAP",
        )

    first_raw = raw[0]
    if isinstance(first_raw, DNASequence):
        first = first_raw
    elif isinstance(first_raw, str):
        first = DNASequence(first_raw, alphabet=_infer_text_alphabet(first_raw))
    else:
        raise ConfigurationError(
            "concat inputs must be DNASequence or normalized DNA text.",
            code="INVALID_CONCAT_INPUT",
            context={"index": 0, "type": type(first_raw).__name__},
        )
    if first.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Circular input concatenation is not supported in the MVP.",
            code="CIRCULAR_CONCAT_NOT_SUPPORTED",
        )
    fragments = [
        coerce_fragment(
            item,
            fallback_alphabet=first.alphabet,
            strandedness=first.strandedness,
        )
        for item in raw
    ]
    separator: DNASequence | Gap | None
    if linker is not None:
        separator = coerce_fragment(
            linker,
            fallback_alphabet=first.alphabet,
            strandedness=first.strandedness,
        )
    else:
        separator = gap

    alphabets = [fragment.alphabet for fragment in fragments]
    if isinstance(separator, DNASequence):
        alphabets.append(separator.alphabet)
    alphabet = promote_alphabet(*alphabets)
    parts: list[str | Gap] = []
    for index, fragment in enumerate(fragments):
        if index and separator is not None:
            if isinstance(separator, DNASequence):
                parts.extend(separator.parts)
            else:
                parts.append(separator)
        parts.extend(fragment.parts)
    return combine_parts(
        parts,
        alphabet=alphabet,
        topology=Topology.LINEAR,
        strandedness=first.strandedness,
    )


def concat_overlap(
    sequences: Iterable[DNASequence | str],
    *,
    min_overlap: int = 1,
    max_overlap: int | None = None,
) -> DNASequence:
    """Concatenate two linear sequences while retaining an exact overlap once.

    The longest exact suffix/prefix overlap is selected automatically.  IUPAC
    symbols are compared literally; ambiguity compatibility is not inferred.
    """

    _validate_overlap_bounds(min_overlap, max_overlap)
    raw = tuple(sequences)
    if len(raw) != 2:
        raise ConfigurationError(
            "concat_overlap requires exactly two sequences.",
            code="INVALID_OVERLAP_CONCAT_INPUTS",
            context={"count": len(raw)},
        )

    first_raw = raw[0]
    if isinstance(first_raw, DNASequence):
        first = first_raw
    elif isinstance(first_raw, str):
        first = DNASequence(first_raw, alphabet=_infer_text_alphabet(first_raw))
    else:
        raise ConfigurationError(
            "concat_overlap inputs must be DNASequence or normalized DNA text.",
            code="INVALID_OVERLAP_CONCAT_INPUT",
            context={"index": 0, "type": type(first_raw).__name__},
        )
    if first.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Circular input overlap concatenation is not supported.",
            code="CIRCULAR_CONCAT_NOT_SUPPORTED",
        )

    fragments = tuple(
        coerce_fragment(
            item,
            fallback_alphabet=first.alphabet,
            strandedness=first.strandedness,
        )
        for item in raw
    )
    if any(fragment.is_gapped for fragment in fragments):
        raise UnsupportedGapOperationError(
            "concat_overlap cannot silently discard explicit Gap objects.",
            code="OVERLAP_CONCAT_GAPPED_INPUT",
            hint="Resolve the gaps first or use concat() when the Gap should be retained.",
        )

    overlap_length = _find_overlap_length(
        fragments[0].symbols,
        fragments[1].symbols,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
    )
    alphabet = promote_alphabet(*(fragment.alphabet for fragment in fragments))
    return combine_parts(
        [fragments[0].symbols + fragments[1].symbols[overlap_length:]],
        alphabet=alphabet,
        topology=Topology.LINEAR,
        strandedness=first.strandedness,
    )


def _validate_overlap_bounds(min_overlap: int, max_overlap: int | None) -> None:
    if isinstance(min_overlap, bool) or not isinstance(min_overlap, int) or min_overlap < 1:
        raise ConfigurationError("min_overlap must be a positive integer.")
    if max_overlap is not None and (
        isinstance(max_overlap, bool) or not isinstance(max_overlap, int) or max_overlap < 1
    ):
        raise ConfigurationError("max_overlap must be a positive integer or None.")
    if max_overlap is not None and min_overlap > max_overlap:
        raise ConfigurationError("min_overlap cannot exceed max_overlap.")


def _find_overlap_length(
    left: str,
    right: str,
    *,
    min_overlap: int,
    max_overlap: int | None,
) -> int:
    upper = min(len(left), len(right))
    if max_overlap is not None:
        upper = min(upper, max_overlap)
    for length in range(upper, min_overlap - 1, -1):
        if left[-length:] == right[:length]:
            return length
    raise ConfigurationError(
        "No exact suffix/prefix overlap satisfies the configured bounds.",
        code="OVERLAP_NOT_FOUND",
        context={
            "left_length": len(left),
            "right_length": len(right),
            "min_overlap": min_overlap,
            "max_overlap": max_overlap,
        },
    )


def _infer_text_alphabet(text: str) -> DNAAlphabet:
    return DNAAlphabet.STRICT if set(text) <= set("ACGT") else DNAAlphabet.IUPAC


__all__ = ["concat", "concat_overlap"]
