"""Bounded palindrome, inverted-repeat, tandem-repeat, and microsatellite scans."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from dnakit.exceptions import ConfigurationError
from dnakit.patterns._shared import (
    SequenceInput,
    build_result,
    iupac_compatible,
    pattern_provenance,
    require_linear,
    resolve_sequence,
    segment_location,
    segments,
    validate_bool,
    validate_positive_int,
)
from dnakit.patterns.results import (
    InvertedRepeatHit,
    PalindromeHit,
    PatternResult,
    TandemRepeatHit,
)

_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)
_MAX_TANDEM_UNIT_LENGTH = 100_000


def _reverse_complement(text: str) -> str:
    return text.translate(_COMPLEMENT)[::-1]


def _compatible(left: str, right: str) -> bool:
    return len(left) == len(right) and all(
        iupac_compatible(left_base, right_base)
        for left_base, right_base in zip(left, right, strict=True)
    )


def find_reverse_complement_palindromes(
    value: SequenceInput,
    *,
    min_length: int = 4,
    max_length: int = 100,
    maximal_per_start: bool = True,
    max_comparisons: int = 5_000_000,
    max_comparison_cells: int = 50_000_000,
    max_matches: int = 100_000,
) -> PatternResult[PalindromeHit]:
    """Find IUPAC-compatible reverse-complement palindromes within fragments."""

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "reverse-complement palindrome detection")
    validate_positive_int(min_length, "min_length")
    validate_positive_int(max_length, "max_length")
    validate_positive_int(max_comparisons, "max_comparisons")
    validate_positive_int(max_comparison_cells, "max_comparison_cells")
    validate_positive_int(max_matches, "max_matches")
    validate_bool(maximal_per_start, "maximal_per_start")
    if min_length > max_length:
        raise ConfigurationError("min_length cannot exceed max_length.")

    hits: list[PalindromeHit] = []
    comparisons = 0
    comparison_cells = 0
    truncated = False
    for item in segments(sequence):
        for start in range(len(item.text)):
            matches: list[tuple[int, str]] = []
            upper = min(max_length, len(item.text) - start)
            for length in range(min_length, upper + 1):
                comparisons += 1
                comparison_cells += length
                if comparisons > max_comparisons:
                    raise ConfigurationError("Palindrome scan exceeded max_comparisons.")
                if comparison_cells > max_comparison_cells:
                    raise ConfigurationError("Palindrome scan exceeded max_comparison_cells.")
                candidate = item.text[start : start + length]
                if _compatible(candidate, _reverse_complement(candidate)):
                    matches.append((length, candidate))
            selected = matches[-1:] if maximal_per_start else matches
            for length, candidate in selected:
                if len(hits) >= max_matches:
                    truncated = True
                    break
                symbol_location, coordinate_location = segment_location(item, start, start + length)
                hits.append(
                    PalindromeHit(
                        sequence=candidate,
                        symbol_location=symbol_location,
                        coordinate_location=coordinate_location,
                        length=length,
                    )
                )
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="reverse_complement_palindrome_scan",
        method="bounded_exhaustive_iupac_compatibility",
        algorithm_version="1.0",
        parameters={
            "min_length": min_length,
            "max_length": max_length,
            "maximal_per_start": maximal_per_start,
            "ambiguity_rule": "IUPAC-set-intersection",
            "comparisons": comparisons,
            "max_comparisons": max_comparisons,
            "comparison_cells": comparison_cells,
            "max_comparison_cells": max_comparison_cells,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(reimplementation=False),
    )
    return cast(PatternResult[PalindromeHit], result)


def find_inverted_repeats(
    value: SequenceInput,
    *,
    min_arm_length: int = 4,
    max_arm_length: int = 50,
    min_loop_length: int = 0,
    max_loop_length: int = 100,
    max_comparisons: int = 5_000_000,
    max_comparison_cells: int = 50_000_000,
    max_matches: int = 100_000,
) -> PatternResult[InvertedRepeatHit]:
    """Find exact/IUPAC-compatible inverted-repeat arms with bounded loops."""

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "inverted-repeat detection")
    validate_positive_int(min_arm_length, "min_arm_length")
    validate_positive_int(max_arm_length, "max_arm_length")
    if (
        isinstance(min_loop_length, bool)
        or not isinstance(min_loop_length, int)
        or min_loop_length < 0
    ):
        raise ConfigurationError("min_loop_length must be a non-negative integer.")
    if (
        isinstance(max_loop_length, bool)
        or not isinstance(max_loop_length, int)
        or max_loop_length < 0
    ):
        raise ConfigurationError("max_loop_length must be a non-negative integer.")
    validate_positive_int(max_comparisons, "max_comparisons")
    validate_positive_int(max_comparison_cells, "max_comparison_cells")
    validate_positive_int(max_matches, "max_matches")
    if min_arm_length > max_arm_length or min_loop_length > max_loop_length:
        raise ConfigurationError("Minimum repeat bounds cannot exceed maximum bounds.")

    hits: list[InvertedRepeatHit] = []
    comparisons = 0
    comparison_cells = 0
    truncated = False
    for item in segments(sequence):
        text = item.text
        for left_start in range(len(text)):
            for arm_length in range(min_arm_length, max_arm_length + 1):
                left_end = left_start + arm_length
                if left_end > len(text):
                    break
                left = text[left_start:left_end]
                expected_right = _reverse_complement(left)
                for loop_length in range(min_loop_length, max_loop_length + 1):
                    right_start = left_end + loop_length
                    right_end = right_start + arm_length
                    if right_end > len(text):
                        break
                    comparisons += 1
                    comparison_cells += arm_length
                    if comparisons > max_comparisons:
                        raise ConfigurationError("Inverted-repeat scan exceeded max_comparisons.")
                    if comparison_cells > max_comparison_cells:
                        raise ConfigurationError(
                            "Inverted-repeat scan exceeded max_comparison_cells."
                        )
                    right = text[right_start:right_end]
                    if not _compatible(right, expected_right):
                        continue
                    if len(hits) >= max_matches:
                        truncated = True
                        break
                    symbol_location, coordinate_location = segment_location(
                        item, left_start, right_end
                    )
                    left_location, _ = segment_location(item, left_start, left_end)
                    loop_location, _ = segment_location(item, left_end, right_start)
                    right_location, _ = segment_location(item, right_start, right_end)
                    hits.append(
                        InvertedRepeatHit(
                            left_arm=left,
                            right_arm=right,
                            arm_length=arm_length,
                            loop_length=loop_length,
                            symbol_location=symbol_location,
                            coordinate_location=coordinate_location,
                            left_symbol_location=left_location,
                            loop_symbol_location=loop_location,
                            right_symbol_location=right_location,
                        )
                    )
                if truncated:
                    break
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="inverted_repeat_scan",
        method="bounded_exact_arm_loop_enumeration",
        algorithm_version="1.0",
        parameters={
            "min_arm_length": min_arm_length,
            "max_arm_length": max_arm_length,
            "min_loop_length": min_loop_length,
            "max_loop_length": max_loop_length,
            "mismatches": 0,
            "indels": 0,
            "ambiguity_rule": "IUPAC-set-intersection",
            "comparisons": comparisons,
            "max_comparisons": max_comparisons,
            "comparison_cells": comparison_cells,
            "max_comparison_cells": max_comparison_cells,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(reimplementation=False),
    )
    return cast(PatternResult[InvertedRepeatHit], result)


def _repeat_thresholds(
    min_unit_length: int,
    max_unit_length: int,
    min_repeats: int,
    values: Mapping[int, int] | None,
) -> dict[int, int]:
    if values is None:
        return {length: min_repeats for length in range(min_unit_length, max_unit_length + 1)}
    resolved: dict[int, int] = {}
    expected_count = max_unit_length - min_unit_length + 1
    if len(values) != expected_count:
        raise ConfigurationError(
            "min_repeats_by_unit must define every accepted unit length.",
            context={"expected_count": expected_count, "actual_count": len(values)},
        )
    for unit_length in range(min_unit_length, max_unit_length + 1):
        try:
            threshold = values[unit_length]
        except KeyError as exc:
            raise ConfigurationError(
                "min_repeats_by_unit must define every accepted unit length.",
                context={"missing_unit_length": unit_length},
            ) from exc
        validate_positive_int(threshold, f"min_repeats_by_unit[{unit_length}]")
        if threshold < 2:
            raise ConfigurationError("Every tandem-repeat threshold must be at least 2.")
        resolved[unit_length] = threshold
    return resolved


def find_tandem_repeats(
    value: SequenceInput,
    *,
    min_unit_length: int = 1,
    max_unit_length: int = 20,
    min_repeats: int = 2,
    min_repeats_by_unit: Mapping[int, int] | None = None,
    overlapping: bool = False,
    max_comparisons: int = 5_000_000,
    max_comparison_cells: int = 50_000_000,
    max_matches: int = 100_000,
) -> PatternResult[TandemRepeatHit]:
    """Find maximal exact repeats using the smallest accepted period at each start."""

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "tandem-repeat detection")
    validate_positive_int(min_unit_length, "min_unit_length")
    validate_positive_int(max_unit_length, "max_unit_length")
    validate_positive_int(min_repeats, "min_repeats")
    validate_positive_int(max_comparisons, "max_comparisons")
    validate_positive_int(max_comparison_cells, "max_comparison_cells")
    validate_positive_int(max_matches, "max_matches")
    validate_bool(overlapping, "overlapping")
    if min_repeats < 2:
        raise ConfigurationError("min_repeats must be at least 2.")
    if min_unit_length > max_unit_length:
        raise ConfigurationError("min_unit_length cannot exceed max_unit_length.")
    if max_unit_length > _MAX_TANDEM_UNIT_LENGTH:
        raise ConfigurationError(
            "max_unit_length exceeds the hard tandem-repeat unit limit.",
            code="TANDEM_REPEAT_UNIT_LIMIT",
            context={"maximum": _MAX_TANDEM_UNIT_LENGTH},
        )
    longest_fragment = max((len(item.text) for item in segments(sequence)), default=0)
    useful_maximum = longest_fragment // min_repeats
    if min_unit_length > useful_maximum:
        effective_max_unit_length = min_unit_length - 1
    else:
        effective_max_unit_length = min(max_unit_length, useful_maximum)
    thresholds = _repeat_thresholds(
        min_unit_length, max_unit_length, min_repeats, min_repeats_by_unit
    )

    hits: list[TandemRepeatHit] = []
    comparisons = 0
    comparison_cells = 0
    truncated = False
    for item in segments(sequence):
        start = 0
        while start < len(item.text):
            found: tuple[str, int, int] | None = None
            for unit_length in range(min_unit_length, effective_max_unit_length + 1):
                unit_end = start + unit_length
                minimum_end = start + unit_length * thresholds[unit_length]
                if minimum_end > len(item.text):
                    continue
                unit = item.text[start:unit_end]
                repeat_count = 1
                cursor = unit_end
                while cursor + unit_length <= len(item.text):
                    comparisons += 1
                    comparison_cells += unit_length
                    if comparisons > max_comparisons:
                        raise ConfigurationError("Tandem-repeat scan exceeded max_comparisons.")
                    if comparison_cells > max_comparison_cells:
                        raise ConfigurationError(
                            "Tandem-repeat scan exceeded max_comparison_cells."
                        )
                    if item.text[cursor : cursor + unit_length] != unit:
                        break
                    repeat_count += 1
                    cursor += unit_length
                if repeat_count >= thresholds[unit_length]:
                    found = unit, repeat_count, cursor
                    break
            if found is None:
                start += 1
                continue
            unit, repeat_count, end = found
            if len(hits) >= max_matches:
                truncated = True
                break
            symbol_location, coordinate_location = segment_location(item, start, end)
            hits.append(
                TandemRepeatHit(
                    unit=unit,
                    unit_length=len(unit),
                    repeat_count=repeat_count,
                    sequence=item.text[start:end],
                    symbol_location=symbol_location,
                    coordinate_location=coordinate_location,
                )
            )
            start += 1 if overlapping else end - start
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="tandem_repeat_scan",
        method="smallest_period_exact_run",
        algorithm_version="1.0",
        parameters={
            "min_unit_length": min_unit_length,
            "max_unit_length": max_unit_length,
            "effective_max_unit_length": effective_max_unit_length,
            "min_repeats": min_repeats,
            "min_repeats_by_unit": {
                str(unit_length): threshold for unit_length, threshold in sorted(thresholds.items())
            },
            "overlapping": overlapping,
            "mismatches": 0,
            "indels": 0,
            "ambiguity_rule": "literal-IUPAC-symbol-equality",
            "period_policy": "smallest-accepted-at-each-start",
            "comparisons": comparisons,
            "max_comparisons": max_comparisons,
            "comparison_cells": comparison_cells,
            "max_comparison_cells": max_comparison_cells,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(reimplementation=False),
    )
    return cast(PatternResult[TandemRepeatHit], result)


def find_microsatellites(
    value: SequenceInput,
    *,
    min_repeats_by_unit: Mapping[int, int] | None = None,
    overlapping: bool = False,
    max_comparisons: int = 5_000_000,
    max_comparison_cells: int = 50_000_000,
    max_matches: int = 100_000,
) -> PatternResult[TandemRepeatHit]:
    """Find exact 1--6 bp STRs with explicit per-unit repeat thresholds."""

    thresholds = (
        {1: 6, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3}
        if min_repeats_by_unit is None
        else dict(min_repeats_by_unit)
    )
    result = find_tandem_repeats(
        value,
        min_unit_length=1,
        max_unit_length=6,
        min_repeats=2,
        min_repeats_by_unit=thresholds,
        overlapping=overlapping,
        max_comparisons=max_comparisons,
        max_comparison_cells=max_comparison_cells,
        max_matches=max_matches,
    )
    parameters = dict(result.parameters)
    parameters.update(
        {
            "definition": "exact-short-tandem-repeat-1-to-6-bp",
            "interruptions_allowed": False,
        }
    )
    return PatternResult(
        name="microsatellite_scan",
        method=result.method,
        algorithm_version=result.algorithm_version,
        sequence_id=result.sequence_id,
        parameters=type(result.parameters)(parameters),
        hits=result.hits,
        inspected_symbol_count=result.inspected_symbol_count,
        gap_count=result.gap_count,
        unknown_gap_count=result.unknown_gap_count,
        max_matches=result.max_matches,
        truncated=result.truncated,
        coordinate_system=result.coordinate_system,
        gap_policy=result.gap_policy,
        topology=result.topology,
        provenance=result.provenance,
        issues=result.issues,
    )


__all__ = [
    "find_inverted_repeats",
    "find_microsatellites",
    "find_reverse_complement_palindromes",
    "find_tandem_repeats",
]
