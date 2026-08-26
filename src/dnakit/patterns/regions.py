"""CpG-island and Shannon low-complexity region detection."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import cast

from dnakit.core import DNASequence, Interval
from dnakit.core.gap import Gap
from dnakit.exceptions import ConfigurationError
from dnakit.patterns._shared import (
    SequenceInput,
    build_result,
    frozen,
    pattern_provenance,
    require_linear,
    resolve_sequence,
    segment_location,
    segments,
    validate_positive_int,
)
from dnakit.patterns.results import LowComplexityResult, PatternResult, RegionHit


@dataclass(frozen=True)
class _CandidateRegion:
    start: int
    end: int


def _validate_fraction(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
        or value > 1
    ):
        raise ConfigurationError(f"{name} must be finite and between 0 and 1.")
    return float(value)


def _merge_window(regions: list[_CandidateRegion], start: int, end: int) -> None:
    if regions and start <= regions[-1].end:
        previous = regions[-1]
        regions[-1] = _CandidateRegion(previous.start, max(previous.end, end))
    else:
        regions.append(_CandidateRegion(start, end))


def _cpg_metrics(text: str) -> tuple[float, float, int, int, int]:
    c_count = text.count("C")
    g_count = text.count("G")
    cpg_count = sum(text[index : index + 2] == "CG" for index in range(len(text) - 1))
    gc_fraction = (c_count + g_count) / len(text) if text else 0.0
    observed_expected = (
        (cpg_count * len(text)) / (c_count * g_count) if c_count and g_count else 0.0
    )
    return gc_fraction, observed_expected, cpg_count, c_count, g_count


def find_cpg_islands(
    value: SequenceInput,
    *,
    window_size: int = 200,
    step: int = 1,
    min_gc: float = 0.5,
    min_observed_expected: float = 0.6,
    min_region_length: int = 200,
    max_windows: int = 2_000_000,
    max_matches: int = 100_000,
) -> PatternResult[RegionHit]:
    """Find merged qualifying windows using a configurable classic CpG rule.

    A window qualifies when it is unambiguous, has ``GC >= min_gc``, and has
    ``observed CpG / expected CpG >= min_observed_expected`` where
    ``O/E = CpG * L / (C * G)``. Overlapping qualifying windows are merged.
    """

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "CpG-island detection")
    validate_positive_int(window_size, "window_size")
    validate_positive_int(step, "step")
    validate_positive_int(min_region_length, "min_region_length")
    validate_positive_int(max_windows, "max_windows")
    validate_positive_int(max_matches, "max_matches")
    resolved_gc = _validate_fraction(min_gc, "min_gc")
    if (
        isinstance(min_observed_expected, bool)
        or not isinstance(min_observed_expected, (int, float))
        or not math.isfinite(min_observed_expected)
        or min_observed_expected < 0
    ):
        raise ConfigurationError("min_observed_expected must be finite and non-negative.")

    hits: list[RegionHit] = []
    windows_checked = 0
    truncated = False
    for item in segments(sequence):
        candidates: list[_CandidateRegion] = []
        for start in range(0, len(item.text) - window_size + 1, step):
            windows_checked += 1
            if windows_checked > max_windows:
                raise ConfigurationError("CpG-island scan exceeded max_windows.")
            window = item.text[start : start + window_size]
            if set(window) - set("ACGT"):
                continue
            gc_fraction, observed_expected, _, _, _ = _cpg_metrics(window)
            if gc_fraction >= resolved_gc and observed_expected >= min_observed_expected:
                _merge_window(candidates, start, start + window_size)
        for candidate in candidates:
            length = candidate.end - candidate.start
            if length < min_region_length:
                continue
            if len(hits) >= max_matches:
                truncated = True
                break
            region = item.text[candidate.start : candidate.end]
            gc_fraction, observed_expected, cpg_count, c_count, g_count = _cpg_metrics(region)
            symbol_location, coordinate_location = segment_location(
                item, candidate.start, candidate.end
            )
            hits.append(
                RegionHit(
                    kind="cpg_island_candidate",
                    symbol_location=symbol_location,
                    coordinate_location=coordinate_location,
                    length=length,
                    score=observed_expected,
                    attributes=frozen(
                        {
                            "gc_fraction": gc_fraction,
                            "observed_expected": observed_expected,
                            "cpg_count": cpg_count,
                            "c_count": c_count,
                            "g_count": g_count,
                        }
                    ),
                )
            )
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="cpg_island_scan",
        method="merged_qualifying_windows",
        algorithm_version="1.0",
        parameters={
            "window_size": window_size,
            "step": step,
            "min_gc": resolved_gc,
            "min_observed_expected": float(min_observed_expected),
            "min_region_length": min_region_length,
            "observed_expected_formula": "CpG * L / (C * G)",
            "threshold_comparison": "greater-than-or-equal",
            "ambiguous_window_policy": "skip",
            "merge_policy": "overlapping-qualifying-windows",
            "windows_checked": windows_checked,
            "max_windows": max_windows,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="Gardiner-Garden and Frommer configurable CpG boundary",
            reference_version="DNAKit-explicit-thresholds-v1",
        ),
    )
    return cast(PatternResult[RegionHit], result)


def _shannon_entropy(text: str) -> float:
    counts = Counter(text)
    total = len(text)
    return -math.fsum(
        (count / total) * math.log2(count / total) for count in counts.values() if count
    )


def _mask_sequence(
    sequence: DNASequence, intervals: tuple[Interval, ...], symbol: str
) -> DNASequence:
    parts: list[str | Gap] = []
    symbol_offset = 0
    for part in sequence.parts:
        if isinstance(part, Gap):
            parts.append(part)
            continue
        masked = list(part)
        part_end = symbol_offset + len(part)
        for interval in intervals:
            overlap_start = max(interval.start, symbol_offset)
            overlap_end = min(interval.end, part_end)
            for position in range(overlap_start, overlap_end):
                masked[position - symbol_offset] = symbol
        parts.append("".join(masked))
        symbol_offset = part_end
    output_alphabet = sequence.alphabet if symbol in "ACGT" else "iupac"
    return DNASequence(
        parts,
        alphabet=output_alphabet,
        topology=sequence.topology,
        strandedness=sequence.strandedness,
    )


def find_low_complexity_regions(
    value: SequenceInput,
    *,
    window_size: int = 32,
    step: int = 1,
    max_entropy: float = 1.2,
    min_canonical_fraction: float = 1.0,
    min_region_length: int = 32,
    mask_symbol: str = "N",
    max_windows: int = 2_000_000,
    max_matches: int = 100_000,
) -> LowComplexityResult:
    """Mark merged low-Shannon-entropy windows and return a masked sequence.

    This implementation is explicitly an entropy-window method, not DUST or
    RepeatMasker. Entropy is calculated in bits over canonical bases after the
    canonical-fraction eligibility check.
    """

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "low-complexity detection")
    validate_positive_int(window_size, "window_size")
    validate_positive_int(step, "step")
    validate_positive_int(min_region_length, "min_region_length")
    validate_positive_int(max_windows, "max_windows")
    validate_positive_int(max_matches, "max_matches")
    if (
        isinstance(max_entropy, bool)
        or not isinstance(max_entropy, (int, float))
        or not math.isfinite(max_entropy)
        or max_entropy < 0
        or max_entropy > 2
    ):
        raise ConfigurationError("max_entropy must be finite and between 0 and 2 bits.")
    resolved_canonical_fraction = _validate_fraction(
        min_canonical_fraction, "min_canonical_fraction"
    )
    if (
        not isinstance(mask_symbol, str)
        or len(mask_symbol) != 1
        or mask_symbol not in "ACGTRYSWKMBDHVN"
    ):
        raise ConfigurationError("mask_symbol must be one uppercase DNA IUPAC symbol.")

    hits: list[RegionHit] = []
    windows_checked = 0
    truncated = False
    for item in segments(sequence):
        candidates: list[_CandidateRegion] = []
        for start in range(0, len(item.text) - window_size + 1, step):
            windows_checked += 1
            if windows_checked > max_windows:
                raise ConfigurationError("Low-complexity scan exceeded max_windows.")
            window = item.text[start : start + window_size]
            canonical = "".join(base for base in window if base in "ACGT")
            canonical_fraction = len(canonical) / len(window)
            if not canonical or canonical_fraction < resolved_canonical_fraction:
                continue
            entropy = _shannon_entropy(canonical)
            if entropy <= max_entropy:
                _merge_window(candidates, start, start + window_size)
        for candidate in candidates:
            length = candidate.end - candidate.start
            if length < min_region_length:
                continue
            if len(hits) >= max_matches:
                truncated = True
                break
            region = item.text[candidate.start : candidate.end]
            canonical = "".join(base for base in region if base in "ACGT")
            entropy = _shannon_entropy(canonical)
            symbol_location, coordinate_location = segment_location(
                item, candidate.start, candidate.end
            )
            hits.append(
                RegionHit(
                    kind="low_complexity",
                    symbol_location=symbol_location,
                    coordinate_location=coordinate_location,
                    length=length,
                    score=entropy,
                    attributes=frozen(
                        {
                            "entropy_bits": entropy,
                            "canonical_fraction": len(canonical) / len(region),
                        }
                    ),
                )
            )
        if truncated:
            break
    analysis = build_result(
        sequence,
        sequence_id,
        name="low_complexity_scan",
        method="shannon_entropy_windows",
        algorithm_version="1.0",
        parameters={
            "window_size": window_size,
            "step": step,
            "max_entropy": float(max_entropy),
            "entropy_unit": "bits",
            "min_canonical_fraction": resolved_canonical_fraction,
            "min_region_length": min_region_length,
            "mask_symbol": mask_symbol,
            "merge_policy": "overlapping-qualifying-windows",
            "windows_checked": windows_checked,
            "max_windows": max_windows,
            "is_dust": False,
            "is_repeatmasker": False,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="Shannon entropy low-complexity boundary",
            reference_version="DNAKit-window-v1",
        ),
    )
    intervals = tuple(
        hit.symbol_location for hit in hits if isinstance(hit.symbol_location, Interval)
    )
    return LowComplexityResult(
        analysis=cast(PatternResult[RegionHit], analysis),
        masked_sequence=_mask_sequence(sequence, intervals, mask_symbol),
        mask_symbol=mask_symbol,
    )


__all__ = ["find_cpg_islands", "find_low_complexity_regions"]
