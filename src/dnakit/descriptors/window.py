"""Native multi-descriptor sliding windows."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Iterator
from itertools import islice

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.descriptors._shared import (
    DescriptorAmbiguityPolicy,
    SequenceInput,
    canonical_runs,
    coerce_ambiguity_policy,
    reject_ambiguity,
    sequence_and_id,
    symbol_view,
    validate_bool,
    validate_positive_int,
)
from dnakit.descriptors.results import WindowDescriptorResult, WindowResult
from dnakit.exceptions import ConfigurationError

_SUPPORTED_DESCRIPTORS = frozenset({"gc", "entropy", "cpg"})


def _ranges(segment_ids: tuple[int, ...], *, cross_gaps: bool) -> Iterator[tuple[int, int]]:
    if not segment_ids:
        return
    if cross_gaps:
        yield 0, len(segment_ids)
        return
    start = 0
    for index in range(1, len(segment_ids)):
        if segment_ids[index] != segment_ids[index - 1]:
            yield start, index
            start = index
    yield start, len(segment_ids)


def _resolve_descriptors(descriptors: Iterable[str]) -> tuple[str, ...]:
    if isinstance(descriptors, (str, bytes)):
        raise ConfigurationError("descriptors must be an iterable of descriptor names.")
    try:
        resolved = tuple(islice(iter(descriptors), len(_SUPPORTED_DESCRIPTORS) + 1))
    except TypeError as exc:
        raise ConfigurationError("descriptors must be an iterable of descriptor names.") from exc
    if not resolved:
        raise ConfigurationError("At least one window descriptor is required.")
    if any(not isinstance(name, str) or name not in _SUPPORTED_DESCRIPTORS for name in resolved):
        raise ConfigurationError(
            "Unsupported window descriptor.",
            context={"descriptors": resolved},
            hint="Choose from: gc, entropy, cpg.",
        )
    if len(set(resolved)) != len(resolved):
        raise ConfigurationError("Window descriptor names must be unique.")
    return resolved


def _window_values(text: str, descriptors: tuple[str, ...], log_base: float) -> FrozenDict:
    canonical = tuple(symbol for symbol in text if symbol in "ACGT")
    values: dict[str, object] = {}
    if "gc" in descriptors:
        denominator = len(canonical)
        gc_count = canonical.count("G") + canonical.count("C")
        values.update(
            {
                "gc_fraction": gc_count / denominator if denominator else None,
                "at_fraction": (denominator - gc_count) / denominator if denominator else None,
                "canonical_base_denominator": denominator,
            }
        )
    if "entropy" in descriptors:
        counts = Counter(canonical)
        total = sum(counts.values())
        entropy = (
            -sum((count / total) * math.log(count / total, log_base) for count in counts.values())
            if total
            else 0.0
        )
        values.update(
            {
                "shannon_entropy": entropy,
                "entropy_observation_count": total,
            }
        )
    if "cpg" in descriptors:
        runs = tuple(canonical_runs(text))
        cpg_count = sum(run.count("CG") for run in runs)
        pair_denominator = sum(max(0, len(run) - 1) for run in runs)
        canonical_length = sum(len(run) for run in runs)
        c_count = sum(run.count("C") for run in runs)
        g_count = sum(run.count("G") for run in runs)
        values.update(
            {
                "cpg_count": cpg_count,
                "cpg_density": cpg_count / pair_denominator if pair_denominator else None,
                "cpg_observed_expected": (
                    cpg_count * canonical_length / (c_count * g_count)
                    if c_count and g_count
                    else None
                ),
                "cpg_pair_denominator": pair_denominator,
            }
        )
    return FrozenDict(values)


def window_descriptors(
    value: SequenceInput,
    descriptors: Iterable[str],
    *,
    window_size: int,
    step: int = 1,
    include_partial: bool = False,
    entropy_log_base: float = 2.0,
    ambiguity_policy: DescriptorAmbiguityPolicy | str = DescriptorAmbiguityPolicy.ERROR,
    cross_gaps: bool = False,
) -> WindowDescriptorResult:
    """Calculate GC, Shannon entropy, and/or CpG values in ordered windows.

    Windows use zero-based, half-open symbol coordinates. By default each gap
    starts a new window series, so no row crosses a gap. ``include_partial``
    emits trailing windows shorter than ``window_size`` at each series boundary.
    """

    sequence, sequence_id = sequence_and_id(value)
    policy = coerce_ambiguity_policy(ambiguity_policy)
    reject_ambiguity(sequence, policy)
    resolved_descriptors = _resolve_descriptors(descriptors)
    validate_positive_int(window_size, "window_size")
    validate_positive_int(step, "step")
    validate_bool(include_partial, "include_partial")
    validate_bool(cross_gaps, "cross_gaps")
    if (
        isinstance(entropy_log_base, bool)
        or not isinstance(entropy_log_base, (int, float))
        or not math.isfinite(entropy_log_base)
        or entropy_log_base <= 0
        or entropy_log_base == 1
    ):
        raise ConfigurationError(
            "entropy_log_base must be finite, positive, and different from one.",
            context={"entropy_log_base": entropy_log_base},
        )

    view = symbol_view(sequence)
    rows: list[WindowResult] = []
    range_segments = view.crossable_segment_ids if cross_gaps else view.segment_ids
    for range_start, range_end in _ranges(range_segments, cross_gaps=False):
        last_start = range_end if include_partial else range_end - window_size + 1
        for start in range(range_start, max(range_start, last_start), step):
            end = min(start + window_size, range_end)
            if end <= start or (not include_partial and end - start < window_size):
                continue
            first_segment = view.segment_ids[start]
            last_segment = view.segment_ids[end - 1]
            crossed_gap_count = last_segment - first_segment
            crossed_unknown = any(
                boundary in view.unknown_gap_boundaries
                for boundary in range(first_segment + 1, last_segment + 1)
            )
            coordinate_start = view.coordinate_positions[start]
            last_coordinate = view.coordinate_positions[end - 1]
            rows.append(
                WindowResult(
                    symbol_start=start,
                    symbol_end=end,
                    coordinate_start=coordinate_start,
                    coordinate_end=None if last_coordinate is None else last_coordinate + 1,
                    is_partial=end - start < window_size,
                    crossed_gap_count=crossed_gap_count,
                    crossed_unknown_gap=crossed_unknown,
                    values=_window_values(
                        view.text[start:end],
                        resolved_descriptors,
                        float(entropy_log_base),
                    ),
                )
            )

    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return WindowDescriptorResult(
        name="window_descriptors",
        method="sliding_symbol_window",
        sequence_id=sequence_id,
        ambiguity_policy=policy,
        cross_gaps=cross_gaps,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        window_size=window_size,
        step=step,
        include_partial=include_partial,
        entropy_log_base=float(entropy_log_base),
        descriptors=resolved_descriptors,
        windows=tuple(rows),
    )


__all__ = ["window_descriptors"]
