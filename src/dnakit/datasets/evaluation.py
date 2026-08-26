"""Exhaustive bounded leakage detection and split-quality metrics."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from dnakit.core import DNA, DNASet
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

from ._advanced_shared import pair_similarity, validate_pairwise_sequences
from ._metadata import metadata_value_key
from .config import LeakageConfig
from .results import (
    LeakageEvent,
    LeakageReport,
    SplitAssignment,
    SplitQualityResult,
)


def _validate_splits(
    splits: Mapping[str, DNA | DNASet],
) -> tuple[tuple[str, DNASet], ...]:
    if not isinstance(splits, Mapping) or len(splits) < 2:
        raise ConfigurationError("splits must map at least two names to DNA values.")
    items = tuple(splits.items())
    if any(
        not isinstance(name, str) or not name.strip() or not isinstance(dataset, (DNA, DNASet))
        for name, dataset in items
    ):
        raise ConfigurationError("Split names must be non-empty and values must be DNA.")
    return tuple(
        (name, dataset.dataset if isinstance(dataset, DNA) else dataset) for name, dataset in items
    )


def detect_leakage(
    splits: Mapping[str, DNA | DNASet],
    *,
    config: LeakageConfig | None = None,
) -> LeakageReport:
    """Exhaustively compare every cross-split pair under hard record/pair/event limits."""

    resolved = LeakageConfig() if config is None else config
    if not isinstance(resolved, LeakageConfig):
        raise ConfigurationError("config must be LeakageConfig.")
    items = _validate_splits(splits)
    total = sum(len(dataset) for _, dataset in items)
    if total > resolved.max_records:
        raise ConfigurationError(
            "Leakage input exceeds max_records.",
            code="LEAKAGE_RECORD_LIMIT",
            context={"record_count": total, "max_records": resolved.max_records},
        )
    pair_total = sum(
        len(left) * len(right)
        for position, (_, left) in enumerate(items)
        for _, right in items[position + 1 :]
    )
    if pair_total > resolved.max_cross_pairs:
        raise ConfigurationError(
            "Leakage cross-pair workload exceeds max_cross_pairs.",
            code="LEAKAGE_PAIR_LIMIT",
            context={"cross_pair_count": pair_total, "max_cross_pairs": resolved.max_cross_pairs},
        )
    all_records = tuple(record for _, dataset in items for record in dataset)
    validate_pairwise_sequences(all_records, operation="detect_leakage")
    events: list[LeakageEvent] = []
    exact_count = 0
    high_count = 0
    for position, (left_name, left_set) in enumerate(items):
        for right_name, right_set in items[position + 1 :]:
            for left_index, left in enumerate(left_set):
                for right_index, right in enumerate(right_set):
                    exact = left.sequence.symbols == right.sequence.symbols
                    similarity = pair_similarity(
                        left,
                        right,
                        method=resolved.method,
                        k=resolved.k,
                        canonical=resolved.canonical,
                        max_alignment_cells=resolved.max_alignment_cells,
                    )
                    if exact or similarity >= resolved.threshold:
                        exact_count += int(exact)
                        high_count += int(not exact and similarity >= resolved.threshold)
                        if len(events) >= resolved.max_events:
                            raise ConfigurationError(
                                "Leakage event count exceeds max_events; "
                                "no partial report returned.",
                                code="LEAKAGE_EVENT_LIMIT",
                                context={"max_events": resolved.max_events},
                            )
                        events.append(
                            LeakageEvent(
                                left_name,
                                right_name,
                                left.id,
                                right.id,
                                left_index,
                                right_index,
                                similarity,
                                exact,
                            )
                        )
    strictness = (
        "exhaustive cross-pair exact dynamic-programming comparison"
        if resolved.method in {"identity", "edit"}
        else "exhaustive cross-pair deterministic k-mer/fingerprint comparison; "
        "may miss biologically similar pairs not captured by configured features"
    )
    return LeakageReport(
        tuple(events),
        resolved.method,
        resolved.threshold,
        strictness,
        exact_count,
        high_count,
        pair_total,
        resolved.max_records,
        resolved.max_cross_pairs,
        resolved.max_events,
        False,
        resolved.k,
        resolved.canonical,
        resolved.max_alignment_cells,
    )


def evaluate_split_quality(
    records: DNA | DNASet,
    assignments: Sequence[SplitAssignment],
    *,
    target_ratios: Mapping[str, float],
    label_key: str | None = None,
    group_keys: Sequence[str] = (),
    leakage_report: LeakageReport | None = None,
) -> SplitQualityResult:
    """Evaluate ratio error, label total variation, group leakage, and pair leakage."""

    if not isinstance(records, (DNA, DNASet)):
        raise ConfigurationError("records must be DNA or DNASet.")
    dataset = records.dataset if isinstance(records, DNA) else records
    assignment_tuple = tuple(assignments)
    if len(assignment_tuple) != len(dataset):
        raise ConfigurationError("assignments must align one-to-one with records.")
    by_index: dict[int, str] = {}
    for item in assignment_tuple:
        if not isinstance(item, SplitAssignment) or item.input_index in by_index:
            raise ConfigurationError("assignments contain invalid or duplicate input indices.")
        if item.input_index >= len(dataset) or dataset[item.input_index].id != item.record_id:
            raise ConfigurationError("assignment index/ID does not align with records.")
        by_index[item.input_index] = item.split
    validated_ratios: dict[str, float] = {}
    for name, value in target_ratios.items():
        if (
            not isinstance(name, str)
            or not name.strip()
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0.0
        ):
            raise ConfigurationError(
                "target_ratios must map non-empty names to finite positive values."
            )
        validated_ratios[name] = float(value)
    split_names = tuple(validated_ratios)
    if set(item.split for item in assignment_tuple) - set(split_names):
        raise ConfigurationError("assignments contain names absent from target_ratios.")
    total_ratio = math.fsum(validated_ratios.values())
    if not math.isclose(total_ratio, 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise ConfigurationError("target_ratios must sum to one.")
    counts = {name: sum(item.split == name for item in assignment_tuple) for name in split_names}
    total = len(dataset)
    achieved = {name: counts[name] / total if total else 0.0 for name in split_names}
    deviations = {name: abs(achieved[name] - validated_ratios[name]) for name in split_names}
    label_tv: dict[str, float] = {name: 0.0 for name in split_names}
    if label_key is not None and total:
        global_counts: dict[object, int] = {}
        split_counts: dict[str, dict[object, int]] = {name: {} for name in split_names}
        for index, record in enumerate(dataset):
            if label_key not in record.metadata:
                raise ConfigurationError("A record is missing label metadata.")
            key = metadata_value_key(record.metadata[label_key])
            global_counts[key] = global_counts.get(key, 0) + 1
            name = by_index[index]
            split_counts[name][key] = split_counts[name].get(key, 0) + 1
        for name in split_names:
            denominator = counts[name]
            if denominator:
                label_tv[name] = 0.5 * math.fsum(
                    abs(split_counts[name].get(key, 0) / denominator - value / total)
                    for key, value in global_counts.items()
                )
    validated_group_keys = tuple(group_keys)
    if any(not isinstance(key, str) or not key.strip() for key in validated_group_keys):
        raise ConfigurationError("group_keys must contain non-empty strings.")
    if len(validated_group_keys) != len(set(validated_group_keys)):
        raise ConfigurationError("group_keys must be unique.")
    if label_key is not None and (not isinstance(label_key, str) or not label_key.strip()):
        raise ConfigurationError("label_key must be non-empty or None.")
    if leakage_report is not None and not isinstance(leakage_report, LeakageReport):
        raise ConfigurationError("leakage_report must be LeakageReport or None.")
    group_leaks = 0
    for key in validated_group_keys:
        seen: dict[object, set[str]] = {}
        for index, record in enumerate(dataset):
            if key not in record.metadata:
                raise ConfigurationError("A record is missing group metadata.")
            group_value = metadata_value_key(record.metadata[key])
            seen.setdefault(group_value, set()).add(by_index[index])
        group_leaks += sum(len(splits) > 1 for splits in seen.values())
    max_ratio = max(deviations.values(), default=0.0)
    max_label = max(label_tv.values(), default=0.0)
    leakage_count = None if leakage_report is None else len(leakage_report.events)
    leakage_penalty = 0.0 if leakage_count is None else min(1.0, leakage_count / max(total, 1))
    score = max(0.0, 1.0 - max_ratio - max_label - min(1.0, group_leaks) - leakage_penalty)
    return SplitQualityResult(
        total,
        FrozenDict(counts),
        FrozenDict(validated_ratios),
        FrozenDict(achieved),
        FrozenDict(deviations),
        max_ratio,
        label_key,
        FrozenDict(label_tv),
        validated_group_keys,
        group_leaks,
        leakage_count,
        score,
        "max(0,1-max_ratio_deviation-max_label_TV-min(1,group_leaks)-leakage_events/max(N,1))",
        FrozenDict(
            {
                "empty_split_label_tv": 0.0,
                "group_leak_definition": "one typed metadata group present in more than one split",
            }
        ),
    )


__all__ = ["detect_leakage", "evaluate_split_quality"]
