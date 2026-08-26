"""Temporal and explicitly heuristic multi-constraint splitting."""

from __future__ import annotations

import math
import random
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from itertools import pairwise
from typing import cast

from dnakit.core import DNARecord, DNASet
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

from ._advanced_shared import (
    UnionFind,
    ensure_pair_limit,
    materialize_limited,
    pair_similarity,
    validate_pairwise_sequences,
)
from ._metadata import metadata_value_key
from .config import JointSplitConfig, TemporalSplitConfig
from .results import (
    JointSplitResult,
    SplitAssignment,
    SplitSubset,
    TemporalSplitResult,
)


def _allocate_counts(size: int, ratios: Mapping[str, object]) -> dict[str, int]:
    numeric = {name: float(cast(int | float, value)) for name, value in ratios.items()}
    exact = {name: size * ratio for name, ratio in numeric.items()}
    counts = {name: math.floor(value) for name, value in exact.items()}
    remainder = size - sum(counts.values())
    order = sorted(numeric, key=lambda name: (-(exact[name] - counts[name]), name))
    for name in order[:remainder]:
        counts[name] += 1
    return counts


def _parse_iso(value: object, *, index: int, key: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            "Temporal metadata must be a non-empty ISO-8601 string.",
            code="INVALID_TEMPORAL_METADATA",
            context={"input_index": index, "metadata_key": key},
        )
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        result = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ConfigurationError(
            "Temporal metadata is not valid ISO-8601.",
            code="INVALID_TEMPORAL_METADATA",
            context={"input_index": index, "value": value},
        ) from exc
    if result.tzinfo is not None:
        result = result.astimezone(timezone.utc).replace(tzinfo=None)
    return result


def temporal_split(
    records: Iterable[DNARecord],
    *,
    config: TemporalSplitConfig | None = None,
) -> TemporalSplitResult:
    """Assign oldest records first using explicit cutoffs or target chronological ratios."""

    resolved = TemporalSplitConfig() if config is None else config
    if not isinstance(resolved, TemporalSplitConfig):
        raise ConfigurationError("config must be TemporalSplitConfig.")
    materialized = materialize_limited(records, max_records=resolved.max_records)
    dated: list[tuple[datetime, int]] = []
    for index, record in enumerate(materialized):
        if resolved.metadata_key not in record.metadata:
            raise ConfigurationError(
                "A record is missing temporal metadata.",
                code="TEMPORAL_METADATA_MISSING",
                context={"input_index": index, "record_id": record.id},
            )
        dated.append(
            (
                _parse_iso(
                    record.metadata[resolved.metadata_key], index=index, key=resolved.metadata_key
                ),
                index,
            )
        )
    dated.sort(key=lambda item: (item[0], item[1]))
    split_names = tuple(resolved.ratios)
    assignment: dict[int, str] = {}
    cutoff_strings: tuple[str, ...]
    if resolved.cutoffs is None:
        counts = _allocate_counts(len(materialized), resolved.ratios)
        cursor = 0
        for name in split_names:
            for _, index in dated[cursor : cursor + counts[name]]:
                assignment[index] = name
            cursor += counts[name]
        cutoff_strings = tuple(
            dated[sum(counts[name] for name in split_names[:position]) - 1][0].isoformat()
            for position in range(1, len(split_names))
            if sum(counts[name] for name in split_names[:position]) > 0
        )
    else:
        parsed_cutoffs = tuple(
            _parse_iso(value, index=index, key="cutoffs")
            for index, value in enumerate(resolved.cutoffs)
        )
        if any(right <= left for left, right in pairwise(parsed_cutoffs)):
            raise ConfigurationError(
                "Temporal cutoffs must be strictly increasing.",
                code="INVALID_TEMPORAL_CUTOFFS",
            )
        for date, index in dated:
            split_position = sum(date > cutoff for cutoff in parsed_cutoffs)
            assignment[index] = split_names[split_position]
        cutoff_strings = resolved.cutoffs
        counts = {name: sum(value == name for value in assignment.values()) for name in split_names}
    indices = {
        name: [index for index in range(len(materialized)) if assignment[index] == name]
        for name in split_names
    }
    if not resolved.preserve_order:
        chronological_rank = {index: rank for rank, (_, index) in enumerate(dated)}
        for values in indices.values():
            values.sort(key=chronological_rank.__getitem__)
    return TemporalSplitResult(
        tuple(
            SplitSubset(name, DNASet(materialized[index] for index in indices[name]))
            for name in split_names
        ),
        tuple(
            SplitAssignment(index, record.id, assignment[index])
            for index, record in enumerate(materialized)
        ),
        FrozenDict(counts),
        cast(FrozenDict, resolved.ratios),
        resolved.metadata_key,
        cutoff_strings,
        "chronological-explicit-cutoffs" if resolved.cutoffs else "chronological-largest-remainder",
        resolved.preserve_order,
        "naive timestamps are interpreted as UTC; aware timestamps are normalized to UTC",
        resolved.max_records,
    )


def _joint_units(
    records: tuple[DNARecord, ...],
    config: JointSplitConfig,
) -> tuple[tuple[tuple[int, ...], ...], int]:
    union = UnionFind(len(records))
    for key in config.group_keys:
        grouped: dict[object, int] = {}
        for index, record in enumerate(records):
            if key not in record.metadata:
                raise ConfigurationError(
                    "Joint split group metadata is missing.",
                    code="JOINT_SPLIT_METADATA_MISSING",
                    context={"input_index": index, "metadata_key": key},
                )
            typed_key = metadata_value_key(record.metadata[key])
            if typed_key in grouped:
                union.union(grouped[typed_key], index)
            else:
                grouped[typed_key] = index
    comparisons = 0
    if config.similarity_threshold is not None:
        validate_pairwise_sequences(records, operation="joint_split similarity constraint")
        comparisons = ensure_pair_limit(len(records), config.max_pairwise_comparisons)
        for left in range(len(records)):
            for right in range(left + 1, len(records)):
                value = pair_similarity(
                    records[left],
                    records[right],
                    method="kmer",
                    k=config.similarity_k,
                    canonical=config.similarity_canonical,
                    max_alignment_cells=1,
                )
                if value >= config.similarity_threshold:
                    union.union(left, right)
    return union.groups(), comparisons


def joint_split(
    records: Iterable[DNARecord],
    *,
    config: JointSplitConfig,
) -> JointSplitResult:
    """Greedily assign atomic constraint components; never claim global optimality."""

    if not isinstance(config, JointSplitConfig):
        raise ConfigurationError("config must be JointSplitConfig.")
    materialized = materialize_limited(records, max_records=config.max_records)
    units, comparisons = _joint_units(materialized, config)
    targets = _allocate_counts(len(materialized), config.ratios)
    split_names = tuple(config.ratios)
    counts = {name: 0 for name in split_names}
    assignment: dict[int, str] = {}
    generator = random.Random(config.seed)
    ordered = list(units)
    generator.shuffle(ordered)
    ordered.sort(key=lambda unit: -len(unit))
    label_totals: dict[object, int] = {}
    if config.label_key is not None:
        for index, record in enumerate(materialized):
            if config.label_key not in record.metadata:
                raise ConfigurationError(
                    "Joint split label metadata is missing.",
                    code="JOINT_SPLIT_METADATA_MISSING",
                    context={"input_index": index, "metadata_key": config.label_key},
                )
            label = metadata_value_key(record.metadata[config.label_key])
            label_totals[label] = label_totals.get(label, 0) + 1
    label_counts: dict[str, dict[object, int]] = {name: {} for name in split_names}
    for unit in ordered:
        unit_labels: dict[object, int] = {}
        if config.label_key is not None:
            for index in unit:
                label = metadata_value_key(materialized[index].metadata[config.label_key])
                unit_labels[label] = unit_labels.get(label, 0) + 1

        def objective(
            name: str,
            unit: tuple[int, ...] = unit,
            unit_labels: dict[object, int] = unit_labels,
        ) -> tuple[float, int]:
            projected = counts[name] + len(unit)
            ratio_penalty = ((projected - targets[name]) / max(targets[name], 1)) ** 2
            label_penalty = 0.0
            for label, total in label_totals.items():
                desired = total * float(cast(int | float, config.ratios[name]))
                current = label_counts[name].get(label, 0) + unit_labels.get(label, 0)
                label_penalty += ((current - desired) / max(desired, 1.0)) ** 2
            return ratio_penalty + label_penalty, split_names.index(name)

        chosen = min(split_names, key=objective)
        counts[chosen] += len(unit)
        for label, value in unit_labels.items():
            label_counts[chosen][label] = label_counts[chosen].get(label, 0) + value
        for index in unit:
            assignment[index] = chosen
    total = len(materialized)
    achieved = {name: (counts[name] / total if total else 0.0) for name in split_names}
    deviations = {
        name: abs(achieved[name] - float(cast(int | float, config.ratios[name])))
        for name in split_names
    }
    max_deviation = max(deviations.values(), default=0.0)
    feasible = max_deviation <= config.ratio_tolerance
    if not feasible and config.infeasible_policy == "error":
        raise ConfigurationError(
            "Atomic constraints make the requested split ratios infeasible within tolerance.",
            code="JOINT_SPLIT_INFEASIBLE",
            context={
                "counts": counts,
                "target_counts": targets,
                "max_ratio_deviation": max_deviation,
                "ratio_tolerance": config.ratio_tolerance,
            },
            hint="Choose infeasible_policy='relax' to return the audited greedy assignment.",
        )
    indices = {
        name: [index for index in range(total) if assignment[index] == name] for name in split_names
    }
    return JointSplitResult(
        tuple(
            SplitSubset(name, DNASet(materialized[index] for index in indices[name]))
            for name in split_names
        ),
        tuple(
            SplitAssignment(index, record.id, assignment[index])
            for index, record in enumerate(materialized)
        ),
        FrozenDict(counts),
        FrozenDict(targets),
        FrozenDict(achieved),
        config.group_keys,
        config.label_key,
        config.similarity_threshold,
        config.similarity_k,
        config.similarity_canonical,
        config.seed,
        "seeded-greedy-atomic-components-ratio-and-label-penalty",
        feasible,
        not feasible,
        ("ratio_tolerance",) if not feasible else (),
        max_deviation,
        config.ratio_tolerance,
        len(units),
        comparisons,
        config.max_records,
        config.max_pairwise_comparisons,
        config.infeasible_policy,
    )


__all__ = ["joint_split", "temporal_split"]
