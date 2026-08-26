"""Exact and reverse-complement-aware record deduplication."""

from __future__ import annotations

from collections.abc import Hashable, Iterable
from typing import TypeAlias

from dnakit.core import DNARecord, DNASequence, DNASet, Gap, Topology
from dnakit.core._json import JSONValue
from dnakit.exceptions import ConfigurationError

from ._metadata import metadata_value_key
from .config import DedupEquivalence, DeduplicationConfig
from .results import DedupAction, DedupGroup, DeduplicationResult, Orientation

SequenceKey: TypeAlias = tuple[tuple[str | Gap, ...], object, object]
EquivalenceKey: TypeAlias = SequenceKey | frozenset[SequenceKey] | tuple[str, object, object]


def _least_rotation(symbols: str) -> tuple[str, int]:
    """Return the lexicographically least rotation and its first offset."""

    if not symbols:
        return symbols, 0
    length = len(symbols)
    doubled = symbols + symbols
    left, right, matched = 0, 1, 0
    while left < length and right < length and matched < length:
        left_symbol = doubled[left + matched]
        right_symbol = doubled[right + matched]
        if left_symbol == right_symbol:
            matched += 1
            continue
        if left_symbol > right_symbol:
            left += matched + 1
            if left <= right:
                left = right + 1
        else:
            right += matched + 1
            if right <= left:
                right = left + 1
        matched = 0
    best_offset = min(left, right)
    return doubled[best_offset : best_offset + length], best_offset


def _sequence_key(sequence: DNASequence) -> SequenceKey:
    # Alphabet is deliberately omitted: ACGT content has the same identity
    # whether its validated declaration is strict or IUPAC.
    return (sequence.parts, sequence.topology, sequence.strandedness)


def _equivalence_key(sequence: DNASequence, equivalence: DedupEquivalence) -> EquivalenceKey:
    forward = _sequence_key(sequence)
    if equivalence == "exact":
        return forward
    if equivalence in {"circular", "circular_reverse_complement"}:
        if sequence.topology is not Topology.CIRCULAR:
            raise ConfigurationError(
                "Circular equivalence requires every sequence to declare circular topology.",
                code="CIRCULAR_DEDUPLICATION_TOPOLOGY_REQUIRED",
            )
        if sequence.is_gapped:
            raise ConfigurationError(
                "Circular rotation equivalence does not support Gap objects.",
                code="CIRCULAR_DEDUPLICATION_GAP_UNSUPPORTED",
            )
        canonical, _ = _least_rotation(sequence.symbols)
        if equivalence == "circular_reverse_complement":
            reverse, _ = _least_rotation(sequence.reverse_complement().symbols)
            canonical = min(canonical, reverse)
        return canonical, sequence.topology, sequence.strandedness
    return frozenset((forward, _sequence_key(sequence.reverse_complement())))


def _materialize(records: Iterable[DNARecord]) -> tuple[DNARecord, ...]:
    try:
        materialized = tuple(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be an iterable of DNARecord objects.",
            code="INVALID_DATASET_RECORDS",
        ) from exc
    for index, record in enumerate(materialized):
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "Every dataset item must be a DNARecord.",
                code="INVALID_DATASET_RECORD",
                context={"input_index": index, "type": type(record).__name__},
            )
    return materialized


def _quality_score(record: DNARecord) -> float:
    quality = record.letter_annotations.get("phred_quality")
    if not quality:
        return float("-inf")
    return sum(quality) / len(quality)


def _representative_index(
    members: tuple[int, ...], records: tuple[DNARecord, ...], policy: str
) -> int:
    if policy == "first":
        return members[0]
    if policy == "last":
        return members[-1]
    return max(members, key=lambda index: (_quality_score(records[index]), -index))


def _conflict(
    members: tuple[int, ...],
    records: tuple[DNARecord, ...],
    field: str | None,
) -> tuple[bool, tuple[JSONValue, ...], int]:
    if field is None:
        return False, (), 0
    values: list[JSONValue] = []
    missing = 0
    for index in members:
        if field not in records[index].metadata:
            missing += 1
        else:
            values.append(records[index].metadata[field])
    unique_by_key = {metadata_value_key(value): value for value in values}
    unique_values = tuple(unique_by_key[key] for key in sorted(unique_by_key, key=repr))
    distinct_count = len(unique_values) + int(missing > 0)
    return distinct_count > 1, unique_values, missing


def _merge_metadata(representative: DNARecord, members: tuple[DNARecord, ...]) -> DNARecord:
    merged: dict[str, object] = dict(representative.metadata)
    keys = {key for member in members for key in member.metadata}
    for key in sorted(keys):
        present = [member.metadata[key] for member in members if key in member.metadata]
        if not present:
            continue
        first_key = metadata_value_key(present[0])
        if all(metadata_value_key(value) == first_key for value in present[1:]):
            merged.setdefault(key, present[0])
    return DNARecord(
        representative.sequence,
        representative.id,
        description=representative.description,
        features=representative.features,
        metadata=merged,
        letter_annotations=representative.letter_annotations,
    )


def _orientations(
    representative: DNARecord,
    members: tuple[int, ...],
    records: tuple[DNARecord, ...],
    equivalence: DedupEquivalence,
) -> tuple[Orientation, ...]:
    representative_key = _sequence_key(representative.sequence)
    reverse_key = _sequence_key(representative.sequence.reverse_complement())
    values: list[Orientation] = []
    for index in members:
        member_key = _sequence_key(records[index].sequence)
        orientation: Orientation = "forward"
        if equivalence == "reverse_complement" and member_key != representative_key:
            if member_key != reverse_key:
                raise AssertionError("Member does not match its reverse-complement group.")
            orientation = "reverse_complement"
        values.append(orientation)
    return tuple(values)


def _circular_transforms(
    representative: DNARecord,
    members: tuple[int, ...],
    records: tuple[DNARecord, ...],
    equivalence: DedupEquivalence,
) -> tuple[tuple[Orientation, ...], tuple[int | None, ...]]:
    if equivalence not in {"circular", "circular_reverse_complement"}:
        return _orientations(representative, members, records, equivalence), ()
    orientations: list[Orientation] = []
    offsets: list[int | None] = []
    for index in members:
        sequence = records[index].sequence
        forward, forward_offset = _least_rotation(sequence.symbols)
        orientation: Orientation = "forward"
        offset = forward_offset
        if equivalence == "circular_reverse_complement":
            reverse, reverse_offset = _least_rotation(sequence.reverse_complement().symbols)
            if reverse < forward:
                orientation = "reverse_complement"
                offset = reverse_offset
        orientations.append(orientation)
        offsets.append(offset)
    return tuple(orientations), tuple(offsets)


def deduplicate(
    records: Iterable[DNARecord],
    *,
    equivalence: DedupEquivalence = "exact",
    config: DeduplicationConfig | None = None,
) -> DeduplicationResult:
    """Deduplicate normalized records while preserving a complete group audit."""

    if equivalence not in {
        "exact",
        "reverse_complement",
        "circular",
        "circular_reverse_complement",
    }:
        raise ConfigurationError(
            "Unknown deduplication equivalence.",
            code="INVALID_DEDUPLICATION_EQUIVALENCE",
            context={"equivalence": equivalence},
        )
    resolved = DeduplicationConfig() if config is None else config
    if not isinstance(resolved, DeduplicationConfig):
        raise ConfigurationError(
            "config must be DeduplicationConfig or None.",
            code="INVALID_DEDUPLICATION_CONFIG",
        )
    materialized = _materialize(records)
    grouped: dict[Hashable, list[int]] = {}
    for index, record in enumerate(materialized):
        grouped.setdefault(_equivalence_key(record.sequence, equivalence), []).append(index)

    output: list[DNARecord] = []
    audit_groups: list[DedupGroup] = []
    conflict_count = 0
    for group_index, raw_members in enumerate(grouped.values()):
        members = tuple(raw_members)
        representative_index = _representative_index(
            members, materialized, resolved.representative_policy
        )
        representative = materialized[representative_index]
        orientations, rotation_offsets = _circular_transforms(
            representative, members, materialized, equivalence
        )
        conflict, conflict_values, missing_count = _conflict(
            members, materialized, resolved.conflict_field
        )
        if conflict:
            conflict_count += 1
            if resolved.conflict_policy == "error":
                raise ConfigurationError(
                    "A duplicate group contains conflicting metadata values.",
                    code="DEDUPLICATION_METADATA_CONFLICT",
                    context={
                        "group_index": group_index,
                        "member_ids": [materialized[index].id for index in members],
                        "conflict_field": resolved.conflict_field,
                    },
                    hint="Choose drop_group, keep_representative, or keep_all explicitly.",
                )
        action: DedupAction = "deduplicated"
        if conflict and resolved.conflict_policy == "drop_group":
            action = "dropped"
        elif conflict and resolved.conflict_policy == "keep_all":
            action = "kept_all"
            output.extend(materialized[index] for index in members)
        else:
            chosen = representative
            if resolved.merge_metadata:
                chosen = _merge_metadata(
                    representative, tuple(materialized[index] for index in members)
                )
            output.append(chosen)
        audit_groups.append(
            DedupGroup(
                group_index=group_index,
                representative_id=representative.id,
                member_ids=tuple(materialized[index].id for index in members),
                orientations=orientations,
                conflict=conflict,
                conflict_values=conflict_values,
                missing_conflict_value_count=missing_count,
                action=action,
                rotation_offsets=rotation_offsets,
                rotation_offset_definition=(
                    "zero-based left-rotation offset to the lexicographically least "
                    "representation after the reported orientation"
                    if equivalence in {"circular", "circular_reverse_complement"}
                    else None
                ),
            )
        )
    dataset = DNASet(output)
    return DeduplicationResult(
        records=dataset,
        groups=tuple(audit_groups),
        equivalence=equivalence,
        representative_policy=resolved.representative_policy,
        conflict_field=resolved.conflict_field,
        conflict_policy=resolved.conflict_policy,
        merge_metadata=resolved.merge_metadata,
        input_count=len(materialized),
        output_count=len(dataset),
        duplicate_count=len(materialized) - len(audit_groups),
        conflicted_group_count=conflict_count,
        removed_count=len(materialized) - len(dataset),
    )


__all__ = ["deduplicate"]
