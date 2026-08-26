"""Feature- and letter-annotation-aware operations on immutable DNA records."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Literal, TypeAlias, cast

from dnakit.core import (
    CompoundLocation,
    DNAFeature,
    DNARecord,
    DNASequence,
    Interval,
    Location,
    Strand,
    UnresolvedLocation,
    reverse_strand_location,
)
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError, UnknownLengthError
from dnakit.ops.circular import canonical_origin, rotate
from dnakit.ops.edit import Edit, delete, insert, mask, substitute, trim

FeatureOverlapPolicy: TypeAlias = Literal["preserve", "truncate", "delete", "split", "unresolved"]
LetterAnnotationPolicy: TypeAlias = Literal["error", "drop"]
FeatureChangeAction: TypeAlias = Literal[
    "preserved",
    "shifted",
    "resized",
    "truncated",
    "split",
    "unresolved",
    "deleted",
    "reverse_complemented",
    "rotated",
]
LetterAnnotationAction: TypeAlias = Literal[
    "absent",
    "preserved",
    "updated",
    "dropped",
    "trimmed",
    "reversed",
    "rotated",
]
RecordOperation: TypeAlias = Literal[
    "insert",
    "delete",
    "substitute",
    "mask",
    "trim",
    "reverse_complement",
    "rotate",
    "canonical_origin",
]

_FEATURE_POLICIES = frozenset({"preserve", "truncate", "delete", "split", "unresolved"})
_ANNOTATION_POLICIES = frozenset({"error", "drop"})
_FEATURE_ACTIONS = frozenset(
    {
        "preserved",
        "shifted",
        "resized",
        "truncated",
        "split",
        "unresolved",
        "deleted",
        "reverse_complemented",
        "rotated",
    }
)
_ANNOTATION_ACTIONS = frozenset(
    {"absent", "preserved", "updated", "dropped", "trimmed", "reversed", "rotated"}
)


@dataclass(frozen=True, slots=True)
class FeatureChange:
    """One feature-coordinate decision made by a record operation."""

    feature_index: int
    feature_id: str | None
    feature_type: str
    action: FeatureChangeAction
    original_location: Location
    new_location: Location | None
    policy: FeatureOverlapPolicy
    affected_edit_indices: tuple[int, ...] = ()
    reason: str | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.feature_index, bool)
            or not isinstance(self.feature_index, int)
            or self.feature_index < 0
        ):
            raise ConfigurationError("FeatureChange feature_index must be non-negative.")
        if not isinstance(self.feature_type, str) or not self.feature_type.strip():
            raise ConfigurationError("FeatureChange feature_type must be non-empty.")
        if self.action not in _FEATURE_ACTIONS:
            raise ConfigurationError("FeatureChange action is invalid.")
        if self.policy not in _FEATURE_POLICIES:
            raise ConfigurationError("FeatureChange policy is invalid.")
        if any(
            isinstance(index, bool) or not isinstance(index, int) or index < 0
            for index in self.affected_edit_indices
        ):
            raise ConfigurationError("FeatureChange edit indices must be non-negative integers.")
        if self.reason is not None and (
            not isinstance(self.reason, str) or not self.reason.strip()
        ):
            raise ConfigurationError("FeatureChange reason must be non-empty or None.")

    def to_dict(self) -> dict[str, Any]:
        """Return this decision as JSON-compatible primitives."""

        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class RecordOperationResult:
    """A transformed record plus complete coordinate and annotation audit."""

    record: DNARecord
    operation: RecordOperation
    edits: tuple[Edit, ...]
    feature_policy: FeatureOverlapPolicy
    feature_changes: tuple[FeatureChange, ...]
    letter_annotation_action: LetterAnnotationAction
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if not isinstance(self.record, DNARecord):
            raise ConfigurationError("RecordOperationResult record must be DNARecord.")
        if self.operation not in {
            "insert",
            "delete",
            "substitute",
            "mask",
            "trim",
            "reverse_complement",
            "rotate",
            "canonical_origin",
        }:
            raise ConfigurationError("RecordOperationResult operation is invalid.")
        if any(not isinstance(edit, Edit) for edit in self.edits):
            raise ConfigurationError("RecordOperationResult edits must contain Edit objects.")
        if self.feature_policy not in _FEATURE_POLICIES:
            raise ConfigurationError("RecordOperationResult feature_policy is invalid.")
        if any(not isinstance(change, FeatureChange) for change in self.feature_changes):
            raise ConfigurationError(
                "RecordOperationResult feature_changes must contain FeatureChange objects."
            )
        if self.letter_annotation_action not in _ANNOTATION_ACTIONS:
            raise ConfigurationError("letter_annotation_action is invalid.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("RecordOperationResult parameters must be FrozenDict.")

    @property
    def sequence(self) -> DNASequence:
        """Convenience alias for the transformed record sequence."""

        return self.record.sequence

    @property
    def deleted_feature_count(self) -> int:
        return sum(change.action == "deleted" for change in self.feature_changes)

    def to_dict(self) -> dict[str, Any]:
        """Return an audit payload without duplicating the full output sequence."""

        payload = {
            "record_id": self.record.id,
            "operation": self.operation,
            "edits": self.edits,
            "feature_policy": self.feature_policy,
            "feature_changes": self.feature_changes,
            "deleted_feature_count": self.deleted_feature_count,
            "letter_annotation_action": self.letter_annotation_action,
            "parameters": self.parameters,
            "output_symbol_length": self.record.sequence.symbol_length,
            "output_coordinate_span": self.record.sequence.coordinate_span,
            "output_topology": self.record.sequence.topology,
        }
        return cast(dict[str, Any], to_json_compatible(payload))


def _require_record(record: DNARecord) -> DNARecord:
    if not isinstance(record, DNARecord):
        raise ConfigurationError(
            "A DNARecord object is required.",
            code="INVALID_RECORD_OPERATION_ARGUMENT",
            context={"type": type(record).__name__},
        )
    return record


def _feature_policy(value: FeatureOverlapPolicy) -> FeatureOverlapPolicy:
    if not isinstance(value, str) or value not in _FEATURE_POLICIES:
        raise ConfigurationError(
            "Unknown feature overlap policy.",
            code="INVALID_FEATURE_OVERLAP_POLICY",
            context={"feature_policy": value},
            hint="Choose preserve, truncate, delete, split, or unresolved.",
        )
    return value


def _annotation_policy(value: LetterAnnotationPolicy) -> LetterAnnotationPolicy:
    if not isinstance(value, str) or value not in _ANNOTATION_POLICIES:
        raise ConfigurationError(
            "Unknown letter annotation policy.",
            code="INVALID_LETTER_ANNOTATION_POLICY",
            context={"letter_annotation_policy": value},
            hint="Choose error or drop.",
        )
    return value


def _copy_feature(feature: DNAFeature, location: Location) -> DNAFeature:
    return DNAFeature(
        feature.type,
        location,
        id=feature.id,
        strand=feature.strand,
        label=feature.label,
        score=feature.score,
        phase=feature.phase,
        qualifiers=feature.qualifiers,
        source=feature.source,
    )


def _location_parts(location: Location) -> tuple[Interval, ...]:
    if isinstance(location, Interval):
        return (location,)
    if isinstance(location, CompoundLocation):
        return location.parts
    return location.anchors


def _resolved_location(parts: Iterable[Interval]) -> Location | None:
    materialized = tuple(parts)
    if not materialized:
        return None
    return materialized[0] if len(materialized) == 1 else CompoundLocation(materialized)


def _edit_delta(edit: Edit) -> int:
    replacement_length = sum(
        len(part) if isinstance(part, str) else cast(int, part.length)
        for part in edit.replacement_parts
    )
    return replacement_length - (edit.end - edit.start)


def _replacement_length(edit: Edit) -> int:
    return edit.end - edit.start + _edit_delta(edit)


def _interval_affected(interval: Interval, edit: Edit) -> bool:
    if edit.start == edit.end:
        return _replacement_length(edit) > 0 and interval.start < edit.start < interval.end
    return interval.start < edit.end and edit.start < interval.end


def _map_start(position: int, edit: Edit) -> int:
    delta = _edit_delta(edit)
    if position < edit.start:
        return position
    if position >= edit.end:
        return position + delta
    return edit.start


def _map_end(position: int, edit: Edit) -> int:
    delta = _edit_delta(edit)
    if position <= edit.start:
        return position
    if position >= edit.end:
        return position + delta
    return edit.start + _replacement_length(edit)


def _preserved_interval(interval: Interval, edit: Edit) -> Interval:
    if edit.start == edit.end:
        delta = _edit_delta(edit)
        if len(interval) == 0:
            mapped = interval.start + delta if interval.start >= edit.start else interval.start
            return Interval(mapped, mapped)
        start = interval.start + delta if interval.start >= edit.start else interval.start
        end = interval.end + delta if interval.end > edit.start else interval.end
        return Interval(start, end)
    if len(interval) == 0:
        position = interval.start
        delta = _edit_delta(edit)
        if position <= edit.start:
            mapped = position
        elif position >= edit.end:
            mapped = position + delta
        else:
            mapped = edit.start
        return Interval(mapped, mapped)
    return Interval(_map_start(interval.start, edit), _map_end(interval.end, edit))


def _surviving_intervals(interval: Interval, edit: Edit) -> tuple[Interval, ...]:
    survivors: list[Interval] = []
    before_end = min(interval.end, edit.start)
    if interval.start < before_end:
        survivors.append(Interval(interval.start, before_end))
    after_start = max(interval.start, edit.end)
    if after_start < interval.end:
        delta = _edit_delta(edit)
        survivors.append(Interval(after_start + delta, interval.end + delta))
    return tuple(survivors)


def _sync_feature_for_edit(
    feature: DNAFeature,
    edit: Edit,
    policy: FeatureOverlapPolicy,
) -> tuple[DNAFeature | None, FeatureChangeAction, str | None]:
    location = feature.location
    parts = _location_parts(location)
    affected = any(_interval_affected(part, edit) for part in parts)
    if isinstance(location, UnresolvedLocation):
        mapped = tuple(_preserved_interval(anchor, edit) for anchor in parts)
        reason = f"{location.reason}; coordinates transformed across {edit.kind}"
        return _copy_feature(feature, UnresolvedLocation(reason, mapped)), "unresolved", reason

    if not affected:
        mapped_location = _resolved_location(_preserved_interval(part, edit) for part in parts)
        assert mapped_location is not None
        action: FeatureChangeAction = "preserved" if mapped_location == location else "shifted"
        return _copy_feature(feature, mapped_location), action, None

    if policy == "delete":
        return None, "deleted", f"feature overlaps {edit.kind}"

    survivors = tuple(survivor for part in parts for survivor in _surviving_intervals(part, edit))
    if policy == "unresolved":
        reason = f"feature overlaps {edit.kind}; exact biological extent is unresolved"
        return _copy_feature(feature, UnresolvedLocation(reason, survivors)), "unresolved", reason
    if policy == "truncate":
        if not survivors:
            return None, "deleted", f"feature is completely replaced by {edit.kind}"
        selected = min(survivors, key=lambda part: (-len(part), part.start, part.end))
        return _copy_feature(feature, selected), "truncated", "longest surviving segment retained"
    if policy == "split":
        mapped_location = _resolved_location(survivors)
        if mapped_location is None:
            return None, "deleted", f"feature is completely replaced by {edit.kind}"
        action = "split" if len(survivors) > 1 else "truncated"
        return _copy_feature(feature, mapped_location), action, "replacement coordinates excluded"

    mapped_location = _resolved_location(_preserved_interval(part, edit) for part in parts)
    assert mapped_location is not None
    action = "preserved" if mapped_location == location else "resized"
    return _copy_feature(feature, mapped_location), action, "replacement coordinates included"


def _sync_features_for_edits(
    features: tuple[DNAFeature, ...],
    edits: tuple[Edit, ...],
    policy: FeatureOverlapPolicy,
) -> tuple[tuple[DNAFeature, ...], tuple[FeatureChange, ...]]:
    output: list[DNAFeature] = []
    changes: list[FeatureChange] = []
    indexed_edits = tuple(
        sorted(enumerate(edits), key=lambda item: (item[1].start, item[1].end), reverse=True)
    )
    for feature_index, original in enumerate(features):
        current: DNAFeature | None = original
        actions: list[FeatureChangeAction] = []
        reasons: list[str] = []
        affected_indices: list[int] = []
        for edit_index, edit in indexed_edits:
            if current is None:
                break
            was_affected = any(
                _interval_affected(part, edit) for part in _location_parts(current.location)
            )
            current, action, reason = _sync_feature_for_edit(current, edit, policy)
            actions.append(action)
            if was_affected:
                affected_indices.append(edit_index)
            if reason is not None:
                reasons.append(reason)
        if current is not None:
            output.append(current)
        significant = next(
            (action for action in reversed(actions) if action not in {"preserved", "shifted"}),
            "shifted" if "shifted" in actions else "preserved",
        )
        changes.append(
            FeatureChange(
                feature_index,
                original.id,
                original.type,
                "deleted" if current is None else significant,
                original.location,
                None if current is None else current.location,
                policy,
                tuple(affected_indices),
                "; ".join(reasons) or None,
            )
        )
    return tuple(output), tuple(changes)


def _symbol_boundary(sequence: DNASequence, coordinate: int) -> int:
    coordinate_cursor = 0
    symbol_cursor = 0
    for part in sequence.parts:
        if isinstance(part, str):
            part_end = coordinate_cursor + len(part)
            if coordinate_cursor <= coordinate <= part_end:
                return symbol_cursor + coordinate - coordinate_cursor
            coordinate_cursor = part_end
            symbol_cursor += len(part)
        else:
            if part.length is None:
                raise UnknownLengthError(
                    "Letter annotations cannot be synchronized across an unknown-length Gap.",
                    code="UNRESOLVED_LETTER_ANNOTATION_COORDINATES",
                )
            part_end = coordinate_cursor + part.length
            if coordinate_cursor <= coordinate <= part_end:
                return symbol_cursor
            coordinate_cursor = part_end
    if coordinate == coordinate_cursor:
        return symbol_cursor
    raise ConfigurationError(
        "Letter annotation coordinate exceeds the sequence span.",
        code="LETTER_ANNOTATION_COORDINATE_OUT_OF_RANGE",
    )


def _replacement_annotation_values(
    record: DNARecord,
    replacement_length: int,
    replacement_annotations: Mapping[str, Iterable[int | float]] | None,
    policy: LetterAnnotationPolicy,
) -> tuple[Mapping[str, tuple[int | float, ...]], str]:
    if policy == "drop":
        if replacement_annotations:
            raise ConfigurationError(
                "replacement_annotations cannot be used when annotations are dropped.",
                code="UNEXPECTED_REPLACEMENT_ANNOTATIONS",
            )
        return {}, "dropped" if record.letter_annotations else "absent"
    if not record.letter_annotations:
        if replacement_annotations:
            raise ConfigurationError(
                "Cannot create partial letter annotations for an unannotated record.",
                code="PARTIAL_REPLACEMENT_ANNOTATIONS",
            )
        return {}, "absent"
    if replacement_length == 0 and replacement_annotations is None:
        return {}, "updated"
    if replacement_annotations is None:
        raise ConfigurationError(
            "Inserted or substituted symbols require replacement letter annotations.",
            code="REPLACEMENT_ANNOTATIONS_REQUIRED",
            context={"replacement_length": replacement_length},
            hint="Provide every existing annotation key or choose letter_annotation_policy='drop'.",
        )
    if set(replacement_annotations) != set(record.letter_annotations):
        raise ConfigurationError(
            "Replacement annotation keys must exactly match existing annotation keys.",
            code="REPLACEMENT_ANNOTATION_KEYS_MISMATCH",
            context={
                "expected": sorted(record.letter_annotations),
                "observed": sorted(replacement_annotations),
            },
        )
    resolved: dict[str, tuple[int | float, ...]] = {}
    for name, values in replacement_annotations.items():
        materialized = tuple(values)
        if len(materialized) != replacement_length:
            raise ConfigurationError(
                "Replacement annotation length must equal replacement symbol length.",
                code="REPLACEMENT_ANNOTATION_LENGTH_MISMATCH",
                context={
                    "name": name,
                    "annotation_length": len(materialized),
                    "replacement_length": replacement_length,
                },
            )
        resolved[name] = materialized
    return resolved, "updated"


def _edit_annotations(
    record: DNARecord,
    edit: Edit,
    replacement_annotations: Mapping[str, Iterable[int | float]] | None,
    policy: LetterAnnotationPolicy,
) -> tuple[Mapping[str, tuple[int | float, ...]], LetterAnnotationAction]:
    replacement_length = len(edit.replacement_symbols)
    replacements, action = _replacement_annotation_values(
        record, replacement_length, replacement_annotations, policy
    )
    if policy == "drop" or not record.letter_annotations:
        return {}, cast(LetterAnnotationAction, action)
    start = _symbol_boundary(record.sequence, edit.start)
    end = _symbol_boundary(record.sequence, edit.end)
    return (
        {
            name: (*values[:start], *replacements.get(name, ()), *values[end:])
            for name, values in record.letter_annotations.items()
        },
        cast(LetterAnnotationAction, action),
    )


def _result(
    record: DNARecord,
    sequence: DNASequence,
    operation: RecordOperation,
    edits: tuple[Edit, ...],
    feature_policy: FeatureOverlapPolicy,
    annotations: Mapping[str, Iterable[int | float]],
    annotation_action: LetterAnnotationAction,
    parameters: Mapping[str, object],
) -> RecordOperationResult:
    features, changes = _sync_features_for_edits(record.features, edits, feature_policy)
    transformed = DNARecord(
        sequence,
        record.id,
        description=record.description,
        features=features,
        metadata=record.metadata,
        letter_annotations=annotations,
    )
    return RecordOperationResult(
        transformed,
        operation,
        edits,
        feature_policy,
        changes,
        annotation_action,
        FrozenDict(
            {
                **parameters,
                "qualifiers": "preserved",
                "feature_phase": "preserved_not_recomputed",
            }
        ),
    )


def insert_record(
    record: DNARecord,
    position: int,
    fragment: DNASequence | str,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult:
    """Insert sequence and synchronize record annotations explicitly."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    annotation_policy = _annotation_policy(letter_annotation_policy)
    edited = insert(source.sequence, position, fragment)
    annotations, action = _edit_annotations(
        source, edited.edits[0], replacement_annotations, annotation_policy
    )
    return _result(
        source,
        edited.sequence,
        "insert",
        edited.edits,
        policy,
        annotations,
        action,
        {"position": position, "letter_annotation_policy": annotation_policy},
    )


def delete_record(
    record: DNARecord,
    start: int,
    end: int,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult:
    """Delete sequence and map or resolve every affected feature."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    annotation_policy = _annotation_policy(letter_annotation_policy)
    edited = delete(source.sequence, start, end)
    annotations, action = _edit_annotations(source, edited.edits[0], None, annotation_policy)
    return _result(
        source,
        edited.sequence,
        "delete",
        edited.edits,
        policy,
        annotations,
        action,
        {"start": start, "end": end, "letter_annotation_policy": annotation_policy},
    )


def substitute_record(
    record: DNARecord,
    start: int,
    end: int,
    fragment: DNASequence | str,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult:
    """Substitute sequence and synchronize features and letter annotations."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    annotation_policy = _annotation_policy(letter_annotation_policy)
    edited = substitute(source.sequence, start, end, fragment)
    annotations, action = _edit_annotations(
        source, edited.edits[0], replacement_annotations, annotation_policy
    )
    return _result(
        source,
        edited.sequence,
        "substitute",
        edited.edits,
        policy,
        annotations,
        action,
        {"start": start, "end": end, "letter_annotation_policy": annotation_policy},
    )


def mask_record(
    record: DNARecord,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult:
    """Mask symbols while preserving aligned letter values unless explicitly dropped."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    annotation_policy = _annotation_policy(letter_annotation_policy)
    edited = mask(source.sequence, intervals, symbol=symbol)
    annotations: Mapping[str, Iterable[int | float]] = (
        {} if annotation_policy == "drop" else source.letter_annotations
    )
    action: LetterAnnotationAction = (
        "dropped"
        if annotation_policy == "drop" and source.letter_annotations
        else "preserved"
        if source.letter_annotations
        else "absent"
    )
    return _result(
        source,
        edited.sequence,
        "mask",
        edited.edits,
        policy,
        annotations,
        action,
        {"symbol": symbol, "letter_annotation_policy": annotation_policy},
    )


def trim_record(
    record: DNARecord,
    *,
    left: int = 0,
    right: int = 0,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult:
    """Trim record ends and synchronize features in original coordinates."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    annotation_policy = _annotation_policy(letter_annotation_policy)
    edited = trim(source.sequence, left=left, right=right)
    if annotation_policy == "drop":
        annotations: Mapping[str, Iterable[int | float]] = {}
        action: LetterAnnotationAction = "dropped" if source.letter_annotations else "absent"
    elif source.letter_annotations:
        span = source.sequence.coordinate_span
        assert span is not None
        symbol_start = _symbol_boundary(source.sequence, left)
        symbol_end = _symbol_boundary(source.sequence, span - right)
        annotations = {
            name: values[symbol_start:symbol_end]
            for name, values in source.letter_annotations.items()
        }
        action = "trimmed" if left or right else "preserved"
    else:
        annotations = {}
        action = "absent"
    return _result(
        source,
        edited.sequence,
        "trim",
        edited.edits,
        policy,
        annotations,
        action,
        {
            "left": left,
            "right": right,
            "letter_annotation_policy": annotation_policy,
            "edit_coordinate_space": "original_record",
        },
    )


def _reverse_strand(strand: Strand) -> Strand:
    if strand is Strand.FORWARD:
        return Strand.REVERSE
    if strand is Strand.REVERSE:
        return Strand.FORWARD
    return strand


def _copy_feature_with_strand(
    feature: DNAFeature, location: Location, strand: Strand
) -> DNAFeature:
    return DNAFeature(
        feature.type,
        location,
        id=feature.id,
        strand=strand,
        label=feature.label,
        score=feature.score,
        phase=feature.phase,
        qualifiers=feature.qualifiers,
        source=feature.source,
    )


def reverse_complement_record(
    record: DNARecord,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
) -> RecordOperationResult:
    """Reverse-complement sequence, locations, strands, and letter arrays."""

    source = _require_record(record)
    policy = _feature_policy(feature_policy)
    span = source.sequence.coordinate_span
    output_features: list[DNAFeature] = []
    changes: list[FeatureChange] = []
    for index, feature in enumerate(source.features):
        if span is None:
            if policy == "delete":
                changes.append(
                    FeatureChange(
                        index,
                        feature.id,
                        feature.type,
                        "deleted",
                        feature.location,
                        None,
                        policy,
                        reason="unknown sequence span prevents reverse coordinate mapping",
                    )
                )
                continue
            if policy != "unresolved":
                raise UnknownLengthError(
                    "Feature coordinates cannot be reverse-complemented across an unknown Gap.",
                    code="UNRESOLVED_REVERSE_FEATURE_COORDINATES",
                    context={"feature_index": index},
                    hint="Choose feature_policy='unresolved' or 'delete'.",
                )
            location: Location = UnresolvedLocation("reverse-complement coordinate span is unknown")
            action: FeatureChangeAction = "unresolved"
        else:
            location = reverse_strand_location(feature.location, sequence_length=span)
            action = "reverse_complemented"
        mapped = _copy_feature_with_strand(feature, location, _reverse_strand(feature.strand))
        output_features.append(mapped)
        changes.append(
            FeatureChange(
                index,
                feature.id,
                feature.type,
                action,
                feature.location,
                location,
                policy,
            )
        )
    transformed = DNARecord(
        source.sequence.reverse_complement(),
        source.id,
        description=source.description,
        features=output_features,
        metadata=source.metadata,
        letter_annotations={
            name: tuple(reversed(values)) for name, values in source.letter_annotations.items()
        },
    )
    return RecordOperationResult(
        transformed,
        "reverse_complement",
        (),
        policy,
        tuple(changes),
        "reversed" if source.letter_annotations else "absent",
        FrozenDict(
            {
                "strand_mapping": "forward<->reverse; both/unknown unchanged",
                "qualifiers": "preserved",
                "feature_phase": "preserved_not_recomputed",
            }
        ),
    )


def _rotate_interval(interval: Interval, offset: int, span: int) -> tuple[Interval, ...]:
    length = len(interval)
    start = (interval.start - offset) % span
    if length == 0:
        return (Interval(start, start),)
    if length == span:
        return (Interval(0, span),)
    end = start + length
    if end <= span:
        return (Interval(start, end),)
    return (Interval(start, span), Interval(0, end - span))


def _rotate_feature(
    feature: DNAFeature,
    offset: int,
    span: int,
    policy: FeatureOverlapPolicy,
) -> tuple[DNAFeature | None, FeatureChangeAction, str | None]:
    if isinstance(feature.location, UnresolvedLocation):
        anchors = tuple(
            mapped
            for anchor in feature.location.anchors
            for mapped in _rotate_interval(anchor, offset, span)
        )
        unresolved_location: Location = UnresolvedLocation(feature.location.reason, anchors)
        return _copy_feature(feature, unresolved_location), "unresolved", None
    source_parts = _location_parts(feature.location)
    mapped_groups = tuple(_rotate_interval(part, offset, span) for part in source_parts)
    mapped = tuple(transformed for group in mapped_groups for transformed in group)
    wrapped = any(len(group) > 1 for group in mapped_groups)
    if wrapped and policy == "delete":
        return None, "deleted", "feature crosses the new circular origin"
    if wrapped and policy == "unresolved":
        reason = "feature crosses the new circular origin"
        return _copy_feature(feature, UnresolvedLocation(reason, mapped)), "unresolved", reason
    if wrapped and policy == "truncate":
        selected = min(mapped, key=lambda part: (-len(part), part.start, part.end))
        return _copy_feature(feature, selected), "truncated", "longest origin-side retained"
    resolved_location = _resolved_location(mapped)
    assert resolved_location is not None
    action: FeatureChangeAction = (
        "split" if wrapped and policy == "split" else "preserved" if wrapped else "rotated"
    )
    return (
        _copy_feature(feature, resolved_location),
        action,
        "origin-spanning feature represented as CompoundLocation" if wrapped else None,
    )


def _rotation_symbol_offset(sequence: DNASequence, coordinate: int) -> int:
    return _symbol_boundary(sequence, coordinate)


def _rotate_record(
    source: DNARecord,
    *,
    offset: int,
    operation: Literal["rotate", "canonical_origin"],
    requested_offset: int | None,
    feature_policy: FeatureOverlapPolicy,
) -> RecordOperationResult:
    rotated = rotate(source.sequence, offset)
    policy = _feature_policy(feature_policy)
    features: list[DNAFeature] = []
    changes: list[FeatureChange] = []
    for index, feature in enumerate(source.features):
        mapped, action, reason = _rotate_feature(
            feature, rotated.effective_offset, rotated.sequence_span, policy
        )
        if mapped is not None:
            features.append(mapped)
        changes.append(
            FeatureChange(
                index,
                feature.id,
                feature.type,
                action,
                feature.location,
                None if mapped is None else mapped.location,
                policy,
                reason=reason,
            )
        )
    symbol_offset = _rotation_symbol_offset(source.sequence, rotated.effective_offset)
    annotations = {
        name: (*values[symbol_offset:], *values[:symbol_offset])
        for name, values in source.letter_annotations.items()
    }
    transformed = DNARecord(
        rotated.sequence,
        source.id,
        description=source.description,
        features=features,
        metadata=source.metadata,
        letter_annotations=annotations,
    )
    origin_rule = (
        "lexicographically_minimal_forward_rotation_then_smallest_offset"
        if operation == "canonical_origin"
        else rotated.rule
    )
    operation_parameters: dict[str, object] = {}
    if operation == "canonical_origin":
        operation_parameters = {
            "comparison": "normalized_symbol_code_point_order",
            "orientation": "forward_only",
            "tie_break": "smallest_offset",
            "algorithm": "booth",
        }
    return RecordOperationResult(
        transformed,
        operation,
        (),
        policy,
        tuple(changes),
        "rotated" if annotations else "absent",
        FrozenDict(
            {
                "requested_offset": requested_offset,
                "effective_offset": rotated.effective_offset,
                "sequence_span": rotated.sequence_span,
                "origin_rule": origin_rule,
                **operation_parameters,
                "qualifiers": "preserved",
                "feature_phase": "preserved_not_recomputed",
            }
        ),
    )


def rotate_record(
    record: DNARecord,
    offset: int,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> RecordOperationResult:
    """Rotate a circular record and map origin-spanning feature locations."""

    return _rotate_record(
        _require_record(record),
        offset=offset,
        operation="rotate",
        requested_offset=offset,
        feature_policy=feature_policy,
    )


def canonical_origin_record(
    record: DNARecord,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> RecordOperationResult:
    """Apply the canonical forward lexicographic origin to a record."""

    source = _require_record(record)
    canonical = canonical_origin(source.sequence)
    return _rotate_record(
        source,
        offset=canonical.effective_offset,
        operation="canonical_origin",
        requested_offset=None,
        feature_policy=feature_policy,
    )


__all__ = [
    "FeatureChange",
    "FeatureChangeAction",
    "FeatureOverlapPolicy",
    "LetterAnnotationAction",
    "LetterAnnotationPolicy",
    "RecordOperationResult",
    "canonical_origin_record",
    "delete_record",
    "insert_record",
    "mask_record",
    "reverse_complement_record",
    "rotate_record",
    "substitute_record",
    "trim_record",
]
