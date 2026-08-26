"""Immutable metadata add, merge, filter, projection, and validation APIs."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Literal, cast, overload

from dnakit.core._json import JSONValue, to_json_compatible
from dnakit.core.collection import DNASet
from dnakit.core.facade import DNA
from dnakit.core.record import DNARecord
from dnakit.exceptions import ConfigurationError, InputFormatError

ConflictPolicy = Literal["error", "replace", "keep"]
MissingPolicy = Literal["error", "keep", "drop"]
ExtraPolicy = Literal["error", "ignore"]
MetadataType = Literal["string", "integer", "number", "boolean", "null", "array", "object"]
MetadataPredicate = Callable[[JSONValue], bool]
MetadataRecordInput = DNA | DNARecord
MetadataDatasetInput = DNA | DNASet


def _copy_record(record: DNARecord, metadata: Mapping[str, object]) -> DNARecord:
    return DNARecord(
        record.sequence,
        record.id,
        description=record.description,
        features=record.features,
        metadata=metadata,
        letter_annotations=record.letter_annotations,
    )


def _copy_set(dataset: DNASet, records: Iterable[DNARecord]) -> DNASet:
    return DNASet(
        records,
        name=dataset.name,
        source=dataset.source,
        version=dataset.version,
        metadata=dataset.metadata,
    )


def _record_value(value: MetadataRecordInput) -> DNARecord:
    if isinstance(value, DNA):
        return value.record
    if isinstance(value, DNARecord):
        return value
    raise TypeError("record must be DNA or DNARecord.")


def _record_result(source: MetadataRecordInput, record: DNARecord) -> DNA | DNARecord:
    return source._derive(((0, record),)) if isinstance(source, DNA) else record


def _dataset_value(value: MetadataDatasetInput) -> DNASet:
    if isinstance(value, DNA):
        return value.dataset
    if isinstance(value, DNASet):
        return value
    raise TypeError("dataset must be DNA or DNASet.")


def _dataset_result(
    source: MetadataDatasetInput,
    indexed_records: Iterable[tuple[int, DNARecord]],
) -> DNA | DNASet:
    selected = tuple(indexed_records)
    if isinstance(source, DNA):
        return source._derive(selected)
    return _copy_set(source, (record for _, record in selected))


@overload
def with_metadata(
    record: DNA,
    values: Mapping[str, object],
    *,
    conflict: ConflictPolicy = "error",
) -> DNA: ...


@overload
def with_metadata(
    record: DNARecord,
    values: Mapping[str, object],
    *,
    conflict: ConflictPolicy = "error",
) -> DNARecord: ...


def with_metadata(
    record: MetadataRecordInput,
    values: Mapping[str, object],
    *,
    conflict: ConflictPolicy = "error",
) -> DNA | DNARecord:
    """Return a new record after a top-level immutable metadata merge."""

    source_record = _record_value(record)
    if not isinstance(values, Mapping):
        raise TypeError("values must be a mapping.")
    if conflict not in {"error", "replace", "keep"}:
        raise ConfigurationError(
            "Unknown metadata conflict policy.", code="INVALID_METADATA_POLICY"
        )
    result: dict[str, object] = dict(source_record.metadata)
    for key, value in values.items():
        if not isinstance(key, str) or not key:
            raise ConfigurationError(
                "Metadata keys must be non-empty strings.", code="INVALID_METADATA_KEY"
            )
        if key in result:
            if conflict == "error" and result[key] != value:
                raise InputFormatError(
                    "Metadata merge encountered a conflicting value.",
                    code="METADATA_CONFLICT",
                    context={"record_id": source_record.id, "field": key},
                )
            if conflict == "keep":
                continue
        result[key] = value
    return _record_result(record, _copy_record(source_record, result))


@overload
def merge_metadata(
    dataset: DNA,
    metadata_by_id: Mapping[str, Mapping[str, object]],
    *,
    conflict: ConflictPolicy = "error",
    missing: MissingPolicy = "error",
    extra: ExtraPolicy = "error",
) -> DNA: ...


@overload
def merge_metadata(
    dataset: DNASet,
    metadata_by_id: Mapping[str, Mapping[str, object]],
    *,
    conflict: ConflictPolicy = "error",
    missing: MissingPolicy = "error",
    extra: ExtraPolicy = "error",
) -> DNASet: ...


def merge_metadata(
    dataset: MetadataDatasetInput,
    metadata_by_id: Mapping[str, Mapping[str, object]],
    *,
    conflict: ConflictPolicy = "error",
    missing: MissingPolicy = "error",
    extra: ExtraPolicy = "error",
) -> DNA | DNASet:
    """Left-join metadata by exact record ID with explicit missing/extra policies."""

    source_dataset = _dataset_value(dataset)
    if not isinstance(metadata_by_id, Mapping):
        raise TypeError("metadata_by_id must be a mapping.")
    if missing not in {"error", "keep", "drop"} or extra not in {"error", "ignore"}:
        raise ConfigurationError("Unknown metadata join policy.", code="INVALID_METADATA_POLICY")
    if len(set(source_dataset.ids)) != len(source_dataset.ids):
        raise InputFormatError(
            "Metadata joins require unique DNASet record IDs.",
            code="DUPLICATE_METADATA_JOIN_ID",
        )
    invalid = tuple(
        key
        for key, value in metadata_by_id.items()
        if not isinstance(key, str) or not key or not isinstance(value, Mapping)
    )
    if invalid:
        raise ConfigurationError(
            "Metadata join keys must be non-empty strings and values must be mappings.",
            code="INVALID_METADATA_TABLE",
            context={"invalid_keys": invalid[:10]},
        )
    input_ids = set(source_dataset.ids)
    extra_ids = tuple(sorted(set(metadata_by_id) - input_ids))
    if extra_ids and extra == "error":
        raise InputFormatError(
            "Metadata table contains IDs absent from the DNASet.",
            code="EXTRA_METADATA_IDS",
            context={"count": len(extra_ids), "sample": extra_ids[:10]},
        )
    output: list[tuple[int, DNARecord]] = []
    for index, record in enumerate(source_dataset):
        values = metadata_by_id.get(record.id)
        if values is None:
            if missing == "error":
                raise InputFormatError(
                    "A DNASet record has no matching metadata row.",
                    code="MISSING_METADATA_ID",
                    context={"record_id": record.id},
                )
            if missing == "keep":
                output.append((index, record))
            continue
        output.append((index, with_metadata(record, values, conflict=conflict)))
    return _dataset_result(dataset, output)


@overload
def filter_by_metadata(
    dataset: DNA,
    criteria: Mapping[str, object | MetadataPredicate],
    *,
    require_all: bool = True,
) -> DNA: ...


@overload
def filter_by_metadata(
    dataset: DNASet,
    criteria: Mapping[str, object | MetadataPredicate],
    *,
    require_all: bool = True,
) -> DNASet: ...


def filter_by_metadata(
    dataset: MetadataDatasetInput,
    criteria: Mapping[str, object | MetadataPredicate],
    *,
    require_all: bool = True,
) -> DNA | DNASet:
    """Filter records by exact values or side-effect-free caller predicates."""

    source_dataset = _dataset_value(dataset)
    if not isinstance(criteria, Mapping) or not criteria:
        raise ConfigurationError(
            "criteria must be a non-empty mapping.", code="INVALID_METADATA_FILTER"
        )
    if any(not isinstance(field, str) or not field for field in criteria):
        raise ConfigurationError(
            "Metadata filter fields must be non-empty strings.",
            code="INVALID_METADATA_FILTER",
        )
    if not isinstance(require_all, bool):
        raise TypeError("require_all must be a boolean.")

    def matches(record: DNARecord) -> bool:
        decisions: list[bool] = []
        for field, expected in criteria.items():
            if field not in record.metadata:
                decisions.append(False)
                continue
            actual = record.metadata[field]
            decisions.append(bool(expected(actual)) if callable(expected) else actual == expected)
        return all(decisions) if require_all else any(decisions)

    return _dataset_result(
        dataset,
        ((index, record) for index, record in enumerate(source_dataset) if matches(record)),
    )


@overload
def select_metadata(
    record: DNA,
    fields: Iterable[str],
    *,
    missing: Literal["error", "ignore"] = "error",
) -> DNA: ...


@overload
def select_metadata(
    record: DNARecord,
    fields: Iterable[str],
    *,
    missing: Literal["error", "ignore"] = "error",
) -> DNARecord: ...


def select_metadata(
    record: MetadataRecordInput,
    fields: Iterable[str],
    *,
    missing: Literal["error", "ignore"] = "error",
) -> DNA | DNARecord:
    """Return a record containing only selected metadata fields."""

    source_record = _record_value(record)
    if missing not in {"error", "ignore"}:
        raise ConfigurationError(
            "Unknown projection missing policy.", code="INVALID_METADATA_POLICY"
        )
    field_tuple = tuple(fields)
    if any(not isinstance(field, str) or not field for field in field_tuple):
        raise ConfigurationError(
            "Metadata projection fields must be non-empty strings.", code="INVALID_METADATA_KEY"
        )
    absent = tuple(field for field in field_tuple if field not in source_record.metadata)
    if absent and missing == "error":
        raise InputFormatError(
            "Metadata projection references missing fields.",
            code="MISSING_METADATA_FIELD",
            context={"record_id": source_record.id, "fields": absent},
        )
    return _record_result(
        record,
        _copy_record(
            source_record,
            {
                field: source_record.metadata[field]
                for field in field_tuple
                if field in source_record.metadata
            },
        ),
    )


def _metadata_type(value: JSONValue) -> MetadataType:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "number"
    if isinstance(value, str):
        return "string"
    if isinstance(value, Mapping):
        return "object"
    return "array"


@dataclass(frozen=True, slots=True)
class MetadataValidationIssue:
    record_id: str
    field: str
    code: Literal["missing", "type_mismatch"]
    expected_type: MetadataType
    actual_type: MetadataType | None

    def __post_init__(self) -> None:
        allowed = {"string", "integer", "number", "boolean", "null", "array", "object"}
        if not isinstance(self.record_id, str) or not self.record_id:
            raise ConfigurationError(
                "Metadata issue record_id must be non-empty.",
                code="INVALID_METADATA_REPORT",
            )
        if not isinstance(self.field, str) or not self.field:
            raise ConfigurationError(
                "Metadata issue field must be non-empty.", code="INVALID_METADATA_REPORT"
            )
        if (
            self.code not in {"missing", "type_mismatch"}
            or self.expected_type not in allowed
            or (self.actual_type is not None and self.actual_type not in allowed)
        ):
            raise ConfigurationError(
                "Metadata issue contains an invalid code or type.",
                code="INVALID_METADATA_REPORT",
            )
        if (self.code == "missing") != (self.actual_type is None):
            raise ConfigurationError(
                "Metadata issue code and actual_type are inconsistent.",
                code="INVALID_METADATA_REPORT",
            )


@dataclass(frozen=True, slots=True)
class MetadataValidationReport:
    record_count: int
    issues: tuple[MetadataValidationIssue, ...]

    def __post_init__(self) -> None:
        if (
            isinstance(self.record_count, bool)
            or not isinstance(self.record_count, int)
            or self.record_count < 0
            or not isinstance(self.issues, tuple)
            or any(not isinstance(issue, MetadataValidationIssue) for issue in self.issues)
        ):
            raise ConfigurationError(
                "Metadata validation report fields are invalid.",
                code="INVALID_METADATA_REPORT",
            )

    @property
    def valid(self) -> bool:
        return not self.issues

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def validate_metadata(
    dataset: MetadataDatasetInput,
    schema: Mapping[str, MetadataType],
    *,
    required: Iterable[str] = (),
) -> MetadataValidationReport:
    """Validate top-level metadata presence and stable JSON scalar/container types."""

    source_dataset = _dataset_value(dataset)
    allowed = {"string", "integer", "number", "boolean", "null", "array", "object"}
    if not isinstance(schema, Mapping) or any(
        not isinstance(field, str) or not field or not isinstance(kind, str) or kind not in allowed
        for field, kind in schema.items()
    ):
        raise ConfigurationError(
            "Metadata schema contains an invalid field or type.", code="INVALID_METADATA_SCHEMA"
        )
    required_tuple = tuple(required)
    if any(not isinstance(field, str) or not field for field in required_tuple):
        raise ConfigurationError(
            "Required metadata fields must be non-empty strings.",
            code="INVALID_METADATA_SCHEMA",
        )
    if any(field not in schema for field in required_tuple):
        raise ConfigurationError(
            "Every required metadata field must appear in schema.", code="INVALID_METADATA_SCHEMA"
        )
    issues: list[MetadataValidationIssue] = []
    for record in source_dataset:
        for field, expected in schema.items():
            value = record.metadata.get(field)
            if field not in record.metadata:
                if field in required_tuple:
                    issues.append(
                        MetadataValidationIssue(record.id, field, "missing", expected, None)
                    )
                continue
            actual = _metadata_type(value)
            compatible = actual == expected or (expected == "number" and actual == "integer")
            if not compatible:
                issues.append(
                    MetadataValidationIssue(record.id, field, "type_mismatch", expected, actual)
                )
    return MetadataValidationReport(len(source_dataset), tuple(issues))


__all__ = [
    "ConflictPolicy",
    "ExtraPolicy",
    "MetadataPredicate",
    "MetadataType",
    "MetadataValidationIssue",
    "MetadataValidationReport",
    "MissingPolicy",
    "filter_by_metadata",
    "merge_metadata",
    "select_metadata",
    "validate_metadata",
    "with_metadata",
]
