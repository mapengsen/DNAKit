from __future__ import annotations

import json

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.exceptions import BackendUnavailableError, InputFormatError
from dnakit.io import (
    filter_by_metadata,
    merge_metadata,
    parquet_backend_status,
    require_parquet_backend,
    select_metadata,
    validate_metadata,
    with_metadata,
)


def _dataset() -> DNASet:
    return DNASet(
        (
            DNARecord(DNASequence("A"), "a", metadata={"species": "human", "batch": 1}),
            DNARecord(DNASequence("C"), "b", metadata={"species": "mouse"}),
        ),
        name="examples",
        metadata={"owner": "test"},
    )


def test_metadata_add_merge_filter_and_projection_are_immutable() -> None:
    dataset = _dataset()
    updated = with_metadata(dataset[0], {"batch": 2, "family": "f1"}, conflict="replace")

    assert dataset[0].metadata["batch"] == 1
    assert updated.metadata["batch"] == 2
    merged = merge_metadata(
        dataset,
        {"a": {"family": "f1"}, "b": {"family": "f2"}},
    )
    assert merged.name == dataset.name
    assert merged.metadata == dataset.metadata
    assert filter_by_metadata(merged, {"species": "human"}).ids == ("a",)
    assert filter_by_metadata(merged, {"family": lambda value: value == "f2"}).ids == ("b",)
    projected = select_metadata(merged[0], ["family"])
    assert dict(projected.metadata) == {"family": "f1"}


def test_metadata_join_policies_are_explicit() -> None:
    dataset = _dataset()

    with pytest.raises(InputFormatError) as conflict:
        with_metadata(dataset[0], {"batch": 2})
    assert conflict.value.code == "METADATA_CONFLICT"

    with pytest.raises(InputFormatError) as missing:
        merge_metadata(dataset, {"a": {}}, extra="ignore")
    assert missing.value.code == "MISSING_METADATA_ID"

    with pytest.raises(InputFormatError) as extra:
        merge_metadata(dataset, {"a": {}, "b": {}, "c": {}})
    assert extra.value.code == "EXTRA_METADATA_IDS"

    dropped = merge_metadata(dataset, {"a": {}}, missing="drop", extra="ignore")
    assert dropped.ids == ("a",)


def test_metadata_validation_report_is_json_compatible() -> None:
    report = validate_metadata(
        _dataset(),
        {"species": "string", "batch": "integer"},
        required=("species", "batch"),
    )

    assert not report.valid
    assert report.issues[0].record_id == "b"
    assert report.issues[0].code == "missing"
    assert json.loads(json.dumps(report.to_dict()))["record_count"] == 2


def test_parquet_probe_never_imports_or_installs_backend() -> None:
    status = parquet_backend_status()

    assert status.implementation == "conditional-adapter"
    assert json.loads(json.dumps(status.to_dict()))["distribution"] == "pyarrow"
    if status.available:
        assert require_parquet_backend().version == status.version
    else:
        with pytest.raises(BackendUnavailableError) as exc_info:
            require_parquet_backend()
        assert getattr(exc_info.value, "code", None) == "PARQUET_BACKEND_UNAVAILABLE"
