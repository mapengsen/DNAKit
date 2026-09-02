from __future__ import annotations

import csv
import json
import platform
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast, overload

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError
from dnakit.io import (
    TableSchema,
    export_result,
    export_table,
    parquet_backend_status,
    read_table,
)
from dnakit.io.tables import TableFormat


def test_csv_tsv_and_json_table_export_are_schema_stable(tmp_path: Path) -> None:
    rows = ({"id": "a", "score": 1.5}, {"id": "b", "score": None})
    schema = TableSchema(("id", "score"), schema_version="example.v1")

    csv_result = export_table(rows, tmp_path / "table.csv", format="csv", schema=schema)
    tsv_result = export_table(rows, tmp_path / "table.tsv", format="tsv", schema=schema)
    json_result = export_table(rows, tmp_path / "table.json", format="json", schema=schema)

    assert csv_result.row_count == tsv_result.row_count == json_result.row_count == 2
    assert csv_result.parameters["backend"] == {
        "name": "python-stdlib",
        "version": platform.python_version(),
        "implementation": "native",
        "license": "PSF-2.0",
    }
    with (tmp_path / "table.csv").open(encoding="utf-8", newline="") as handle:
        assert next(iter(csv.DictReader(handle))) == {"id": "a", "score": "1.5"}
    assert (tmp_path / "table.tsv").read_text(encoding="utf-8").startswith("id\tscore\n")
    payload = json.loads((tmp_path / "table.json").read_text(encoding="utf-8"))
    assert payload == {
        "columns": ["id", "score"],
        "schema_version": "example.v1",
        "rows": [{"id": "a", "score": 1.5}, {"id": "b", "score": None}],
    }


def test_table_export_accepts_sequence_rows_and_is_atomic_on_failure(tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    path.write_text("keep", encoding="utf-8")

    with pytest.raises(InputFormatError) as exc_info:
        export_table(
            (("a", 1), ("b", {"nested": True})),
            path,
            format="csv",
            schema=TableSchema(("id", "value")),
            overwrite=True,
        )

    assert exc_info.value.code == "INVALID_TABLE_CELL"
    assert path.read_text(encoding="utf-8") == "keep"


class _InfiniteColumns:
    def __init__(self) -> None:
        self.consumed = 0

    def __iter__(self) -> Iterator[str]:
        while True:
            self.consumed += 1
            yield f"column-{self.consumed}"


def test_table_schema_enforces_a_hard_bound_without_exhausting_iterables() -> None:
    columns = _InfiniteColumns()

    with pytest.raises(ConfigurationError) as exc_info:
        TableSchema(cast(Any, columns))

    assert exc_info.value.code == "INVALID_TABLE_SCHEMA"
    assert exc_info.value.context == {"max_columns": 10_000}
    assert columns.consumed == 10_001


class _InfiniteMapping(Mapping[str, object]):
    def __init__(self) -> None:
        self.consumed = 0

    def __getitem__(self, key: str) -> object:
        return key

    def __iter__(self) -> Iterator[str]:
        while True:
            self.consumed += 1
            yield f"key-{self.consumed}"

    def __len__(self) -> int:
        return 1


class _InfiniteSequence(Sequence[object]):
    def __init__(self) -> None:
        self.consumed = 0

    @overload
    def __getitem__(self, index: int) -> object: ...

    @overload
    def __getitem__(self, index: slice) -> Sequence[object]: ...

    def __getitem__(self, index: int | slice) -> object | Sequence[object]:
        return () if isinstance(index, slice) else index

    def __len__(self) -> int:
        return 1

    def __iter__(self) -> Iterator[object]:
        while True:
            self.consumed += 1
            yield self.consumed


@pytest.mark.parametrize("row", [_InfiniteMapping(), _InfiniteSequence()])
def test_table_rows_consume_at_most_schema_width_plus_one(
    tmp_path: Path, row: _InfiniteMapping | _InfiniteSequence
) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        export_table(
            (row,),
            tmp_path / "unbounded.json",
            format="json",
            schema=TableSchema(("id", "value")),
        )

    assert exc_info.value.code == "TABLE_SCHEMA_MISMATCH"
    assert row.consumed == 3
    assert not (tmp_path / "unbounded.json").exists()


def test_table_export_enforces_schema_rows_and_cell_limits(tmp_path: Path) -> None:
    schema = TableSchema(("id",))

    with pytest.raises(InputFormatError) as schema_error:
        export_table(({"other": "x"},), tmp_path / "schema.json", format="json", schema=schema)
    assert schema_error.value.code == "TABLE_SCHEMA_MISMATCH"

    with pytest.raises(InputFormatError) as row_error:
        export_table(
            ({"id": "a"}, {"id": "b"}),
            tmp_path / "rows.json",
            format="json",
            schema=schema,
            max_rows=1,
        )
    assert row_error.value.code == "TABLE_ROW_LIMIT_EXCEEDED"

    with pytest.raises(InputFormatError) as cell_error:
        export_table(
            ({"id": "long"},),
            tmp_path / "cell.json",
            format="json",
            schema=schema,
            max_cell_characters=3,
        )
    assert cell_error.value.code == "TABLE_CELL_LIMIT_EXCEEDED"

    with pytest.raises(ConfigurationError) as column_error:
        export_table(
            ({"id": "a"},),
            tmp_path / "columns.json",
            format="json",
            schema=schema,
            max_columns=0,
        )
    assert column_error.value.code == "INVALID_TABLE_LIMIT"

    with pytest.raises(InputFormatError) as column_count_error:
        export_table(
            ({"id": "a", "value": 1},),
            tmp_path / "too-many-columns.json",
            format="json",
            schema=TableSchema(("id", "value")),
            max_columns=1,
        )
    assert column_count_error.value.code == "TABLE_COLUMN_LIMIT_EXCEEDED"


@pytest.mark.parametrize("format", ["csv", "tsv", "json"])
def test_text_table_output_limit_is_atomic_and_preserves_existing_target(
    tmp_path: Path, format: str
) -> None:
    target = tmp_path / f"limited.{format}"
    target.write_text("keep", encoding="utf-8")

    with pytest.raises(InputFormatError) as exc_info:
        export_table(
            ({"id": "a" * 100},),
            target,
            format=cast(TableFormat, format),
            schema=TableSchema(("id",)),
            overwrite=True,
            max_output_bytes=5,
        )

    assert exc_info.value.code == "TABLE_OUTPUT_LIMIT_EXCEEDED"
    assert target.read_text(encoding="utf-8") == "keep"
    assert list(tmp_path.glob(f".{target.name}.*.tmp")) == []


@dataclass(frozen=True)
class ExampleResult:
    def to_dict(self) -> dict[str, Any]:
        return {"score": 1.0, "details": {"method": "native", "values": [1, 2]}}


def test_arbitrary_result_object_nested_values_export_as_canonical_json(tmp_path: Path) -> None:
    result = export_result(ExampleResult(), tmp_path / "result.csv", format="csv")

    assert result.row_count == 1
    with (tmp_path / "result.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["score"] == "1.0"
    assert row["details"] == '{"method":"native","values":[1,2]}'


def test_parquet_export_is_conditional_without_installation(tmp_path: Path) -> None:
    if parquet_backend_status().available:
        pytest.skip(
            "PyArrow is available; local optional-backend integration is environment-specific."
        )
    with pytest.raises(BackendUnavailableError) as exc_info:
        export_table(
            ({"id": "a"},),
            tmp_path / "table.parquet",
            format="parquet",
            schema=TableSchema(("id",)),
        )

    assert exc_info.value.code == "PARQUET_BACKEND_UNAVAILABLE"
    assert not (tmp_path / "table.parquet").exists()


def test_parquet_export_round_trips_when_optional_backend_is_available(tmp_path: Path) -> None:
    status = parquet_backend_status()
    if not status.available:
        pytest.skip("optional pyarrow backend is not installed")
    import pyarrow.parquet as pq  # type: ignore[import-untyped]

    target = tmp_path / "table.parquet"
    result = export_table(
        [{"id": "one", "score": 1.25}, {"id": "two", "score": None}],
        target,
        format="parquet",
        schema=TableSchema(("id", "score"), "dnakit.parquet-test.v1"),
    )

    table = pq.read_table(target)
    assert result.row_count == 2
    assert result.byte_count == target.stat().st_size
    assert table.to_pylist() == [
        {"id": "one", "score": 1.25},
        {"id": "two", "score": None},
    ]
    assert table.schema.metadata[b"dnakit_schema_version"] == b"dnakit.parquet-test.v1"


@dataclass(frozen=True)
class InvalidNestedResult:
    def to_dict(self) -> dict[str, Any]:
        return {"details": {"invalid": float("nan")}}


def test_result_export_rejects_nonstandard_nested_nan(tmp_path: Path) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        export_result(InvalidNestedResult(), tmp_path / "bad.json")

    assert exc_info.value.code == "INVALID_RESULT_EXPORT"
    assert not (tmp_path / "bad.json").exists()


@dataclass(frozen=True)
class StructurallyInvalidResult:
    value: object

    def to_dict(self) -> dict[str, Any]:
        return {"details": self.value}


@pytest.mark.parametrize("kind", ["cycle", "deep"])
def test_result_export_wraps_core_json_structure_errors(tmp_path: Path, kind: str) -> None:
    nested: list[object] = []
    if kind == "cycle":
        nested.append(nested)
    else:
        current = nested
        for _ in range(65):
            child: list[object] = []
            current.append(child)
            current = child

    target = tmp_path / f"{kind}.json"
    with pytest.raises(InputFormatError) as exc_info:
        export_result(StructurallyInvalidResult(nested), target)

    assert exc_info.value.code == "INVALID_RESULT_EXPORT"
    assert exc_info.value.context == {"field": "details"}
    assert not target.exists()


@pytest.mark.parametrize("format", ["csv", "tsv", "json"])
def test_generic_table_round_trip_preserves_schema_and_types(tmp_path: Path, format: str) -> None:
    schema = TableSchema(
        ("id", "count", "score", "enabled"),
        "typed.v1",
        {"id": "string", "count": "integer", "score": "number", "enabled": "boolean"},
        ("score",),
    )
    target = tmp_path / f"table.{format}"
    table_format = cast(TableFormat, format)
    export_table(
        (
            {"id": "a", "count": 1, "score": 1.25, "enabled": True},
            {"id": "b", "count": 2, "score": None, "enabled": False},
        ),
        target,
        format=table_format,
        schema=schema,
    )

    result = read_table(target, format=table_format, schema=schema)

    assert result.row_count == 2
    assert dict(result.rows[0]) == {"id": "a", "count": 1, "score": 1.25, "enabled": True}
    assert dict(result.rows[1]) == {"id": "b", "count": 2, "score": None, "enabled": False}


@pytest.mark.parametrize("format", ["csv", "tsv"])
def test_delimited_table_distinguishes_empty_string_from_null(tmp_path: Path, format: str) -> None:
    schema = TableSchema(("id", "value"), column_types={"value": "string"})
    target = tmp_path / f"nulls.{format}"
    export_table(
        ({"id": "empty", "value": ""}, {"id": "null", "value": None}),
        target,
        format=cast(TableFormat, format),
        schema=schema,
    )

    result = read_table(target, format=cast(TableFormat, format), schema=schema)

    assert [dict(row) for row in result.rows] == [
        {"id": "empty", "value": ""},
        {"id": "null", "value": None},
    ]

    with pytest.raises(InputFormatError) as collision_error:
        export_table(
            ({"id": "literal", "value": r"\N"},),
            tmp_path / f"collision.{format}",
            format=cast(TableFormat, format),
            schema=schema,
        )
    assert collision_error.value.code == "TABLE_NULL_SENTINEL_COLLISION"


def test_generic_table_reader_enforces_schema_missing_and_resource_limits(tmp_path: Path) -> None:
    target = tmp_path / "table.csv"
    target.write_text("id,count\na,\\N\n", encoding="utf-8")
    required = TableSchema(("id", "count"), column_types={"count": "integer"}, nullable=())

    with pytest.raises(InputFormatError) as null_error:
        read_table(target, format="csv", schema=required)
    assert null_error.value.code == "TABLE_NULL_NOT_ALLOWED"

    with pytest.raises(InputFormatError) as size_error:
        read_table(target, format="csv", schema=required, max_file_bytes=2)
    assert size_error.value.code == "TABLE_FILE_LIMIT_EXCEEDED"

    target.write_text("id,count\na,1\nb,2\n", encoding="utf-8")
    with pytest.raises(InputFormatError) as row_error:
        read_table(target, format="csv", schema=required, max_rows=1)
    assert row_error.value.code == "TABLE_ROW_LIMIT_EXCEEDED"

    with pytest.raises(InputFormatError) as decoded_error:
        read_table(
            target,
            format="csv",
            schema=required,
            max_decoded_bytes=1,
        )
    assert decoded_error.value.code == "TABLE_DECODED_LIMIT_EXCEEDED"


def test_parquet_generic_table_round_trip_with_installed_backend(tmp_path: Path) -> None:
    if not parquet_backend_status().available:
        pytest.skip("optional pyarrow backend is not installed")
    schema = TableSchema(
        ("id", "score"),
        "typed.parquet.v1",
        {"id": "string", "score": "number"},
        ("score",),
    )
    target = tmp_path / "typed.parquet"
    exported = export_table(
        ({"id": "a", "score": 1.5}, {"id": "b", "score": None}),
        target,
        format="parquet",
        schema=schema,
        parquet_compression="gzip",
    )

    result = read_table(target, format="parquet", schema=schema)

    assert result.row_count == 2
    assert exported.parameters["backend"] == {
        "name": "pyarrow",
        "version": parquet_backend_status().version,
        "implementation": "conditional-adapter",
        "license": "Apache-2.0",
        "compression": "gzip",
    }
    assert result.parameters["backend"] == {
        "name": "pyarrow",
        "version": parquet_backend_status().version,
        "implementation": "conditional-adapter",
        "license": "Apache-2.0",
        "compression": ("gzip",),
    }
    assert [dict(row) for row in result.rows] == [
        {"id": "a", "score": 1.5},
        {"id": "b", "score": None},
    ]

    with pytest.raises(InputFormatError) as decoded_error:
        read_table(
            target,
            format="parquet",
            schema=schema,
            max_decoded_bytes=1,
        )
    assert decoded_error.value.code == "TABLE_DECODED_LIMIT_EXCEEDED"

    limited = tmp_path / "limited.parquet"
    with pytest.raises(InputFormatError) as output_error:
        export_table(
            ({"id": "a", "score": 1.5},),
            limited,
            format="parquet",
            schema=schema,
            max_output_bytes=10,
        )
    assert output_error.value.code == "TABLE_OUTPUT_LIMIT_EXCEEDED"
    assert not limited.exists()
    assert list(tmp_path.glob(f".{limited.name}.*.tmp")) == []


def test_json_table_reader_rejects_duplicate_keys_and_nonstandard_constants(
    tmp_path: Path,
) -> None:
    schema = TableSchema(("id",))
    target = tmp_path / "invalid.json"
    target.write_text(
        '{"schema_version":"dnakit.table.v1","columns":["id"],"rows":[{"id":"a","id":"b"}]}',
        encoding="utf-8",
    )

    with pytest.raises(InputFormatError) as duplicate_error:
        read_table(target, format="json", schema=schema)
    assert duplicate_error.value.code == "INVALID_TABLE_JSON"

    target.write_text(
        '{"schema_version":"dnakit.table.v1","columns":["id"],"rows":[{"id":NaN}]}',
        encoding="utf-8",
    )
    with pytest.raises(InputFormatError) as constant_error:
        read_table(target, format="json", schema=schema)
    assert constant_error.value.code == "INVALID_TABLE_JSON"
