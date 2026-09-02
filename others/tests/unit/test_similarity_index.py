"""Tests for sketch index construction, search, and persistence."""

import json
from collections.abc import Iterator

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.similarity import (
    build_sketch_index,
    load_sketch_index,
    nearest_neighbors,
    save_sketch_index,
)
from dnakit.similarity.index import SketchIndex


def _records() -> DNASet:
    return DNASet.from_records(
        [
            DNARecord(DNASequence("ACGTACGT"), "a"),
            DNARecord(DNASequence("ACGTTCGT"), "b"),
            DNARecord(DNASequence("TTTTTTTT"), "c"),
        ]
    )


def test_index_nearest_neighbors_are_stable_and_auditable() -> None:
    index = build_sketch_index(_records(), k=3, num_hashes=100, seed=9)
    result = nearest_neighbors(_records()[0], index, top_k=2)

    assert [hit.record_id for hit in result.hits] == ["a", "b"]
    assert result.hits[0].similarity == 1
    assert result.parameters["seed"] == 9


def test_index_json_roundtrip_and_checksum(tmp_path: object) -> None:
    from pathlib import Path

    target = Path(str(tmp_path)) / "index.json"
    index = build_sketch_index(_records(), k=2, num_hashes=10)
    digest = save_sketch_index(index, target)

    assert len(digest) == 64
    assert load_sketch_index(target) == index
    envelope = json.loads(target.read_text(encoding="utf-8"))
    envelope["payload"]["ids"][0] = "tampered"
    target.write_text(json.dumps(envelope), encoding="utf-8")
    with pytest.raises(InputFormatError) as error:
        load_sketch_index(target)
    assert error.value.code == "SKETCH_INDEX_CHECKSUM"


def test_index_bounds_generator_and_rejects_duplicates() -> None:
    consumed = 0

    def records() -> Iterator[DNARecord]:
        nonlocal consumed
        while True:
            consumed += 1
            yield DNARecord(DNASequence("AC"), f"r{consumed}")

    with pytest.raises(ConfigurationError) as limit:
        build_sketch_index(records(), k=1, max_records=2)
    assert limit.value.code == "INDEX_RECORD_LIMIT_EXCEEDED"
    assert consumed == 3
    with pytest.raises(ConfigurationError) as duplicate:
        build_sketch_index(
            [
                DNARecord(DNASequence("AC"), "same"),
                DNARecord(DNASequence("GT"), "same"),
            ],
            k=1,
        )
    assert duplicate.value.code == "DUPLICATE_INDEX_ID"


def test_index_save_refuses_overwrite(tmp_path: object) -> None:
    from pathlib import Path

    target = Path(str(tmp_path)) / "index.json"
    index = build_sketch_index(_records(), k=2)
    save_sketch_index(index, target)
    with pytest.raises(ConfigurationError) as error:
        save_sketch_index(index, target)
    assert error.value.code == "INDEX_EXISTS"


def test_empty_index_still_validates_schema_parameters() -> None:
    with pytest.raises(ConfigurationError):
        SketchIndex((), (), 0, 10, True, 0)
    with pytest.raises(ConfigurationError):
        SketchIndex((), (), 2, 10, True, 2**64)


def test_index_loader_rejects_tampered_envelope_schema(tmp_path: object) -> None:
    from pathlib import Path

    target = Path(str(tmp_path)) / "index.json"
    save_sketch_index(build_sketch_index(_records(), k=2), target)
    envelope = json.loads(target.read_text(encoding="utf-8"))
    envelope["schema_version"] = "unknown"
    target.write_text(json.dumps(envelope), encoding="utf-8")

    with pytest.raises(InputFormatError) as error:
        load_sketch_index(target)
    assert error.value.code == "INVALID_SKETCH_INDEX"


def test_index_loader_rejects_duplicate_keys_and_deep_json(tmp_path: object) -> None:
    from pathlib import Path

    duplicate = Path(str(tmp_path)) / "duplicate.json"
    duplicate.write_text('{"schema_version":"a","schema_version":"b"}', encoding="utf-8")
    deep = Path(str(tmp_path)) / "deep.json"
    deep.write_text("[" * 1_200 + "]" * 1_200, encoding="utf-8")

    for target in (duplicate, deep):
        with pytest.raises(InputFormatError) as error:
            load_sketch_index(target)
        assert error.value.code == "INVALID_SKETCH_INDEX"
