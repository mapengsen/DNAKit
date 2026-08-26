from __future__ import annotations

import gzip
import json
import os
from collections.abc import Iterator, Mapping
from pathlib import Path

import pytest

from dnakit.exceptions import DuplicateIDError, InputFormatError
from dnakit.io import (
    build_fasta_index,
    build_fastq_index,
    iter_chunks,
    load_fasta_index,
    load_fastq_index,
)


def test_iter_chunks_is_lazy_bounded_and_stable() -> None:
    consumed: list[int] = []

    def values() -> Iterator[int]:
        for value in range(5):
            consumed.append(value)
            yield value

    chunks = iter_chunks(values(), chunk_size=2)

    assert next(chunks) == (0, 1)
    assert consumed == [0, 1]
    assert list(chunks) == [(2, 3), (4,)]


def test_fasta_index_persistence_id_and_coordinate_fetch(tmp_path: Path) -> None:
    fasta = tmp_path / "records.fa"
    index_path = tmp_path / "records.index.json"
    fasta.write_text(">alpha first\nACGT\nAC\n>beta\nNNRY\n", encoding="utf-8")

    built = build_fasta_index(fasta, index_path)
    loaded = load_fasta_index(index_path)

    assert built.ids == loaded.ids == ("alpha", "beta")
    assert loaded.fetch("alpha").sequence.to_string() == "ACGTAC"
    selected = loaded.fetch("alpha", start=1, end=5, strand="reverse")
    assert selected.sequence.to_string() == "TACG"
    query = selected.metadata["fasta_index_query"]
    assert isinstance(query, Mapping)
    assert query["coordinate_system"] == "0-based-half-open"
    assert "_by_id" not in json.loads(index_path.read_text(encoding="utf-8"))


def test_fasta_index_detects_stale_checksum_and_query_bounds(tmp_path: Path) -> None:
    fasta = tmp_path / "records.fa"
    index_path = tmp_path / "records.index.json"
    fasta.write_text(">alpha\nAAAA\n", encoding="utf-8")
    index = build_fasta_index(fasta, index_path)

    with pytest.raises(InputFormatError) as bounds_error:
        index.fetch("alpha", start=0, end=5)
    assert bounds_error.value.code == "FASTA_QUERY_OUT_OF_BOUNDS"

    fasta.write_text(">alpha\nCCCC\n", encoding="utf-8")
    with pytest.raises(InputFormatError) as stale_error:
        load_fasta_index(index_path)
    assert stale_error.value.code == "STALE_FASTA_INDEX"


def test_fasta_index_rejects_duplicate_ids_and_gzip(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.fa"
    duplicate.write_text(">x\nA\n>x\nC\n", encoding="utf-8")
    with pytest.raises(DuplicateIDError):
        build_fasta_index(duplicate)

    compressed = tmp_path / "records.fa.gz"
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write(">x\nA\n")
    with pytest.raises(InputFormatError) as exc_info:
        build_fasta_index(compressed)
    assert exc_info.value.code == "COMPRESSED_FASTA_INDEX_UNSUPPORTED"


def test_fasta_fetch_enforces_record_byte_limit(tmp_path: Path) -> None:
    fasta = tmp_path / "record.fa"
    fasta.write_text(">x\nACGT\n", encoding="utf-8")
    index = build_fasta_index(fasta)

    with pytest.raises(InputFormatError) as exc_info:
        index.fetch("x", max_record_bytes=2)

    assert exc_info.value.code == "FASTA_RECORD_LIMIT_EXCEEDED"


def test_fasta_index_json_types_are_not_silently_coerced(tmp_path: Path) -> None:
    fasta = tmp_path / "record.fa"
    index_path = tmp_path / "record.index.json"
    fasta.write_text(">x\nA\n", encoding="utf-8")
    build_fasta_index(fasta, index_path)
    payload = json.loads(index_path.read_text(encoding="utf-8"))
    payload["entries"][0]["byte_start"] = True
    index_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(InputFormatError) as exc_info:
        load_fasta_index(index_path)

    assert exc_info.value.code == "INVALID_FASTA_INDEX"


def test_fastq_index_persistence_fetch_and_quality_synchronized_reverse_slice(
    tmp_path: Path,
) -> None:
    fastq = tmp_path / "records.fastq"
    index_path = tmp_path / "records.index.json"
    fastq.write_text(
        "@alpha first\nACGTAC\n+\n!+5?IS\n@beta\nNNRY\n+beta\nIIII\n", encoding="utf-8"
    )

    built = build_fastq_index(fastq, index_path)
    loaded = load_fastq_index(index_path)

    assert built.ids == loaded.ids == ("alpha", "beta")
    assert loaded.fetch("alpha").letter_annotations["phred_quality"] == (0, 10, 20, 30, 40, 50)
    selected = loaded.fetch("alpha", start=1, end=5, strand="reverse")
    assert selected.sequence.to_string() == "TACG"
    assert selected.letter_annotations["phred_quality"] == (40, 30, 20, 10)
    assert selected.metadata["fastq_index_query"] == {
        "start": 1,
        "end": 5,
        "strand": "reverse",
        "coordinate_system": "0-based-half-open",
    }


def test_fastq_index_rejects_duplicate_gzip_malformed_and_resource_limits(
    tmp_path: Path,
) -> None:
    duplicate = tmp_path / "duplicate.fastq"
    duplicate.write_text("@x\nA\n+\nI\n@x\nC\n+\nI\n", encoding="utf-8")
    with pytest.raises(DuplicateIDError):
        build_fastq_index(duplicate)

    compressed = tmp_path / "records.fastq.gz"
    with gzip.open(compressed, "wt", encoding="utf-8") as handle:
        handle.write("@x\nA\n+\nI\n")
    with pytest.raises(InputFormatError) as gzip_error:
        build_fastq_index(compressed)
    assert gzip_error.value.code == "COMPRESSED_FASTQ_INDEX_UNSUPPORTED"

    malformed = tmp_path / "malformed.fastq"
    malformed.write_text("@x\nAC\n+\nI\n", encoding="utf-8")
    with pytest.raises(InputFormatError) as format_error:
        build_fastq_index(malformed)
    assert format_error.value.code == "FASTQ_QUALITY_LENGTH_MISMATCH"

    with pytest.raises(InputFormatError) as source_error:
        build_fastq_index(duplicate, max_source_bytes=1)
    assert source_error.value.code == "FASTQ_SOURCE_LIMIT_EXCEEDED"


def test_fastq_index_detects_stale_source_bounds_and_strict_index_schema(tmp_path: Path) -> None:
    fastq = tmp_path / "record.fastq"
    index_path = tmp_path / "record.index.json"
    fastq.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
    index = build_fastq_index(fastq, index_path)

    with pytest.raises(InputFormatError) as bounds_error:
        index.fetch("x", end=5)
    assert bounds_error.value.code == "FASTQ_QUERY_OUT_OF_BOUNDS"

    payload = json.loads(index_path.read_text(encoding="utf-8"))
    payload["unexpected"] = True
    index_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(InputFormatError) as schema_error:
        load_fastq_index(index_path)
    assert schema_error.value.code == "INVALID_FASTQ_INDEX"

    build_fastq_index(fastq, index_path, overwrite=True)
    fastq.write_text("@x\nTGCA\n+\nIIII\n", encoding="utf-8")
    with pytest.raises(InputFormatError) as stale_error:
        load_fastq_index(index_path)
    assert stale_error.value.code == "STALE_FASTQ_INDEX"


def test_fastq_index_loader_rejects_duplicate_json_keys_and_bounded_file(tmp_path: Path) -> None:
    fastq = tmp_path / "record.fastq"
    index_path = tmp_path / "record.index.json"
    fastq.write_text("@x\nA\n+\nI\n", encoding="utf-8")
    build_fastq_index(fastq, index_path)

    with pytest.raises(InputFormatError) as size_error:
        load_fastq_index(index_path, max_index_bytes=1)
    assert size_error.value.code == "FASTQ_INDEX_FILE_LIMIT_EXCEEDED"

    payload = index_path.read_text(encoding="utf-8").rstrip()
    index_path.write_text(
        payload[:-1] + ',"schema_version":"dnakit.fastq-index.v1"}', encoding="utf-8"
    )
    with pytest.raises(InputFormatError) as duplicate_error:
        load_fastq_index(index_path)
    assert duplicate_error.value.code == "INVALID_FASTQ_INDEX"


def test_fastq_index_loader_requires_original_mtime_even_when_bytes_match(tmp_path: Path) -> None:
    fastq = tmp_path / "record.fastq"
    index_path = tmp_path / "record.index.json"
    fastq.write_text("@x\nA\n+\nI\n", encoding="utf-8")
    build_fastq_index(fastq, index_path)
    stat = fastq.stat()
    os.utime(fastq, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

    with pytest.raises(InputFormatError) as stale_error:
        load_fastq_index(index_path)

    assert stale_error.value.code == "STALE_FASTQ_INDEX"
