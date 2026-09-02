from __future__ import annotations

import gzip
import io

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import InputFormatError, RecordSourceClosedError
from dnakit.io import ReadConfig, RecordSource, WriteConfig, read, read_one, write


def _record(record_id: str, sequence: str = "ACGT") -> DNARecord:
    return DNARecord(DNASequence(sequence), record_id)


def test_record_source_is_single_use_and_releases_callback_once() -> None:
    calls: list[str] = []
    source = RecordSource(iter((_record("one"),)), close_callback=lambda: calls.append("close"))

    assert iter(source) is source
    assert [record.id for record in source] == ["one"]
    assert source.closed
    assert source.exhausted
    assert calls == ["close"]
    with pytest.raises(StopIteration):
        next(source)

    source.close()
    assert calls == ["close"]


def test_manual_close_is_idempotent_and_later_consumption_fails() -> None:
    calls: list[str] = []
    source = RecordSource(
        iter((_record("one"), _record("two"))), close_callback=lambda: calls.append("close")
    )
    assert next(source).id == "one"

    source.close()
    source.close()

    assert calls == ["close"]
    assert source.closed
    assert not source.exhausted
    with pytest.raises(RecordSourceClosedError) as exc_info:
        next(source)
    assert exc_info.value.code == "RECORD_SOURCE_CLOSED"
    with pytest.raises(RecordSourceClosedError):
        source.collect()


def test_context_manager_closes_active_source() -> None:
    calls: list[str] = []
    source = RecordSource(
        iter((_record("one"), _record("two"))), close_callback=lambda: calls.append("close")
    )

    with source as records:
        assert next(records).id == "one"

    assert calls == ["close"]
    with pytest.raises(RecordSourceClosedError):
        next(source)


def test_collect_materializes_only_remaining_records() -> None:
    source = RecordSource(iter((_record("one"), _record("two"), _record("three"))))
    assert next(source).id == "one"

    dataset = source.collect()

    assert dataset.ids == ("two", "three")
    assert source.exhausted


def test_borrowed_text_stream_stays_open_by_default() -> None:
    stream = io.StringIO(">one\nACGT\n")

    assert read_one(stream, format="fasta").id == "one"

    assert not stream.closed


def test_explicitly_owned_text_stream_is_closed() -> None:
    stream = io.StringIO(">one\nACGT\n")

    assert read_one(stream, format="fasta", config=ReadConfig(close_source=True)).id == "one"

    assert stream.closed


def test_borrowed_binary_stream_stays_open_and_is_detached() -> None:
    stream = io.BytesIO(b">one\nACGT\n")

    assert read_one(stream, format="fasta").sequence.to_string() == "ACGT"

    assert not stream.closed
    assert stream.tell() == len(stream.getvalue())


def test_borrowed_gzip_binary_stream_roundtrip_stays_open() -> None:
    target = io.BytesIO()
    write(
        DNASequence("ACGT"),
        target,
        format="fasta",
        config=WriteConfig(compression="gzip"),
    )

    assert not target.closed
    assert gzip.decompress(target.getvalue()) == b">sequence_1\nACGT\n"
    target.seek(0)
    record = read_one(
        target,
        format="fasta",
        config=ReadConfig(compression="gzip"),
    )
    assert record.sequence.to_string() == "ACGT"
    assert not target.closed


def test_owned_output_stream_is_closed_after_write() -> None:
    stream = io.StringIO()

    write(
        DNASequence("A"),
        stream,
        format="fasta",
        config=WriteConfig(close_target=True),
    )

    assert stream.closed


def test_invalid_utf8_and_invalid_gzip_have_structured_errors() -> None:
    with pytest.raises(InputFormatError) as decode_error:
        source = read(io.BytesIO(b">one\n\xff\n"), format="fasta")
        tuple(source)
    assert decode_error.value.code == "INPUT_DECODE_FAILED"

    with pytest.raises(InputFormatError) as gzip_error:
        compressed = read(
            io.BytesIO(b"not gzip"),
            format="fasta",
            config=ReadConfig(compression="gzip"),
        )
        tuple(compressed)
    assert gzip_error.value.code == "INPUT_READ_FAILED"


@pytest.mark.parametrize(
    ("text", "record_count"),
    [("", 0), (">one\nA\n>two\nC\n", "at_least_2")],
)
def test_read_one_rejects_any_count_other_than_one(text: str, record_count: object) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_one(io.StringIO(text), format="fasta")

    assert exc_info.value.code == "EXPECTED_ONE_RECORD"
    assert exc_info.value.context["record_count"] == record_count


def test_stream_without_filename_requires_explicit_format() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read(io.StringIO(">one\nACGT\n"))

    assert exc_info.value.code == "FORMAT_REQUIRED"
