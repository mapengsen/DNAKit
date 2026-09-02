from __future__ import annotations

import gzip
import hashlib
import io
import os
from pathlib import Path
from typing import Literal

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence
from dnakit.exceptions import InputFormatError
from dnakit.io import ReadConfig, WriteConfig, read, read_one, read_set, write


def test_read_multiline_fasta_with_description_and_iupac() -> None:
    stream = io.StringIO(">alpha first sample\nACGT\nNRY\n>beta\nTTAA\n")

    records = read_set(stream, format="fasta")

    assert records.ids == ("alpha", "beta")
    assert records[0].description == "first sample"
    assert records[0].sequence.to_string() == "ACGTNRY"
    assert records[0].sequence.alphabet is DNAAlphabet.IUPAC


def test_fasta_does_not_silently_strip_sequence_spaces() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(">alpha\n ACGT\n"), format="fasta")

    assert exc_info.value.code == "INVALID_SEQUENCE_CONTENT"


def test_streaming_sequence_line_and_record_limits_are_enforced() -> None:
    with pytest.raises(InputFormatError) as line_error:
        read_set(
            io.StringIO(">alpha\nACGT\n"),
            format="fasta",
            config=ReadConfig(max_field_size=3),
        )
    assert line_error.value.code == "SEQUENCE_LINE_TOO_LONG"

    with pytest.raises(InputFormatError) as record_error:
        read_set(
            io.StringIO(">alpha\nA\nC\n"),
            format="fasta",
            config=ReadConfig(max_record_lines=2),
        )
    assert record_error.value.code == "SEQUENCE_RECORD_LINE_LIMIT_EXCEEDED"

    with pytest.raises(InputFormatError) as fastq_error:
        read_set(
            io.StringIO("@alpha\nACGT\n+\n!!!!\n"),
            format="fastq",
            config=ReadConfig(max_field_size=3),
        )
    assert fastq_error.value.code == "SEQUENCE_LINE_TOO_LONG"


@pytest.mark.parametrize(
    ("format", "text"),
    [
        ("fasta", ">alpha\nACGT\n"),
        ("fastq", "@alpha\nACGT\n+\n!!!!\n"),
        ("jsonl", '{"id":"alpha","sequence":"ACGT"}\n'),
        ("json", '[{"id":"alpha","sequence":"ACGT"}]'),
        ("csv", "id,sequence\nalpha,ACGT\n"),
        ("tsv", "id\tsequence\nalpha\tACGT\n"),
    ],
)
def test_all_sequence_readers_enforce_symbol_limit(format: str, text: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format=format, config=ReadConfig(max_sequence_symbols=3))

    assert exc_info.value.code == "SEQUENCE_SYMBOL_LIMIT_EXCEEDED"


def test_input_byte_limit_counts_decoded_gzip_payload(tmp_path: Path) -> None:
    path = tmp_path / "large.fa.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(">alpha\n" + "A" * 100 + "\n")

    with pytest.raises(InputFormatError) as exc_info:
        read_set(path, config=ReadConfig(max_input_bytes=20))

    assert exc_info.value.code == "INPUT_BYTE_LIMIT_EXCEEDED"


@pytest.mark.parametrize(("suffix", "compression"), [(".fa", "none"), (".fa.gz", "gzip")])
def test_path_output_byte_limit_is_atomic_for_plain_and_gzip(
    tmp_path: Path, suffix: str, compression: Literal["none", "gzip"]
) -> None:
    path = tmp_path / f"limited{suffix}"

    with pytest.raises(InputFormatError) as exc_info:
        write(
            DNARecord(DNASequence("A" * 100), "alpha"),
            path,
            config=WriteConfig(compression=compression, max_output_bytes=10),
        )

    assert exc_info.value.code == "OUTPUT_BYTE_LIMIT_EXCEEDED"
    assert not path.exists()
    assert list(tmp_path.glob(f".{path.name}.*.tmp")) == []


def test_text_and_binary_stream_output_byte_limits_raise_before_exceeding() -> None:
    text_target = io.StringIO()
    with pytest.raises(InputFormatError) as text_error:
        write(
            DNARecord(DNASequence("A" * 20), "alpha"),
            text_target,
            format="fasta",
            config=WriteConfig(max_output_bytes=10),
        )
    assert text_error.value.code == "OUTPUT_BYTE_LIMIT_EXCEEDED"
    assert len(text_target.getvalue().encode()) <= 10

    binary_target = io.BytesIO()
    with pytest.raises(InputFormatError) as binary_error:
        write(
            DNARecord(DNASequence("A" * 100), "alpha"),
            binary_target,
            format="fasta",
            config=WriteConfig(compression="gzip", max_output_bytes=10),
        )
    assert binary_error.value.code == "OUTPUT_BYTE_LIMIT_EXCEEDED"
    assert len(binary_target.getvalue()) <= 10


@pytest.mark.parametrize("format", ["fasta", "fastq", "csv", "tsv", "json", "jsonl", "genbank"])
def test_every_record_format_honors_output_byte_limit(tmp_path: Path, format: str) -> None:
    path = tmp_path / f"limited-{format}.out"
    record = DNARecord(
        DNASequence("ACGT"),
        "alpha",
        letter_annotations={"phred_quality": (0, 1, 2, 3)},
    )

    with pytest.raises(InputFormatError) as exc_info:
        write(
            record,
            path,
            format=format,
            config=WriteConfig(max_output_bytes=5),
        )

    assert exc_info.value.code == "OUTPUT_BYTE_LIMIT_EXCEEDED"
    assert not path.exists()


def test_fasta_path_format_inference_and_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "records.fasta"
    original = [
        DNARecord(DNASequence("ACGT"), "alpha", description="first"),
        DNARecord(DNASequence("NNRY", alphabet=DNAAlphabet.IUPAC), "beta"),
    ]

    result = write(original, path, config=WriteConfig(line_width=2))
    restored = read_set(path)

    assert result.format == "fasta"
    assert result.record_count == 2
    assert result.byte_count == path.stat().st_size
    assert result.target_artifact is not None
    assert result.target_artifact.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert not os.path.isabs(result.target_artifact.relative_path)
    assert [record.sequence.to_string() for record in restored] == ["ACGT", "NNRY"]
    assert restored[0].description == "first"
    assert path.read_text(encoding="utf-8").splitlines()[1:3] == ["AC", "GT"]


def test_gzip_fasta_write_and_lazy_read(tmp_path: Path) -> None:
    path = tmp_path / "records.fasta.gz"
    records = [DNARecord(DNASequence("ACGT"), "alpha")]

    result = write(records, path)
    source = read(path)
    restored = next(source)

    assert result.byte_count == path.stat().st_size
    assert path.read_bytes()[:2] == b"\x1f\x8b"
    assert gzip.decompress(path.read_bytes()).decode() == ">alpha\nACGT\n"
    assert restored.sequence.to_string() == "ACGT"
    with pytest.raises(StopIteration):
        next(source)


def test_gzip_output_is_byte_for_byte_reproducible(tmp_path: Path) -> None:
    first = tmp_path / "first.fa.gz"
    second = tmp_path / "second.fa.gz"
    record = DNARecord(DNASequence("ACGT"), "alpha")

    write(record, first)
    write(record, second)

    assert first.read_bytes() == second.read_bytes()


def test_anonymous_sequences_receive_stable_input_order_ids() -> None:
    stream = io.StringIO()

    result = write(
        [DNASequence("A"), DNARecord(DNASequence("C"), "named"), DNASequence("G")],
        stream,
        format="fasta",
    )

    assert [(item.input_index, item.generated_id) for item in result.generated_ids] == [
        (0, "sequence_1"),
        (2, "sequence_3"),
    ]
    assert stream.getvalue() == ">sequence_1\nA\n>named\nC\n>sequence_3\nG\n"


def test_anonymous_sequence_can_be_rejected() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        write(
            DNASequence("ACGT"),
            io.StringIO(),
            format="fasta",
            config=WriteConfig(anonymous_id_policy="error"),
        )

    assert exc_info.value.code == "ANONYMOUS_SEQUENCE_ID_REQUIRED"


def test_writer_refuses_overwrite_and_can_explicitly_replace(tmp_path: Path) -> None:
    path = tmp_path / "record.fa"
    path.write_text("original", encoding="utf-8")

    with pytest.raises(FileExistsError):
        write(DNASequence("A"), path)
    assert path.read_text(encoding="utf-8") == "original"

    write(DNASequence("C"), path, config=WriteConfig(overwrite=True))
    assert path.read_text(encoding="utf-8") == ">sequence_1\nC\n"


def test_failed_atomic_write_leaves_no_output_or_temporary_file(tmp_path: Path) -> None:
    path = tmp_path / "bad.fastq"
    record_without_quality = DNARecord(DNASequence("ACGT"), "alpha")

    with pytest.raises(InputFormatError) as exc_info:
        write(record_without_quality, path)

    assert exc_info.value.code == "FASTQ_QUALITY_MISSING"
    assert not path.exists()
    assert list(tmp_path.glob(f".{path.name}.*.tmp")) == []


def test_failed_overwrite_preserves_existing_file(tmp_path: Path) -> None:
    path = tmp_path / "bad.fastq"
    path.write_bytes(b"keep-me")

    with pytest.raises(InputFormatError):
        write(
            DNARecord(DNASequence("ACGT"), "alpha"),
            path,
            config=WriteConfig(overwrite=True),
        )

    assert path.read_bytes() == b"keep-me"


def test_fastq_reads_iupac_and_boundary_phred_values() -> None:
    record = read_one(io.StringIO("@alpha sample\nANV\n+alpha\n!I~\n"), format="fastq")

    assert record.sequence.to_string() == "ANV"
    assert record.description == "sample"
    assert record.letter_annotations["phred_quality"] == (0, 40, 93)


def test_fastq_empty_sequence_and_quality_are_valid() -> None:
    record = read_one(io.StringIO("@empty\n\n+\n\n"), format="fastq")

    assert record.sequence.symbol_length == 0
    assert record.letter_annotations["phred_quality"] == ()


@pytest.mark.parametrize(
    ("text", "error_code"),
    [
        ("@alpha\nAN\n+\n!\n", "FASTQ_QUALITY_LENGTH_MISMATCH"),
        ("@alpha\nAN\n+\n !\n", "FASTQ_QUALITY_OUT_OF_RANGE"),
        ("@alpha\nAN\n+beta\n!!\n", "FASTQ_HEADER_MISMATCH"),
        ("@alpha\nAN\nnot-plus\n!!\n", "FASTQ_MISSING_SEPARATOR"),
        ("@alpha\nAN\n+\n", "TRUNCATED_FASTQ_RECORD"),
        ("@alpha\nA\n+\n\x7f\n", "FASTQ_QUALITY_OUT_OF_RANGE"),
        ("@alpha\nA\n+\n\u0080\n", "FASTQ_QUALITY_OUT_OF_RANGE"),
    ],
)
def test_fastq_rejects_malformed_records(text: str, error_code: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format="fastq")

    assert exc_info.value.code == error_code


def test_fastq_strict_alphabet_rejects_iupac_ambiguity() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_one(
            io.StringIO("@alpha\nAN\n+\n!!\n"),
            format="fastq",
            config=ReadConfig(alphabet=DNAAlphabet.STRICT),
        )

    assert exc_info.value.code == "INVALID_SEQUENCE_CONTENT"


def test_fastq_custom_phred_offset_roundtrip() -> None:
    original = DNARecord(
        DNASequence("ANV", alphabet=DNAAlphabet.IUPAC),
        "alpha",
        letter_annotations={"phred_quality": (0, 40, 62)},
    )
    stream = io.StringIO()

    write(original, stream, format="fastq", config=WriteConfig(phred_offset=64))
    stream.seek(0)
    restored = read_one(stream, format="fastq", config=ReadConfig(phred_offset=64))

    assert stream.getvalue() == "@alpha\nANV\n+\n@h~\n"
    assert restored.letter_annotations["phred_quality"] == (0, 40, 62)


def test_fastq_writer_rejects_non_integer_or_unencodable_quality() -> None:
    float_quality = DNARecord(
        DNASequence("A"), "float", letter_annotations={"phred_quality": (1.5,)}
    )
    high_quality = DNARecord(DNASequence("A"), "high", letter_annotations={"phred_quality": (93,)})

    with pytest.raises(InputFormatError):
        write(float_quality, io.StringIO(), format="fastq")
    with pytest.raises(InputFormatError):
        write(high_quality, io.StringIO(), format="fastq", config=WriteConfig(phred_offset=64))


@pytest.mark.parametrize(
    "record",
    [
        DNARecord(DNASequence("A"), "two words"),
        DNARecord(DNASequence("A"), "one", description="line one\nline two"),
    ],
)
def test_fasta_and_fastq_reject_unsafe_headers(record: DNARecord) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        write(record, io.StringIO(), format="fasta")

    assert exc_info.value.code == "INVALID_SEQUENCE_HEADER"


@pytest.mark.parametrize("format", ["parquet", "fai"])
def test_advanced_formats_fail_with_explicit_unsupported_error(format: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read(io.StringIO(""), format=format)

    assert exc_info.value.code == "UNSUPPORTED_FORMAT"
    assert exc_info.value.hint is not None
