from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from dnakit.core import (
    DNAAlphabet,
    DNAFeature,
    DNARecord,
    DNASequence,
    Gap,
    GapKind,
    Interval,
    Strand,
    Strandedness,
    Topology,
)
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io import ReadConfig, WriteConfig, read_one, read_set, write


def _records() -> list[DNARecord]:
    return [
        DNARecord(
            DNASequence(
                "ACGN",
                alphabet=DNAAlphabet.IUPAC,
                topology=Topology.CIRCULAR,
                strandedness=Strandedness.DOUBLE,
            ),
            "alpha",
            description="first sample",
            features=[
                DNAFeature(
                    "motif",
                    Interval(1, 3),
                    id="feature-1",
                    strand=Strand.FORWARD,
                    label="CG",
                    score=0.75,
                    qualifiers={"source": ["manual"]},
                    source="test",
                )
            ],
            metadata={"species": "human", "batch": 2},
            letter_annotations={"phred_quality": (10, 20, 30, 40)},
        ),
        DNARecord(DNASequence("TTAA"), "beta", metadata={"species": "mouse"}),
    ]


@pytest.mark.parametrize("format", ["csv", "tsv", "json", "jsonl"])
def test_structured_format_roundtrip(format: str) -> None:
    stream = io.StringIO()

    result = write(_records(), stream, format=format)
    stream.seek(0)
    restored = read_set(stream, format=format)

    assert result.record_count == 2
    assert restored.ids == ("alpha", "beta")
    assert restored[0].sequence.to_string() == "ACGN"
    assert restored[0].sequence.topology is Topology.CIRCULAR
    assert restored[0].sequence.strandedness is Strandedness.DOUBLE
    assert restored[0].description == "first sample"
    assert restored[0].features == _records()[0].features
    assert restored[0].metadata["species"] == "human"
    assert restored[0].metadata["batch"] == 2
    assert restored[0].letter_annotations["phred_quality"] == (10, 20, 30, 40)


def test_json_format_is_standard_array_not_json_lines() -> None:
    stream = io.StringIO()

    write(_records(), stream, format="json")

    payload = json.loads(stream.getvalue())
    assert isinstance(payload, list)
    assert [item["id"] for item in payload] == ["alpha", "beta"]


def test_csv_and_tsv_delimiters_are_distinct() -> None:
    csv_stream = io.StringIO()
    tsv_stream = io.StringIO()

    write(_records()[:1], csv_stream, format="csv")
    write(_records()[:1], tsv_stream, format="tsv")

    assert csv_stream.getvalue().splitlines()[0].count(",") == 10
    assert tsv_stream.getvalue().splitlines()[0].count("\t") == 10


def test_csv_path_inference_and_gzip_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "records.csv.gz"

    write(_records(), path)
    restored = read_set(path)

    assert restored.ids == ("alpha", "beta")
    assert restored[0].metadata["species"] == "human"


def test_structured_input_generates_id_for_missing_value() -> None:
    record = read_one(io.StringIO('[{"sequence":"ACGT","metadata":{"batch":1}}]'), format="json")

    assert record.id == "sequence_1"
    assert record.metadata["batch"] == 1


def test_json_requires_array_root() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO('{"id":"one","sequence":"ACGT"}'), format="json")

    assert exc_info.value.code == "INVALID_JSON_ROOT"


def test_csv_requires_sequence_column() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO("id,value\none,ACGT\n"), format="csv")

    assert exc_info.value.code == "MISSING_SEQUENCE_COLUMN"


def test_csv_rejects_duplicate_headers_and_extra_fields() -> None:
    with pytest.raises(InputFormatError) as duplicate_error:
        read_set(io.StringIO("id,sequence,sequence\nx,AC,GT\n"), format="csv")
    assert duplicate_error.value.code == "INVALID_TABLE_HEADER"

    with pytest.raises(InputFormatError) as extra_error:
        read_set(io.StringIO("id,sequence\nx,AC,SECRET\n"), format="csv")
    assert extra_error.value.code == "EXTRA_TABLE_FIELDS"


def test_invalid_json_metadata_cell_is_rejected() -> None:
    text = 'id,sequence,metadata\none,ACGT,"{not-json}"\n'

    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format="csv")

    assert exc_info.value.code == "INVALID_JSON_COLUMN"


@pytest.mark.parametrize("format", ["json", "jsonl"])
def test_structured_json_rejects_duplicate_object_keys(format: str) -> None:
    row = '{"id":"one","id":"two","sequence":"A"}'
    text = f"[{row}]" if format == "json" else row + "\n"

    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format=format)

    assert exc_info.value.code == "INVALID_JSON"


def test_delimited_embedded_json_rejects_duplicate_keys() -> None:
    text = 'id,sequence,metadata\none,A,"{ ""tag"": 1, ""tag"": 2 }"\n'

    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format="csv")

    assert exc_info.value.code == "INVALID_JSON_COLUMN"


@pytest.mark.parametrize("format", ["json", "jsonl"])
def test_structured_json_converts_excessive_nesting_to_stable_error(format: str) -> None:
    nested = "[" * 1_200 + "0" + "]" * 1_200
    row = f'{{"id":"one","sequence":"A","metadata":{nested}}}'
    text = f"[{row}]" if format == "json" else row + "\n"

    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(text), format=format)

    assert exc_info.value.code == "JSON_STRUCTURE_LIMIT_EXCEEDED"


def test_structured_json_enforces_configured_depth_and_node_limits() -> None:
    nested = {"level": {"level": {"level": "value"}}}
    depth_text = json.dumps([{"id": "one", "sequence": "A", "metadata": nested}])
    with pytest.raises(InputFormatError) as depth_error:
        read_set(
            io.StringIO(depth_text),
            format="json",
            config=ReadConfig(max_json_depth=4),
        )
    assert depth_error.value.code == "JSON_STRUCTURE_LIMIT_EXCEEDED"

    nodes_text = json.dumps({"id": "one", "sequence": "A", "metadata": {"values": list(range(10))}})
    with pytest.raises(InputFormatError) as node_error:
        read_set(
            io.StringIO(nodes_text + "\n"),
            format="jsonl",
            config=ReadConfig(max_json_nodes=8),
        )
    assert node_error.value.code == "JSON_STRUCTURE_LIMIT_EXCEEDED"


def test_binary_output_stream_is_borrowed_by_default() -> None:
    stream = io.BytesIO()

    write(_records()[:1], stream, format="json")

    assert not stream.closed
    assert json.loads(stream.getvalue())[0]["id"] == "alpha"


def test_json_lines_anonymous_ids_follow_record_order_not_blank_lines() -> None:
    stream = io.StringIO('\n\n{"sequence":"A"}\n\n{"sequence":"C"}\n')

    restored = read_set(stream, format="jsonl")

    assert restored.ids == ("sequence_1", "sequence_2")


def test_json_lines_respects_configured_line_limit() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(
            io.StringIO('{"sequence":"ACGT"}\n'),
            format="jsonl",
            config=ReadConfig(max_field_size=10),
        )

    assert exc_info.value.code == "STRUCTURED_LINE_TOO_LONG"


def test_json_lines_rejects_multiline_pretty_print_configuration() -> None:
    stream = io.StringIO()

    with pytest.raises(ConfigurationError) as exc_info:
        write(
            DNASequence("AC"),
            stream,
            format="jsonl",
            config=WriteConfig(json_indent=2),
        )

    assert exc_info.value.code == "JSONL_INDENT_NOT_ALLOWED"
    assert stream.getvalue() == ""


def test_invalid_falsy_sequence_type_is_not_silently_defaulted() -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO('[{"sequence":"A","alphabet":0}]'), format="json")

    assert exc_info.value.code == "INVALID_SEQUENCE_TYPE_COLUMN"


def test_invalid_structured_feature_is_rejected() -> None:
    payload = '[{"id":"one","sequence":"AC","features":[{"type":"x"}]}]'

    with pytest.raises(InputFormatError) as exc_info:
        read_set(io.StringIO(payload), format="json")

    assert exc_info.value.code == "INVALID_FEATURES_COLUMN"


@pytest.mark.parametrize("format", ["csv", "tsv", "json", "jsonl"])
def test_structured_formats_losslessly_roundtrip_explicit_gap_parts(format: str) -> None:
    gap = Gap(
        None,
        kind=GapKind.SCAFFOLD,
        crossable=False,
        evidence=("assembly",),
        metadata={"source": "test"},
    )
    record = DNARecord(DNASequence(["AC", gap, "GT"]), "gapped")
    stream = io.StringIO()

    write(record, stream, format=format)
    stream.seek(0)
    restored = read_one(stream, format=format)

    assert restored.sequence == record.sequence
    assert restored.sequence.parts == ("AC", gap, "GT")


def test_structured_parts_mismatch_and_unknown_schema_are_rejected() -> None:
    mismatch = (
        '[{"schema_version":"dnakit.record.v1","id":"one","sequence":"AA",'
        '"parts":[{"kind":"symbols","symbols":"AC"}]}]'
    )
    with pytest.raises(InputFormatError) as mismatch_error:
        read_set(io.StringIO(mismatch), format="json")
    assert mismatch_error.value.code == "SEQUENCE_PARTS_MISMATCH"

    with pytest.raises(InputFormatError) as schema_error:
        read_set(
            io.StringIO('[{"schema_version":"future","id":"one","sequence":"AC"}]'),
            format="json",
        )
    assert schema_error.value.code == "UNSUPPORTED_RECORD_SCHEMA"


def test_long_csv_roundtrip_and_configured_field_limit() -> None:
    record = DNARecord(DNASequence("A" * 200_000), "long")
    stream = io.StringIO()
    write(record, stream, format="csv")
    stream.seek(0)

    assert read_one(stream, format="csv").sequence.symbol_length == 200_000

    stream.seek(0)
    with pytest.raises(InputFormatError) as exc_info:
        read_set(stream, format="csv", config=ReadConfig(max_field_size=1_000))
    assert exc_info.value.code == "CSV_FIELD_ERROR"


@pytest.mark.parametrize("format", ["fasta", "fastq"])
def test_lossy_formats_reject_features_unless_drop_is_explicit(format: str, tmp_path: Path) -> None:
    record = _records()[0]
    target = tmp_path / "record.fastq" if format == "fastq" else tmp_path / "record.fasta"

    with pytest.raises(InputFormatError) as exc_info:
        write(record, target, format=format)

    assert exc_info.value.code == "FEATURE_LOSS_NOT_ALLOWED"
    assert not target.exists()
    assert list(tmp_path.glob(".*.tmp")) == []

    write(record, target, format=format, config=WriteConfig(feature_policy="drop"))
    assert target.exists()


@pytest.mark.parametrize("format", ["fasta", "fastq"])
def test_sequence_only_formats_reject_explicit_gap_loss(format: str, tmp_path: Path) -> None:
    record = DNARecord(DNASequence(["A", Gap(2), "T"]), "gapped")
    target = tmp_path / f"record.{format}"

    with pytest.raises(InputFormatError) as exc_info:
        write(record, target, format=format)

    assert exc_info.value.code == "GAP_LOSS_NOT_ALLOWED"
    assert not target.exists()
