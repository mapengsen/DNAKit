from __future__ import annotations

import io
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from dnakit.core import DNAAlphabet, DNASequence, Strandedness, Topology
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io import ReadConfig, WriteConfig, WriteResult, read_set, write


def test_read_config_coerces_sequence_type_enums() -> None:
    config = ReadConfig(alphabet="strict", topology="circular", strandedness="double")  # type: ignore[arg-type]

    assert config.alphabet is DNAAlphabet.STRICT
    assert config.topology is Topology.CIRCULAR
    assert config.strandedness is Strandedness.DOUBLE


@pytest.mark.parametrize(
    "factory",
    [
        lambda: ReadConfig(encoding="not-a-real-codec"),
        lambda: ReadConfig(phred_offset=32),
        lambda: ReadConfig(delimiter="::"),
        lambda: ReadConfig(uppercase=1),  # type: ignore[arg-type]
        lambda: ReadConfig(close_source=1),  # type: ignore[arg-type]
        lambda: ReadConfig(id_column="sequence", sequence_column="sequence"),
        lambda: ReadConfig(features_column="metadata", metadata_column="metadata"),
        lambda: ReadConfig(parts_column="sequence", sequence_column="sequence"),
        lambda: ReadConfig(max_field_size=True),
        lambda: ReadConfig(delimiter="\n"),
        lambda: ReadConfig(delimiter="\x00"),
        lambda: ReadConfig(delimiter='"'),
        lambda: ReadConfig(max_field_size=sys.maxsize + 1),
        lambda: ReadConfig(max_records=0),
        lambda: ReadConfig(max_record_lines=True),
        lambda: ReadConfig(max_sequence_symbols=0),
        lambda: ReadConfig(max_input_bytes=True),
        lambda: ReadConfig(max_json_depth=0),
        lambda: ReadConfig(max_json_nodes=True),
        lambda: ReadConfig(anonymous_id_prefix="two words"),
        lambda: WriteConfig(phred_offset=65),
        lambda: WriteConfig(compression_level=10),
        lambda: WriteConfig(line_width=0),
        lambda: WriteConfig(overwrite=None),  # type: ignore[arg-type]
        lambda: WriteConfig(close_target=1),  # type: ignore[arg-type]
        lambda: WriteConfig(json_indent=-1),
        lambda: WriteConfig(max_output_bytes=0),
        lambda: WriteConfig(feature_policy="silently-ignore"),  # type: ignore[arg-type]
        lambda: WriteConfig(delimiter="\r"),
        lambda: WriteConfig(delimiter='"'),
    ],
)
def test_invalid_io_configuration_is_rejected(factory: Callable[[], object]) -> None:
    with pytest.raises(ConfigurationError):
        factory()


def test_custom_delimiter_roundtrip() -> None:
    stream = io.StringIO()
    write(DNASequence("ACGT"), stream, format="csv", config=WriteConfig(delimiter=";"))
    stream.seek(0)

    restored = read_set(stream, format="csv", config=ReadConfig(delimiter=";"))

    assert restored.ids == ("sequence_1",)


def test_file_object_name_supports_format_inference(tmp_path: Path) -> None:
    path = tmp_path / "one.fasta"
    path.write_text(">one\nACGT\n", encoding="utf-8")

    with path.open(encoding="utf-8") as handle:
        restored = read_set(handle)
        assert not handle.closed

    assert restored.ids == ("one",)


def test_create_parents_is_explicit(tmp_path: Path) -> None:
    path = tmp_path / "nested" / "one.fasta"

    with pytest.raises(FileNotFoundError):
        write(DNASequence("A"), path)

    write(DNASequence("A"), path, config=WriteConfig(create_parents=True))
    assert path.read_text(encoding="utf-8") == ">sequence_1\nA\n"


def test_text_stream_compression_conflict_honors_explicit_ownership() -> None:
    output = io.StringIO()
    with pytest.raises(ConfigurationError):
        write(
            DNASequence("A"),
            output,
            format="fasta",
            config=WriteConfig(compression="gzip", close_target=True),
        )
    assert output.closed

    source = io.StringIO(">one\nA\n")
    with pytest.raises(ConfigurationError):
        read_set(
            source,
            format="fasta",
            config=ReadConfig(compression="gzip", close_source=True),
        )
    assert source.closed


def test_write_result_is_json_compatible() -> None:
    result = write(
        DNASequence("A"),
        io.StringIO(),
        format="fasta",
        config=WriteConfig(line_width=7, feature_policy="drop"),
    )

    payload = result.to_dict()

    assert payload["format"] == "fasta"
    assert result.provenance.implementation.label.value == "reimplementation"
    assert payload["record_count"] == 1
    assert payload["generated_ids"] == [{"input_index": 0, "generated_id": "sequence_1"}]
    assert payload["target_artifact"] is None
    assert payload["parameters"]["line_width"] == 7
    assert payload["parameters"]["feature_policy"] == "drop"
    assert payload["parameters"]["compression"] == "auto"
    assert payload["provenance"]["implementation"]["label"] == "reimplementation"


def test_write_result_rejects_invalid_provenance_type() -> None:
    with pytest.raises(ConfigurationError):
        WriteResult("fasta", 1, provenance=object())  # type: ignore[arg-type]


def test_read_record_limit_is_enforced_lazily() -> None:
    stream = io.StringIO(">one\nA\n>two\nC\n")

    with pytest.raises(InputFormatError) as exc_info:
        read_set(stream, format="fasta", config=ReadConfig(max_records=1))

    assert getattr(exc_info.value, "code", None) == "INPUT_RECORD_LIMIT_EXCEEDED"
