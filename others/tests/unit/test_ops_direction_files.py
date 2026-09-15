"""Text and streaming file contracts for OPS-001.3."""

from __future__ import annotations

import csv
import gzip
import io
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import ConfigurationError, InputFormatError, InvalidAlphabetError
from dnakit.io import ReadConfig, WriteConfig, WriteResult, read, read_set, write
from dnakit.ops import reverse_complement


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [("AATCGC", "GCGATT"), ("aa CGt\nN", "NACGTT"), ("ARYKMBVDHSWN", "NWSDHBVKMRYT")],
)
def test_direct_dna_text(sequence: str, expected: str) -> None:
    assert reverse_complement(sequence) == expected


def test_dna_text_does_not_guess_an_existing_filename(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.chdir(tmp_path)
    Path("ACGT").write_text("different contents")
    assert reverse_complement("ACGT") == "ACGT"
    with pytest.raises(InvalidAlphabetError):
        reverse_complement("ACGX")


@pytest.mark.parametrize("suffix", [".fna", ".fa", ".fasta", ".FNA", ".fna.gz"])
def test_file_inference_and_default_output(tmp_path: Path, suffix: str) -> None:
    source = tmp_path / ("assembly.v1" + suffix)
    text = ">chr1 first\nAaCGTN\n>contig\nARY\n"
    if suffix.endswith(".gz"):
        with gzip.open(source, "wt") as handle:
            handle.write(text)
    else:
        source.write_text(text)
    original = source.read_bytes()
    result = reverse_complement(str(source), show_progress=False)
    assert isinstance(result, WriteResult)
    assert result.record_count == 2
    assert result.target_artifact is not None
    output = Path(result.target_artifact.relative_path)
    assert output.name == "assembly.v1_reverse_complement" + suffix
    with read(output) as records:
        first, second = records
    assert (first.id, first.description, first.sequence.to_string()) == ("chr1", "first", "NACGTT")
    assert (second.id, second.sequence.to_string()) == ("contig", "RYT")
    assert source.read_bytes() == original


@pytest.mark.parametrize(
    "text_mode,expected", [("lines", "CGTT\nAAGC\n"), ("sequence", "AAGCCGTT\n")]
)
def test_txt_records_and_wrapped_sequence(tmp_path: Path, text_mode: Any, expected: str) -> None:
    source = tmp_path / "seqs.txt"
    source.write_text("aaCG\n\nGCTT\n")
    output = tmp_path / "rc.txt"
    reverse_complement(
        source, output, read_config=ReadConfig(text_mode=text_mode), show_progress=False
    )
    assert output.read_text() == expected


def test_csv_named_column_and_metadata(tmp_path: Path) -> None:
    source = tmp_path / "seqs.csv"
    with source.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "dna", "metadata"])
        writer.writerow(["sample1", "AAGC", '{"label": 1}'])
    output = tmp_path / "rc.csv"
    result = reverse_complement(source, output, sequence_column="dna", show_progress=False)
    assert result.record_count == 1
    with read(output, config=ReadConfig(sequence_column="dna")) as records:
        record = next(records)
    assert record.sequence.to_string() == "GCTT"
    assert record.id == "sample1"
    assert record.metadata["label"] == 1


def test_csv_missing_sequence_is_rejected_without_output(tmp_path: Path) -> None:
    source = tmp_path / "seqs.csv"
    source.write_text("id,label\na,1\n")
    output = tmp_path / "rc.csv"
    with pytest.raises(InputFormatError, match="sequence"):
        reverse_complement(source, output, show_progress=False)
    assert not output.exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_large_local_file_limits_are_derived_without_reading_whole_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import dnakit._processing.inputs as file_direction

    source = tmp_path / "genome.fna"
    with source.open("wb") as handle:
        handle.truncate(3_200_000_000)  # Sparse fixture: no genome allocation or full read.

    @contextmanager
    def fake_read(path: Path, *, config: ReadConfig) -> Any:
        assert path == source
        assert config.max_input_bytes >= 3_200_000_000
        assert config.max_sequence_symbols >= 248_956_422
        assert config.max_record_lines > 3_000_000
        yield iter([DNARecord(DNASequence("AAGC"), "chr1")])

    monkeypatch.setattr(file_direction, "read", fake_read)
    result = reverse_complement(source, show_progress=False)
    assert result.record_count == 1


def test_explicit_limits_still_apply(tmp_path: Path) -> None:
    source = tmp_path / "seq.fna"
    source.write_text(">chr1\nAAGCT\n")
    with pytest.raises(InputFormatError):
        reverse_complement(
            source, read_config=ReadConfig(max_sequence_symbols=4), show_progress=False
        )
    assert not (tmp_path / "seq_reverse_complement.fna").exists()


def test_no_input_overwrite_and_explicit_output_replacement(tmp_path: Path) -> None:
    source = tmp_path / "seq.txt"
    source.write_text("AAGC\n")
    with pytest.raises(ConfigurationError, match="differ"):
        reverse_complement(source, source, write_config=WriteConfig(overwrite=True))
    output = tmp_path / "rc.txt"
    output.write_text("old\n")
    with pytest.raises(FileExistsError):
        reverse_complement(source, output, show_progress=False)
    assert output.read_text() == "old\n"
    reverse_complement(
        source, output, write_config=WriteConfig(overwrite=True), show_progress=False
    )
    assert output.read_text() == "GCTT\n"
    assert source.read_text() == "AAGC\n"


def test_malformed_later_record_preserves_existing_output(tmp_path: Path) -> None:
    source = tmp_path / "seq.txt"
    source.write_text("AAGC\nINVALID\n")
    output = tmp_path / "rc.txt"
    output.write_text("old\n")
    with pytest.raises(InputFormatError):
        reverse_complement(
            source, output, write_config=WriteConfig(overwrite=True), show_progress=False
        )
    assert output.read_text() == "old\n"


@pytest.mark.parametrize("show_progress", [True, False])
def test_file_console_messages(tmp_path: Path, capsys: Any, show_progress: bool) -> None:
    source = tmp_path / "seq.txt"
    source.write_text("AAGC\n")
    reverse_complement(source, show_progress=show_progress)
    output = capsys.readouterr()
    assert output.out == ""
    if show_progress:
        assert "Starting reverse complement:" in output.err
        assert "Reverse complement complete:" in output.err
        assert "1 records" in output.err
    else:
        assert output.err == ""


def test_missing_and_unsupported_files_are_not_dna_text(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        reverse_complement(str(tmp_path / "missing.fna"))
    source = tmp_path / "seq.unknown"
    source.write_text("ACGT")
    with pytest.raises(InputFormatError, match="format"):
        reverse_complement(source, show_progress=False)


def test_plain_text_io_and_limits() -> None:
    records = read_set(io.StringIO("AAGC\n\nNRY\n"), format="txt")
    assert records.ids == ("sequence_1", "sequence_2")
    target = io.StringIO()
    write(records, target, format="txt")
    assert target.getvalue() == "AAGC\nNRY\n"
    with pytest.raises(InputFormatError):
        read_set(io.StringIO("AAGC\n"), format="txt", config=ReadConfig(max_sequence_symbols=3))
    with pytest.raises(InputFormatError):
        read_set(
            io.StringIO("AA\n\n"),
            format="txt",
            config=ReadConfig(text_mode="sequence", max_record_lines=1),
        )


def test_sequence_input_rejects_file_only_options() -> None:
    with pytest.raises(ConfigurationError, match="File options"):
        reverse_complement("ACGT", "out.fna")
