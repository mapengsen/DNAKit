"""Tests for dnakit.fasta_fastq.screed_gap_functions (screed index build and FASTA/FASTQ dumps).

Inputs are small files written to tmp_path; nothing is downloaded.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pytest.importorskip("screed")

from dnakit.fasta_fastq import screed_gap_functions as facades  # noqa: E402

FASTA = ">s1 desc\nACGTACGT\n>s2\nGGGCCCAAA\n"
FASTQ = "@r1 x\nACGTNN\n+\nIIIIII\n"


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    path = tmp_path / "seqs.fa"
    path.write_text(FASTA, encoding="utf-8")
    return path


@pytest.fixture()
def fastq_file(tmp_path: Path) -> Path:
    path = tmp_path / "reads.fq"
    path.write_text(FASTQ, encoding="utf-8")
    return path


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_create_database_from_fasta(fasta_file, tmp_path):
    result = facades.create_screed_database(str(fasta_file), str(tmp_path / "seqs"), verbose=False)
    assert result["input_format"] == "FASTA"
    assert result["record_count"] == 2
    assert result["database"].endswith("_screed")
    assert Path(result["database"]).is_file()


def test_create_database_from_fastq(fastq_file, tmp_path):
    result = facades.create_screed_database(str(fastq_file), str(tmp_path / "reads.db"), verbose=False)
    assert result["input_format"] == "FASTQ"
    assert result["record_count"] == 1


def test_dump_fasta_round_trips_sequences(fasta_file, tmp_path):
    database = facades.create_screed_database(str(fasta_file), str(tmp_path / "seqs"), verbose=False)["database"]
    out = tmp_path / "back.fa"
    result = facades.dump_screed_fasta(database, str(out), verbose=False)
    assert result["record_count"] == 2
    text = out.read_text(encoding="utf-8")
    assert "ACGTACGT" in text and "GGGCCCAAA" in text
    assert text.count(">") == 2


def test_dump_fastq_from_fastq_database(fastq_file, tmp_path):
    database = facades.create_screed_database(str(fastq_file), str(tmp_path / "reads"), verbose=False)["database"]
    out = tmp_path / "back.fq"
    result = facades.dump_screed_fastq(database, str(out), verbose=False)
    assert result["record_count"] == 1
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("@r1")
    assert lines[1] == "ACGTNN"
    assert lines[3] == "IIIIII"


def test_dump_fastq_needs_quality_values(fasta_file, tmp_path):
    database = facades.create_screed_database(str(fasta_file), str(tmp_path / "seqs"), verbose=False)["database"]
    with pytest.raises(ConfigurationError):
        facades.dump_screed_fastq(database, str(tmp_path / "out.fq"), verbose=False)


def test_dump_refuses_empty_database(tmp_path):
    from screed import createscreed, fasta

    empty = str(tmp_path / "empty_screed")
    createscreed.create_db(empty, fasta.FieldTypes, iter([]))
    with pytest.raises(ConfigurationError):
        facades.dump_screed_fasta(empty, str(tmp_path / "out.fa"), verbose=False)


def test_existing_outputs_are_refused(fasta_file, tmp_path):
    existing = tmp_path / "taken"
    existing.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        facades.create_screed_database(str(fasta_file), str(existing), verbose=False)
    database = facades.create_screed_database(str(fasta_file), str(tmp_path / "seqs"), verbose=False)["database"]
    with pytest.raises(ConfigurationError):
        facades.dump_screed_fasta(database, str(existing), verbose=False)


def test_missing_input_is_reported(tmp_path):
    with pytest.raises(ConfigurationError):
        facades.create_screed_database(str(tmp_path / "none.fa"), str(tmp_path / "x"), verbose=False)
    with pytest.raises(ConfigurationError):
        facades.dump_screed_fasta(str(tmp_path / "none_screed"), str(tmp_path / "x.fa"), verbose=False)


def test_missing_screed_gives_install_hint(monkeypatch, fasta_file, tmp_path):
    import importlib

    real = importlib.import_module

    def blocked(name, package=None):
        if name.startswith("screed"):
            raise ImportError("blocked for test")
        return real(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", blocked)
    with pytest.raises(BackendUnavailableError) as info:
        facades.create_screed_database(str(fasta_file), str(tmp_path / "x"), verbose=False)
    assert "bio-screed" in (info.value.hint or "")
