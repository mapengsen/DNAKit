"""Tests for dnakit.fasta_fastq.cutadapt_gap_functions (cutadapt command-line facades)."""

from __future__ import annotations

import inspect
import shutil
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import cutadapt_gap_functions as facade

_LOCAL_CUTADAPT = Path(sys.executable).parent / "cutadapt"
_CUTADAPT = shutil.which("cutadapt") or (str(_LOCAL_CUTADAPT) if _LOCAL_CUTADAPT.is_file() else None)
needs_cutadapt = pytest.mark.skipif(_CUTADAPT is None, reason="cutadapt executable is not installed")

_READS = [
    ("r1", "ACGTACGTAGATCGGAAGAGCTTTTTTTT"),
    ("r2", "GGGGCCCCAAAATTTTGGGGCCCC"),
]


def _write_fastq(path: Path) -> Path:
    path.write_text("".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in _READS), encoding="utf-8")
    return path


def _fake_cutadapt(directory: Path) -> Path:
    """A stand-in executable that records its argument list and writes a one-read FASTQ output."""
    script = directory / "fake_cutadapt"
    script.write_text(
        f"#!{sys.executable}\n"
        "import pathlib, sys\n"
        "args = sys.argv[1:]\n"
        "if args == ['--version']:\n"
        "    print('9.9-fake')\n"
        "    raise SystemExit(0)\n"
        "out = args[args.index('-o') + 1]\n"
        "pathlib.Path(out).write_text('@r\\nA\\n+\\nI\\n', encoding='utf-8')\n"
        "print('fake report')\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def test_public_functions_take_keyword_only_verbose() -> None:
    names = [n for n in facade.__all__]
    assert names
    for name in names:
        signature = inspect.signature(getattr(facade, name))
        verbose = signature.parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True


def test_invalid_arguments_are_rejected_before_running(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(tmp_path / "a.fq"), adapter_3p="ACGTX!", cutadapt_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(tmp_path / "b.fq"), error_rate=1.5, cutadapt_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(tmp_path / "c.fq"), quality_cutoff="abc", cutadapt_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(
            str(source), str(tmp_path / "d.fq"), minimum_length=30, maximum_length=10, cutadapt_path="/nope"
        )
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(tmp_path / "e.fq"), cores=0, cutadapt_path="/nope")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(tmp_path / "f.fq"), report="long", cutadapt_path="/nope")
    existing = tmp_path / "exists.fq"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim(str(source), str(existing), cutadapt_path="/nope")


def test_missing_executable_raises_with_install_hint(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    with pytest.raises(BackendUnavailableError) as info:
        facade.cutadapt_trim(str(source), str(tmp_path / "out.fq"), cutadapt_path=str(tmp_path / "missing"))
    assert "conda install -c bioconda cutadapt" in (info.value.hint or "")


def test_argument_list_maps_fixed_flags(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    fake = _fake_cutadapt(tmp_path)
    result = facade.cutadapt_trim(
        str(source),
        str(tmp_path / "out.fq"),
        adapter_3p="agatcggaagagc",
        quality_cutoff="10,20",
        minimum_length=20,
        maximum_length=200,
        discard_untrimmed=True,
        cores=2,
        cutadapt_path=str(fake),
    )
    assert result["status"] == "completed"
    assert result["backend"] == "cutadapt"
    assert result["version"] == "9.9-fake"
    assert result["report"].strip() == "fake report"
    command = result["command"]
    assert command[0] == "cutadapt"
    assert command[command.index("-a") + 1] == "AGATCGGAAGAGC"
    assert command[command.index("-q") + 1] == "10,20"
    assert command[command.index("-m") + 1] == "20"
    assert command[command.index("-M") + 1] == "200"
    assert "--discard-untrimmed" in command
    assert command[command.index("-j") + 1] == "2"
    assert command[command.index("--report") + 1] == "full"
    assert command[-1] == str(source.resolve())
    assert result["outputs"] == [str((tmp_path / "out.fq").resolve())]


@needs_cutadapt
def test_trim_removes_3p_adapter(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    output = tmp_path / "trimmed.fq"
    result = facade.cutadapt_trim(str(source), str(output), adapter_3p="AGATCGGAAGAGC", cutadapt_path=_CUTADAPT)
    assert result["status"] == "completed"
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[1] == "ACGTACGT"
    assert lines[5] == _READS[1][1]
    assert "Total reads processed" in result["report"]


@needs_cutadapt
def test_full_and_minimal_reports(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    full = facade.cutadapt_full_report(str(source), str(tmp_path / "full.fq"), adapter_3p="AGATCGGAAGAGC",
                                       cutadapt_path=_CUTADAPT)
    assert full["status"] == "completed"
    assert "=== Summary ===" in full["report"]
    minimal = facade.cutadapt_minimal_report(str(source), str(tmp_path / "min.fq"), adapter_3p="AGATCGGAAGAGC",
                                             cutadapt_path=_CUTADAPT)
    assert minimal["status"] == "completed"
    header, row = minimal["report"].strip().splitlines()[:2]
    assert header.startswith("status\tin_reads")
    assert row.split("\t")[1] == "2"


@needs_cutadapt
def test_json_report_is_written_and_parsed(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "in.fq")
    result = facade.cutadapt_json_report(
        str(source),
        str(tmp_path / "out.fq"),
        str(tmp_path / "report.json"),
        adapter_3p="AGATCGGAAGAGC",
        cutadapt_path=_CUTADAPT,
    )
    assert result["status"] == "completed"
    assert result["json"]["tag"] == "Cutadapt report"
    assert (tmp_path / "report.json").is_file()
