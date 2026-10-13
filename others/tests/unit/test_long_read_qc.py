"""Contract tests for long-read NanoFilt, Filtlong and NanoPlot facades."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc import (
    filter_long_reads,
    filter_nanopore_reads,
    summarize_long_reads,
    trim_porechop,
)
from dnakit.tools import default_tool_registry


def _fake_filter(path: Path, label: str) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "source = Path(sys.argv[-1])\n"
        "print(source.read_text(encoding='ascii'), end='')\n"
        f"print('{label} diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_nanoplot(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "directory = Path(sys.argv[sys.argv.index('--outdir') + 1])\n"
        "(directory / 'NanoStats.txt').write_text('N50\t4\\n', encoding='ascii')\n"
        "(directory / 'NanoPlot-report.html').write_text('<html></html>', encoding='ascii')\n"
        "print('nanoplot diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_porechop(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "arguments = sys.argv\n"
        "source = Path(arguments[arguments.index('-i') + 1])\n"
        "output = Path(arguments[arguments.index('-o') + 1])\n"
        "output.write_text(source.read_text(encoding='ascii'), encoding='ascii')\n"
        "print('porechop diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_filter_nanopore_reads_captures_fastq_and_diagnostics(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_filter(tmp_path / "fake-nanofilt", "nanofilt")

    result = filter_nanopore_reads(
        source,
        tmp_path / "nanofilt-output",
        min_quality=10,
        min_length=4,
        headcrop=1,
        tailcrop=1,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "nanofilt"
    assert "filtered.fastq" in result["result"]["artifacts"]
    assert "filtered.fastq.stderr" in result["result"]["artifacts"]


def test_filter_long_reads_captures_fastq(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_filter(tmp_path / "fake-filtlong", "filtlong")

    result = filter_long_reads(
        source,
        tmp_path / "filtlong-output",
        min_length=4,
        keep_percent=90,
        target_bases=4,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "filtlong"
    assert "filtered.fastq" in result["result"]["artifacts"]


def test_summarize_long_reads_checks_nanoplot_reports(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_nanoplot(tmp_path / "fake-nanoplot")

    result = summarize_long_reads(
        source,
        tmp_path / "nanoplot-output",
        threads=2,
        min_length=4,
        min_quality=10,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "nanoplot"
    assert "NanoStats.txt" in result["result"]["artifacts"]
    assert "NanoPlot-report.html" in result["result"]["artifacts"]


def test_trim_porechop_writes_trimmed_fastq(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_porechop(tmp_path / "fake-porechop")

    result = trim_porechop(
        source,
        tmp_path / "porechop-output",
        threads=2,
        no_split=True,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "porechop"
    assert "trimmed.fastq" in result["result"]["artifacts"]


def test_long_read_facades_validate_numeric_limits(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    source.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    with pytest.raises(ConfigurationError, match="keep_percent"):
        filter_long_reads(source, tmp_path / "invalid", keep_percent=101, verbose=False)
    with pytest.raises(ConfigurationError, match="threads"):
        summarize_long_reads(source, tmp_path / "invalid-plot", threads=0, verbose=False)


def test_long_read_facades_are_registered_as_external_tools() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.sequencing_qc.filter_nanopore_reads",
        "dnakit.sequencing_qc.filter_long_reads",
        "dnakit.sequencing_qc.summarize_long_reads",
        "dnakit.sequencing_qc.trim_porechop",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True
