"""Contract tests for the bounded Trimmomatic facade."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc import run_afterqc, trim_trimmomatic
from dnakit.tools import default_tool_registry


def _fake_trimmomatic(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "arguments = sys.argv\n"
        "mode = arguments[1]\n"
        "if mode == 'SE':\n"
        "    source = Path(arguments[arguments.index('-phred33') + 1])\n"
        "    output = Path(arguments[arguments.index('-phred33') + 2])\n"
        "    output.write_text(source.read_text(encoding='ascii'), encoding='ascii')\n"
        "else:\n"
        "    source1 = Path(arguments[arguments.index('-phred33') + 1])\n"
        "    source2 = Path(arguments[arguments.index('-phred33') + 2])\n"
        "    start = arguments.index('-phred33') + 3\n"
        "    for source, output in ((source1, arguments[start]), (source2, arguments[start + 2])):\n"
        "        Path(output).write_text(source.read_text(encoding='ascii'), encoding='ascii')\n"
        "    Path(arguments[start + 1]).write_text('', encoding='ascii')\n"
        "    Path(arguments[start + 3]).write_text('', encoding='ascii')\n"
        "print('trimmomatic diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_trim_trimmomatic_single_end(tmp_path: Path) -> None:
    read1 = tmp_path / "reads.fastq"
    adapter = tmp_path / "adapters.fa"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    adapter.write_text(">adapter\nAGAT\n", encoding="ascii")
    executable = _fake_trimmomatic(tmp_path / "fake-trimmomatic")

    result = trim_trimmomatic(
        read1,
        tmp_path / "single-output",
        adapter_fasta=adapter,
        min_length=4,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "trimmomatic"
    assert "trimmed.fastq.gz" in result["result"]["artifacts"]


def test_trim_trimmomatic_paired_end(tmp_path: Path) -> None:
    read1 = tmp_path / "reads_R1.fastq"
    read2 = tmp_path / "reads_R2.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    read2.write_text("@r1\nTGCA\n+\nIIII\n", encoding="ascii")
    executable = _fake_trimmomatic(tmp_path / "fake-trimmomatic")

    result = trim_trimmomatic(
        read1,
        tmp_path / "paired-output",
        read2=read2,
        leading=3,
        trailing=3,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "trimmed_R1.paired.fastq.gz" in result["result"]["artifacts"]
    assert "trimmed_R2.paired.fastq.gz" in result["result"]["artifacts"]


def test_trim_trimmomatic_requires_a_step(tmp_path: Path) -> None:
    read1 = tmp_path / "reads.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    with pytest.raises(ConfigurationError, match="trimming step"):
        trim_trimmomatic(read1, tmp_path / "invalid", verbose=False)


def test_trim_trimmomatic_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.sequencing_qc.trim_trimmomatic")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True


def _fake_afterqc(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "arguments = sys.argv[1:]\n"
        "def value(flag):\n"
        "    return Path(arguments[arguments.index(flag) + 1])\n"
        "value('-g').mkdir()\n"
        "value('-b').mkdir()\n"
        "Path('QC').mkdir()\n"
        "Path('QC/report.html').write_text('<html>ok</html>', encoding='utf-8')\n"
        "print('afterqc diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_run_afterqc_single_end(tmp_path: Path) -> None:
    read1 = tmp_path / "reads_R1.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_afterqc(tmp_path / "after.py")

    result = run_afterqc(
        read1,
        tmp_path / "afterqc-output",
        qualified_quality_phred=20,
        seq_len_req=4,
        gzip_output=True,
        compression=5,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "afterqc"
    assert "QC/report.html" in result["result"]["artifacts"]


def test_run_afterqc_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.sequencing_qc.run_afterqc")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True
