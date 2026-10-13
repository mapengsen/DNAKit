"""Contract tests for the bounded Cutadapt and fastp adapters."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.sequencing_qc import run_cutadapt, run_fastp, run_fastqc, run_multiqc


def _backend(tmp_path: Path, *, name: str) -> Path:
    script = tmp_path / name
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "def value(flag):\n"
        "    return args[args.index(flag) + 1]\n"
        "for flag in ('-o', '-O', '-p', '--json', '-j', '-h'):\n"
        "    if flag in args:\n"
        "        path = Path(value(flag))\n"
        "        path.parent.mkdir(parents=True, exist_ok=True)\n"
        "        content = 'synthetic backend output\\n'\n"
        "        if path.suffix == '.json':\n"
        "            content = '{\\\"ok\\\": true}\\n'\n"
        "        path.write_text(content)\n"
        "print('synthetic backend completed')\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def _report_backend(tmp_path: Path, *, name: str) -> Path:
    script = tmp_path / name
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "outdir = Path(args[args.index('--outdir') + 1])\n"
        "outdir.joinpath('sample_fastqc.html').write_text('html\\n')\n"
        "outdir.joinpath('sample_fastqc.zip').write_text('zip\\n')\n"
        "if '--filename' in args:\n"
        "    outdir.joinpath(args[args.index('--filename') + 1]).write_text('multiqc\\n')\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def test_cutadapt_adapter_executes_paired_reads_and_report(tmp_path: Path) -> None:
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    read2.write_text("@r1\nTGCA\n+\nIIII\n", encoding="utf-8")
    output1 = tmp_path / "trimmed.R1.fastq"
    output2 = tmp_path / "trimmed.R2.fastq"
    report = tmp_path / "cutadapt.json"

    result = run_cutadapt(
        {"read1": read1, "read2": read2},
        options={
            "executable_path": _backend(tmp_path, name="cutadapt"),
            "output1": output1,
            "output2": output2,
            "adapters": ["AGATCGGAAG"],
            "json_report": report,
            "cores": 2,
        },
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["citation_doi"] == "10.14806/ej.17.1.200"
    assert output1.is_file() and output2.is_file() and report.is_file()
    assert "-p" in result["arguments"]


def test_fastp_adapter_executes_single_end_with_report(tmp_path: Path) -> None:
    read1 = tmp_path / "reads.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    output = tmp_path / "filtered.fastq"
    report = tmp_path / "fastp.json"

    result = run_fastp(
        read1,
        options={
            "executable_path": _backend(tmp_path, name="fastp"),
            "output": output,
            "json_report": report,
            "length_required": 3,
            "qualified_quality_phred": 20,
        },
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["citation_doi"] == "10.1002/imt2.70078"
    assert output.is_file() and report.is_file()
    assert "-l" in result["arguments"]


def test_external_adapter_rejects_unallowlisted_options(tmp_path: Path) -> None:
    read1 = tmp_path / "reads.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Unsupported fastp options"):
        run_fastp(
            read1,
            options={"output": tmp_path / "out.fastq", "arbitrary_flag": "--unsafe"},
            verbose=False,
        )


def test_fastqc_and_multiqc_adapters_monitor_reports(tmp_path: Path) -> None:
    reads = tmp_path / "sample.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    output_directory = tmp_path / "reports"
    output_directory.mkdir()
    fastqc = run_fastqc(
        reads,
        options={
            "executable_path": _report_backend(tmp_path, name="fastqc"),
            "output_directory": output_directory,
        },
        verbose=False,
    )
    multiqc_input = tmp_path / "multiqc-input"
    multiqc_input.mkdir()
    (multiqc_input / "fastqc.html").write_text("html\n", encoding="utf-8")
    multiqc_directory = tmp_path / "multiqc-reports"
    multiqc_directory.mkdir()
    multiqc = run_multiqc(
        multiqc_input,
        options={
            "executable_path": _report_backend(tmp_path, name="multiqc"),
            "output_directory": multiqc_directory,
        },
        verbose=False,
    )
    assert fastqc["status"] == "completed"
    assert multiqc["status"] == "completed"
    assert fastqc["license"] == "GPL-3.0"
    assert multiqc["citation_doi"] == "10.1093/bioinformatics/btw354"


def test_external_adapter_rejects_success_without_declared_artifacts(tmp_path: Path) -> None:
    read1 = tmp_path / "reads.fastq"
    read1.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    backend = tmp_path / "silent-fastp"
    backend.write_text(f"#!{sys.executable}\n", encoding="utf-8")
    backend.chmod(0o755)

    with pytest.raises(BackendExecutionError, match="did not produce"):
        run_fastp(
            read1,
            options={"executable_path": backend, "output": tmp_path / "out.fastq"},
            verbose=False,
        )


def test_fastqc_passes_kmer_option_and_uses_expected_report_names(tmp_path: Path) -> None:
    reads = tmp_path / "sample.fastq.gz"
    reads.write_text("synthetic\n", encoding="utf-8")
    output_directory = tmp_path / "reports"
    output_directory.mkdir()
    backend = tmp_path / "fastqc"
    backend.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "assert '--kmers' in args and args[args.index('--kmers') + 1] == '7'\n"
        "outdir = Path(args[args.index('--outdir') + 1])\n"
        "outdir.joinpath('sample_fastqc.html').write_text('html\\n')\n"
        "outdir.joinpath('sample_fastqc.zip').write_text('zip\\n')\n",
        encoding="utf-8",
    )
    backend.chmod(0o755)

    result = run_fastqc(
        reads,
        options={
            "executable_path": backend,
            "output_directory": output_directory,
            "kmer": 7,
        },
        verbose=False,
    )

    assert result["status"] == "completed"
    assert {Path(path).name for path in result["outputs"]} == {
        "sample_fastqc.html",
        "sample_fastqc.zip",
    }
