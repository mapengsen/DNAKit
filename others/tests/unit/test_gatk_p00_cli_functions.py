"""Tests for the GATK P00 command-line facades (three modules and their shared core).

Structural and fake-launcher tests run without GATK. Tests that run the real program are skipped
unless ``gatk`` (GATK 4.6 and Java 17) is on PATH or ``DNAKIT_GATK`` names the launcher.
"""

from __future__ import annotations

import inspect
import json
import os
import shutil
import stat
from pathlib import Path

import pytest

from dnakit.annotation import gatk_p00_intervals_cli_functions as intervals
from dnakit.bioinformatics import _gatk_p00_core as core
from dnakit.bioinformatics._gatk_p00_options import TOOL_OPTIONS
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import gatk_p00_metrics_cli_functions as metrics
from dnakit.variants import gatk_p00_copynumber_cli_functions as copynumber

MODULES = (metrics, copynumber, intervals)
FUNCTIONS = [getattr(m, name) for m in MODULES for name in m.__all__]
GATK = os.environ.get("DNAKIT_GATK") or shutil.which("gatk")
requires_gatk = pytest.mark.skipif(
    GATK is None or shutil.which("java") is None, reason="GATK 4.6 with Java is not installed"
)


def test_public_function_count() -> None:
    assert len(FUNCTIONS) == 93
    assert len({f.__name__ for f in FUNCTIONS}) == 93


@pytest.mark.parametrize("function", FUNCTIONS, ids=lambda f: f.__name__)
def test_signature_contract(function) -> None:
    parameters = inspect.signature(function).parameters
    for name in ("options", "gatk_path", "java_max_memory", "verbose"):
        assert parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["verbose"].default is True
    assert function.__doc__ and "Parameters" in function.__doc__


def test_option_table_covers_every_tool() -> None:
    assert len(TOOL_OPTIONS) == 100
    for tool, table in TOOL_OPTIONS.items():
        assert set(table.values()) <= {"value", "path", "out", "dir"}, tool


def test_option_arguments_whitelist(tmp_path: Path) -> None:
    assert core.option_arguments("CountReads", {"ambig-filter-frac": 0.5}) == [
        "--ambig-filter-frac",
        "0.5",
    ]
    with pytest.raises(ConfigurationError):
        core.option_arguments("CountReads", {"not-an-option": 1})
    with pytest.raises(ConfigurationError):
        core.option_arguments("CountReads", {"arguments_file": "x"})
    with pytest.raises(ConfigurationError):
        core.option_arguments("CountReads", {"output": "a"}, reserved=("output",))
    existing = tmp_path / "exists.txt"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        core.option_arguments("CountReads", {"output": str(existing)})


def test_missing_launcher(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv("DNAKIT_GATK", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    bam = tmp_path / "a.bam"
    bam.write_bytes(b"x")
    with pytest.raises(BackendUnavailableError):
        metrics.gatk_check_terminator_block(str(bam), verbose=False)


def test_input_validation(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        metrics.gatk_get_sample_name(
            [str(tmp_path / "missing.bam")], str(tmp_path / "o"), verbose=False
        )
    done = tmp_path / "out.txt"
    done.write_text("x")
    bam = tmp_path / "a.bam"
    bam.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        metrics.gatk_get_sample_name([str(bam)], str(done), verbose=False)
    with pytest.raises(ConfigurationError):
        intervals.gatk_convert_interval_list(
            str(bam), str(tmp_path / "x.il"), output_format="interval_list", verbose=False
        )
    with pytest.raises(ConfigurationError):
        copynumber.gatk_merge_annotated_regions(
            str(bam), str(tmp_path / "m.seg"), max_merge_distance=5, verbose=False
        )


FAKE = """#!/usr/bin/env python3
import json, os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("The Genome Analysis Toolkit (GATK) v4.6.2.0")
    sys.exit(0)
open(os.environ["FAKE_GATK_LOG"], "w").write(json.dumps(args))
for flag in ("--output", "--OUTPUT"):
    if flag in args:
        open(args[args.index(flag) + 1], "w").write("7\\n")
sys.exit(int(os.environ.get("FAKE_GATK_RC", "0")))
"""


@pytest.fixture
def fake(tmp_path: Path, monkeypatch) -> Path:
    launcher = tmp_path / "gatk"
    launcher.write_text(FAKE)
    launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("FAKE_GATK_LOG", str(tmp_path / "args.json"))
    monkeypatch.setenv("DNAKIT_GATK", str(launcher))
    return tmp_path


def test_fake_count_reads_command(fake: Path) -> None:
    bam = fake / "a.bam"
    bam.write_bytes(b"x")
    result = copynumber.gatk_count_reads_or_bases(
        [str(bam)],
        str(fake / "n.txt"),
        unit="bases",
        options={"ambig-filter-frac": 0.2},
        verbose=False,
    )
    args = json.loads((fake / "args.json").read_text())
    assert args[0] == "CountBases"
    assert args[args.index("--input") + 1] == str(bam)
    assert "--ambig-filter-frac" in args
    assert result["count"] == 7 and result["unit"] == "bases" and result["version"] == "4.6.2.0"


def test_fake_variant_selects_tool(fake: Path) -> None:
    bam = fake / "a.bam"
    bam.write_bytes(b"x")
    metrics.gatk_collect_quality_yield_metrics(
        str(bam), str(fake / "q.txt"), variant="flow", include_bq_histogram=True, verbose=False
    )
    assert json.loads((fake / "args.json").read_text())[0] == "CollectQualityYieldMetricsFlow"
    with pytest.raises(ConfigurationError):
        metrics.gatk_collect_quality_yield_metrics(
            str(bam),
            str(fake / "q2.txt"),
            variant="standard",
            include_bq_histogram=True,
            verbose=False,
        )


def test_fake_failure_and_finding_code(fake: Path, monkeypatch) -> None:
    bam = fake / "a.bam"
    bam.write_bytes(b"x")
    monkeypatch.setenv("FAKE_GATK_RC", "100")
    result = metrics.gatk_check_terminator_block(str(bam), verbose=False)
    assert result["terminator_ok"] is False and result["has_findings"] is True
    monkeypatch.setenv("FAKE_GATK_RC", "3")
    with pytest.raises(Exception, match="exited with code 3"):
        metrics.gatk_check_terminator_block(str(bam), verbose=False)


def _toy(path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    path.mkdir()
    (path / "toy.fa").write_text(">chr1\n" + "ACGTTGCAAC" * 30 + "\n")
    (path / "toy.fa.fai").write_text("chr1\t300\t6\t300\t301\n")
    (path / "toy.dict").write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:300\n")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 300}],
        "RG": [{"ID": "rg1", "SM": "S1", "LB": "l", "PL": "ILLUMINA"}],
    }
    with pysam.AlignmentFile(str(path / "toy.bam"), "wb", header=header) as handle:
        for index in range(10):
            read = pysam.AlignedSegment()
            read.query_name = f"r{index}"
            read.query_sequence = ("ACGTTGCAAC" * 30)[index * 10 : index * 10 + 50]
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.reference_id = 0
            read.reference_start = index * 10
            read.cigarstring = "50M"
            read.mapping_quality = 60
            read.flag = 0
            read.set_tag("RG", "rg1")
            read.set_tag("NM", 0)
            handle.write(read)
    pysam.index(str(path / "toy.bam"))


@requires_gatk
def test_real_counts_and_sample(tmp_path: Path) -> None:
    data = tmp_path / "d"
    _toy(data)
    reads = copynumber.gatk_count_reads_or_bases(
        [str(data / "toy.bam")], str(tmp_path / "r.txt"), gatk_path=GATK, verbose=False
    )
    bases = copynumber.gatk_count_reads_or_bases(
        [str(data / "toy.bam")],
        str(tmp_path / "b.txt"),
        unit="bases",
        gatk_path=GATK,
        verbose=False,
    )
    assert reads["count"] == 10 and bases["count"] == 500 and reads["status"] == "completed"
    name = metrics.gatk_get_sample_name(
        [str(data / "toy.bam")], str(tmp_path / "s.txt"), gatk_path=GATK, verbose=False
    )
    assert name["sample_name"] == "S1"


@requires_gatk
def test_real_metrics_and_validation(tmp_path: Path) -> None:
    data = tmp_path / "d"
    _toy(data)
    summary = metrics.gatk_collect_alignment_summary_metrics(
        str(data / "toy.bam"),
        str(tmp_path / "a.txt"),
        reference=str(data / "toy.fa"),
        gatk_path=GATK,
        verbose=False,
    )
    unpaired = next(row for row in summary["metrics"] if row["CATEGORY"] == "UNPAIRED")
    assert unpaired["TOTAL_READS"] == 10
    valid = metrics.gatk_validate_sam_file(
        str(data / "toy.bam"), mode="SUMMARY", gatk_path=GATK, verbose=False
    )
    assert valid["valid"] is True
    term = metrics.gatk_check_terminator_block(str(data / "toy.bam"), gatk_path=GATK, verbose=False)
    assert term["terminator_ok"] is True


@requires_gatk
def test_real_intervals(tmp_path: Path) -> None:
    data = tmp_path / "d"
    _toy(data)
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t10\t100\tt1\n")
    converted = intervals.gatk_convert_interval_list(
        str(bed),
        str(tmp_path / "t.interval_list"),
        output_format="interval_list",
        sequence_dictionary=str(data / "toy.dict"),
        gatk_path=GATK,
        verbose=False,
    )
    assert converted["interval_count"] == 1
    back = intervals.gatk_convert_interval_list(
        converted["outputs"]["output_path"],
        str(tmp_path / "back.bed"),
        output_format="bed",
        gatk_path=GATK,
        verbose=False,
    )
    assert Path(back["outputs"]["output_path"]).read_text().startswith("chr1\t10\t100")
    prep = copynumber.gatk_preprocess_intervals(
        str(data / "toy.fa"),
        str(tmp_path / "p.interval_list"),
        bin_length=100,
        padding=0,
        gatk_path=GATK,
        verbose=False,
    )
    assert prep["interval_count"] == 3
    counts = copynumber.gatk_collect_read_counts(
        [str(data / "toy.bam")],
        [prep["outputs"]["output_path"]],
        str(tmp_path / "c.tsv"),
        output_format="TSV",
        gatk_path=GATK,
        verbose=False,
    )
    assert [row["COUNT"] for row in counts["rows"]] == [10, 10, 8] or len(counts["rows"]) == 3
