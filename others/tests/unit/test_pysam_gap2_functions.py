"""Tests for dnakit.bioinformatics.pysam_gap2_functions (pysam wrappers added for Fig. 2 coverage)."""

from __future__ import annotations

from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from dnakit.bioinformatics import pysam_gap2_functions as gap2
from dnakit.exceptions import ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)


def test_quality_string_round_trip() -> None:
    values = [0, 10, 20, 30, 40]
    text = gap2.array_to_qualitystring(values)
    assert text == "!+5?I"
    assert gap2.qualitystring_to_array(text) == values
    assert gap2.qualities_to_qualitystring(values) == text


def test_quality_offset_changes_encoding() -> None:
    assert gap2.array_to_qualitystring([0, 10], offset=64) == "@J"
    assert gap2.qualitystring_to_array("@J", offset=64) == [0, 10]


def test_quality_rejects_non_string() -> None:
    with pytest.raises(ConfigurationError):
        gap2.qualitystring_to_array(123)  # type: ignore[arg-type]


def test_bam_to_fastq_runs_on_shared_bam(tmp_path: Path) -> None:
    bam = DATA / "reads.bam"
    if not bam.exists():
        pytest.skip("shared test BAM not present")
    text = gap2.bam_to_fastq(["-n", str(bam)], str(tmp_path / "out.fq"))
    assert text.startswith("@") and text.count("\n") >= 4
    with pytest.raises(ConfigurationError):
        gap2.bam_to_fastq(["-n", str(bam)], str(tmp_path / "out.fq"))


def test_arguments_must_be_strings() -> None:
    with pytest.raises(ConfigurationError):
        gap2.bam_to_fastq("not-a-list", "/tmp/unused.fq")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        gap2.shuffle_bam(["ok", 3])  # type: ignore[list-item]
