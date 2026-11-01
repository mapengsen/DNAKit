"""Tests for deepTools BAM mapping statistics (batch-0 coverage-gap review)."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

pytest.importorskip("deeptools")
pytest.importorskip("pysam")

import pysam

from dnakit.epigenome import deeptools_batch0_review_functions as review
from dnakit.exceptions import ConfigurationError

SHARED_BAM = Path(
    "/zengxiangxiang/mps/DNAKit/paper/write_paper/figure/fig2/task_coverage/consistency/data/reads.bam"
)


def _write_bam(path: Path, *, index: bool) -> None:
    """Two references; four mapped reads (3 on chr1, 1 on chr2) and two unplaced unmapped reads."""

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [
        {"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}]}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for name, tid, pos in [("r1", 0, 10), ("r2", 0, 100), ("r3", 0, 300), ("r4", 1, 50)]:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "ACGTACGTAC"
            read.flag = 0
            read.reference_id = tid
            read.reference_start = pos
            read.mapping_quality = 60
            read.cigar = ((0, 10),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            out.write(read)
        for name in ("u1", "u2"):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "ACGTACGTAC"
            read.flag = 4
            read.reference_id = -1
            read.reference_start = -1
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            out.write(read)
    pysam.sort("-o", str(path), str(unsorted))
    unsorted.unlink()
    if index:
        pysam.index(str(path))


def test_public_functions_are_verbose_keyword_only() -> None:
    for name in review.__all__:
        function = getattr(review, name)
        verbose = inspect.signature(function).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True
        assert function.__doc__


def test_counts_from_an_indexed_bam(tmp_path: Path) -> None:
    bam = tmp_path / "small.bam"
    _write_bam(bam, index=True)
    result = review.bam_mapping_statistics(bam, verbose=False)
    assert result["status"] == "completed"
    assert result["operation"] == "bam_mapping_statistics"
    assert result["format"] == "BAM"
    assert result["mapped_reads"] == 4
    assert result["unmapped_reads"] == 2
    assert result["references"] == [
        {"name": "chr1", "mapped": 3, "unmapped": 0},
        {"name": "chr2", "mapped": 1, "unmapped": 0},
    ]
    assert result["provenance"]["package"].startswith("deeptools ")


def test_threads_are_checked_and_do_not_change_bam_counts(tmp_path: Path) -> None:
    bam = tmp_path / "threads.bam"
    _write_bam(bam, index=True)
    one = review.bam_mapping_statistics(bam, verbose=False)
    four = review.bam_mapping_statistics(bam, threads=4, verbose=False)
    assert one["mapped_reads"] == four["mapped_reads"]
    assert one["references"] == four["references"]
    with pytest.raises(ConfigurationError, match="threads"):
        review.bam_mapping_statistics(bam, threads=0, verbose=False)
    with pytest.raises(ConfigurationError, match="threads"):
        review.bam_mapping_statistics(bam, threads=True, verbose=False)


def test_unindexed_and_missing_files_are_rejected(tmp_path: Path) -> None:
    bam = tmp_path / "noindex.bam"
    _write_bam(bam, index=False)
    with pytest.raises(ConfigurationError, match="index"):
        review.bam_mapping_statistics(bam, verbose=False)
    with pytest.raises(ConfigurationError):
        review.bam_mapping_statistics(tmp_path / "absent.bam", verbose=False)


@pytest.mark.skipif(not SHARED_BAM.is_file(), reason="shared consistency data not present")
def test_shared_reads_bam_totals() -> None:
    result = review.bam_mapping_statistics(SHARED_BAM, verbose=False)
    assert result["mapped_reads"] == 5000
    assert result["unmapped_reads"] == 0
    assert sum(item["mapped"] for item in result["references"]) == 5000
    assert [item["name"] for item in result["references"]] == ["chr1", "chr2"]
