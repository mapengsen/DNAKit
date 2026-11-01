"""Tests for dnakit.alignment.samtools_p06_cli_functions (samtools fqidx, coverage -m and checksum facades)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from _p06_htslib_data import DATA, requires_data, requires_samtools

from dnakit.alignment import samtools_p06_cli_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pytestmark = [requires_samtools, requires_data]


def _samtools(*arguments: str) -> str:
    return subprocess.run(["samtools", *arguments], check=True, capture_output=True, text=True).stdout


def _fastq(tmp_path: Path) -> Path:
    target = tmp_path / "reads.fastq"
    shutil.copyfile(DATA / "reads.fastq", target)
    return target


def _records(path: Path) -> dict[str, tuple[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return {lines[i][1:]: (lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)}


def test_fastq_region_matches_the_source_record(tmp_path: Path) -> None:
    source = _fastq(tmp_path)
    full = _records(DATA / "reads.fastq")["read00000"]
    out = tmp_path / "region.fq"
    result = m.samtools_extract_fastq_regions(
        source, out, regions=["read00000:5-24", "read00001"], verbose=False
    )
    assert result["status"] == "completed" and result["backend"] == "samtools"
    got = _records(out)
    assert got["read00000:5-24"] == (full[0][4:24], full[1][4:24])
    assert got["read00001"] == _records(DATA / "reads.fastq")["read00001"]
    assert result["records"] == [
        {"name": "read00000:5-24", "length": 20},
        {"name": "read00001", "length": len(got["read00001"][0])},
    ]
    assert (tmp_path / "reads.fastq.fai").is_file()


def test_fastq_region_equals_fqidx_and_faidx_dash_f(tmp_path: Path) -> None:
    source = _fastq(tmp_path)
    out = tmp_path / "a.fq"
    m.samtools_extract_fastq_regions(source, out, regions=["read00002:10-30"], verbose=False)
    assert out.read_text(encoding="utf-8") == _samtools("fqidx", str(source), "read00002:10-30")
    assert out.read_text(encoding="utf-8") == _samtools("faidx", "-f", str(source), "read00002:10-30")


def test_fastq_reverse_complement_and_strand_mark(tmp_path: Path) -> None:
    source = _fastq(tmp_path)
    forward = tmp_path / "f.fq"
    reverse = tmp_path / "r.fq"
    m.samtools_extract_fastq_regions(source, forward, regions=["read00003:1-12"], verbose=False)
    m.samtools_extract_fastq_regions(
        source, reverse, regions=["read00003:1-12"], reverse_complement=True, mark_strand="sign", verbose=False
    )
    seq, qual = next(iter(_records(forward).values()))
    name, (rseq, rqual) = next(iter(_records(reverse).items()))
    assert name == "read00003:1-12(-)"
    assert rseq == seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    assert rqual == qual[::-1]


def test_fastq_region_file_index_path_and_missing_region(tmp_path: Path) -> None:
    source = _fastq(tmp_path)
    regions = tmp_path / "regions.txt"
    regions.write_text("read00000:1-5\nread00004:2-6\n", encoding="utf-8")
    out = tmp_path / "from_file.fq"
    index = tmp_path / "custom.fai"
    result = m.samtools_extract_fastq_regions(
        source, out, region_file=regions, index_path=index, verbose=False
    )
    assert [row["name"] for row in result["records"]] == ["read00000:1-5", "read00004:2-6"]
    assert index.is_file() and not (tmp_path / "reads.fastq.fai").exists()
    with pytest.raises(Exception, match="exited with code"):
        m.samtools_extract_fastq_regions(
            source, tmp_path / "bad.fq", regions=["nosuchread:1-5"], index_path=index, verbose=False
        )
    assert not (tmp_path / "bad.fq").exists()
    skipped = m.samtools_extract_fastq_regions(
        source, tmp_path / "skip.fq", regions=["nosuchread:1-5", "read00000:1-5"],
        continue_on_missing=True, index_path=index, verbose=False,
    )
    assert skipped["records"][-1]["name"] == "read00000:1-5"


def test_fastq_validation_and_missing_program(tmp_path: Path) -> None:
    source = _fastq(tmp_path)
    out = tmp_path / "x.fq"
    with pytest.raises(ConfigurationError):
        m.samtools_extract_fastq_regions(source, out, verbose=False)
    with pytest.raises(ConfigurationError):
        m.samtools_extract_fastq_regions(source, out, regions=["-bad"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.samtools_extract_fastq_regions(source, out, regions=["a"], mark_strand="weird", verbose=False)
    out.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.samtools_extract_fastq_regions(source, out, regions=["read00000"], verbose=False)
    with pytest.raises(BackendUnavailableError):
        m.samtools_extract_fastq_regions(
            source, tmp_path / "y.fq", regions=["read00000"], samtools_path=str(tmp_path / "nope"), verbose=False
        )


def test_coverage_histogram_statistics_match_the_table(tmp_path: Path) -> None:
    bam = DATA / "reads.bam"
    out = tmp_path / "hist.txt"
    result = m.samtools_coverage_histogram([bam], out, ascii_only=True, bins=30, verbose=False)
    table = [line.split("\t") for line in _samtools("coverage", str(bam)).splitlines()[1:]]
    assert [c["name"] for c in result["contigs"]] == [row[0] for row in table]
    for contig, row in zip(result["contigs"], table, strict=True):
        assert contig["statistics"]["Number of reads"] == int(row[3])
        assert contig["statistics"]["Percent covered"] == pytest.approx(float(row[5]), abs=0.01)
        assert contig["statistics"]["Mean coverage"] == pytest.approx(float(row[6]), abs=0.01)
    text = out.read_text(encoding="utf-8")
    assert text == _samtools("coverage", "-m", "-A", "-w", "30", str(bam))
    assert all(ord(ch) < 128 for ch in text)


def test_coverage_histogram_region_and_filters(tmp_path: Path) -> None:
    bam = DATA / "reads.bam"
    result = m.samtools_coverage_histogram(
        [bam], tmp_path / "region.txt", region="chr2:1-20000", min_mapping_quality=10, verbose=False
    )
    assert [c["name"] for c in result["contigs"]] == ["chr2"]
    direct = _samtools("coverage", "-m", "-q", "10", "-r", "chr2:1-20000", "-w", "40", str(bam))
    assert (tmp_path / "region.txt").read_text(encoding="utf-8") == direct
    with pytest.raises(ConfigurationError):
        m.samtools_coverage_histogram([bam], tmp_path / "bad.txt", required_flags="not a flag", verbose=False)
    with pytest.raises(ConfigurationError):
        m.samtools_coverage_histogram([], tmp_path / "bad2.txt", verbose=False)


def test_checksum_is_order_agnostic_and_merges(tmp_path: Path) -> None:
    bam = DATA / "reads.bam"
    shuffled = tmp_path / "byname.bam"
    subprocess.run(["samtools", "sort", "-n", "-o", str(shuffled), str(bam)], check=True, capture_output=True)
    first = m.samtools_checksum([bam], tmp_path / "a.chk", verbose=False)
    second = m.samtools_checksum([shuffled], tmp_path / "b.chk", verbose=False)
    assert first["groups"] == second["groups"]
    assert first["groups"][0]["count"] == 5000 and first["groups"][0]["Group"] == "all"
    assert (tmp_path / "a.chk").read_text(encoding="utf-8").splitlines()[3:] == _samtools(
        "checksum", str(bam)
    ).splitlines()[3:]
    ordered = m.samtools_checksum([bam], tmp_path / "o.chk", in_order=True, check_position=True, verbose=False)
    assert ordered["groups"][0]["combined"] != first["groups"][0]["combined"]
    merged = m.samtools_checksum(
        [tmp_path / "a.chk", tmp_path / "b.chk"], tmp_path / "m.chk", merge=True, verbose=False
    )
    assert merged["merged"] is True and merged["groups"][0]["count"] == 10000


def test_checksum_tabs_and_limits(tmp_path: Path) -> None:
    bam = DATA / "reads.bam"
    tabbed = m.samtools_checksum([bam], tmp_path / "t.chk", tab_separated=True, max_records=100, verbose=False)
    assert tabbed["groups"][0]["count"] == 100
    assert "\t" in (tmp_path / "t.chk").read_text(encoding="utf-8").splitlines()[-1]
    with pytest.raises(ConfigurationError):
        m.samtools_checksum([bam], tmp_path / "x.chk", exclude_flags=";rm", verbose=False)
    with pytest.raises(ConfigurationError):
        m.samtools_checksum([], tmp_path / "y.chk", verbose=False)
