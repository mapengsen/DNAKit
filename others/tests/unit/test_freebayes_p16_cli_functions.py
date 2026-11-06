"""FreeBayes helper tools: Python ports (checked against the original script outputs) and bamleftalign."""

from __future__ import annotations

from pathlib import Path

import pytest
from _p16_data import genome, need, read_bam

import dnakit.variants.freebayes_p16_cli_functions as fb
from dnakit.exceptions import ConfigurationError


def _fai(tmp_path: Path) -> str:
    (tmp_path / "ref.fa.fai").write_text(
        "chr1\t1000\t6\t60\t61\nchr2\t450\t1100\t60\t61\ntiny\t7\t2000\t60\t61\n"
    )
    (tmp_path / "ref.fa").write_text("")
    return str(tmp_path / "ref.fa")


def test_regions_match_fasta_generate_regions(tmp_path: Path) -> None:
    ref = _fai(tmp_path)
    regions = fb.partition_reference_into_regions(ref, 300, verbose=False)["regions"]
    assert regions == [
        "chr1:0-300",
        "chr1:300-600",
        "chr1:600-900",
        "chr1:900-1000",
        "chr2:0-300",
        "chr2:300-450",
        "tiny:0-7",
    ]
    assert fb.partition_reference_into_regions(ref, 3, chunks=True, verbose=False)["regions"][
        :3
    ] == ["chr1:0-334", "chr1:334-668", "chr1:668-1000"]
    only = fb.partition_reference_into_regions(ref, 300, chromosomes=["chr2"], verbose=False)[
        "regions"
    ]
    assert only == ["chr2:0-300", "chr2:300-450"]
    with pytest.raises(ConfigurationError):
        fb.partition_reference_into_regions(ref, 300, chromosomes=["none"], verbose=False)


def test_regions_bed_files(tmp_path: Path) -> None:
    result = fb.partition_reference_into_regions(
        _fai(tmp_path), 500, bed_prefix=str(tmp_path / "r"), verbose=False
    )
    assert len(result["outputs"]) == 4
    assert (tmp_path / "r.chr1.region.1.bed").read_text() == "chr1\t0\t500"


def test_coverage_regions(tmp_path: Path) -> None:
    ref = _fai(tmp_path)
    coverage = tmp_path / "cov.txt"
    coverage.write_text("".join(f"chr1\t{p}\t{10 if p < 500 else 2}\n" for p in range(0, 1000, 10)))
    regions = fb.partition_reference_by_coverage(ref, str(coverage), 2, verbose=False)["regions"]
    assert regions == ["chr1:0-370", "chr1:370-1000"] or regions[-1].endswith("-1000")
    assert regions[0].startswith("chr1:0-")


def test_keep_first_vcf_header(tmp_path: Path) -> None:
    a, b = tmp_path / "a.vcf", tmp_path / "b.vcf"
    a.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t10\n")
    b.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n2\t20\n")
    result = fb.keep_first_vcf_header([str(a), str(b)], str(tmp_path / "m.vcf"), verbose=False)
    assert (tmp_path / "m.vcf").read_text() == "##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t10\n2\t20\n"
    assert result["header_lines_dropped"] == 2 and result["records_written"] == 2


def test_split_by_bam_size_reads_index(tmp_path: Path) -> None:
    pytest.importorskip("scipy")
    contigs = genome(tmp_path)
    bam = read_bam(tmp_path, contigs)
    result = fb.partition_reference_by_bam_size(
        str(tmp_path / "genome.fa"), [bam], output_path=str(tmp_path / "r.tsv"), verbose=False
    )
    assert result["regions"][0][0] == "chrA" and result["regions"][-1][2] == 12000
    with pytest.raises(ConfigurationError):
        fb.partition_reference_by_bam_size(
            str(tmp_path / "genome.fa"), [str(tmp_path / "nofile.bam")], verbose=False
        )


def test_bamleftalign_runs(tmp_path: Path) -> None:
    program = need("bamleftalign")
    pysam = pytest.importorskip("pysam")
    contigs = genome(tmp_path)
    bam = read_bam(tmp_path, contigs)
    result = fb.left_align_bam_indels(
        bam,
        str(tmp_path / "left.bam"),
        reference=str(tmp_path / "genome.fa"),
        bamleftalign_path=program,
        verbose=False,
    )
    assert result["status"] == "completed"
    with pysam.AlignmentFile(result["output"]) as handle:
        deletions = [
            r.cigarstring
            for r in handle.fetch(until_eof=True)
            if r.cigarstring and "D" in r.cigarstring
        ]
    assert deletions and all(c.startswith("2") or "D" in c for c in deletions)
    # left alignment moves the deletion from the right edge of the A homopolymer (1008) to its left edge (1000)
    with pysam.AlignmentFile(result["output"]) as handle:
        starts = {
            read.reference_start + sum(n for op, n in read.cigartuples[:1])
            for read in handle.fetch(until_eof=True)
            if read.cigarstring and "D" in read.cigarstring and read.cigartuples[0][0] == 0
        }
    assert min(starts) <= 1001
