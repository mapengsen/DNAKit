"""Tests for dnakit.annotation.augustus_p04_hints_functions (AUGUSTUS hints and alignment converters)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.annotation import augustus_p04_hints_functions as h
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
GENOME = DATA / "genome.fa"
BAM = DATA / "reads.bam"
requires_augustus = pytest.mark.skipif(
    shutil.which("bam2hints") is None or shutil.which("blat2hints.pl") is None,
    reason="AUGUSTUS (conda install -c bioconda augustus) is not on PATH",
)
PSL = "100\t0\t0\t0\t0\t0\t1\t300\t+\tq1\t100\t0\t100\tchr1\t120000\t1000\t1400\t2\t50,50,\t0,50,\t1000,1350,\n"
MAF = (
    "##maf version=1\n\na score=100\ns hg.chr1 100 20 + 1000 ACGT-ACGTACGTACGTACGT\n"
    "s mm.chr2 200 20 + 1000 ACGTTACGTACGTACGTACGT\n\n"
)


def test_missing_program_and_validation(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        h.augustus_bam_to_hints(str(BAM), str(tmp_path / "o"), augustus_bin_dir=str(tmp_path), verbose=False)
    with pytest.raises(ConfigurationError):
        h.augustus_polya_to_hints(str(BAM), str(tmp_path / "o"), input_format="x", verbose=False)
    with pytest.raises(ConfigurationError):
        h.augustus_filter_splice_hints(str(GENOME), str(BAM), str(tmp_path / "o"), splice_pattern="GT", verbose=False)
    with pytest.raises(ConfigurationError):
        h.augustus_bam_to_wiggle(str(BAM), str(tmp_path / "o"), method="bamToWig.py", verbose=False)
    with pytest.raises(ConfigurationError):
        h.augustus_filter_maf(str(BAM), str(tmp_path / "o"), interval="chr1", verbose=False)


@requires_augustus
def test_hints_from_psl_and_bam(tmp_path: Path) -> None:
    psl = tmp_path / "a.psl"
    psl.write_text(PSL)
    result = h.augustus_psl_to_hints(str(psl), str(tmp_path / "psl.hints"), verbose=False)
    lines = Path(result["outputs"]["output"]).read_text().splitlines()
    assert any("\tintron\t1051\t1350\t" in line for line in lines)
    assert result["version"]
    bam = h.augustus_bam_to_hints(str(BAM), str(tmp_path / "bam.hints"), verbose=False)
    assert bam["status"] == "completed" and bam["lines"] == 0  # unspliced reads give no intron hints
    with pytest.raises(ConfigurationError):
        h.augustus_psl_to_hints(str(psl), str(tmp_path / "psl.hints"), verbose=False)


@requires_augustus
def test_coverage_converters(tmp_path: Path) -> None:
    wig = h.augustus_bam_to_wiggle(str(BAM), str(tmp_path / "cov.wig"), verbose=False)
    assert "variableStep chrom=chr1" in Path(wig["outputs"]["output"]).read_text()
    hints = h.augustus_wiggle_to_hints(wig["outputs"]["output"], str(tmp_path / "wig.hints"), priority=2, verbose=False)
    assert hints["lines"] > 0 and "pri=2" in hints["preview"][0]
    chosen = h.augustus_choose_wiggle_intervals(wig["outputs"]["output"], str(tmp_path / "choose.bed"), verbose=False)
    assert chosen["lines"] >= 1
    bedgraph = tmp_path / "c.bedgraph"
    bedgraph.write_text("chr1\t0\t100\t3\nchr1\t100\t200\t5\n")
    converted = h.augustus_bedgraph_to_wiggle(str(bedgraph), str(tmp_path / "c.wig"), verbose=False)
    assert converted["lines"] == 2
    psl = tmp_path / "a.psl"
    psl.write_text(PSL)
    aln = h.augustus_alignments_to_wiggle(str(psl), str(tmp_path / "aln.wig"), verbose=False)
    assert aln["lines"] > 100


@requires_augustus
def test_filters_and_maf(tmp_path: Path) -> None:
    psl = tmp_path / "a.psl"
    psl.write_text(PSL)
    kept = h.augustus_filter_psl(str(psl), str(tmp_path / "f.psl"), verbose=False)
    assert kept["lines"] == 1
    dropped = h.augustus_filter_psl(str(psl), str(tmp_path / "f2.psl"), min_identity=100, min_coverage=100, no_introns=True, verbose=False)
    assert dropped["lines"] == 0
    filtered = h.augustus_filter_bam(str(BAM), str(tmp_path / "f.bam"), verbose=False)
    assert filtered["bytes"] > 0
    maf = tmp_path / "t.maf"
    maf.write_text(MAF)
    inside = h.augustus_filter_maf(str(maf), str(tmp_path / "in.maf"), interval="hg.chr1:50-300", verbose=False)
    assert inside["lines"] == 3
    outside = h.augustus_filter_maf(str(maf), str(tmp_path / "out.maf"), interval="hg.chr1:500-600", verbose=False)
    assert outside["lines"] == 0
    weeded = h.augustus_weed_maf(str(maf), str(tmp_path / "weed.maf"), verbose=False)
    assert weeded["status"] == "completed"
    cons = h.augustus_maf_to_conservation_wiggle(str(maf), str(tmp_path / "cons.wig"), verbose=False)
    assert cons["preview"][0].startswith("fixedStep chrom=chr1")


@requires_augustus
def test_hint_tools(tmp_path: Path) -> None:
    hints = tmp_path / "h.gff"
    hints.write_text("chr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=a;pri=4;src=E\nchr1\tb2h\tintron\t1051\t1350\t0\t+\t.\tgrp=b;pri=4;src=E\n")
    joined = h.augustus_join_multiplicity_hints(str(hints), str(tmp_path / "j.gff"), verbose=False)
    assert joined["lines"] == 1 and "mult=2" in joined["preview"][0]
    redundant = h.augustus_remove_redundant_hints(str(hints), str(hints), str(tmp_path / "r.gff"), verbose=False)
    assert redundant["status"] == "completed"
    filtered = h.augustus_filter_splice_hints(str(GENOME), str(hints), str(tmp_path / "s.gff"), verbose=False)
    assert filtered["status"] == "completed"
