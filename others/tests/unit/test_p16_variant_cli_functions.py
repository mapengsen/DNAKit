"""mosdepth, Sniffles2, cuteSV and vcfanno facades (each part skips when its program is missing)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from _p16_data import genome, need, read_bam, sv_bam

import dnakit.variants.p16_variant_cli_functions as pv
from dnakit.exceptions import ConfigurationError


def test_mosdepth_depth_regions_and_quantized(tmp_path: Path) -> None:
    program = need("mosdepth")
    contigs = genome(tmp_path)
    bam = read_bam(tmp_path, contigs)
    depth = pv.compute_mosdepth_depth(bam, str(tmp_path / "d"), mosdepth_path=program, verbose=False)
    totals = {row["chrom"]: row for row in depth["summary"]}
    assert totals["total"]["length"] == 32000 and 5.0 < totals["chrA"]["mean"] < 6.5
    with gzip.open(tmp_path / "d.per-base.bed.gz", "rt") as handle:
        assert handle.readline().split("\t")[0] == "chrA"
    regions = pv.compute_mosdepth_regions(bam, str(tmp_path / "r"), by="1000", thresholds=[1, 10], mosdepth_path=program, verbose=False)
    assert len(regions["regions"]) == 32 and not (tmp_path / "r.per-base.bed.gz").exists()
    assert (tmp_path / "r.thresholds.bed.gz").exists() and (tmp_path / "r.mosdepth.region.dist.txt").exists()
    bed = tmp_path / "t.bed"
    bed.write_text("chrA\t0\t5000\nchrB\t100\t3000\n")
    by_bed = pv.compute_mosdepth_regions(bam, str(tmp_path / "b"), by=str(bed), use_median=True, mosdepth_path=program, verbose=False)
    assert [r["end"] for r in by_bed["regions"]] == [5000, 3000]
    quantized = pv.compute_mosdepth_quantized_coverage(bam, str(tmp_path / "q"), quantize="0:1:4:", mosdepth_path=program, verbose=False)
    assert quantized["segments"][0]["chrom"] == "chrA" and quantized["segments"][0]["label"]
    with pytest.raises(ConfigurationError):
        pv.compute_mosdepth_depth(bam, str(tmp_path / "d"), mosdepth_path=program, verbose=False)


def test_sniffles_cutesv_and_vcfanno(tmp_path: Path) -> None:
    sniffles, cutesv, vcfanno = need("sniffles"), need("cuteSV"), need("vcfanno")
    contigs = genome(tmp_path)
    bam = sv_bam(tmp_path, contigs)
    one = pv.call_structural_variants_sniffles2(bam, str(tmp_path / "s1.vcf"), snf_path=str(tmp_path / "s1.snf"), min_support=2,
                                                sniffles_path=sniffles, verbose=False)
    assert one["records"] >= 2 and (tmp_path / "s1.snf").exists()
    types = {line.split("SVTYPE=")[1].split(";")[0] for line in (tmp_path / "s1.vcf").read_text().splitlines() if not line.startswith("#") and "SVTYPE=" in line}
    assert {"DEL", "INS"} <= types
    mosaic = pv.call_structural_variants_sniffles2(bam, str(tmp_path / "s2.vcf"), mosaic=True, min_support=2, sniffles_path=sniffles, verbose=False)
    assert mosaic["mosaic"] is True
    two = pv.call_structural_variants_sniffles2(bam, str(tmp_path / "s3.vcf"), snf_path=str(tmp_path / "s3.snf"), min_support=2, sniffles_path=sniffles, verbose=False)
    joint = pv.merge_structural_variants_sniffles2_snf([str(tmp_path / "s1.snf"), str(tmp_path / "s3.snf")], str(tmp_path / "joint.vcf"),
                                                       min_support=1, sniffles_path=sniffles, verbose=False)
    assert joint["samples"] == 2 and two["records"] >= 2
    typed = pv.genotype_structural_variants_cutesv(bam, str(tmp_path / "genome.fa"), str(tmp_path / "cs.vcf"), str(tmp_path / "cs_work"),
                                                   platform="ont", min_support=2, cutesv_path=cutesv, verbose=False)
    body = [l for l in (tmp_path / "cs.vcf").read_text().splitlines() if not l.startswith("#")]
    assert typed["records"] == len(body) >= 2 and all(":" in l.split("\t")[8] and "GT" in l.split("\t")[8] for l in body)
    annotated = pv.annotate_vcf_from_bam_vcfanno(str(tmp_path / "s1.vcf"), bam, str(tmp_path / "anno.vcf"), vcfanno_path=vcfanno, verbose=False)
    info = [l for l in (tmp_path / "anno.vcf").read_text().splitlines() if not l.startswith("#")][0].split("\t")[7]
    assert "sample_depth=" in info and "sample_mapq=" in info and annotated["records"] >= 2


def test_sniffles_multi_sample_needs_two_files(tmp_path: Path) -> None:
    one = tmp_path / "a.snf"
    one.write_text("x")
    with pytest.raises(ConfigurationError):
        pv.merge_structural_variants_sniffles2_snf([str(one)], str(tmp_path / "o.vcf"), verbose=False)
    with pytest.raises(ConfigurationError, match="operations"):
        pv.annotate_vcf_from_bam_vcfanno(str(one), str(one), str(tmp_path / "x.vcf"), operations=["gc"], verbose=False)
