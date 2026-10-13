"""Unit tests for the G3 annotation wrappers.

Covers bioframe, pybedtools, gffutils, gtfparse, pyBigWig, and pyranges1. Each wrapped function is
checked against a hand calculation or a direct upstream call. Missing optional packages skip their
tests through pytest.importorskip. Fixtures are written inside tmp_path.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError

BED_A = "chr1\t10\t20\ta\t5\t+\nchr1\t40\t50\tc\t1\t-\nchr2\t5\t9\td\t2\t+\n"
BED_A_FULL = (
    "chr1\t10\t20\ta\t5\t+\nchr1\t15\t30\tb\t7\t+\nchr1\t40\t50\tc\t1\t-\nchr2\t5\t9\td\t2\t+\n"
)
BED_B = "chr1\t15\t18\tx\nchr1\t45\t60\ty\nchr2\t100\t120\tz\n"
GFF3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=GENE1\n"
    "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1;Name=GENE1-201\n"
    "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=tx1\n"
    "chr1\tsrc\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=tx1\n"
    "chr2\tsrc\tgene\t50\t90\t.\t-\t.\tID=gene2\n"
)
GTF = (
    'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "g1"; gene_name "GENE1";\n'
    'chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number 1;\n'
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_overlap_bioframe_intervals(tmp_path):
    pytest.importorskip("bioframe")
    import bioframe

    from dnakit.annotation.bioframe_functions import overlap_bioframe_intervals

    a = _write(tmp_path / "a.bed", BED_A_FULL)
    b = _write(tmp_path / "b.bed", BED_B)
    result = overlap_bioframe_intervals(a, b)
    spans = sorted((r["chrom"], r["start"], r["end"]) for r in result["records"])
    assert spans == [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 40, 50)]
    direct = bioframe.overlap(
        bioframe.read_table(str(a), schema="bed6"),
        bioframe.read_table(str(b), schema="bed4"),
        how="inner",
    )
    assert result["record_count"] == len(direct) == 3
    with pytest.raises(ConfigurationError):
        overlap_bioframe_intervals(a, b, how="cross")


def test_subtract_bioframe_intervals(tmp_path):
    pytest.importorskip("bioframe")
    from dnakit.annotation.bioframe_functions import subtract_bioframe_intervals

    a = _write(tmp_path / "a.bed", BED_A_FULL)
    b = _write(tmp_path / "b.bed", BED_B)
    result = subtract_bioframe_intervals(a, b)
    spans = sorted((r["chrom"], r["start"], r["end"]) for r in result["intervals"])
    # chr1 10-20 minus 15-18 -> 10-15, 18-20; chr1 15-30 minus 15-18 -> 18-30;
    # chr1 40-50 minus 45-60 -> 40-45; chr2 5-9 is untouched.
    assert spans == [
        ("chr1", 10, 15),
        ("chr1", 18, 20),
        ("chr1", 18, 30),
        ("chr1", 40, 45),
        ("chr2", 5, 9),
    ]
    assert result["interval_count"] == 5


def test_expand_bioframe_intervals(tmp_path):
    pytest.importorskip("bioframe")
    from dnakit.annotation.bioframe_functions import expand_bioframe_intervals

    a = _write(tmp_path / "a.bed", BED_A)
    result = expand_bioframe_intervals(a, pad=5)
    spans = sorted((r["chrom"], r["start"], r["end"]) for r in result["intervals"])
    assert spans == [("chr1", 5, 25), ("chr1", 35, 55), ("chr2", 0, 14)]
    assert result["pad"] == 5


def test_cluster_bioframe_intervals(tmp_path):
    pytest.importorskip("bioframe")
    from dnakit.annotation.bioframe_functions import cluster_bioframe_intervals

    a = _write(tmp_path / "a.bed", BED_A_FULL)
    result = cluster_bioframe_intervals(a)
    by_span = {(r["chrom"], r["start"], r["end"]): r["cluster"] for r in result["records"]}
    # chr1 10-20 and chr1 15-30 overlap, so they share a cluster; 40-50 and chr2 are separate.
    assert by_span[("chr1", 10, 20)] == by_span[("chr1", 15, 30)]
    assert len({by_span[("chr1", 40, 50)], by_span[("chr2", 5, 9)], by_span[("chr1", 10, 20)]}) == 3


def test_binnify_bioframe_chromosomes():
    pytest.importorskip("bioframe")
    from dnakit.annotation.bioframe_functions import binnify_bioframe_chromosomes

    result = binnify_bioframe_chromosomes({"chr1": 25, "chr2": 10}, binsize=10)
    spans = [(r["chrom"], r["start"], r["end"]) for r in result["bins"]]
    assert spans == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 25), ("chr2", 0, 10)]
    assert result["bin_count"] == 4


def test_parse_bioframe_region():
    pytest.importorskip("bioframe")
    from dnakit.annotation.bioframe_functions import parse_bioframe_region

    result = parse_bioframe_region("chr1:1,101-2,000")
    assert (result["chrom"], result["start"], result["end"]) == ("chr1", 1101, 2000)
    with pytest.raises(ConfigurationError):
        parse_bioframe_region("not a region")


def test_five_prime_pybedtools_coordinates(tmp_path):
    pytest.importorskip("pybedtools")
    from dnakit.annotation.pybedtools_functions import five_prime_pybedtools_coordinates

    a = _write(tmp_path / "a.bed", BED_A)
    result = five_prime_pybedtools_coordinates(a, upstream=5, downstream=5)
    spans = [(r["start"], r["end"]) for r in result["records"]]
    # + strand 10-20 -> 5' end 10 -> 5-15; - strand 40-50 -> 5' end 50 -> 45-55; 5-9 -> 0-10.
    assert sorted(spans) == [(0, 10), (5, 15), (45, 55)]


def test_midpoint_pybedtools_coordinates(tmp_path):
    pytest.importorskip("pybedtools")
    from dnakit.annotation.pybedtools_functions import midpoint_pybedtools_coordinates

    a = _write(tmp_path / "a.bed", BED_A)
    result = midpoint_pybedtools_coordinates(a)
    assert sorted((r["start"], r["end"]) for r in result["records"]) == [(7, 8), (15, 16), (45, 46)]


def test_tag_pybedtools_bedpe(tmp_path):
    pytest.importorskip("pybedtools")
    from dnakit.annotation.pybedtools_functions import tag_pybedtools_bedpe

    pairs = _write(tmp_path / "pairs.bedpe", "chr1\t10\t20\tchr1\t100\t110\tp1\t0\t+\t-\n")
    queries = _write(tmp_path / "q.bed", "chr1\t12\t18\tq1\nchr1\t101\t105\tq2\n")
    result = tag_pybedtools_bedpe(pairs, queries)
    # q1 covers 12-18 inside 10-20 (6 bp); q2 covers 101-105 inside 100-110 (4 bp).
    assert result["pairs"] == [
        {
            "name": "p1",
            "left": [{"query": "q1", "overlap_bp": 6}],
            "right": [{"query": "q2", "overlap_bp": 4}],
        }
    ]


def test_write_pybedtools_chromsizes(tmp_path):
    pytest.importorskip("pybedtools")
    from dnakit.annotation.pybedtools_functions import write_pybedtools_chromsizes

    target = tmp_path / "sizes.txt"
    result = write_pybedtools_chromsizes({"chr2": 50, "chr1": 100}, target)
    assert target.read_text(encoding="utf-8").split() == ["chr1", "100", "chr2", "50"]
    assert result["chromosome_count"] == 2
    with pytest.raises(ConfigurationError):
        write_pybedtools_chromsizes({"chr1": 100}, target)


def _gffutils_db(tmp_path: Path) -> Path:
    gffutils = pytest.importorskip("gffutils")
    gff = _write(tmp_path / "genes.gff3", GFF3)
    db = tmp_path / "genes.db"
    gffutils.create_db(str(gff), str(db), force=True, verbose=False, merge_strategy="create_unique")
    return db


def test_to_bed12_gffutils_transcript(tmp_path):
    db = _gffutils_db(tmp_path)
    from dnakit.annotation.gffutils_functions import to_bed12_gffutils_transcript

    fields = to_bed12_gffutils_transcript(db, "tx1")["bed12"].split("\t")
    # exon1 1-based 100-200 -> BED 99-200 (size 101); exon2 300-500 -> BED 299-500 (size 201).
    assert fields[1:3] == ["99", "500"]
    assert fields[9:] == ["2", "101,201", "0,200"]
    with pytest.raises(ConfigurationError):
        to_bed12_gffutils_transcript(db, "missing")


def test_tss_gffutils_intervals(tmp_path):
    pytest.importorskip("pybedtools")
    db = _gffutils_db(tmp_path)
    from dnakit.annotation.gffutils_functions import tss_gffutils_intervals

    records = tss_gffutils_intervals(db)["records"]
    # The only transcript tx1 starts at 1-based 100, which is BED start 99.
    assert records == [{"chrom": "chr1", "start": 99, "end": 100, "name": "tx1", "strand": "+"}]


def test_parse_gffutils_feature_line():
    pytest.importorskip("gffutils")
    from dnakit.annotation.gffutils_functions import parse_gffutils_feature_line

    line = "chr1\tsrc\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=tx1"
    result = parse_gffutils_feature_line(line)
    assert (result["featuretype"], result["start"], result["end"]) == ("exon", 300, 500)
    assert result["attributes"] == {"ID": ["exon2"], "Parent": ["tx1"]}
    with pytest.raises(ConfigurationError):
        parse_gffutils_feature_line("chr1 only three fields")


def test_inspect_gffutils_source(tmp_path):
    pytest.importorskip("gffutils")
    from dnakit.annotation.gffutils_functions import inspect_gffutils_source

    gff = _write(tmp_path / "genes.gff3", GFF3)
    result = inspect_gffutils_source(gff)
    assert result["featuretype"] == {"exon": 2, "gene": 2, "mRNA": 1}
    assert result["feature_count"] == 5
    assert result["chrom"] == {"chr1": 4, "chr2": 1}


def test_parse_gtfparse_fixed_columns(tmp_path):
    gtfparse = pytest.importorskip("gtfparse")
    from dnakit.annotation.gtfparse_functions import parse_gtfparse_fixed_columns

    gtf = _write(tmp_path / "genes.gtf", GTF)
    result = parse_gtfparse_fixed_columns(gtf)
    assert result["record_count"] == len(gtfparse.parse_gtf(str(gtf), split_attributes=False)) == 3
    exon = result["records"][2]
    assert (exon["feature"], exon["start"], exon["end"]) == ("exon", 100, 200)
    assert "exon_number 1" in exon["attribute"]


def test_expand_gtfparse_attributes(tmp_path):
    pytest.importorskip("gtfparse")
    from dnakit.annotation.gtfparse_functions import expand_gtfparse_attributes

    gtf = _write(tmp_path / "genes.gtf", GTF)
    result = expand_gtfparse_attributes(gtf)
    assert {"gene_id", "transcript_id", "exon_number"} <= set(result["columns"])
    assert str(result["records"][2]["exon_number"]) == "1"


def test_write_gtfparse_table(tmp_path):
    pytest.importorskip("gtfparse")
    from dnakit.annotation.gtfparse_functions import write_gtfparse_table

    gtf = _write(tmp_path / "genes.gtf", GTF)
    target = tmp_path / "round.gtf"
    result = write_gtfparse_table(gtf, target)
    lines = target.read_text(encoding="utf-8").splitlines()
    assert result["line_count"] == len(lines) == 3
    assert lines[2].split("\t")[3:5] == ["100", "200"]


def _bigwig_fixture(tmp_path: Path) -> tuple[Path, Path]:
    pytest.importorskip("pyBigWig")
    bedgraph = _write(
        tmp_path / "cov.bedgraph", "chr1\t0\t10\t1.5\nchr1\t10\t20\t3\nchr2\t5\t9\t2\n"
    )
    return bedgraph, tmp_path / "cov.bw"


def test_read_bigwig_header(tmp_path):
    bedgraph, target = _bigwig_fixture(tmp_path)
    from dnakit.annotation.pyBigWig_functions import read_bigwig_header, write_bigwig_intervals

    write_bigwig_intervals(bedgraph, target, chromsizes={"chr1": 100, "chr2": 50})
    result = read_bigwig_header(target)
    assert result["chromosomes"] == [
        {"chrom": "chr1", "length": 100},
        {"chrom": "chr2", "length": 50},
    ]
    # Covered bases: 10 + 10 + 4.
    assert result["header"]["nBasesCovered"] == 24


def test_write_bigwig_intervals(tmp_path):
    bedgraph, target = _bigwig_fixture(tmp_path)
    from dnakit.annotation.pyBigWig_functions import write_bigwig_intervals

    result = write_bigwig_intervals(bedgraph, target, chromsizes={"chr1": 100, "chr2": 50})
    assert result["interval_count"] == 3
    with pytest.raises(ConfigurationError):
        write_bigwig_intervals(bedgraph, target, chromsizes={"chr1": 100, "chr2": 50})


def test_summarize_bigwig_bins(tmp_path):
    bedgraph, target = _bigwig_fixture(tmp_path)
    from dnakit.annotation.pyBigWig_functions import summarize_bigwig_bins, write_bigwig_intervals

    write_bigwig_intervals(bedgraph, target, chromsizes={"chr1": 100, "chr2": 50})
    result = summarize_bigwig_bins(target, chrom="chr1", start=0, end=20, bins=2)
    # [0, 10) has value 1.5 and [10, 20) has value 3.
    assert result["values"] == [1.5, 3.0]


def test_tile_pyranges1_genome():
    pytest.importorskip("pyranges1")
    from dnakit.annotation.pyranges1_functions import tile_pyranges1_genome

    result = tile_pyranges1_genome({"chr1": 25, "chr2": 10}, tile_size=10)
    spans = [(r["chrom"], r["start"], r["end"]) for r in result["tiles"]]
    assert spans == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 25), ("chr2", 0, 10)]


def test_count_pyranges1_overlaps(tmp_path):
    pytest.importorskip("pyranges1")
    from dnakit.annotation.pyranges1_functions import count_pyranges1_overlaps

    a = _write(tmp_path / "a.bed", BED_A)
    b = _write(tmp_path / "b.bed", BED_B)
    result = count_pyranges1_overlaps(a, b)
    counts = [(r["chrom"], r["start"], r["overlaps"]) for r in result["records"]]
    # chr1 10-20 meets chr1 15-18 and chr1 40-50 meets chr1 45-60; chr2 5-9 has no partner.
    assert counts == [("chr1", 10, 1), ("chr1", 40, 1), ("chr2", 5, 0)]


def test_parse_pyranges1_string():
    pytest.importorskip("pyranges1")
    from dnakit.annotation.pyranges1_functions import parse_pyranges1_string

    result = parse_pyranges1_string("Chromosome Start End Strand\nchr1 10 20 +\n")
    assert result["records"] == [{"chrom": "chr1", "start": 10, "end": 20, "strand": "+"}]
    with pytest.raises(ConfigurationError):
        parse_pyranges1_string("not a table")


def test_jaccard_pyranges1_intervals(tmp_path):
    pyranges1 = pytest.importorskip("pyranges1")
    from dnakit.annotation.pyranges1_functions import jaccard_pyranges1_intervals

    a = _write(tmp_path / "a.bed", BED_A)
    b = _write(tmp_path / "b.bed", BED_B)
    result = jaccard_pyranges1_intervals(a, b)
    # Shared bases: 15-18 (3 bp) and 45-50 (5 bp) = 8 bp; union: 30 bp on chr1 + 24 bp on chr2 = 54.
    assert result["jaccard"] == pytest.approx(8 / 54)
    direct = pyranges1.stats.jaccard(pyranges1.read_bed(str(a)), pyranges1.read_bed(str(b)))
    assert result["jaccard"] == pytest.approx(direct)
