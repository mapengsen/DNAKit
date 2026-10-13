"""Contract tests for optional interval, annotation, and region-set backends."""

from __future__ import annotations

import random
import shutil
import sys
from collections import Counter
from pathlib import Path

import pytest

from dnakit.annotation import backends
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

A_BED = "chr1\t10\t20\ta\t5\t+\nchr1\t15\t30\tb\t7\t+\nchr1\t40\t50\tc\t1\t-\nchr2\t5\t9\td\t2\t+\n"
B_BED = "chr1\t12\t18\tx\nchr1\t45\t60\ty\nchr2\t100\t120\tz\n"
TOUCH_BED = "chr1\t0\t10\nchr1\t10\t20\nchr1\t30\t40\n"
SEQ_BED = "chr1\t10\t20\tp1\t0\t+\nchr1\t0\t5\tp2\t0\t-\nchr2\t5\t9\tp3\t0\t+\n"
FASTA = ">chr1\nACGTACGTAC\nGGGGTTTTCC\n>chr2\nATATATATAT\n"
UNIVERSE_BED = "chr1\t12\t18\nchr1\t45\t60\nchr2\t100\t120\n"
PARTIAL_BED = "chr1\t15\t50\nchr1\t13\t14\n"
GFF3 = "\n".join(
    [
        "##gff-version 3",
        "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=GENE1",
        "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1;Name=GENE1-201",
        "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=tx1",
        "chr1\tsrc\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=tx1",
        "chr1\tsrc\tgene\t800\t900\t.\t-\t.\tID=gene2;Name=GENE2",
        "chr1\tsrc\tmRNA\t800\t900\t.\t-\t.\tID=tx2;Parent=gene2",
        "chr1\tsrc\texon\t800\t900\t.\t-\t.\tID=exon3;Parent=tx2",
        "",
    ]
)
GTF = "\n".join(
    [
        'chr1\tens\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1"; '
        'gene_biotype "protein_coding";',
        'chr1\tens\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; '
        'gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; '
        'transcript_version "2";',
        'chr1\tens\texon\t100\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; '
        'exon_number "1"; gene_name "GENE1"; gene_biotype "protein_coding"; '
        'transcript_biotype "protein_coding";',
        'chr1\tens\texon\t300\t500\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; '
        'exon_number "2"; gene_name "GENE1"; gene_biotype "protein_coding"; '
        'transcript_biotype "protein_coding";',
        'chr1\tens\tgene\t800\t900\t.\t-\t.\tgene_id "ENSG2"; gene_name "GENE2"; '
        'gene_biotype "lncRNA";',
        'chr1\tens\ttranscript\t800\t900\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; '
        'gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";',
        'chr1\tens\texon\t800\t900\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; '
        'exon_number "1"; gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";',
        "",
    ]
)
GENCODE_GTF = (
    'chr1\thavana\tgene\t10\t90\t.\t+\t.\tgene_id "G1"; gene_name "X"; gene_type "lncRNA";\n'
)

requires_bedtools = pytest.mark.skipif(
    shutil.which("bedtools") is None, reason="the bedtools executable is not installed"
)


def _write(tmp_path: Path, text: str, name: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _intervals(records: list[dict]) -> list[tuple]:
    return [(row["chrom"], row["start"], row["end"]) for row in records]


def _write_bigwig(path: Path) -> Path:
    pybigwig = pytest.importorskip("pyBigWig")
    handle = pybigwig.open(str(path), "w")
    handle.addHeader([("chr1", 100), ("chr2", 50)])
    handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
    handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
    handle.close()
    return path


# bioframe


def test_merge_bioframe_intervals_joins_overlaps_and_matches_upstream(tmp_path: Path) -> None:
    bioframe = pytest.importorskip("bioframe")
    path = _write(tmp_path, A_BED, "a.bed")

    result = backends.merge_bioframe_intervals(path, verbose=False)

    # chr1 [10,20) and [15,30) overlap into [10,30); [40,50) and chr2 [5,9) stay separate.
    assert _intervals(result["intervals"]) == [("chr1", 10, 30), ("chr1", 40, 50), ("chr2", 5, 9)]
    assert [row["n_intervals"] for row in result["intervals"]] == [2, 1, 1]
    upstream = bioframe.merge(bioframe.read_table(path, schema="bed6"), min_dist=0)
    assert result["interval_count"] == len(upstream)


def test_merge_bioframe_intervals_joins_bookended_at_zero_distance(tmp_path: Path) -> None:
    pytest.importorskip("bioframe")
    path = _write(tmp_path, TOUCH_BED, "touch.bed")

    result = backends.merge_bioframe_intervals(path, min_distance=0)

    assert _intervals(result["intervals"]) == [("chr1", 0, 20), ("chr1", 30, 40)]


def test_closest_bioframe_intervals_reports_zero_and_gap_distances(tmp_path: Path) -> None:
    pytest.importorskip("bioframe")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.closest_bioframe_intervals(left, right, k=1)

    distances = {(row["chrom"], row["start"]): row["distance"] for row in result["records"]}
    assert distances[("chr1", 10)] == 0  # [10,20) overlaps [12,18)
    assert distances[("chr1", 40)] == 0  # [40,50) overlaps [45,60)
    assert distances[("chr2", 5)] == 91  # 100 - 9 for separated half-open intervals
    assert result["record_count"] == 4


def test_coverage_bioframe_intervals_counts_covered_bases(tmp_path: Path) -> None:
    pytest.importorskip("bioframe")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.coverage_bioframe_intervals(left, right)

    # [10,20) covers 6 bases of [12,18); [15,30) covers 3; [40,50) covers 5; chr2 covers 0.
    assert [row["coverage"] for row in result["records"]] == [6, 3, 5, 0]


# pybedtools


@requires_bedtools
def test_merge_pybedtools_intervals_follows_distance(tmp_path: Path) -> None:
    pytest.importorskip("pybedtools")
    path = _write(tmp_path, A_BED, "a.bed")

    touching = backends.merge_pybedtools_intervals(path, distance=0)
    gapped = backends.merge_pybedtools_intervals(path, distance=10)

    assert _intervals(touching["intervals"]) == [
        ("chr1", 10, 30),
        ("chr1", 40, 50),
        ("chr2", 5, 9),
    ]
    # A gap of exactly 10 bases between [30,40) and [40,50) joins at distance 10.
    assert _intervals(gapped["intervals"]) == [("chr1", 10, 50), ("chr2", 5, 9)]


@requires_bedtools
def test_closest_pybedtools_intervals_matches_bedtools_distance(tmp_path: Path) -> None:
    pybedtools = pytest.importorskip("pybedtools")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.closest_pybedtools_intervals(left, right, k=1)

    direct = (
        pybedtools.BedTool(str(left))
        .sort()
        .closest(pybedtools.BedTool(str(right)).sort(), d=True, t="first")
    )
    direct_distances = [int(line.split("\t")[-1]) for line in str(direct).splitlines() if line]
    assert [row["distance"] for row in result["records"]] == direct_distances
    # BEDTools counts one extra base for separated intervals: 100 - 9 + 1.
    assert result["records"][-1]["distance"] == 92


@requires_bedtools
def test_coverage_pybedtools_intervals_reports_bases_and_fraction(tmp_path: Path) -> None:
    pytest.importorskip("pybedtools")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    records = backends.coverage_pybedtools_intervals(left, right)["records"]

    first = records[0]
    assert (first["overlaps"], first["bases_covered"], first["length"]) == (1, 6, 10)
    assert first["fraction"] == pytest.approx(0.6)
    assert [row["bases_covered"] for row in records] == [6, 3, 5, 0]


@requires_bedtools
def test_extract_pybedtools_sequences_matches_fasta_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("pybedtools")
    bed = _write(tmp_path, SEQ_BED, "seqs.bed")
    fasta = _write(tmp_path, FASTA, "ref.fa")

    result = backends.extract_pybedtools_sequences(bed, fasta)

    # chr1 is "ACGTACGTAC" + "GGGGTTTTCC", so [10,20) is GGGGTTTTCC; chr2 [5,9) is TATA.
    assert [row["sequence"] for row in result["sequences"]] == ["GGGGTTTTCC", "ACGTA", "TATA"]
    assert result["sequences"][0]["name"] == "chr1:10-20"


@requires_bedtools
def test_extract_pybedtools_sequences_reverse_complements_minus_strand(tmp_path: Path) -> None:
    pytest.importorskip("pybedtools")
    bed = _write(tmp_path, SEQ_BED, "seqs.bed")
    fasta = _write(tmp_path, FASTA, "ref.fa")

    result = backends.extract_pybedtools_sequences(bed, fasta, use_strand=True)

    assert result["sequences"][1]["sequence"] == "TACGT"  # reverse complement of ACGTA


@requires_bedtools
def test_extract_pybedtools_sequences_rejects_intervals_outside_fasta(tmp_path: Path) -> None:
    pytest.importorskip("pybedtools")
    bed = _write(tmp_path, A_BED, "a.bed")  # chr1 intervals run past the 20 bp sequence
    fasta = _write(tmp_path, FASTA, "ref.fa")

    with pytest.raises(ConfigurationError, match="fewer sequences than intervals"):
        backends.extract_pybedtools_sequences(bed, fasta)


# gtfparse


def test_read_gtfparse_table_expands_attributes(tmp_path: Path) -> None:
    pytest.importorskip("gtfparse")
    path = _write(tmp_path, GTF, "genes.gtf")

    result = backends.read_gtfparse_table(path)

    assert result["record_count"] == 7
    assert {"gene_name", "transcript_id", "transcript_version"} <= set(result["columns"])
    transcript = next(row for row in result["records"] if row["feature"] == "transcript")
    assert transcript["transcript_id"] == "ENST1"
    assert transcript["transcript_version"] == 2
    gene = next(row for row in result["records"] if row["feature"] == "gene")
    assert gene["transcript_id"] is None  # absent attributes are reported as None


def test_read_gtfparse_table_filters_features_and_normalizes_biotypes(tmp_path: Path) -> None:
    pytest.importorskip("gtfparse")
    exons = backends.read_gtfparse_table(_write(tmp_path, GTF, "genes.gtf"), features=["exon"])
    assert exons["record_count"] == 3
    assert {row["feature"] for row in exons["records"]} == {"exon"}

    gencode = _write(tmp_path, GENCODE_GTF, "gencode.gtf")
    raw = backends.read_gtfparse_table(gencode)["records"][0]
    normalized = backends.read_gtfparse_table(gencode, normalize_gencode_biotypes=True)
    assert "gene_type" in raw and "gene_biotype" not in raw
    assert normalized["records"][0]["gene_biotype"] == "lncRNA"
    assert normalized["provenance"]["gencode_biotype_aliases"] is True


def test_select_gtfparse_records_by_gene_and_transcript(tmp_path: Path) -> None:
    pytest.importorskip("gtfparse")
    path = _write(tmp_path, GTF, "genes.gtf")

    by_name = backends.select_gtfparse_records(path, gene_names=["GENE2"])
    by_transcript = backends.select_gtfparse_records(path, transcript_ids=["ENST1"])
    by_gene = backends.select_gtfparse_records(path, gene_ids=["ENSG1"])

    assert by_name["record_count"] == 3  # gene, transcript, and one exon of GENE2
    assert by_transcript["record_count"] == 3  # transcript and its two exons
    assert by_gene["record_count"] == 4  # gene, transcript, and two exons of ENSG1
    assert by_gene["gene_ids"] == ["ENSG1"]
    assert by_gene["transcript_ids"] == ["ENST1"]


def test_summarize_gtfparse_annotation_counts_features_and_biotypes(tmp_path: Path) -> None:
    pytest.importorskip("gtfparse")
    result = backends.summarize_gtfparse_annotation(_write(tmp_path, GTF, "genes.gtf"))

    assert result["feature_counts"] == {"gene": 2, "transcript": 2, "exon": 3}
    assert result["gene_count"] == 2
    assert result["transcript_count"] == 2
    assert result["gene_biotype_counts"] == {"protein_coding": 1, "lncRNA": 1}
    assert result["chromosome_counts"] == {"chr1": 7}


# gffutils


def test_create_gffutils_database_counts_features(tmp_path: Path) -> None:
    pytest.importorskip("gffutils")
    source = _write(tmp_path, GFF3, "genes.gff3")
    target = tmp_path / "genes.db"

    result = backends.create_gffutils_database(source, target)

    assert result["input_format"] == "gff3"
    assert result["featuretype_counts"] == {"exon": 3, "gene": 2, "mRNA": 2}
    assert result["feature_count"] == 7
    assert Path(result["database"]) == target.resolve()


def test_create_gffutils_database_refuses_existing_path(tmp_path: Path) -> None:
    pytest.importorskip("gffutils")
    source = _write(tmp_path, GFF3, "genes.gff3")
    existing = _write(tmp_path, "not a database", "existing.db")

    with pytest.raises(ConfigurationError, match="already exists"):
        backends.create_gffutils_database(source, existing)


def test_list_gffutils_children_by_level_and_featuretype(tmp_path: Path) -> None:
    pytest.importorskip("gffutils")
    database = tmp_path / "genes.db"
    backends.create_gffutils_database(_write(tmp_path, GFF3, "genes.gff3"), database)

    direct = backends.list_gffutils_children(database, "gene1", level=1)
    grandchildren = backends.list_gffutils_children(database, "gene1", level=2)
    exons = backends.list_gffutils_children(database, "tx1", featuretype="exon")

    assert [row["id"] for row in direct["records"]] == ["tx1"]
    assert [row["id"] for row in grandchildren["records"]] == ["exon1", "exon2"]
    assert [(row["start"], row["end"]) for row in exons["records"]] == [(100, 200), (300, 500)]


def test_list_gffutils_children_drops_gtf_self_links(tmp_path: Path) -> None:
    gffutils = pytest.importorskip("gffutils")
    database = tmp_path / "genes_gtf.db"
    backends.create_gffutils_database(_write(tmp_path, GTF, "genes.gtf"), database)

    # gffutils stores a GTF-inferred transcript as its own level-1 child; the wrapper drops it.
    raw = [feature.id for feature in gffutils.FeatureDB(str(database)).children("ENST1", level=1)]
    result = backends.list_gffutils_children(database, "ENST1")

    assert "ENST1" in raw
    assert [row["id"] for row in result["records"]] == ["exon_1", "exon_2"]


def test_count_gffutils_featuretypes_matches_upstream(tmp_path: Path) -> None:
    gffutils = pytest.importorskip("gffutils")
    database = tmp_path / "genes.db"
    backends.create_gffutils_database(_write(tmp_path, GFF3, "genes.gff3"), database)

    result = backends.count_gffutils_featuretypes(database)

    upstream = gffutils.FeatureDB(str(database))
    expected = {kind: upstream.count_features_of_type(kind) for kind in upstream.featuretypes()}
    assert result["featuretype_counts"] == expected
    assert result["feature_count"] == 7


def test_list_gffutils_children_rejects_unknown_parent_and_nul(tmp_path: Path) -> None:
    pytest.importorskip("gffutils")
    database = tmp_path / "genes.db"
    backends.create_gffutils_database(_write(tmp_path, GFF3, "genes.gff3"), database)

    with pytest.raises(ConfigurationError, match="is not in"):
        backends.list_gffutils_children(database, "missing-feature")
    with pytest.raises(ConfigurationError, match="NUL"):
        backends.list_gffutils_children(database, "gene\x001")


# pyBigWig


def test_read_pybigwig_values_per_base_with_missing_as_none(tmp_path: Path) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = _write_bigwig(tmp_path / "signal.bw")

    result = backends.read_pybigwig_values(bigwig, "chr1", 8, 30)

    # Bases 10-14 hold 1.0 and bases 20-24 hold 3.0; the other 12 bases of [8,30) are missing.
    values = result["values"]
    assert result["length"] == 22
    assert result["covered_bases"] == 10
    assert values[:2] == [None, None]
    assert values[2:7] == [1.0] * 5
    assert sum(value for value in values if value is not None) == pytest.approx(20.0)


def test_read_pybigwig_intervals_returns_stored_intervals(tmp_path: Path) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = _write_bigwig(tmp_path / "signal.bw")

    result = backends.read_pybigwig_intervals(bigwig, "chr1", 8, 30)

    assert [(row["start"], row["end"], row["value"]) for row in result["records"]] == [
        (10, 15, 1.0),
        (20, 25, 3.0),
    ]


def test_summarize_pybigwig_bed_intervals_statistics(tmp_path: Path) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = _write_bigwig(tmp_path / "signal.bw")
    bed = _write(tmp_path, SEQ_BED, "seqs.bed")

    mean = backends.summarize_pybigwig_bed_intervals(bigwig, bed, statistic="mean")
    total = backends.summarize_pybigwig_bed_intervals(bigwig, bed, statistic="sum")

    values = {row["name"]: row["value"] for row in mean["records"]}
    assert values == {"p1": 1.0, "p2": None, "p3": 2.0}  # p2 lies where the track has no data
    sums = {row["name"]: row["value"] for row in total["records"]}
    assert sums["p1"] == pytest.approx(5.0)  # five bases of value 1.0
    assert sums["p3"] == pytest.approx(8.0)  # four bases of value 2.0


def test_pybigwig_rejects_bad_regions_and_statistics(tmp_path: Path) -> None:
    pytest.importorskip("pyBigWig")
    bigwig = _write_bigwig(tmp_path / "signal.bw")
    bed = _write(tmp_path, SEQ_BED, "seqs.bed")

    with pytest.raises(ConfigurationError, match="no chromosome"):
        backends.read_pybigwig_values(bigwig, "chrX", 0, 10)
    with pytest.raises(ConfigurationError, match="must satisfy 0 <= start < end"):
        backends.read_pybigwig_values(bigwig, "chr1", 0, 1000)
    with pytest.raises(ConfigurationError, match="max_bases"):
        backends.read_pybigwig_values(bigwig, "chr1", 0, 50, max_bases=10)
    with pytest.raises(ConfigurationError, match="statistic must be one of"):
        backends.summarize_pybigwig_bed_intervals(bigwig, bed, statistic="median")


# pyranges1


def test_intersect_pyranges1_intervals_clips_to_overlap(tmp_path: Path) -> None:
    pytest.importorskip("pyranges1")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.intersect_pyranges1_intervals(left, right)

    assert [(row["Start"], row["End"]) for row in result["records"]] == [
        (12, 18),
        (15, 18),
        (45, 50),
    ]


def test_nearest_pyranges1_intervals_uses_bedtools_distance_convention(tmp_path: Path) -> None:
    pytest.importorskip("pyranges1")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.nearest_pyranges1_intervals(left, right, k=1)

    distances = [row["Distance"] for row in result["records"]]
    assert distances == [0, 0, 0, 92]  # 100 - 9 + 1 for the chr2 pair


def test_complement_pyranges1_intervals_returns_gaps(tmp_path: Path) -> None:
    pytest.importorskip("pyranges1")
    path = _write(tmp_path, A_BED, "a.bed")

    result = backends.complement_pyranges1_intervals(path, chrom_sizes={"chr1": 100, "chr2": 50})

    assert [(row["Chromosome"], row["Start"], row["End"]) for row in result["records"]] == [
        ("chr1", 0, 10),
        ("chr1", 30, 40),
        ("chr1", 50, 100),
        ("chr2", 0, 5),
        ("chr2", 9, 50),
    ]


# polars-bio


def test_merge_polars_bio_intervals_respects_min_distance(tmp_path: Path) -> None:
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, A_BED, "a.bed")
    touch = _write(tmp_path, TOUCH_BED, "touch.bed")

    overlaps = backends.merge_polars_bio_intervals(path, min_distance=0)
    gapped = backends.merge_polars_bio_intervals(path, min_distance=10)
    bookended = backends.merge_polars_bio_intervals(touch, min_distance=0)

    assert _intervals(overlaps["intervals"]) == [
        ("chr1", 10, 30),
        ("chr1", 40, 50),
        ("chr2", 5, 9),
    ]
    assert _intervals(gapped["intervals"]) == [("chr1", 10, 50), ("chr2", 5, 9)]
    assert _intervals(bookended["intervals"]) == [("chr1", 0, 20), ("chr1", 30, 40)]


def test_nearest_polars_bio_intervals_zero_based_distance(tmp_path: Path) -> None:
    pytest.importorskip("polars_bio")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.nearest_polars_bio_intervals(left, right, k=1)

    assert [row["distance"] for row in result["records"]] == [0, 0, 0, 91]  # 100 - 9


def test_coverage_polars_bio_intervals_counts_covered_bases(tmp_path: Path) -> None:
    pytest.importorskip("polars_bio")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.coverage_polars_bio_intervals(left, right)

    assert [row["coverage"] for row in result["records"]] == [6, 3, 5, 0]


# gtars


def test_jaccard_gtars_regions_base_pair_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("gtars")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.jaccard_gtars_regions(left, right)

    # Covered bases: A = 20 + 10 + 4 = 34, B = 6 + 15 + 20 = 41, shared = 6 + 5 = 11.
    # Jaccard = 11 / (34 + 41 - 11) = 11 / 64.
    assert result["jaccard"] == pytest.approx(11 / 64)
    assert (result["source_regions"], result["other_regions"]) == (4, 3)


def test_gaps_gtars_regions_includes_regionless_chromosomes(tmp_path: Path) -> None:
    pytest.importorskip("gtars")
    path = _write(tmp_path, A_BED, "a.bed")

    result = backends.gaps_gtars_regions(path, {"chr1": 100, "chr2": 50, "chr3": 7})

    assert sorted(_intervals(result["gaps"])) == [
        ("chr1", 0, 10),
        ("chr1", 30, 40),
        ("chr1", 50, 100),
        ("chr2", 0, 5),
        ("chr2", 9, 50),
        ("chr3", 0, 7),
    ]


def test_tokenize_gtars_regions_keeps_one_token_per_overlap(tmp_path: Path) -> None:
    pytest.importorskip("gtars")
    partial = _write(tmp_path, PARTIAL_BED, "partial.bed")
    universe = _write(tmp_path, UNIVERSE_BED, "universe.bed")

    result = backends.tokenize_gtars_regions(partial, universe)

    # chr1:15-50 overlaps two universe regions and chr1:13-14 overlaps one.
    assert Counter(result["tokens"]) == Counter(["chr1:12-18", "chr1:12-18", "chr1:45-60"])
    assert result["token_count"] == 3
    assert len(result["token_ids"]) == 3


def test_tokenize_gtars_regions_requires_bed3_universe(tmp_path: Path) -> None:
    pytest.importorskip("gtars")
    partial = _write(tmp_path, PARTIAL_BED, "partial.bed")
    wide = _write(tmp_path, A_BED, "universe.bed")  # BED6, not BED3

    with pytest.raises(ConfigurationError, match="BED3"):
        backends.tokenize_gtars_regions(partial, wide)


# ncls


def _brute_pairs(
    targets: list[tuple[str, int, int]],
    queries: list[tuple[str, int, int]],
    containment: bool,
) -> list[tuple[int, int]]:
    pairs = []
    for query_index, (q_chrom, q_start, q_end) in enumerate(queries):
        for target_index, (t_chrom, t_start, t_end) in enumerate(targets):
            if t_chrom != q_chrom:
                continue
            if containment:
                hit = t_start <= q_start and q_end <= t_end
            else:
                hit = t_start < q_end and q_start < t_end
            if hit:
                pairs.append((query_index, target_index))
    return sorted(pairs)


def _random_bed(rng: random.Random, count: int) -> tuple[str, list[tuple[str, int, int]]]:
    rows = []
    for _ in range(count):
        chrom = rng.choice(["chr1", "chr2"])
        start = rng.randrange(0, 500)
        rows.append((chrom, start, start + rng.randrange(1, 60)))
    text = "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in rows)
    return text, rows


def test_ncls_overlaps_and_containment_match_brute_force(tmp_path: Path) -> None:
    pytest.importorskip("ncls")
    rng = random.Random(7)
    target_text, targets = _random_bed(rng, 60)
    query_text, queries = _random_bed(rng, 40)
    target_path = _write(tmp_path, target_text, "targets.bed")
    query_path = _write(tmp_path, query_text, "queries.bed")

    overlaps = backends.query_ncls_overlaps(target_path, query_path)
    containing = backends.find_ncls_containing_intervals(target_path, query_path)

    found = sorted((row["query_index"], row["target_index"]) for row in overlaps["records"])
    assert found == _brute_pairs(targets, queries, containment=False)
    found_containing = sorted(
        (row["query_index"], row["target_index"]) for row in containing["records"]
    )
    assert found_containing == _brute_pairs(targets, queries, containment=True)


def test_count_ncls_overlaps_matches_hand_counts(tmp_path: Path) -> None:
    pytest.importorskip("ncls")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.count_ncls_overlaps(left, right)

    # chr1 [12,18) overlaps two source intervals, chr1 [45,60) one, chr2 [100,120) none.
    assert [row["overlap_count"] for row in result["records"]] == [2, 1, 0]


def test_ncls_containment_requires_full_containment(tmp_path: Path) -> None:
    pytest.importorskip("ncls")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.find_ncls_containing_intervals(left, right)

    # [12,18) lies inside [10,20) but not inside [15,30).
    assert [(row["query_index"], row["target_index"]) for row in result["records"]] == [(0, 0)]


# genomicranges


def test_reduce_genomicranges_intervals_merges_overlaps(tmp_path: Path) -> None:
    pytest.importorskip("genomicranges")
    pytest.importorskip("pandas")
    path = _write(tmp_path, A_BED, "a.bed")

    result = backends.reduce_genomicranges_intervals(path)

    assert _intervals(result["intervals"]) == [("chr1", 10, 30), ("chr1", 40, 50), ("chr2", 5, 9)]


def test_reduce_genomicranges_min_gap_width_controls_bookended_ranges(tmp_path: Path) -> None:
    pytest.importorskip("genomicranges")
    pytest.importorskip("pandas")
    path = _write(tmp_path, TOUCH_BED, "touch.bed")

    joined = backends.reduce_genomicranges_intervals(path, min_gap_width=1)
    separate = backends.reduce_genomicranges_intervals(path, min_gap_width=0)

    assert _intervals(joined["intervals"]) == [("chr1", 0, 20), ("chr1", 30, 40)]
    assert _intervals(separate["intervals"]) == [
        ("chr1", 0, 10),
        ("chr1", 10, 20),
        ("chr1", 30, 40),
    ]


def test_intersect_genomicranges_intervals_keeps_shared_positions(tmp_path: Path) -> None:
    pytest.importorskip("genomicranges")
    pytest.importorskip("pandas")
    left = _write(tmp_path, A_BED, "a.bed")
    right = _write(tmp_path, B_BED, "b.bed")

    result = backends.intersect_genomicranges_intervals(left, right)

    # A covers [10,30) and [40,50) on chr1, B covers [12,18) and [45,60); shared: [12,18), [45,50).
    assert _intervals(result["intervals"]) == [("chr1", 12, 18), ("chr1", 45, 50)]


def test_flank_genomicranges_intervals_start_and_end_sides(tmp_path: Path) -> None:
    pytest.importorskip("genomicranges")
    pytest.importorskip("pandas")
    path = _write(tmp_path, A_BED, "a.bed")

    upstream = backends.flank_genomicranges_intervals(path, 5, side="start")
    downstream = backends.flank_genomicranges_intervals(path, 5, side="end")

    assert _intervals(upstream["intervals"]) == [
        ("chr1", 5, 10),
        ("chr1", 10, 15),
        ("chr1", 35, 40),
        ("chr2", 0, 5),
    ]
    assert _intervals(downstream["intervals"]) == [
        ("chr1", 20, 25),
        ("chr1", 30, 35),
        ("chr1", 50, 55),
        ("chr2", 9, 14),
    ]


def test_flank_genomicranges_intervals_rejects_flanks_before_position_zero(
    tmp_path: Path,
) -> None:
    pytest.importorskip("genomicranges")
    pytest.importorskip("pandas")
    path = _write(tmp_path, A_BED, "a.bed")

    with pytest.raises(ConfigurationError, match="before position 0"):
        backends.flank_genomicranges_intervals(path, 50)


# Validation and module contract (no upstream package needed).


def test_public_names_are_exported_and_callable() -> None:
    assert len(backends.__all__) == len(set(backends.__all__))
    for name in backends.__all__:
        assert callable(getattr(backends, name)), name


def test_integer_options_reject_bools_and_out_of_range(tmp_path: Path) -> None:
    path = _write(tmp_path, A_BED, "a.bed")

    with pytest.raises(ConfigurationError, match="min_distance"):
        backends.merge_bioframe_intervals(path, min_distance=True)
    with pytest.raises(ConfigurationError, match="k must be an integer"):
        backends.closest_pybedtools_intervals(path, path, k=0)
    with pytest.raises(ConfigurationError, match="max_rows"):
        backends.merge_pybedtools_intervals(path, max_rows=2)


def test_bed_reader_rejects_malformed_and_empty_files(tmp_path: Path) -> None:
    empty = _write(tmp_path, "# only a comment\n", "empty.bed")
    reversed_bed = _write(tmp_path, "chr1\t20\t10\n", "reversed.bed")
    narrow = _write(tmp_path, "chr1\t10\n", "narrow.bed")

    with pytest.raises(ConfigurationError, match="no BED intervals"):
        backends.merge_pybedtools_intervals(empty)
    with pytest.raises(ConfigurationError, match="0 <= start < end"):
        backends.summarize_pybigwig_bed_intervals(empty, reversed_bed)
    with pytest.raises(ConfigurationError, match="BED3 to BED6"):
        backends.merge_pybedtools_intervals(narrow)


def test_select_gtfparse_records_requires_a_filter(tmp_path: Path) -> None:
    path = _write(tmp_path, GTF, "genes.gtf")

    with pytest.raises(ConfigurationError, match="at least one"):
        backends.select_gtfparse_records(path)


def test_gtf_reader_enforces_max_rows(tmp_path: Path) -> None:
    path = _write(tmp_path, GTF, "genes.gtf")

    with pytest.raises(ConfigurationError, match="max_rows=3"):
        backends.read_gtfparse_table(path, max_rows=3)


@pytest.mark.parametrize(
    ("module", "extra", "call"),
    [
        (
            "bioframe",
            "bio-bioframe",
            lambda p: backends.merge_bioframe_intervals(p),
        ),
        (
            "pybedtools",
            "bio-pybedtools",
            lambda p: backends.merge_pybedtools_intervals(p),
        ),
        (
            "gtfparse",
            "bio-gtfparse",
            lambda p: backends.read_gtfparse_table(p),
        ),
        (
            "gffutils",
            "bio-gffutils",
            lambda p: backends.count_gffutils_featuretypes(p),
        ),
        (
            "pyBigWig",
            "bio-pybigwig",
            lambda p: backends.read_pybigwig_values(p, "chr1", 0, 5),
        ),
        (
            "pyranges1",
            "bio-pyranges1",
            lambda p: backends.intersect_pyranges1_intervals(p, p),
        ),
        (
            "polars_bio",
            "bio-polarsbio",
            lambda p: backends.merge_polars_bio_intervals(p),
        ),
        (
            "gtars.models",
            "bio-gtars",
            lambda p: backends.jaccard_gtars_regions(p, p),
        ),
        (
            "ncls",
            "bio-ncls",
            lambda p: backends.count_ncls_overlaps(p, p),
        ),
        (
            "genomicranges",
            "bio-genomicranges",
            lambda p: backends.reduce_genomicranges_intervals(p),
        ),
    ],
)
def test_missing_optional_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, module: str, extra: str, call
) -> None:
    path = _write(tmp_path, A_BED, "a.bed")
    monkeypatch.setitem(sys.modules, module, None)

    with pytest.raises(BackendUnavailableError, match=extra):
        call(path)
