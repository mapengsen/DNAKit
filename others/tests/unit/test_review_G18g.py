"""Upstream-checked and hand-calculated wrappers for pybedtools, gffutils, gtfparse (G18g)."""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

from dnakit.annotation import gffutils_review_functions as gu_api
from dnakit.annotation import gtfparse_review_functions as gt_api
from dnakit.annotation import pybedtools_review_functions as pb_api
from dnakit.exceptions import ConfigurationError

ALL_WRAPPERS = [
    pb_api.bed_to_gff_pybedtools,
    pb_api.center_pybedtools_coordinates,
    pb_api.cis_trans_pybedtools_interactions,
    pb_api.filter_pybedtools_intervals_by_length,
    pb_api.gff_to_bed_pybedtools,
    pb_api.normalize_pybedtools_scores_to_length,
    pb_api.overlap_pybedtools_bp,
    pb_api.random_jaccard_pybedtools_intervals,
    pb_api.scale_pybedtools_bedgraph_values,
    pb_api.three_prime_pybedtools_coordinates,
    gu_api.asinterval_gffutils_feature,
    gu_api.bin_gffutils_interval,
    gu_api.canonical_gffutils_transcript_sequences,
    gu_api.features_gffutils_to_bedtool,
    gu_api.gffutils_line_to_seqfeature,
    gu_api.infer_gffutils_attribute_dialect,
    gu_api.merge_gffutils_features,
    gu_api.read_gffutils_records,
    gu_api.sanitize_gffutils_file,
    gu_api.seqfeature_to_gffutils_line,
    gt_api.create_missing_gtfparse_features,
    gt_api.expand_gtfparse_attribute_strings,
]

SMALL_GFF3 = """##gff-version 3
chr1\tens\tgene\t101\t400\t.\t+\t.\tID=g1;Name=GENE1;gene_biotype=protein_coding
chr1\tens\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1;Name=T1
chr1\tens\texon\t101\t200\t.\t+\t.\tID=e1;Parent=t1
chr1\tens\texon\t301\t400\t.\t+\t.\tID=e2;Parent=t1
chr1\tens\tCDS\t151\t200\t.\t+\t0\tID=c1;Parent=t1
chr1\tens\tgene\t501\t700\t.\t-\t.\tID=g2;Name=GENE2
chr1\tens\tmRNA\t501\t700\t.\t-\t.\tID=t2;Parent=g2
chr1\tens\texon\t501\t700\t.\t-\t.\tID=e3;Parent=t2
chr2\tens\tgene\t11\t40\t.\t+\t.\tID=g3;Name=GENE3
chr2\tens\tmRNA\t11\t40\t.\t+\t.\tID=t3;Parent=g3
chr2\tens\texon\t11\t40\t.\t+\t.\tID=e4;Parent=t3
"""
EXONS_GTF = (
    'chr1\tHAVANA\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENEA";\n'
    'chr1\tHAVANA\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENEA";\n'
    'chr1\tHAVANA\texon\t501\t600\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENEB";\n'
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _pattern(length: int) -> str:
    """Reference bases with base i (0-based) equal to 'ACGT'[i % 4]."""

    return ("ACGT" * (length // 4 + 1))[:length]


def _revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fasta(tmp_path: Path) -> Path:
    text = f">chr1\n{_pattern(1000)}\n>chr2\n{'G' * 100}\n"
    return _write(tmp_path / "ref.fa", text)


def test_every_wrapper_has_native_keyword_only_verbose() -> None:
    for function in ALL_WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_module_exports_match_wrappers() -> None:
    assert set(pb_api.__all__) == {
        f.__name__ for f in ALL_WRAPPERS if f.__module__.endswith("pybedtools_review_functions")
    }
    assert set(gu_api.__all__) == {
        f.__name__ for f in ALL_WRAPPERS if f.__module__.endswith("gffutils_review_functions")
    }
    assert set(gt_api.__all__) == {
        f.__name__ for f in ALL_WRAPPERS if f.__module__.endswith("gtfparse_review_functions")
    }


# pybedtools


def test_overlap_matches_hand_calculation_and_upstream_argument_order() -> None:
    pybedtools = importlib.import_module("pybedtools")
    result = pb_api.overlap_pybedtools_bp(0, 10, 5, 15, verbose=False)
    assert result["overlap_bp"] == 5 == pybedtools.overlap(0, 5, 10, 15)
    assert result["overlaps"] is True
    touching = pb_api.overlap_pybedtools_bp(0, 10, 10, 20, verbose=False)
    assert touching["overlap_bp"] == 0 and touching["overlaps"] is False
    gap = pb_api.overlap_pybedtools_bp(0, 10, 20, 30, verbose=False)
    assert gap["overlap_bp"] == -10 and gap["overlaps"] is False
    with pytest.raises(ConfigurationError):
        pb_api.overlap_pybedtools_bp(10, 10, 0, 5, verbose=False)


def test_normalize_scores_matches_hand_value_and_rejects_non_numeric(tmp_path: Path) -> None:
    bedgraph = _write(tmp_path / "bg.bedgraph", "chr1\t0\t100\t10\nchr1\t100\t200\t20\n")
    result = pb_api.normalize_pybedtools_scores_to_length(
        bedgraph, index=3, scalar=1000, verbose=False
    )
    # value * scalar / length: 10 * 1000 / 100 = 100 and 20 * 1000 / 100 = 200
    assert [r["name"] for r in result["records"]] == ["100.0", "200.0"]
    assert result["record_count"] == 2
    bed3 = _write(tmp_path / "b3.bed", "chr1\t0\t100\n")
    with pytest.raises(ConfigurationError):
        pb_api.normalize_pybedtools_scores_to_length(bed3, index=4, verbose=False)
    with pytest.raises(ConfigurationError):
        pb_api.normalize_pybedtools_scores_to_length(bedgraph, index=0, verbose=False)


def test_scale_bedgraph_values_matches_hand_value(tmp_path: Path) -> None:
    bedgraph = _write(tmp_path / "bg.bedgraph", "chr1\t0\t100\t10\nchr1\t100\t200\t20\n")
    result = pb_api.scale_pybedtools_bedgraph_values(bedgraph, scalar=0.5, verbose=False)
    assert [r["name"] for r in result["records"]] == ["5.0", "10.0"]
    with pytest.raises(ConfigurationError):
        pb_api.scale_pybedtools_bedgraph_values(
            _write(tmp_path / "x.bed", "chr1\t0\t100\tp\t5\n"), scalar=2, verbose=False
        )


def test_center_windows_match_hand_calculation(tmp_path: Path) -> None:
    peaks = _write(
        tmp_path / "peaks.bed",
        "chr1\t100\t200\tpeakA\t5\t+\nchr1\t300\t340\tpeakB\t9\t-\nchr2\t50\t60\tpeakC\t2\t+\n",
    )
    forty = pb_api.center_pybedtools_coordinates(peaks, width=40, verbose=False)
    by_name = {r["name"]: (r["start"], r["end"]) for r in forty["records"]}
    # midpoint of 100-200 is 150; a 40 bp window is 130-170; the 10 bp feature is unchanged
    assert by_name == {"peakA": (130, 170), "peakB": (300, 340), "peakC": (50, 60)}
    single = pb_api.center_pybedtools_coordinates(peaks, width=1, verbose=False)
    assert {r["name"]: (r["start"], r["end"]) for r in single["records"]}["peakA"] == (150, 151)
    odd = pb_api.center_pybedtools_coordinates(peaks, width=41, verbose=False)
    assert {r["name"]: (r["start"], r["end"]) for r in odd["records"]}["peakA"] == (130, 170)


def test_three_prime_windows_follow_strand(tmp_path: Path) -> None:
    peaks = _write(
        tmp_path / "peaks.bed",
        "chr1\t100\t200\tpeakA\t5\t+\nchr1\t300\t340\tpeakB\t9\t-\nchr2\t50\t60\tpeakC\t2\t+\n",
    )
    result = pb_api.three_prime_pybedtools_coordinates(
        peaks, upstream=10, downstream=20, verbose=False
    )
    windows = {r["name"]: (r["start"], r["end"]) for r in result["records"]}
    # + strand 3' end is stop: [200 - 10, 200 + 20); - strand 3' end is start: [300 - 20, 300 + 10)
    assert windows == {"peakA": (190, 220), "peakB": (280, 310), "peakC": (50, 80)}


def test_filter_by_length_uses_strict_thresholds(tmp_path: Path) -> None:
    peaks = _write(
        tmp_path / "peaks.bed",
        "chr1\t100\t200\tpeakA\t5\t+\nchr1\t300\t340\tpeakB\t9\t-\nchr2\t50\t60\tpeakC\t2\t+\n",
    )

    def names(result: dict) -> list[str]:
        return [record["name"] for record in result["records"]]

    assert names(
        pb_api.filter_pybedtools_intervals_by_length(peaks, greater_than=30, verbose=False)
    ) == [
        "peakA",
        "peakB",
    ]
    assert names(
        pb_api.filter_pybedtools_intervals_by_length(peaks, less_than=100, verbose=False)
    ) == [
        "peakB",
        "peakC",
    ]
    both = pb_api.filter_pybedtools_intervals_by_length(
        peaks, greater_than=30, less_than=100, verbose=False
    )
    assert names(both) == ["peakB"]
    with pytest.raises(ConfigurationError):
        pb_api.filter_pybedtools_intervals_by_length(peaks, verbose=False)


def test_gff_to_bed_uses_name_field_and_matches_upstream(tmp_path: Path) -> None:
    pybedtools = importlib.import_module("pybedtools")
    featurefuncs = importlib.import_module("pybedtools.featurefuncs")
    gff = _write(
        tmp_path / "g.gff",
        "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1;Name=GENE1\n"
        "chr1\tsrc\tmRNA\t101\t180\t5\t+\t.\tID=t1;Parent=g1\n",
    )
    named = pb_api.gff_to_bed_pybedtools(gff, name_field="Name", verbose=False)
    assert [(r["start"], r["end"], r["name"]) for r in named["records"]] == [
        (100, 180, "."),
        (100, 200, "GENE1"),
    ]
    upstream = sorted(
        str(x).rstrip("\n")
        for x in pybedtools.BedTool(str(gff)).each(featurefuncs.gff2bed, name_field="Name")
    )
    ours = sorted(
        "\t".join([r["chrom"], str(r["start"]), str(r["end"]), r["name"], r["score"], r["strand"]])
        for r in named["records"]
    )
    assert ours == sorted("\t".join(line.split("\t")[:6]) for line in upstream)
    default = pb_api.gff_to_bed_pybedtools(gff, verbose=False)
    assert sorted(r["name"] for r in default["records"]) == ["g1", "t1"]
    with pytest.raises(ConfigurationError):
        pb_api.gff_to_bed_pybedtools(_write(tmp_path / "b.bed", "chr1\t1\t2\n"), verbose=False)


def test_bed_to_gff_shifts_start_to_one_based(tmp_path: Path) -> None:
    peaks = _write(tmp_path / "peaks.bed", "chr1\t100\t200\tpeakA\t5\t+\n")
    result = pb_api.bed_to_gff_pybedtools(peaks, verbose=False)
    record = result["records"][0]
    assert (record["seqid"], record["start"], record["end"]) == ("chr1", 101, 200)
    assert record["strand"] == "+" and record["attributes"] == 'Name="peakA";'


def test_cis_trans_matches_upstream_table(tmp_path: Path) -> None:
    pybedtools = importlib.import_module("pybedtools")
    long_range = importlib.import_module("pybedtools.contrib.long_range_interaction")
    bedpe = _write(
        tmp_path / "pairs.bedpe",
        "chr1\t1\t10\tchr1\t50\t90\tpair1\t5\t+\t-\nchr1\t2\t15\tchr1\t200\t210\tpair2\t1\t+\t+\n",
    )
    tss = _write(
        tmp_path / "tss.bed", "chr1\t5\t6\tgene1\nchr1\t60\t61\tgene2\nchr1\t88\t89\tgene3\n"
    )
    peaks = _write(tmp_path / "pk.bed", "chr1\t3\t4\tpeak1\t50\t.\n")
    result = pb_api.cis_trans_pybedtools_interactions(
        bedpe, {"tss": tss, "pk": peaks}, verbose=False
    )
    iterator, count, extra = long_range.tag_bedpe(
        pybedtools.BedTool(str(bedpe)),
        {"tss": pybedtools.BedTool(str(tss)), "pk": pybedtools.BedTool(str(peaks))},
        verbose=False,
    )
    frame = long_range.cis_trans_interactions(iterator, count, extra, verbose=False)
    upstream = sorted(
        tuple(str(row[column]) for column in sorted(row)) for row in frame.to_dict("records")
    )
    ours = sorted(tuple(str(row[column]) for column in sorted(row)) for row in result["records"])
    assert ours == upstream
    assert {
        "target_label": "tss",
        "target_name": "gene1",
        "cis_label": "pk",
        "cis_name": "peak1",
        "distal_label": "tss",
        "distal_name": "gene2",
        "label": "pair1",
    } in result["records"]
    with pytest.raises(ConfigurationError):
        pb_api.cis_trans_pybedtools_interactions(bedpe, {"tss": tss}, verbose=False)
    with pytest.raises(ConfigurationError):
        pb_api.cis_trans_pybedtools_interactions(bedpe, {"t s": tss, "pk": peaks}, verbose=False)


def test_random_jaccard_observed_hand_value_and_seeded_upstream_match(tmp_path: Path) -> None:
    stats = importlib.import_module("pybedtools.stats")
    pybedtools = importlib.import_module("pybedtools")
    a = _write(tmp_path / "a.bed", "chr1\t100\t300\tpeakA\n")
    b = _write(tmp_path / "b.bed", "chr1\t150\t350\tb1\nchr1\t600\t700\tb2\n")
    genome = _write(tmp_path / "genome.tsv", "chr1\t1000\nchr2\t500\n")
    result = pb_api.random_jaccard_pybedtools_intervals(a, b, genome, seed=7, verbose=False)
    # observed: intersection [150, 300) = 150 bp; union of both sets = 350 bp
    assert result["observed"]["intersection_bp"] == 150
    assert result["observed"]["union_bp"] == 350
    assert result["observed"]["jaccard"] == pytest.approx(150 / 350, abs=1e-6)
    direct = stats.random_jaccard(
        pybedtools.BedTool(str(a)),
        pybedtools.BedTool(str(b)),
        str(genome),
        {"seed": 7},
        {},
    )
    assert result["shuffled"]["intersection_bp"] == int(direct["intersection"])
    assert result["shuffled"]["union_bp"] == int(direct["union-intersection"])
    again = pb_api.random_jaccard_pybedtools_intervals(a, b, genome, seed=7, verbose=False)
    assert again["shuffled"] == result["shuffled"]
    with pytest.raises(ConfigurationError):
        pb_api.random_jaccard_pybedtools_intervals(
            a, b, _write(tmp_path / "bad.tsv", "chr1\n"), verbose=False
        )


# gffutils


def test_bin_matches_ucsc_hand_values() -> None:
    assert gu_api.bin_gffutils_interval(100, 200, verbose=False)["bins"] == [4681]
    assert gu_api.bin_gffutils_interval(200_000, 200_100, verbose=False)["bins"] == [4682]
    assert gu_api.bin_gffutils_interval(131_000, 131_200, verbose=False)["bins"] == [585]
    assert gu_api.bin_gffutils_interval(0, 100, fmt="bed", verbose=False)["bins"] == [4681]
    many = gu_api.bin_gffutils_interval(1, 70_000, one=False, verbose=False)
    assert many["bins"] == [1, 9, 73, 585, 4681]
    with pytest.raises(ConfigurationError):
        gu_api.bin_gffutils_interval(0, 100, fmt="gff", verbose=False)
    with pytest.raises(ConfigurationError):
        gu_api.bin_gffutils_interval(10, 100, fmt="other", verbose=False)


def test_read_records_matches_datairerator_count(tmp_path: Path) -> None:
    gffutils = importlib.import_module("gffutils")
    gff = _write(tmp_path / "small.gff3", SMALL_GFF3)
    result = gu_api.read_gffutils_records(gff, verbose=False)
    assert result["record_count"] == len(list(gffutils.DataIterator(str(gff)))) == 11
    first = result["records"][0]
    assert first["featuretype"] == "gene" and (first["start"], first["end"]) == (101, 400)
    assert first["attributes"]["Name"] == ["GENE1"]
    with pytest.raises(ConfigurationError):
        gu_api.read_gffutils_records(gff, max_rows=3, verbose=False)


def test_asinterval_matches_hand_coordinates(tmp_path: Path) -> None:
    gff = _write(tmp_path / "small.gff3", SMALL_GFF3)
    result = gu_api.asinterval_gffutils_feature(gff, "g1", verbose=False)
    assert (result["chrom"], result["start"], result["end"]) == ("chr1", 100, 400)
    assert result["name"] == "g1" and result["strand"] == "+"
    with pytest.raises(ConfigurationError):
        gu_api.asinterval_gffutils_feature(gff, "missing", verbose=False)


def test_features_to_bedtool_is_sorted_and_filters_type(tmp_path: Path) -> None:
    gff = _write(tmp_path / "small.gff3", SMALL_GFF3)
    genes = gu_api.features_gffutils_to_bedtool(gff, featuretype="gene", verbose=False)
    assert [(r["chrom"], r["start"], r["end"], r["name"]) for r in genes["records"]] == [
        ("chr1", 100, 400, "g1"),
        ("chr1", 500, 700, "g2"),
        ("chr2", 10, 40, "g3"),
    ]
    everything = gu_api.features_gffutils_to_bedtool(gff, verbose=False)
    assert everything["record_count"] == 11


def test_infer_dialect_recognises_gff3_and_gtf() -> None:
    gff3 = gu_api.infer_gffutils_attribute_dialect("ID=g1;Name=x", verbose=False)
    assert gff3["fmt"] == "gff3" and gff3["dialect"]["keyval separator"] == "="
    gtf = gu_api.infer_gffutils_attribute_dialect('gene_id "A"; transcript_id "B";', verbose=False)
    assert gtf["fmt"] == "gtf" and gtf["dialect"]["order"] == ["gene_id", "transcript_id"]


def test_canonical_transcripts_return_hand_sequences(tmp_path: Path) -> None:
    gff = _write(tmp_path / "small.gff3", SMALL_GFF3)
    fasta = _fasta(tmp_path)
    result = gu_api.canonical_gffutils_transcript_sequences(gff, fasta, verbose=False)
    by_id = {r["transcript_id"]: r for r in result["records"]}
    reference = _pattern(1000)
    # t1 has a CDS (151-200), so only its coding exon is used on the + strand
    assert by_id["t1"]["sequence"] == reference[150:200]
    assert by_id["t1"]["gene_ids"] == ["g1"]
    # t2 has no CDS: its exon 501-700 is taken and reverse-complemented on the - strand
    assert by_id["t2"]["sequence"] == _revcomp(reference[500:700])
    assert by_id["t3"]["sequence"] == "G" * 30
    assert [r["transcript_id"] for r in result["records"]] == ["t1", "t2", "t3"]


def test_gffutils_line_to_seqfeature_uses_zero_based_location() -> None:
    line = "chr1\tens\tgene\t101\t400\t.\t+\t.\tID=g1;Name=GENE1"
    result = gu_api.gffutils_line_to_seqfeature(line, verbose=False)
    assert result["type"] == "gene"
    assert result["location"] == {"start": 100, "end": 400, "strand": 1}
    assert result["qualifiers"]["Name"] == ["GENE1"]
    with pytest.raises(ConfigurationError):
        gu_api.gffutils_line_to_seqfeature("chr1 gene 1 2", verbose=False)


def test_seqfeature_to_line_writes_one_based_gff() -> None:
    result = gu_api.seqfeature_to_gffutils_line(
        "chr1", "gene", 99, 400, feature_id="g9", attributes={"Name": "Z"}, verbose=False
    )
    assert result["line"] == "chr1\t.\tgene\t100\t400\t.\t+\t.\tID=g9;Name=Z"
    assert (result["start"], result["end"]) == (100, 400)
    with pytest.raises(ConfigurationError):
        gu_api.seqfeature_to_gffutils_line("chr1", "gene", 5, 5, verbose=False)
    with pytest.raises(ConfigurationError):
        gu_api.seqfeature_to_gffutils_line(
            "chr1", "gene", 0, 5, attributes={"seqid": "x"}, verbose=False
        )


def test_sanitize_writes_new_file_and_leaves_source_untouched(tmp_path: Path) -> None:
    source = _write(tmp_path / "small.gff3", SMALL_GFF3)
    before = source.read_bytes()
    output = tmp_path / "sanitised.gff3"
    result = gu_api.sanitize_gffutils_file(source, output, verbose=False)
    assert source.read_bytes() == before
    assert result["gene_count"] == 3
    lines = [line for line in output.read_text(encoding="utf-8").splitlines() if line[:1] != "#"]
    assert len(lines) == result["feature_count"] == 11
    assert all("gid=" in line for line in lines)
    with pytest.raises(ConfigurationError):
        gu_api.sanitize_gffutils_file(source, output, verbose=False)
    with pytest.raises(ConfigurationError):
        gu_api.sanitize_gffutils_file(source, source, verbose=False)


def test_merge_matches_upstream_groups(tmp_path: Path) -> None:
    gffutils = importlib.import_module("gffutils")
    merge_criteria = importlib.import_module("gffutils.merge_criteria")
    gff = _write(tmp_path / "small.gff3", SMALL_GFF3)
    db = gffutils.create_db(str(gff), ":memory:", verbose=False)
    upstream = list(
        db.merge(
            db.features_of_type("exon", order_by=("seqid", "start")),
            merge_criteria=(
                merge_criteria.seqid,
                merge_criteria.strand,
                merge_criteria.feature_type,
            ),
            multiline=False,
        )
    )
    result = gu_api.merge_gffutils_features(
        gff,
        featuretype="exon",
        criteria=["seqid", "strand", "feature_type"],
        verbose=False,
    )
    assert result["merged_count"] == len(upstream) == 3
    grouped = [r for r in result["merged"] if r["member_ids"]]
    assert [(r["start"], r["end"], r["member_ids"]) for r in grouped] == [(101, 400, ["e1", "e2"])]
    thresholded = gu_api.merge_gffutils_features(
        gff,
        featuretype="exon",
        criteria=["seqid", "overlap_any_threshold"],
        threshold=5,
        verbose=False,
    )
    assert thresholded["threshold"] == 5 and thresholded["merged_count"] >= 1
    with pytest.raises(ConfigurationError):
        gu_api.merge_gffutils_features(gff, featuretype="exon", criteria=["nope"], verbose=False)
    with pytest.raises(ConfigurationError):
        gu_api.merge_gffutils_features(
            gff, featuretype="exon", criteria=["seqid", "seqid"], verbose=False
        )


# gtfparse


def test_expand_attribute_strings_matches_upstream_and_hand_values() -> None:
    attribute_parsing = importlib.import_module("gtfparse.attribute_parsing")
    rows = [
        'gene_id "G1"; transcript_id "T1"; gene_name "A";',
        None,
        'gene_id "G2"; tag "x"; tag "y";',
    ]
    result = gt_api.expand_gtfparse_attribute_strings(rows, verbose=False)
    assert result["columns"]["gene_id"] == ["G1", "", "G2"]
    assert result["columns"]["tag"] == ["", "", "x,y"]
    assert result["row_count"] == 3
    upstream = attribute_parsing.expand_attribute_strings(rows, quote_char="'", usecols=None)
    assert result["columns"] == {key: list(values) for key, values in sorted(upstream.items())}
    only = gt_api.expand_gtfparse_attribute_strings(rows, usecols=["gene_id"], verbose=False)
    assert list(only["columns"]) == ["gene_id"]
    with pytest.raises(ConfigurationError):
        gt_api.expand_gtfparse_attribute_strings("gene_id G1", verbose=False)


def test_create_missing_transcripts_from_exons(tmp_path: Path) -> None:
    gtf = _write(tmp_path / "exons.gtf", EXONS_GTF)
    result = gt_api.create_missing_gtfparse_features(
        gtf,
        feature="transcript",
        key="transcript_id",
        extra_columns=["gene_id", "gene_name"],
        verbose=False,
    )
    assert result["created_count"] == 2
    created = {r["transcript_id"]: (r["start"], r["end"], r["strand"]) for r in result["records"]}
    assert created == {"T1": (101, 400, "+"), "T2": (501, 600, "-")}
    assert result["records"][0]["gene_name"] == "GENEA"


def test_create_missing_genes_and_existing_feature(tmp_path: Path) -> None:
    gtf = _write(tmp_path / "exons.gtf", EXONS_GTF)
    genes = gt_api.create_missing_gtfparse_features(
        gtf, feature="gene", key="gene_id", extra_columns=["gene_name"], verbose=False
    )
    assert [(r["gene_id"], r["start"], r["end"]) for r in genes["records"]] == [
        ("G1", 101, 400),
        ("G2", 501, 600),
    ]
    with_transcript = _write(
        tmp_path / "with_transcript.gtf",
        EXONS_GTF
        + 'chr1\tHAVANA\ttranscript\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n',
    )
    # upstream skips a feature type that already has rows, so no transcript is added
    existing = gt_api.create_missing_gtfparse_features(with_transcript, verbose=False)
    assert existing["created_count"] == 0 and existing["records"] == []
    with pytest.raises(ConfigurationError):
        gt_api.create_missing_gtfparse_features(gtf, extra_columns=["nope"], verbose=False)
