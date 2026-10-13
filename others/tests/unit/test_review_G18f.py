"""Unit tests for the G18f review wrappers: polars-bio and gtars, plus the mygene and myvariant
review modules that export no functions.

Each wrapped function is checked against a hand calculation or a direct call to the upstream
package. Inputs are written into pytest's ``tmp_path``; no network access is used.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dnakit.annotation import gtars_review_functions as gt_review
from dnakit.annotation import polars_bio_review_functions as pb_review
from dnakit.databases import mygene_review_functions, myvariant_review_functions
from dnakit.exceptions import ConfigurationError

SAM = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:73",
        "@RG\tID:rg1\tSM:s1",
        "r1\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTTA\tIIIIIIIIII\tRG:Z:rg1",
        "r2\t16\tchr1\t20\t40\t8M\t*\t0\t0\tGCCGATTA\tIIIIIIII\tRG:Z:rg1",
        "r3\t0\tchr1\t40\t30\t5M2I3M\t*\t0\t0\tACGTAGCTAG\tIIIIIIIIII\tNM:i:2",
        "",
    ]
)
REFERENCE = ">chr1\nACGTACGTTAGCCGATTAGCATGCAATGCCGTAGCTAGCTAGGATCGATCGTAGCTAGCTAGCTACGATCGAT\n"
VCF = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=1000>",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
        "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=12\tGT\t0/1\t1/1",
        "chr1\t25\t.\tC\tT\t40\tPASS\tDP=8\tGT\t0/0\t0/1",
        "",
    ]
)
BED = "chr1\t0\t10\tfeatA\t500\t+\nchr1\t20\t35\tfeatB\t300\t-\nchr2\t5\t15\tfeatC\t100\t+\n"
A2M = ">seq1\nACGT-ACGT\n>seq2\nAC-TACGTT\n"
A3M = ">seq1\nACGTACGT\n>seq2\nacgtACGT\n"
STO = "# STOCKHOLM 1.0\nseq1  ACGT-ACGT\nseq2  AC-TACGTT\n//\n"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture()
def alignments(tmp_path: Path) -> dict[str, Path]:
    pysam = pytest.importorskip("pysam")
    sam = _write(tmp_path / "reads.sam", SAM)
    reference = _write(tmp_path / "ref.fa", REFERENCE)
    pysam.faidx(str(reference))
    bam = tmp_path / "reads.bam"
    cram = tmp_path / "reads.cram"
    pysam.view("-b", "-o", str(bam), str(sam), catch_stdout=False)
    pysam.index(str(bam))
    pysam.view("-C", "-T", str(reference), "-o", str(cram), str(sam), catch_stdout=False)
    return {"sam": sam, "bam": bam, "cram": cram, "reference": reference}


@pytest.fixture()
def variants(tmp_path: Path) -> dict[str, Path]:
    pysam = pytest.importorskip("pysam")
    vcf = _write(tmp_path / "sample.vcf", VCF)
    bcf = tmp_path / "sample.bcf"
    source = pysam.VariantFile(str(vcf))
    target = pysam.VariantFile(str(bcf), "wb", header=source.header)
    for record in source:
        target.write(record)
    target.close()
    source.close()
    return {"vcf": vcf, "bcf": bcf}


# polars-bio: headers and schemas.


def test_read_file_header_reports_vcf_and_sam_header_fields(
    variants: dict[str, Path], alignments: dict[str, Path]
) -> None:
    vcf = pb_review.read_file_header(variants["vcf"], file_format="vcf", verbose=False)
    assert vcf["status"] == "completed"
    assert vcf["header"]["sample_names"] == ["s1", "s2"]
    assert vcf["header"]["contigs"] == [{"id": "chr1", "length": 1000}]
    assert "_datafusion_table_name" not in vcf["header"]
    sam = pb_review.read_file_header(alignments["sam"], file_format="sam", verbose=False)
    assert sam["header"]["sort_order"] == "coordinate"
    assert '"length":73' in sam["header"]["reference_sequences"]
    with pytest.raises(ConfigurationError):
        pb_review.read_file_header(variants["vcf"], file_format="fasta", verbose=False)


def test_describe_genomic_file_vcf_bcf_cool_and_cram(
    variants: dict[str, Path], alignments: dict[str, Path], tmp_path: Path
) -> None:
    vcf = pb_review.describe_genomic_file(variants["vcf"], file_format="vcf", verbose=False)
    assert [entry["name"] for entry in vcf["entries"]] == ["DP", "genotypes"]
    bcf = pb_review.describe_genomic_file(variants["bcf"], file_format="bcf", verbose=False)
    assert bcf["entry_count"] == vcf["entry_count"] == 2
    cram = pb_review.describe_genomic_file(
        alignments["cram"],
        file_format="cram",
        reference_path=alignments["reference"],
        verbose=False,
    )
    names = {entry["column_name"] for entry in cram["entries"]}
    assert {"name", "chrom", "start", "end", "NM", "RG"} <= names
    with pytest.raises(ConfigurationError):
        pb_review.describe_genomic_file(alignments["cram"], file_format="cram", verbose=False)
    with pytest.raises(ConfigurationError):
        pb_review.describe_genomic_file(
            variants["vcf"],
            file_format="vcf",
            reference_path=alignments["reference"],
            verbose=False,
        )
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    chromsizes = pd.Series({"chr1": 100, "chr2": 100}, name="length")
    bins = cooler.util.binnify(chromsizes, 20)
    pixels = pd.DataFrame({"bin1_id": [0, 0, 1, 5], "bin2_id": [0, 1, 2, 6], "count": [3, 1, 2, 4]})
    cool = tmp_path / "contacts.cool"
    cooler.create_cooler(str(cool), bins, pixels, dtypes={"count": "int"}, ordered=True)
    summary = pb_review.describe_genomic_file(cool, file_format="cool", verbose=False)
    assert summary["entries"][0]["resolution"] == 20
    assert summary["entries"][0]["nchroms"] == 2


# polars-bio: alignments, intervals, and sequences.


def test_read_cram_records_matches_hand_coordinates(alignments: dict[str, Path]) -> None:
    result = pb_review.read_cram_records(alignments["cram"], alignments["reference"], verbose=False)
    assert result["record_count"] == 3
    by_name = {row["name"]: row for row in result["alignments"]}
    # SAM POS 5 is 0-based start 4; the 10M alignment ends at 14.
    assert (by_name["r1"]["start"], by_name["r1"]["end"], by_name["r1"]["cigar"]) == (4, 14, "10M")
    # SAM POS 40 with 5M2I3M consumes 8 reference bases, so the end is 47.
    assert (by_name["r3"]["start"], by_name["r3"]["end"]) == (39, 47)


def test_read_bigbed_intervals_splits_stored_fields(tmp_path: Path) -> None:
    pybigtools = pytest.importorskip("pybigtools")
    path = tmp_path / "features.bb"
    writer = pybigtools.open(str(path), "w")
    writer.write(
        {"chr1": 100, "chr2": 100},
        [
            ("chr1", 0, 10, "featA\t500\t+"),
            ("chr1", 20, 35, "featB\t300\t-"),
            ("chr2", 5, 15, "featC\t100\t+"),
        ],
    )
    result = pb_review.read_bigbed_intervals(path, verbose=False)
    assert result["interval_count"] == 3
    first = result["intervals"][0]
    assert (first["chrom"], first["start"], first["end"]) == ("chr1", 0, 10)
    assert first["fields"] == ["featA", "500", "+"]


def test_read_alignment_sequences_a2m_and_a3m(tmp_path: Path) -> None:
    a2m = pb_review.read_alignment_sequences(_write(tmp_path / "a.a2m", A2M), verbose=False)
    assert a2m["sequence_count"] == 2
    assert a2m["alignment_columns"] == 9
    assert a2m["sequences"][0]["sequence"] == "ACGT-ACGT"
    a3m = pb_review.read_alignment_sequences(
        _write(tmp_path / "a.a3m", A3M), file_format="a3m", verbose=False
    )
    # The A3M lowercase insertion in seq2 leaves both rows at 8 characters.
    assert a3m["alignment_columns"] == 8
    with pytest.raises(ConfigurationError):
        pb_review.read_alignment_sequences(tmp_path / "a.a2m", file_format="sto", verbose=False)


def test_write_fastq_reads_round_trip_and_validation(tmp_path: Path) -> None:
    reads = [
        {"name": "r1", "description": "first", "sequence": "acgtn", "quality_scores": "IIII#"},
        {"name": "r2", "description": None, "sequence": "GG", "quality_scores": "II"},
    ]
    target = tmp_path / "out.fastq"
    written = pb_review.write_fastq_reads(reads, target, verbose=False)
    assert written["records_written"] == 2
    assert target.read_text(encoding="utf-8") == ("@r1 first\nACGTN\n+\nIIII#\n@r2\nGG\n+\nII\n")
    with pytest.raises(ConfigurationError):
        pb_review.write_fastq_reads(
            [{"name": "bad", "sequence": "ACGT", "quality_scores": "II"}],
            tmp_path / "bad.fastq",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        pb_review.write_fastq_reads(
            [{"name": "two words", "sequence": "AC", "quality_scores": "II"}],
            tmp_path / "space.fastq",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        pb_review.write_fastq_reads(reads, target, verbose=False)


def test_write_alignment_file_keeps_record_positions(alignments: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    out_dir = alignments["sam"].parent
    to_bam = pb_review.write_alignment_file(
        alignments["sam"],
        out_dir / "converted.bam",
        source_format="sam",
        output_format="bam",
        verbose=False,
    )
    assert to_bam["records_written"] == 3
    with pysam.AlignmentFile(str(out_dir / "converted.bam"), "rb", check_sq=False) as handle:
        positions = [record.reference_start + 1 for record in handle.fetch(until_eof=True)]
    assert positions == [5, 20, 40]
    pb_review.write_alignment_file(
        alignments["cram"],
        out_dir / "converted.sam",
        source_format="cram",
        output_format="sam",
        reference_path=alignments["reference"],
        verbose=False,
    )
    with pytest.raises(ConfigurationError):
        pb_review.write_alignment_file(
            alignments["sam"], out_dir / "needs_ref.cram", output_format="cram", verbose=False
        )


def test_query_genomic_tables_runs_sql_and_removes_tables(tmp_path: Path) -> None:
    pb = pytest.importorskip("polars_bio")
    bed = _write(tmp_path / "features.bed", BED)
    vcf = _write(tmp_path / "sample.vcf", VCF)
    result = pb_review.query_genomic_tables(
        {"peaks_t1": (bed, "bed"), "vars_t1": (vcf, "vcf")},
        'SELECT chrom, start, "end", name FROM peaks_t1 ORDER BY start',
        verbose=False,
    )
    # BED starts are 0-based; ORDER BY keeps the query order.
    assert [row["start"] for row in result["rows"]] == [0, 5, 20]
    assert result["columns"] == ["chrom", "start", "end", "name"]
    variants = pb_review.query_genomic_tables(
        {"vars_t2": (vcf, "vcf")}, "SELECT chrom, start, id FROM vars_t2", verbose=False
    )
    # VCF POS 10 is 0-based start 9.
    assert {row["start"] for row in variants["rows"]} == {9, 24}
    with pytest.raises(KeyError):
        pb.ctx.table("peaks_t1")
    with pytest.raises(ConfigurationError):
        pb_review.query_genomic_tables({"peaks_t3": (bed, "bed")}, "DROP TABLE peaks_t3")
    pb.register_bed(str(bed), name="taken_t4")
    try:
        with pytest.raises(ConfigurationError):
            pb_review.query_genomic_tables({"taken_t4": (bed, "bed")}, "SELECT * FROM taken_t4")
        assert pb.ctx.table("taken_t4") is not None
    finally:
        pb.ctx.deregister_table("taken_t4")


# gtars: LOLA universe helpers, region tokens, and VRS locations.


def _bed(path: Path, rows: list[tuple[str, int, int]]) -> Path:
    return _write(path, "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in rows))


def test_lola_restricted_universe_splits_overlaps(tmp_path: Path) -> None:
    set_a = _bed(tmp_path / "a.bed", [("chr1", 0, 10), ("chr1", 5, 20)])
    set_b = _bed(tmp_path / "b.bed", [("chr1", 0, 10), ("chr2", 1, 4)])
    result = gt_review.lola_restricted_universe([set_a, set_b], verbose=False)
    # chr1 [0, 20) is cut at 5 and 10 into three pieces; chr2 [1, 4) is kept whole.
    assert [(r["chrom"], r["start"], r["end"]) for r in result["regions"]] == [
        ("chr1", 0, 5),
        ("chr1", 5, 10),
        ("chr1", 10, 20),
        ("chr2", 1, 4),
    ]
    with pytest.raises(ConfigurationError):
        gt_review.lola_restricted_universe([], verbose=False)


def test_lola_check_universe_reports_coverage_and_many_to_many(tmp_path: Path) -> None:
    user = _bed(tmp_path / "user.bed", [("chr1", 0, 10), ("chr1", 5, 20)])
    universe = _bed(
        tmp_path / "uni.bed",
        [("chr1", 0, 10), ("chr1", 5, 20), ("chr2", 0, 50), ("chr1", 100, 110)],
    )
    result = gt_review.lola_check_universe([user], universe, verbose=False)
    (only,) = result["user_sets"]
    assert only["total_regions"] == 2
    assert only["regions_in_universe"] == 2
    assert only["coverage"] == 1.0
    # Each user region overlaps two universe regions, so both count as many-to-many.
    assert only["many_to_many"] == 2
    assert result["warnings"]


def test_lola_redefine_user_sets_maps_onto_universe(tmp_path: Path) -> None:
    set_a = _bed(tmp_path / "a.bed", [("chr1", 0, 10), ("chr1", 5, 20)])
    set_b = _bed(tmp_path / "b.bed", [("chr2", 1, 4)])
    universe = _bed(tmp_path / "uni.bed", [("chr1", 0, 10), ("chr1", 5, 20), ("chr2", 0, 50)])
    result = gt_review.lola_redefine_user_sets([set_a, set_b], universe, verbose=False)
    first, second = result["user_sets"]
    assert [(r["start"], r["end"]) for r in first["regions"]] == [(0, 10), (5, 20)]
    assert [(r["chrom"], r["start"], r["end"]) for r in second["regions"]] == [("chr2", 0, 50)]


def test_region_tokens_round_trip_and_validation(tmp_path: Path) -> None:
    target = tmp_path / "tokens.gtok"
    written = gt_review.write_region_tokens([3, 1, 2], target, verbose=False)
    assert written["token_count"] == 3
    read = gt_review.read_region_tokens(target, verbose=False)
    assert read["tokens"] == [3, 1, 2]
    with pytest.raises(ConfigurationError):
        gt_review.write_region_tokens([3, -1], tmp_path / "neg.gtok", verbose=False)
    with pytest.raises(ConfigurationError):
        gt_review.write_region_tokens([1], target, verbose=False)


def test_vrs_location_digest_matches_upstream(tmp_path: Path) -> None:
    refget = pytest.importorskip("gtars.refget")
    vrs = pytest.importorskip("gtars.vrs")
    digest = str(refget.sha512t24u_digest("ACGTACGTAC"))
    result = gt_review.vrs_location_digest(digest, 0, 4, verbose=False)
    assert result["digest"] == str(vrs.location_digest(digest, 0, 4))
    assert result["identifier"] == f"ga4gh:SL.{result['digest']}"
    bare = gt_review.vrs_location_digest(digest, 0, 4, prefixed=False, verbose=False)
    assert bare["identifier"] == result["digest"]
    with pytest.raises(ConfigurationError):
        gt_review.vrs_location_digest(digest, 4, 4, verbose=False)


# Review modules without functions.


def test_mygene_and_myvariant_review_modules_export_nothing() -> None:
    assert mygene_review_functions.__all__ == []
    assert myvariant_review_functions.__all__ == []
