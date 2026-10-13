"""Tests for the H4 family-variant check: new wrappers and the claims kept or rejected.

New wrappers (pyensembl filenames and source check, gtars alphabet records, bioframe row
selection) are checked against upstream calls or hand calculations. For the ledger rows kept as
``family-variant`` the tests call the existing dnakit wrapper with the claimed parameter and
compare with the upstream function. Rows marked ``duplicate-of`` are checked by object identity
or by reading the upstream source. Inputs are written into pytest's ``tmp_path``; there is no
network access.
"""

from __future__ import annotations

import hashlib
import importlib
import inspect
import sys
from pathlib import Path

import pytest

from dnakit.annotation import bioframe_family_functions as bf_family
from dnakit.annotation import gtars_family_functions as gt_family
from dnakit.annotation import gtars_functions as gt_api
from dnakit.annotation import gtars_review_functions as gt_review
from dnakit.annotation import polars_bio_functions as pb_api
from dnakit.annotation import polars_bio_review_functions as pb_review
from dnakit.databases import genomepy_functions as gp_api
from dnakit.databases import pyensembl_family_functions as pe_family
from dnakit.databases import pyensembl_functions as pe_api
from dnakit.databases import pyensembl_review_functions as pe_review
from dnakit.exceptions import ConfigurationError


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _resolve(qualname: str):
    parts = qualname.split(".")
    for cut in range(len(parts), 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:cut]))
        except ImportError:
            continue
        for part in parts[cut:]:
            obj = getattr(obj, part)
        return obj
    raise AssertionError(qualname)


# ----------------------------------------------------------------------------------------
# New wrapper: pyensembl filenames and canonical-source check.
# ----------------------------------------------------------------------------------------


def test_gtf_filename_matches_upstream_and_url_tail() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import ensembl_url_templates as templates

    for release, species in [(75, "homo_sapiens"), (109, "homo_sapiens"), (110, "mus_musculus")]:
        result = pe_family.ensembl_gtf_filename(release, species, verbose=False)
        assert result["filename"] == templates.make_gtf_filename(release, species)
        url = pe_api.ensembl_gtf_url(release, species, verbose=False)["url"]
        assert url.endswith("/" + result["filename"])
    # Hand value: release 75 is GRCh37 and has no patch variant.
    assert (
        pe_family.ensembl_gtf_filename(75, "homo_sapiens", verbose=False)["filename"]
        == "Homo_sapiens.GRCh37.75.gtf.gz"
    )


def test_sequence_fasta_filename_matches_upstream_for_each_type() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import ensembl_url_templates as templates

    expected = {
        "cdna": "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "ncrna": "Homo_sapiens.GRCh38.ncrna.fa.gz",
        "pep": "Homo_sapiens.GRCh38.pep.all.fa.gz",
    }
    for kind, filename in expected.items():
        result = pe_family.ensembl_sequence_fasta_filename(109, "homo_sapiens", kind, verbose=False)
        assert result["filename"] == filename
        assert filename == templates.make_fasta_filename(109, "homo_sapiens", kind)
        url = pe_api.ensembl_sequence_fasta_url(109, "homo_sapiens", kind, verbose=False)["url"]
        assert url.endswith("/" + filename)
    old = pe_family.ensembl_sequence_fasta_filename(75, "homo_sapiens", "cdna", verbose=False)
    assert old["filename"] == "Homo_sapiens.GRCh37.75.cdna.all.fa.gz"


def test_filename_wrappers_reject_bad_input() -> None:
    pytest.importorskip("pyensembl")
    with pytest.raises(ConfigurationError):
        pe_family.ensembl_gtf_filename(0, "homo_sapiens", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_family.ensembl_gtf_filename(109, "not_a_species", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_family.ensembl_sequence_fasta_filename(109, "homo_sapiens", "dna", verbose=False)
    with pytest.raises(ConfigurationError):
        pe_family.check_ensembl_canonical_source("", verbose=False)


def test_canonical_source_check_differs_from_the_dna_options_flag() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import genome_fasta_cache

    base = "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/"
    cases = {
        base + "Homo_sapiens.GRCh38.dna.toplevel.fa.gz": (True, True),
        # Official directory, but a cDNA file name: official path, not a DNA source.
        base + "Homo_sapiens.GRCh38.cdna.all.fa.gz": (True, False),
        "https://example.org/pub/release-109/fasta/homo_sapiens/dna/x.fa.gz": (False, False),
        "http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/x.dna.toplevel.fa.gz": (
            False,
            False,
        ),
    }
    for url, (official, dna_flag) in cases.items():
        checked = pe_family.check_ensembl_canonical_source(url, verbose=False)
        assert checked["official_source"] is official
        assert checked["official_source"] == genome_fasta_cache.is_canonical_source(url)
        options = pe_review.ensembl_genome_fasta_source_options(url, verbose=False)
        assert options["canonical"] is dna_flag


# ----------------------------------------------------------------------------------------
# New wrappers: gtars alphabet records.
# ----------------------------------------------------------------------------------------

FASTA = ">a first record\nACGTAC\nGGTT\n>b\nNNACGU\n"


def test_sequence_digest_record_matches_digests_and_gtars() -> None:
    gtars = pytest.importorskip("gtars")
    result = gt_family.sequence_digest_record("acgtacgtnn", name="s1", verbose=False)
    assert result["length"] == 10 and result["name"] == "s1" and result["description"] is None
    assert result["md5"] == hashlib.md5(b"ACGTACGTNN").hexdigest()
    sha = gt_api.sequence_digest("ACGTACGTNN", verbose=False)["digest"]
    assert result["sha512t24u"] == sha
    assert (
        gt_api.sequence_digest("ACGTACGTNN", algorithm="md5", verbose=False)["digest"]
        == (result["md5"])
    )
    direct = gtars.refget.digest_sequence(b"ACGTACGTNN", name="s1")
    assert result["alphabet"] == str(direct.metadata.alphabet) == "dna3bit"
    assert gt_family.sequence_digest_record("ACGT", verbose=False)["alphabet"] == "dna2bit"
    described = gt_family.sequence_digest_record("ACGT", description="x", verbose=False)
    assert described["description"] == "x" and described["name"] == ""


def test_sequence_digest_record_rejects_non_dna() -> None:
    pytest.importorskip("gtars")
    with pytest.raises(ConfigurationError):
        gt_family.sequence_digest_record("ACGT!", verbose=False)
    with pytest.raises(ConfigurationError):
        gt_family.sequence_digest_record("", verbose=False)


def test_load_fasta_collection_matches_digest_wrapper(tmp_path: Path) -> None:
    gtars = pytest.importorskip("gtars")
    fasta = _write(tmp_path / "t.fa", FASTA)
    loaded = gt_family.load_fasta_sequence_collection(fasta, verbose=False)
    digested = gt_api.digest_fasta_sequences(fasta, verbose=False)
    assert loaded["collection_digest"] == digested["collection_digest"]
    by_name = {row["name"]: row for row in loaded["sequences"]}
    assert by_name["a"]["sequence"] == "ACGTACGGTT" and by_name["a"]["length"] == 10
    assert by_name["a"]["description"] == "first record" and by_name["b"]["description"] is None
    assert by_name["b"]["sequence"] == "NNACGU"
    assert by_name["a"]["alphabet"] == "dna2bit" and by_name["b"]["alphabet"] == "dnaio"
    assert by_name["a"]["md5"] == hashlib.md5(b"ACGTACGGTT").hexdigest()
    for row in digested["sequences"]:
        assert by_name[row["name"]]["sha512t24u"] == row["sha512t24u"]
        assert by_name[row["name"]]["md5"] == row["md5"]
    direct = gtars.refget.load_fasta(str(fasta))
    assert [str(item.metadata.alphabet) for item in direct.sequences] == ["dna2bit", "dnaio"]
    assert loaded["total_bases"] == 16
    assert [row["name"] for row in loaded["sequences"]] == ["a", "b"]


def test_load_fasta_collection_options(tmp_path: Path) -> None:
    pytest.importorskip("gtars")
    fasta = _write(tmp_path / "t.fa", FASTA)
    lean = gt_family.load_fasta_sequence_collection(fasta, include_sequences=False, verbose=False)
    assert all("sequence" not in row for row in lean["sequences"])
    with pytest.raises(ConfigurationError):
        gt_family.load_fasta_sequence_collection(fasta, max_bases=5, verbose=False)
    with pytest.raises(ConfigurationError):
        gt_family.load_fasta_sequence_collection(fasta, max_rows=1, verbose=False)
    with pytest.raises(ConfigurationError):
        gt_family.load_fasta_sequence_collection(fasta, include_sequences="yes", verbose=False)


# ----------------------------------------------------------------------------------------
# New wrapper: bioframe row selection.
# ----------------------------------------------------------------------------------------

BED = "chr1\t10\t30\nchr1\t50\t80\nchr2\t5\t15\nchr1\t100\t120\n"


def test_select_overlap_rows_mask_indices_and_labels(tmp_path: Path) -> None:
    bioframe = pytest.importorskip("bioframe")
    bed = _write(tmp_path / "a.bed", BED)
    result = bf_family.select_bioframe_overlap_rows(bed, "chr1:25-55", verbose=False)
    # Rows 0 (10-30) and 1 (50-80) overlap 25-55; row 2 is on chr2; row 3 starts at 100.
    assert result["mask"] == [True, True, False, False]
    assert result["indices"] == [0, 1] and result["labels"] == [0, 1]
    assert result["row_count"] == 4 and result["selected_count"] == 2
    frame = bioframe.read_table(str(bed), schema="bed3", comment="#")
    assert result["mask"] == list(bioframe.select_mask(frame, "chr1:25-55"))
    assert result["indices"] == list(bioframe.select_indices(frame, "chr1:25-55"))
    assert result["labels"] == list(bioframe.select_labels(frame, "chr1:25-55"))
    whole = bf_family.select_bioframe_overlap_rows(bed, "chr1", verbose=False)
    assert whole["indices"] == [0, 1, 3]
    # Half-open: an interval ending exactly at the region start is not selected.
    edge = bf_family.select_bioframe_overlap_rows(bed, "chr1:30-50", verbose=False)
    assert edge["indices"] == []


def test_select_overlap_rows_agrees_with_region_wrapper(tmp_path: Path) -> None:
    pytest.importorskip("bioframe")
    from dnakit.annotation import bioframe_review_functions as bf_review

    bed = _write(tmp_path / "a.bed", BED)
    rows = bf_family.select_bioframe_overlap_rows(bed, "chr1:25-55", verbose=False)
    records = bf_review.select_bioframe_region(bed, "chr1:25-55", verbose=False)
    assert records["interval_count"] == rows["selected_count"]


def test_select_overlap_rows_rejects_bad_input(tmp_path: Path) -> None:
    pytest.importorskip("bioframe")
    bed = _write(tmp_path / "a.bed", BED)
    with pytest.raises(ConfigurationError):
        bf_family.select_bioframe_overlap_rows(bed, "", verbose=False)
    with pytest.raises(ConfigurationError):
        bf_family.select_bioframe_overlap_rows(bed, "chr1:5-1", verbose=False)
    with pytest.raises(ConfigurationError):
        bf_family.select_bioframe_overlap_rows(bed, "chr1", max_rows=2, verbose=False)


def test_bioframe_ops_aliases_are_the_same_objects() -> None:
    bioframe = pytest.importorskip("bioframe")
    from bioframe import ops

    for name in ("select_mask", "select_indices", "select_labels"):
        assert getattr(ops, name) is getattr(bioframe, name)


# ----------------------------------------------------------------------------------------
# family-variant rows kept: polars-bio.
# ----------------------------------------------------------------------------------------

SAM = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:73",
        "r1\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTTA\tIIIIIIIIII",
        "r2\t16\tchr1\t20\t40\t8M\t*\t0\t0\tGCCGATTA\tIIIIIIII",
        "",
    ]
)
REFERENCE = ">chr1\nACGTACGTTAGCCGATTAGCATGCAATGCCGTAGCTAGCTAGGATCGATCGTAGCTAGCTAGCTACGATCGAT\n"
VCF = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=1000>",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=12",
        "chr1\t25\t.\tC\tT\t40\tPASS\tDP=8",
        "",
    ]
)
GTF = (
    'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_id "G1";\n'
    'chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
)
GFF = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1;Name=A\n"
    "chr1\tsrc\texon\t101\t150\t.\t+\t.\tID=e1;Parent=g1\n"
)
BEDFILE = "chr1\t0\t10\tfeatA\t500\t+\nchr1\t20\t35\tfeatB\t300\t-\n"
FASTQ = "@r1\nACGT\n+\nIIII\n@r2\nTTGA\n+\nIIII\n"
A3M = ">seq1\nACGTACGT\n>seq2\nacgtACGT\n"


@pytest.fixture()
def pysam_files(tmp_path: Path) -> dict[str, Path]:
    pysam = pytest.importorskip("pysam")
    sam = _write(tmp_path / "reads.sam", SAM)
    reference = _write(tmp_path / "ref.fa", REFERENCE)
    pysam.faidx(str(reference))
    bam = tmp_path / "reads.bam"
    cram = tmp_path / "reads.cram"
    pysam.view("-b", "-o", str(bam), str(sam), catch_stdout=False)
    pysam.index(str(bam))
    pysam.view("-C", "-T", str(reference), "-o", str(cram), str(sam), catch_stdout=False)
    vcf = _write(tmp_path / "s.vcf", VCF)
    bcf = tmp_path / "s.bcf"
    source = pysam.VariantFile(str(vcf))
    target = pysam.VariantFile(str(bcf), "wb", header=source.header)
    for record in source:
        target.write(record)
    target.close()
    source.close()
    return {"sam": sam, "bam": bam, "cram": cram, "reference": reference, "vcf": vcf, "bcf": bcf}


def test_sql_registrars_cover_every_register_row() -> None:
    pb = pytest.importorskip("polars_bio")
    claimed = "gff gtf vcf bcf fastq fasta a2m a3m sto bam sam bed bigwig cool bigbed pairs".split()  # noqa: SIM905
    mapping = pb_review._SQL_REGISTRARS
    for fmt in claimed:
        assert mapping[fmt] == f"register_{fmt}"
        assert callable(getattr(pb, f"register_{fmt}"))


@pytest.mark.parametrize("fmt", ["bed", "gtf", "gff", "vcf", "fasta", "fastq", "sam", "bam"])
def test_query_genomic_tables_reaches_each_register_function(
    fmt: str, tmp_path: Path, pysam_files: dict[str, Path]
) -> None:
    pytest.importorskip("polars_bio")
    files = {
        "bed": _write(tmp_path / "x.bed", BEDFILE),
        "gtf": _write(tmp_path / "x.gtf", GTF),
        "gff": _write(tmp_path / "x.gff", GFF),
        "vcf": pysam_files["vcf"],
        "fasta": _write(tmp_path / "x.fasta", ">s1\nACGT\n>s2\nTTGA\n"),
        "fastq": _write(tmp_path / "x.fastq", FASTQ),
        "sam": pysam_files["sam"],
        "bam": pysam_files["bam"],
    }
    expected = {"bed": 2, "gtf": 2, "gff": 2, "vcf": 2, "fasta": 2, "fastq": 2, "sam": 2, "bam": 2}
    result = pb_review.query_genomic_tables(
        {"t": (files[fmt], fmt)}, "SELECT COUNT(*) AS n FROM t", verbose=False
    )
    assert result["rows"] == [{"n": expected[fmt]}]


@pytest.mark.parametrize("fmt", ["sam", "bam", "cram"])
def test_describe_genomic_file_reaches_each_describe_function(
    fmt: str, pysam_files: dict[str, Path]
) -> None:
    pb = pytest.importorskip("polars_bio")
    path = pysam_files[fmt]
    if fmt == "cram":
        result = pb_review.describe_genomic_file(
            path, file_format=fmt, reference_path=pysam_files["reference"], verbose=False
        )
        direct = pb.describe_cram(str(path), reference_path=str(pysam_files["reference"]))
    else:
        result = pb_review.describe_genomic_file(path, file_format=fmt, verbose=False)
        direct = getattr(pb, f"describe_{fmt}")(str(path))
    assert result["entry_count"] == direct.height > 0


def test_read_wrappers_reach_bam_gff_bcf_a3m(tmp_path: Path, pysam_files: dict[str, Path]) -> None:
    pb = pytest.importorskip("polars_bio")
    bam = pb_api.read_alignment_records(pysam_files["bam"], file_format="bam", verbose=False)
    assert bam["file_format"] == "bam" and bam["record_count"] == 2
    direct = pb.read_bam(str(pysam_files["bam"]), use_zero_based=True)
    assert sorted(r["start"] for r in bam["alignments"]) == sorted(direct["start"].to_list())

    gff = _write(tmp_path / "x.gff", GFF)
    annotation = pb_api.read_gene_annotation(gff, file_format="gff", verbose=False)
    assert annotation["file_format"] == "gff" and annotation["feature_count"] == 2
    assert sorted(f["type"] for f in annotation["features"]) == ["exon", "gene"]

    bcf = pb_api.read_variant_calls(pysam_files["bcf"], file_format="bcf", verbose=False)
    vcf = pb_api.read_variant_calls(pysam_files["vcf"], verbose=False)
    assert bcf["variant_count"] == 2
    assert [v["start"] for v in bcf["variants"]] == [v["start"] for v in vcf["variants"]]

    a3m = _write(tmp_path / "x.a3m", A3M)
    msa = pb_review.read_alignment_sequences(a3m, file_format="a3m", verbose=False)
    assert msa["file_format"] == "a3m" and msa["sequence_count"] == 2
    assert pb.read_a3m(str(a3m)).height == 2


@pytest.mark.parametrize("output_format", ["bam", "cram"])
def test_write_alignment_file_reaches_bam_and_cram_writers(
    output_format: str, tmp_path: Path, pysam_files: dict[str, Path]
) -> None:
    pytest.importorskip("polars_bio")
    pysam = pytest.importorskip("pysam")
    target = tmp_path / f"out.{output_format}"
    result = pb_review.write_alignment_file(
        pysam_files["sam"],
        target,
        output_format=output_format,
        reference_path=pysam_files["reference"],
        verbose=False,
    )
    assert result["output_format"] == output_format and result["records_written"] == 2
    with pysam.AlignmentFile(
        str(target),
        "rc" if output_format == "cram" else "rb",
        reference_filename=str(pysam_files["reference"]),
    ) as handle:
        assert [read.query_name for read in handle] == ["r1", "r2"]


def test_scan_and_sink_rows_are_the_lazy_forms_of_the_wrapped_calls() -> None:
    pb = pytest.importorskip("polars_bio")
    from polars_bio.io import IOOperations

    formats = "bam sam cram bed cool bigbed fasta a2m a3m sto fastq gff gtf pairs bcf vcf".split()  # noqa: SIM905
    for fmt in formats:
        reader = IOOperations.__dict__[f"read_{fmt}"]
        source = inspect.getsource(getattr(reader, "__func__", reader))
        assert f"scan_{fmt}(" in source and ".collect()" in source, fmt
        assert callable(getattr(pb, f"scan_{fmt}"))
    for fmt in ("vcf", "fasta", "fastq", "bam", "sam", "cram"):
        writer = IOOperations.__dict__[f"write_{fmt}"]
        sink = IOOperations.__dict__[f"sink_{fmt}"]
        for function in (writer, sink):
            code = inspect.getsource(getattr(function, "__func__", function))
            assert "_write_file(" in code or "_write_bam_file(" in code, (fmt, function)


# ----------------------------------------------------------------------------------------
# family-variant rows kept: gffutils merge criteria.
# ----------------------------------------------------------------------------------------

MERGE_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\tsrc\texon\t1\t100\t.\t+\t.\tID=e1",
        "chr1\tsrc\texon\t50\t150\t.\t+\t.\tID=e2",
        "chr1\tsrc\texon\t151\t200\t.\t+\t.\tID=e3",
        "chr1\tsrc\texon\t300\t400\t.\t+\t.\tID=e4",
        "chr1\tsrc\texon\t320\t380\t.\t-\t.\tID=e5",
        "chr2\tsrc\texon\t10\t20\t.\t+\t.\tID=e6",
        "chr2\tsrc\texon\t20\t30\t.\t+\t.\tID=e7",
        "",
    ]
)
MERGE_CRITERIA = {
    "seqid": [("chr1", 1, 400), ("chr2", 10, 30)],
    "strand": [("chr1", 1, 400), ("chr1", 320, 380), ("chr2", 10, 30)],
    "feature_type": [("chr1,chr2", 1, 400)],
    "exact_coordinates_only": [
        ("chr1", 1, 100),
        ("chr1", 50, 150),
        ("chr1", 151, 200),
        ("chr1", 300, 400),
        ("chr1", 320, 380),
        ("chr2", 10, 20),
        ("chr2", 20, 30),
    ],
    "overlap_end_inclusive": [("chr1", 1, 200), ("chr1", 300, 400), ("chr2", 10, 30)],
    "overlap_start_inclusive": [
        ("chr1", 1, 100),
        ("chr1", 50, 150),
        ("chr1", 151, 200),
        ("chr1", 300, 400),
        ("chr2", 10, 20),
        ("chr2", 20, 30),
    ],
    "overlap_any_inclusive": [("chr1", 1, 200), ("chr1", 300, 400), ("chr2", 10, 30)],
}


@pytest.mark.parametrize("name", sorted(MERGE_CRITERIA))
def test_merge_criteria_selected_by_name_match_gffutils(name: str, tmp_path: Path) -> None:
    gffutils = pytest.importorskip("gffutils")
    from gffutils import merge_criteria

    from dnakit.annotation import gffutils_review_functions as gu_review

    path = _write(tmp_path / "m.gff", MERGE_GFF)
    result = gu_review.merge_gffutils_features(
        path, featuretype="exon", criteria=[name], verbose=False
    )
    assert [(m["seqid"], m["start"], m["end"]) for m in result["merged"]] == MERGE_CRITERIA[name]
    db = gffutils.create_db(str(path), ":memory:", verbose=False)
    direct = list(
        db.merge(
            db.features_of_type("exon", order_by=("seqid", "start")),
            merge_criteria=[getattr(merge_criteria, name)],
            multiline=False,
        )
    )
    assert sorted((f.seqid, f.start, f.end) for f in direct) == MERGE_CRITERIA[name]


@pytest.mark.parametrize(
    "name", ["overlap_end_threshold", "overlap_start_threshold", "overlap_any_threshold"]
)
def test_threshold_merge_criteria_use_the_threshold_argument(name: str, tmp_path: Path) -> None:
    gffutils = pytest.importorskip("gffutils")
    from gffutils import merge_criteria

    from dnakit.annotation import gffutils_review_functions as gu_review

    path = _write(tmp_path / "m.gff", MERGE_GFF)
    db = gffutils.create_db(str(path), ":memory:", verbose=False)
    for threshold in (0, 10, 60):
        result = gu_review.merge_gffutils_features(
            path, featuretype="exon", criteria=[name], threshold=threshold, verbose=False
        )
        direct = list(
            db.merge(
                db.features_of_type("exon", order_by=("seqid", "start")),
                merge_criteria=[getattr(merge_criteria, name)(threshold)],
                multiline=False,
            )
        )
        assert sorted((m["seqid"], m["start"], m["end"]) for m in result["merged"]) == sorted(
            (f.seqid, f.start, f.end) for f in direct
        ), (name, threshold)

    # Hand value: the threshold is the gap tolerated between features. e1 and e2 merge to 1-150,
    # which ends one base before e3 (151-200); a threshold of 0 keeps them apart and 10 joins them.
    def spans(threshold: int) -> set[tuple[str, int, int]]:
        merged = gu_review.merge_gffutils_features(
            path, featuretype="exon", criteria=[name], threshold=threshold, verbose=False
        )["merged"]
        return {(m["seqid"], m["start"], m["end"]) for m in merged}

    if name != "overlap_start_threshold":
        assert ("chr1", 151, 200) in spans(0) and ("chr1", 1, 200) not in spans(0)
        assert ("chr1", 1, 200) in spans(10)


# ----------------------------------------------------------------------------------------
# family-variant rows kept: gtars, pybedtools, pyensembl, genomepy.
# ----------------------------------------------------------------------------------------

GENES = "\n".join(
    [
        'chr1\tens\tgene\t1001\t2000\t.\t+\t.\tgene_id "G1"; gene_name "ABC"; '
        'gene_type "protein_coding";',
        'chr1\tens\texon\t1001\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_type "protein_coding";',
        'chr1\tens\texon\t1801\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_type "protein_coding";',
        "",
    ]
)
REGIONS = "chr1\t950\t1010\tp1\nchr1\t1500\t1600\tp2\nchr1\t9000\t9100\tp4\n"


def test_partition_counts_chrom_sizes_reaches_expected_partitions(tmp_path: Path) -> None:
    gtars = pytest.importorskip("gtars")
    regions = _write(tmp_path / "r.bed", REGIONS)
    genes = _write(tmp_path / "g.gtf", GENES)
    sizes = {"chr1": 10000}
    wrapped = gt_api.partition_region_counts(regions, genes, chrom_sizes=sizes, verbose=False)
    assert wrapped["mode"] == "expected"
    models = gtars.models
    partitions = models.PartitionList.from_gtf(
        str(genes),
        core_prom=100,
        prox_prom=2000,
        filter_protein_coding=True,
        convert_ensembl_ucsc=True,
        chrom_sizes=sizes,
    )
    direct = gtars.genomic_distributions.calc_expected_partitions(
        models.RegionSet(str(regions)), partitions, sizes
    )
    assert [row["partition"] for row in wrapped["partitions"]] == [
        str(x) for x in direct["partition"]
    ]
    assert [row["expected"] for row in wrapped["partitions"]] == pytest.approx(
        [float(x) for x in direct["expected"]]
    )
    plain = gt_api.partition_region_counts(regions, genes, verbose=False)
    assert plain["mode"] == "observed" and "expected" not in plain["partitions"][0]


def test_sequence_digest_md5_algorithm_reaches_md5_digest() -> None:
    gtars = pytest.importorskip("gtars")
    wrapped = gt_api.sequence_digest("acgtnn", algorithm="md5", verbose=False)
    assert (
        wrapped["digest"] == gtars.refget.md5_digest("ACGTNN") == hashlib.md5(b"ACGTNN").hexdigest()
    )


def test_vrs_allele_id_unprefixed_reaches_vrs_digest() -> None:
    gtars = pytest.importorskip("gtars")
    seq_digest = gt_api.sequence_digest("ACGTACGTAC", verbose=False)["digest"]
    bare = gt_api.vrs_allele_id(seq_digest, 3, 4, "G", prefixed=False, verbose=False)
    full = gt_api.vrs_allele_id(seq_digest, 3, 4, "G", verbose=False)
    assert bare["identifier"] == gtars.vrs.vrs_digest(seq_digest, 3, 4, "G")
    assert full["identifier"] == "ga4gh:VA." + bare["identifier"]


def test_gtok_string_tokens_are_the_same_ids_as_text(tmp_path: Path) -> None:
    utils = pytest.importorskip("gtars.utils")
    path = tmp_path / "t.gtok"
    gt_review.write_region_tokens([5, 1, 9], path, verbose=False)
    ids = gt_review.read_region_tokens(path, verbose=False)["tokens"]
    assert ids == [5, 1, 9]
    assert utils.read_tokens_from_gtok_as_strings(str(path)) == [str(token) for token in ids]


def test_pybedtools_length_filters_are_selected_by_the_arguments(tmp_path: Path) -> None:
    pybedtools = pytest.importorskip("pybedtools")
    from pybedtools import featurefuncs

    from dnakit.annotation import pybedtools_review_functions as pbt_review

    bed = _write(tmp_path / "x.bed", "chr1\t0\t10\nchr1\t20\t50\nchr1\t60\t65\n")
    longer = pbt_review.filter_pybedtools_intervals_by_length(bed, greater_than=8, verbose=False)
    shorter = pbt_review.filter_pybedtools_intervals_by_length(bed, less_than=8, verbose=False)
    sizes = lambda r: sorted(i["end"] - i["start"] for i in r["records"])  # noqa: E731
    assert sizes(longer) == [10, 30] and sizes(shorter) == [5]
    direct = pybedtools.BedTool(str(bed)).filter(featurefuncs.greater_than, size=8)
    assert sorted(len(i) for i in direct) == [10, 30]
    direct = pybedtools.BedTool(str(bed)).filter(featurefuncs.less_than, size=8)
    assert sorted(len(i) for i in direct) == [5]


def test_pybedtools_random_intersection_is_broken_upstream(tmp_path: Path) -> None:
    pybedtools = pytest.importorskip("pybedtools")
    from pybedtools import stats

    first = _write(tmp_path / "a.bed", "chr1\t10\t30\nchr1\t50\t80\n")
    second = _write(tmp_path / "b.bed", "chr1\t20\t60\n")
    genome = _write(tmp_path / "g.txt", "chr1\t1000\n")
    # close_or_delete throws StopIteration into an already closed generator, so every call
    # raises. If a pybedtools release fixes this, wrap both functions and revisit the ledger.
    for function in (stats.random_intersection, stats.random_intersection_bp):
        with pytest.raises(StopIteration):
            function(
                pybedtools.BedTool(str(first)),
                pybedtools.BedTool(str(second)),
                str(genome),
                {"seed": 1},
                {},
            )


def test_gtfparse_pandas_variant_delegates_to_parse_gtf(tmp_path: Path) -> None:
    pytest.importorskip("gtfparse")
    import gtfparse

    read_module = sys.modules["gtfparse.read_gtf"]
    assert gtfparse.parse_gtf_pandas is read_module.parse_gtf_pandas
    source = inspect.getsource(gtfparse.parse_gtf_pandas)
    assert "df = parse_gtf(*args, **kwargs)" in source
    path = _write(tmp_path / "x.gtf", GENES)
    plain = gtfparse.parse_gtf(str(path), split_attributes=False)
    again = gtfparse.parse_gtf_pandas(str(path), split_attributes=False)
    assert type(plain) is type(again) and plain.equals(again)
    from dnakit.annotation import gtfparse_functions as gtf_api

    wrapped = gtf_api.parse_gtfparse_fixed_columns(path, verbose=False)
    assert wrapped["record_count"] == len(plain)


def test_genomicranges_parse_gtf_is_the_first_step_of_read_gtf(tmp_path: Path) -> None:
    genomicranges = pytest.importorskip("genomicranges")
    from genomicranges.io import gtf

    from dnakit.annotation import genomicranges_functions as gr_api

    source = inspect.getsource(gtf.read_gtf)
    assert "parse_gtf(file" in source and "GenomicRanges.from_pandas(data)" in source
    path = _write(tmp_path / "x.gtf", GENES)
    frame = gtf.parse_gtf(str(path), compressed=False)
    wrapped = gr_api.read_gtf_ranges(path, verbose=False)
    # parse_gtf ends are the GTF ends minus one; the wrapper restores the GTF end.
    assert sorted(frame["ends"] + 1) == sorted(f["end"] for f in wrapped["features"])
    assert sorted(frame["starts"] - 1) == sorted(f["start"] for f in wrapped["features"])
    assert genomicranges.read_gtf(str(path)) is not None


def test_species_record_by_reference_reaches_find_species_by_reference() -> None:
    pytest.importorskip("pyensembl")
    from pyensembl import reference_name

    wrapped = pe_review.ensembl_species_record("grch37", by="reference", verbose=False)
    direct = reference_name.find_species_by_reference("grch37")
    assert wrapped["by"] == "reference" and wrapped["latin_name"] == direct.latin_name
    assert "by" in inspect.signature(pe_review.ensembl_species_record).parameters


def test_bed_sequences_reach_the_region_and_bed_helpers(tmp_path: Path) -> None:
    genomepy = pytest.importorskip("genomepy")
    from genomepy.genome import sequences as helpers

    from dnakit.databases import write_bed_sequences

    fasta = _write(
        tmp_path / "g.fa", ">chr1\nACGTACGTACGTAAACCCGGGTTT\n>chr2\nTTTTGGGGCCCCAAAATTTT\n"
    )
    genome = genomepy.Genome(str(fasta))
    assert "bed" in inspect.signature(write_bed_sequences).parameters

    # region list file (regions_to_seqs and region_to_seq)
    regions = _write(tmp_path / "r.regions", "chr1:2-10\nchr2:0-5\n")
    out = tmp_path / "regions.fa"
    write_bed_sequences(fasta, regions, out, extend_down=2, verbose=False)
    direct = {s.name: s.seq for s in helpers.regions_to_seqs(genome, str(regions), 0, 2)}
    lines = out.read_text().split()
    assert dict(zip(lines[0::2], lines[1::2], strict=True)) == {
        f">{k}": v for k, v in direct.items()
    }
    assert helpers.region_to_seq(genome, "chr1:2-10", 0, 2) == direct["chr1:2-10"]

    # BED6 file, minus strand (bed_to_seqs and bed_to_seq)
    bed = _write(tmp_path / "r.bed", "chr1\t0\t6\tfeatA\t0\t-\nchr2\t0\t5\tfeatB\t0\t+\n")
    out = tmp_path / "bed.fa"
    write_bed_sequences(fasta, bed, out, stranded=True, verbose=False)
    text = out.read_text().split("\n")
    # chr1:0-6 is ACGTAC; its reverse complement is GTACGT.
    assert text[0] == ">chr1:0-6 featA" and text[1] == "GTACGT"
    one = helpers.bed_to_seq(genome, ["chr1", "0", "6", "featA", "0", "-"], stranded=True)
    assert one.seq == "GTACGT" and one.name == "chr1:0-6 featA"
    assert [s.seq for s in helpers.bed_to_seqs(genome, str(bed), stranded=True)] == [
        "GTACGT",
        "TTTTG",
    ]


# ----------------------------------------------------------------------------------------
# duplicate-of rows: re-exports are the same objects.
# ----------------------------------------------------------------------------------------

PYENSEMBL_REEXPORTS = {
    "pyensembl.normalize_reference_name": "pyensembl.reference_name.normalize_reference_name",
    "pyensembl.find_species_by_reference": "pyensembl.reference_name.find_species_by_reference",
    "pyensembl.which_reference": "pyensembl.reference_name.which_reference",
    "pyensembl.find_nearest_locus": "pyensembl.search.find_nearest_locus",
    "pyensembl.find_species_by_name": "pyensembl.species.find_species_by_name",
    "pyensembl.database.normalize_chromosome": "pyensembl.normalization.normalize_chromosome",
    "pyensembl.database.normalize_strand": "pyensembl.normalization.normalize_strand",
    "pyensembl.database.match_version": "pyensembl.versioned_ids.match_version",
    "pyensembl.dated_releases.find_species_by_name": "pyensembl.species.find_species_by_name",
    "pyensembl.ensembl_annotation.make_dated_release_urls": (
        "pyensembl.ensembl_url_templates.make_dated_release_urls"
    ),
    "pyensembl.ensembl_annotation.find_species_by_name": "pyensembl.species.find_species_by_name",
    "pyensembl.ensembl_release.canonical_source_options": (
        "pyensembl.genome_fasta_cache.canonical_source_options"
    ),
    "pyensembl.ensembl_release.dated_genome_fasta_mask": (
        "pyensembl.ensembl_url_templates.dated_genome_fasta_mask"
    ),
    "pyensembl.ensembl_release.make_dated_release_urls": (
        "pyensembl.ensembl_url_templates.make_dated_release_urls"
    ),
    "pyensembl.ensembl_release.make_gtf_url": "pyensembl.ensembl_url_templates.make_gtf_url",
    "pyensembl.ensembl_release.make_fasta_url": "pyensembl.ensembl_url_templates.make_fasta_url",
    "pyensembl.ensembl_release.make_genome_fasta_url": (
        "pyensembl.ensembl_url_templates.make_genome_fasta_url"
    ),
    "pyensembl.ensembl_url_templates.find_species_by_name": (
        "pyensembl.species.find_species_by_name"
    ),
    "pyensembl.genome.merge_intervals": "pyensembl.common.merge_intervals",
    "pyensembl.genome.normalize_aliases": "pyensembl.gene_name_aliases.normalize_aliases",
    "pyensembl.genome.stable_ensembl_gene_id": "pyensembl.gene_name_aliases.stable_ensembl_gene_id",
    "pyensembl.genome.normalize_chromosome": "pyensembl.normalization.normalize_chromosome",
    "pyensembl.genome.normalize_strand": "pyensembl.normalization.normalize_strand",
    "pyensembl.genome.find_nearest_locus": "pyensembl.search.find_nearest_locus",
    "pyensembl.genome.match_version": "pyensembl.versioned_ids.match_version",
    "pyensembl.genome_fasta.normalize_chromosome": "pyensembl.normalization.normalize_chromosome",
    "pyensembl.locus.normalize_chromosome": "pyensembl.normalization.normalize_chromosome",
    "pyensembl.locus.normalize_strand": "pyensembl.normalization.normalize_strand",
    "pyensembl.reference_name.find_species_by_name": "pyensembl.species.find_species_by_name",
    "pyensembl.sequence_data.first_fasta_id": "pyensembl.fasta.first_fasta_id",
    "pyensembl.sequence_data.parse_fasta_dictionary": "pyensembl.fasta.parse_fasta_dictionary",
    "pyensembl.sequence_data.match_version": "pyensembl.versioned_ids.match_version",
    "pyensembl.shell.canonical_source_options": (
        "pyensembl.genome_fasta_cache.canonical_source_options"
    ),
    "pyensembl.shell.find_species_by_reference": (
        "pyensembl.reference_name.find_species_by_reference"
    ),
    "pyensembl.shell.normalize_reference_name": (
        "pyensembl.reference_name.normalize_reference_name"
    ),
    "pyensembl.shell.find_species_by_name": "pyensembl.species.find_species_by_name",
    "pyensembl.transcript.merge_intervals": "pyensembl.common.merge_intervals",
    "pyensembl.transcript.lookup_sequence_with_version_fallback": (
        "pyensembl.sequence_data.lookup_sequence_with_version_fallback"
    ),
}


def test_pyensembl_reexports_are_the_wrapped_objects() -> None:
    pytest.importorskip("pyensembl")
    for alias, target in PYENSEMBL_REEXPORTS.items():
        assert _resolve(alias) is _resolve(target), alias


def test_pyensembl_wrappers_call_the_re_exported_functions() -> None:
    pytest.importorskip("pyensembl")
    calls = {
        pe_api.normalize_genome_reference_name: "normalize_reference_name",
        pe_api.ensembl_reference_name: "which_reference",
        pe_api.find_nearest_genomic_locus: "find_nearest_locus",
        pe_api.normalize_chromosome_name: "normalize_chromosome",
        pe_api.normalize_strand_symbol: "normalize_strand",
        pe_api.match_identifier_version: "match_version",
        pe_api.ensembl_gtf_url: "make_gtf_url",
        pe_api.ensembl_sequence_fasta_url: "make_fasta_url",
        pe_api.ensembl_genome_fasta_url: "make_genome_fasta_url",
        pe_api.merge_genomic_intervals: "merge_intervals",
        pe_api.normalize_gene_aliases: "normalize_aliases",
        pe_api.stable_ensembl_id: "stable_ensembl_gene_id",
        pe_api.read_first_fasta_id: "first_fasta_id",
        pe_api.read_fasta_records: "parse_fasta_dictionary",
        pe_api.lookup_sequence_version_fallback: "lookup_sequence_with_version_fallback",
        pe_review.ensembl_dated_release_urls: "make_dated_release_urls",
        pe_review.ensembl_genome_fasta_source_options: "canonical_source_options",
        pe_review.ensembl_dated_genome_fasta_mask: "dated_genome_fasta_mask",
    }
    for function, upstream in calls.items():
        assert upstream in inspect.getsource(function), function.__name__
    species = inspect.getsource(pe_review.ensembl_species_record)
    assert "find_species_by_name" in species and "find_species_by_reference" in species


def test_genomepy_reexports_are_the_wrapped_objects() -> None:
    pytest.importorskip("genomepy")
    pairs = {
        "genomepy.annotation.read_annot": "genomepy.annotation.utils.read_annot",
        "genomepy.providers.base.read_annot": "genomepy.annotation.utils.read_annot",
        "genomepy.annotation.sanitize.write_annot": "genomepy.annotation.utils.write_annot",
        "genomepy.providers.base.write_annot": "genomepy.annotation.utils.write_annot",
        "genomepy.annotation.map_locations": "genomepy.providers.map_locations",
    }
    for alias, target in pairs.items():
        assert _resolve(alias) is _resolve(target), alias
    assert "annotation.read_annot" in inspect.getsource(gp_api.read_annotation_records)
    assert "annotation.write_annot" in inspect.getsource(gp_api.write_annotation_records)
