"""Tests for the G5 wrappers: pysam, cyvcf2 and vcfpy variant and interval helpers."""

from __future__ import annotations

import gzip
import sys
import warnings
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import cyvcf2_functions, pysam_functions, vcfpy_functions

REFERENCE = "ACGTACGTACGTACGTTTGACCAAGGTTACGATCGATCGAACGTACGTAC"

VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=50>
##FILTER=<ID=q10,Description="Quality below 10">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t3\trs1\tG\tT\t50\tPASS\tDP=20;AF=0.5\tGT:AD:DP\t0/1:10,10:20\t1/1:0,12:12\t0/0:15,0:15
chr1\t10\t.\tC\tCA\t40\tPASS\tDP=10\tGT:AD:DP\t0/1:5,5:10\t0/0:10,0:10\t0/1:4,6:10
chr1\t22\trs3\tC\tG,A\t60\tq10\tDP=30;AF=0.3,0.1\tGT:AD:DP\t1|2:0,10,5:15\t0/0:20,0,0:20\t0/2:10,0,10:20
chr1\t35\trs4\tG\tC\t70\tPASS\tDP=12;AF=0.4\tGT:AD:DP\t0/1:6,6:12\t0/1:5,5:10\t./.:0,0:0
"""

ANNOTATION_TEXT = "chr1\t3\tG\tT\tcommon\nchr1\t22\tC\tG\trare\n"

INFO_HEADER = '##INFO=<ID=NOTE,Number=1,Type=String,Description="note">'


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _records(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


@pytest.fixture()
def workdir(tmp_path: Path) -> dict[str, Path]:
    """Reference FASTA with index, plain and indexed VCF, and an indexed annotation table."""

    pysam = pytest.importorskip("pysam")
    fasta = _write(
        tmp_path / "ref.fa",
        ">chr1\n" + "\n".join(REFERENCE[i : i + 10] for i in range(0, len(REFERENCE), 10)) + "\n",
    )
    pysam.faidx(str(fasta))
    plain = _write(tmp_path / "sample.vcf", VCF_TEXT)
    pysam.tabix_index(str(plain), preset="vcf", force=True, keep_original=True)
    annotation = _write(tmp_path / "annotation.tsv", ANNOTATION_TEXT)
    pysam.tabix_compress(str(annotation), str(tmp_path / "annotation.tsv.gz"), force=True)
    pysam.tabix_index(
        str(tmp_path / "annotation.tsv.gz"), seq_col=0, start_col=1, end_col=1, force=True
    )
    return {
        "dir": tmp_path,
        "fasta": fasta,
        "plain": plain,
        "indexed": tmp_path / "sample.vcf.gz",
        "annotation": tmp_path / "annotation.tsv.gz",
    }


# ---------------------------------------------------------------- pysam


def test_index_tabix_file_bgzips_and_indexes_for_region_queries(workdir: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    result = pysam_functions.index_tabix_file(workdir["plain"], force=True, verbose=False)
    assert result["compressed"] == str(workdir["indexed"])
    assert Path(result["index"]).is_file() and workdir["plain"].is_file()
    # Direct upstream region query: POS 10 and 22 fall in chr1:9-22 (0-based 8..22).
    with pysam.VariantFile(str(workdir["indexed"])) as handle:
        assert [record.pos for record in handle.fetch("chr1", 8, 22)] == [10, 22]
    with pytest.raises(ConfigurationError, match="already exists"):
        pysam_functions.index_tabix_file(workdir["plain"], verbose=False)


def test_compress_bgzip_file_round_trips_text_and_checks_suffix(workdir: dict[str, Path]) -> None:
    pytest.importorskip("pysam")
    target = workdir["dir"] / "copy.vcf.gz"
    result = pysam_functions.compress_bgzip_file(workdir["plain"], target, verbose=False)
    assert result["uncompressed_bytes"] == len(VCF_TEXT.encode())
    with gzip.open(target, "rt", encoding="utf-8") as handle:
        assert handle.read() == VCF_TEXT
    with pytest.raises(ConfigurationError, match=r"must end with \.gz"):
        pysam_functions.compress_bgzip_file(
            workdir["plain"], workdir["dir"] / "x.vcf", verbose=False
        )


def test_parse_genomic_region_converts_to_zero_based_half_open(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    result = pysam_functions.parse_genomic_region("chr1:11-20", verbose=False)
    assert (result["reference"], result["start"], result["end"], result["length"]) == (
        "chr1",
        10,
        20,
        10,
    )
    direct = pysam.libcutils.parse_region(region="chr1:11-20")
    assert direct == ("chr1", 10, 20)
    whole = pysam_functions.parse_genomic_region("chr2", verbose=False)
    assert (whole["start"], whole["end"]) == (0, None)
    comma = pysam_functions.parse_genomic_region("chr1:1,001-1,020", verbose=False)
    assert (comma["start"], comma["end"]) == (1000, 1020)
    for bad in ("chr1:11", "chr1:0-5", "chr1:20-11", "chr1:x-5"):
        with pytest.raises(ConfigurationError):
            pysam_functions.parse_genomic_region(bad, verbose=False)


def test_parse_genomic_regions_is_strict_where_upstream_is_not() -> None:
    pytest.importorskip("pysam")
    result = pysam_functions.parse_genomic_regions("chr1:10-20,chr2", verbose=False)
    assert result["regions"] == [
        {"reference": "chr1", "start": 9, "end": 20},
        {"reference": "chr2", "start": 0, "end": None},
    ]
    # pysam.libcvcf.parse_regions silently keeps "chr1:10-20" from "chr1:10-20-30".
    with pytest.raises(ConfigurationError, match="region must look like"):
        pysam_functions.parse_genomic_regions("chr1:10-20-30", verbose=False)


def test_fetch_padded_sequence_pads_outside_the_contig(workdir: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    right = pysam_functions.fetch_padded_sequence(workdir["fasta"], "chr1", 45, 60, verbose=False)
    assert right["sequence"] == "CGTACNNNNNNNNNN"
    assert (right["left_padding"], right["right_padding"]) == (0, 10)
    left = pysam_functions.fetch_padded_sequence(workdir["fasta"], "chr1", -3, 4, verbose=False)
    assert left["sequence"] == "NNNACGT"
    assert (left["left_padding"], left["right_padding"]) == (3, 0)
    with pysam.FastaFile(str(workdir["fasta"])) as handle:
        assert handle.fetch("chr1", 0, 5).upper() == "ACGTA"
    with pytest.raises(ConfigurationError, match="is not in"):
        pysam_functions.fetch_padded_sequence(workdir["fasta"], "chrX", 0, 4, verbose=False)


def test_query_vcf_fields_returns_text_rows_and_genotypes(workdir: dict[str, Path]) -> None:
    pytest.importorskip("pysam")
    result = pysam_functions.query_vcf_fields(
        workdir["plain"], ["POS", "REF", "ALT", "FILTER", "INFO/DP"], verbose=False
    )
    assert result["rows"][:2] == [["3", "G", "T", "PASS", "20"], ["10", "C", "CA", "PASS", "10"]]
    assert result["record_count"] == 4
    indexed = pysam_functions.query_vcf_fields(
        workdir["indexed"], ["POS"], region="chr1:9-22", genotypes=True, verbose=False
    )
    assert indexed["columns"] == ["POS", "GT:S1", "GT:S2", "GT:S3"]
    assert indexed["rows"] == [["10", "0/1", "0/0", "0/1"], ["22", "1|2", "0/0", "0/2"]]
    with pytest.raises(ConfigurationError, match="unsupported field"):
        pysam_functions.query_vcf_fields(workdir["plain"], ["FORMAT/DP"], verbose=False)
    with pytest.raises(ConfigurationError, match="needs a tabix"):
        pysam_functions.query_vcf_fields(
            workdir["plain"], ["POS"], region="chr1:9-22", verbose=False
        )


def test_select_vcf_records_subsets_samples_region_and_pass(workdir: dict[str, Path]) -> None:
    pytest.importorskip("pysam")
    target = workdir["dir"] / "subset.vcf"
    result = pysam_functions.select_vcf_records(
        workdir["indexed"],
        target,
        samples=["S2", "S1"],
        region="chr1:9-22",
        pass_only=True,
        verbose=False,
    )
    # POS 10 is PASS; POS 22 carries FILTER q10 and is dropped.
    assert result["records"] == 1
    header = [line for line in target.read_text().splitlines() if line.startswith("#CHROM")]
    assert header[0].split("\t")[9:] == ["S2", "S1"]
    assert _records(target.read_text())[0][1] == "10"
    with pytest.raises(ConfigurationError, match="already exists"):
        pysam_functions.select_vcf_records(workdir["indexed"], target, verbose=False)


def test_concatenate_vcf_files_restores_split_records(workdir: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    lines = VCF_TEXT.splitlines(keepends=True)
    header = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    first = _write(workdir["dir"] / "part1.vcf", "".join(header + body[:2]))
    second = _write(workdir["dir"] / "part2.vcf", "".join(header + body[2:]))
    for part in (first, second):
        pysam.tabix_index(str(part), preset="vcf", force=True, keep_original=True)
    target = workdir["dir"] / "joined.vcf"
    result = pysam_functions.concatenate_vcf_files(
        [workdir["dir"] / "part1.vcf.gz", workdir["dir"] / "part2.vcf.gz"], target, verbose=False
    )
    assert result["records"] == 4
    assert [row[1] for row in _records(target.read_text())] == ["3", "10", "22", "35"]


def test_intersect_vcf_sites_counts_private_and_shared(workdir: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    extra = VCF_TEXT + "chr1\t40\t.\tA\tT\t30\tPASS\t.\tGT\t0/1\t0/1\t0/1\n"
    second_plain = _write(workdir["dir"] / "more.vcf", extra)
    pysam.tabix_index(str(second_plain), preset="vcf", force=True, keep_original=True)
    result = pysam_functions.intersect_vcf_sites(
        workdir["indexed"], workdir["dir"] / "more.vcf.gz", workdir["dir"] / "isec", verbose=False
    )
    assert (result["private_first"], result["private_second"], result["shared"]) == (0, 1, 4)
    assert "sites.txt" in result["files"]
    with pytest.raises(ConfigurationError, match="already exists"):
        pysam_functions.intersect_vcf_sites(
            workdir["indexed"],
            workdir["dir"] / "more.vcf.gz",
            workdir["dir"] / "isec",
            verbose=False,
        )


def test_sort_vcf_records_orders_unsorted_input(workdir: dict[str, Path]) -> None:
    pytest.importorskip("pysam")
    lines = VCF_TEXT.splitlines(keepends=True)
    header = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    shuffled = _write(
        workdir["dir"] / "shuffled.vcf", "".join([*header, body[3], body[0], body[2], body[1]])
    )
    target = workdir["dir"] / "sorted.vcf"
    result = pysam_functions.sort_vcf_records(shuffled, target, verbose=False)
    assert result["records"] == 4
    assert [row[1] for row in _records(target.read_text())] == ["3", "10", "22", "35"]


def test_annotate_vcf_from_tabix_copies_info_tag(workdir: dict[str, Path]) -> None:
    pytest.importorskip("pysam")
    target = workdir["dir"] / "annotated.vcf"
    result = pysam_functions.annotate_vcf_from_tabix(
        workdir["plain"],
        workdir["annotation"],
        ["CHROM", "POS", "REF", "ALT", "INFO/NOTE"],
        target,
        header_lines=[INFO_HEADER],
        verbose=False,
    )
    assert result["annotated"] == {"NOTE": 2}
    notes = {row[1]: row[7] for row in _records(target.read_text())}
    assert "NOTE=common" in notes["3"] and "NOTE=rare" in notes["22"]
    assert "NOTE" not in notes["10"]
    with pytest.raises(ConfigurationError, match="NOTE"):
        pysam_functions.annotate_vcf_from_tabix(
            workdir["plain"],
            workdir["annotation"],
            ["CHROM", "POS", "REF", "ALT", "INFO/NOTE"],
            workdir["dir"] / "no_header.vcf",
            verbose=False,
        )


def _edited(reference: str, edits: dict[int, str]) -> str:
    """Replace 1-based single-base positions with the given bases."""

    chars = list(reference)
    for position, base in edits.items():
        chars[position - 1] = base
    return "".join(chars)


def test_build_vcf_consensus_applies_one_haplotype(workdir: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    # S2 haplotype 1: 1/1 at 3 gives T; 0/0 elsewhere; 0/1 at 35 keeps the reference on hap 1.
    s2_hap1 = pysam_functions.build_vcf_consensus(
        workdir["indexed"],
        workdir["fasta"],
        "S2",
        workdir["dir"] / "s2h1.fa",
        haplotype=1,
        verbose=False,
    )
    assert s2_hap1["sequences"] == {"chr1": 50}
    expected_s2_h1 = _edited(REFERENCE, {3: "T"})
    assert _sequence_of(workdir["dir"] / "s2h1.fa") == expected_s2_h1
    # S2 haplotype 2 also carries the 0/1 allele C at 35.
    pysam_functions.build_vcf_consensus(
        workdir["indexed"],
        workdir["fasta"],
        "S2",
        workdir["dir"] / "s2h2.fa",
        haplotype=2,
        verbose=False,
    )
    assert _sequence_of(workdir["dir"] / "s2h2.fa") == _edited(REFERENCE, {3: "T", 35: "C"})
    # S1 haplotype 2: T at 3, CA insertion after 10, A at 22 (allele 2 of 1|2), C at 35.
    s1_hap2 = pysam_functions.build_vcf_consensus(
        workdir["indexed"],
        workdir["fasta"],
        "S1",
        workdir["dir"] / "s1h2.fa",
        haplotype=2,
        verbose=False,
    )
    snv_only = _edited(REFERENCE, {3: "T", 22: "A", 35: "C"})
    # The CA allele inserts A after position 10, so every later base shifts by one.
    expected_s1 = snv_only[:10] + "A" + snv_only[10:]
    assert _sequence_of(workdir["dir"] / "s1h2.fa") == expected_s1
    assert s1_hap2["sequences"] == {"chr1": 51}
    # Direct upstream call gives the same sequence for S1 haplotype 2.
    upstream = pysam.bcftools.consensus(
        "-f", str(workdir["fasta"]), "-s", "S1", "-H", "2", str(workdir["indexed"])
    )
    assert "".join(upstream.splitlines()[1:]) == expected_s1
    with pytest.raises(ConfigurationError, match="bcftools consensus failed"):
        pysam_functions.build_vcf_consensus(
            workdir["indexed"], workdir["fasta"], "ZZ", workdir["dir"] / "bad.fa", verbose=False
        )


def _sequence_of(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def test_pysam_module_requires_bio_reads_extra(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    source = _write(tmp_path / "plain.vcf", VCF_TEXT)
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(BackendUnavailableError, match="bio-reads"):
        pysam_functions.index_tabix_file(source, verbose=False)


# ---------------------------------------------------------------- cyvcf2


def _indexed_vcf(tmp_path: Path) -> Path:
    pysam = pytest.importorskip("pysam")
    plain = _write(tmp_path / "sample.vcf", VCF_TEXT)
    pysam.tabix_index(str(plain), preset="vcf", force=True, keep_original=True)
    return tmp_path / "sample.vcf.gz"


def test_sample_heterozygosity_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("cyvcf2")
    path = _indexed_vcf(tmp_path)
    result = cyvcf2_functions.summarize_sample_heterozygosity(
        path, ["S1", "S2", "S3"], [("chr1", 3), ("chr1", 10)], verbose=False
    )
    s1, s2, s3 = (result["samples"][name] for name in ("S1", "S2", "S3"))
    # S1 is heterozygous at both sites; S2 at neither; S3 only at POS 10.
    assert (s1["het_count"], s1["het_ratio"], s1["mean_depth"]) == (2, 1.0, 15.0)
    assert (s2["het_count"], s2["het_ratio"]) == (0, 0.0)
    assert (s3["het_count"], s3["het_ratio"]) == (1, 0.5)
    assert s1["p10"] == s1["p90"] == 0.5  # alt fraction 10/20 and 5/10
    assert s3["p90"] == pytest.approx(0.6)  # alt fraction 6/10 at POS 10
    assert result["sites_used"] == ["chr1:10:C:CA", "chr1:3:G:T"]


def test_sample_heterozygosity_agrees_with_direct_het_check(tmp_path: Path) -> None:
    cyvcf2 = pytest.importorskip("cyvcf2")
    path = _indexed_vcf(tmp_path)
    sites = [("chr1", 3), ("chr1", 10)]
    result = cyvcf2_functions.summarize_sample_heterozygosity(
        path, ["S1", "S2", "S3"], sites, verbose=False
    )
    direct = cyvcf2.VCF(str(path), samples=["S1", "S2", "S3"], gts012=True).het_check(
        min_depth=8, percentiles=(10, 90), sites=[list(site) for site in sites]
    )[0]
    for name in ("S1", "S2", "S3"):
        assert result["samples"][name]["het_count"] == int(direct[name]["het_count"])
        assert result["samples"][name]["het_ratio"] == pytest.approx(direct[name]["het_ratio"])


def test_pairwise_relatedness_agrees_with_direct_site_relatedness(tmp_path: Path) -> None:
    cyvcf2 = pytest.importorskip("cyvcf2")
    path = _indexed_vcf(tmp_path)
    sites = [("chr1", 3, "G", "T"), ("chr1", 10, "C", "CA")]
    result = cyvcf2_functions.summarize_pairwise_relatedness(
        path, ["S1", "S2", "S3"], sites, verbose=False
    )
    with warnings.catch_warnings():
        # The direct call warns about S2 having no heterozygous sites, as the wrapper reports.
        warnings.simplefilter("ignore", UserWarning)
        direct = cyvcf2.VCF(str(path), samples=["S1", "S2", "S3"], gts012=True).site_relatedness(
            sites=[list(site) for site in sites], min_depth=5
        )
    assert [(p["sample_a"], p["sample_b"]) for p in result["pairs"]] == [
        ("S1", "S2"),
        ("S1", "S3"),
        ("S2", "S3"),
    ]
    assert [p["ibs2"] for p in result["pairs"]] == [int(v) for v in direct["ibs2"]]
    assert [p["relatedness"] for p in result["pairs"]] == pytest.approx(
        [float(v) for v in direct["rel"]]
    )
    assert all(p["informative_sites"] == 2 for p in result["pairs"])


def test_pearson_genotype_correlation_hand_values() -> None:
    cyvcf2 = pytest.importorskip("cyvcf2")
    np = pytest.importorskip("numpy")
    assert cyvcf2_functions.pearson_genotype_correlation([0, 1, 2, 1], [0, 2, 1, 1], verbose=False)[
        "correlation"
    ] == pytest.approx(0.5)
    assert cyvcf2_functions.pearson_genotype_correlation([0, 1, 2], [2, 1, 0], verbose=False)[
        "correlation"
    ] == pytest.approx(-1.0)
    # Pair 3 is skipped; the remaining pairs give r = sqrt(3) / 2.
    with_missing = cyvcf2_functions.pearson_genotype_correlation(
        [0, 1, None, 1], [0, 2, 1, 1], verbose=False
    )
    assert with_missing["correlation"] == pytest.approx(3**0.5 / 2)
    assert with_missing["missing_pairs"] == 1
    assert (
        cyvcf2_functions.pearson_genotype_correlation([0, 0, 0], [0, 1, 2], verbose=False)[
            "correlation"
        ]
        is None
    )
    direct = cyvcf2.r_unphased(
        np.array([0, 1, 2, 1], dtype=np.int32), np.array([0, 2, 1, 1], dtype=np.int32), 0.0, 4
    )
    assert direct == pytest.approx(0.5)


def test_cyvcf2_validation_and_missing_package(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    with pytest.raises(ConfigurationError, match="genotype dosages"):
        cyvcf2_functions.pearson_genotype_correlation([0, 3], [0, 1], verbose=False)
    plain = _write(tmp_path / "plain.vcf", VCF_TEXT)
    with pytest.raises(ConfigurationError, match="bgzip-compressed"):
        cyvcf2_functions.summarize_sample_heterozygosity(
            plain, ["S1"], [("chr1", 3)], verbose=False
        )
    monkeypatch.setitem(sys.modules, "cyvcf2", None)
    with pytest.raises(BackendUnavailableError, match="bio-cyvcf2"):
        cyvcf2_functions.pearson_genotype_correlation([0, 1], [0, 1], verbose=False)


# ---------------------------------------------------------------- vcfpy


def test_classify_vcf_alt_covers_allele_classes() -> None:
    pytest.importorskip("vcfpy")
    cases = {
        ("G", "T"): ("substitution", "SNV", "T"),
        ("CAT", "C"): ("substitution", "DEL", "C"),
        ("C", "CA"): ("substitution", "INS", "CA"),
        ("ACG", "TTT"): ("substitution", "MNV", "TTT"),
        ("C", "<DEL>"): ("symbolic", None, "DEL"),
    }
    for (ref, alt), (kind, subtype, value) in cases.items():
        result = vcfpy_functions.classify_vcf_alt(ref, alt, verbose=False)
        assert (result["kind"], result["subtype"], result["value"]) == (kind, subtype, value)
    breakend = vcfpy_functions.classify_vcf_alt("C", "G]17:198982]", verbose=False)
    assert breakend["kind"] == "breakend"
    assert breakend["breakend"]["mate_chrom"] == "17"
    assert breakend["breakend"]["mate_pos"] == 198982
    single = vcfpy_functions.classify_vcf_alt("C", ".A", verbose=False)
    assert single["kind"] == "single_breakend"
    with pytest.raises(ConfigurationError, match="VCF ALT"):
        vcfpy_functions.classify_vcf_alt("C", ".", verbose=False)


def test_classify_vcf_substitution_matches_upstream_process_sub() -> None:
    vcfpy = pytest.importorskip("vcfpy")
    for ref, alt in [("G", "T"), ("AT", "A"), ("C", "CA"), ("AC", "G"), ("ACG", "TTT")]:
        result = vcfpy_functions.classify_vcf_substitution(ref, alt, verbose=False)
        upstream = vcfpy.parser.process_sub(ref, alt)
        assert (result["subtype"], result["value"]) == (upstream.type, upstream.value)
    assert vcfpy_functions.classify_vcf_substitution("G", "T", verbose=False)["subtype"] == "SNV"


def test_parse_vcf_breakend_hand_values() -> None:
    pytest.importorskip("vcfpy")
    result = vcfpy_functions.parse_vcf_breakend("]13:123456]T", verbose=False)
    assert (result["mate_chrom"], result["mate_pos"]) == ("13", 123456)
    assert (result["orientation"], result["mate_orientation"]) == ("+", "-")
    assert (result["sequence"], result["within_main_assembly"]) == ("T", True)
    assert vcfpy_functions.parse_vcf_breakend("G]17:198982]", verbose=False)["orientation"] == "-"
    with pytest.raises(ConfigurationError, match="bracket notation"):
        vcfpy_functions.parse_vcf_breakend("A", verbose=False)
    with pytest.raises(ConfigurationError, match="failed"):
        vcfpy_functions.parse_vcf_breakend("G]17]", verbose=False)


def test_vcf_field_values_parse_and_format_round_trip() -> None:
    pytest.importorskip("vcfpy")
    assert vcfpy_functions.parse_vcf_field_value(
        "0.5,.", value_type="Float", number="A", verbose=False
    )["value"] == [0.5, None]
    assert (
        vcfpy_functions.parse_vcf_field_value("12", value_type="Integer", number=1, verbose=False)[
            "value"
        ]
        == 12
    )
    assert (
        vcfpy_functions.parse_vcf_field_value(".", value_type="Integer", number=1, verbose=False)[
            "value"
        ]
        is None
    )
    assert (
        vcfpy_functions.format_vcf_field_value(
            [0.5, 0.25], value_type="Float", number="A", verbose=False
        )["text"]
        == "0.5,0.25"
    )
    assert (
        vcfpy_functions.format_vcf_field_value(12, value_type="Integer", number=1, verbose=False)[
            "text"
        ]
        == "12"
    )
    assert (
        vcfpy_functions.format_vcf_field_value(None, value_type="Integer", number=1, verbose=False)[
            "text"
        ]
        == "."
    )
    with pytest.raises(ConfigurationError, match="cannot be converted"):
        vcfpy_functions.parse_vcf_field_value("x", value_type="Integer", number=1, verbose=False)
    with pytest.raises(ConfigurationError, match="does not match"):
        vcfpy_functions.format_vcf_field_value("x", value_type="Integer", number=1, verbose=False)


def test_remove_vcf_header_lines_matches_upstream(tmp_path: Path) -> None:
    vcfpy = pytest.importorskip("vcfpy")
    path = _write(tmp_path / "sample.vcf", VCF_TEXT)
    result = vcfpy_functions.remove_vcf_header_lines(path, [("FILTER", "q10")], verbose=False)
    assert result["removed_lines"] == 1
    assert not any(line.startswith("##FILTER=<ID=q10") for line in result["header_lines"])
    assert result["samples"] == ["S1", "S2", "S3"]
    reader = vcfpy.Reader.from_path(str(path))
    try:
        upstream = vcfpy.header_without_lines(reader.header, [("FILTER", "q10")])
    finally:
        reader.close()
    assert [line.serialize() for line in upstream.lines] == result["header_lines"]


def test_vcfpy_module_requires_bio_vcfpy_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "vcfpy", None)
    with pytest.raises(BackendUnavailableError, match="bio-vcfpy"):
        vcfpy_functions.classify_vcf_alt("G", "T", verbose=False)
