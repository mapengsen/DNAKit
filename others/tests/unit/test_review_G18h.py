"""Offline tests for the G18h review wrappers: pysam (samtools, bcftools) and hgvs.

Fixtures are small and built inside tmp_path. Expected values are hand calculations from the
fixture reads (see the comments), or come from the upstream call the wrapper makes.
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError

pysam = pytest.importorskip("pysam")

from dnakit.variants import hgvs_review_functions as hgvs_mod  # noqa: E402
from dnakit.variants import pysam_review_functions as pysam_mod  # noqa: E402

REF_CHR1 = "ACGT" * 30  # 120 bases
REF_CHR2 = "GGCC" * 20  # 80 bases

SAM_TEXT = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:120\n"
    "@SQ\tSN:chr2\tLN:80\n"
    "@RG\tID:rg1\tSM:sampleA\n"
    "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTTA\tIIIIIIIIII\tRG:Z:rg1\n"
    "r2\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTTTAGCC\tIIIIIIIIII\tRG:Z:rg1\n"
    "r3\t16\tchr1\t20\t30\t8M2I\t*\t0\t0\tGATCGATCGG\tIIIIIIIIII\tRG:Z:rg1\n"
    "r4\t0\tchr2\t3\t60\t10M\t*\t0\t0\tGGCCTTAAGG\tIIIIIIIIII\tRG:Z:rg1\n"
    "r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\n"
)


def _write_ref(tmp_path: Path) -> Path:
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr1\n{REF_CHR1}\n>chr2\n{REF_CHR2}\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    return fasta


def _sam_fixture(tmp_path: Path) -> Path:
    """Write the reads as BAM with upstream pysam; sorted by coordinate."""

    sam = tmp_path / "reads.sam"
    sam.write_text(SAM_TEXT, encoding="utf-8")
    _write_ref(tmp_path)
    bam = tmp_path / "reads.bam"
    bam.write_bytes(pysam.samtools.view("-b", str(sam), catch_stdout=True))
    return bam


def _ref(tmp_path: Path) -> Path:
    return tmp_path / "ref.fa"


def _indexed_bam(tmp_path: Path) -> Path:
    bam = _sam_fixture(tmp_path)
    pysam.samtools.index(str(bam), catch_stdout=True)
    return bam


def _vcf(
    tmp_path: Path, name: str, samples: tuple[str, str], rows: list[tuple[int, str, str, str]]
) -> Path:
    """Write a small VCF; rows hold (position, alt, filter:quality, genotypes)."""

    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=120>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples[0]}\t{samples[1]}",
    ]
    for position, alt, quality_filter, genotypes in rows:
        ref = REF_CHR1[position - 1]
        if alt == ref:
            alt = "C" if ref != "C" else "G"
        quality, filt = quality_filter.split(":")
        first, second = genotypes.split(",")
        lines.append(
            f"chr1\t{position}\t.\t{ref}\t{alt}\t{quality}\t{filt}\tAF=0.5\tGT:DP\t{first}:20\t{second}:18"
        )
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _gz_vcf(
    tmp_path: Path, name: str, samples: tuple[str, str], rows: list[tuple[int, str, str, str]]
) -> Path:
    plain = _vcf(tmp_path, name, samples, rows)
    pysam.tabix_index(str(plain), preset="vcf", force=True, keep_original=True)
    return Path(str(plain) + ".gz")


SITES = [
    (5, "G", "50:PASS", "0/1,1/1"),
    (30, "T", "60:PASS", "0/0,0/1"),
    (60, "C", "40:LowQual", "0/1,0/0"),
]


def test_all_public_functions_take_keyword_verbose() -> None:
    for module in (pysam_mod, hgvs_mod):
        for name in module.__all__:
            function = getattr(module, name)
            parameter = inspect.signature(function).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is True


def test_add_bam_read_group_sets_the_read_group(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    result = pysam_mod.add_bam_read_group(bam, tmp_path / "rg.bam", "@RG\tID:rgX\tSM:s9")
    assert result["read_group_id"] == "rgX"
    assert result["records"] == 5
    with pytest.raises(ConfigurationError):
        pysam_mod.add_bam_read_group(bam, tmp_path / "bad.bam", "@PG\tID:x")


def test_clip_amplicon_primers_keeps_every_read(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    bed = tmp_path / "primers.bed"
    bed.write_text("chr1\t0\t10\tamp1_LEFT\t60\t+\nchr1\t100\t110\tamp1_RIGHT\t60\t-\n")
    result = pysam_mod.clip_amplicon_primers(bam, bed, tmp_path / "clipped.bam")
    assert result["records"] == 5


def test_summarize_amplicon_stats_reads_the_amplicon_layout(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    bed = tmp_path / "primers.bed"
    bed.write_text("chr1\t0\t10\tamp1_LEFT\t60\t+\nchr1\t100\t110\tamp1_RIGHT\t60\t-\n")
    result = pysam_mod.summarize_amplicon_stats(bam, bed)
    assert result["count"] == 1
    assert result["amplicons"][0]["reference"] == "chr1"
    # The FSS rows are per reference that carries amplicons; chr1 holds r1, r2 and r3.
    metrics = {
        row["metric"]: row["value"] for row in result["summary"] if row["reference"] == "chr1"
    }
    assert metrics["raw total sequences"] == 3


def test_extract_fastq_and_fasta_split_reads_by_role(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    fastq = pysam_mod.extract_fastq_from_bam(bam, tmp_path / "fq")
    # Reads carry no READ1/READ2 flag, so all five go to the unpaired file.
    assert fastq["records"]["other"] == 5
    fasta = pysam_mod.extract_fasta_from_bam(bam, tmp_path / "fa")
    assert fasta["records"]["other"] == 5
    text = Path(fastq["files"]["other"]).read_text(encoding="utf-8")
    assert text.splitlines()[1] == "ACGTACGTTA"


def test_bam2fq_matches_fastq_output(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    upstream = pysam.samtools.bam2fq("-0", str(tmp_path / "b0.fq"), str(bam), catch_stdout=True)
    assert upstream == ""
    pysam_mod.extract_fastq_from_bam(bam, tmp_path / "fq2")
    assert (tmp_path / "b0.fq").read_bytes() == (tmp_path / "fq2" / "other.fastq").read_bytes()


def test_collate_and_concatenate_and_merge(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    assert pysam_mod.collate_bam_reads(bam, tmp_path / "col.bam")["records"] == 5
    assert pysam_mod.concatenate_bam_files([bam, bam], tmp_path / "cat.bam")["records"] == 10
    assert pysam_mod.merge_bam_files([bam, bam], tmp_path / "merged.bam")["records"] == 10


def test_bamshuf_order_matches_collate(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    pysam.samtools.bamshuf(str(bam), str(tmp_path / "shuf"), catch_stdout=True)
    pysam_mod.collate_bam_reads(bam, tmp_path / "col2.bam")
    shuffled = pysam.samtools.view("-h", str(tmp_path / "shuf.bam"), catch_stdout=True)
    collated = pysam.samtools.view("-h", str(tmp_path / "col2.bam"), catch_stdout=True)
    records = lambda text: [line for line in text.splitlines() if not line.startswith("@")]  # noqa: E731
    assert records(shuffled) == records(collated)


def test_bedcov_counts_overlapping_bases(tmp_path: Path) -> None:
    bam = _indexed_bam(tmp_path)
    bed = tmp_path / "intervals.bed"
    bed.write_text("chr1\t0\t10\nchr1\t0\t120\n")
    rows = pysam_mod.count_bam_reads_per_interval(bed, bam)["intervals"]
    # r1 covers 10 bases of [0,10); r2 covers 6 bases (5-10); total 16.
    # Over [0,120): 10 + 10 + 8 aligned reference bases = 28.
    assert [row["covered_bases"] for row in rows] == [16, 28]


def test_coverage_and_depth_hand_calculation(tmp_path: Path) -> None:
    bam = _indexed_bam(tmp_path)
    chr1 = next(
        row for row in pysam_mod.summarize_bam_coverage(bam)["contigs"] if row["rname"] == "chr1"
    )
    # Distinct covered bases: 1-14 (14) plus 20-27 (8) = 22.
    assert chr1["covbases"] == 22
    assert chr1["numreads"] == 3
    depth = pysam_mod.calculate_bam_read_depth(bam, region="chr1:1-10")["positions"]
    assert len(depth) == 10
    assert depth[0]["depth"] == 1  # only r1 at position 1
    assert depth[4]["depth"] == 2  # r1 and r2 at position 5


def test_dict_and_flags(tmp_path: Path) -> None:
    _sam_fixture(tmp_path)
    result = pysam_mod.write_fasta_sequence_dictionary(_ref(tmp_path), tmp_path / "ref.dict")
    assert [(item["name"], item["length"]) for item in result["contigs"]] == [
        ("chr1", 120),
        ("chr2", 80),
    ]
    flags = pysam_mod.decode_sam_flags([99, 0, 4])["flags"]
    assert flags[0]["names"] == ["MREVERSE", "PAIRED", "PROPER_PAIR", "READ1"]  # 99 = 1+2+32+64
    assert flags[1]["names"] == []
    assert flags[2]["names"] == ["UNMAP"]


def test_head_idxstats_index_and_stats(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    head = pysam_mod.read_bam_header_and_records(bam, records=2)
    assert head["header"][0].startswith("@HD")
    assert len(head["records"]) == 2
    index = pysam_mod.index_alignment_file(bam)
    assert index["format"] == "BAI"
    stats = pysam_mod.summarize_bam_index_stats(bam)
    chr1 = next(item for item in stats["contigs"] if item["contig"] == "chr1")
    assert chr1["mapped"] == 3 and stats["unplaced_unmapped"] == 1
    assert pysam_mod.index_alignment_file(bam, csi=True)["format"] == "CSI"
    assert pysam_mod.summarize_alignment_stats(bam)["summary"]["raw total sequences"] == 5


def test_index_refuses_to_overwrite_without_force(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    pysam_mod.index_alignment_file(bam)
    with pytest.raises(ConfigurationError):
        pysam_mod.index_alignment_file(bam)


def test_markdup_sort_and_list_samples(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    marked = pysam_mod.mark_bam_duplicates(bam, tmp_path / "md.bam")
    assert marked["duplicates"] == 0 and marked["records"] == 5
    assert (
        pysam_mod.sort_alignment_file(bam, tmp_path / "byname.bam", by_name=True)["order"]
        == "queryname"
    )
    assert pysam_mod.sort_alignment_file(bam, tmp_path / "bycoord.bam")["order"] == "coordinate"
    samples = pysam_mod.list_read_group_samples([bam])["samples"]
    assert samples[0]["sample"] == "sampleA"


def test_mpileup_positions(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    result = pysam_mod.pileup_alignment_positions(
        bam, tmp_path / "pile.txt", reference=_ref(tmp_path)
    )
    # chr1: 1-14 and 20-27 = 22 covered positions; chr2: 3-12 = 10.
    assert result["positions"] == 32


def test_fix_mates_needs_name_sorted_input(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    with pytest.raises(ConfigurationError):
        pysam_mod.fix_bam_mate_fields(bam, tmp_path / "fixed_coord.bam")
    byname = pysam_mod.sort_alignment_file(bam, tmp_path / "byname.bam", by_name=True)
    assert pysam_mod.fix_bam_mate_fields(byname["output"], tmp_path / "fixed.bam")["records"] == 5


def test_integrity_check_flags_bad_file(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    junk = tmp_path / "junk.bam"
    junk.write_text("not a bam", encoding="utf-8")
    result = pysam_mod.check_alignment_file_integrity([bam, junk])
    assert [item["valid"] for item in result["files"]] == [True, False]
    assert result["invalid"] == 1


def test_replace_header_and_reset_unmaps(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    header = tmp_path / "header.sam"
    header.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:120\n@SQ\tSN:chr2\tLN:80\n@RG\tID:new\tSM:sampleB\n"
    )
    replaced = pysam_mod.replace_bam_header(header, bam, tmp_path / "rh.bam")
    assert replaced["records"] == 5
    view = pysam.samtools.view("-h", replaced["output"], catch_stdout=True)
    assert "ID:new" in view
    reset = pysam_mod.reset_alignments_to_unmapped(bam, tmp_path / "reset.bam")
    assert reset["records"] == 5
    flags = {
        line.split("\t")[1]
        for line in pysam.samtools.view(reset["output"], catch_stdout=True).splitlines()
    }
    assert flags == {"4"}


def test_split_by_read_group_and_select_records(tmp_path: Path) -> None:
    bam = _indexed_bam(tmp_path)
    split = pysam_mod.split_bam_by_read_group(bam, tmp_path / "split")
    assert [item["file"] for item in split["files"]] == ["rg1.bam"]
    assert split["files"][0]["records"] == 5
    # Region chr1:1-10 with MAPQ >= 30: r1 and r2 overlap; r3 starts at 20.
    region = pysam_mod.select_alignment_records(
        bam, tmp_path / "sel.bam", region="chr1:1-10", min_mapping_quality=30
    )
    assert region["records"] == 2
    # Unmapped reads (flag 4) removed.
    no_unmapped = pysam_mod.select_alignment_records(bam, tmp_path / "sel.sam", excluded_flags=4)
    assert no_unmapped["records"] == 4


def test_region_query_needs_an_index(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    with pytest.raises(ConfigurationError):
        pysam_mod.select_alignment_records(bam, tmp_path / "x.bam", region="chr1:1-10")


def test_reconstruct_reference_from_md_tags(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    calmd = pysam_mod.recompute_nm_md_tags(bam, _ref(tmp_path), tmp_path / "calmd.bam")
    assert calmd["records"] == 5
    # Positions covered by aligned bases: chr1 has 22 (1-14, 20-27), chr2 has 10 (3-12).
    result = pysam_mod.reconstruct_reference_from_alignments(
        calmd["output"], tmp_path / "rebuilt.fa"
    )
    determined = {item["name"]: item["determined_bases"] for item in result["contigs"]}
    assert determined == {"chr1": 22, "chr2": 10}


def test_phase_writes_haplotype_files(tmp_path: Path) -> None:
    bam = _sam_fixture(tmp_path)
    result = pysam_mod.phase_bam_haplotypes(bam, tmp_path / "phased")
    names = {Path(item["file"]).name for item in result["files"]}
    assert names <= {"phased.0.bam", "phased.1.bam", "phased.chimera.bam"} and names
    with pytest.raises(ConfigurationError):
        pysam_mod.phase_bam_haplotypes(bam, tmp_path / "phased")


def test_fastq_index_and_import(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fq"
    fastq.write_text(
        "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nGGGGCCCCAA\n+\nIIIIIIIIII\n", encoding="utf-8"
    )
    assert pysam_mod.index_fastq_file(fastq)["entries"] == 2
    imported = pysam_mod.import_fastq_to_bam([fastq], tmp_path / "imported.bam")
    assert imported["records"] == 2 and imported["layout"] == "single"


def test_bcftools_filter_csq_head_and_gtcheck(tmp_path: Path) -> None:
    _write_ref(tmp_path)
    plain = _vcf(tmp_path, "calls.vcf", ("S1", "S2"), SITES)
    gz = _gz_vcf(tmp_path, "gzcalls.vcf", ("S1", "S2"), SITES)
    included = pysam_mod.filter_vcf_by_expression(gz, "QUAL>45", tmp_path / "f.vcf")
    assert included["records"] == 2
    excluded = pysam_mod.filter_vcf_by_expression(
        gz, 'FILTER="PASS"', tmp_path / "e.vcf", mode="exclude"
    )
    assert excluded["records"] == 1
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t120\t.\t+\t.\tID=gene1;Name=G1\n"
        "chr1\ttest\tmRNA\t1\t120\t.\t+\t.\tID=tx1;Parent=gene1;Name=T1\n"
        "chr1\ttest\tCDS\t1\t120\t.\t+\t0\tID=cds1;Parent=tx1\n",
        encoding="utf-8",
    )
    # bcftools csq stops on unphased heterozygous genotypes, so the test uses a phased copy.
    phased = tmp_path / "phased.vcf"
    phased.write_text(plain.read_text(encoding="utf-8").replace("/", "|"), encoding="utf-8")
    csq = pysam_mod.annotate_vcf_consequence(phased, gff, _ref(tmp_path), tmp_path / "csq.vcf")
    assert csq["annotated_records"] == 3
    header = pysam_mod.read_vcf_header(plain)
    assert header["samples"] == ["S1", "S2"]
    pairs = pysam_mod.check_vcf_sample_concordance(gz)["pairs"]
    self_pairs = [item for item in pairs if item["query"] == item["genotyped"]]
    assert all(item["discordance"] == 0 for item in self_pairs)
    assert importlib.import_module("pysam.bcftools")  # the submodule the wrappers call


def test_bcftools_index_merge_mpileup_rename_and_roh(tmp_path: Path) -> None:
    _write_ref(tmp_path)
    gz1 = _gz_vcf(tmp_path, "one.vcf", ("S1", "S2"), SITES)
    gz2 = _gz_vcf(tmp_path, "two.vcf", ("S3", "S4"), SITES)
    csi = pysam_mod.index_vcf_csi(gz1)
    assert csi["index"].endswith(".csi") and Path(csi["index"]).is_file()
    with pytest.raises(ConfigurationError):
        pysam_mod.index_vcf_csi(gz1)
    # Without an index the merge is refused; a CSI index on the second file satisfies it.
    Path(str(gz2) + ".tbi").unlink()
    with pytest.raises(ConfigurationError):
        pysam_mod.merge_indexed_vcfs([gz1, gz2], tmp_path / "unindexed.vcf")
    pysam_mod.index_vcf_csi(gz2)
    merged = pysam_mod.merge_indexed_vcfs([gz1, gz2], tmp_path / "merged.vcf")
    assert merged["samples"] == ["S1", "S2", "S3", "S4"]
    bam = _sam_fixture(tmp_path)
    pileup = pysam_mod.pileup_bam_to_vcf(bam, _ref(tmp_path), tmp_path / "pile.vcf")
    assert pileup["records"] > 0
    names = tmp_path / "names.txt"
    names.write_text("A1\nA2\n", encoding="utf-8")
    renamed = pysam_mod.rename_vcf_samples(
        _vcf(tmp_path, "rn.vcf", ("S1", "S2"), SITES), names, tmp_path / "rn2.vcf"
    )
    assert renamed["samples"] == ["A1", "A2"]
    assert "A1\tA2" in Path(renamed["output"]).read_text(encoding="utf-8")
    roh_rows = [(position, "A", "50:PASS", "1/1,0/0") for position in range(1, 40)]
    roh_vcf = _gz_vcf(tmp_path, "roh.vcf", ("S1", "S2"), roh_rows)
    # On this small fixture bcftools reports no runs; the call itself must still succeed.
    roh = pysam_mod.detect_bcftools_roh(roh_vcf, sample="S1")
    assert roh["sample"] == "S1" and roh["count"] == 0
    with pytest.raises(ConfigurationError):
        pysam_mod.detect_bcftools_roh(roh_vcf, sample="S1", genotype_quality=0)


def test_roh_parser_reads_region_rows(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Parse an RG row in bcftools roh's column layout (the upstream call is mocked here)."""

    text = (
        "# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\n"
        "RG\tS1\tchr1\t1\t60\t60\t12\t27.5\n"
    )

    def fake_command(program, subcommand, arguments):
        Path(arguments[arguments.index("-o") + 1]).write_text(text, encoding="utf-8")
        return ""

    monkeypatch.setattr(pysam_mod, "_command", fake_command)
    source = tmp_path / "x.vcf"
    source.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    result = pysam_mod.detect_bcftools_roh(source, sample="S1")
    assert result["runs"] == [
        {
            "sample": "S1",
            "chromosome": "chr1",
            "start": 1,
            "end": 60,
            "length": 60,
            "markers": 12,
            "quality": 27.5,
        }
    ]


def test_hgvs_trim_interbase_and_repeat_notation() -> None:
    trimmed = hgvs_mod.trim_common_suffix_strings(["GGCTA", "AGCTA"])
    assert trimmed["trimmed_length"] == 4 and trimmed["strings"] == ["G", "A"]
    assert hgvs_mod.hgvs_interbase_span("NC_000001.11:g.100_102del")["start"] == 99
    span = hgvs_mod.hgvs_interbase_span("NC_000001.11:g.100_102del")
    assert (span["start"], span["end"]) == (99, 102)
    block = [{"repeat_unit": "AC", "repeat_count": 3, "block_size": 6, "block": "ACACAC"}]
    assert hgvs_mod.repeat_notation_for_alleles("ACACAC", "ACAC", block)["formatted"] == "AC[3]"
    assert hgvs_mod.repeat_notation_for_alleles("ACGT", "ACGG", [])["formatted"] is None
    with pytest.raises(ConfigurationError):
        hgvs_mod.trim_common_suffix_strings(["A"])
    with pytest.raises(ConfigurationError):
        hgvs_mod.hgvs_interbase_span("not hgvs")
