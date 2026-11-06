"""Tests for dnakit.variants.bcftools_plugins_p06_cli_functions (bcftools plugin facades).

Each function is checked against the same plugin run directly on the command line and, where the answer is
easy to derive, against a value computed in Python from the input.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from _p06_htslib_data import (
    DATA,
    PED_TEXT,
    body,
    data_vcf_gz,
    requires_bcftools,
    requires_data,
    requires_plugins,
    trio_vcf_text,
    write_vcf,
)

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import bcftools_plugins_p06_cli_functions as m

pytestmark = [requires_bcftools, requires_plugins, requires_data]


def _direct(*arguments: str) -> str:
    return subprocess.run(
        ["bcftools", *arguments], check=True, capture_output=True, text=True
    ).stdout


def _records(text: str) -> list[str]:
    return [line for line in text.splitlines() if line and not line.startswith("#")]


def _same_records(path: Path, direct: str) -> None:
    assert body(path) == _records(direct)


@pytest.fixture
def variants(tmp_path: Path) -> Path:
    return data_vcf_gz(tmp_path)


@pytest.fixture
def trio(tmp_path: Path) -> Path:
    return write_vcf(tmp_path / "trio.vcf", trio_vcf_text(), compress=True, index=True)


def _plain_rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in body(path)]


def test_ad_bias_reports_the_fisher_test_hit(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_ad_bias(
        trio, tmp_path / "bias.txt", sample_pairs=[["F", "M"]], verbose=False
    )
    hits = [row for row in result["report"]["rows"] if row[0] == "FT"]
    assert hits and hits[0][1:7] == ["F", "M", "chr1", 1000, "A", "G"]
    assert hits[0][-1] == pytest.approx(5.832134e-06, rel=1e-4)
    pairs = tmp_path / "pairs.txt"
    pairs.write_text("F\tM\n", encoding="utf-8")
    direct = [
        line
        for line in _direct("+ad-bias", str(trio), "--", "-s", str(pairs)).splitlines()
        if line.startswith("FT")
    ]
    assert [
        row
        for row in (tmp_path / "bias.txt").read_text(encoding="utf-8").splitlines()
        if row.startswith("FT")
    ] == direct


def test_ad_bias_clean_vcf_and_validation(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_ad_bias(
        trio,
        tmp_path / "clean.vcf",
        sample_pairs=[["F", "M"]],
        clean_vcf=True,
        variant_type="snp",
        verbose=False,
    )
    pairs = tmp_path / "pairs.txt"
    pairs.write_text("F\tM\n", encoding="utf-8")
    _same_records(
        tmp_path / "clean.vcf",
        _direct("+ad-bias", str(trio), "--", "-s", str(pairs), "-v", "snp", "-c"),
    )
    assert result["records"] == len(body(tmp_path / "clean.vcf"))
    with pytest.raises(ConfigurationError):
        m.bcftools_ad_bias(trio, tmp_path / "x.txt", sample_pairs=[["F"]], verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_ad_bias(
            trio, tmp_path / "x.txt", sample_pairs=[["F", "M"]], variant_type="mnp", verbose=False
        )


def test_af_distribution_matches_the_command_line(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_af_distribution(
        variants,
        tmp_path / "af.txt",
        deviation_bins=[0, 0.5, 1],
        probability_bins=[0, 0.5, 1],
        verbose=False,
    )
    rows = result["report"]["rows"]
    assert [row[0] for row in rows] == ["PROB_DIST", "PROB_DIST", "DEV_DIST", "DEV_DIST"]
    direct = _direct(
        "+af-dist", str(variants), "--", "-d", "0.0,0.5,1.0", "-p", "0.0,0.5,1.0", "-t", "AF"
    )
    assert [
        r
        for r in (tmp_path / "af.txt").read_text(encoding="utf-8").splitlines()
        if not r.startswith("#")
    ] == [r for r in direct.splitlines() if not r.startswith("#")]
    with pytest.raises(ConfigurationError):
        m.bcftools_af_distribution(
            variants, tmp_path / "bad.txt", deviation_bins=[0.5], verbose=False
        )


def test_allele_length_counts_agree_with_python(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_allele_length_counts(variants, tmp_path / "len.txt", verbose=False)
    assert result["report"]["columns"][:3] == ["LENGTH", "REF", "ALT"]
    rows = (
        _plain_rows(DATA / "variants.vcf")
        if False
        else [
            line.split("\t")
            for line in (DATA / "variants.vcf").read_text(encoding="utf-8").splitlines()
            if not line.startswith("#")
        ]
    )
    ref_counts: dict[int, int] = {}
    alt_counts: dict[int, int] = {}
    for row in rows:
        ref_counts[len(row[3])] = ref_counts.get(len(row[3]), 0) + 1
        alt_counts[len(row[4])] = alt_counts.get(len(row[4]), 0) + 1
    table = {row[0]: row for row in result["report"]["rows"]}
    for length, count in ref_counts.items():
        assert table[length][1] == count
    for length, count in alt_counts.items():
        assert table[length][2] == count


def test_check_ploidy_reports_every_sample(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_check_ploidy(variants, tmp_path / "ploidy.txt", verbose=False)
    rows = result["report"]["rows"]
    assert len({row[0] for row in rows}) == 24 and {row[4] for row in rows} == {2}
    assert (tmp_path / "ploidy.txt").read_text(encoding="utf-8") == _direct(
        "+check-ploidy", str(variants)
    )
    with_missing = m.bcftools_check_ploidy(
        variants, tmp_path / "ploidy_m.txt", use_missing=True, verbose=False
    )
    assert with_missing["status"] == "completed"


def test_check_sparsity_lists_samples_without_genotypes(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFULL\tEMPTY\n"
        "chr1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t./.\nchr1\t20\t.\tA\tG\t.\t.\t.\tGT\t1/1\t./.\n"
    )
    source = write_vcf(tmp_path / "sparse.vcf", text, compress=True, index=True)
    result = m.bcftools_check_sparsity(source, tmp_path / "sparse.txt", verbose=False)
    assert "EMPTY" in (tmp_path / "sparse.txt").read_text(encoding="utf-8")
    assert "FULL" not in (tmp_path / "sparse.txt").read_text(encoding="utf-8")
    assert result["status"] == "completed"
    with pytest.raises(ConfigurationError):
        m.bcftools_check_sparsity(source, tmp_path / "bad.txt", min_markers=0, verbose=False)


def test_genotype_dosage_from_gt(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_genotype_dosage(
        variants, tmp_path / "dosage.txt", tags=["GT"], verbose=False
    )
    source_rows = [
        line.split("\t")
        for line in (DATA / "variants.vcf").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    expected = []
    for row in source_rows:
        values = []
        for call in row[9:]:
            alleles = call.split(":")[0].replace("|", "/").split("/")
            values.append(-1.0 if "." in alleles else float(sum(int(a) > 0 for a in alleles)))
        expected.append(values)
    got = [row[4:] for row in result["report"]["rows"]]
    assert got == expected
    with pytest.raises(ConfigurationError):
        m.bcftools_genotype_dosage(variants, tmp_path / "bad.txt", tags=[], verbose=False)


def test_genotype_intersections_counts_all_sites(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_genotype_intersections(
        trio, tmp_path / "gtisec.txt", annotate_subsets=True, verbose=False
    )
    assert (tmp_path / "gtisec.txt").read_text(encoding="utf-8").splitlines()[3:] == _direct(
        "+GTisec", str(trio), "--", "-v"
    ).splitlines()[3:]
    assert result["report"]["rows"]
    human = m.bcftools_genotype_intersections(
        trio, tmp_path / "human.txt", human_readable=True, verbose=False
    )
    assert human["status"] == "completed"


def test_indel_stats_counts_indels(tmp_path: Path) -> None:
    source = write_vcf(
        tmp_path / "indel.vcf", trio_vcf_text(indels=True), compress=True, index=True
    )
    result = m.bcftools_indel_stats(source, tmp_path / "indel.txt", verbose=False)
    sn = [row for row in result["report"]["rows"] if str(row[0]).startswith("SN")]
    assert sn and sn[0][1] == 4  # four samples
    direct = _direct("+indel-stats", "-c", "CSQ", "--max-len", "20", str(source))
    assert [
        r
        for r in (tmp_path / "indel.txt").read_text(encoding="utf-8").splitlines()
        if not r.startswith("#")
    ] == [r for r in direct.splitlines() if not r.startswith("#")]
    ped = tmp_path / "trio.ped"
    ped.write_text(PED_TEXT, encoding="utf-8")
    de_novo = m.bcftools_indel_stats(source, tmp_path / "dnm.txt", ped_path=ped, verbose=False)
    assert de_novo["status"] == "completed"
    with pytest.raises(ConfigurationError):
        m.bcftools_indel_stats(source, tmp_path / "bad.txt", max_length=0, verbose=False)


def test_vcf_to_table_prints_the_variant_tables(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_vcf_to_table(
        variants,
        tmp_path / "table.txt",
        hide=["HOM_REF", "NO_CALL"],
        regions="chr1:392-392",
        verbose=False,
    )
    assert result["report"]["variants"] == 1
    text = (tmp_path / "table.txt").read_text(encoding="utf-8")
    assert "<<< chr1:392:T" in text and "| ID     | var0024 |" in text
    assert text == _direct(
        "+vcf2table", "-r", "chr1:392-392", str(variants), "--", "--hide", "HOM_REF,NO_CALL"
    )


def test_vrfs_scores_sites_and_validates(tmp_path: Path) -> None:
    reads = DATA / "reads.bam"
    reference = DATA / "genome.fa"
    if not reads.exists() or not reference.exists():
        pytest.skip("shared reads.bam / genome.fa missing")
    listing = tmp_path / "bams.txt"
    listing.write_text(f"{reads}\n", encoding="utf-8")
    sites = tmp_path / "sites.txt"
    sites.write_text("chr1\t392\tT\tC\nchr1\t859\tG\tT\n", encoding="utf-8")
    result = m.bcftools_variant_read_frequency_score(
        tmp_path / "scores.txt",
        alignment_list=listing,
        reference_fasta=reference,
        sites_file=sites,
        min_depth=1,
        verbose=False,
    )
    assert result["status"] == "completed"
    direct = _direct(
        "+vrfs", "-a", str(listing), "-f", str(reference), "-s", str(sites), "-d", "1", "-n", "20"
    )
    assert (tmp_path / "scores.txt").read_text(encoding="utf-8").splitlines()[
        3:
    ] == direct.splitlines()[3:]
    with pytest.raises(ConfigurationError):
        m.bcftools_variant_read_frequency_score(tmp_path / "x.txt", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_variant_read_frequency_score(
            tmp_path / "y.txt", alignment_list=listing, merge_batch_files=[sites], verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_variant_read_frequency_score(
            tmp_path / "z.txt",
            alignment_list=listing,
            reference_fasta=reference,
            sites_file=sites,
            batch="x",
            verbose=False,
        )


def test_add_variantkey_and_tables(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_add_variantkey(variants, tmp_path / "vk.vcf", verbose=False)
    assert result["records"] == 400
    first = _plain_rows(tmp_path / "vk.vcf")[0]
    assert "VKX=" in first[7] and "RSX=" in first[7]
    _same_records(tmp_path / "vk.vcf", _direct("+add-variantkey", str(variants)))
    tables = m.bcftools_variantkey_tables(variants, tmp_path / "keys", verbose=False)
    assert tables["files"] == ["nrvk.unsorted.tsv", "rsvk.unsorted.hex", "vkrs.unsorted.hex"]
    assert tables["counters"]["VariantKeys"] == 400
    with pytest.raises(ConfigurationError):
        m.bcftools_variantkey_tables(variants, tmp_path / "keys", verbose=False)


def test_contrast_groups_adds_association_tags(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_contrast_groups(
        trio,
        tmp_path / "contrast.vcf",
        control_samples=["F", "M"],
        case_samples=["C"],
        verbose=False,
    )
    assert result["records"] == 40
    _same_records(
        tmp_path / "contrast.vcf", _direct("+contrast", "-0", "F,M", "-1", "C", str(trio))
    )
    assert "PASSOC=" in _plain_rows(tmp_path / "contrast.vcf")[0][7]
    with pytest.raises(ConfigurationError):
        m.bcftools_contrast_groups(
            trio,
            tmp_path / "bad.vcf",
            control_samples=["F"],
            case_samples=["C"],
            annotations=["NOPE"],
            verbose=False,
        )


def test_fill_from_fasta_sets_the_ancestral_allele(tmp_path: Path, variants: Path) -> None:
    header = tmp_path / "aa.hdr"
    header.write_text(
        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n', encoding="utf-8"
    )
    result = m.bcftools_fill_from_fasta(
        variants,
        tmp_path / "aa.vcf",
        fasta_path=DATA / "genome.fa",
        column="AA",
        header_lines_path=header,
        include='TYPE="snp"',
        verbose=False,
    )
    assert result["records"] == 400
    rows = _plain_rows(tmp_path / "aa.vcf")
    genome: dict[str, str] = {}
    name = ""
    for line in (DATA / "genome.fa").read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            genome[name] = ""
        else:
            genome[name] += line.strip()
    snps = [row for row in rows if len(row[3]) == 1 and len(row[4]) == 1]
    assert snps
    for row in snps[:25]:
        assert f"AA={genome[row[0]][int(row[1]) - 1]}" in row[7]
    with pytest.raises(ConfigurationError):
        m.bcftools_fill_from_fasta(
            variants,
            tmp_path / "bad.vcf",
            fasta_path=DATA / "genome.fa",
            column="AA",
            include="a",
            exclude="b",
            verbose=False,
        )


def test_fill_tags_computes_an_and_ac(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_fill_tags(
        variants, tmp_path / "tags.vcf", tags=["AN", "AC", "NS"], verbose=False
    )
    assert result["records"] == 400
    source_rows = [
        line.split("\t")
        for line in (DATA / "variants.vcf").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    for got, original in list(zip(_plain_rows(tmp_path / "tags.vcf"), source_rows, strict=True))[
        :30
    ]:
        alleles = [
            a
            for call in original[9:]
            for a in call.split(":")[0].replace("|", "/").split("/")
            if a != "."
        ]
        info = dict(item.split("=") for item in got[7].split(";") if "=" in item)
        assert int(info["AN"]) == len(alleles)
        assert int(info["AC"].split(",")[0]) == sum(a == "1" for a in alleles)
    _same_records(
        tmp_path / "tags.vcf", _direct("+fill-tags", str(variants), "--", "-t", "AN,AC,NS")
    )
    with pytest.raises(ConfigurationError):
        m.bcftools_fill_tags(variants, tmp_path / "bad.vcf", tags=[], verbose=False)


def test_exclusive_genotype_sites(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_exclusive_genotype_sites(
        variants, tmp_path / "gts.vcf", samples=["S00", "S01"], verbose=False
    )
    _same_records(tmp_path / "gts.vcf", _direct("+GTsubset", str(variants), "--", "-s", "S00,S01"))
    assert result["records"] == 2
    for row in _plain_rows(tmp_path / "gts.vcf"):
        calls = [c.split(":")[0] for c in row[9:]]
        assert calls[0] == calls[1]


def test_isec_genotypes_sets_differences_to_missing(tmp_path: Path, trio: Path) -> None:
    changed = trio_vcf_text().replace("0/1:60,0,60:15,15", "1/1:60,0,60:15,15", 3)
    other = write_vcf(tmp_path / "other.vcf", changed, compress=True, index=True)
    result = m.bcftools_isec_genotypes(trio, other, tmp_path / "isec.vcf", verbose=False)
    _same_records(tmp_path / "isec.vcf", _direct("+isecGT", str(trio), str(other)))
    missing = sum(
        call.startswith("./.") for row in _plain_rows(tmp_path / "isec.vcf") for call in row[9:]
    )
    assert result["records"] == 40 and missing > 0


def test_prune_sites_by_r2_and_annotation(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_prune_sites(
        variants, tmp_path / "pruned.vcf", max_value=0.6, window="1000", verbose=False
    )
    _same_records(
        tmp_path / "pruned.vcf", _direct("+prune", "-w", "1000", "-m", "r2=0.6", str(variants))
    )
    assert 0 < result["records"] < 400
    annotated = m.bcftools_prune_sites(
        variants,
        tmp_path / "annotated.vcf",
        annotate=["r2", "LD"],
        window="1",
        set_filter=".",
        verbose=False,
    )
    assert annotated["records"] == 400 and "R2=" in _plain_rows(tmp_path / "annotated.vcf")[1][7]
    thinned = m.bcftools_prune_sites(
        variants,
        tmp_path / "thin.vcf",
        sites_per_window=1,
        window="50kb",
        sites_per_window_mode="1st",
        verbose=False,
    )
    assert thinned["records"] < 400
    with pytest.raises(ConfigurationError):
        m.bcftools_prune_sites(variants, tmp_path / "bad.vcf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_prune_sites(
            variants, tmp_path / "bad.vcf", max_value=0.5, statistic="r3", verbose=False
        )


def test_remove_overlaps_removes_duplicates_and_lists_sites(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\nchr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\nchr1\t50\t.\tA\tT\t.\t.\t.\tGT\t1/1\n"
    )
    source = write_vcf(tmp_path / "dup.vcf", text)
    kept = m.bcftools_remove_overlaps(source, tmp_path / "kept.vcf", mark="dup", verbose=False)
    assert kept["records"] == 1 and _plain_rows(tmp_path / "kept.vcf")[0][1] == "50"
    reverse = m.bcftools_remove_overlaps(
        source, tmp_path / "rev.vcf", mark="dup", reverse=True, verbose=False
    )
    assert reverse["records"] == 2
    marked = m.bcftools_remove_overlaps(
        source, tmp_path / "marked.vcf", mark="dup", mark_tag="DUP", verbose=False
    )
    assert (
        marked["records"] == 3
        and sum("DUP" in row[7] for row in _plain_rows(tmp_path / "marked.vcf")) == 2
    )
    listing = m.bcftools_remove_overlaps(
        source, tmp_path / "sites.txt", mark="dup", output_type="t", verbose=False
    )
    assert (tmp_path / "sites.txt").read_text(encoding="utf-8").split() == [
        "chr1",
        "50",
    ] and listing["status"] == "completed"
    with pytest.raises(ConfigurationError):
        m.bcftools_remove_overlaps(source, tmp_path / "bad.vcf", mark="whatever", verbose=False)


def test_scatter_vcf_by_chunk_and_region(tmp_path: Path, variants: Path) -> None:
    chunks = m.bcftools_scatter_vcf(
        variants, tmp_path / "chunks", sites_per_chunk=150, verbose=False
    )
    assert chunks["files"] == ["0.vcf", "1.vcf", "2.vcf"]
    assert sum(len(body(tmp_path / "chunks" / name)) for name in chunks["files"]) == 400
    regions = m.bcftools_scatter_vcf(
        variants,
        tmp_path / "regions",
        regions_to_scatter=["chr1", "chr2"],
        prefix="part-",
        verbose=False,
    )
    assert regions["files"] == ["part-chr1.vcf", "part-chr2.vcf"]
    assert len(body(tmp_path / "regions" / "part-chr1.vcf")) == 194
    with pytest.raises(ConfigurationError):
        m.bcftools_scatter_vcf(variants, tmp_path / "none", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_scatter_vcf(
            variants,
            tmp_path / "two",
            sites_per_chunk=5,
            regions_to_scatter=["chr1"],
            verbose=False,
        )
    assert not (tmp_path / "none").exists()


def test_set_genotypes_and_missing2ref_equivalence(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_set_genotypes(
        variants, tmp_path / "ref.vcf", target_genotypes=".", new_genotype="0p", verbose=False
    )
    assert result["records"] == 400
    assert all("." not in call for row in _plain_rows(tmp_path / "ref.vcf") for call in row[9:])
    _same_records(tmp_path / "ref.vcf", _direct("+missing2ref", str(variants), "--", "-p"))
    _same_records(
        tmp_path / "ref.vcf", _direct("+setGT", str(variants), "--", "-t", ".", "-n", "0p")
    )
    major = m.bcftools_set_genotypes(
        variants, tmp_path / "major.vcf", new_genotype="M", verbose=False
    )
    _same_records(tmp_path / "major.vcf", _direct("+missing2ref", str(variants), "--", "-m"))
    assert major["records"] == 400


def test_set_genotypes_by_expression_and_validation(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_set_genotypes(
        trio,
        tmp_path / "q.vcf",
        target_genotypes="q",
        new_genotype=".",
        genotype_include="FMT/DP<20",
        verbose=False,
    )
    assert result["records"] == 40
    custom = m.bcftools_set_genotypes(
        trio, tmp_path / "c.vcf", target_genotypes="a", new_genotype="c:0/1", verbose=False
    )
    assert custom["records"] == 40
    with pytest.raises(ConfigurationError):
        m.bcftools_set_genotypes(
            trio, tmp_path / "bad.vcf", target_genotypes="q", new_genotype=".", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_set_genotypes(
            trio, tmp_path / "bad.vcf", target_genotypes="zz", new_genotype=".", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_set_genotypes(trio, tmp_path / "bad.vcf", new_genotype="; rm", verbose=False)


def test_split_by_sample_and_groups(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_split_vcf_by_sample(variants, tmp_path / "bysample", verbose=False)
    assert len(result["files"]) == 24 and result["files"][0] == "S00.vcf"
    header = [
        line
        for line in (tmp_path / "bysample" / "S05.vcf").read_text(encoding="utf-8").splitlines()
        if line.startswith("#CHROM")
    ]
    assert header[0].split("\t")[9:] == ["S05"]
    groups = tmp_path / "samples.txt"
    groups.write_text("S00,S01\tA,B\tduo\nS02\n", encoding="utf-8")
    grouped = m.bcftools_split_vcf_by_sample(
        variants, tmp_path / "grouped", samples_file=groups, verbose=False
    )
    assert grouped["files"] == ["S02.vcf", "duo.vcf"] or sorted(grouped["files"]) == [
        "S02.vcf",
        "duo.vcf",
    ]
    with pytest.raises(ConfigurationError):
        m.bcftools_split_vcf_by_sample(
            variants, tmp_path / "both", samples_file=groups, groups_file=groups, verbose=False
        )


CSQ_VCF = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    'Format: Allele|Consequence|IMPACT|SYMBOL|Feature">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t10\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|MODERATE|GENE1|T1,G|synonymous_variant|LOW|GENE1|T2\n"
    "chr1\t20\t.\tC\tT\t.\t.\tCSQ=T|stop_gained|HIGH|GENE2|T3\n"
)


def test_split_vep_columns_text_and_field_list(tmp_path: Path) -> None:
    source = write_vcf(tmp_path / "csq.vcf", CSQ_VCF)
    extracted = m.bcftools_split_vep(
        source,
        tmp_path / "vep.vcf",
        columns="Consequence,IMPACT,SYMBOL",
        select="worst",
        verbose=False,
    )
    assert extracted["records"] == 2
    assert "IMPACT=MODERATE" in _plain_rows(tmp_path / "vep.vcf")[0][7]
    table = m.bcftools_split_vep(
        source,
        tmp_path / "vep.txt",
        select="worst",
        text_format="%CHROM\t%POS\t%Consequence\t%IMPACT\n",
        verbose=False,
    )
    assert table["report"]["rows"] == [
        ["chr1", 10, "missense_variant", "MODERATE"],
        ["chr1", 20, "stop_gained", "HIGH"],
    ]
    fields = m.bcftools_split_vep(source, tmp_path / "fields.txt", list_fields=True, verbose=False)
    assert [f["name"] for f in fields["fields"]] == [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Feature",
    ]
    duplicated = m.bcftools_split_vep(
        source,
        tmp_path / "dup.txt",
        text_format="%CHROM:%POS %SYMBOL %Feature\n",
        duplicate=True,
        verbose=False,
    )
    assert len(duplicated["report"]["rows"]) == 3
    with pytest.raises(ConfigurationError):
        m.bcftools_split_vep(source, tmp_path / "bad.vcf", columns="-bad", verbose=False)


def test_convert_format_tags_pl_to_gl_and_validation(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_convert_format_tags(
        trio, tmp_path / "gl.vcf", conversion="PL-to-GL", verbose=False
    )
    assert result["records"] == 40
    _same_records(
        tmp_path / "gl.vcf",
        _direct("+tag2tag", str(trio), "--", "--PL-to-GL", "-t", "0.1", "-s", "0"),
    )
    first = _plain_rows(tmp_path / "gl.vcf")[0]
    assert "GL" in first[8].split(":")
    gt = m.bcftools_convert_format_tags(
        trio, tmp_path / "gt.vcf", conversion="PL-to-GT", verbose=False
    )
    assert gt["status"] == "completed"
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_format_tags(
            trio, tmp_path / "bad.vcf", conversion="GT-to-PL", verbose=False
        )


def test_convert_format_tags_localized_and_qs(tmp_path: Path) -> None:
    localized = write_vcf(
        tmp_path / "loc.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=LAA,Number=.,Type=Integer,Description="Localized alleles">\n'
        '##FORMAT=<ID=LPL,Number=.,Type=Integer,Description="Localized PL">\n'
        '##FORMAT=<ID=LAD,Number=.,Type=Integer,Description="Localized AD">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tC,G,T\t.\t.\t.\tGT:LAA:LPL:LAD\t1/2:1,2:30,0,40,50,60,70:2,10,12\n",
        compress=True,
        index=True,
    )
    expanded = m.bcftools_convert_format_tags(
        localized, tmp_path / "xx.vcf", conversion="LXX-to-XX", verbose=False
    )
    row = _plain_rows(tmp_path / "xx.vcf")[0]
    assert row[8].split(":")[-2:] == ["AD", "PL"]
    assert (
        row[9].split(":")[-2:] == ["2,10,12,.", "30,0,40,50,60,70,.,.,.,."]
        and expanded["records"] == 1
    )
    qs = write_vcf(
        tmp_path / "qr.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=QR,Number=1,Type=Integer,Description="QR">\n'
        '##FORMAT=<ID=QA,Number=A,Type=Integer,Description="QA">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tC\t.\t.\t.\tGT:QR:QA\t0/1:100:200\n",
    )
    m.bcftools_convert_format_tags(qs, tmp_path / "qs.vcf", conversion="QR-QA-to-QS", verbose=False)
    assert _plain_rows(tmp_path / "qs.vcf")[0][9] == "0/1:100:200:100,200"
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_format_tags(
            localized, tmp_path / "bad.vcf", conversion="XX-to-LXX", verbose=False
        )


def test_variant_distance_matches_neighbours(tmp_path: Path, variants: Path) -> None:
    result = m.bcftools_variant_distance(
        variants, tmp_path / "dist.vcf", direction="both", verbose=False
    )
    assert result["records"] == 400
    rows = _plain_rows(tmp_path / "dist.vcf")
    positions = [(row[0], int(row[1])) for row in rows]
    for index in range(1, 30):
        info = dict(item.split("=") for item in rows[index][7].split(";") if "=" in item)
        before, after = (int(v) for v in info["DIST"].split(","))
        assert before == positions[index][1] - positions[index - 1][1]
        assert after == positions[index + 1][1] - positions[index][1]
    nearest = m.bcftools_variant_distance(
        variants, tmp_path / "near.vcf", tag_name="NEAR", verbose=False
    )
    assert "NEAR=" in _plain_rows(tmp_path / "near.vcf")[0][7] and nearest["records"] == 400
    with pytest.raises(ConfigurationError):
        m.bcftools_variant_distance(
            variants, tmp_path / "bad.vcf", direction="sideways", verbose=False
        )


def test_fixref_statistics_and_conversion(tmp_path: Path, variants: Path) -> None:
    stats = m.bcftools_fixref(
        variants, tmp_path / "fix.txt", reference_fasta=DATA / "genome.fa", verbose=False
    )
    counts = {str(row[0]).strip(): row[1:] for row in stats["sections"]["NS"]}
    assert stats["mode"] == "stats" and counts["total"][0] == 400
    parts = ("ref match", "ref mismatch", "errors", "skipped")
    assert sum(counts[name][0] for name in parts) == counts["total"][0]
    assert sum(row[1] for row in stats["sections"]["ST"]) > 0
    fixed = m.bcftools_fixref(
        variants,
        tmp_path / "fixed.vcf",
        reference_fasta=DATA / "genome.fa",
        mode="flip",
        verbose=False,
    )
    assert fixed["records"] <= 400 and fixed["mode"] == "flip"
    with pytest.raises(ConfigurationError):
        m.bcftools_fixref(
            variants,
            tmp_path / "bad.vcf",
            reference_fasta=DATA / "genome.fa",
            mode="id",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_fixref(
            variants,
            tmp_path / "bad.vcf",
            reference_fasta=DATA / "genome.fa",
            mode="weird",
            verbose=False,
        )


def test_color_chromosomes_finds_shared_segments(tmp_path: Path) -> None:
    source = write_vcf(
        tmp_path / "phased.vcf", trio_vcf_text(phased=True), compress=True, index=True
    )
    result = m.bcftools_color_chromosomes(
        source, tmp_path / "colors", trio=["M", "F", "C"], verbose=False
    )
    assert result["files"] == ["color-chrs.dat"]
    assert result["segments"] and result["segments"][0][:3] == ["chr1", 1, 4799]
    assert {row[0] for row in result["switches"]} == {"M", "F"}
    with pytest.raises(ConfigurationError):
        m.bcftools_color_chromosomes(
            source, tmp_path / "both", trio=["M", "F", "C"], unrelated=["M", "F"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_color_chromosomes(source, tmp_path / "short", trio=["M", "F"], verbose=False)
    assert not (tmp_path / "short").exists()


def test_missing_plugin_program_has_an_install_hint(tmp_path: Path, variants: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        m.bcftools_allele_length_counts(
            variants, tmp_path / "x.txt", bcftools_path=str(tmp_path / "absent"), verbose=False
        )
    existing = tmp_path / "exists.txt"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.bcftools_allele_length_counts(variants, existing, verbose=False)
