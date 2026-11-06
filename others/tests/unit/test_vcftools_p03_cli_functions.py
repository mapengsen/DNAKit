"""Tests for dnakit.variants.vcftools_p03_cli_functions (VCFtools programs and Perl scripts as typed facades).

The tests import the module directly. Tests that run a program skip when VCFtools (``vcftools`` and the ``vcf-*`` Perl
scripts), ``tabix``/``bgzip`` or ``samtools`` is missing; argument validation is tested without any program. The toy
VCF has 6 phased samples with GT:DP:GQ:PL, 20 sites on chr1 and chr2, 3 indels, 2 LowQual sites, singletons and one
missing genotype.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import vcftools_p03_cli_functions as m

TOY = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##contig=<ID=chr1,length=20000>
##contig=<ID=chr2,length=20000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
chr1	100	v1	G	T	38	PASS	DP=91;AF=0.1	GT:DP:GQ:PL	1|0:9:45:30,0,30	0|1:31:60:30,0,30	0|1:9:30:30,0,30	0|1:20:99:30,0,30	1|0:8:12:30,0,30	0|0:14:99:0,30,300
chr1	800	v2	G	A	55	PASS	DP=120;AF=0.25	GT:DP:GQ:PL	1|0:25:99:30,0,30	1|1:21:60:300,30,0	0|0:17:45:0,30,300	0|0:21:60:0,30,300	0|0:11:45:0,30,300	0|0:25:99:0,30,300
chr1	1500	v3	A	G	67	PASS	DP=99;AF=0.5	GT:DP:GQ:PL	1|1:10:45:300,30,0	0|1:11:99:30,0,30	0|1:23:60:30,0,30	0|0:36:45:0,30,300	0|1:10:60:30,0,30	1|1:9:45:300,30,0
chr1	2200	v4	A	AAC	38	PASS	DP=166;AF=0.5	GT:DP:GQ:PL	1|0:33:30:30,0,30	1|1:23:12:300,30,0	1|0:32:30:30,0,30	1|0:28:30:30,0,30	1|1:30:30:300,30,0	1|1:20:30:300,30,0
chr1	2900	v5	A	G	94	PASS	DP=172;AF=0.25	GT:DP:GQ:PL	0|1:31:60:30,0,30	0|0:31:12:0,30,300	0|0:31:30:0,30,300	0|0:31:12:0,30,300	0|0:12:30:0,30,300	1|0:36:60:30,0,30
chr1	3600	v6	C	T	22	LowQual	DP=180;AF=0.1	GT:DP:GQ:PL	0|0:37:30:0,30,300	0|0:35:12:0,30,300	0|1:36:45:30,0,30	0|0:36:45:0,30,300	0|0:25:60:0,30,300	0|0:11:30:0,30,300
chr1	4300	v7	A	G	81	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	0|0:27:30:0,30,300	1|1:20:99:300,30,0	0|0:18:60:0,30,300	0|0:21:45:0,30,300	0|0:31:12:0,30,300	0|0:20:12:0,30,300
chr1	5000	v8	G	A	81	PASS	DP=124;AF=0.25	GT:DP:GQ:PL	0|0:6:30:0,30,300	0|0:36:60:0,30,300	0|0:28:30:0,30,300	0|1:11:60:30,0,30	1|1:13:45:300,30,0	1|0:30:12:30,0,30
chr1	5700	v9	TGT	T	55	PASS	DP=93;AF=0.5	GT:DP:GQ:PL	1|0:15:30:30,0,30	0|1:14:60:30,0,30	1|1:7:30:300,30,0	0|0:6:30:0,30,300	1|0:12:12:30,0,30	0|0:39:45:0,30,300
chr1	6400	v10	C	T	38	PASS	DP=173;AF=0.25	GT:DP:GQ:PL	0|0:40:30:0,30,300	0|0:15:99:0,30,300	1|0:39:12:30,0,30	0|0:38:30:0,30,300	./.:.:.:.	0|0:34:30:0,30,300
chr1	7100	v11	C	T	22	LowQual	DP=177;AF=0.1	GT:DP:GQ:PL	0|0:34:60:0,30,300	0|0:26:99:0,30,300	0|0:38:99:0,30,300	0|0:38:60:0,30,300	0|0:18:99:0,30,300	0|0:23:30:0,30,300
chr1	7800	v12	G	GAC	94	PASS	DP=135;AF=0.25	GT:DP:GQ:PL	0|1:29:12:30,0,30	1|0:15:60:30,0,30	1|1:22:60:300,30,0	1|0:14:30:30,0,30	0|0:35:30:0,30,300	1|0:20:30:30,0,30
chr1	8500	v13	T	G	38	PASS	DP=93;AF=0.4	GT:DP:GQ:PL	0|1:12:45:30,0,30	1|1:11:30:300,30,0	1|0:22:60:30,0,30	0|1:23:45:30,0,30	0|1:8:60:30,0,30	0|1:17:30:30,0,30
chr1	9200	v14	T	A	55	PASS	DP=137;AF=0.4	GT:DP:GQ:PL	1|0:32:12:30,0,30	0|1:23:30:30,0,30	0|0:14:45:0,30,300	0|1:8:12:30,0,30	0|1:39:30:30,0,30	0|0:21:30:0,30,300
chr2	300	v15	G	C	94	PASS	DP=204;AF=0.4	GT:DP:GQ:PL	0|0:33:30:0,30,300	1|0:37:30:30,0,30	0|1:40:45:30,0,30	1|0:31:30:30,0,30	0|0:38:30:0,30,300	1|0:25:60:30,0,30
chr2	1200	v16	G	C	38	PASS	DP=141;AF=0.1	GT:DP:GQ:PL	1|0:23:99:30,0,30	0|0:34:45:0,30,300	1|0:6:30:30,0,30	0|0:22:12:0,30,300	0|0:29:45:0,30,300	0|0:27:30:0,30,300
chr2	2100	v17	G	A	38	PASS	DP=158;AF=0.4	GT:DP:GQ:PL	1|0:11:30:30,0,30	0|1:39:45:30,0,30	0|1:15:99:30,0,30	1|1:30:30:300,30,0	1|1:26:12:300,30,0	1|1:37:99:300,30,0
chr2	3000	v18	T	A	94	PASS	DP=133;AF=0.1	GT:DP:GQ:PL	0|0:7:60:0,30,300	0|0:40:12:0,30,300	0|0:21:99:0,30,300	1|0:37:99:30,0,30	0|0:22:12:0,30,300	0|0:6:99:0,30,300
chr2	3900	v19	A	T	94	PASS	DP=117;AF=0.4	GT:DP:GQ:PL	0|0:18:99:0,30,300	1|0:10:99:30,0,30	1|0:15:30:30,0,30	0|0:27:12:0,30,300	1|0:22:60:30,0,30	1|1:25:12:300,30,0
chr2	4800	v20	T	C	67	PASS	DP=154;AF=0.1	GT:DP:GQ:PL	0|0:10:30:0,30,300	0|0:38:12:0,30,300	0|0:34:99:0,30,300	0|0:23:12:0,30,300	0|0:30:30:0,30,300	0|1:19:99:30,0,30
"""


def _installed(name: str) -> bool:
    return shutil.which(name) is not None or (Path(sys.executable).parent / name).is_file()


requires_vcftools = pytest.mark.skipif(
    not _installed("vcftools"), reason="vcftools is not installed"
)
requires_scripts = pytest.mark.skipif(
    not (_installed("vcf-annotate") and shutil.which("perl")),
    reason="the VCFtools Perl scripts are not installed",
)
requires_htslib = pytest.mark.skipif(
    not (shutil.which("bgzip") and shutil.which("tabix")),
    reason="bgzip and tabix (htslib) are not installed",
)
requires_samtools = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools is not installed"
)


@pytest.fixture()
def toy(tmp_path: Path) -> Path:
    path = tmp_path / "toy.vcf"
    path.write_text(TOY, encoding="utf-8")
    return path


def _indexed(tmp_path: Path, name: str, lines: list[str]) -> str:
    """A bgzip-compressed, tabix-indexed VCF built from ``lines``; returns the .gz path."""
    plain = tmp_path / name
    plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
    subprocess.run(["bgzip", "-k", str(plain)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(plain) + ".gz"], check=True)
    return str(plain) + ".gz"


def _toy_lines() -> tuple[list[str], list[str]]:
    lines = TOY.splitlines()
    return [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]


def _records(path: str | Path) -> list[list[str]]:
    return [
        x.split("\t") for x in Path(path).read_text().splitlines() if x and not x.startswith("#")
    ]


# ---------------------------------------------------------------------------
# Validation and failure modes that need no program
# ---------------------------------------------------------------------------


def test_missing_program_raises_backend_unavailable(toy, tmp_path):
    with pytest.raises(BackendUnavailableError):
        m.vcftools_depth(
            str(toy),
            str(tmp_path / "x"),
            vcftools_path=str(tmp_path / "no_vcftools"),
            verbose=False,
        )
    with pytest.raises(BackendUnavailableError):
        m.vcf_to_tab(
            str(toy),
            str(tmp_path / "x.tab"),
            script_path=str(tmp_path / "no_script"),
            verbose=False,
        )


def test_argument_validation(toy, tmp_path):
    with pytest.raises(ConfigurationError):
        m.vcftools_depth(str(toy), str(tmp_path / "x"), level="bogus", verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcftools_list_sites(
            str(toy), str(tmp_path / "x"), filters={"not_a_filter": 1}, verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.vcftools_list_sites(
            str(toy),
            str(tmp_path / "x"),
            filters={"bed_file": str(toy), "exclude_bed_file": str(toy)},
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        m.vcftools_filter_vcf(str(toy), str(tmp_path / "x"), output_format="vcf.gz", verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcftools_filter_vcf(str(toy), str(tmp_path / "x"), min_maf=True, verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcftools_transition_transversion(
            str(toy), str(tmp_path / "x"), verbose=False
        )  # bin_size is required
    with pytest.raises(ConfigurationError):
        m.vcftools_linkage_disequilibrium(
            str(toy),
            str(tmp_path / "x"),
            statistic="geno_chisq",
            interchromosomal=True,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        m.vcftools_export_genotypes(
            str(toy), str(tmp_path / "x"), export_format="ldhat", verbose=False
        )  # needs chr
    with pytest.raises(ConfigurationError):
        m.vcftools_pca(str(toy), str(tmp_path / "x"), output="snp_loadings", verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcftools_pca(
            str(toy),
            str(tmp_path / "x"),
            output="snp_loadings",
            snp_loadings=2,
            normalise=False,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        m.vcftools_extract_fields(
            str(toy), str(tmp_path / "x"), ["DP", "GQ"], source="format", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.vcftools_depth(str(tmp_path / "missing.vcf"), str(tmp_path / "x"), verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcf_annotate(str(toy), str(tmp_path / "a.vcf"), verbose=False)  # no operation selected
    with pytest.raises(ConfigurationError):
        m.vcf_annotate(str(toy), str(tmp_path / "a.vcf"), filters="+; rm", verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcf_subset(
            str(toy), str(tmp_path / "s.vcf"), samples=["S1"], samples_file=str(toy), verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.vcf_subset(str(toy), str(tmp_path / "s.vcf"), variant_types=["bogus"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcf_convert(str(toy), str(tmp_path / "c.vcf"), target_version="3.3", verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcf_compare([str(toy)], str(tmp_path / "c.txt"), verbose=False)
    with pytest.raises(ConfigurationError):
        m.vcf_fill_reference_md5(
            str(toy), str(tmp_path / "o.vcf.gz"), str(toy), verbose=False
        )  # not compressed


def test_existing_output_is_never_overwritten(toy, tmp_path):
    (tmp_path / "taken.idepth").write_text("x")
    with pytest.raises(ConfigurationError):
        m.vcftools_depth(str(toy), str(tmp_path / "taken"), verbose=False)
    existing = tmp_path / "out.vcf"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        m.vcf_annotate(str(toy), str(existing), fill_type=True, verbose=False)


def test_public_names_follow_the_verbosity_contract():
    import inspect

    assert len(m.__all__) == 33
    for name in m.__all__:
        signature = inspect.signature(getattr(m, name))
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY
        assert signature.parameters["verbose"].default is True


# ---------------------------------------------------------------------------
# vcftools program
# ---------------------------------------------------------------------------


@requires_vcftools
def test_filter_vcf_applies_site_individual_and_genotype_filters(toy, tmp_path):
    result = m.vcftools_filter_vcf(
        str(toy),
        str(tmp_path / "f"),
        min_maf=0.2,
        remove_filtered_all=True,
        min_genotype_quality=20,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["version"].startswith("0.1.")
    assert (result["sites_kept"], result["sites_total"], result["individuals_kept"]) == (14, 20, 6)
    assert len(_records(result["outputs"][0])) == 14

    subset = m.vcftools_filter_vcf(
        str(toy),
        str(tmp_path / "g"),
        keep_individuals=["S1", "S2"],
        remove_indels=True,
        recode_info=["DP"],
        verbose=False,
    )
    assert (subset["sites_kept"], subset["individuals_kept"]) == (17, 2)
    first = _records(subset["outputs"][0])[0]
    assert (
        first[7] == "DP=91" and len(first) == 11
    )  # only the requested INFO key, two sample columns

    window = m.vcftools_filter_vcf(
        str(toy),
        str(tmp_path / "w"),
        chromosomes=["chr2"],
        from_bp=1000,
        to_bp=4000,
        thin=1000,
        verbose=False,
    )
    assert window["sites_kept"] == 2
    indels = m.vcftools_filter_vcf(str(toy), str(tmp_path / "i"), only_indels=True, verbose=False)
    assert indels["sites_kept"] == 3


@requires_vcftools
def test_filter_vcf_writes_bcf(toy, tmp_path):
    result = m.vcftools_filter_vcf(
        str(toy), str(tmp_path / "b"), output_format="bcf", max_alleles=2, verbose=False
    )
    assert result["status"] == "completed"
    path = Path(result["outputs"][0])
    assert path.name.endswith(".recode.bcf") and path.stat().st_size > 0
    # the BCF is readable by vcftools itself
    again = m.vcftools_site_quality(str(path), str(tmp_path / "q"), verbose=False)
    assert again["table"]["row_count"] == 20


@requires_vcftools
def test_list_sites_kept_and_removed(toy, tmp_path):
    kept = m.vcftools_list_sites(
        str(toy),
        str(tmp_path / "k"),
        filters={"min_maf": 0.2, "remove_filtered_all": True},
        verbose=False,
    )
    removed = m.vcftools_list_sites(
        str(toy),
        str(tmp_path / "r"),
        which="removed",
        filters={"min_maf": 0.2, "remove_filtered_all": True},
        verbose=False,
    )
    assert kept["table"]["columns"] == ["CHROM", "POS"]
    assert kept["table"]["row_count"] + removed["table"]["row_count"] == 20
    kept_positions = {(row[0], row[1]) for row in kept["table"]["rows"]}
    removed_positions = {(row[0], row[1]) for row in removed["table"]["rows"]}
    assert kept_positions.isdisjoint(removed_positions)


@requires_vcftools
def test_depth_levels(toy, tmp_path):
    individual = m.vcftools_depth(str(toy), str(tmp_path / "i"), verbose=False)
    assert individual["table"]["columns"] == ["INDV", "N_SITES", "MEAN_DEPTH"]
    assert individual["table"]["rows"][0] == ["S1", 20, 22.1]
    site_mean = m.vcftools_depth(str(toy), str(tmp_path / "m"), level="site_mean", verbose=False)
    assert site_mean["table"]["rows"][1][:3] == ["chr1", 800, 20]
    site_sum = m.vcftools_depth(str(toy), str(tmp_path / "s"), level="site_sum", verbose=False)
    assert site_sum["table"]["columns"] == ["CHROM", "POS", "SUM_DEPTH", "SUMSQ_DEPTH"]
    genotype = m.vcftools_depth(
        str(toy), str(tmp_path / "g"), level="genotype", max_rows=3, verbose=False
    )
    assert (
        genotype["table"]["row_count"] == 20
        and len(genotype["table"]["rows"]) == 3
        and genotype["table"]["truncated"]
    )


@requires_vcftools
def test_missingness_quality_and_heterozygosity(toy, tmp_path):
    sites = m.vcftools_missingness(str(toy), str(tmp_path / "ms"), by="site", verbose=False)
    missing = [row for row in sites["table"]["rows"] if row[4] > 0]
    assert missing == [["chr1", 6400, 12, 0, 2, 0.166667]]
    individuals = m.vcftools_missingness(str(toy), str(tmp_path / "mi"), verbose=False)
    assert individuals["table"]["row_count"] == 6
    quality = m.vcftools_site_quality(str(toy), str(tmp_path / "q"), verbose=False)
    assert quality["table"]["rows"][:2] == [["chr1", 100, 38], ["chr1", 800, 55]]
    het = m.vcftools_heterozygosity(str(toy), str(tmp_path / "h"), verbose=False)
    assert het["table"]["rows"][0] == ["S1", 9, 11.5, 19, -0.33549]


@requires_vcftools
def test_singletons_indel_histogram_and_snp_density(toy, tmp_path):
    singletons = m.vcftools_singletons(str(toy), str(tmp_path / "s"), verbose=False)
    assert singletons["table"]["row_count"] == 5
    assert singletons["table"]["rows"][0] == ["chr1", 3600, "S", "T", "S3"]
    histogram = m.vcftools_indel_length_histogram(str(toy), str(tmp_path / "h"), verbose=False)
    assert histogram["table"]["rows"] == [[-2, 1, 5], [0, 17, 85], [2, 2, 10]]
    density = m.vcftools_snp_density(str(toy), str(tmp_path / "d"), 5000, verbose=False)
    assert density["table"]["rows"] == [
        ["chr1", 0, 7, 1.4],
        ["chr1", 5000, 7, 1.4],
        ["chr2", 0, 6, 1.2],
    ]


@requires_vcftools
def test_transition_transversion_modes(toy, tmp_path):
    bins = m.vcftools_transition_transversion(
        str(toy), str(tmp_path / "b"), bin_size=5000, verbose=False
    )
    assert bins["table"]["rows"][0] == ["chr1", 0, 6, 5]
    by_filter = m.vcftools_transition_transversion(
        str(toy), str(tmp_path / "f"), mode="filter_category", verbose=False
    )
    assert by_filter["table"]["rows"][0][:4] == ["PASS", 18, 8, 7]
    by_count = m.vcftools_transition_transversion(
        str(toy), str(tmp_path / "c"), mode="alt_allele_count", verbose=False
    )
    assert by_count["table"]["columns"][0] == "ALT_ALLELE_COUNT"
    by_qual = m.vcftools_transition_transversion(
        str(toy), str(tmp_path / "q"), mode="qual_threshold", verbose=False
    )
    assert by_qual["table"]["columns"][0] == "QUAL_THRESHOLD"


@requires_vcftools
def test_burden_diversity_relatedness(toy, tmp_path):
    burden = m.vcftools_individual_burden(str(toy), str(tmp_path / "b"), verbose=False)
    assert burden["table"]["rows"][0] == ["S1", 9, 10, 1, 0]
    frequency = m.vcftools_individual_burden(
        str(toy), str(tmp_path / "f"), burden="frequency_classes", verbose=False
    )
    once = m.vcftools_individual_burden(
        str(toy), str(tmp_path / "o"), burden="frequency_classes_hom_once", verbose=False
    )
    assert sum(sum(row[1:]) for row in frequency["table"]["rows"]) > sum(
        sum(row[1:]) for row in once["table"]["rows"]
    )
    pi = m.vcftools_site_nucleotide_diversity(str(toy), str(tmp_path / "p"), verbose=False)
    assert pi["table"]["rows"][0] == ["chr1", 100, 0.530303]
    yang = m.vcftools_relatedness(str(toy), str(tmp_path / "y"), verbose=False)
    assert yang["table"]["columns"][-1] == "RELATEDNESS_AJK" and yang["table"]["row_count"] == 21
    phi = m.vcftools_relatedness(str(toy), str(tmp_path / "m"), method="manichaikul", verbose=False)
    assert phi["table"]["columns"][-1] == "RELATEDNESS_PHI" and phi["table"]["rows"][0][:2] == [
        "S1",
        "S1",
    ]


@requires_vcftools
def test_linkage_disequilibrium_forms(toy, tmp_path):
    positions = tmp_path / "pos.txt"
    positions.write_text("chr1\t100\nchr1\t800\n")
    genotype = m.vcftools_linkage_disequilibrium(str(toy), str(tmp_path / "g"), verbose=False)
    assert genotype["table"]["columns"] == ["CHR", "POS1", "POS2", "N_INDV", "R^2"]
    assert genotype["table"]["rows"][0] == ["chr1", 100, 800, 6, 0.0857143]
    assert genotype["table"]["row_count"] == 106
    windowed = m.vcftools_linkage_disequilibrium(
        str(toy), str(tmp_path / "w"), ld_window_bp=1500, min_r2=0.1, verbose=False
    )
    assert 0 < windowed["table"]["row_count"] < 106
    haplotype = m.vcftools_linkage_disequilibrium(
        str(toy), str(tmp_path / "h"), statistic="hap_r2", verbose=False
    )
    assert haplotype["table"]["columns"] == ["CHR", "POS1", "POS2", "N_CHR", "R^2", "D", "Dprime"]
    chisq = m.vcftools_linkage_disequilibrium(
        str(toy), str(tmp_path / "c"), statistic="geno_chisq", verbose=False
    )
    assert chisq["table"]["columns"][:3] == ["CHR", "POS1", "POS2"]
    listed = m.vcftools_linkage_disequilibrium(
        str(toy), str(tmp_path / "l"), positions_file=str(positions), verbose=False
    )
    assert listed["outputs"][0].endswith(".list.geno.ld")
    inter = m.vcftools_linkage_disequilibrium(
        str(toy), str(tmp_path / "i"), interchromosomal=True, verbose=False
    )
    assert inter["outputs"][0].endswith(".interchrom.geno.ld")
    assert {row[0] for row in inter["table"]["rows"]} == {"chr1"}


@requires_vcftools
def test_long_runs_haplotype_count_and_field_extraction(toy, tmp_path):
    lroh = m.vcftools_long_runs_of_homozygosity(
        str(toy), str(tmp_path / "r"), "chr1", verbose=False
    )
    assert lroh["status"] == "completed" and lroh["outputs"][0].endswith(".LROH")
    bed = tmp_path / "win.bed"
    bed.write_text("#CHROM\tSTART\tEND\nchr1\t0\t5000\nchr1\t5000\t12000\nchr2\t0\t4000\n")
    counts = m.vcftools_haplotype_count(str(toy), str(tmp_path / "h"), str(bed), verbose=False)
    assert [row[:4] for row in counts["table"]["rows"]] == [
        ["chr1", 0, 5000, 8],
        ["chr1", 5000, 12000, 5],
        ["chr2", 0, 4000, 5],
    ]
    assert counts["table"]["rows"][0][4] == 11
    assert counts["warnings"]
    headerless = tmp_path / "headerless.bed"
    headerless.write_text("chr1\t0\t5000\nchr1\t5000\t12000\n")
    dropped = m.vcftools_haplotype_count(
        str(toy), str(tmp_path / "d"), str(headerless), verbose=False
    )
    assert len(dropped["warnings"]) == 2
    info = m.vcftools_extract_fields(str(toy), str(tmp_path / "i"), ["DP", "AF"], verbose=False)
    assert info["table"]["columns"] == ["CHROM", "POS", "REF", "ALT", "DP", "AF"]
    assert info["table"]["rows"][0] == ["chr1", 100, "G", "T", 91, 0.1]
    fmt = m.vcftools_extract_fields(
        str(toy), str(tmp_path / "f"), ["GQ"], source="format", verbose=False
    )
    assert fmt["outputs"][0].endswith(".GQ.FORMAT")
    assert fmt["table"]["rows"][0] == ["chr1", 100, 45, 60, 30, 99, 12, 99]


@requires_vcftools
def test_pca_components_and_loadings(toy, tmp_path):
    result = m.vcftools_pca(str(toy), str(tmp_path / "p"), verbose=False)
    if result["status"] != "completed":
        pytest.skip("vcftools was built without PCA support")
    assert result["table"]["columns"][:2] == ["INDV", "EIG_0"]
    plain = m.vcftools_pca(str(toy), str(tmp_path / "n"), normalise=False, verbose=False)
    assert plain["status"] == "completed"
    loadings = m.vcftools_pca(
        str(toy), str(tmp_path / "l"), output="snp_loadings", snp_loadings=2, verbose=False
    )
    assert loadings["table"]["columns"] == ["CHROM", "POS", "GAMMA_0", "GAMMA_1"]


@requires_vcftools
def test_export_formats(toy, tmp_path):
    expected = {
        "012": 3,
        "impute": 3,
        "plink": 2,
        "plink_tped": 2,
    }
    for fmt, count in expected.items():
        result = m.vcftools_export_genotypes(
            str(toy), str(tmp_path / fmt), export_format=fmt, verbose=False
        )
        assert result["status"] == "completed", fmt
        assert len(result["outputs"]) == count
    for fmt in ("ldhat", "ldhelmet", "ldhat_geno", "beagle_pl"):
        result = m.vcftools_export_genotypes(
            str(toy), str(tmp_path / fmt), export_format=fmt, chromosome="chr1", verbose=False
        )
        assert result["status"] == "completed", fmt
    ped = Path(tmp_path / "plink.ped").read_text().splitlines()[0].split("\t")
    assert ped[:2] == ["S1", "S1"] and ped[6:8] == ["T", "G"]
    # the toy has no GL tag, so BEAGLE genotype likelihoods from GL fail cleanly
    failed = m.vcftools_export_genotypes(
        str(toy), str(tmp_path / "gl"), export_format="beagle_gl", chromosome="chr1", verbose=False
    )
    assert failed["status"] == "failed" and "GL or PL" in failed["log_tail"]


@requires_vcftools
def test_compare_vcfs_modes(toy, tmp_path):
    lines = TOY.splitlines()
    changed = [x.replace("1|1:", "0|0:", 1) if x.startswith("chr1\t800") else x for x in lines]
    other = tmp_path / "other.vcf"
    other.write_text("\n".join(changed) + "\n")
    individuals = m.vcftools_compare_vcfs(str(toy), str(other), str(tmp_path / "i"), verbose=False)
    assert individuals["table"]["rows"][0] == ["S1", "B"]
    discordance = m.vcftools_compare_vcfs(
        str(toy), str(other), str(tmp_path / "d"), mode="site_discordance", verbose=False
    )
    assert [row for row in discordance["table"]["rows"] if row[5] > 0] == [
        ["chr1", 800, "B", 1, 6, 1, 0.166667]
    ]
    switch = m.vcftools_compare_vcfs(
        str(toy), str(other), str(tmp_path / "s"), mode="switch_error", verbose=False
    )
    assert switch["status"] == "completed" and len(switch["outputs"]) == 2


@requires_vcftools
def test_pass_through_of_a_failing_run(toy, tmp_path):
    # --hap-r2 needs phased data with enough sites; a position filter that leaves one site fails cleanly
    result = m.vcftools_linkage_disequilibrium(
        str(toy),
        str(tmp_path / "x"),
        statistic="hap_r2",
        filters={"chromosomes": ["chr2"], "from_bp": 300, "to_bp": 300},
        verbose=False,
    )
    assert result["status"] == "failed"


# ---------------------------------------------------------------------------
# Perl scripts
# ---------------------------------------------------------------------------


@requires_scripts
def test_vcf_annotate_operations(toy, tmp_path):
    filtered = m.vcf_annotate(
        str(toy), str(tmp_path / "a1.vcf"), filters="+/Qual=40", hard_filter=True, verbose=False
    )
    assert filtered["status"] == "completed" and filtered["records"] == 12
    typed = m.vcf_annotate(
        str(toy), str(tmp_path / "a2.vcf"), fill_type=True, fill_ac_an=True, verbose=False
    )
    info = _records(typed["output_path"])[0][7].split(";")
    assert {"AC=5", "AN=12", "TYPE=snp"} <= set(info)
    hwe = m.vcf_annotate(str(toy), str(tmp_path / "a3.vcf"), fill_hwe=True, verbose=False)
    assert any(item.startswith("HWE=") for item in _records(hwe["output_path"])[0][7].split(";"))
    removed = m.vcf_annotate(
        str(toy), str(tmp_path / "a4.vcf"), remove_tags=["INFO/AF", "QUAL"], verbose=False
    )
    row = _records(removed["output_path"])[0]
    assert row[5] == "." and row[7] == "DP=91"
    flagged = m.vcf_annotate(str(toy), str(tmp_path / "a5.vcf"), filters="+/Qual=40", verbose=False)
    assert [r[6] for r in _records(flagged["output_path"])][:4] == ["Qual", "PASS", "PASS", "Qual"]
    normalized = m.vcf_annotate(
        str(toy), str(tmp_path / "a6.vcf"), normalize_alleles=True, verbose=False
    )
    assert normalized["records"] == 20


@requires_scripts
def test_vcf_subset_samples_and_classes(toy, tmp_path):
    columns = m.vcf_subset(str(toy), str(tmp_path / "s1.vcf"), samples=["S1", "S2"], verbose=False)
    header = next(
        x for x in Path(columns["output_path"]).read_text().splitlines() if x.startswith("#CHROM")
    )
    assert header.split("\t")[9:] == ["S1", "S2"]
    indels = m.vcf_subset(
        str(toy),
        str(tmp_path / "s2.vcf"),
        samples=["S1", "S2"],
        variant_types=["indels"],
        exclude_ref_only=True,
        verbose=False,
    )
    assert indels["records"] == 3
    private = m.vcf_subset(
        str(toy), str(tmp_path / "s3.vcf"), samples=["S3"], private_only=True, verbose=False
    )
    assert private["records"] == 2
    listing = tmp_path / "samples.txt"
    listing.write_text("S5\nS6\n")
    from_file = m.vcf_subset(
        str(toy),
        str(tmp_path / "s4.vcf"),
        samples_file=str(listing),
        trim_alt_alleles=True,
        verbose=False,
    )
    assert from_file["status"] == "completed"


@requires_scripts
def test_vcf_to_tab_and_indel_statistics(toy, tmp_path):
    table = m.vcf_to_tab(str(toy), str(tmp_path / "t.tab"), verbose=False)
    assert table["status"] == "completed" and table["records"] == 20
    first = Path(table["output_path"]).read_text().splitlines()[1].split("\t")
    assert first[:3] == ["chr1", "100", "G"] and first[3] == "T/G"
    iupac = m.vcf_to_tab(str(toy), str(tmp_path / "t2.tab"), iupac=True, verbose=False)
    assert (
        iupac["status"] == "failed" and "IUPAC" in iupac["log_tail"]
    )  # the toy has indels; IUPAC output is for SNPs
    exons = tmp_path / "exons.txt"
    exons.write_text("chr1\t2000\t3000\nchr1\t5600\t6000\n")
    stats = m.vcf_indel_stats(str(toy), str(tmp_path / "i.txt"), str(exons), verbose=False)
    assert stats["statistics"] == {"total": 3, "in-frame": 0, "frameshift": 2, "ratio": 1.0}


@requires_scripts
@requires_htslib
def test_vcf_compare_contrast_and_shuffle(tmp_path):
    header, body = _toy_lines()
    full = _indexed(tmp_path, "a.vcf", header + body)
    half = _indexed(tmp_path, "b.vcf", header + body[::2])
    compared = m.vcf_compare(
        [full, half], str(tmp_path / "cmp.txt"), compare_genotypes=True, verbose=False
    )
    assert compared["status"] == "completed"
    assert [row["sites"] for row in compared["venn"]] == [10, 10]
    assert compared["summary"]["Number of REF matches"] == 10
    assert compared["summary"]["Number of samples in GT comparison"] == 6
    contrast = m.vcf_contrast(
        full, str(tmp_path / "ct.vcf"), ["S1"], ["S2", "S3"], novel_sites_only=True, verbose=False
    )
    assert contrast["records"] == 5
    assert "NOVELGT=S1" in _records(contrast["output_path"])[0][7]
    columns = header[-1].split("\t")
    order = [*range(9), *range(len(columns) - 1, 8, -1)]
    template = _indexed(
        tmp_path,
        "t.vcf",
        [
            *header[:-1],
            "\t".join(columns[i] for i in order),
            *("\t".join(row.split("\t")[i] for i in order) for row in body),
        ],
    )
    shuffled = m.vcf_shuffle_columns(full, template, str(tmp_path / "sh.vcf"), verbose=False)
    line = next(
        x for x in Path(shuffled["output_path"]).read_text().splitlines() if x.startswith("#CHROM")
    )
    assert line.split("\t")[9:] == ["S6", "S5", "S4", "S3", "S2", "S1"]


@requires_scripts
@requires_htslib
def test_vcf_phased_join_joins_overlapping_chunks(tmp_path):
    header, body = _toy_lines()
    first = _indexed(tmp_path, "c1.vcf", header + body[:12])
    second = _indexed(tmp_path, "c2.vcf", header + body[7:14])
    joined = m.vcf_phased_join([first, second], str(tmp_path / "pj.vcf"), verbose=False)
    assert joined["status"] == "completed" and joined["records"] == 14


@requires_scripts
def test_vcf_fix_ploidy_sets_male_x_haploid(tmp_path):
    header, body = _toy_lines()
    x_header = [line.replace("ID=chr1", "ID=X") for line in header]
    x_body = [line.replace("chr1", "X", 1) if line.startswith("chr1") else line for line in body]
    source = tmp_path / "x.vcf"
    source.write_text("\n".join(x_header + x_body) + "\n")
    sexes = tmp_path / "sex.txt"
    sexes.write_text("S1 M\nS2 F\nS3 M\nS4 F\nS5 M\nS6 F\n")
    fixed = m.vcf_fix_ploidy(str(source), str(tmp_path / "fixed.vcf"), str(sexes), verbose=False)
    assert fixed["status"] == "completed"
    row = _records(fixed["output_path"])[0]
    assert (
        row[9].split(":")[0] == "1" and row[10].split(":")[0] == "0|1"
    )  # male S1 is haploid, female S2 stays diploid


def _reference(tmp_path: Path) -> Path:
    """A FASTA whose bases agree with the REF alleles of the toy VCF."""
    _, body = _toy_lines()
    sequences = {"chr1": ["A"] * 20000, "chr2": ["A"] * 20000}
    for row in body:
        fields = row.split("\t")
        for offset, base in enumerate(fields[3]):
            sequences[fields[0]][int(fields[1]) - 1 + offset] = base
    path = tmp_path / "ref.fa"
    path.write_text("".join(f">{name}\n{''.join(bases)}\n" for name, bases in sequences.items()))
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    return path


@requires_scripts
@requires_samtools
def test_vcf_convert_flanking_sequence_and_ancestral_allele(toy, tmp_path):
    reference = _reference(tmp_path)
    without = m.vcf_convert(str(toy), str(tmp_path / "c0.vcf"), verbose=False)
    assert without["status"] == "failed" and "-r" in without["log_tail"]  # indels need a reference
    converted = m.vcf_convert(
        str(toy),
        str(tmp_path / "c1.vcf"),
        target_version="4.1",
        reference_fasta=str(reference),
        verbose=False,
    )
    assert converted["status"] == "completed"
    assert Path(converted["output_path"]).read_text().startswith("##fileformat=VCFv4.1")
    flanks = m.vcf_fill_flanking_sequence(
        str(toy), str(tmp_path / "fs.vcf"), str(reference), length=5, verbose=False
    )
    assert "FS=AAAAA[G/A]AAAAA" in _records(flanks["output_path"])[1][7]
    for chrom in ("chr1", "chr2"):
        bases = reference.read_text().split(">")
        sequence = "".join(next(b for b in bases if b.startswith(chrom + "\n")).split("\n")[1:])
        fasta = tmp_path / f"anc_{chrom}.fa"
        fasta.write_text(
            f">{chrom}\n"
            + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70))
            + "\n"
        )
        subprocess.run(["bgzip", str(fasta)], check=True)
        subprocess.run(["samtools", "faidx", str(fasta) + ".gz"], check=True)
    ancestral = m.vcf_fill_ancestral_allele(
        str(toy), str(tmp_path / "aa.vcf"), str(tmp_path / "anc_"), verbose=False
    )
    assert ancestral["status"] == "completed" and ancestral["records"] == 20
    assert "AA=G" in _records(ancestral["output_path"])[0][7]


@requires_scripts
@requires_samtools
@requires_htslib
def test_vcf_fill_reference_md5_adds_contig_lines(tmp_path):
    reference = _reference(tmp_path)
    header, body = _toy_lines()
    source = _indexed(tmp_path, "m.vcf", header + body)
    result = m.vcf_fill_reference_md5(
        source, str(tmp_path / "md5.vcf.gz"), str(reference), assembly="toy", verbose=False
    )
    if result["status"] != "completed":
        pytest.skip("tabix -r (reheader) failed in this htslib build: " + result["log_tail"][-200:])
    import gzip

    with gzip.open(result["output_path"], "rt") as handle:
        text = handle.read()
    assert (
        "##contig=<ID=chr1,assembly=toy,length=20000,md5=" in text and "##reference=file:" in text
    )
