"""Contract tests for optional VCF backends and scikit-allel population statistics."""

from __future__ import annotations

from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import (
    calculate_allele_frequency,
    calculate_fst,
    calculate_nucleotide_diversity,
    calculate_observed_heterozygosity,
    calculate_pca,
    calculate_site_frequency_spectrum,
    calculate_tajimas_d,
    calculate_watterson_theta,
    calculate_windowed_ld,
    fetch_vcf_region,
    read_vcf,
)

VCF_TEXT = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t10\trs1\tA\tG\t50\tPASS\tDP=20;AF=0.5;DB\tGT\t0/1\t1/1\t0/0\t./.
chr1\t20\t.\tC\tT,A\t.\tq10\tDP=7;AF=0.2,0.1\tGT\t0|1\t0/0\t0/2\t1/2
chr1\t30\trs3\tG\tC\t99\t.\tDP=12\tGT\t0/0\t0/1\t1/1\t0/1
chr1\t40\trs4\tT\tA\t60\tPASS\tDP=9;AF=0.4\tGT\t1/1\t0/1\t0/0\t0/0
chr1\t50\trs5\tG\tA\t70\tPASS\tDP=11;AF=0.3\tGT\t0/1\t0/1\t./.\t0/0
"""


def _write(tmp_path: Path, text: str, name: str = "sample.vcf") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _site(pos: int, *calls: str, alt: str = "G") -> dict:
    return {"chrom": "chr1", "pos": pos, "ref": "A", "alt": [alt], "genotypes": list(calls)}


def _two_sample_dataset() -> dict:
    """Two diploid samples, four chromosomes; sites at 10 (alt=1), 30 (alt=3), 50 (monomorphic)."""

    return {
        "samples": ["S1", "S2"],
        "variants": [
            _site(10, "0/1", "0/0"),
            _site(30, "1/1", "0/1"),
            _site(50, "0/0", "0/0"),
        ],
    }


@pytest.mark.parametrize("backend", ["cyvcf2", "pysam", "vcfpy"])
def test_read_vcf_backends_return_identical_records(tmp_path: Path, backend: str) -> None:
    reference = read_vcf(_write(tmp_path, VCF_TEXT), backend="cyvcf2", verbose=False)
    module = {"cyvcf2": "cyvcf2", "pysam": "pysam", "vcfpy": "vcfpy"}[backend]
    pytest.importorskip(module)

    result = read_vcf(_write(tmp_path, VCF_TEXT), backend=backend, verbose=False)

    assert result["samples"] == ["S1", "S2", "S3", "S4"]
    assert result["variants"] == reference["variants"]


def test_read_vcf_normalizes_filters_missing_calls_and_phasing(tmp_path: Path) -> None:
    pytest.importorskip("cyvcf2")
    result = read_vcf(_write(tmp_path, VCF_TEXT), backend="cyvcf2", verbose=False)
    first, second, third = result["variants"][0], result["variants"][1], result["variants"][2]

    assert first["genotypes"] == ["0/1", "1/1", "0/0", "./."]
    assert first["filter"] == []  # PASS is represented by an empty list
    assert second["genotypes"][0] == "0|1"  # phased calls keep '|'
    assert second["filter"] == ["q10"]
    assert second["id"] is None and first["id"] == "rs1"
    assert first["info"] == {"DP": 20, "AF": 0.5, "DB": True}  # single-value INFO is a scalar
    assert second["info"]["AF"] == [0.2, 0.1]
    assert third["qual"] is None or third["qual"] == 99.0


def test_read_vcf_rejects_unknown_backend_and_excess_records(tmp_path: Path) -> None:
    path = _write(tmp_path, VCF_TEXT)
    with pytest.raises(ConfigurationError, match="backend must be one of"):
        read_vcf(path, backend="bcftools", verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=2"):
        read_vcf(path, backend="cyvcf2", max_records=2, verbose=False)


def test_missing_optional_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, VCF_TEXT)
    monkeypatch.setitem(__import__("sys").modules, "cyvcf2", None)
    with pytest.raises(BackendUnavailableError, match="bio-cyvcf2"):
        read_vcf(path, backend="cyvcf2", verbose=False)


def test_fetch_vcf_region_reads_only_requested_positions(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    plain = _write(tmp_path, VCF_TEXT)
    compressed = tmp_path / "sample.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)

    for backend in ("pysam", "cyvcf2"):
        if backend == "cyvcf2":
            pytest.importorskip("cyvcf2")
        result = fetch_vcf_region(compressed, "chr1:15-35", backend=backend, verbose=False)
        assert [variant["pos"] for variant in result["variants"]] == [20, 30]
        assert result["region"] == "chr1:15-35"


def test_fetch_vcf_region_requires_an_index_and_valid_region(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    path = _write(tmp_path, VCF_TEXT)
    with pytest.raises(ConfigurationError, match="region must look like"):
        fetch_vcf_region(path, "chr1:-5", verbose=False)
    with pytest.raises(ConfigurationError):
        fetch_vcf_region(path, "chr1:15-35", verbose=False)  # plain text VCF has no index


def test_observed_heterozygosity_matches_hand_calculation() -> None:
    # Site heterozygosity: 1/2, 1/2 and 0 -> mean 1/3.
    result = calculate_observed_heterozygosity(_two_sample_dataset(), verbose=False)
    population = result["populations"]["all"]
    assert population["variants"] == 3
    assert population["mean_observed_heterozygosity"] == pytest.approx(1 / 3)


def test_observed_heterozygosity_by_population() -> None:
    data = {
        "samples": ["S1", "S2"],
        "variants": [_site(10, "0/1", "0/0"), _site(20, "0/1", "0/1")],
    }
    result = calculate_observed_heterozygosity(
        data, populations={"A": ["S1"], "B": ["S2"]}, verbose=False
    )["populations"]
    assert result["A"]["mean_observed_heterozygosity"] == pytest.approx(1.0)
    assert result["B"]["mean_observed_heterozygosity"] == pytest.approx(0.5)


def test_watterson_theta_matches_closed_form() -> None:
    # S = 2 segregating sites, n = 4 chromosomes, a_n = 1 + 1/2 + 1/3 = 11/6.
    result = calculate_watterson_theta(_two_sample_dataset(), sequence_length=100, verbose=False)
    expected = (2 / (11 / 6)) / 100
    assert result["theta_per_base"] == pytest.approx(expected)


def test_nucleotide_diversity_matches_pairwise_count() -> None:
    # Site 10: 1 * 3 / C(4,2) = 0.5; site 30: 3 * 1 / 6 = 0.5; site 50 contributes 0.
    result = calculate_nucleotide_diversity(
        _two_sample_dataset(), sequence_length=100, verbose=False
    )
    assert result["pi_per_base"] == pytest.approx(1.0 / 100)


def _three_segregating_sites_dataset() -> dict:
    """Four chromosomes, three segregating sites: the minimum scikit-allel needs for Tajima's D."""

    return {
        "samples": ["S1", "S2"],
        "variants": [
            _site(10, "0/1", "0/0"),
            _site(30, "1/1", "0/1"),
            _site(40, "0/1", "0/1"),
        ],
    }


def test_tajimas_d_matches_scikit_allel_directly() -> None:
    allel = pytest.importorskip("allel")
    counts = allel.GenotypeArray(
        [[[0, 1], [0, 0]], [[1, 1], [0, 1]], [[0, 1], [0, 1]]]
    ).count_alleles()
    expected = allel.tajima_d(counts, pos=[10, 30, 40])

    result = calculate_tajimas_d(_three_segregating_sites_dataset(), verbose=False)

    assert result["tajimas_d"] == pytest.approx(float(expected))


def test_tajimas_d_is_none_below_three_segregating_sites() -> None:
    # scikit-allel returns NaN below min_sites=3 segregating sites; the wrapper reports None.
    result = calculate_tajimas_d(_two_sample_dataset(), verbose=False)
    assert result["tajimas_d"] is None


def test_site_frequency_spectrum_counts_alt_alleles() -> None:
    # alt counts 1, 3 and 0 across n = 4 chromosomes.
    result = calculate_site_frequency_spectrum(_two_sample_dataset(), verbose=False)
    assert result["chromosomes"] == 4
    assert result["spectrum"] == [1, 1, 0, 1, 0]


def test_fst_is_one_for_fixed_differences() -> None:
    data = {
        "samples": ["S1", "S2", "S3", "S4"],
        "variants": [
            _site(10, "0/0", "0/0", "1/1", "1/1"),
            _site(20, "0/0", "0/0", "1/1", "1/1"),
            _site(30, "0/0", "0/0", "1/1", "1/1"),
        ],
    }
    populations = {"A": ["S1", "S2"], "B": ["S3", "S4"]}
    hudson = calculate_fst(data, populations=populations, method="hudson", verbose=False)
    wc = calculate_fst(data, populations=populations, method="weir_cockerham", verbose=False)
    assert hudson["fst"] == pytest.approx(1.0)
    assert wc["fst"] == pytest.approx(1.0)


def test_weir_cockerham_fst_needs_two_samples_per_population() -> None:
    with pytest.raises(ConfigurationError, match="at least two samples"):
        calculate_fst(
            _two_sample_dataset(),
            populations={"A": ["S1"], "B": ["S2"]},
            method="weir_cockerham",
            verbose=False,
        )


def test_fst_matches_scikit_allel_directly() -> None:
    allel = pytest.importorskip("allel")
    data = {
        "samples": ["S1", "S2", "S3", "S4"],
        "variants": [
            _site(10, "0/1", "0/0", "1/1", "0/1"),
            _site(20, "0/0", "0/1", "0/1", "1/1"),
            _site(30, "1/1", "1/1", "0/0", "0/1"),
        ],
    }
    genotypes = allel.GenotypeArray(
        [
            [[0, 1], [0, 0], [1, 1], [0, 1]],
            [[0, 0], [0, 1], [0, 1], [1, 1]],
            [[1, 1], [1, 1], [0, 0], [0, 1]],
        ]
    )
    expected, _se, _blocks, _jack = allel.average_hudson_fst(
        genotypes.take([0, 1], axis=1).count_alleles(),
        genotypes.take([2, 3], axis=1).count_alleles(),
        1,
    )
    result = calculate_fst(
        data, populations={"A": ["S1", "S2"], "B": ["S3", "S4"]}, method="hudson", verbose=False
    )
    assert result["fst"] == pytest.approx(float(expected))


def test_weir_cockerham_fst_matches_scikit_allel_directly() -> None:
    allel = pytest.importorskip("allel")
    data = {
        "samples": ["S1", "S2", "S3", "S4"],
        "variants": [
            _site(10, "0/1", "0/0", "1/1", "0/1"),
            _site(20, "0/0", "0/1", "0/1", "1/1"),
            _site(30, "1/1", "1/1", "0/0", "0/1"),
        ],
    }
    genotypes = allel.GenotypeArray(
        [
            [[0, 1], [0, 0], [1, 1], [0, 1]],
            [[0, 0], [0, 1], [0, 1], [1, 1]],
            [[1, 1], [1, 1], [0, 0], [0, 1]],
        ]
    )
    expected, _se, _blocks, _jack = allel.average_weir_cockerham_fst(genotypes, [[0, 1], [2, 3]], 1)
    result = calculate_fst(
        data,
        populations={"A": ["S1", "S2"], "B": ["S3", "S4"]},
        method="weir_cockerham",
        verbose=False,
    )
    assert result["fst"] == pytest.approx(float(expected))


def test_fst_rejects_wrong_population_count() -> None:
    with pytest.raises(ConfigurationError, match="exactly two populations"):
        calculate_fst(
            _two_sample_dataset(),
            populations={"A": ["S1"], "B": ["S2"], "C": ["S1"]},
            verbose=False,
        )


def test_pca_matches_scikit_allel_directly() -> None:
    allel = pytest.importorskip("allel")
    rows = [
        ("0/0", "0/1", "1/1", "0/0"),
        ("0/1", "0/1", "1/1", "0/0"),
        ("1/1", "0/0", "0/1", "1/1"),
        ("0/1", "1/1", "0/0", "0/1"),
        ("1/1", "1/1", "0/1", "0/0"),
        ("0/0", "0/1", "0/1", "1/1"),
    ]
    data = {
        "samples": ["S1", "S2", "S3", "S4"],
        "variants": [_site(10 * (i + 1), *row) for i, row in enumerate(rows)],
    }
    result = calculate_pca(data, n_components=2, verbose=False)

    alt = allel.GenotypeArray([[[int(c[0]), int(c[2])] for c in row] for row in rows]).to_n_alt()
    coords, model = allel.pca(alt, n_components=2, scaler="patterson")

    assert len(result["coordinates"]) == 4
    assert result["explained_variance_ratio"] == pytest.approx(
        [float(v) for v in model.explained_variance_ratio_]
    )
    assert result["coordinates"]["S1"] == pytest.approx([float(v) for v in coords[0]])


def test_windowed_ld_is_one_for_identical_sites() -> None:
    calls = ("0/1", "1/1", "0/0", "0/1")
    data = {"samples": ["S1", "S2", "S3", "S4"], "variants": [_site(10, *calls), _site(20, *calls)]}
    result = calculate_windowed_ld(data, window_size=100, verbose=False)
    populated = [window for window in result["windows"] if window["variants"] == 2]
    assert populated and populated[0]["median_r_squared"] == pytest.approx(1.0)


def test_population_statistics_reject_bad_calls() -> None:
    data = {"samples": ["S1", "S2"], "variants": [_site(10, "0/x", "0/0")]}
    with pytest.raises(ConfigurationError, match="unsupported genotype call"):
        calculate_observed_heterozygosity(data, verbose=False)


def test_unknown_sample_name_is_rejected() -> None:
    with pytest.raises(ConfigurationError, match="unknown sample names"):
        calculate_observed_heterozygosity(
            _two_sample_dataset(), populations={"A": ["nope"]}, verbose=False
        )


def test_read_vcf_output_feeds_native_statistics(tmp_path: Path) -> None:
    pytest.importorskip("cyvcf2")
    records = read_vcf(_write(tmp_path, VCF_TEXT), backend="cyvcf2", verbose=False)
    result = calculate_allele_frequency(records, verbose=False)["result"]
    assert len(result["variants"]) == 5
    assert result["total_called_alleles"] > 0


def test_fetch_vcf_region_rejects_unknown_contig_on_every_backend(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    plain = _write(tmp_path, VCF_TEXT)
    compressed = tmp_path / "sample.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    for backend in ("pysam", "cyvcf2"):
        if backend == "cyvcf2":
            pytest.importorskip("cyvcf2")
        with pytest.raises(ConfigurationError, match="invalid contig"):
            fetch_vcf_region(compressed, "chrX:1-100", backend=backend, verbose=False)
