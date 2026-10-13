"""Contract tests for optional Hail MatrixTable, QC, LD, PCA and association adapters."""

from __future__ import annotations

import sys
from fractions import Fraction
from math import factorial, sqrt
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants.hail_backends import (
    calculate_association_tests,
    calculate_hardy_weinberg_pvalues,
    calculate_hwe_normalized_pca,
    calculate_sample_qc,
    calculate_variant_qc,
    import_vcf_matrix_table,
    prune_linkage_disequilibrium,
    summarize_vcf_matrix_table,
)

# Four samples, three biallelic chr1 sites and one multi-allelic chr2 site.
VCF_TEXT = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=248956422>\n"
    "##contig=<ID=chr2,length=242193529>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
    "chr1\t10\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0\t./.\n"
    "chr1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\t0/0\t0/1\t1/1\n"
    "chr1\t30\trs3\tG\tC\t99\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n"
    "chr2\t5\trs4\tA\tC,T\t60\tPASS\t.\tGT\t0/1\t1/2\t./.\t0/0\n"
)
PHENOTYPE_TEXT = "sample\ty\tcase\nS1\t1.0\t0\nS2\t2.5\t1\nS3\t3.0\t0\nS4\t2.0\t1\n"

# Alternate-allele dosages per site (None = missing), in sample order S1..S4.
DOSAGES = {10: [1, 2, 0, None], 20: [1, 0, 1, 2], 30: [0, 1, 2, 1]}
QUANTITATIVE = [1.0, 2.5, 3.0, 2.0]
BINARY = [0.0, 1.0, 0.0, 1.0]


@pytest.fixture()
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    vcf = tmp_path / "cohort.vcf"
    vcf.write_text(VCF_TEXT, encoding="utf-8")
    phenotype = tmp_path / "phenotype.tsv"
    phenotype.write_text(PHENOTYPE_TEXT, encoding="utf-8")
    return vcf, phenotype


def _hail():
    return pytest.importorskip("hail")


def _levene_haldane_midp(n_hom_ref: int, n_het: int, n_hom_var: int, one_sided: bool) -> Fraction:
    """Exact Hardy-Weinberg mid-p computed from the Levene-Haldane distribution."""
    n = n_hom_ref + n_het + n_hom_var
    n_a = 2 * n_hom_ref + n_het
    n_b = 2 * n_hom_var + n_het

    def pmf(h: int) -> Fraction:
        aa, bb = (n_a - h) // 2, (n_b - h) // 2
        return Fraction(
            factorial(n) * factorial(n_a) * factorial(n_b) * 2**h,
            factorial(2 * n) * factorial(aa) * factorial(h) * factorial(bb),
        )

    support = {h: pmf(h) for h in range(min(n_a, n_b) + 1) if (n_a - h) % 2 == 0}
    assert sum(support.values()) == 1
    observed = support[n_het]
    if one_sided:
        upper = [p for h, p in support.items() if h > n_het]
        return sum(upper, Fraction(0)) + observed / 2
    less = sum((p for p in support.values() if p < observed), Fraction(0))
    return less + observed / 2


def _impute(values: list[int | None]) -> list[float]:
    observed = [value for value in values if value is not None]
    mean = sum(observed) / len(observed)
    return [float(mean if value is None else value) for value in values]


def _r_squared(x: list[float], y: list[float]) -> float:
    n = len(x)
    mean_x, mean_y = sum(x) / n, sum(y) / n
    sxy = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y, strict=True))
    sxx = sum((a - mean_x) ** 2 for a in x)
    syy = sum((b - mean_y) ** 2 for b in y)
    return sxy * sxy / (sxx * syy)


def _linear_hand(x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Intercept-only OLS with four samples: beta, standard error and two-sided p-value."""
    n = len(x)
    assert n == 4
    mean_x, mean_y = sum(x) / n, sum(y) / n
    sxx = sum((a - mean_x) ** 2 for a in x)
    beta = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y, strict=True)) / sxx
    intercept = mean_y - beta * mean_x
    ssr = sum((b - (intercept + beta * a)) ** 2 for a, b in zip(x, y, strict=True))
    standard_error = sqrt(ssr / (n - 2) / sxx)
    t_stat = beta / standard_error
    # Two-sided p-value of Student's t with 2 degrees of freedom.
    p_value = 1 - abs(t_stat) / sqrt(t_stat * t_stat + 2)
    return beta, standard_error, p_value


def test_missing_hail_reports_optional_extra(monkeypatch: pytest.MonkeyPatch, inputs) -> None:
    vcf, _ = inputs
    monkeypatch.setitem(sys.modules, "hail", None)
    with pytest.raises(BackendUnavailableError, match="bio-hail"):
        summarize_vcf_matrix_table(vcf)


def test_invalid_arguments_are_rejected_before_hail_starts(tmp_path: Path, inputs) -> None:
    vcf, phenotype = inputs
    cases = [
        (lambda: calculate_variant_qc(vcf, max_variants=True), "max_variants"),
        (lambda: calculate_variant_qc(vcf, max_variants=0), "max_variants"),
        (lambda: calculate_variant_qc(vcf, reference_genome="hg19"), "reference_genome"),
        (lambda: summarize_vcf_matrix_table(vcf, max_contigs=0), "max_contigs"),
        (lambda: calculate_hardy_weinberg_pvalues(vcf, one_sided=1), "one_sided"),
        (lambda: prune_linkage_disequilibrium(vcf, r2=0), "r2"),
        (lambda: prune_linkage_disequilibrium(vcf, r2=True), "r2"),
        (lambda: prune_linkage_disequilibrium(vcf, bp_window_size=0), "bp_window_size"),
        (lambda: calculate_hwe_normalized_pca(vcf, n_components=0), "n_components"),
        (lambda: import_vcf_matrix_table(vcf, tmp_path), "already exists"),
        (
            lambda: calculate_association_tests(vcf, phenotype, phenotype_column="y", model="x"),
            "model",
        ),
    ]
    for call, message in cases:
        with pytest.raises(ConfigurationError, match=message):
            call()


@pytest.mark.parametrize(
    ("text", "kwargs", "message"),
    [
        ("sample\ty\nS1\t1\nS1\t2\n", {"phenotype_column": "y"}, "Duplicate"),
        ("sample\ty\nS1\tabc\n", {"phenotype_column": "y"}, "not numeric"),
        ("sample\tcase\nS1\t2\n", {"phenotype_column": "case", "model": "logistic"}, "0, 1"),
        ("sample\ty\nS1\t1\n", {"phenotype_column": "missing"}, "must contain"),
        ("sample\ty\nS1\t1\nS2\t2\n", {"phenotype_column": "y", "max_samples": 1}, "max_samples"),
    ],
)
def test_invalid_phenotype_tables_are_rejected(
    tmp_path: Path, inputs, text: str, kwargs: dict, message: str
) -> None:
    vcf, _ = inputs
    phenotype = tmp_path / "bad.tsv"
    phenotype.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigurationError, match=message):
        calculate_association_tests(vcf, phenotype, **kwargs)


def test_import_summary_and_matrix_table_round_trip(tmp_path: Path, inputs) -> None:
    _hail()
    vcf, _ = inputs
    written = import_vcf_matrix_table(vcf, tmp_path / "cohort.mt")
    assert written["n_variants"] == 4
    assert written["n_samples"] == 4

    from_vcf = summarize_vcf_matrix_table(vcf)
    from_matrix = summarize_vcf_matrix_table(tmp_path / "cohort.mt")
    assert from_vcf == from_matrix
    assert from_vcf["biallelic_variants"] == 3
    assert from_vcf["multiallelic_variants"] == 1
    assert from_vcf["snv_variants"] == 3
    assert from_vcf["contigs"] == [
        {"contig": "chr1", "variants": 3},
        {"contig": "chr2", "variants": 1},
    ]
    assert from_vcf["samples"] == ["S1", "S2", "S3", "S4"]
    with pytest.raises(ConfigurationError, match="max_contigs"):
        summarize_vcf_matrix_table(vcf, max_contigs=1)


def test_variant_and_sample_qc_match_hand_counts(inputs) -> None:
    _hail()
    vcf, _ = inputs
    variants = {row["pos"]: row for row in calculate_variant_qc(vcf)["variants"]}
    assert sorted(variants) == [10, 20, 30]
    # Site 10: S4 is ./., so 3 of 4 samples are called; alt copies 1+2+0 out of 6 alleles.
    assert variants[10]["n_called"] == 3
    assert variants[10]["call_rate"] == pytest.approx(0.75)
    assert (variants[10]["allele_count"], variants[10]["allele_number"]) == (3, 6)
    assert variants[10]["allele_frequency"] == pytest.approx(0.5)
    assert (variants[10]["n_het"], variants[10]["n_hom_alt"]) == (1, 1)
    # Site 20: all called; alt copies 1+0+1+2 out of 8 alleles, two heterozygotes.
    assert variants[20]["call_rate"] == pytest.approx(1.0)
    assert variants[20]["allele_frequency"] == pytest.approx(0.5)
    assert (variants[20]["n_het"], variants[20]["n_hom_alt"]) == (2, 1)

    samples = {row["sample"]: row for row in calculate_sample_qc(vcf)["samples"]}
    assert samples["S3"]["call_rate"] == pytest.approx(0.75)
    assert (samples["S1"]["n_het"], samples["S1"]["n_hom_ref"]) == (3, 1)
    # Transition and transversion counts use alternate-allele copies (hom-alt counts twice).
    # S2: A>G hom-alt (2 transitions), G>C het (1 transversion), A>C/A>T at 1/2 (2 transversions).
    assert (samples["S2"]["n_transition"], samples["S2"]["n_transversion"]) == (2, 3)
    assert samples["S2"]["r_ti_tv"] == pytest.approx(2 / 3)

    with pytest.raises(ConfigurationError, match="max_variants"):
        calculate_variant_qc(vcf, max_variants=2)


def test_hardy_weinberg_matches_levene_haldane_mid_p(inputs) -> None:
    _hail()
    vcf, _ = inputs
    result = calculate_hardy_weinberg_pvalues(vcf)
    by_pos = {row["pos"]: row for row in result["variants"]}
    assert result["skipped_multiallelic"] == 1
    # Genotype counts (hom_ref, het, hom_var): site 10 is (1, 1, 1); sites 20 and 30 are (1, 2, 1).
    assert (by_pos[10]["n_hom_ref"], by_pos[10]["n_het"], by_pos[10]["n_hom_var"]) == (1, 1, 1)
    assert by_pos[10]["p_value"] == pytest.approx(float(_levene_haldane_midp(1, 1, 1, False)))
    assert by_pos[20]["p_value"] == pytest.approx(float(_levene_haldane_midp(1, 2, 1, False)))
    assert by_pos[30]["p_value"] == pytest.approx(23 / 35)
    assert by_pos[20]["het_freq_hwe"] == pytest.approx(4 / 7)

    one_sided = {
        row["pos"]: row for row in calculate_hardy_weinberg_pvalues(vcf, one_sided=True)["variants"]
    }
    assert one_sided[20]["p_value"] == pytest.approx(float(_levene_haldane_midp(1, 2, 1, True)))


def test_ld_pruning_keeps_an_independent_set(inputs) -> None:
    hl = _hail()
    vcf, _ = inputs
    result = prune_linkage_disequilibrium(vcf, r2=0.2)
    kept = {row["pos"] for row in result["retained_variants"]}
    assert result["n_input_variants"] == 3
    # Hand r-squared on mean-imputed dosages: r2(10,20)=0.25, r2(10,30)=0.25, r2(20,30)=0.
    r2 = {
        (a, b): _r_squared(_impute(DOSAGES[a]), _impute(DOSAGES[b]))
        for a, b in [(10, 20), (10, 30), (20, 30)]
    }
    assert r2[(10, 20)] == pytest.approx(0.25)
    assert r2[(20, 30)] == pytest.approx(0.0)
    assert len(kept) >= 1
    for a, b in [(10, 20), (10, 30), (20, 30)]:
        if a in kept and b in kept:
            assert r2[(a, b)] < 0.2
    for dropped in {10, 20, 30} - kept:
        assert any(r2[tuple(sorted((dropped, other)))] >= 0.2 for other in kept), (
            "a dropped site must be correlated with a retained site"
        )

    matrix = hl.import_vcf(str(vcf), reference_genome="GRCh38")
    matrix = matrix.filter_rows(hl.len(matrix.alleles) == 2)
    direct = hl.ld_prune(matrix.GT, r2=0.2, bp_window_size=1000000).collect()
    assert sorted(row.locus.position for row in direct) == sorted(kept)


def test_hwe_pca_matches_direct_hail_call(inputs) -> None:
    hl = _hail()
    vcf, _ = inputs
    result = calculate_hwe_normalized_pca(vcf, n_components=2)
    matrix = hl.import_vcf(str(vcf), reference_genome="GRCh38")
    matrix = matrix.filter_rows(hl.len(matrix.alleles) == 2)
    eigenvalues, scores, _ = hl.hwe_normalized_pca(matrix.GT, k=2, compute_loadings=False)
    assert result["eigenvalues"] == pytest.approx(list(eigenvalues))
    direct = {row.s: list(row.scores) for row in scores.collect()}
    assert [item["sample"] for item in result["samples"]] == sorted(direct)
    for item in result["samples"]:
        assert item["scores"] == pytest.approx(direct[item["sample"]])


def test_linear_association_matches_closed_form(inputs) -> None:
    _hail()
    vcf, phenotype = inputs
    result = calculate_association_tests(vcf, phenotype, phenotype_column="y")
    assert result["n_samples_with_phenotype"] == 4
    assert result["covariates"] == []
    by_pos = {row["pos"]: row for row in result["results"]}
    for pos in (10, 20, 30):
        beta, standard_error, p_value = _linear_hand(_impute(DOSAGES[pos]), QUANTITATIVE)
        row = by_pos[pos]
        assert row["beta"] == pytest.approx(beta, abs=1e-9)
        assert row["standard_error"] == pytest.approx(standard_error, abs=1e-9)
        assert row["p_value"] == pytest.approx(p_value, abs=1e-9)


def test_logistic_association_reports_separation_and_zero_effect(inputs) -> None:
    _hail()
    vcf, phenotype = inputs
    result = calculate_association_tests(vcf, phenotype, phenotype_column="case", model="logistic")
    by_pos = {row["pos"]: row for row in result["results"]}
    # Site 10: cases at dosage 1 and 2 and controls at 0 and 1 give quasi-complete separation,
    # so no finite maximum-likelihood estimate exists and Hail reports non-convergence.
    assert by_pos[10]["converged"] is False
    assert by_pos[10]["beta"] is None
    # Sites 20 and 30: the score at beta=0 is zero, so beta=0, p=1 and SE=sqrt(2) by hand.
    for pos in (20, 30):
        assert by_pos[pos]["converged"] is True
        assert by_pos[pos]["beta"] == pytest.approx(0.0, abs=1e-8)
        assert by_pos[pos]["standard_error"] == pytest.approx(sqrt(2), abs=1e-6)
        assert by_pos[pos]["p_value"] == pytest.approx(1.0, abs=1e-8)
