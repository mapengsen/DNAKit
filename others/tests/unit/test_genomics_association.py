"""Known-value checks for the PLINK .ped reader, the GCTA GRM and the P3D LMM GWAS scan."""

from __future__ import annotations

import json
import math

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.genomics import genomic_relationship_matrix, mixed_model_gwas, read_plink_ped
from dnakit.genomics.association import _brent_fmin

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

PED = """F1 s1 0 0 0 -9 A A 0 0 C T
F2 s2 0 0 0 -9 A G G G C C
F3 s3 0 0 0 -9 G G T T T T
F4 s4 0 0 0 -9 0 0 G G C C
"""
MAP = "1\trs1\t0\t100\n1\trs2\t0\t200\n2\trs3\t0\t50\n"


def _write(tmp_path):
    (tmp_path / "toy.ped").write_text(PED)
    (tmp_path / "toy.map").write_text(MAP)
    return tmp_path / "toy.ped"


def test_read_plink_ped_minor_coding(tmp_path) -> None:
    result = read_plink_ped(_write(tmp_path), verbose=False)
    assert result["sample_ids"] == ["s1", "s2", "s3", "s4"]
    assert result["variant_ids"] == ["rs1", "rs2", "rs3"]
    assert result["chromosomes"] == ["1", "1", "2"] and result["positions"] == [100, 200, 50]
    # rs1: A x3, G x3 (tie -> 'G' sorts last); rs2: G x4, T x2; rs3: C x5, T x3
    assert result["counted_alleles"] == ["G", "T", "T"]
    assert result["genotypes"] == [[0.0, None, 1.0], [1.0, 0.0, 0.0], [2.0, 2.0, 2.0], [None, 0.0, 0.0]]
    assert result["missing_counts"] == [1, 1, 0]
    json.dumps(result, allow_nan=False)


def test_read_plink_ped_bglr_coding_reproduces_first_sample_reference(tmp_path) -> None:
    result = read_plink_ped(_write(tmp_path), coding="bglr", verbose=False)
    # rs2 is missing in the first sample: BGLR codes every homozygote as 2 (looks monomorphic)
    assert result["genotypes"] == [[0.0, None, 1.0], [1.0, 2.0, 0.0], [2.0, 2.0, 2.0], [None, 2.0, 0.0]]
    assert result["bglr_first_sample_missing_indices"] == [1]
    assert result["counted_alleles"][0] == "G" and result["other_alleles"][1] is None
    summary = read_plink_ped(_write(tmp_path), coding="bglr", include_genotypes=False, verbose=False)
    assert "genotypes" not in summary and summary["n_missing_genotypes"] == 2


def test_read_plink_ped_rejects_mismatched_map(tmp_path) -> None:
    path = _write(tmp_path)
    (tmp_path / "toy.map").write_text("1\trs1\t0\t100\n")
    with pytest.raises(ConfigurationError):
        read_plink_ped(path, verbose=False)
    with pytest.raises(ConfigurationError):
        read_plink_ped(path, coding="plink2", verbose=False)


def test_gcta_grm_matches_pairwise_definition() -> None:
    rng = np.random.default_rng(3)
    g = rng.integers(0, 3, (6, 40)).astype(float)
    g[rng.random(g.shape) < 0.1] = np.nan
    result = genomic_relationship_matrix(g, method="gcta", eigenvectors=2, verbose=False)
    keep = result["marker_indices"]
    x = g[:, keep]
    p = np.nanmean(x, axis=0) / 2
    z = (x - 2 * p) / np.sqrt(p * (1 - p))
    expected = np.zeros((6, 6))
    for j in range(6):
        for k in range(6):
            both = ~np.isnan(z[j]) & ~np.isnan(z[k])
            expected[j, k] = np.sum(z[j, both] * z[k, both]) / (2 * both.sum())
    np.testing.assert_allclose(result["relationship_matrix"], expected, rtol=1e-12)
    vectors = np.asarray(result["eigenvectors"])
    np.testing.assert_allclose(np.asarray(result["relationship_matrix"]) @ vectors[:, 0],
                               result["eigenvalues"][0] * vectors[:, 0], atol=1e-10)
    assert abs(sum(result["explained_variance_ratio"]) - sum(result["eigenvalues"][:2]) / sum(result["eigenvalues"])) < 1e-12
    json.dumps(result, allow_nan=False)


def test_brent_fmin_reaches_minimum() -> None:
    assert abs(_brent_fmin(lambda v: (math.log(v) - 1.3) ** 2, 1e-9, 1e9, 1e-10) - math.exp(1.3)) < 1e-5


def _toy_lmm(seed: int = 7, n: int = 60, m: int = 25):
    rng = np.random.default_rng(seed)
    g = rng.integers(0, 3, (n, m)).astype(float)
    z = (g - g.mean(0)) / g.std(0)
    k = z @ z.T / m + 1e-6 * np.eye(n)
    y = 2.0 * g[:, 0] + rng.multivariate_normal(np.zeros(n), k) + rng.normal(size=n)
    return g, y, k


def _reference_f(y, x, h):
    """Dense GLS on [1, x] with H, residual variance re-estimated: two-sided t(n-2) P value."""
    stats = pytest.importorskip("scipy.stats")
    hi = np.linalg.inv(h)
    design = np.column_stack([np.ones(len(y)), x])
    w = design.T @ hi @ design
    beta = np.linalg.solve(w, design.T @ hi @ y)
    resid = y - design @ beta
    s2 = resid @ hi @ resid / (len(y) - 2)
    t = beta[1] / math.sqrt(s2 * np.linalg.inv(w)[1, 1])
    return 2 * stats.t.sf(abs(t), len(y) - 2)


def test_mixed_model_gwas_matches_dense_gls_and_missing_modes() -> None:
    g, y, k = _toy_lmm()
    g[3, 1] = g[10, 1] = np.nan
    y[5] = np.nan
    results = {
        mode: mixed_model_gwas(g, y, k, missing_genotypes=mode, verbose=False)
        for mode in ("exact", "rrblup", "mean")
    }
    keep = ~np.isnan(y)
    lam = results["exact"]["lambda"]
    h = k[np.ix_(keep, keep)] + lam * np.eye(int(keep.sum()))
    # complete marker: all three modes equal the dense GLS F test
    for mode in results:
        assert results[mode]["records"][0]["p"] == pytest.approx(_reference_f(y[keep], g[keep, 0], h), rel=1e-9)
    # marker with missing genotypes: 'exact' re-inverts H on the observed samples
    obs = ~np.isnan(g[keep, 1])
    exact = _reference_f(y[keep][obs], g[keep, 1][obs], h[np.ix_(obs, obs)])
    assert results["exact"]["records"][1]["p"] == pytest.approx(exact, rel=1e-9)
    assert results["exact"]["records"][1]["n"] == int(obs.sum())
    assert results["rrblup"]["records"][1]["p"] != pytest.approx(exact, rel=1e-9)
    filled = g[keep, 1].copy()
    filled[~obs] = np.nanmean(filled)
    assert results["mean"]["records"][1]["p"] == pytest.approx(_reference_f(y[keep], filled, h), rel=1e-9)
    assert results["exact"]["n_phenotyped"] == int(keep.sum())
    assert results["exact"]["significant_indices"][0] == 0
    json.dumps(results["exact"], allow_nan=False)


def test_mixed_model_gwas_reml_maximises_restricted_likelihood() -> None:
    g, y, k = _toy_lmm(seed=11)
    result = mixed_model_gwas(g, y, k, verbose=False)
    n = len(y)
    s = np.eye(n) - 1.0 / n

    def neg_reml(lam: float) -> float:
        h = k + lam * np.eye(n)
        hi = np.linalg.inv(h)
        p = hi - hi @ np.ones((n, 1)) @ np.ones((1, n)) @ hi / hi.sum()
        return (n - 1) * math.log(y @ p @ y) + np.linalg.slogdet(h)[1] + math.log(hi.sum())

    lam = result["lambda"]
    assert neg_reml(lam) <= min(neg_reml(lam * 0.98), neg_reml(lam * 1.02)) + 1e-9
    assert 0 < result["heritability"] < 1
    assert s.shape == (n, n)


def test_mixed_model_gwas_reports_untestable_markers_and_validates() -> None:
    g, y, k = _toy_lmm(seed=5, m=4)
    g[:, 2] = 1.0
    result = mixed_model_gwas(g, y, k, min_maf=0.0, variant_ids=["a", "b", "c", "d"], verbose=False)
    assert result["records"][2]["status"] == "monomorphic_or_all_missing"
    assert result["records"][2]["p"] is None and result["n_tested"] == 3
    assert result["bonferroni_threshold"] == pytest.approx(0.05 / 3)
    with pytest.raises(ConfigurationError):
        mixed_model_gwas(g, y, k, missing_genotypes="impute", verbose=False)
    with pytest.raises(ConfigurationError):
        mixed_model_gwas(g, y[:-1], k, verbose=False)
    with pytest.raises(ConfigurationError):
        mixed_model_gwas(g, y, -np.eye(len(y)), verbose=False)
