"""Tests for dnakit.microbiome.skbio_p20_composition_functions (compositional data analysis)."""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("skbio")

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import skbio_p20_composition_functions as comp

X = [[1, 2, 3], [4, 5, 6]]


@pytest.fixture(scope="module")
def table():
    rng = np.random.default_rng(1)
    groups = ["a"] * 4 + ["b"] * 4 + ["c"] * 4
    counts = rng.poisson(20, (12, 6)) + 1
    counts[4:8, 0] += 40
    return counts.tolist(), groups, {"grp": groups, "cont": rng.normal(size=12).tolist()}


def test_basic_operations_match_scikit_bio():
    import skbio.stats.composition as ref

    assert np.allclose(comp.closure_compositions(X, verbose=False)["compositions"], ref.closure(X))
    assert np.allclose(
        comp.replace_zeros_multiplicative([[1, 0, 3], [4, 5, 6]], verbose=False)["compositions"],
        ref.multi_replace([[1, 0, 3], [4, 5, 6]]),
    )
    assert np.allclose(
        comp.perturb_compositions(X[0], X[1], verbose=False)["composition"], ref.perturb(X[0], X[1])
    )
    assert np.allclose(
        comp.perturb_compositions(X[0], X[1], inverse=True, verbose=False)["composition"],
        ref.perturb_inv(X[0], X[1]),
    )
    assert np.allclose(
        comp.power_composition(X[0], 2, verbose=False)["composition"], ref.power(X[0], 2)
    )
    assert comp.aitchison_inner_product(X[0], X[1], verbose=False)[
        "inner_product"
    ] == pytest.approx(float(ref.inner(X[0], X[1])))
    assert np.allclose(
        comp.centralize_compositions(X, verbose=False)["compositions"], ref.centralize(X)
    )


@pytest.mark.parametrize("method", ["clr", "rclr", "ilr", "alr"])
def test_logratio_round_trip(method):
    import skbio.stats.composition as ref

    forward = comp.logratio_transform(X, method, verbose=False)
    assert np.allclose(forward["values"], getattr(ref, method)(X))
    back = comp.logratio_transform(forward["values"], method, inverse=True, verbose=False)
    assert np.allclose(back["values"], ref.closure(X))


def test_logratio_options_and_errors():
    alr = comp.logratio_transform(X, "alr", reference_index=1, verbose=False)["values"]
    assert alr[0] == pytest.approx([float(np.log(1 / 2)), float(np.log(3 / 2))])
    with pytest.raises(ConfigurationError):
        comp.logratio_transform([[0, 1, 2]], "clr", verbose=False)
    with pytest.raises(ConfigurationError):
        comp.logratio_transform(X, "alr", reference_index=9, verbose=False)


def test_vlr_and_bases(table):
    counts, _, _ = table
    pairwise = comp.variance_log_ratio(counts, feature_ids=list("abcdef"), verbose=False)
    single = comp.variance_log_ratio(
        counts, feature_ids=list("abcdef"), feature_x="a", feature_y="b", verbose=False
    )
    assert pairwise["ids"] == list("abcdef") and pairwise["matrix"][0][1] == pytest.approx(
        single["vlr"]
    )
    sbp = comp.ilr_basis(sbp=[[1, 1, -1], [1, -1, 0]], verbose=False)["basis"]
    assert np.allclose(np.array(sbp) @ np.array(sbp).T, np.eye(2))
    tree = comp.ilr_basis(newick="((a,b)x,c)r;", verbose=False)
    assert tree["balance_nodes"] == ["r", "x"]
    with pytest.raises(ConfigurationError):
        comp.ilr_basis(verbose=False)


def test_differential_abundance(table):
    counts, groups, meta = table
    ancom = comp.ancom_test(counts, groups, verbose=False)
    assert ancom["significant"][0] is True and not any(ancom["significant"][1:])
    bc = comp.differential_abundance_ancombc(counts, meta, "grp", verbose=False)
    bc2 = comp.differential_abundance_ancombc2(counts, meta, "grp", verbose=False)
    assert bc["columns"][:2] == ["FeatureID", "Covariate"]
    row = next(r for r in bc2["rows"] if r[0] == "F1" and r[1] == "grp[T.b]")
    assert row[2] == pytest.approx(1.1044655960160505, rel=1e-6) and row[-1] is True
    for test in ("global", "pairwise", "dunnett", "trend"):
        out = comp.ancombc_posthoc_test(
            counts, meta, "grp", "grp", test, bootstraps=20, seed=1, verbose=False
        )
        assert out["test"] == test and out["rows"][0][0] == "F1"
    zeros = comp.detect_structural_zeros(counts, groups, verbose=False)
    assert zeros["groups"] == ["a", "b", "c"] and not any(any(r) for r in zeros["structural_zero"])
    ttest = comp.dirmult_ttest_test(
        counts, groups, treatment="b", reference="a", draws=16, seed=1, verbose=False
    )
    assert ttest["rows"][0][-1] is True
    lme = comp.dirmult_mixed_model(counts, meta, "cont", "grp", draws=8, seed=1, verbose=False)
    assert "Covariate" in lme["columns"]
    with pytest.raises(ConfigurationError):
        comp.ancom_test([[0, 1], [1, 2]], ["a", "b"], verbose=False)
    with pytest.raises(ConfigurationError):
        comp.ancombc_posthoc_test(counts, meta, "grp", "grp", "bogus", verbose=False)
