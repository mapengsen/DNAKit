"""Tests for dnakit.microbiome.skbio_p20_stats_functions (distance tests, ordination, alpha diversity)."""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("skbio")

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import skbio_p20_stats_functions as stats

IDS = [f"s{i}" for i in range(8)]
GROUPS = ["a"] * 4 + ["b"] * 4


@pytest.fixture(scope="module")
def data():
    rng = np.random.default_rng(0)
    counts = rng.poisson(5, (8, 6))
    from skbio.diversity import beta_diversity

    distance = beta_diversity("braycurtis", counts, ids=IDS).data.tolist()
    env = {"p": rng.normal(size=8).tolist(), "q": rng.normal(size=8).tolist()}
    return counts, distance, env


def test_anosim_permdisp_match_scikit_bio(data):
    from skbio import DistanceMatrix
    from skbio.stats.distance import anosim, permdisp

    _, distance, _ = data
    ours = stats.anosim_test(distance, GROUPS, ids=IDS, permutations=99, seed=1, verbose=False)
    ref = anosim(DistanceMatrix(distance, ids=IDS), GROUPS, permutations=99, seed=1)
    assert ours["test_statistic"] == pytest.approx(float(ref["test statistic"]))
    assert ours["p_value"] == pytest.approx(float(ref["p-value"]))
    disp = stats.permdisp_test(distance, GROUPS, ids=IDS, permutations=99, seed=1, verbose=False)
    ref = permdisp(DistanceMatrix(distance, ids=IDS), GROUPS, permutations=99, seed=1, dimensions=7)
    assert (
        disp["test_statistic"] == pytest.approx(float(ref["test statistic"]))
        and disp["dimensions"] == 7
    )
    with pytest.raises(ConfigurationError):
        stats.anosim_test(distance, GROUPS[:-1], ids=IDS, verbose=False)


def test_mantel_and_bioenv(data):
    _, distance, env = data
    single = stats.mantel_test(
        [distance, distance], ids=IDS, permutations=49, seed=1, verbose=False
    )
    assert single["coefficient"] == pytest.approx(1.0) and single["n"] == 8
    pairs = stats.mantel_test(
        [distance] * 3, ids=IDS, permutations=9, seed=1, labels=["x", "y", "z"], verbose=False
    )
    assert [r[:2] for r in pairs["rows"]] == [["x", "y"], ["x", "z"], ["y", "z"]]
    search = stats.bioenv_search(distance, env, ids=IDS, verbose=False)
    assert search["columns"] == ["vars", "size", "correlation"] and len(search["rows"]) == 2
    with pytest.raises(ConfigurationError):
        stats.mantel_test([distance], ids=IDS, verbose=False)


def test_ordination_functions(data, tmp_path):
    counts, distance, env = data
    biplot = stats.pcoa_biplot_projection(distance, env, ids=IDS, verbose=False)
    assert set(biplot["ordination"]["features"]["index"]) == {"p", "q"}
    ca = stats.correspondence_analysis((counts + 1).tolist(), verbose=False)
    assert len(ca["ordination"]["eigvals"]) == 5
    rda = stats.constrained_ordination((counts + 1).tolist(), env, verbose=False)
    cca = stats.constrained_ordination((counts + 1).tolist(), env, method="cca", verbose=False)
    assert rda["ordination"]["method"] == "RDA" and cca["ordination"]["method"] == "CCA"
    plot = stats.plot_ordination(
        distance, str(tmp_path / "o.png"), ids=IDS, grouping=GROUPS, verbose=False
    )
    assert (tmp_path / "o.png").stat().st_size > 1000 and plot["axes"] == [0, 1]
    with pytest.raises(ConfigurationError):
        stats.plot_ordination(distance, str(tmp_path / "o.txt"), ids=IDS, verbose=False)


def test_mmvec_fit_predict_score(data):
    pytest.importorskip("torch")
    counts, _, _ = data
    other = np.random.default_rng(2).poisson(5, (8, 4))
    result = stats.mmvec_analysis(
        counts.tolist(),
        other.tolist(),
        dimensions=2,
        max_iter=30,
        seed=1,
        predict_x=counts[:2].tolist(),
        score_x=counts.tolist(),
        score_y=other.tolist(),
        verbose=False,
    )
    assert len(result["ranks"]["data"]) == 6 and len(result["ranks"]["data"][0]) == 4
    assert len(result["predicted_y_proportions"]["data"]) == 2
    assert result["q_squared"] is None or isinstance(result["q_squared"], float)
    with pytest.raises(ConfigurationError):
        stats.mmvec_analysis(
            counts.tolist(), other.tolist(), score_x=counts.tolist(), verbose=False
        )


def test_alpha_metrics_match_scikit_bio(data):
    import skbio.diversity.alpha as alpha

    counts = [1, 1, 2, 3, 5, 8, 0, 13, 1, 21]
    table = [counts, counts[::-1]]
    expected = {
        "ace": alpha.ace(counts),
        "chao1": alpha.chao1(counts),
        "berger_parker_d": alpha.berger_parker_d(counts),
        "brillouin_d": alpha.brillouin_d(counts),
        "dominance": alpha.dominance(counts),
        "doubles": alpha.doubles(counts),
        "inv_simpson": alpha.inv_simpson(counts),
        "fisher_alpha": alpha.fisher_alpha(counts),
        "gini_index": alpha.gini_index(counts),
        "goods_coverage": alpha.goods_coverage(counts),
        "heip_e": alpha.heip_e(counts),
        "hill": alpha.hill(counts),
        "kempton_taylor_q": alpha.kempton_taylor_q(counts),
    }
    for metric, value in expected.items():
        got = stats.alpha_diversity_metric(table, metric, verbose=False)["values"]
        assert got[0] == pytest.approx(float(value)), metric
    assert stats.alpha_diversity_metric(table, "enspie", verbose=False)["values"][
        0
    ] == pytest.approx(alpha.enspie(counts))
    assert stats.alpha_diversity_metric(table, "hill", order=0, verbose=False)["values"][0] == 9
    assert stats.alpha_diversity_metric(table, "chao1_ci", verbose=False)["values"][
        0
    ] == pytest.approx(list(alpha.chao1_ci(counts)))
    assert stats.alpha_diversity_metric(table, "esty_ci", verbose=False)["values"][
        0
    ] == pytest.approx(list(alpha.esty_ci(counts)))
    lladser = stats.alpha_diversity_metric(table, "lladser_pe", seed=1, verbose=False)["values"]
    assert all(0.0 <= value <= 1.0 for value in lladser)
    faith = stats.alpha_diversity_metric(
        [[1, 0, 2, 3]],
        "faith_pd",
        feature_ids=list("abcd"),
        tree_newick="((a:1,b:2):1,(c:3,d:1):2);",
        verbose=False,
    )
    assert faith["values"] == [8.0]
    with pytest.raises(ConfigurationError):
        stats.alpha_diversity_metric(table, "faith_pd", verbose=False)
    with pytest.raises(ConfigurationError):
        stats.alpha_diversity_metric([[1.5, 2]], "chao1", verbose=False)
    with pytest.raises(ConfigurationError):
        stats.alpha_diversity_metric(table, "bogus", verbose=False)
