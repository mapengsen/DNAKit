from __future__ import annotations

import math

import pytest

pytest.importorskip("skbio")
pytest.importorskip("pandas")

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import skbio_p21_functions as m

TREE = "(((a:1,b:1):1,c:2):1,(d:1,e:1):2);"
COUNTS = [[1, 2, 0, 5, 1], [3, 0, 4, 1, 0], [0, 0, 1, 2, 9]]
TAXA = list("abcde")


def test_alpha_matches_skbio_for_every_metric():
    from skbio.diversity import alpha

    vector = [1, 2, 0, 5, 1, 7]
    tree = "((((a:1,b:1):1,c:2):1,(d:1,e:1):2):1,f:1);"
    for metric in m.ALPHA_METRICS:
        result = m.calculate_alpha_diversity(
            vector, metric=metric, taxa=list("abcdef"), tree_newick=tree, verbose=False
        )
        assert result["status"] == "completed"
        if metric == "osd":
            assert result["value"] == {"observed": 5, "singles": 2, "doubles": 1}
        elif metric == "phydiv":
            assert result["value"] == 10.0
        elif metric == "michaelis_menten_fit":
            assert (
                math.isfinite(result["value"]) and result["value"] >= 5
            )  # random subsampling, unseeded
        else:
            assert result["value"] == pytest.approx(float(getattr(alpha, metric)(vector)))


def test_alpha_options_and_edge_cases():
    assert m.calculate_alpha_diversity([5, 5], metric="shannon", base=2, verbose=False)[
        "value"
    ] == pytest.approx(1.0)
    assert m.calculate_alpha_diversity([1, 3], metric="simpson", finite=True, verbose=False)[
        "value"
    ] == pytest.approx(0.5)
    assert m.calculate_alpha_diversity([0, 0, 0], metric="sobs", verbose=False)["value"] == 0
    with pytest.raises(ConfigurationError):
        m.calculate_alpha_diversity([1, 2], metric="nope", verbose=False)
    with pytest.raises(ConfigurationError):
        m.calculate_alpha_diversity([1.5, 2], metric="shannon", verbose=False)
    with pytest.raises(ConfigurationError):
        m.calculate_alpha_diversity([1, 2], metric="phydiv", verbose=False)


def test_unifrac_matches_skbio():
    from skbio import TreeNode
    from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac

    tree = TreeNode.read([TREE])
    un = m.calculate_unifrac_distance(COUNTS, TAXA, TREE, verbose=False)["distance_matrix"]
    wt = m.calculate_unifrac_distance(
        COUNTS, TAXA, TREE, weighted=True, normalized=True, verbose=False
    )["distance_matrix"]
    assert un[0][1] == pytest.approx(unweighted_unifrac(COUNTS[0], COUNTS[1], TAXA, tree))
    assert wt[0][2] == pytest.approx(
        weighted_unifrac(COUNTS[0], COUNTS[2], TAXA, tree, normalized=True)
    )
    assert un[1][1] == 0.0
    with pytest.raises(ConfigurationError):
        m.calculate_unifrac_distance(COUNTS, TAXA[:3], TREE, verbose=False)


def test_rarefy_and_isubsample():
    a = m.rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)
    b = m.rarefy_count_vector([5, 3, 0, 2], 4, seed=1, verbose=False)
    assert a["counts"] == b["counts"] and sum(a["counts"]) == 4 and a["counts"][2] == 0
    assert (
        sum(m.rarefy_count_vector([1, 1], 10, replace=True, seed=2, verbose=False)["counts"]) == 10
    )
    with pytest.raises(ConfigurationError):
        m.rarefy_count_vector([1, 1], 10, verbose=False)
    chosen = m.subsample_items_per_bin(
        list("abcdefgh"), 2, bins=list("xxxyyyzz"), seed=1, verbose=False
    )["selected"]
    assert sorted(b for b, _ in chosen) == ["x", "x", "y", "y", "z", "z"]
    dropped = m.subsample_items_per_bin(
        list("abc"), 2, minimum=2, bins=["x", "y", "y"], seed=1, verbose=False
    )["selected"]
    assert {b for b, _ in dropped} == {"y"}
    with pytest.raises(ConfigurationError):
        m.subsample_items_per_bin(["a"], 1, minimum=2, verbose=False)


def test_hommola():
    host = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
    para = [[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 3], [1, 2, 3, 0]]
    link = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]]
    result = m.run_hommola_cospeciation_test(host, para, link, seed=1, verbose=False)
    again = m.run_hommola_cospeciation_test(host, para, link, seed=1, verbose=False)
    assert result["interaction_edges"] == 4 and result == again
    assert -1 <= result["correlation"] <= 1 and 0 < result["p_value"] <= 1
    with pytest.raises(ConfigurationError):
        m.run_hommola_cospeciation_test(host, para, [[1, 0, 0]], verbose=False)


def _gradient_inputs():
    coords = {f"s{i}": [0.1 * i, (i * 7 % 5) / 5, (i * 3 % 4) / 4] for i in range(1, 13)}
    meta = {k: {"g": "A" if i % 2 else "B", "t": i // 2} for i, k in enumerate(coords)}
    return coords, meta


@pytest.mark.parametrize("algorithm", m._GRADIENT_ALGORITHMS)
def test_gradient_anova_algorithms(algorithm):
    coords, meta = _gradient_inputs()
    result = m.calculate_gradient_anova(
        coords,
        [0.5, 0.3, 0.2],
        meta,
        algorithm=algorithm,
        sort_category="t",
        trajectory_categories=["g"],
        verbose=False,
    )
    groups = result["categories"][0]["groups"]
    assert set(groups) == {"A", "B"}
    assert all(math.isfinite(g["mean"]) for g in groups.values())


def test_gradient_anova_errors():
    coords, meta = _gradient_inputs()
    with pytest.raises(ConfigurationError):
        m.calculate_gradient_anova(coords, [0.5, 0.3, 0.2], meta, algorithm="x", verbose=False)
    with pytest.raises(ConfigurationError):
        m.calculate_gradient_anova(coords, [0.5, 0.3], meta, verbose=False)
    with pytest.raises(ConfigurationError):
        m.calculate_gradient_anova(coords, [0.5, 0.3, 0.2], meta, weighted=True, verbose=False)


TABLE = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]


def test_table_augmentation_shapes_and_seed():
    mix = m.mixup_table(TABLE, 3, labels=[0, 1, 0, 1], seed=1, verbose=False)
    assert len(mix["matrix"]) == 3 and len(mix["labels"][0]) == 2
    assert (
        mix["matrix"]
        == m.mixup_table(TABLE, 3, labels=[0, 1, 0, 1], seed=1, verbose=False)["matrix"]
    )
    assert len(m.mixup_table(TABLE, 3, append=True, seed=1, verbose=False)["matrix"]) == 7
    ait = m.aitchison_mixup_table(TABLE, 2, seed=1, verbose=False)
    assert all(sum(row) == pytest.approx(1.0) for row in ait["matrix"])
    cut = m.compositional_cutmix_table(TABLE, 2, labels=[0, 1, 0, 1], seed=1, verbose=False)
    assert all(sum(row) == pytest.approx(1.0) for row in cut["matrix"])
    phylo = m.phylomix_table(TABLE, 2, TREE, TAXA, seed=1, verbose=False)
    assert len(phylo["matrix"]) == 2 and len(phylo["matrix"][0]) == 5
    with pytest.raises(ConfigurationError):
        m.mixup_table(TABLE, 0, verbose=False)
    with pytest.raises(ConfigurationError):
        m.phylomix_table(TABLE, 2, TREE, TAXA[:3], verbose=False)


def test_sample_metadata():
    a = {"s1": {"x": 1, "c": "u"}, "s2": {"x": 5, "c": "v"}, "s3": {"x": 7, "c": "u"}}
    b = {"s1": {"y": 2.5}, "s2": {"y": 3.5}, "s4": {"y": 1.0}}
    merged = m.merge_sample_metadata([a, b], verbose=False)
    assert merged["ids"] == ["s1", "s2"] and merged["columns"] == ["x", "c", "y"]
    assert m.query_sample_ids(a, where='"x" > 2', verbose=False)["ids"] == ["s2", "s3"]
    assert m.query_sample_ids(a, verbose=False)["ids"] == ["s1", "s2", "s3"]
    with pytest.raises(ConfigurationError):
        m.query_sample_ids(a, where="nonsense ((", verbose=False)
    with pytest.raises(ConfigurationError):
        m.merge_sample_metadata([a], verbose=False)


def test_interval_metadata():
    first = {"upper_bound": 3, "intervals": [{"bounds": [[0, 2]], "metadata": {"gene": "sagA"}}]}
    second = {"upper_bound": 4, "intervals": [{"bounds": [[1, 4]], "metadata": {"gene": "sagB"}}]}
    hit = m.query_interval_metadata(
        first, bounds=[[1, 3]], attributes={"gene": "sagA"}, verbose=False
    )
    assert hit["n_matches"] == 1
    assert (
        m.query_interval_metadata(first, attributes={"gene": "recA"}, verbose=False)["n_matches"]
        == 0
    )
    joined = m.concat_interval_metadata([first, second], verbose=False)
    assert joined["upper_bound"] == 7 and joined["intervals"][1]["bounds"] == [[4, 7]]
    with pytest.raises(ConfigurationError):
        m.query_interval_metadata(
            {"upper_bound": 3, "intervals": [{"bounds": [[0, 9]]}]}, verbose=False
        )


def test_embedding_ordination():
    result = m.ordinate_embedding_vectors(
        {"a": [1, 2, 3], "b": [2, 1, 0], "c": [0, 5, 1], "d": [3, 3, 3]}, verbose=False
    )
    assert result["ids"] == ["a", "b", "c", "d"] and len(result["coordinates"]) == 4
    assert sum(result["proportion_explained"]) == pytest.approx(1.0)
    with pytest.raises(ConfigurationError):
        m.ordinate_embedding_vectors({"a": [1, 2], "b": [1]}, verbose=False)
