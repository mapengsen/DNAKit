"""Tests for the P19 GenePop editing and Bio.Cluster functions (no network)."""

from __future__ import annotations

import pytest

pytest.importorskip("Bio")
np = pytest.importorskip("numpy")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.population import biopython_p19_cluster_functions as m  # noqa: E402

GENEPOP = """Test genepop
Loc1
Loc2
Loc3
Pop
Ind1, 0101 0202 0303
Ind2, 0102 0203 0000
Pop
Ind3, 0201 0202 0303
Ind4, 0101 0102 0103
Pop
Ind5, 0202 0101 0101
"""
DATA = [[1.0, 2, 3], [1.1, 2.1, 3.2], [8, 9, 10], [8.2, 9.1, 9.9], [0.9, 2, 3.1]]


@pytest.fixture()
def genepop(tmp_path):
    path = tmp_path / "a.gen"
    path.write_text(GENEPOP)
    return path


def _read(path):
    from Bio.PopGen import GenePop

    with open(path) as handle:
        return GenePop.read(handle)


def test_split_by_population(genepop, tmp_path) -> None:
    result = m.split_genepop_by_population(genepop, tmp_path / "pops", verbose=False)
    assert set(result["files"]) == {"pop1", "pop2", "pop3"}
    second = _read(result["files"]["pop2"])
    assert [name for name, _ in second.populations[0]] == ["Ind3", "Ind4"]
    assert second.loci_list == ["Loc1", "Loc2", "Loc3"]
    named = m.split_genepop_by_population(genepop, tmp_path / "named", population_names=["x", "y", "z"], verbose=False)
    assert set(named["files"]) == {"x", "y", "z"}
    with pytest.raises(ConfigurationError):
        m.split_genepop_by_population(genepop, tmp_path / "bad", population_names=["x"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.split_genepop_by_population(genepop, tmp_path / "pops", verbose=False)  # exists


def test_split_by_locus(genepop, tmp_path) -> None:
    result = m.split_genepop_by_locus(genepop, tmp_path / "loci", verbose=False)
    assert set(result["files"]) == {"Loc1", "Loc2", "Loc3"}
    record = _read(result["files"]["Loc2"])
    assert record.loci_list == ["Loc2"]
    assert len(record.populations) == 3
    assert record.populations[0][1] == ("Ind2", [(2, 3)])


def test_remove_populations(genepop, tmp_path) -> None:
    result = m.remove_genepop_populations(genepop, tmp_path / "out.gen", [0, 2], verbose=False)
    assert result["population_count"] == 1
    record = _read(result["output"])
    assert [name for name, _ in record.populations[0]] == ["Ind3", "Ind4"]
    assert record.marker_len == 2
    with pytest.raises(ConfigurationError):
        m.remove_genepop_populations(genepop, tmp_path / "o2.gen", [0, 1, 2], verbose=False)
    with pytest.raises(ConfigurationError):
        m.remove_genepop_populations(genepop, tmp_path / "o3.gen", [7], verbose=False)


def test_remove_loci_by_name_and_position(genepop, tmp_path) -> None:
    by_name = m.remove_genepop_loci(genepop, tmp_path / "n.gen", names=["Loc2"], verbose=False)
    assert by_name["remaining_loci"] == ["Loc1", "Loc3"]
    by_position = m.remove_genepop_loci(genepop, tmp_path / "p.gen", positions=[0, 2], verbose=False)
    assert by_position["remaining_loci"] == ["Loc2"]
    record = _read(by_position["output"])
    assert record.populations[0][0] == ("Ind1", [(2, 2)])
    with pytest.raises(ConfigurationError):
        m.remove_genepop_loci(genepop, tmp_path / "x.gen", verbose=False)  # neither given
    with pytest.raises(ConfigurationError):
        m.remove_genepop_loci(genepop, tmp_path / "y.gen", names=["Nope"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.remove_genepop_loci(genepop, tmp_path / "z.gen", positions=[0, 1, 2], verbose=False)


def test_kmeans_and_kmedians_match_biopython() -> None:
    from Bio import Cluster

    for method, code, dist, name in (("mean", "a", "euclidean", "e"), ("median", "m", "cityblock", "b")):
        result = m.kmeans_cluster(
            DATA, 2, method=method, distance=dist, initial_clusters=[0, 0, 1, 1, 0], verbose=False
        )
        ids, error, found = Cluster.kcluster(
            np.array(DATA), 2, method=code, dist=name, initialid=np.array([0, 0, 1, 1, 0], dtype=np.int32)
        )
        assert result["cluster_ids"] == [int(i) for i in ids]
        assert result["error"] == pytest.approx(float(error))
    random_start = m.kmeans_cluster(DATA, 2, n_passes=5, verbose=False)
    assert sorted(set(random_start["cluster_ids"])) == [0, 1]
    assert random_start["cluster_ids"][0] == random_start["cluster_ids"][1] == random_start["cluster_ids"][4]
    with pytest.raises(ConfigurationError):
        m.kmeans_cluster(DATA, 9, verbose=False)
    with pytest.raises(ConfigurationError):
        m.kmeans_cluster([[1, 2], [3]], 1, verbose=False)
    with pytest.raises(ConfigurationError):
        m.kmeans_cluster([[1, float("nan")]], 1, verbose=False)


def test_kmedoids_matches_biopython() -> None:
    from Bio import Cluster

    matrix = m.profile_distance_matrix(DATA, verbose=False)["matrix"]
    result = m.kmedoids_cluster(matrix, 2, initial_clusters=[0, 0, 1, 1, 0], verbose=False)
    ids, error, _ = Cluster.kmedoids(np.array(matrix), 2, initialid=np.array([0, 0, 1, 1, 0], dtype=np.int32))
    assert result["cluster_ids"] == [int(i) for i in ids]
    assert result["error"] == pytest.approx(float(error))
    with pytest.raises(ConfigurationError):
        m.kmedoids_cluster([[0, 1], [1, 0], [2, 2]], 2, verbose=False)


@pytest.mark.parametrize("linkage,code", [("single", "s"), ("complete", "m"), ("average", "a"), ("centroid", "c")])
def test_hierarchical_matches_biopython(linkage, code) -> None:
    from Bio import Cluster

    result = m.hierarchical_cluster(DATA, linkage=linkage, verbose=False)
    tree = Cluster.treecluster(np.array(DATA), method=code, dist="e")
    assert [(n["left"], n["right"]) for n in result["nodes"]] == [(tree[i].left, tree[i].right) for i in range(len(tree))]
    assert [n["distance"] for n in result["nodes"]] == pytest.approx([tree[i].distance for i in range(len(tree))])
    assert result["item_count"] == 5


def test_hierarchical_from_distance_matrix_and_errors() -> None:
    matrix = m.profile_distance_matrix(DATA, verbose=False)["matrix"]
    result = m.hierarchical_cluster(distance_matrix=matrix, linkage="average", verbose=False)
    assert len(result["nodes"]) == 4
    with pytest.raises(ConfigurationError):
        m.hierarchical_cluster(distance_matrix=matrix, linkage="centroid", verbose=False)
    with pytest.raises(ConfigurationError):
        m.hierarchical_cluster(DATA, distance_matrix=matrix, verbose=False)
    with pytest.raises(ConfigurationError):
        m.hierarchical_cluster(verbose=False)


def test_cut_and_sort_tree() -> None:
    from Bio import Cluster

    nodes = m.hierarchical_cluster(DATA, linkage="average", verbose=False)["nodes"]
    cut = m.cut_cluster_tree(nodes, 2, verbose=False)
    tree = Cluster.treecluster(np.array(DATA), method="a", dist="e")
    assert cut["cluster_ids"] == [int(i) for i in tree.cut(2)]
    assert len(set(cut["cluster_ids"])) == 2
    ordered = m.sort_cluster_tree(nodes, verbose=False)
    assert sorted(ordered["leaf_order"]) == [0, 1, 2, 3, 4]
    assert ordered["leaf_order"] == [int(i) for i in tree.sort()]
    custom = m.sort_cluster_tree(nodes, order=[4, 3, 2, 1, 0], verbose=False)
    assert sorted(custom["leaf_order"]) == [0, 1, 2, 3, 4]
    with pytest.raises(ConfigurationError):
        m.cut_cluster_tree(nodes, 99, verbose=False)
    with pytest.raises(ConfigurationError):
        m.cut_cluster_tree([{"left": 0}], 1, verbose=False)


def test_som_shapes() -> None:
    result = m.som_cluster(DATA, nx_grid=2, ny_grid=1, n_iterations=30, verbose=False)
    assert len(result["cell_ids"]) == 5 and len(result["cell_data"]) == 2
    assert len(result["cell_data"][0][0]) == 3
    assert all(0 <= x < 2 and y == 0 for x, y in result["cell_ids"])


def test_cluster_distance_methods_match_biopython() -> None:
    from Bio import Cluster

    table = {"mean": "a", "median": "m", "single": "s", "complete": "x", "average": "v"}
    for method, code in table.items():
        result = m.cluster_distance(DATA, [0, 1], [2, 3], method=method, verbose=False)
        expected = Cluster.clusterdistance(np.array(DATA), index1=np.array([0, 1], dtype=np.int32), index2=np.array([2, 3], dtype=np.int32), method=code, dist="e")
        assert result["distance"] == pytest.approx(float(expected))
    with pytest.raises(ConfigurationError):
        m.cluster_distance(DATA, [0, 9], [2], verbose=False)


def test_cluster_centroids() -> None:
    mean = m.cluster_centroids(DATA, [0, 0, 1, 1, 0], verbose=False)
    assert mean["centroids"][0] == pytest.approx([(1 + 1.1 + 0.9) / 3, 2.0333333333, 3.1])
    median = m.cluster_centroids(DATA, [0, 0, 1, 1, 0], method="median", verbose=False)
    assert median["centroids"][0] == pytest.approx([1.0, 2.0, 3.1])
    assert mean["centroid_mask"] == [[1, 1, 1], [1, 1, 1]]
    with pytest.raises(ConfigurationError):
        m.cluster_centroids(DATA, [0, 1], verbose=False)


@pytest.mark.parametrize(
    "distance,code",
    [("euclidean", "e"), ("cityblock", "b"), ("pearson", "c"), ("abs_pearson", "a"), ("uncentered", "u"),
     ("abs_uncentered", "x"), ("spearman", "s"), ("kendall", "k")],
)
def test_distance_matrix_matches_biopython(distance, code) -> None:
    from Bio import Cluster

    result = m.profile_distance_matrix(DATA, distance=distance, verbose=False)
    lower = Cluster.distancematrix(np.array(DATA), dist=code)
    for i in range(5):
        for j in range(i):
            assert result["matrix"][i][j] == pytest.approx(float(lower[i][j]))
            assert result["matrix"][j][i] == result["matrix"][i][j]
    transposed = m.profile_distance_matrix(DATA, distance=distance, transpose=True, verbose=False)
    assert transposed["item_count"] == 3


def test_distance_matrix_with_mask() -> None:
    mask = [[1, 1, 1], [1, 0, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]
    masked = m.profile_distance_matrix(DATA, mask=mask, verbose=False)["matrix"]
    plain = m.profile_distance_matrix(DATA, verbose=False)["matrix"]
    assert masked[1][0] != plain[1][0]
    with pytest.raises(ConfigurationError):
        m.profile_distance_matrix(DATA, mask=[[1]], verbose=False)


def test_pca_reconstructs_data() -> None:
    result = m.cluster_pca(DATA, verbose=False)
    rebuilt = np.array(result["coordinates"]) @ np.array(result["components"]) + np.array(result["column_means"])
    assert np.allclose(rebuilt, np.array(DATA))
    assert result["eigenvalues"] == sorted(result["eigenvalues"], reverse=True)
    with pytest.raises(ConfigurationError):
        m.cluster_pca([[1.0, 2.0]], verbose=False)
