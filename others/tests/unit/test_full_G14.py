"""Contract tests for the G14 DendroPy, toytree and pyvolve phylogeny wrappers.

Expected values are hand calculations where they are short (for example the four-sequence
alignment below has 3 segregating sites), and direct upstream calls elsewhere.
"""

from __future__ import annotations

import inspect
import math
import random
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import dendropy_functions as dp
from dnakit.phylogeny import pyvolve_functions as pv
from dnakit.phylogeny import toytree_functions as tt

# Columns (1-based): 1,2,3,6,7,9,10 are monomorphic; 4, 5 and 8 are segregating.
ALIGNMENT = {
    "A": "ACGTACGTAC",
    "B": "ACGTACGAAC",
    "C": "ACGTTCGAAC",
    "D": "ACGAACGAAC",
}
TREE_T1 = "((A:1,B:1):1,(C:2,D:2):0.5);"
TREE_T2 = "((A:1,C:1):1,(B:2,D:2):0.5);"
TREE_FIVE = "((A,B),(C,(D,E)));"
TREE_UNROOTED = "(A:1,B:1,(C:2,D:2):1.5);"


def _dendropy():
    return pytest.importorskip("dendropy")


def _toytree():
    return pytest.importorskip("toytree")


def _write(tmp_path: Path, text: str, name: str = "tree.nwk") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _upstream_pair(dendropy, text_a: str, text_b: str):
    namespace = dendropy.TaxonNamespace()
    first = dendropy.Tree.get(data=text_a, schema="newick", taxon_namespace=namespace)
    second = dendropy.Tree.get(data=text_b, schema="newick", taxon_namespace=namespace)
    return first, second


# ---------------------------------------------------------------- DendroPy


def test_dendropy_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "dendropy", None)
    monkeypatch.setitem(sys.modules, "dendropy.calculate.popgenstat", None)
    with pytest.raises(BackendUnavailableError, match="bio-dendropy"):
        dp.count_segregating_sites(ALIGNMENT)


def test_segregating_sites_matches_hand_count_and_upstream() -> None:
    dendropy = _dendropy()
    from dendropy.calculate import popgenstat

    result = dp.count_segregating_sites(ALIGNMENT)
    assert result["status"] == "completed"
    assert result["operation"] == "count_segregating_sites"
    assert result["result"]["segregating_sites"] == 3  # columns 4, 5 and 8
    assert result["result"]["sequence_count"] == 4
    assert result["result"]["site_count"] == 10
    upstream = popgenstat.num_segregating_sites(dendropy.DnaCharacterMatrix.from_dict(ALIGNMENT))
    assert result["result"]["segregating_sites"] == upstream
    assert result["provenance"]["package"].startswith("dendropy==")


def test_mean_pairwise_differences_hand_value() -> None:
    # Pairwise differences: AB 1, AC 2, AD 2, BC 1, BD 1, CD 2 -> sum 9 over 6 pairs.
    result = dp.calculate_mean_pairwise_differences(ALIGNMENT)
    assert result["result"]["mean_pairwise_differences"] == pytest.approx(1.5)


def test_nucleotide_diversity_hand_value() -> None:
    # pi = k / L = 1.5 / 10 sites.
    result = dp.calculate_nucleotide_diversity(ALIGNMENT)
    assert result["result"]["nucleotide_diversity"] == pytest.approx(0.15)


def test_tajimas_d_matches_formula_and_is_none_without_segregating_sites() -> None:
    # Tajima (1989) with n=4, S=3, k=1.5, computed independently here.
    n, seg, k = 4, 3, 1.5
    a1 = sum(1 / i for i in range(1, n))
    a2 = sum(1 / i**2 for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / a1**2
    e1, e2 = c1 / a1, c2 / (a1**2 + a2)
    expected = (k - seg / a1) / math.sqrt(e1 * seg + e2 * seg * (seg - 1))
    result = dp.calculate_tajimas_d(ALIGNMENT)
    assert result["result"]["tajimas_d"] == pytest.approx(expected, abs=1e-9)
    assert result["result"]["segregating_sites"] == 3

    flat = dp.calculate_tajimas_d({"a": "AAA", "b": "AAA"})
    assert flat["result"]["tajimas_d"] is None
    assert flat["result"]["segregating_sites"] == 0


def test_wattersons_theta_hand_value() -> None:
    # theta_W = S / a_n with a_n = 1 + 1/2 + 1/3 = 11/6, so 3 / (11/6) = 18/11.
    result = dp.calculate_wattersons_theta(ALIGNMENT)
    assert result["result"]["wattersons_theta"] == pytest.approx(18 / 11)


def test_derived_state_matrix_rows_follow_hand_states_and_input_order() -> None:
    # Ancestral sequence A; derived (1) where a sample differs from A.
    expected = [
        {"id": "B", "states": "0000000100"},
        {"id": "C", "states": "0000100100"},
        {"id": "D", "states": "0001000100"},
    ]
    result = dp.build_derived_state_matrix(ALIGNMENT, ancestral_id="A")
    assert result["result"]["rows"] == expected
    assert result["result"]["ancestral_id"] == "A"

    shuffled = dict(random.Random(7).sample(list(ALIGNMENT.items()), k=4))
    assert dp.build_derived_state_matrix(shuffled, ancestral_id="A") == result


def test_unfolded_site_frequency_spectrum_hand_counts() -> None:
    # Derived-allele counts among B, C, D: 7 sites with 0, 2 with 1, 0 with 2, 1 with 3.
    result = dp.calculate_unfolded_site_frequency_spectrum(ALIGNMENT, ancestral_id="A")
    assert result["result"]["site_frequency_spectrum"] == [7, 2, 0, 1]
    assert result["result"]["ingroup_count"] == 3


def test_popgen_validation_errors() -> None:
    with pytest.raises(ConfigurationError, match="At least two"):
        dp.count_segregating_sites({"only": "ACGT"})
    with pytest.raises(ConfigurationError, match="same length"):
        dp.count_segregating_sites({"a": "ACGT", "b": "ACG"})
    with pytest.raises(ConfigurationError, match="unsupported DNA"):
        dp.count_segregating_sites({"a": "ACGT", "b": "ACGX"})
    with pytest.raises(ConfigurationError, match="ancestral_id"):
        dp.build_derived_state_matrix(ALIGNMENT, ancestral_id="Z")
    with pytest.raises(ConfigurationError, match="ignore_uncertain"):
        dp.count_segregating_sites(ALIGNMENT, ignore_uncertain="yes")
    with pytest.raises(ConfigurationError, match="max_sequences"):
        dp.count_segregating_sites(ALIGNMENT, max_sequences=3)


def test_weighted_rf_matches_upstream_and_is_zero_for_identical_trees() -> None:
    dendropy = _dendropy()
    from dendropy.calculate import treecompare

    first, second = _upstream_pair(dendropy, TREE_T1, TREE_T2)
    first.encode_bipartitions()
    second.encode_bipartitions()
    upstream = treecompare.weighted_robinson_foulds_distance(first, second)
    result = dp.calculate_weighted_robinson_foulds_distance(TREE_T1, TREE_T2)
    assert result["result"]["distance"] == pytest.approx(upstream)
    assert result["result"]["tip_count"] == 4

    same = dp.calculate_weighted_robinson_foulds_distance(TREE_T1, TREE_T1)
    assert same["result"]["distance"] == pytest.approx(0.0)


def test_false_positive_and_negative_splits_hand_counts(tmp_path: Path) -> None:
    # Splits: T1 has {A,B}|{C,D}; T2 has {A,C}|{B,D}. Each tree has one split the other lacks.
    reference = _write(tmp_path, TREE_T1, "ref.nwk")
    result = dp.count_false_positive_and_negative_splits(reference, TREE_T2)
    assert result["result"]["false_positives"] == 1
    assert result["result"]["false_negatives"] == 1
    same = dp.count_false_positive_and_negative_splits(TREE_T1, TREE_T1)
    assert same["result"]["false_positives"] == 0
    assert same["result"]["false_negatives"] == 0


def test_colless_imbalance_hand_values() -> None:
    # Raw Colless I: root 2-vs-3 gives 1, the (C,(D,E)) node 1-vs-2 gives 1, so I = 2.
    # Maximum for n=5 is (n-1)(n-2)/2 = 6, so the normalized value is 2/6.
    raw = dp.calculate_colless_tree_imbalance(TREE_FIVE, normalize=None)
    assert raw["result"]["colless_imbalance"] == pytest.approx(2.0)
    assert raw["result"]["normalize"] == "none"
    normalized = dp.calculate_colless_tree_imbalance(TREE_FIVE)
    assert normalized["result"]["colless_imbalance"] == pytest.approx(2 / 6)


def test_colless_rejects_non_bifurcating_tree() -> None:
    with pytest.raises(ConfigurationError, match="bifurcating"):
        dp.calculate_colless_tree_imbalance("(A,B,C,D);")


def test_sackin_index_hand_values() -> None:
    # In ((A,B),(C,D)) every tip has two ancestors (its cherry and the root): 4 * 2 = 8.
    raw = dp.calculate_sackin_index("((A,B),(C,D));", normalize=None)
    assert raw["result"]["sackin_index"] == pytest.approx(8.0)
    per_tip = dp.calculate_sackin_index("((A,B),(C,D));")
    assert per_tip["result"]["sackin_index"] == pytest.approx(2.0)
    with pytest.raises(ConfigurationError, match="normalize"):
        dp.calculate_sackin_index("((A,B),(C,D));", normalize="max")


def test_b1_index_hand_value() -> None:
    # Non-root internal nodes: (A,B) and (D,E) contribute 1 each; (C,(D,E)) contributes 1/2.
    result = dp.calculate_b1_index(TREE_FIVE)
    assert result["result"]["b1"] == pytest.approx(2.5)


def test_pybus_harvey_gamma_matches_upstream_and_rejects_non_ultrametric() -> None:
    dendropy = _dendropy()
    from dendropy.calculate import treemeasure

    ultrametric = "((A:1,B:1):1,(C:1,D:1):1);"
    upstream = treemeasure.pybus_harvey_gamma(dendropy.Tree.get(data=ultrametric, schema="newick"))
    result = dp.calculate_pybus_harvey_gamma(ultrametric)
    assert result["result"]["gamma"] == pytest.approx(upstream)
    with pytest.raises(ConfigurationError, match="ultrametric"):
        dp.calculate_pybus_harvey_gamma(TREE_T1)
    with pytest.raises(ConfigurationError, match="precision"):
        dp.calculate_pybus_harvey_gamma(ultrametric, precision=0)


def test_treeness_hand_value() -> None:
    # Internal branches 1 + 0.5 = 1.5 of total length 1 + 1 + 2 + 2 + 1 + 0.5 = 7.5.
    result = dp.calculate_treeness(TREE_T1)
    assert result["result"]["treeness"] == pytest.approx(1.5 / 7.5)


def test_parsimony_score_hand_value_and_id_check() -> None:
    # Each of the three informative columns needs exactly one change on ((A,B),(C,D)).
    chars = {"A": "AAAC", "B": "AAGC", "C": "TAGC", "D": "TCGC"}
    result = dp.calculate_parsimony_score("((A,B),(C,D));", chars)
    assert result["result"]["parsimony_score"] == 3
    with pytest.raises(ConfigurationError, match="match the tree tips"):
        dp.calculate_parsimony_score("((A,B),(C,D));", {"A": "AA", "B": "AA", "C": "AA"})


def test_reconciliation_discordance_matches_upstream_and_is_zero_for_same_tree() -> None:
    dendropy = _dendropy()
    from dendropy.model import reconcile

    species = "((A,B),(C,(D,E)));"
    gene = "((A,C),(B,(D,E)));"
    gene_tree, species_tree = _upstream_pair(dendropy, gene, species)
    gene_tree.is_rooted = True
    species_tree.is_rooted = True
    gene_tree.encode_bipartitions()
    species_tree.encode_bipartitions()
    upstream = reconcile.reconciliation_discordance(gene_tree, species_tree)
    result = dp.calculate_reconciliation_discordance(gene, species)
    assert result["result"]["duplications"] == upstream
    assert dp.calculate_reconciliation_discordance(species, species)["result"]["duplications"] == 0


def test_monophyletic_partition_discordance_hand_and_upstream() -> None:
    dendropy = _dendropy()
    from dendropy.model import reconcile

    groups = [["A", "B"], ["C"], ["D", "E"]]
    # Species tree already has A,B and D,E as clades, so forcing those groups costs nothing.
    assert (
        dp.calculate_monophyletic_partition_discordance(TREE_FIVE, groups)["result"][
            "deep_coalescences"
        ]
        == 0
    )

    gene = "((A,C),(B,(D,E)));"
    namespace = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(data=gene, schema="newick", taxon_namespace=namespace)
    partition = dendropy.TaxonNamespacePartition(
        taxon_namespace=namespace,
        membership_lists=[[namespace.get_taxon(label) for label in group] for group in groups],
    )
    tree.encode_bipartitions()
    upstream = reconcile.monophyletic_partition_discordance(tree, partition)
    result = dp.calculate_monophyletic_partition_discordance(gene, groups)
    assert result["result"]["deep_coalescences"] == upstream == 1
    assert result["result"]["group_count"] == 3


def test_tree_input_from_file_matches_text(tmp_path: Path) -> None:
    path = _write(tmp_path, TREE_FIVE)
    from_file = dp.calculate_b1_index(path)
    from_text = dp.calculate_b1_index(TREE_FIVE)
    assert from_file["result"] == from_text["result"]


def test_tree_input_validation(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="exactly one tree"):
        dp.calculate_treeness("((A,B),C);((A,C),B);")
    with pytest.raises(ConfigurationError, match="max_tips"):
        dp.calculate_treeness(TREE_FIVE, max_tips=4)
    with pytest.raises(ConfigurationError, match="does not exist"):
        dp.calculate_treeness(tmp_path / "missing.nwk")
    with pytest.raises(ConfigurationError, match="same tip labels"):
        dp.calculate_weighted_robinson_foulds_distance(TREE_T1, "((A:1,B:1):1,(C:2,E:2):0.5);")
    with pytest.raises(ConfigurationError, match="partition"):
        dp.calculate_monophyletic_partition_discordance(TREE_FIVE, [["A", "B"], ["C"]])


def test_dendropy_verbose_is_keyword_only_and_results_do_not_depend_on_it() -> None:
    result = dp.count_segregating_sites(ALIGNMENT, verbose=False)
    assert result["result"]["segregating_sites"] == 3
    for name in dp.__all__:
        parameter = inspect.signature(getattr(dp, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


# ---------------------------------------------------------------- toytree


def test_toytree_missing_package_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "toytree", None)
    with pytest.raises(BackendUnavailableError, match="bio-toytree"):
        tt.unroot_tree(TREE_T1)


def test_midpoint_rooting_matches_hand_newick() -> None:
    # Longest path C..A = 2 + 1.5 + 1 = 4.5; the midpoint is 2.25 from C, which is 0.25 past
    # the internal node of (C,D) on the 1.5 edge. Root-to-A and root-to-C are both 2.25.
    _toytree()
    result = tt.root_tree_on_midpoint(TREE_UNROOTED)
    assert result["result"]["newick"] == "((C:2,D:2):0.25,(A:1,B:1):1.25);"
    assert result["result"]["tips"] == ["A", "B", "C", "D"]


def test_balanced_and_mad_rooting_match_upstream() -> None:
    toytree = _toytree()
    from toytree import mod

    balanced = tt.root_tree_on_balanced_midpoint(TREE_UNROOTED, tolerance=1e-9)
    upstream = mod.root_on_balanced_midpoint(toytree.tree(TREE_UNROOTED), tolerance=1e-9)
    assert balanced["result"]["newick"] == upstream.write()

    mad = tt.root_tree_on_minimal_ancestor_deviation(TREE_UNROOTED)
    upstream_mad = mod.root_on_minimal_ancestor_deviation(toytree.tree(TREE_UNROOTED))
    assert mad["result"]["newick"] == upstream_mad.write()


def test_unroot_merges_the_root_edges_by_hand() -> None:
    # Root edges 1 (to A,B) and 0.5 (to C,D): the C,D node keeps 0.5 + 1 = 1.5 above it.
    result = tt.unroot_tree(TREE_T1)
    assert result["result"]["newick"] == "(A:1,B:1,(C:2,D:2):1.5);"


def test_ladderize_matches_upstream_and_direction_reverses_order() -> None:
    toytree = _toytree()
    from toytree import mod

    text = "((A:1,B:1):1,(C:2,(D:2,E:2):1):0.5);"
    default = tt.ladderize_tree(text)
    assert default["result"]["newick"] == mod.ladderize(toytree.tree(text)).write()
    reversed_order = tt.ladderize_tree(text, direction=True)
    assert reversed_order["result"]["newick"] == "(((D:2,E:2):1,C:2):0.5,(A:1,B:1):1);"
    assert reversed_order["result"]["newick"] != default["result"]["newick"]


def test_extract_subtree_hand_clade_and_missing_tip_error() -> None:
    result = tt.extract_subtree(TREE_T1, ["C", "D"])
    assert result["result"]["newick"] == "(C:2,D:2):0.5;"
    assert result["result"]["tips"] == ["C", "D"]
    with pytest.raises(ConfigurationError, match="not in the tree"):
        tt.extract_subtree(TREE_T1, ["C", "Z"])


def test_collapse_weak_nodes_hand_by_support_label() -> None:
    # The (C,D) node has support 30 < 50, so it collapses and its children gain 0.5 each.
    text = "((A:1,B:1)90:1,(C:2,D:2)30:0.5);"
    result = tt.collapse_weak_nodes(text, min_support=50)
    assert result["result"]["newick"] == "((A:1,B:1)90:1,C:2.5,D:2.5);"
    untouched = tt.collapse_weak_nodes(text, min_support=0)
    assert untouched["result"]["newick"] == "((A:1,B:1)90:1,(C:2,D:2)30:0.5);"


def test_remove_unary_nodes_hand_length_merge() -> None:
    # A hangs from a single-child node of length 1; merging gives A a length of 1 + 1 = 2.
    result = tt.remove_unary_nodes("((A:1):1,(B:2,C:2):1);")
    assert result["result"]["newick"] == "((B:2,C:2):1,A:2);"


def test_scale_tree_to_root_height_hand_factor() -> None:
    # Longest root-to-tip is 2.5 (C), so the factor for height 2 is 0.8 without the stem.
    result = tt.scale_tree_to_root_height(TREE_T1, height=2.0)
    assert result["result"]["newick"] == "((A:0.8,B:0.8):0.8,(C:1.6,D:1.6):0.4);"


def test_tree_vcv_and_correlation_matrices_hand_values() -> None:
    vcv = tt.calculate_tree_vcv_matrix(TREE_T1)
    assert vcv["result"]["tips"] == ["A", "B", "C", "D"]
    assert vcv["result"]["matrix"] == [
        [2.0, 1.0, 0.0, 0.0],
        [1.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 2.5, 0.5],
        [0.0, 0.0, 0.5, 2.5],
    ]
    corr = tt.calculate_tree_correlation_matrix(TREE_T1)["result"]["matrix"]
    assert corr[0][1] == pytest.approx(0.5)  # 1 / sqrt(2 * 2)
    assert corr[2][3] == pytest.approx(0.2)  # 0.5 / sqrt(2.5 * 2.5)


@pytest.mark.parametrize(
    ("wrapper", "upstream_name"),
    [
        (tt.calculate_information_rf_distance, "get_treedist_rfi"),
        (tt.calculate_matching_split_distance, "get_treedist_rfg_ms"),
        (tt.calculate_matching_split_information_distance, "get_treedist_rfg_msi"),
        (tt.calculate_shared_phylogenetic_information_distance, "get_treedist_rfg_spi"),
        (tt.calculate_mutual_clustering_information_distance, "get_treedist_rfg_mci"),
    ],
)
def test_tree_distance_wrappers_match_upstream(wrapper, upstream_name: str) -> None:
    toytree = _toytree()
    from toytree import distance

    upstream = getattr(distance, upstream_name)(
        toytree.tree(TREE_T1), toytree.tree(TREE_T2), normalize=False
    )
    result = wrapper(TREE_T1, TREE_T2)
    assert result["result"]["distance"] == pytest.approx(float(upstream))
    assert result["result"]["normalize"] is False


@pytest.mark.parametrize(
    "wrapper",
    [
        tt.calculate_information_rf_distance,
        tt.calculate_matching_split_distance,
        tt.calculate_matching_split_information_distance,
        tt.calculate_mutual_clustering_information_distance,
    ],
)
def test_tree_distance_is_zero_for_identical_trees(wrapper) -> None:
    assert wrapper(TREE_T1, TREE_T1)["result"]["distance"] == pytest.approx(0.0, abs=1e-9)


def test_shared_phylogenetic_information_rejects_identical_split_sets() -> None:
    # toytree's SPI similarity is infinite for identical splits and the assignment fails.
    with pytest.raises(ConfigurationError, match="Identical split"):
        tt.calculate_shared_phylogenetic_information_distance(TREE_T1, TREE_T1)


def test_matching_split_distance_hand_value() -> None:
    # The two trees share no nontrivial split, so each contributes one mismatched split: 2.
    assert tt.calculate_matching_split_distance(TREE_T1, TREE_T2)["result"]["distance"] == 2.0


def test_tree_distance_requires_identical_tips() -> None:
    with pytest.raises(ConfigurationError, match="same tip labels"):
        tt.calculate_information_rf_distance(TREE_T1, "((A:1,C:1):1,(B:2,E:2):0.5);")


def test_list_bipartitions_hand_splits() -> None:
    result = tt.list_tree_bipartitions(TREE_T1)
    assert result["result"]["splits"] == [[["A", "B"], ["C", "D"]]]
    with_singletons = tt.list_tree_bipartitions(TREE_T1, include_singletons=True)
    assert with_singletons["result"]["split_count"] == 5
    assert [["A"], ["B", "C", "D"]] in with_singletons["result"]["splits"]


def test_relative_evolutionary_divergence_matches_upstream() -> None:
    toytree = _toytree()
    from toytree import pcm

    upstream = pcm.get_relative_evolutionary_divergence(toytree.tree(TREE_T1))
    rows = tt.calculate_relative_evolutionary_divergence(TREE_T1)["result"]["nodes"]
    assert len(rows) == len(upstream)
    for row in rows:
        assert row["red"] == pytest.approx(float(upstream[row["node_index"]]))
    assert [row["node_index"] for row in rows] == sorted(row["node_index"] for row in rows)


def test_equal_splits_matches_upstream_for_tips_only() -> None:
    toytree = _toytree()
    from toytree import pcm

    series = pcm.get_equal_splits(toytree.tree(TREE_T1))
    tips = tt.calculate_equal_splits(TREE_T1)["result"]["tips"]
    assert [row["tip"] for row in tips] == ["A", "B", "C", "D"]
    for row in tips:
        assert row["equal_splits"] == pytest.approx(float(series[row["tip"]]))


def test_toytree_validation_errors() -> None:
    with pytest.raises(ConfigurationError, match="min_support"):
        tt.collapse_weak_nodes(TREE_T1, min_support=-1)
    with pytest.raises(ConfigurationError, match="at least two tip labels"):
        tt.extract_subtree(TREE_T1, ["A"])
    with pytest.raises(ConfigurationError, match="max_tips"):
        tt.unroot_tree(TREE_FIVE, max_tips=4)
    with pytest.raises(ConfigurationError, match="positive"):
        tt.scale_tree_to_root_height(TREE_T1, height=0)


def test_toytree_verbose_is_keyword_only_and_results_do_not_depend_on_it() -> None:
    quiet = tt.unroot_tree(TREE_T1, verbose=False)
    assert quiet["result"]["newick"] == "(A:1,B:1,(C:2,D:2):1.5);"
    for name in tt.__all__:
        parameter = inspect.signature(getattr(tt, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


# ---------------------------------------------------------------- pyvolve


def test_pyvolve_module_exports_no_new_wrappers() -> None:
    pyvolve = pytest.importorskip("pyvolve")
    assert pv.__all__ == []
    # The inventory's only tree reader is still available for the existing simulate_* backends.
    assert callable(pyvolve.read_tree)
