"""Tests for dnakit.phylogeny.skbio_p20_tree_functions (scikit-bio tree and taxonomy operations)."""

from __future__ import annotations

import pytest

pytest.importorskip("skbio")

from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny import skbio_p20_tree_functions as tree_fn

TREE = "((a:1,b:2)90:1,(c:3,d:1)0.5:2,e:1);"


def test_shuffle_is_reproducible_and_keeps_tips():
    first = tree_fn.shuffle_tree_tips(TREE, seed=3, n=2, verbose=False)["trees"]
    again = tree_fn.shuffle_tree_tips(TREE, seed=3, n=2, verbose=False)["trees"]
    assert first == again and len(first) == 2
    assert all(sorted(c for c in t if c in "abcde") == list("abcde") for t in first)
    with pytest.raises(ConfigurationError):
        tree_fn.shuffle_tree_tips("not a tree", verbose=False)


def test_reroot_methods_agree():
    expected = "((d:1.0,((a:1.0,b:2.0)90:1.0,e:1.0):2.0)0.5:3.0)c;"
    for method in ("root_at", "unrooted_copy", "unrooted_move"):
        assert tree_fn.reroot_tree(TREE, "c", method=method, verbose=False)["newick"] == expected
    above = tree_fn.reroot_tree(TREE, ["c", "d"], above=True, verbose=False)
    assert above["tip_count"] == 5
    with pytest.raises(ConfigurationError):
        tree_fn.reroot_tree(TREE, "zzz", verbose=False)


def test_subsets_supports_counts_and_diameter():
    clades = tree_fn.tree_tip_subsets(TREE, verbose=False)["clades"]
    assert [c["tips"] for c in clades] == [["a", "b"], ["c", "d"]]
    assert tree_fn.tree_tip_subsets(TREE, node=["c", "d"], verbose=False)["tips"] == ["c", "d"]
    lengths = tree_fn.tree_tip_subsets(TREE, include_tips=True, map_to_length=True, verbose=False)[
        "clades"
    ]
    assert {"tips": ["c"], "length": 3.0} in lengths
    supports = tree_fn.extract_tree_supports(TREE, verbose=False)["supports"]
    assert supports == [{"clade": "a,b", "support": 90.0}, {"clade": "c,d", "support": 0.5}]
    counts = tree_fn.propagate_tip_counts(TREE, {"a": 1, "c": 2}, verbose=False)["nodes"]
    assert {"node": "a,b,c,d,e", "count": 3} in counts and {"node": "c,d", "count": 2} in counts
    diameter = tree_fn.tree_diameter(TREE, verbose=False)
    assert diameter["distance"] == 8.0 and diameter["tips"] == ["b", "c"]
    assert tree_fn.tree_diameter(TREE, use_length=False, verbose=False)["distance"] == 4.0


def test_linkage_and_taxonomy_trees(tmp_path):
    linkage = [[0, 1, 1.0, 2], [2, 3, 1.0, 2], [4, 5, 5.0, 4]]
    built = tree_fn.tree_from_linkage_matrix(linkage, list("abcd"), verbose=False)
    assert built["newick"] == "((a:0.5,b:0.5):2.0,(c:0.5,d:0.5):2.0);"
    lineages = {"g1": ["k__B", "p__X"], "g3": ["k__A"]}
    taxonomy = tree_fn.build_taxonomy_tree(lineages, verbose=False)
    exported = tree_fn.export_tree_lineages(taxonomy["newick"], verbose=False)["lineages"]
    assert exported == [
        {"id": "g1", "lineage": ["k__B", "p__X"]},
        {"id": "g3", "lineage": ["k__A"]},
    ]
    assert (
        "'k__B'"
        not in tree_fn.build_taxonomy_tree(lineages, extract_rank=True, verbose=False)["newick"]
    )
    row = "\t|\t".join(["{}", "{}", "{}", "", "8", "0", "1", "0", "0", "0", "0", "0", ""]) + "\t|\n"
    (tmp_path / "nodes.dmp").write_text(
        row.format(1, 1, "no rank") + row.format(2, 1, "superkingdom")
    )
    (tmp_path / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
    )
    result = tree_fn.build_tree_from_taxdump(
        str(tmp_path / "nodes.dmp"), names_path=str(tmp_path / "names.dmp"), verbose=False
    )
    assert result["newick"] == "(Bacteria)root;"
    with pytest.raises(ConfigurationError):
        tree_fn.tree_from_linkage_matrix([[0, 1, 1.0]], ["a", "b"], verbose=False)
