"""Contract tests for dendropy, toytree and pyvolve phylogeny backends."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny.backends import (
    calculate_branch_length_distance,
    calculate_node_coordinates,
    calculate_patristic_distances,
    calculate_quartet_distance,
    read_tree_summary,
    root_tree_on_outgroup,
    simulate_nucleotide_alignment,
    simulate_nucleotide_sequences,
    simulate_substitution_counts,
)

TREE_4 = "((a:1,b:2):0.5,(c:0.25,d:3):1.5);\n"
NEXUS_2 = """#NEXUS
begin taxa;
  dimensions ntax=3;
  taxlabels X Y Z;
end;
begin trees;
  tree one = [&R] ((X:1,Y:2):3,Z:4);
  tree two = ((X:1,Z:2):3,Y:4);
end;
"""


def _write(tmp_path: Path, text: str, name: str = "tree.nwk") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


# ---------------------------------------------------------------- dendropy


def test_read_tree_summary_reports_hand_computed_tips_and_lengths(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    path = _write(tmp_path, TREE_4)
    result = read_tree_summary(path)
    assert result["status"] == "completed"
    assert result["operation"] == "read_tree_summary"
    body = result["result"]
    # Tips sorted by label with their pendant lengths; total = 1 + 2 + 0.5 + 0.25 + 3 + 1.5.
    assert [(t["label"], t["branch_length"]) for t in body["tips"]] == [
        ("a", 1.0),
        ("b", 2.0),
        ("c", 0.25),
        ("d", 3.0),
    ]
    assert body["branch_lengths"] == [0.5, 1.0, 2.0, 1.5, 0.25, 3.0]
    assert body["total_branch_length"] == pytest.approx(8.25)
    assert body["missing_branch_lengths"] == 0
    assert body["rooting"] == "undefined"
    assert body["tree_count"] == 1
    assert result["provenance"]["package"].startswith("dendropy==")


def test_read_tree_summary_selects_nexus_tree_and_keeps_underscores(tmp_path: Path) -> None:
    dendropy = pytest.importorskip("dendropy")
    nexus = _write(tmp_path, NEXUS_2, "trees.nex")
    first = read_tree_summary(nexus, schema="nexus", tree_index=0)["result"]
    second = read_tree_summary(nexus, schema="nexus", tree_index=1)["result"]
    assert first["tree_count"] == second["tree_count"] == 2
    assert first["rooting"] == "rooted"
    assert second["rooting"] == "undefined"
    assert [t["label"] for t in first["tips"]] == ["X", "Y", "Z"]
    # Every tree in the Nexus file has the same total length: 1 + 2 + 3 + 4.
    assert first["total_branch_length"] == pytest.approx(10.0)
    assert second["total_branch_length"] == pytest.approx(
        dendropy.TreeList.get(path=str(nexus), schema="nexus")[1].length()
    )

    underscored = _write(tmp_path, "((Homo_sapiens:1,Pan:2):1,(Gorilla:3,Pongo:4):1);\n", "u.nwk")
    labels = [t["label"] for t in read_tree_summary(underscored)["result"]["tips"]]
    assert labels == ["Gorilla", "Homo_sapiens", "Pan", "Pongo"]


def test_read_tree_summary_matches_direct_dendropy(tmp_path: Path) -> None:
    dendropy = pytest.importorskip("dendropy")
    path = _write(tmp_path, TREE_4)
    direct = dendropy.Tree.get(path=str(path), schema="newick")
    summary = read_tree_summary(path)["result"]
    assert summary["total_branch_length"] == pytest.approx(direct.length())
    assert summary["tip_count"] == len(list(direct.leaf_node_iter())) == 4


def test_patristic_matrix_matches_hand_sums(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    path = _write(tmp_path, TREE_4)
    body = calculate_patristic_distances(path)["result"]
    assert body["taxa"] == ["a", "b", "c", "d"]
    expected = [
        [0.0, 3.0, 3.25, 6.0],  # a-b = 1+2; a-c = 1+0.5+1.5+0.25; a-d = 1+0.5+1.5+3
        [3.0, 0.0, 4.25, 7.0],  # b-c = 2+0.5+1.5+0.25; b-d = 2+0.5+1.5+3
        [3.25, 4.25, 0.0, 3.25],  # c-d = 0.25+3
        [6.0, 7.0, 3.25, 0.0],
    ]
    for row, expected_row in zip(body["matrix"], expected, strict=True):
        assert row == pytest.approx(expected_row)


def test_patristic_matrix_matches_root_distance_identity(tmp_path: Path) -> None:
    dendropy = pytest.importorskip("dendropy")
    newick = "(((a:0.3,b:0.1):0.2,c:0.5):0.4,(d:0.25,(e:0.15,f:0.6):0.35):0.1);\n"
    path = _write(tmp_path, newick, "six.nwk")
    tree = dendropy.Tree.get(path=str(path), schema="newick")

    def ancestors(node):
        chain = []
        while node is not None:
            chain.append(node)
            node = node.parent_node
        return chain

    def root_distance(node):
        return sum(item.edge.length for item in ancestors(node)[:-1])

    leaves = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter()}
    body = calculate_patristic_distances(path, max_taxa=10)["result"]
    labels = body["taxa"]
    for i, first in enumerate(labels):
        for j, second in enumerate(labels):
            if i == j:
                continue
            # Path length via the MRCA: d(a,b) = depth(a) + depth(b) - 2 * depth(MRCA).
            second_chain = set(map(id, ancestors(leaves[second])))
            mrca = next(n for n in ancestors(leaves[first]) if id(n) in second_chain)
            expected = (
                root_distance(leaves[first])
                + root_distance(leaves[second])
                - 2 * root_distance(mrca)
            )
            assert body["matrix"][i][j] == pytest.approx(expected)


def test_patristic_rejects_missing_branch_lengths(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    path = _write(tmp_path, "((a,b):1,c);\n", "missing.nwk")
    with pytest.raises(ConfigurationError, match="branch length"):
        calculate_patristic_distances(path)


def test_branch_length_distance_matches_hand_value_and_dendropy(tmp_path: Path) -> None:
    dendropy = pytest.importorskip("dendropy")
    treecompare = pytest.importorskip("dendropy.calculate.treecompare")
    first = _write(tmp_path, "((A:1,B:2):1,(C:3,D:4):1);\n", "x.nwk")
    second = _write(tmp_path, "((A:1,B:1):1,(C:2,D:4):1);\n", "y.nwk")
    result = calculate_branch_length_distance(first, second)["result"]
    # Only B (2 vs 1) and C (3 vs 2) differ, so the distance is sqrt(1 + 1).
    assert result["distance"] == pytest.approx(math.sqrt(2.0))
    assert result["tip_count"] == 4
    namespace = dendropy.TaxonNamespace()
    tree_a = dendropy.Tree.get(path=str(first), schema="newick", taxon_namespace=namespace)
    tree_b = dendropy.Tree.get(path=str(second), schema="newick", taxon_namespace=namespace)
    assert result["distance"] == pytest.approx(treecompare.euclidean_distance(tree_a, tree_b))


def test_branch_length_distance_ignores_root_placement(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    # Same unrooted splits: the two root-adjacent edges sum to 4 on both trees.
    first = _write(tmp_path, "((A:1,B:1):1,(C:1,D:1):3);\n", "a.nwk")
    second = _write(tmp_path, "((A:1,B:1):2,(C:1,D:1):2);\n", "b.nwk")
    assert calculate_branch_length_distance(first, second)["result"]["distance"] == 0.0


def test_branch_length_distance_rejects_different_tips(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    first = _write(tmp_path, "((A:1,B:1):1,(C:1,D:1):1);\n", "a.nwk")
    second = _write(tmp_path, "((A:1,B:1):1,(C:1,E:1):1);\n", "b.nwk")
    with pytest.raises(ConfigurationError, match="same tip labels"):
        calculate_branch_length_distance(first, second)


def test_tree_readers_reject_bad_input(tmp_path: Path) -> None:
    pytest.importorskip("dendropy")
    path = _write(tmp_path, TREE_4)
    with pytest.raises(ConfigurationError, match="schema"):
        read_tree_summary(path, schema="phyloxml")
    with pytest.raises(ConfigurationError, match="tree_index"):
        read_tree_summary(path, tree_index=True)
    with pytest.raises(ConfigurationError, match="out of range"):
        read_tree_summary(path, tree_index=1)
    with pytest.raises(ConfigurationError, match="max_tips"):
        read_tree_summary(path, max_tips=3)
    with pytest.raises(ConfigurationError, match="max_file_bytes"):
        read_tree_summary(path, max_file_bytes=5)
    with pytest.raises(ConfigurationError, match="does not exist"):
        read_tree_summary(tmp_path / "missing.nwk")
    with pytest.raises(ConfigurationError, match="could not parse"):
        read_tree_summary(_write(tmp_path, "((a:1,b:1);\n", "bad.nwk"))
    with pytest.raises(ConfigurationError, match="Duplicate taxon labels"):
        read_tree_summary(_write(tmp_path, "((a:1,a:1):1,b:1);\n", "dup.nwk"))


def test_dendropy_missing_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, TREE_4)
    monkeypatch.setitem(sys.modules, "dendropy", None)
    with pytest.raises(BackendUnavailableError, match="bio-dendropy"):
        read_tree_summary(path)


# ----------------------------------------------------------------- toytree

LINEAR_TIPS_AND_NODES = {
    # index: (x, y) from toytree's linear layout of TREE_4, derived by hand.
    0: (-3.0, 0.0),  # a: root distance 1.5 minus max root distance 4.5
    1: (-2.0, 1.0),  # b: 2.5 - 4.5
    2: (-2.75, 2.0),  # c: 1.75 - 4.5
    3: (0.0, 3.0),  # d: 4.5 - 4.5
    4: (-4.0, 0.5),  # X: 0.5 - 4.5, y is the mean of a and b
    5: (-3.0, 2.5),  # Y: 1.5 - 4.5, y is the mean of c and d
    6: (-4.5, 1.5),  # root: y is the mean of X and Y
}


def test_node_coordinates_match_hand_linear_layout(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    path = _write(tmp_path, TREE_4)
    result = calculate_node_coordinates(path, layout="r")["result"]
    assert result["tip_count"] == 4
    assert result["node_count"] == 7
    rows = {row["index"]: row for row in result["nodes"]}
    for index, (x_value, y_value) in LINEAR_TIPS_AND_NODES.items():
        assert rows[index]["x"] == pytest.approx(x_value)
        assert rows[index]["y"] == pytest.approx(y_value)
    assert [row["name"] for row in result["nodes"][:4]] == ["a", "b", "c", "d"]
    assert rows[6]["parent"] is None


def test_node_coordinates_match_direct_toytree_layout(tmp_path: Path) -> None:
    toytree = pytest.importorskip("toytree")
    layout_module = pytest.importorskip("toytree.layout")
    style_module = pytest.importorskip("toytree.style")
    path = _write(tmp_path, TREE_4)
    tree = toytree.tree(path)
    direct = layout_module.LinearLayout(tree, style_module.TreeStyle(layout="r")).coords
    rows = calculate_node_coordinates(path, layout="r")["result"]["nodes"]
    for row in rows:
        assert row["x"] == pytest.approx(float(direct[row["index"]][0]))
        assert row["y"] == pytest.approx(float(direct[row["index"]][1]))


def test_node_coordinates_without_edge_lengths_align_tips(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    path = _write(tmp_path, TREE_4)
    rows = calculate_node_coordinates(path, layout="r", use_edge_lengths=False)["result"]["nodes"]
    tip_x = [row["x"] for row in rows if row["is_tip"]]
    assert tip_x == [0.0, 0.0, 0.0, 0.0]


def test_circular_and_unrooted_layouts_return_finite_rows(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    path = _write(tmp_path, TREE_4)
    for layout in ("c", "unrooted"):
        result = calculate_node_coordinates(path, layout=layout)["result"]
        assert result["node_count"] == 7
        assert all(row["x"] is not None and row["y"] is not None for row in result["nodes"])


def test_quartet_distance_four_tips_hand_value(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    first = _write(tmp_path, "((a:1,b:1):1,(c:1,d:1):1);\n", "a.nwk")
    second = _write(tmp_path, "((a:1,c:1):1,(b:1,d:1):1);\n", "b.nwk")
    result = calculate_quartet_distance(first, second)["result"]
    # One 4-tip subset; ab|cd and ac|bd disagree.
    assert result["total_quartets"] == 1
    assert result["different_quartets"] == 1
    assert result["shared_quartets"] == 0
    assert result["different_fraction"] == pytest.approx(1.0)
    same = calculate_quartet_distance(first, first)["result"]
    assert same["shared_quartets"] == 1 and same["different_quartets"] == 0


def test_quartet_distance_five_tips_matches_hand_count(tmp_path: Path) -> None:
    toytree_distance = pytest.importorskip("toytree.distance")
    toytree = pytest.importorskip("toytree")
    # Quartets omitting e, d differ; omitting a, b, c agree -> 3 shared, 2 different.
    first = _write(tmp_path, "((a:1,b:1):1,c:1,(d:1,e:1):1);\n", "x.nwk")
    second = _write(tmp_path, "((a:1,c:1):1,b:1,(d:1,e:1):1);\n", "y.nwk")
    result = calculate_quartet_distance(first, second)["result"]
    assert (result["total_quartets"], result["shared_quartets"], result["different_quartets"]) == (
        5,
        3,
        2,
    )
    direct = toytree_distance.get_treedist_quartets(toytree.tree(first), toytree.tree(second))
    assert result["different_quartets"] == round(float(direct["D"]))
    assert result["shared_quartets"] == round(float(direct["S"]))


def test_root_tree_on_outgroup_hand_rooting(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    path = _write(tmp_path, "((a:1,b:1):1,(c:1,d:1):1);\n")
    target = tmp_path / "rooted.nwk"
    result = root_tree_on_outgroup(path, ["c"], target)["result"]
    # Root placed on the 1.0 edge above c: each side gets 0.5, and the old root merges into the
    # path to a and b, so that edge becomes 1 + 1 = 2. Total length stays 6.
    assert result["newick"] == "(c:0.5,(d:1,(a:1,b:1):2):0.5);\n"
    assert target.read_text(encoding="utf-8") == result["newick"]
    assert result["tip_count"] == 4


def test_root_tree_on_outgroup_preserves_total_length(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    dendropy = pytest.importorskip("dendropy")
    path = _write(tmp_path, TREE_4)
    result = root_tree_on_outgroup(path, "c", tmp_path / "c_rooted.nwk")["result"]
    assert result["newick"] == "(c:0.125,(d:3,(a:1,b:2):2):0.125);\n"
    rooted = dendropy.Tree.get(path=str(tmp_path / "c_rooted.nwk"), schema="newick")
    assert rooted.length() == pytest.approx(8.25)


def test_root_tree_on_outgroup_rejects_bad_outgroup(tmp_path: Path) -> None:
    pytest.importorskip("toytree")
    path = _write(tmp_path, TREE_4)
    with pytest.raises(ConfigurationError, match="not in the tree"):
        root_tree_on_outgroup(path, ["zzz"], tmp_path / "x.nwk")
    with pytest.raises(ConfigurationError, match="ingroup"):
        root_tree_on_outgroup(path, ["a", "b", "c", "d"], tmp_path / "x.nwk")
    with pytest.raises(ConfigurationError, match="could not root"):
        root_tree_on_outgroup(path, ["a", "c"], tmp_path / "x.nwk")
    with pytest.raises(ConfigurationError, match="unique"):
        root_tree_on_outgroup(path, ["a", "a"], tmp_path / "x.nwk")
    existing = tmp_path / "exists.nwk"
    existing.write_text("old\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        root_tree_on_outgroup(path, ["c"], existing)


def test_toytree_missing_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, TREE_4)
    monkeypatch.setitem(sys.modules, "toytree", None)
    with pytest.raises(BackendUnavailableError, match="bio-toytree"):
        calculate_node_coordinates(path)


# ----------------------------------------------------------------- pyvolve

PAIR_TREE = "(a:0.1,b:0.1);\n"


def test_alignment_is_reproducible_for_a_fixed_seed(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, TREE_4)
    first = simulate_nucleotide_alignment(path, tmp_path / "one.fasta", length=60, seed=7)
    simulate_nucleotide_alignment(path, tmp_path / "two.fasta", length=60, seed=7)
    simulate_nucleotide_alignment(path, tmp_path / "three.fasta", length=60, seed=8)
    one = _fasta(tmp_path / "one.fasta")
    assert one == _fasta(tmp_path / "two.fasta")
    assert one != _fasta(tmp_path / "three.fasta")
    assert first["result"]["taxa"] == ["a", "b", "c", "d"]
    assert first["result"]["sequence_count"] == 4
    assert first["result"]["site_count"] == 60
    assert set("".join(one.values())) <= set("ACGT")
    assert all(len(seq) == 60 for seq in one.values())


def test_file_and_in_memory_sequences_agree(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, TREE_4)
    simulate_nucleotide_alignment(path, tmp_path / "file.fasta", length=50, seed=3)
    in_memory = simulate_nucleotide_sequences(path, length=50, seed=3)["result"]["sequences"]
    assert in_memory == _fasta(tmp_path / "file.fasta")


def test_zero_length_branches_give_identical_tips(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, "(a:0,b:0);\n")
    sequences = simulate_nucleotide_sequences(path, length=500, seed=1)["result"]["sequences"]
    assert sequences["a"] == sequences["b"]


def test_base_composition_follows_state_freqs(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, "(a:0,b:0);\n")
    freqs = (0.4, 0.1, 0.1, 0.4)
    seq = simulate_nucleotide_sequences(path, length=20000, state_freqs=freqs, seed=3)["result"][
        "sequences"
    ]["a"]
    for base, expected in zip("ACGT", freqs, strict=True):
        assert seq.count(base) / len(seq) == pytest.approx(expected, abs=0.01)


def test_jukes_cantor_divergence_matches_formula(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, PAIR_TREE)
    sequences = simulate_nucleotide_sequences(path, length=20000, seed=11)["result"]["sequences"]
    observed = sum(p != q for p, q in zip(sequences["a"], sequences["b"], strict=True)) / 20000
    # Two tips, each 0.1 from the root; JC69 divergence over 0.2 expected substitutions per site.
    expected = 0.75 * (1 - math.exp(-4.0 / 3.0 * 0.2))
    assert observed == pytest.approx(expected, abs=0.01)


def test_kappa_sets_transition_fraction(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, "(a:0.05,b:0.05);\n")
    kappa = 10.0
    sequences = simulate_nucleotide_sequences(path, length=20000, kappa=kappa, seed=2)["result"][
        "sequences"
    ]
    differing = [(p, q) for p, q in zip(sequences["a"], sequences["b"], strict=True) if p != q]
    transitions = sum(1 for p, q in differing if {p, q} in ({"A", "G"}, {"C", "T"}))
    # Equal frequencies: a single substitution is a transition with probability kappa/(kappa+2).
    assert transitions / len(differing) == pytest.approx(kappa / (kappa + 2), abs=0.05)


def test_substitution_counts_match_expected_per_branch(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, PAIR_TREE)
    counts = simulate_substitution_counts(path, length=50000, seed=5)["result"]
    per_branch = counts["branch_substitution_counts"]
    # Expected substitutions per branch = branch length (0.1) times sites (50000) = 5000.
    assert per_branch["a"] == pytest.approx(5000, abs=300)
    assert per_branch["b"] == pytest.approx(5000, abs=300)
    assert counts["total_substitutions"] == sum(per_branch.values())


def test_zero_length_branches_have_no_substitutions(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, "(a:0,b:0);\n")
    counts = simulate_substitution_counts(path, length=300, seed=4)["result"]
    assert counts["total_substitutions"] == 0
    assert set(counts["branch_substitution_counts"].values()) == {0}


def test_pyvolve_rejects_bad_inputs(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, TREE_4)
    with pytest.raises(ConfigurationError, match="kappa"):
        simulate_nucleotide_sequences(path, kappa=0)
    with pytest.raises(ConfigurationError, match="sum to 1"):
        simulate_nucleotide_sequences(path, state_freqs=(0.5, 0.5, 0.5, 0.5))
    with pytest.raises(ConfigurationError, match="four frequencies"):
        simulate_nucleotide_sequences(path, state_freqs=(0.5, 0.5))
    with pytest.raises(ConfigurationError, match="seed"):
        simulate_nucleotide_sequences(path, seed=True)
    with pytest.raises(ConfigurationError, match="seed"):
        simulate_nucleotide_sequences(path, seed=-1)
    with pytest.raises(ConfigurationError, match="length"):
        simulate_nucleotide_sequences(path, length=0)
    with pytest.raises(ConfigurationError, match="max_bases"):
        simulate_nucleotide_sequences(path, length=50, max_bases=100)
    with pytest.raises(ConfigurationError, match="branch length"):
        simulate_nucleotide_sequences(_write(tmp_path, "(a,b);\n", "nolen.nwk"))
    with pytest.raises(ConfigurationError, match="unique"):
        simulate_nucleotide_sequences(_write(tmp_path, "(a:1,a:1);\n", "dup.nwk"))
    with pytest.raises(ConfigurationError, match="at least two tips"):
        simulate_nucleotide_sequences(_write(tmp_path, "(a:1);\n", "one.nwk"))


def test_alignment_refuses_to_overwrite_output(tmp_path: Path) -> None:
    pytest.importorskip("pyvolve")
    path = _write(tmp_path, TREE_4)
    existing = tmp_path / "existing.fasta"
    existing.write_text(">old\nACGT\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        simulate_nucleotide_alignment(path, existing, length=10, seed=1)
    assert existing.read_text(encoding="utf-8") == ">old\nACGT\n"


def test_pyvolve_missing_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, TREE_4)
    monkeypatch.setitem(sys.modules, "pyvolve", None)
    with pytest.raises(BackendUnavailableError, match="bio-pyvolve"):
        simulate_nucleotide_sequences(path, length=10, seed=1)
