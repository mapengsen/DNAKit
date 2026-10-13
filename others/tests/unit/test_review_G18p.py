"""Tests for the G18p review wrappers: toytree, dendropy and pyranges1 functions.

Expectations are hand calculations on tiny trees and intervals, or direct calls to the upstream
package. Tests skip when an upstream package is missing and never use the network.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from dnakit.annotation import pyranges1_review_functions as pr_review
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import dendropy_review_functions as dp_review
from dnakit.phylogeny import toytree_review_functions as tt_review

TREE = "((A:1,B:2):1,(C:3,D:1):2);"
NETWORK = "((A:1,(B:1,#H1:0::0.3):1):1,(#H1:0::0.7,C:2):1);"


def _matrix_close(actual: list, expected: list) -> None:
    """Compare two matrices row by row with approximate floats."""

    assert len(actual) == len(expected)
    for row, wanted in zip(actual, expected, strict=True):
        assert list(row) == pytest.approx(list(wanted))


# ---------------------------------------------------------------- toytree: node queries


def test_toytree_node_path_and_distances_match_hand_values():
    pytest.importorskip("toytree")
    path = tt_review.list_tree_node_path(TREE, "A", "C")["result"]
    assert [node["idx"] for node in path["nodes"]] == [0, 4, 6, 5, 2]
    assert path["edge_count"] == 4
    patristic = tt_review.calculate_tree_node_distance(TREE, "A", "C")["result"]
    assert patristic["distance"] == pytest.approx(7.0)  # 1 + 1 + 2 + 3
    topo = tt_review.calculate_tree_node_distance(TREE, "A", "C", topology_only=True)
    assert topo["result"]["distance"] == pytest.approx(4.0)


def test_toytree_all_node_matrix_is_symmetric_and_matches_upstream():
    toytree = pytest.importorskip("toytree")
    from toytree import distance

    result = tt_review.calculate_all_node_distance_matrix(TREE)["result"]
    assert result["labels"] == ["A", "B", "C", "D", "node_4", "node_5", "node_6"]
    assert result["matrix"][0] == pytest.approx([0, 3, 7, 5, 1, 4, 2])  # A to every node
    upstream = distance.get_node_distance_matrix(toytree.tree(TREE))
    _matrix_close(result["matrix"], upstream.tolist())
    for i, row in enumerate(result["matrix"]):
        for j, value in enumerate(row):
            assert value == pytest.approx(result["matrix"][j][i])


def test_toytree_internal_matrix_and_descendants_hand_values():
    pytest.importorskip("toytree")
    internal = tt_review.calculate_internal_node_distance_matrix(TREE)["result"]
    assert internal["labels"] == ["node_4", "node_5", "node_6"]
    _matrix_close(internal["matrix"], [[0, 3, 1], [3, 0, 2], [1, 2, 0]])
    desc = tt_review.calculate_descendant_distances(TREE)["result"]["descendants"]
    by_name = {row["name"]: row["distance"] for row in desc if row["name"]}
    assert by_name == pytest.approx({"A": 2.0, "B": 3.0, "C": 5.0, "D": 3.0})
    under_a = tt_review.calculate_descendant_distances(TREE, "A")["result"]["descendants"]
    assert [row["distance"] for row in under_a] == [0.0]


def test_toytree_farthest_node_hand_value():
    pytest.importorskip("toytree")
    far = tt_review.find_tree_farthest_node(TREE, "A")["result"]
    assert far["farthest"]["name"] == "C"
    assert far["distance"] == pytest.approx(7.0)  # A to C: 1 + 1 + 2 + 3
    value = tt_review.calculate_farthest_node_distance(TREE, "A")["result"]["distance"]
    assert value == pytest.approx(7.0)


# ---------------------------------------------------------------- toytree: enumeration


def test_toytree_quadripartitions_and_quartets_hand_values():
    pytest.importorskip("toytree")
    quad4 = tt_review.list_tree_quadripartitions(TREE)["result"]
    assert quad4["quadripartitions"] == [[["A"], ["B"], ["C"], ["D"]]]
    quad5 = tt_review.list_tree_quadripartitions("((A,B),(C,D),E);")["result"]
    assert quad5["count"] == 2
    assert [["A"], ["B"], ["C", "D"], ["E"]] in quad5["quadripartitions"]
    quartets = tt_review.list_tree_quartets("((A,B),(C,D),E);")["result"]
    pairs = {(tuple(row["left"]), tuple(row["right"])) for row in quartets["quartets"]}
    expected = {
        (("A", "B"), ("C", "D")),
        (("A", "B"), ("C", "E")),
        (("A", "B"), ("D", "E")),
        (("C", "D"), ("A", "E")),
        (("C", "D"), ("B", "E")),
    }
    assert pairs == expected
    assert quartets["count"] == 5


def test_toytree_consensus_support_is_fraction_of_trees():
    pytest.importorskip("toytree")
    target = "((A,B),(C,D));"
    trees = ["((A,B),(C,D));", "((A,B),(C,D));", "((A,C),(B,D));"]
    result = tt_review.map_consensus_tree_features(target, trees)["result"]
    node = next(row for row in result["nodes"] if row["tips"] == ["C", "D"])
    assert node["values"]["support"] == pytest.approx(2 / 3)
    assert result["tree_count"] == 3
    with pytest.raises(ConfigurationError, match="same tips"):
        tt_review.map_consensus_tree_features(target, ["((A,B),(C,E));"])


# ---------------------------------------------------------------- toytree: editing


def test_toytree_alignment_and_height_edits_hand_values():
    pytest.importorskip("toytree")
    aligned = tt_review.extend_tree_tips_to_align(TREE)["result"]
    assert aligned["newick"] == "((A:4,B:4):1,(C:3,D:3):2);"  # all tips end at depth 5
    import toytree
    from toytree import mod

    heights = tt_review.set_tree_node_heights(TREE, {"A": 2.0})["result"]
    assert heights["newick"] == mod.edges_set_node_heights(toytree.tree(TREE), {"A": 2.0}).write()


def test_toytree_rotate_and_drop_hand_values():
    pytest.importorskip("toytree")
    rotated = tt_review.rotate_tree_node(TREE, 4)["result"]
    assert rotated["newick"] == "((B:2,A:1):1,(C:3,D:1):2);"
    dropped = tt_review.drop_tree_tips(TREE, ["A"])["result"]
    assert dropped["newick"] == "((C:3,D:1):2,B:3);"  # B's edge becomes 2 + 1
    assert dropped["tips"] == ["B", "C", "D"]
    with pytest.raises(ConfigurationError, match="Unknown tip"):
        tt_review.drop_tree_tips(TREE, ["Z"])
    with pytest.raises(ConfigurationError, match="at least two tips"):
        tt_review.drop_tree_tips(TREE, ["A", "B", "C"])


def test_toytree_bisect_splits_clade_from_remainder():
    pytest.importorskip("toytree")
    split = tt_review.bisect_tree(TREE, 4)["result"]
    assert split["subtree_newick"] == "(A:1,B:2):1;"
    assert split["subtree_tips"] == ["A", "B"]
    assert split["remainder_newick"] == "((C:3,D:1):2);"
    assert split["remainder_tips"] == ["C", "D"]


def test_toytree_resolve_polytomies_is_seeded_and_binary():
    toytree = pytest.importorskip("toytree")
    from toytree import mod

    first = tt_review.resolve_tree_polytomies("(A,B,C,D);", seed=3)["result"]
    second = tt_review.resolve_tree_polytomies("(A,B,C,D);", seed=3)["result"]
    assert first["newick"] == second["newick"]
    assert first["newick"] == mod.resolve_polytomies(toytree.tree("(A,B,C,D);"), seed=3).write()
    assert toytree.tree(first["newick"]).nnodes == 7  # bifurcating 4-tip tree
    with pytest.raises(ConfigurationError, match="seed"):
        tt_review.resolve_tree_polytomies("(A,B,C,D);", seed=True)


def test_toytree_add_node_functions_match_upstream_and_hand_counts():
    toytree = pytest.importorskip("toytree")
    from toytree import mod

    source = toytree.tree(TREE)
    internal = tt_review.add_tree_internal_node(TREE, "A", name="X", dist=0.5)["result"]
    assert internal["newick"] == mod.add_internal_node(source, "A", dist=0.5, name="X").write()
    assert internal["name"] == "X"
    child = tt_review.add_tree_child_node(TREE, 4, name="Z", dist=0.5)["result"]
    assert child["newick"] == "((A:1,B:2,Z:0.5):1,(C:3,D:1):2);"
    sister = tt_review.add_tree_sister_node(TREE, "A", name="Z", dist=0.5)["result"]
    assert sister["tips"] == ["A", "B", "C", "D", "Z"]
    both = tt_review.add_tree_internal_node_and_child(
        TREE, "A", name="Z", dist=0.5, parent_name="P", parent_dist=0.5
    )["result"]
    assert both["tips"] == ["A", "B", "C", "D", "Z"]
    graft = tt_review.add_tree_internal_node_and_subtree(
        TREE, "A", "(X:1,Y:1);", subtree_stem_dist=0.5, parent_name="P", parent_dist=0.5
    )["result"]
    assert graft["tips"] == ["A", "B", "C", "D", "X", "Y"]
    assert graft["subtree_tips"] == ["X", "Y"]
    with pytest.raises(ConfigurationError, match="must not repeat"):
        tt_review.add_tree_internal_node_and_subtree(TREE, "A", "(A:1,Y:1);")


def test_toytree_remove_and_merge_nodes():
    pytest.importorskip("toytree")
    removed = tt_review.remove_tree_nodes(TREE, [4])["result"]
    assert removed["newick"] == "((C:3,D:1):2,A:2,B:3);"  # preserved distances
    with pytest.raises(ConfigurationError, match="at least two tips"):
        tt_review.remove_tree_nodes(TREE, [0, 1, 2])
    merged = tt_review.merge_tree_nodes("((A:1,B:1):1,(C:1,D:1):1);", "dist")["result"]
    assert merged["newick"] == "(A:2,C:2);"
    assert merged["removed_node_count"] == 4
    with pytest.raises(ConfigurationError, match="identifier"):
        tt_review.merge_tree_nodes(TREE, "bad feature")


def test_toytree_network_parsing_matches_upstream():
    toytree = pytest.importorskip("toytree")
    from toytree.network import parse_network_to_tree_and_admix

    result = tt_review.parse_tree_network(NETWORK)["result"]
    major, admix = parse_network_to_tree_and_admix(NETWORK)
    assert result["newick"] == major.write()
    assert result["tips"] == ["A", "B", "C"]
    assert result["admixture_count"] == len(admix) == 1
    edge = result["admixture_edges"][0]
    assert edge["source_tips"] == ["B"]
    assert edge["gamma"] == pytest.approx(0.3)
    assert toytree is not None


# ---------------------------------------------------------------- toytree: tables


def test_vcv_distance_matrix_is_v_ii_plus_v_jj_minus_twice_v_ij():
    pytest.importorskip("toytree")
    result = tt_review.calculate_tip_distance_matrix_from_vcv([[2, 1], [1, 2]], ["x", "y"])
    _matrix_close(result["result"]["matrix"], [[0.0, 2.0], [2.0, 0.0]])
    with pytest.raises(ConfigurationError, match="symmetric"):
        tt_review.calculate_tip_distance_matrix_from_vcv([[2, 1], [0, 2]], ["x", "y"])
    with pytest.raises(ConfigurationError, match=r"square|finite"):
        tt_review.calculate_tip_distance_matrix_from_vcv([[1, 0, 0], [0, 1, 0]], ["x", "y"])


def test_toytree_tip_diversification_matches_upstream_for_tips_only():
    toytree = pytest.importorskip("toytree")
    from toytree import pcm

    rates = tt_review.calculate_tip_diversification_rates(TREE)["result"]["rates"]
    upstream = pcm.get_tip_level_diversification(toytree.tree(TREE))
    expected = {str(label): float(upstream[label]) for label in ["A", "B", "C", "D"]}
    assert {row["tip"]: row["diversification_rate"] for row in rates} == pytest.approx(expected)


def test_toytree_module_reports_missing_package(monkeypatch):
    monkeypatch.setitem(sys.modules, "toytree", None)
    with pytest.raises(BackendUnavailableError) as info:
        tt_review.list_tree_node_path(TREE, "A", "C")
    assert "bio-toytree" in (info.value.hint or "")


# ---------------------------------------------------------------- dendropy


def test_missing_bipartitions_hand_values_for_five_tips():
    pytest.importorskip("dendropy")
    result = dp_review.list_missing_bipartitions("((A,B),(C,D),E);", "((A,C),(B,D),E);")["result"]
    assert result["bipartitions"] == [["A", "B"], ["C", "D"]]
    assert result["missing_count"] == 2
    with pytest.raises(ConfigurationError):
        dp_review.list_missing_bipartitions("((A,B),C);", "((A,B),D);")


def test_n_bar_hand_values():
    dendropy = pytest.importorskip("dendropy")
    from dendropy.calculate import treemeasure

    # Tips A and B have 2 nodes above them, C has 1: mean (2 + 2 + 1) / 3.
    assert dp_review.calculate_n_bar_index("((A,B),C);")["result"]["n_bar"] == pytest.approx(5 / 3)
    balanced = dp_review.calculate_n_bar_index(TREE)["result"]["n_bar"]
    upstream = treemeasure.N_bar(dendropy.Tree.get(data=TREE, schema="newick"))
    assert balanced == pytest.approx(2.0) == pytest.approx(upstream)


def test_time_vectors_hand_values_on_weighted_tree():
    pytest.importorskip("dendropy")
    # Forward depths: root 0, cherry AB 1, cherry CD 2, tips A 2, B 3, C 5, D 3.
    depths = dp_review.calculate_node_depths(TREE)["result"]["values"]
    assert depths == pytest.approx([0, 1, 2, 2, 3, 3, 5])
    # Backward ages are 5 minus the depths.
    ages = dp_review.calculate_node_ages(TREE)["result"]["values"]
    assert ages == pytest.approx([0, 2, 2, 3, 3, 4, 5])
    assert dp_review.calculate_coalescence_ages(TREE)["result"]["values"] == pytest.approx(
        [3, 4, 5]
    )
    assert dp_review.calculate_divergence_times(TREE)["result"]["values"] == pytest.approx(
        [0, 1, 2]
    )
    internal_depths = dp_review.calculate_node_depths(TREE, internal_only=True)["result"]
    assert internal_depths["values"] == pytest.approx([0, 1, 2])
    internal_ages = dp_review.calculate_node_ages(TREE, internal_only=True)["result"]
    assert internal_ages["values"] == pytest.approx([3, 4, 5])


def test_time_measures_need_branch_lengths():
    pytest.importorskip("dendropy")
    with pytest.raises(ConfigurationError, match="branch length"):
        dp_review.calculate_divergence_times("((A,B),(C,D));")


def test_dendropy_module_reports_missing_package(monkeypatch):
    monkeypatch.setitem(sys.modules, "dendropy", None)
    with pytest.raises(BackendUnavailableError) as info:
        dp_review.calculate_n_bar_index("((A,B),C);")
    assert "bio-dendropy" in (info.value.hint or "")


# ---------------------------------------------------------------- pyranges1 fixtures


GTF = (
    'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\tCDS\t100\t130\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\tCDS\t200\t260\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\tCDS\t400\t450\t.\t-\t0\tgene_id "g2"; transcript_id "t2";\n'
)
GFF3 = (
    "##gff-version 3\n"
    "chr1\tsrc\tgene\t101\t500\t.\t+\t.\tID=g1;Name=G1\n"
    "chr1\tsrc\tCDS\t101\t130\t.\t+\t0\tID=c1;Parent=g1\n"
)
BED_X = "chr1\t0\t100\tx1\t0\t+\nchr1\t200\t300\tx2\t0\t-\nchr2\t50\t80\tx3\t0\t+\n"
BED_Y = "chr1\t50\t150\ty1\t0\t+\nchr2\t60\t70\ty2\t0\t+\n"
FASTA = ">chr1\nGGGATGGCTGCAAAGTGATTTCCCATGCCCGCGTAAGGGCCC\n"


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def files(tmp_path: Path) -> dict[str, Path]:
    return {
        "gtf": _write(tmp_path, "a.gtf", GTF),
        "gff3": _write(tmp_path, "a.gff3", GFF3),
        "x": _write(tmp_path, "x.bed", BED_X),
        "y": _write(tmp_path, "y.bed", BED_Y),
        "cds": _write(tmp_path, "cds.bed", "chr1\t9\t21\tg\t0\t+\n"),
        "fasta": _write(tmp_path, "g.fa", FASTA),
    }


# ---------------------------------------------------------------- pyranges1 readers


def test_gff_readers_keep_zero_based_starts_and_attributes(files):
    pytest.importorskip("pyranges1")
    gtf = pr_review.read_pyranges1_gff(files["gtf"])
    assert gtf["record_count"] == 4
    gene = next(r for r in gtf["records"] if r["Feature"] == "gene")
    assert (gene["chrom"], gene["start"], gene["end"]) == ("chr1", 99, 500)
    assert gene["transcript_id"] == "t1"
    gff3 = pr_review.read_pyranges1_gff3(files["gff3"])
    gene3 = next(r for r in gff3["records"] if r["Feature"] == "gene")
    assert gene3["start"] == 100 and gene3["ID"] == "g1"
    full = pr_review.read_pyranges1_gtf_full(files["gtf"])
    assert full["record_count"] == 4


def test_reader_row_limit_is_enforced(files):
    pytest.importorskip("pyranges1")
    with pytest.raises(ConfigurationError, match="max_rows"):
        pr_review.read_pyranges1_gff(files["gtf"], max_rows=2)


def test_bigwig_reader_returns_track_intervals(tmp_path):
    pyBigWig = pytest.importorskip("pyBigWig")
    pytest.importorskip("pyranges1")
    path = tmp_path / "s.bw"
    handle = pyBigWig.open(str(path), "w")
    handle.addHeader([("chr1", 1000), ("chr2", 500)])
    handle.addEntries(
        ["chr1", "chr1", "chr2"], [0, 200, 50], ends=[100, 260, 80], values=[1.5, 2.0, 0.25]
    )
    handle.close()
    result = pr_review.read_pyranges1_bigwig(path)
    rows = [(r["chrom"], r["start"], r["end"], r["value"]) for r in result["records"]]
    assert rows == [("chr1", 0, 100, 1.5), ("chr1", 200, 260, 2.0), ("chr2", 50, 80, 0.25)]
    with pytest.raises(ConfigurationError, match="max_file_bytes"):
        pr_review.read_pyranges1_bigwig(path, max_file_bytes=10)


def test_concatenate_sorts_intervals_from_both_files(files):
    pytest.importorskip("pyranges1")
    result = pr_review.concatenate_pyranges1_intervals([files["x"], files["y"]])
    keys = [(r["chrom"], r["start"], r["end"]) for r in result["records"]]
    assert keys == [
        ("chr1", 0, 100),
        ("chr1", 50, 150),
        ("chr1", 200, 300),
        ("chr2", 50, 80),
        ("chr2", 60, 70),
    ]
    assert result["source_count"] == 2


# ---------------------------------------------------------------- pyranges1 ORFs


def test_cds_frame_uses_cumulative_cds_length_modulo_three(tmp_path):
    pytest.importorskip("pyranges1")
    text = (
        'chr1\ts\tCDS\t1\t10\t.\t+\t0\ttranscript_id "t1";\n'
        'chr1\ts\tCDS\t101\t110\t.\t+\t0\ttranscript_id "t1";\n'
        'chr1\ts\tCDS\t201\t210\t.\t+\t0\ttranscript_id "t1";\n'
    )
    path = _write(tmp_path, "three.gtf", text)
    frames = pr_review.calculate_pyranges1_cds_frame(path)["records"]
    assert [row["frame"] for row in frames] == [0, 1, 2]  # 0, 10 mod 3, 20 mod 3
    with pytest.raises(ConfigurationError, match="No feature"):
        pr_review.calculate_pyranges1_cds_frame(path, feature="exon")


def test_extend_orfs_hand_coordinates(files):
    pytest.importorskip("pyranges1")
    pytest.importorskip("pyfaidx")
    # The 42-base sequence has ATG at 3 and 25, and TAA at 33; the CDS is [9, 21).
    up = pr_review.extend_pyranges1_orfs(files["cds"], files["fasta"], direction="up")
    assert (up["records"][0]["start"], up["records"][0]["end"]) == (3, 21)
    down = pr_review.extend_pyranges1_orfs(files["cds"], files["fasta"], direction="down")
    assert (down["records"][0]["start"], down["records"][0]["end"]) == (9, 36)
    both = pr_review.extend_pyranges1_orfs(files["cds"], files["fasta"])
    assert (both["records"][0]["start"], both["records"][0]["end"]) == (3, 36)
    bed3 = _write(files["cds"].parent, "cds3.bed", "chr1\t9\t21\n")
    with pytest.raises(ConfigurationError, match="BED6"):
        pr_review.extend_pyranges1_orfs(bed3, files["fasta"])


# ---------------------------------------------------------------- pyranges1 statistics


def test_forbes_matches_upstream(files):
    pytest.importorskip("pyranges1")
    from pyranges1 import read_bed
    from pyranges1.stats import forbes

    sizes = {"chr1": 1000, "chr2": 500}
    value = pr_review.calculate_pyranges1_forbes(files["x"], files["y"], chromsizes=sizes)
    upstream = forbes(read_bed(str(files["x"])), read_bed(str(files["y"])), chromsizes=sizes)
    assert value["forbes"] == pytest.approx(float(upstream))


def test_mcc_counts_match_hand_base_pair_overlaps(files):
    pytest.importorskip("pyranges1")
    result = pr_review.calculate_pyranges1_mcc([files["x"], files["y"]], ["x", "y"])
    rows = {(r["truth_label"], r["prediction_label"]): r for r in result["records"]}
    # x covers 230 bp, y covers 110 bp, they share 60 bp; genome extent is 380 bp.
    xy = rows[("x", "y")]
    assert (xy["tp"], xy["fp"], xy["fn"], xy["tn"]) == (60, 50, 170, 100)
    expected = (60 * 100 - 50 * 170) / math.sqrt(110 * 230 * 150 * 270)
    assert xy["mcc"] == pytest.approx(expected)
    assert rows[("x", "x")]["mcc"] == pytest.approx(1.0)
    with pytest.raises(ConfigurationError, match="unique"):
        pr_review.calculate_pyranges1_mcc([files["x"], files["y"]], ["a", "a"])


def test_relative_distance_matches_upstream(files):
    pytest.importorskip("pyranges1")
    from pyranges1 import read_bed
    from pyranges1.stats import relative_distance

    result = pr_review.calculate_pyranges1_relative_distance(files["x"], files["y"])
    upstream = relative_distance(read_bed(str(files["x"])), read_bed(str(files["y"])))
    assert result["record_count"] == len(upstream)
    assert result["records"][0]["fraction"] == pytest.approx(float(upstream["fraction"].iloc[0]))


def test_pyranges1_module_reports_missing_package(monkeypatch, files):
    monkeypatch.setitem(sys.modules, "pyranges1", None)
    with pytest.raises(BackendUnavailableError) as info:
        pr_review.read_pyranges1_gff(files["gtf"])
    assert "bio-pyranges1" in (info.value.hint or "")
