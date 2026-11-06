"""Tests for the P19 Biopython phylogeny, NEXUS and PAML functions (no network)."""

from __future__ import annotations

import io
import os
import random
import shutil
import warnings

import pytest

pytest.importorskip("Bio")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import biopython_p19_phylo_functions as m

GAPPED = {"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"}
ALN = {"a": "ACGTACGTAC", "b": "ACGTACGTTC", "c": "ACGAACGTTC", "d": "TCGAACGATC"}
NEWICK = "((A_x:1,A_y:2)0.9:1,(C:1,D:1)0.8:2,E:3);"


def _codon_alignment() -> dict[str, str]:
    rnd = random.Random(5)
    stops = {"TAA", "TAG", "TGA"}
    base = "".join(rnd.choice("ACGT") for _ in range(600))
    base = "".join(
        base[i : i + 3] if base[i : i + 3] not in stops else "GCT" for i in range(0, 600, 3)
    )

    def mutate(text: str, count: int) -> str:
        letters = list(text)
        for _ in range(count):
            letters[rnd.randrange(len(letters))] = rnd.choice("ACGT")
        mutated = "".join(letters)
        return "".join(
            mutated[i : i + 3] if mutated[i : i + 3] not in stops else "GCT"
            for i in range(0, len(mutated), 3)
        )

    return {"a": base, "b": mutate(base, 25), "c": mutate(base, 45), "d": mutate(base, 70)}


def _nexus(alignment: dict[str, str]):
    from Bio.Nexus import Nexus

    names = list(alignment)
    text = (
        "#NEXUS\nbegin data;\n"
        f"dimensions ntax={len(names)} nchar={len(alignment[names[0]])};\n"
        "format datatype=dna missing=? gap=-;\nmatrix\n"
        + "\n".join(f"{n} {alignment[n]}" for n in names)
        + "\n;\nend;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Nexus.Nexus(text)


def test_terminal_gaps_to_missing_matches_biopython() -> None:
    result = m.nexus_terminal_gaps_to_missing(GAPPED, verbose=False)
    reference = _nexus(GAPPED)
    reference.terminal_gap_to_missing()
    assert result["alignment"] == {k: str(v) for k, v in reference.matrix.items()}
    assert result["alignment"]["A"] == "??ACGT--A"
    assert (
        m.nexus_terminal_gaps_to_missing({"x": "-----"}, verbose=False)["alignment"]["x"] == "?????"
    )
    with pytest.raises(ConfigurationError):
        m.nexus_terminal_gaps_to_missing({"A": "AC", "B": "ACG"}, verbose=False)


def test_combine_matrices_fills_missing_taxa() -> None:
    result = m.nexus_combine_matrices(
        {"g1": GAPPED, "g2": {"A": "ACGTAA", "D": "TTTTGG"}}, verbose=False
    )
    assert result["taxa"] == ["A", "B", "C", "D"]
    assert result["character_count"] == 15
    assert result["alignment"]["D"] == "?" * 9 + "TTTTGG"
    assert result["alignment"]["B"].endswith("??????")
    assert result["charsets"]["g2"] == list(range(9, 15))
    with pytest.raises(ConfigurationError):
        m.nexus_combine_matrices({}, verbose=False)


def test_write_partition_files(tmp_path) -> None:
    out = tmp_path / "parts"
    result = m.nexus_write_partition_files(
        GAPPED, {"p1": [0, 1, 2], "p2": [3, 4, 5, 6, 7, 8]}, out, exclude=[1], verbose=False
    )
    assert set(result["files"]) == {"p1", "p2"}
    with open(result["files"]["p1"]) as handle:
        text = handle.read()
    assert "nchar=2" in text and "begin data;" in text
    assert result["partition_sizes"] == {"p1": 2, "p2": 6}
    with pytest.raises(ConfigurationError):
        m.nexus_write_partition_files(GAPPED, {"p1": [0]}, out, verbose=False)  # directory exists
    with pytest.raises(ConfigurationError):
        m.nexus_write_partition_files(GAPPED, {"p1": [99]}, tmp_path / "x", verbose=False)


def test_weighted_stepmatrix_matches_biopython() -> None:
    result = m.nexus_weighted_stepmatrix(GAPPED, name="steps", verbose=False)
    assert result["stepmatrix"] == _nexus(GAPPED).weighted_stepmatrix(name="steps")
    assert result["symbols"] == ["A", "C", "G", "T", "-"]
    assert result["weights"]["A"]["C"] == pytest.approx(1.39, abs=0.01)
    assert result["weights"]["A"]["G"] is None


def test_bootstrap_alignment_is_seeded_and_restores_state() -> None:
    random.seed(99)
    expected_next = random.random()
    random.seed(99)
    first = m.nexus_bootstrap_alignment(GAPPED, seed=3, verbose=False)["alignment"]
    assert random.random() == expected_next  # global state untouched
    assert first == m.nexus_bootstrap_alignment(GAPPED, seed=3, verbose=False)["alignment"]
    reference = _nexus(GAPPED)
    random.seed(3)
    assert first == {k: str(v) for k, v in reference.bootstrap().items()}
    assert {len(v) for v in first.values()} == {9}


def test_distance_matrix_identity_and_blastn() -> None:
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    msa = MultipleSeqAlignment([SeqRecord(Seq(s), id=n) for n, s in ALN.items()])
    for model in ("identity", "blastn"):
        result = m.biopython_alignment_distance_matrix(ALN, model=model, verbose=False)
        reference = DistanceCalculator(model).get_distance(msa)
        assert result["ids"] == list(ALN)
        assert result["matrix"][3][1] == pytest.approx(reference["d", "b"])
        assert result["matrix"][1][3] == result["matrix"][3][1]
    assert m.biopython_alignment_distance_matrix(ALN, verbose=False)["matrix"][0][
        1
    ] == pytest.approx(0.1)
    with pytest.raises(ConfigurationError):
        m.biopython_alignment_distance_matrix(ALN, model="blosum62", verbose=False)


def test_parsimony_tree_search() -> None:
    result = m.parsimony_tree_search(ALN, verbose=False)
    assert result["parsimony_score"] == 4.0
    for name in ALN:
        assert name in result["newick"]
    seeded = m.parsimony_tree_search(ALN, starting_tree="((a:1,b:1):1,(c:1,d:1):1);", verbose=False)
    assert seeded["starting_tree"] == "supplied"
    with pytest.raises(ConfigurationError):
        m.parsimony_tree_search(ALN, starting_tree="((a:1,b:1):1,(c:1,z:1):1);", verbose=False)
    with pytest.raises(ConfigurationError):
        m.parsimony_tree_search(
            ALN, starting_tree="((a,b),(c,d));", verbose=False
        )  # no branch lengths


def test_strict_consensus_matches_biopython_and_ignores_tip_order() -> None:
    from Bio import Phylo
    from Bio.Phylo import Consensus

    trees = ["(((A,B),C),(D,E));", "(((A,B),C),(D,E));", "((A,B),(C,(D,E)));"]
    result = m.strict_consensus_tree(trees, verbose=False)
    reference = Consensus.strict_consensus([Phylo.read(io.StringIO(t), "newick") for t in trees])
    ref_sets = {
        frozenset(t.name for t in c.get_terminals()) for c in reference.find_clades(terminal=False)
    }
    mine = Phylo.read(io.StringIO(result["newick"]), "newick")
    assert {
        frozenset(t.name for t in c.get_terminals()) for c in mine.find_clades(terminal=False)
    } == ref_sets
    # different tip order: (C,D) vs (C,E) must not agree
    clash = m.strict_consensus_tree(["((A,B),(C,D),E);", "((A,B),(C,E),D);"], verbose=False)
    assert clash["shared_clade_count"] == 1
    with pytest.raises(ConfigurationError):
        m.strict_consensus_tree(["((A,B),C);", "((A,B),D);"], verbose=False)


def test_collapse_genera() -> None:
    result = m.collapse_tree_genera(NEWICK, verbose=False)
    assert "A <collapsed>" in result["newick"]
    assert result["terminal_count_before"] == 5 and result["terminal_count_after"] == 4
    plain = m.collapse_tree_genera("((A_x,B_y),(A_z,B_w));", verbose=False)
    assert plain["terminal_count_after"] == 4


def test_split_clade_and_errors() -> None:
    result = m.split_tree_clade(NEWICK, "C", count=3, branch_length=0.5, verbose=False)
    assert result["new_clades"] == ["C0", "C1", "C2"]
    assert "C2:0.5" in result["newick"]
    with pytest.raises(ConfigurationError):
        m.split_tree_clade(NEWICK, "nope", verbose=False)


def test_random_tree_is_seeded() -> None:
    taxa = list("ABCDEFG")
    first = m.random_tree(taxa, seed=3, verbose=False)["newick"]
    assert first == m.random_tree(taxa, seed=3, verbose=False)["newick"]
    assert first != m.random_tree(taxa, seed=4, verbose=False)["newick"]
    for name in taxa:
        assert name in first
    assert m.tree_is_bifurcating(first, verbose=False)["is_bifurcating"]
    with pytest.raises(ConfigurationError):
        m.random_tree(["A", "A"], verbose=False)


def test_bifurcation_descendant_and_find() -> None:
    assert m.tree_is_bifurcating("((A,B),(C,D));", verbose=False)["is_bifurcating"]
    star = m.tree_is_bifurcating("(A,B,C,D);", verbose=False)
    assert star["is_bifurcating"] is False and star["internal_child_counts"] == [4]
    tree = "((A,B)AB,(C,D)CD)root;"
    assert m.tree_is_parent_of(tree, "AB", "A", verbose=False)["is_parent_of"] is True
    assert m.tree_is_parent_of(tree, "AB", "C", verbose=False)["is_parent_of"] is False
    found = m.find_tree_clades(NEWICK, name="A.*", verbose=False)
    assert [c["name"] for c in found["clades"]] == ["A_x", "A_y"]
    assert m.find_tree_clades(NEWICK, terminal=False, verbose=False)["match_count"] == 3
    assert m.find_tree_clades(NEWICK, name="C", first_only=True, verbose=False)["match_count"] == 1
    assert m.find_tree_clades(NEWICK, branch_length=3, verbose=False)["clades"][0]["name"] == "E"
    with pytest.raises(ConfigurationError):
        m.find_tree_clades(NEWICK, name="(", verbose=False)


def test_identical_and_compatible() -> None:
    assert (
        m.trees_are_identical("((A,B),(C,D));", "((D,C),(B,A));", verbose=False)["identical"]
        is True
    )
    assert (
        m.trees_are_identical("((A,B),(C,D));", "((A,C),(B,D));", verbose=False)["identical"]
        is False
    )
    ok = m.trees_are_compatible(
        "((A,B)90:1,(C,D)90:1,E);", "((A,B)95:1,(C,D)80:1,E);", threshold=50, verbose=False
    )
    assert ok["compatible"] is True
    clash = m.trees_are_compatible(
        "((A,B)90:1,(C,D)90:1,E);", "((A,C)90:1,(B,D)90:1,E);", threshold=50, verbose=False
    )
    assert clash["compatible"] is False and clash["conflict_count"] >= 1
    with pytest.raises(ConfigurationError):
        m.trees_are_compatible("((A,B)90:1,C);", "((A,B)90:1,D);", threshold=50, verbose=False)


def test_convert_support() -> None:
    result = m.convert_tree_support(
        "((A:1,B:2)90:1,(C:1,D:1)80:2,E:3);", replicates=100, verbose=False
    )
    supports = [n["support"] for n in result["nodes"] if n["support"] is not None]
    assert supports == [0.9, 0.8]
    moved = m.convert_tree_support(
        "((A:0.9,B:2)0.9:1,(C:1,D:1)0.8:2,E:3);", mode="branchlength_to_support", verbose=False
    )
    assert all(n["branch_length"] == 0.0 for n in moved["nodes"])
    with pytest.raises(ConfigurationError):
        m.convert_tree_support(NEWICK, verbose=False)  # replicates missing


@pytest.mark.parametrize("method", ["NG86", "LWL85", "YN00", "ML"])
def test_codon_dn_ds_matches_biopython(method) -> None:
    import numpy as np
    from Bio.Align import Alignment, analysis
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aln = _codon_alignment()
    if method == "ML":  # maximum likelihood is slow: three taxa, 60 bases
        aln = {name: text[:60] for name, text in list(aln.items())[:3]}
    length = len(next(iter(aln.values())))
    result = m.codon_dn_ds_estimates(aln, method=method, verbose=False)
    records = [SeqRecord(Seq(s), id=n) for n, s in aln.items()]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dn, ds = analysis.calculate_dn_ds_matrix(
            Alignment(records, np.array([[0, length]] * len(aln))),
            method=method,
            codon_table=CodonTable.generic_by_id[1],
        )
    assert result["dn"][2][1] == pytest.approx(float(dn.matrix[2][1]))
    assert result["ds"][2][0] == pytest.approx(float(ds.matrix[2][0]))
    assert result["dn"][0][2] == result["dn"][2][0]


def test_codon_dn_ds_edge_cases() -> None:
    with pytest.raises(ConfigurationError):
        m.codon_dn_ds_estimates({"a": "ATGA", "b": "ATGC"}, verbose=False)  # not codons
    with pytest.raises(ConfigurationError):
        m.codon_dn_ds_estimates({"a": "ATG---", "b": "ATGAAA"}, verbose=False)  # gaps refused
    with pytest.raises(ConfigurationError):
        m.codon_dn_ds_estimates(_codon_alignment(), method="XX", verbose=False)


def test_codon_dn_ds_trees() -> None:
    result = m.codon_dn_ds_trees(
        _codon_alignment(), dn_ds_method="NG86", tree_method="NJ", verbose=False
    )
    for key in ("dn_tree_newick", "ds_tree_newick"):
        for name in "abcd":
            assert name in result[key]
    upgma = m.codon_dn_ds_trees(_codon_alignment(), verbose=False)
    assert upgma["tree_method"] == "UPGMA"


def test_read_paml_results_errors(tmp_path) -> None:
    empty = tmp_path / "empty.out"
    empty.write_text("")
    with pytest.raises(ConfigurationError):
        m.read_paml_results("codeml", empty, verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_paml_results("codeml", tmp_path / "missing.out", verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_paml_results("paup", empty, verbose=False)


def test_run_paml_validation_and_missing_program(tmp_path) -> None:
    phy = tmp_path / "a.phy"
    phy.write_text(" 2 6\nA   ATGAAA\nB   ATGAAG\n")
    with pytest.raises(ConfigurationError):
        m.run_paml("codeml", phy, tmp_path / "w1", options={"not_an_option": 1}, verbose=False)
    with pytest.raises(ConfigurationError):
        m.run_paml("yn00", phy, tmp_path / "w2", tree_path=phy, verbose=False)
    if shutil.which("yn00") is None:
        with pytest.raises(BackendUnavailableError):
            m.run_paml("yn00", phy, tmp_path / "w3", verbose=False)


def _paml_inputs(tmp_path):
    aln = _codon_alignment()
    phy = tmp_path / "a.phy"
    phy.write_text(" 4 600\n" + "\n".join(f"{n}  {s}" for n, s in aln.items()) + "\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("((a:0.1,b:0.1):0.05,(c:0.1,d:0.1):0.05);\n")
    return phy, tree


@pytest.mark.skipif(shutil.which("yn00") is None, reason="PAML is not installed")
def test_run_paml_yn00_end_to_end(tmp_path) -> None:
    phy, _ = _paml_inputs(tmp_path)
    result = m.run_paml(
        "yn00", phy, tmp_path / "w", options={"icode": 0, "verbose": 0}, verbose=False
    )
    assert result["status"] == "completed"
    assert set(result["results"]) == {"a", "b", "c", "d"}
    assert 0 < result["results"]["a"]["b"]["YN00"]["dN"] < 1
    parsed = m.read_paml_results("yn00", result["results_file"], verbose=False)
    assert parsed["results"] == result["results"]
    assert os.path.exists(result["control_file"])


@pytest.mark.skipif(shutil.which("codeml") is None, reason="PAML is not installed")
def test_run_paml_codeml_end_to_end(tmp_path) -> None:
    phy, tree = _paml_inputs(tmp_path)
    options = {"runmode": 0, "seqtype": 1, "CodonFreq": 2, "model": 0, "NSsites": [0], "verbose": 0}
    result = m.run_paml(
        "codeml", phy, tmp_path / "w", tree_path=tree, options=options, verbose=False
    )
    assert result["status"] == "completed"
    site = result["results"]["NSsites"]["0"]
    assert site["lnL"] < 0 and site["parameters"]["omega"] > 0
    assert (
        m.read_paml_results("codeml", result["results_file"], verbose=False)["results"]
        == result["results"]
    )


@pytest.mark.skipif(shutil.which("baseml") is None, reason="PAML is not installed")
def test_run_paml_baseml_end_to_end(tmp_path) -> None:
    phy, tree = _paml_inputs(tmp_path)
    result = m.run_paml(
        "baseml", phy, tmp_path / "w", tree_path=tree, options={"model": 4}, verbose=False
    )
    assert result["status"] == "completed"
    assert result["results"]["parameters"]["kappa"] > 0
    assert (
        m.read_paml_results("baseml", result["results_file"], verbose=False)["results"]
        == result["results"]
    )
