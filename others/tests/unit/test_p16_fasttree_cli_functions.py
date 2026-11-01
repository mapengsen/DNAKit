"""FastTree facades (skip without FastTree)."""

from __future__ import annotations

from pathlib import Path

import pytest
from _p16_data import alignment_fasta, need

import dnakit.phylogeny.p16_fasttree_cli_functions as ft
from dnakit.exceptions import ConfigurationError


def test_distance_matrix_and_tree_modes(tmp_path: Path) -> None:
    program = need("FastTree")
    alignment = alignment_fasta(tmp_path)
    matrix = ft.compute_fasttree_distance_matrix(alignment, str(tmp_path / "d.phy"), fasttree_path=program, verbose=False)
    assert matrix["labels"] == ["t1", "t2", "t3", "t4", "t5"]
    assert matrix["matrix"][0][0] == 0.0 and matrix["matrix"][0][1] == matrix["matrix"][1][0] > 0
    tree = ft.infer_fasttree_tree(alignment, str(tmp_path / "t.nwk"), fasttree_path=program, verbose=False)
    assert tree["tips"] == 5 and tree["version"]
    noml = ft.infer_fasttree_tree(alignment, str(tmp_path / "noml.nwk"), refine="noml", fasttree_path=program, verbose=False)
    assert noml["newick"] != tree["newick"]
    refined = ft.infer_fasttree_tree(
        alignment, str(tmp_path / "ref.nwk"), starting_tree=str(tmp_path / "t.nwk"), refine="mllen", bootstrap_replicates=10,
        fasttree_path=program, verbose=False,
    )
    assert refined["tips"] == 5


def test_option_conflicts(tmp_path: Path) -> None:
    alignment = alignment_fasta(tmp_path)
    with pytest.raises(ConfigurationError):
        ft.infer_fasttree_tree(alignment, str(tmp_path / "a"), refine="mllen_fixed_topology", verbose=False)
    with pytest.raises(ConfigurationError):
        ft.infer_fasttree_tree(alignment, str(tmp_path / "b"), bootstrap_replicates=5, no_support=True, verbose=False)
    with pytest.raises(ConfigurationError):
        ft.infer_fasttree_tree(alignment, str(tmp_path / "c"), constraint_weight=10.0, verbose=False)
