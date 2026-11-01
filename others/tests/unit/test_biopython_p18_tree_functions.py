"""Tests for dnakit.phylogeny.biopython_p18_tree_functions (tree files, PHYLIP distances)."""

from __future__ import annotations

import io

import pytest

pytest.importorskip("Bio")

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix

from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny import biopython_p18_tree_functions as m

QUIET = {"verbose": False}
TREES = ["((A:0.1,B:0.2)90:0.3,(C:0.4,D:0.5):0.1);", "(A:1,(B:1,C:1):1);"]


@pytest.mark.parametrize("fmt", ["newick", "nexus", "phyloxml", "nexml", "cdao"])
def test_write_tree_file_round_trip(tmp_path, fmt):
    out = tmp_path / f"trees.{fmt}"
    result = m.write_tree_file(TREES, out, file_format=fmt, **QUIET)
    assert result["tree_count"] == 2 and result["bytes_written"] == out.stat().st_size
    back = list(Phylo.parse(str(out), fmt))
    tips = [sorted(tip.name for tip in tree.get_terminals()) for tree in back]
    assert sorted(tips, key=len) == [["A", "B", "C"], ["A", "B", "C", "D"]]  # CDAO may reorder
    four = next(tree for tree in back if len(tree.get_terminals()) == 4)
    lengths = sorted(round(tip.branch_length, 6) for tip in four.get_terminals())
    assert lengths == [0.1, 0.2, 0.4, 0.5]


def test_write_tree_file_matches_biopython_text(tmp_path):
    out = tmp_path / "t.nwk"
    m.write_tree_file(TREES[:1], out, file_format="newick", **QUIET)
    tree = Phylo.read(io.StringIO(TREES[0]), "newick")
    buffer = io.StringIO()
    Phylo.write([tree], buffer, "newick")
    assert out.read_text(encoding="utf-8") == buffer.getvalue()


def test_write_tree_file_errors(tmp_path):
    out = tmp_path / "t.nwk"
    out.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        m.write_tree_file(TREES, out, **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_tree_file([], tmp_path / "a.nwk", **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_tree_file("(A,B);", tmp_path / "b.nwk", **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_tree_file(["((A,B"], tmp_path / "c.nwk", **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_tree_file(TREES, tmp_path / "d.nwk", file_format="svg", **QUIET)
    assert not (tmp_path / "c.nwk").exists()


def test_format_phylip_distance_matrix_matches_biopython(tmp_path):
    names = ["A", "B", "C"]
    full = [[0, 0.2, 0.3], [0.2, 0, 0.4], [0.3, 0.4, 0]]
    out = tmp_path / "d.phy"
    result = m.format_phylip_distance_matrix(names, full, output_path=out, **QUIET)
    buffer = io.StringIO()
    DistanceMatrix(names, [[0], [0.2, 0], [0.3, 0.4, 0]]).format_phylip(buffer)
    assert result["text"] == buffer.getvalue()
    assert out.read_text(encoding="utf-8") == result["text"]
    assert result["text"].splitlines()[0].strip() == "3"


def test_format_phylip_distance_matrix_errors():
    with pytest.raises(ConfigurationError, match="symmetric"):
        m.format_phylip_distance_matrix(["A", "B"], [[0, 1], [2, 0]], **QUIET)
    with pytest.raises(ConfigurationError, match="diagonal"):
        m.format_phylip_distance_matrix(["A", "B"], [[1, 1], [1, 0]], **QUIET)
    with pytest.raises(ConfigurationError, match="matrix"):
        m.format_phylip_distance_matrix(["A", "B"], [[0, 1]], **QUIET)
    with pytest.raises(ConfigurationError, match="unique"):
        m.format_phylip_distance_matrix(["A", "A"], [[0, 1], [1, 0]], **QUIET)


def test_public_functions_have_verbose_keyword_and_docstrings():
    import inspect

    for name in m.__all__:
        function = getattr(m, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
