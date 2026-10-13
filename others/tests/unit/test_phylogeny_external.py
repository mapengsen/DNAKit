"""Contract tests for the formal IQ-TREE 2 phylogeny facade."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny import infer_fast_tree, infer_ml_tree
from dnakit.tools import default_tool_registry


def _fake_iqtree(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "prefix = Path(sys.argv[sys.argv.index('-pre') + 1])\n"
        "prefix.with_suffix('.treefile').write_text('(a:0.1,b:0.2);', encoding='ascii')\n"
        "prefix.with_suffix('.iqtree').write_text('fake report', encoding='ascii')\n"
        "print('iqtree diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _fake_fasttree(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path('tree.newick').write_text('(a:0.1,b:0.2);', encoding='ascii')\n"
        "print('fasttree diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_infer_ml_tree_serializes_iqtree_artifacts(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(
        ">a\nACGT\n>b\nACGA\n",
        encoding="ascii",
    )
    executable = _fake_iqtree(tmp_path / "fake-iqtree2")

    result = infer_ml_tree(
        alignment,
        tmp_path / "tree",
        model="GTR+G",
        bootstrap=100,
        threads=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "iqtree2"
    treefile = result["result"]["treefile"]
    assert Path(treefile).read_text(encoding="ascii") == "(a:0.1,b:0.2);"
    assert "tree.iqtree" in result["result"]["artifacts"]


def test_infer_ml_tree_validates_bootstrap(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">a\nACGT\n>b\nACGA\n", encoding="ascii")
    with pytest.raises(ConfigurationError, match="bootstrap"):
        infer_ml_tree(alignment, tmp_path / "tree", bootstrap=-1, verbose=False)


def test_infer_fast_tree_captures_newick_output(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">a\nACGT\n>b\nACGA\n", encoding="ascii")
    executable = _fake_fasttree(tmp_path / "fake-fasttree")

    result = infer_fast_tree(
        alignment,
        tmp_path / "fasttree-output",
        nucleotide=True,
        gtr=True,
        gamma=True,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "fasttree"
    assert "tree.newick" in result["result"]["artifacts"]
    assert Path(result["result"]["treefile"]).read_text(encoding="ascii") == "(a:0.1,b:0.2);"


def test_iqtree_facade_is_registered_as_external_tool() -> None:
    registry = default_tool_registry(verbose=False)
    for name in ("dnakit.phylogeny.infer_ml_tree", "dnakit.phylogeny.infer_fast_tree"):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True
