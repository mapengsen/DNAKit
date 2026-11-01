"""Tests for the tskit gap modules: dnakit.population.tskit_gap_functions and tskit_cli_functions.

Offline only: no test uses the network. The tests need the optional tskit package
(dnakit[bio-tskit]) and are
skipped when it is not installed. Tree files are written into pytest's temporary directory.
"""

from __future__ import annotations

import importlib
import inspect
import io
import json
from pathlib import Path

import pytest

pytest.importorskip("tskit")
pytest.importorskip("numpy")

import tskit

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import tskit_cli_functions as cli
from dnakit.population import tskit_gap_functions as gap

TABLES = {
    "sequence_length": 100.0,
    "nodes": [
        {"time": 0.0, "sample": True},
        {"time": 0.0, "sample": True},
        {"time": 0.0, "sample": True},
        {"time": 1.0},
        {"time": 2.0},
    ],
    "edges": [
        {"left": 0.0, "right": 100.0, "parent": 3, "child": 0},
        {"left": 0.0, "right": 100.0, "parent": 3, "child": 1},
        {"left": 0.0, "right": 100.0, "parent": 4, "child": 2},
        {"left": 0.0, "right": 100.0, "parent": 4, "child": 3},
    ],
    "sites": [
        {"position": 10.0, "ancestral_state": "A"},
        {"position": 60.0, "ancestral_state": "C"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 2, "derived_state": "T"},
    ],
}


@pytest.fixture()
def trees_file(tmp_path: Path) -> Path:
    """A .trees file with the same content as TABLES, plus one provenance record."""

    tables = tskit.TableCollection(100.0)
    for time, sample in [(0, 1), (0, 1), (0, 1), (1, 0), (2, 0)]:
        tables.nodes.add_row(time=time, flags=tskit.NODE_IS_SAMPLE if sample else 0)
    tables.edges.add_row(0, 100, 3, 0)
    tables.edges.add_row(0, 100, 3, 1)
    tables.edges.add_row(0, 100, 4, 2)
    tables.edges.add_row(0, 100, 4, 3)
    tables.sites.add_row(10.0, "A")
    tables.sites.add_row(60.0, "C")
    tables.mutations.add_row(0, 0, "G")
    tables.mutations.add_row(1, 2, "T")
    tables.provenances.add_row(
        record=json.dumps({"software": {"name": "dnakit-test"}}), timestamp="2020"
    )
    path = tmp_path / "toy.trees"
    tables.tree_sequence().dump(str(path))
    return path


def test_public_functions_are_verbose_keyword_only() -> None:
    for module in (gap, cli):
        for name in module.__all__:
            function = getattr(module, name)
            verbose = inspect.signature(function).parameters["verbose"]
            assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
            assert verbose.default is True
            assert function.__doc__


# --- random nucleotides ------------------------------------------------------------------------


def test_random_nucleotides_match_tskit_and_repeat_for_a_seed() -> None:
    first = gap.tskit_random_nucleotides(length=40, seed=7, verbose=False)
    second = gap.tskit_random_nucleotides(length=40, seed=7, verbose=False)
    assert first["sequence"] == second["sequence"] == tskit.random_nucleotides(40, seed=7)
    assert first["length"] == 40
    assert set(first["sequence"]) <= set("ACGT")
    assert gap.tskit_random_nucleotides(length=0, seed=1, verbose=False)["sequence"] == ""


def test_random_nucleotides_checks_length_and_seed() -> None:
    with pytest.raises(ConfigurationError):
        gap.tskit_random_nucleotides(length=-1, verbose=False)
    with pytest.raises(ConfigurationError):
        gap.tskit_random_nucleotides(length=5, seed=-3, verbose=False)


# --- SVG drawings ------------------------------------------------------------------------------


def test_draw_tree_sequence_svg_whole_sequence_and_one_tree(tmp_path: Path) -> None:
    whole = gap.draw_tree_sequence_svg(TABLES, tmp_path / "whole.svg", verbose=False)
    assert whole["scope"] == "tree_sequence"
    assert whole["num_trees"] == 1
    assert whole["intervals"] == [[0.0, 100.0]]
    assert whole["num_sites"] == 2
    assert "<svg" in (tmp_path / "whole.svg").read_text(encoding="utf-8")

    one = gap.draw_tree_sequence_svg(TABLES, tmp_path / "tree.svg", tree_index=0, verbose=False)
    assert one["scope"] == "tree"
    assert one["tree_index"] == 0
    assert (tmp_path / "tree.svg").read_text(encoding="utf-8") != (
        tmp_path / "whole.svg"
    ).read_text(encoding="utf-8")


def test_draw_tree_sequence_svg_x_lim_and_checks(tmp_path: Path) -> None:
    limited = gap.draw_tree_sequence_svg(
        TABLES, tmp_path / "xlim.svg", x_lim=[0, 50], verbose=False
    )
    assert limited["x_lim"] == [0.0, 50.0]
    with pytest.raises(ConfigurationError, match="0 <= left < right"):
        gap.draw_tree_sequence_svg(TABLES, tmp_path / "bad.svg", x_lim=[10, 500], verbose=False)
    with pytest.raises(ConfigurationError, match="whole tree sequence"):
        gap.draw_tree_sequence_svg(
            TABLES, tmp_path / "bad2.svg", tree_index=0, x_lim=[0, 5], verbose=False
        )
    with pytest.raises(ConfigurationError, match="tree_index"):
        gap.draw_tree_sequence_svg(TABLES, tmp_path / "bad3.svg", tree_index=3, verbose=False)
    existing = tmp_path / "exists.svg"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        gap.draw_tree_sequence_svg(TABLES, existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep"


# --- HTML views --------------------------------------------------------------------------------


def test_draw_tree_sequence_html_writes_tskit_summary(tmp_path: Path) -> None:
    result = gap.draw_tree_sequence_html(TABLES, tmp_path / "summary.html", verbose=False)
    text = (tmp_path / "summary.html").read_text(encoding="utf-8")
    assert result["format"] == "html"
    assert result["num_trees"] == 1
    assert result["num_samples"] == 3
    assert result["bytes"] == len(text.encode("utf-8"))
    assert "<table" in text
    assert "<div" in text


def test_draw_tree_html_and_variant_html(tmp_path: Path) -> None:
    tree = gap.draw_tree_html(TABLES, 0, tmp_path / "tree.html", verbose=False)
    assert tree["interval"] == [0.0, 100.0]
    assert tree["num_roots"] == 1
    assert "<" in (tmp_path / "tree.html").read_text(encoding="utf-8")
    with pytest.raises(ConfigurationError, match="tree_index"):
        gap.draw_tree_html(TABLES, 1, tmp_path / "bad.html", verbose=False)

    second = gap.draw_variant_html(TABLES, 1, tmp_path / "variant.html", verbose=False)
    assert second["position"] == 60.0
    assert second["alleles"] == ["C", "T"]
    assert second["site_index"] == 1
    first = gap.draw_variant_html(TABLES, 0, tmp_path / "variant0.html", verbose=False)
    assert first["alleles"] == ["A", "G"]
    with pytest.raises(ConfigurationError, match="site_index"):
        gap.draw_variant_html(TABLES, 2, tmp_path / "bad_variant.html", verbose=False)


# --- tree file facades -------------------------------------------------------------------------


def test_tree_file_info_matches_the_table_counts(trees_file: Path) -> None:
    result = cli.tree_file_info(trees_file, verbose=False)
    assert result["num_trees"] == 1
    assert result["num_samples"] == 3
    assert result["num_nodes"] == 5
    assert result["num_edges"] == 4
    assert result["num_sites"] == 2
    assert result["num_mutations"] == 2
    assert result["num_provenances"] == 1
    assert result["sequence_length"] == 100.0
    assert "TreeSequence" in result["summary_text"]


def test_tree_file_trees_lists_intervals_and_draws_on_request(trees_file: Path) -> None:
    result = cli.tree_file_trees(trees_file, verbose=False)
    assert result["num_trees"] == 1
    assert result["trees"][0] == {"index": 0, "num_sites": 2, "interval": [0.0, 100.0]}
    assert "drawing" not in result["trees"][0]
    drawn = cli.tree_file_trees(trees_file, draw=True, precision=2, verbose=False)
    assert "drawing" in drawn["trees"][0]
    assert "0" in drawn["trees"][0]["drawing"]
    with pytest.raises(ConfigurationError, match="precision"):
        cli.tree_file_trees(trees_file, precision=40, verbose=False)


@pytest.mark.parametrize(
    ("function", "table", "rows"),
    [
        (cli.export_tree_file_nodes_table, "nodes", 5),
        (cli.export_tree_file_edges_table, "edges", 4),
        (cli.export_tree_file_sites_table, "sites", 2),
        (cli.export_tree_file_mutations_table, "mutations", 2),
        (cli.export_tree_file_individuals_table, "individuals", 0),
        (cli.export_tree_file_populations_table, "populations", 0),
        (cli.export_tree_file_migrations_table, "migrations", 0),
        (cli.export_tree_file_provenances_table, "provenances", 1),
    ],
)
def test_table_exports_write_tskit_text(
    trees_file: Path, tmp_path: Path, function, table: str, rows: int
) -> None:
    target = tmp_path / f"{table}.txt"
    result = function(trees_file, target, verbose=False)
    expected = io.StringIO()
    # tskit's populations and provenances dumps take no precision argument (as in tskit's own CLI)
    dump_options = {} if table in ("populations", "provenances") else {"precision": 6}
    tskit.load(str(trees_file)).dump_text(**{table: expected}, **dump_options)
    assert target.read_text(encoding="utf-8") == expected.getvalue()
    assert result["row_count"] == rows
    assert result["bytes"] == target.stat().st_size
    assert result["destination"] == str(target)


def test_nodes_table_has_the_expected_rows(trees_file: Path, tmp_path: Path) -> None:
    target = tmp_path / "nodes.txt"
    cli.export_tree_file_nodes_table(trees_file, target, precision=3, verbose=False)
    lines = target.read_text(encoding="utf-8").splitlines()
    assert lines[1].split("\t")[:3] == ["0", "1", "0.000"]
    assert lines[4].split("\t")[:3] == ["3", "0", "1.000"]


def test_provenance_human_readable_block(trees_file: Path, tmp_path: Path) -> None:
    target = tmp_path / "prov.txt"
    result = cli.export_tree_file_provenances_table(
        trees_file, target, human_readable=True, verbose=False
    )
    text = target.read_text(encoding="utf-8")
    assert result["row_count"] == 1
    assert text.startswith("id=0, timestamp=2020, record={")
    assert '"name": "dnakit-test"' in text


def test_fasta_and_vcf_exports(trees_file: Path, tmp_path: Path) -> None:
    fasta = cli.export_tree_file_fasta(trees_file, tmp_path / "out.fa", verbose=False)
    lines = (tmp_path / "out.fa").read_text(encoding="utf-8").splitlines()
    assert fasta["num_samples"] == 3
    assert sum(line.startswith(">") for line in lines) == 3
    one_line = cli.export_tree_file_fasta(
        trees_file, tmp_path / "flat.fa", wrap_width=0, verbose=False
    )
    flat = (tmp_path / "flat.fa").read_text(encoding="utf-8").splitlines()
    assert len(flat) == 6
    assert one_line["wrap_width"] == 0

    vcf = cli.export_tree_file_vcf(trees_file, tmp_path / "out.vcf", verbose=False)
    text = (tmp_path / "out.vcf").read_text(encoding="utf-8")
    assert text.startswith("##fileformat=VCFv4.2")
    assert vcf["num_records"] == 2
    with pytest.raises(ConfigurationError, match="multiple of ploidy"):
        cli.export_tree_file_vcf(trees_file, tmp_path / "bad.vcf", ploidy=2, verbose=False)


def test_tree_file_functions_check_their_input(trees_file: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="does not exist"):
        cli.tree_file_info(tmp_path / "missing.trees", verbose=False)
    text_file = tmp_path / "not_a_tree.trees"
    text_file.write_text("not a tree sequence\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="could not read"):
        cli.tree_file_info(text_file, verbose=False)
    with pytest.raises(ConfigurationError, match="above max_bytes"):
        cli.tree_file_info(trees_file, max_bytes=10, verbose=False)
    existing = tmp_path / "nodes.txt"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        cli.export_tree_file_nodes_table(trees_file, existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep"


def test_missing_tskit_raises_backend_unavailable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    real_import = importlib.import_module

    def fake_import(name: str, *args: object, **kwargs: object) -> object:
        if name == "tskit" or name.startswith("tskit."):
            raise ImportError("simulated missing tskit")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    with pytest.raises(BackendUnavailableError, match=r"dnakit\[bio-tskit\]"):
        gap.draw_tree_sequence_svg(TABLES, tmp_path / "x.svg", verbose=False)
    with pytest.raises(BackendUnavailableError, match=r"dnakit\[bio-tskit\]"):
        cli.tree_file_info(tmp_path / "x.trees", verbose=False)
