"""Tests for dnakit.phylogeny.toytree_gap_functions (toytree drawings, node data, prune and root).

Offline only: no test uses the network. The tests need the optional toytree package
(dnakit[bio-toytree]) and
are skipped when it is not installed.
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

pytest.importorskip("toytree")
pytest.importorskip("toyplot")
pytest.importorskip("pandas")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.phylogeny import toytree_gap_functions as tg

NEWICK = "((A:1.0,B:2.0)n1:0.5,(C:1.5,D:0.7)n2:0.8,E:3.0);"
PLAIN_NEWICK = "((A:1,B:2):0.5,(C:1.5,D:0.7):0.8,E:3);"
PIES = [[0.5, 0.5]] * 7 + [[0.2, 0.8]]  # eight node rows (nnodes = 8), each summing to 1

ANNOTATIONS = [
    "tip_labels",
    "tip_markers",
    "node_markers",
    "edge_markers",
    "node_labels",
    "edge_labels",
    "node_bars",
    "scale_bar",
    "box_outline",
    "box_outline_spines",
    "external_ticks",
]


def _svg(result: dict) -> str:
    return Path(result["result"]["image_path"]).read_text(encoding="utf-8")


def test_public_functions_are_verbose_keyword_only() -> None:
    for name in tg.__all__:
        function = getattr(tg, name)
        verbose = inspect.signature(function).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True
        assert function.__doc__


def test_draw_tree_figure_writes_each_image_format(tmp_path: Path) -> None:
    for fmt in ("svg", "png", "pdf", "html"):
        target = tmp_path / f"tree.{fmt}"
        result = tg.draw_tree_figure(NEWICK, target, image_format=fmt, verbose=False)
        payload = result["result"]
        assert result["status"] == "completed"
        assert Path(payload["image_path"]) == target
        assert target.stat().st_size == payload["image_bytes"] > 0
        assert payload["newick"] == PLAIN_NEWICK
        assert payload["tip_labels"] == ["A", "B", "C", "D", "E"]
        assert payload["tip_count"] == 5
        assert len(payload["nodes"]) == 8
        assert sum(node["is_tip"] for node in payload["nodes"]) == 5
        assert len(payload["tip_label_angles"]) == 5


def test_draw_tree_figure_reads_a_newick_file(tmp_path: Path) -> None:
    source = tmp_path / "tree.nwk"
    source.write_text(NEWICK + "\n", encoding="utf-8")
    result = tg.draw_tree_figure(source, tmp_path / "from_file.svg", verbose=False)
    assert result["result"]["tip_labels"] == ["A", "B", "C", "D", "E"]


@pytest.mark.parametrize("annotation", ANNOTATIONS)
def test_each_annotation_changes_the_drawing(tmp_path: Path, annotation: str) -> None:
    base = tg.draw_tree_figure(NEWICK, tmp_path / "base.svg", verbose=False)
    kwargs: dict[str, object] = {"annotations": [annotation]}
    if annotation == "node_bars":
        kwargs["feature"] = "dist"
    drawn = tg.draw_tree_figure(NEWICK, tmp_path / f"{annotation}.svg", verbose=False, **kwargs)
    assert drawn["result"]["annotations"] == [annotation]
    assert _svg(drawn) != _svg(base), f"annotation {annotation} left the drawing unchanged"


def test_tip_labels_annotation_uses_the_feature_as_label(tmp_path: Path) -> None:
    result = tg.draw_tree_figure(
        NEWICK, tmp_path / "labels.svg", annotations=["tip_labels"], feature="name", verbose=False
    )
    text = _svg(result)
    for label in ("A", "B", "C", "D", "E"):
        assert f">{label}<" in text


def test_pie_charts_draw_proportions_on_nodes_and_edges(tmp_path: Path) -> None:
    base = tg.draw_tree_figure(NEWICK, tmp_path / "base.svg", verbose=False)
    for annotation in ("node_pie_charts", "edge_pie_charts"):
        result = tg.draw_tree_figure(
            NEWICK,
            tmp_path / f"{annotation}.svg",
            annotations=[annotation],
            pie_proportions=PIES,
            verbose=False,
        )
        assert _svg(result) != _svg(base)
        assert result["result"]["annotations"] == [annotation]


@pytest.mark.parametrize(
    ("annotations", "pie", "message"),
    [
        (["node_pie_charts"], None, "need pie_proportions"),
        ([], PIES, "only used with node_pie_charts"),
        (["node_pie_charts"], [[0.5, 0.5]] * 3, "needs 8 rows"),
        (["node_pie_charts"], [[0.5, 0.4]] * 7 + [[0.2, 0.8]], "must sum to 1"),
        (["node_pie_charts"], [[1.5, -0.5]] * 8, "finite number in \\[0.0, 1.0\\]"),
        (["node_pie_charts"], [[1.0]] * 8, "at least two numbers"),
    ],
)
def test_pie_proportions_are_checked(
    tmp_path: Path, annotations: list[str], pie, message: str
) -> None:
    with pytest.raises(ConfigurationError, match=message):
        tg.draw_tree_figure(
            NEWICK,
            tmp_path / "bad.svg",
            annotations=annotations,
            pie_proportions=pie,
            verbose=False,
        )


def test_node_bars_need_a_feature(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="node_bars needs feature"):
        tg.draw_tree_figure(NEWICK, tmp_path / "bars.svg", annotations=["node_bars"], verbose=False)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"annotations": ["unknown_mark"]}, "annotations item must be one of"),
        ({"annotations": ["scale_bar", "scale_bar"]}, "must not repeat"),
        ({"layout": "sideways"}, "layout must be one of"),
        ({"tree_style": "z"}, "tree_style must be one of"),
        ({"image_format": "gif"}, "image_format must be one of"),
        ({"width": 10}, "width must be an integer"),
        ({"max_tips": 3}, "max_tips"),
    ],
)
def test_draw_tree_figure_rejects_bad_options(tmp_path: Path, kwargs: dict, message: str) -> None:
    with pytest.raises(ConfigurationError, match=message):
        tg.draw_tree_figure(NEWICK, tmp_path / "bad.svg", verbose=False, **kwargs)


def test_draw_tree_figure_never_overwrites_or_mismatches_suffix(tmp_path: Path) -> None:
    existing = tmp_path / "exists.svg"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        tg.draw_tree_figure(NEWICK, existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep"
    with pytest.raises(ConfigurationError, match=r"must end with \.png"):
        tg.draw_tree_figure(NEWICK, tmp_path / "tree.svg", image_format="png", verbose=False)
    with pytest.raises(ConfigurationError, match="Output directory does not exist"):
        tg.draw_tree_figure(NEWICK, tmp_path / "missing" / "tree.svg", verbose=False)


def test_draw_tree_node_feature_figure_colours_and_sizes(tmp_path: Path) -> None:
    colour = tg.draw_tree_node_feature_figure(
        NEWICK, "dist", tmp_path / "colour.svg", verbose=False
    )
    payload = colour["result"]
    assert payload["encoding"] == "color"
    assert payload["feature"] == "dist"
    assert len(payload["mapped"]) == len(payload["nodes"]) == 8
    assert all(value.startswith("#") and len(value) == 7 for value in payload["mapped"])
    assert payload["feature_values"][0] == pytest.approx(1.0)  # tip A has dist 1.0
    assert payload["feature_values"][4] == pytest.approx(3.0)  # tip E has dist 3.0
    assert colour["result"]["image_path"].endswith("colour.svg")

    size = tg.draw_tree_node_feature_figure(
        NEWICK,
        "dist",
        tmp_path / "size.svg",
        encoding="size",
        min_size=2.0,
        max_size=9.0,
        verbose=False,
    )
    sizes = size["result"]["mapped"]
    assert len(sizes) == 8
    assert min(sizes) == pytest.approx(2.0)
    assert max(sizes) == pytest.approx(9.0)


def test_draw_tree_node_feature_figure_rejects_unknown_feature(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="is not a node feature"):
        tg.draw_tree_node_feature_figure(
            NEWICK, "no_such_feature", tmp_path / "x.svg", verbose=False
        )
    with pytest.raises(ConfigurationError, match="encoding must be one of"):
        tg.draw_tree_node_feature_figure(
            NEWICK, "dist", tmp_path / "y.svg", encoding="shape", verbose=False
        )


def test_map_node_feature_values_maps_colours_and_sizes() -> None:
    colours = tg.map_node_feature_values([0.0, 0.5, 1.0, None], verbose=False)
    assert colours["result"]["encoding"] == "color"
    assert colours["result"]["count"] == 4
    assert all(value.startswith("#") for value in colours["result"]["mapped"])
    assert colours["result"]["mapped"][0] != colours["result"]["mapped"][2]

    sizes = tg.map_node_feature_values(
        [1, 2, 3], encoding="size", min_size=5, max_size=15, verbose=False
    )
    assert sizes["result"]["mapped"] == pytest.approx([5.0, 10.0, 15.0])


@pytest.mark.parametrize(
    ("values", "kwargs", "message"),
    [
        ([], {}, "non-empty list"),
        ("1,2", {}, "non-empty list"),
        ([1.0, "x"], {}, "values item must be a number"),
        ([1.0, 2.0], {"encoding": "shape"}, "encoding must be one of"),
        ([1.0, 2.0], {"min_size": 9.0, "max_size": 1.0}, "max_size must be at least min_size"),
    ],
)
def test_map_node_feature_values_rejects_bad_input(values, kwargs: dict, message: str) -> None:
    with pytest.raises(ConfigurationError, match=message):
        tg.map_node_feature_values(values, verbose=False, **kwargs)


def test_draw_tree_cloud_and_grid_figures(tmp_path: Path) -> None:
    cloud = tg.draw_tree_cloud_figure([NEWICK, NEWICK], tmp_path / "cloud.svg", verbose=False)
    assert cloud["result"]["tree_count"] == 2
    assert len(cloud["result"]["newick"]) == 2
    assert cloud["result"]["tip_count"] == 5
    assert (tmp_path / "cloud.svg").stat().st_size > 0

    grid = tg.draw_tree_grid_figure(
        [NEWICK, NEWICK], tmp_path / "grid.svg", columns=2, verbose=False
    )
    assert grid["result"]["shape"] == [1, 2]
    assert grid["result"]["tree_count"] == 2
    assert [tree["tip_labels"] for tree in grid["result"]["trees"]] == [
        ["A", "B", "C", "D", "E"]
    ] * 2


def test_multitree_figures_reject_mismatched_tips_and_small_grids(tmp_path: Path) -> None:
    other = "((A:1,B:2):0.5,(C:1.5,X:0.7):0.8,E:3);"
    with pytest.raises(ConfigurationError, match="same tip labels"):
        tg.draw_tree_cloud_figure([NEWICK, other], tmp_path / "bad_cloud.svg", verbose=False)
    with pytest.raises(ConfigurationError, match="do not fit"):
        tg.draw_tree_grid_figure(
            [NEWICK] * 3, tmp_path / "bad_grid.svg", rows=1, columns=2, verbose=False
        )


def test_export_tree_node_data_returns_the_node_table(tmp_path: Path) -> None:
    result = tg.export_tree_node_data(NEWICK, verbose=False)["result"]
    assert result["columns"] == ["idx", "name", "height", "dist", "support"]
    assert result["row_count"] == 8
    assert [row["idx"] for row in result["rows"]] == list(range(8))
    assert result["rows"][0]["name"] == "A"
    assert result["rows"][0]["support"] is None  # tips have no support value

    tips = tg.export_tree_node_data(NEWICK, tips_only=True, verbose=False)["result"]
    assert tips["row_count"] == 5
    assert {row["name"] for row in tips["rows"]} == {"A", "B", "C", "D", "E"}


def test_export_tree_node_data_selects_features_and_nodes_in_order(tmp_path: Path) -> None:
    target = tmp_path / "data.tsv"
    result = tg.export_tree_node_data(
        NEWICK,
        features=["dist"],
        nodes=["E", "n1"],
        output_path=target,
        verbose=False,
    )["result"]
    assert result["columns"] == ["idx", "dist"]
    assert [row["idx"] for row in result["rows"]] == [4, 5]
    assert result["rows"][0]["dist"] == pytest.approx(3.0)
    lines = target.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "idx\tdist"
    assert lines[1] == "4\t3.0"
    assert result["output_bytes"] == target.stat().st_size


def test_export_tree_node_data_checks_its_arguments(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="not found among the node names"):
        tg.export_tree_node_data(NEWICK, nodes=["ZZ"], verbose=False)
    with pytest.raises(ConfigurationError, match="must be tips"):
        tg.export_tree_node_data(NEWICK, tips_only=True, nodes=["n1"], verbose=False)
    with pytest.raises(ConfigurationError, match="single character"):
        tg.export_tree_node_data(NEWICK, separator=";;", verbose=False)
    existing = tmp_path / "exists.tsv"
    existing.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        tg.export_tree_node_data(NEWICK, output_path=existing, verbose=False)


def test_prune_tree_tips_keeps_the_connecting_subtree() -> None:
    result = tg.prune_tree_tips(NEWICK, ["A", "C", "E"], verbose=False)["result"]
    assert result["tips"] == ["A", "C", "E"]
    assert result["tip_count"] == 3
    assert result["retained_names"] == ["A", "C", "E"]
    # with preserve_dists the path length A-n1-root is kept: 1.0 + 0.5
    assert "A:1.5" in result["newick"]


def test_prune_tree_tips_checks_names_and_count() -> None:
    with pytest.raises(ConfigurationError, match="at least 2 node name"):
        tg.prune_tree_tips(NEWICK, ["A"], verbose=False)
    with pytest.raises(ConfigurationError, match="not found among the node names"):
        tg.prune_tree_tips(NEWICK, ["A", "ZZ"], verbose=False)


def test_root_tree_on_outgroup_places_the_root_on_the_outgroup_edge() -> None:
    toytree = pytest.importorskip("toytree")
    result = tg.root_tree_on_outgroup(NEWICK, ["E"], verbose=False)["result"]
    assert result["tips"] == ["A", "B", "C", "D", "E"]
    rooted = toytree.tree(result["newick"])
    children = [child.name for child in rooted.treenode.children]
    assert "E" in children
    assert len(children) == 2
    assert result["outgroup"] == ["E"]


def test_root_tree_on_outgroup_with_several_names_and_checks() -> None:
    result = tg.root_tree_on_outgroup(NEWICK, ["A", "B"], edge_features=["dist"], verbose=False)[
        "result"
    ]
    assert result["tips"] == ["A", "B", "C", "D", "E"]
    assert result["edge_features"] == ["dist"]
    with pytest.raises(ConfigurationError, match="not found among the node names"):
        tg.root_tree_on_outgroup(NEWICK, ["ZZ"], verbose=False)
    with pytest.raises(ConfigurationError, match="at least 1 node name"):
        tg.root_tree_on_outgroup(NEWICK, [], verbose=False)


def test_missing_toytree_raises_backend_unavailable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    real_import = importlib.import_module

    def fake_import(name: str, *args: object, **kwargs: object) -> object:
        if name == "toytree" or name.startswith("toytree."):
            raise ImportError("simulated missing toytree")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    with pytest.raises(BackendUnavailableError, match=r"dnakit\[bio-toytree\]"):
        tg.draw_tree_figure(NEWICK, tmp_path / "x.svg", verbose=False)
    with pytest.raises(BackendUnavailableError, match=r"dnakit\[bio-toytree\]"):
        tg.prune_tree_tips(NEWICK, ["A", "B"], verbose=False)
