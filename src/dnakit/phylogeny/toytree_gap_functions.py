"""Toytree drawing, node-data and tree-editing functions for the coverage gap.

- ``draw_tree_figure()`` draws one tree with ``toytree.ToyTree.draw`` (the routine behind
  ``toytree.draw``, ``toytree.drawing.draw_toytree`` and ``toytree.core.tree.draw_toytree``). It
  adds the optional annotations (tip, node and edge labels and markers, node bars, pie charts, scale
  bar, axes outlines, external ticks) with ``toytree.annotate`` and ``toytree.utils``, and writes
  the figure with ``toytree.save``.
- ``draw_tree_node_feature_figure()`` colours or sizes nodes by a node feature with
  ``toytree.style.get_color_mapped_feature()`` or ``get_range_mapped_feature()``.
- ``map_node_feature_values()`` applies the same colour or size mapping to a plain list of values.
- ``draw_tree_cloud_figure()`` and ``draw_tree_grid_figure()`` draw several trees with toytree's
  multi-tree API (``MultiTree.draw_cloud_tree`` and ``MultiTree.draw``).
- ``export_tree_node_data()`` is the operation of toytree's ``get-node-data`` command: the node (or
  tip) feature table, optionally written as delimited text.
- ``prune_tree_tips()`` and ``root_tree_on_outgroup()`` are the operations of toytree's ``prune``
  and ``root`` commands (``toytree.mod.prune`` and ``toytree.mod.root``).

Every drawing returns the image path and the data it was drawn from: the Newick string, tip labels
and per-node table with layout coordinates. No plotting object is returned. Output is SVG, PNG
(toyplot's pure-Python renderer, no browser), HTML or PDF, written to a new file at an explicit
path. Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-toytree]"`` hint.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.phylogeny.toytree_functions import _newick_result, _read_tree, _tip_labels

_EXTRA = "bio-toytree"
_LAYOUTS = ("r", "l", "u", "d", "c", "c0-180", "un")
_TREE_STYLES = ("n", "r", "c", "s", "o", "b")
_EDGE_TYPES = ("p", "c", "b")
_IMAGE_FORMATS = ("svg", "png", "html", "pdf")
_ENCODINGS = ("color", "size")
_ANNOTATIONS = (
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
    "node_pie_charts",
    "edge_pie_charts",
)
_PIE_ANNOTATIONS = ("node_pie_charts", "edge_pie_charts")
_NODE_TABLE_COLUMNS = ("idx", "name", "height", "dist", "support")
_DRAW_CONVENTION = (
    "Coordinates are toytree's layout coordinates (mark.ntable, one row per node index). "
    "Tip labels are the sorted tip names. Images are written by toytree.save from toyplot; "
    "the file must not exist yet."
)


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def _int(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _number(value: object, name: str, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{low}, {high}].")
    number = float(value)
    if not math.isfinite(number) or not low <= number <= high:
        raise ConfigurationError(f"{name} must be a finite number in [{low}, {high}].")
    return number


def _optional_number(value: object, name: str) -> float | None:
    if value is None:
        return None
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        raise ConfigurationError(f"{name} must be a finite number or null.")
    return float(value)


def _feature(value: object, name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip() or len(value) > 200 or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string of at most 200 characters.")
    return value


def _annotations(values: object) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise ConfigurationError("annotations must be a list of annotation names.")
    names = tuple(_choice(item, "annotations item", _ANNOTATIONS) for item in values)
    if len(set(names)) != len(names):
        raise ConfigurationError("annotations must not repeat a name.")
    return names


def _package() -> str:
    return f"toytree=={metadata.version('toytree')}"


def _new_file(output_path: object) -> Path:
    if not isinstance(output_path, (str, Path)) or "\x00" in str(output_path):
        raise ConfigurationError("output_path must be a file path.")
    target = Path(output_path).expanduser().resolve()
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    if target.exists():
        raise ConfigurationError(
            f"Output file already exists and will not be overwritten: {target}"
        )
    return target


def _new_output(output_path: object, image_format: str) -> Path:
    if (
        isinstance(output_path, (str, Path))
        and Path(output_path).suffix.lower() != f".{image_format}"
    ):
        raise ConfigurationError(
            f"output_path must end with .{image_format} for image_format={image_format}."
        )
    return _new_file(output_path)


def _names(values: object, name: str, *, minimum: int = 1) -> list[str]:
    """A list of distinct, non-empty node names (tip or internal names)."""
    if isinstance(values, str) or not isinstance(values, (list, tuple)) or len(values) < minimum:
        raise ConfigurationError(f"{name} must be a list of at least {minimum} node name(s).")
    names = [_feature(item, f"{name} item") for item in values]
    if any(item is None for item in names):
        raise ConfigurationError(f"{name} items must be non-empty strings.")
    unique = [str(item) for item in names if item is not None]
    if len(set(unique)) != len(unique):
        raise ConfigurationError(f"{name} must not repeat a node name.")
    return unique


def _feature_list(values: object) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, str) or not isinstance(values, (list, tuple)):
        raise ConfigurationError("features must be a list of node feature names or null.")
    result = [_feature(item, "features item") for item in values]
    if any(item is None for item in result):
        raise ConfigurationError("features items must be non-empty strings.")
    return [str(item) for item in result if item is not None]


def _separator(value: object) -> str:
    if not isinstance(value, str) or len(value) != 1:
        raise ConfigurationError("separator must be a single character.")
    return value


def _node_frame(table: Any, requested: list[str] | None) -> Any:
    """The node table as a DataFrame with an idx column first (toytree returns a Series for one
    feature).
    """
    if hasattr(table, "to_frame"):
        label = table.name if table.name is not None else (requested or ["value"])[0]
        frame = table.to_frame(name=label)
    else:
        frame = table.copy()
    if "idx" not in frame.columns:
        frame.insert(0, "idx", frame.index.to_numpy())
    return frame.reset_index(drop=True)


def _known_node_names(tree: Any) -> set[str]:
    return {str(name) for name in tree.get_node_data()["name"].tolist() if str(name)}


def _check_names_in_tree(tree: Any, names: list[str], name: str) -> None:
    missing = sorted(set(names) - _known_node_names(tree))
    if missing:
        raise ConfigurationError(
            f"{name} not found among the node names of the tree: {', '.join(missing)}"
        )


def _pie_array(tree: Any, proportions: object) -> Any:
    """Pie proportions as an (nnodes, categories) array; each row sums to 1 (row = node index)."""
    numpy = require("numpy", _EXTRA)
    if isinstance(proportions, str) or not isinstance(proportions, (list, tuple)):
        raise ConfigurationError(
            "pie_proportions must be a list with one row of proportions per node."
        )
    if len(proportions) != tree.nnodes:
        raise ConfigurationError(
            f"pie_proportions needs {tree.nnodes} rows (one per node index), "
            f"got {len(proportions)}."
        )
    rows: list[list[float]] = []
    for row in proportions:
        if isinstance(row, str) or not isinstance(row, (list, tuple)) or len(row) < 2:
            raise ConfigurationError(
                "each pie_proportions row must be a list of at least two numbers."
            )
        rows.append([_number(value, "pie_proportions value", 0.0, 1.0) for value in row])
    if any(len(row) != len(rows[0]) for row in rows):
        raise ConfigurationError(
            "all pie_proportions rows must have the same number of categories."
        )
    array = numpy.array(rows, dtype=float)
    if not numpy.allclose(array.sum(axis=1), 1.0, atol=1e-6):
        raise ConfigurationError("each pie_proportions row must sum to 1.")
    return array


def _save(toytree: Any, canvas: Any, target: Path) -> int:
    try:
        toytree.save(canvas, str(target))
    except Exception as exc:
        raise BackendExecutionError(f"toytree could not write {target.name}: {exc}") from exc
    if not target.is_file() or target.stat().st_size == 0:
        raise BackendExecutionError(f"toytree did not write a non-empty image to {target}.")
    return target.stat().st_size


def _float_or_none(value: object) -> float | None:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _hex(color: Any) -> str:
    if getattr(color, "dtype", None) is not None and color.dtype.names:
        channels = (color["r"], color["g"], color["b"])
    else:
        channels = tuple(color[:3])
    red, green, blue = (round(float(channel) * 255) for channel in channels)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _node_rows(tree: Any, mark: Any) -> list[dict[str, object]]:
    table = tree.get_node_data()
    coords = mark.ntable
    rows: list[dict[str, object]] = []
    for position, record in enumerate(table.to_dict("records")):
        name = record["name"]
        support = _float_or_none(record["support"])
        rows.append(
            {
                "idx": int(record["idx"]),
                "name": name if isinstance(name, str) else "",
                "is_tip": position < tree.ntips,
                "height": _float_or_none(record["height"]),
                "dist": _float_or_none(record["dist"]),
                "support": support,
                "x": _float_or_none(coords[position][0]),
                "y": _float_or_none(coords[position][1]),
            }
        )
    return rows


def _draw_result(
    operation: str, target: Path, image_format: str, size: int, fields: dict[str, object]
) -> dict[str, object]:
    return {
        "status": "completed",
        "operation": operation,
        "result": {
            "image_path": str(target),
            "image_format": image_format,
            "image_bytes": size,
            **fields,
        },
        "provenance": {"package": _package(), "convention": _DRAW_CONVENTION},
    }


def _add_annotation(name: str, tree: Any, axes: Any, feature: str | None, pie_data: Any) -> None:
    annotate = require("toytree.annotate", _EXTRA)
    utils = require("toytree.utils", _EXTRA)
    if name == "node_pie_charts":
        annotate.add_node_pie_charts(tree, axes, data=pie_data)
    elif name == "edge_pie_charts":
        annotate.add_edge_pie_charts(tree, axes, data=pie_data)
    elif name == "tip_labels":
        annotate.add_tip_labels(tree, axes, labels=feature or "name")
    elif name == "tip_markers":
        annotate.add_tip_markers(tree, axes)
    elif name == "node_markers":
        annotate.add_node_markers(tree, axes)
    elif name == "edge_markers":
        annotate.add_edge_markers(tree, axes)
    elif name == "node_labels":
        annotate.add_node_labels(tree, axes, labels=feature or "idx")
    elif name == "edge_labels":
        annotate.add_edge_labels(tree, axes, labels=feature or "idx")
    elif name == "node_bars":
        annotate.add_node_bars(tree, axes, bar_min=0, bar_max=feature)
    elif name == "scale_bar":
        annotate.add_axes_scale_bar(tree, axes)
    elif name == "box_outline":
        annotate.add_axes_box_outline(tree, axes)
    elif name == "box_outline_spines":
        utils.set_axes_box_outline(axes)
    else:
        utils.set_axes_ticks_external(axes)


def draw_tree_figure(
    tree: object,
    output_path: str | Path,
    *,
    image_format: str = "svg",
    layout: str = "r",
    tree_style: str = "n",
    edge_type: str | None = None,
    width: int | None = None,
    height: int | None = None,
    tip_labels_align: bool = False,
    annotations: Sequence[str] = (),
    feature: str | None = None,
    pie_proportions: Sequence[Sequence[float]] | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Draw one tree with toytree to a new SVG, PNG, HTML or PDF file and return the image path and
    the drawn data.

    tree: Newick text or a path to a one-tree Newick file. layout: 'r', 'l', 'u', 'd', 'c', 'c0-180'
    or 'un'. tree_style: toytree base style 'n', 'r', 'c', 's', 'o' or 'b'. edge_type: None (toytree
    default, phylogram), 'p', 'c' (cladogram) or 'b' (bezier). annotations: any of tip_labels,
    tip_markers, node_markers, edge_markers, node_labels, edge_labels, node_bars (needs feature),
    scale_bar, box_outline, box_outline_spines, external_ticks, node_pie_charts and edge_pie_charts
    (both need pie_proportions); when 'tip_labels' is requested, the tip labels come from
    toytree.annotate.add_tip_labels (labelled by feature or the tip name) instead of tree.draw.
    feature: optional node feature used as label source (node/edge labels, tip labels) or as the bar
    height (node_bars). pie_proportions: one row per node index (the node table order, root last) of
    category proportions that sum to 1; all rows have the same number of categories. Result:
    image_path, image_format, image_bytes, the Newick string, the sorted tip labels, a per-node
    table (idx, name, height, dist, support, x, y) and tip label angles.
    """

    toytree = require("toytree", _EXTRA)
    fmt = _choice(image_format, "image_format", _IMAGE_FORMATS)
    layout_name = _choice(layout, "layout", _LAYOUTS)
    style_name = _choice(tree_style, "tree_style", _TREE_STYLES)
    edge = None if edge_type is None else _choice(edge_type, "edge_type", _EDGE_TYPES)
    names = _annotations(annotations)
    label_feature = _feature(feature, "feature")
    if "node_bars" in names and label_feature is None:
        raise ConfigurationError(
            "node_bars needs feature: the node feature that sets each bar height."
        )
    pies = [name for name in names if name in _PIE_ANNOTATIONS]
    if pies and pie_proportions is None:
        raise ConfigurationError(f"{' and '.join(pies)} need pie_proportions.")
    if pie_proportions is not None and not pies:
        raise ConfigurationError(
            "pie_proportions is only used with node_pie_charts or edge_pie_charts."
        )
    width_value = None if width is None else _int(width, "width", 50, 10_000)
    height_value = None if height is None else _int(height, "height", 50, 10_000)
    align = _flag(tip_labels_align, "tip_labels_align")
    target = _new_output(output_path, fmt)
    tree_obj = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    pie_data = None if not pies else _pie_array(tree_obj, pie_proportions)
    try:
        canvas, axes, mark = tree_obj.draw(
            layout=layout_name,
            tree_style=style_name,
            edge_type=edge,
            width=width_value,
            height=height_value,
            tip_labels_align=align,
            tip_labels="tip_labels" not in names,
        )
    except Exception as exc:
        raise ConfigurationError(f"toytree could not draw the tree: {exc}") from exc
    for name in names:
        _add_annotation(name, tree_obj, axes, label_feature, pie_data)
    size = _save(toytree, canvas, target)
    angles = require("toytree.layout", _EXTRA).get_tip_labels_angles(tree_obj, mark.ntable)
    return _draw_result(
        "draw_tree_figure",
        target,
        fmt,
        size,
        {
            "layout": layout_name,
            "tree_style": style_name,
            "edge_type": edge,
            "annotations": list(names),
            "tip_count": len(_tip_labels(tree_obj)),
            "tip_labels": _tip_labels(tree_obj),
            "newick": str(tree_obj.write()),
            "nodes": _node_rows(tree_obj, mark),
            "tip_label_angles": [_float_or_none(value) for value in angles],
        },
    )


def _feature_values(tree: Any, feature: str) -> list[float | None]:
    table = tree.get_node_data()
    if feature not in table.columns:
        raise ConfigurationError(f"feature {feature!r} is not a node feature of the tree.")
    return [_float_or_none(value) for value in table[feature].tolist()]


def draw_tree_node_feature_figure(
    tree: object,
    feature: str,
    output_path: str | Path,
    *,
    encoding: str = "color",
    image_format: str = "svg",
    layout: str = "r",
    tree_style: str = "n",
    tips_only: bool = False,
    cmap: str | None = None,
    domain_min: float | None = None,
    domain_max: float | None = None,
    nan_value: float | None = None,
    min_size: float = 5.0,
    max_size: float = 15.0,
    reverse: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Draw a tree with its nodes coloured or sized by a numeric node feature, and return the mapped
    values.

    encoding 'color' colours the node markers with toytree.style.get_color_mapped_feature (cmap,
    domain_min, domain_max, nan_value, tips_only, reverse). encoding 'size' sizes them with
    get_range_mapped_feature (min_size..max_size; nan_value defaults to 0 in toytree when not set).
    The feature must be a node data column such as 'dist', 'height' or 'support'. Result: the image
    path, the node table, the raw feature values and the mapped colours (hex) or sizes, one entry
    per node index.
    """

    toytree = require("toytree", _EXTRA)
    style = require("toytree.style", _EXTRA)
    fmt = _choice(image_format, "image_format", _IMAGE_FORMATS)
    mode = _choice(encoding, "encoding", _ENCODINGS)
    layout_name = _choice(layout, "layout", _LAYOUTS)
    style_name = _choice(tree_style, "tree_style", _TREE_STYLES)
    name = _feature(feature, "feature")
    if name is None:
        raise ConfigurationError("feature is required.")
    tips = _flag(tips_only, "tips_only")
    colormap = _feature(cmap, "cmap")
    low = _optional_number(domain_min, "domain_min")
    high = _optional_number(domain_max, "domain_max")
    fill = _optional_number(nan_value, "nan_value")
    lo_size = _number(min_size, "min_size", 0.0, 1000.0)
    hi_size = _number(max_size, "max_size", 0.0, 1000.0)
    if hi_size < lo_size:
        raise ConfigurationError("max_size must be at least min_size.")
    reversed_map = _flag(reverse, "reverse")
    target = _new_output(output_path, fmt)
    tree_obj = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    values = _feature_values(tree_obj, name)
    try:
        if mode == "color":
            mapped = style.get_color_mapped_feature(
                tree_obj,
                name,
                cmap=colormap,
                domain_min=low,
                domain_max=high,
                nan_value=fill,
                tips_only=tips,
                reverse=reversed_map,
            )
            mapped_out: list[object] = [_hex(color) for color in mapped]
            draw_kwargs: dict[str, object] = {"node_sizes": 8, "node_colors": list(mapped)}
        else:
            mapped = style.get_range_mapped_feature(
                tree_obj,
                name,
                min_value=lo_size,
                max_value=hi_size,
                nan_value=0.0 if fill is None else fill,
                tips_only=tips,
            )
            mapped_out = [float(value) for value in mapped]
            draw_kwargs = {"node_sizes": list(mapped)}
        canvas, _axes, mark = tree_obj.draw(
            layout=layout_name, tree_style=style_name, **draw_kwargs
        )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"toytree could not map feature {name!r}: {exc}") from exc
    size = _save(toytree, canvas, target)
    return _draw_result(
        "draw_tree_node_feature_figure",
        target,
        fmt,
        size,
        {
            "feature": name,
            "encoding": mode,
            "tips_only": tips,
            "tip_labels": _tip_labels(tree_obj),
            "newick": str(tree_obj.write()),
            "nodes": _node_rows(tree_obj, mark),
            "feature_values": values,
            "mapped": mapped_out,
        },
    )


def map_node_feature_values(
    values: Sequence[float | None],
    *,
    encoding: str = "color",
    cmap: str | None = None,
    domain_min: float | None = None,
    domain_max: float | None = None,
    nan_value: float | None = None,
    min_size: float = 5.0,
    max_size: float = 15.0,
    reverse: bool = False,
    verbose: bool = True,
) -> dict[str, object]:
    """Map a list of numbers (null allowed) to colours or sizes with toytree's value mappers, no
    drawing.

    encoding 'color' calls toytree.style.get_color_mapped_values and returns hex colours; 'size'
    calls get_range_mapped_values and returns sizes between min_size and max_size. null entries
    become NaN and are filled with nan_value (for 'size' the default fill is 0.0, as in toytree).
    Result: encoding, mapped values in input order.
    """

    style = require("toytree.style", _EXTRA)
    mode = _choice(encoding, "encoding", _ENCODINGS)
    if isinstance(values, str) or not isinstance(values, (list, tuple)) or not values:
        raise ConfigurationError("values must be a non-empty list of numbers or nulls.")
    if len(values) > 1_000_000:
        raise ConfigurationError("values can hold at most 1000000 entries.")
    numbers = [
        math.nan if item is None else _number(item, "values item", -1e300, 1e300) for item in values
    ]
    colormap = _feature(cmap, "cmap")
    low = _optional_number(domain_min, "domain_min")
    high = _optional_number(domain_max, "domain_max")
    fill = _optional_number(nan_value, "nan_value")
    lo_size = _number(min_size, "min_size", 0.0, 1000.0)
    hi_size = _number(max_size, "max_size", 0.0, 1000.0)
    if hi_size < lo_size:
        raise ConfigurationError("max_size must be at least min_size.")
    try:
        if mode == "color":
            mapped = style.get_color_mapped_values(
                numbers,
                cmap=colormap,
                domain_min=low,
                domain_max=high,
                nan_value=fill,
                reverse=_flag(reverse, "reverse"),
            )
            out: list[object] = [_hex(color) for color in mapped]
        else:
            mapped = style.get_range_mapped_values(
                numbers,
                min_value=lo_size,
                max_value=hi_size,
                nan_value=0.0 if fill is None else fill,
            )
            out = [float(value) for value in mapped]
    except Exception as exc:
        raise ConfigurationError(f"toytree could not map the values: {exc}") from exc
    return {
        "status": "completed",
        "operation": "map_node_feature_values",
        "result": {"encoding": mode, "count": len(out), "mapped": out},
        "provenance": {
            "package": _package(),
            "convention": (
                "Colours are the RGB values of toytree's mapper as #rrggbb; "
                "sizes are toytree's range map."
            ),
        },
    }


def _trees_for_multitree(
    toytree: Any, trees: object, max_tips: int, max_file_bytes: int
) -> tuple[list[Any], list[str]]:
    if isinstance(trees, str) or not isinstance(trees, (list, tuple)) or not 1 <= len(trees) <= 100:
        raise ConfigurationError(
            "trees must be a list of 1 to 100 Newick strings or tree file paths."
        )
    parsed = [
        _read_tree(toytree, item, max_tips=max_tips, max_file_bytes=max_file_bytes)
        for item in trees
    ]
    tips = _tip_labels(parsed[0])
    if any(_tip_labels(tree) != tips for tree in parsed[1:]):
        raise ConfigurationError("All trees must have the same tip labels.")
    return parsed, tips


def draw_tree_cloud_figure(
    trees: Sequence[str],
    output_path: str | Path,
    *,
    image_format: str = "svg",
    width: int | None = None,
    height: int | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Draw several trees with the same tips overlaid as one cloud (toytree
    MultiTree.draw_cloud_tree).

    trees: a list of Newick strings or tree file paths sharing tip labels. Only the first tree shows
    tip labels; the tip order follows the consensus tree (toytree's default). Result: image_path,
    tree_count, tip labels, the Newick string of each input tree and the consensus Newick used for
    the tip order.
    """

    toytree = require("toytree", _EXTRA)
    fmt = _choice(image_format, "image_format", _IMAGE_FORMATS)
    width_value = None if width is None else _int(width, "width", 50, 10_000)
    height_value = None if height is None else _int(height, "height", 50, 10_000)
    target = _new_output(output_path, fmt)
    parsed, tips = _trees_for_multitree(toytree, trees, max_tips, max_file_bytes)
    multi = toytree.mtree(parsed)
    try:
        canvas, _axes, _marks = multi.draw_cloud_tree(width=width_value, height=height_value)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not draw the tree cloud: {exc}") from exc
    size = _save(toytree, canvas, target)
    return _draw_result(
        "draw_tree_cloud_figure",
        target,
        fmt,
        size,
        {
            "tree_count": len(parsed),
            "tip_count": len(tips),
            "tip_labels": tips,
            "newick": [str(tree.write()) for tree in parsed],
            "consensus_newick": str(multi.get_consensus_tree().write()),
        },
    )


def draw_tree_grid_figure(
    trees: Sequence[str],
    output_path: str | Path,
    *,
    image_format: str = "svg",
    rows: int = 1,
    columns: int = 4,
    width: int | None = None,
    height: int | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Draw several trees side by side in a rows x columns grid (toytree MultiTree.draw, toytree's
    multitree layout).

    trees: a list of Newick strings or tree file paths; at most rows * columns trees. Each tree gets
    its own cell. Result: image_path, grid shape, tree_count, and for each tree its Newick string
    and tip labels.
    """

    toytree = require("toytree", _EXTRA)
    fmt = _choice(image_format, "image_format", _IMAGE_FORMATS)
    grid_rows = _int(rows, "rows", 1, 20)
    grid_columns = _int(columns, "columns", 1, 20)
    width_value = None if width is None else _int(width, "width", 50, 10_000)
    height_value = None if height is None else _int(height, "height", 50, 10_000)
    target = _new_output(output_path, fmt)
    parsed, _ = _trees_for_multitree(toytree, trees, max_tips, max_file_bytes)
    if len(parsed) > grid_rows * grid_columns:
        raise ConfigurationError(
            f"{len(parsed)} trees do not fit a {grid_rows}x{grid_columns} grid."
        )
    multi = toytree.mtree(parsed)
    try:
        canvas = multi.draw(
            shape=(grid_rows, grid_columns), width=width_value, height=height_value
        )[0]
    except Exception as exc:
        raise ConfigurationError(f"toytree could not draw the tree grid: {exc}") from exc
    size = _save(toytree, canvas, target)
    return _draw_result(
        "draw_tree_grid_figure",
        target,
        fmt,
        size,
        {
            "shape": [grid_rows, grid_columns],
            "tree_count": len(parsed),
            "trees": [
                {"newick": str(tree.write()), "tip_labels": _tip_labels(tree)} for tree in parsed
            ],
        },
    )


def export_tree_node_data(
    tree: object,
    *,
    features: Sequence[str] | None = None,
    nodes: Sequence[str] | None = None,
    tips_only: bool = False,
    separator: str = "\t",
    output_path: str | Path | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Return the node data table of a tree (toytree get-node-data), optionally writing it to a new
    delimited file.

    tree: Newick text or a path to a one-tree Newick file. features: node feature names to keep
    (null keeps all of idx, name, height, dist, support and any other feature). nodes: node names to
    select, returned in the given order (null keeps every node); a name that is not a node of the
    tree is an error. tips_only: keep only tips, as toytree's get_tip_data does. separator: one
    character for the optional file (tab by default, as toytree's command). output_path: an optional
    new file; it has the idx column first and no pandas index. Result: the columns, the rows as
    dicts (NaN becomes null), the row count and the file path if one was written.
    """

    feature_names = _feature_list(features)
    selected = None if nodes is None else _names(nodes, "nodes")
    tips = _flag(tips_only, "tips_only")
    sep = _separator(separator)
    target = None if output_path is None else _new_file(output_path)
    toytree = require("toytree", _EXTRA)
    tree_obj = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    if selected is not None:
        _check_names_in_tree(tree_obj, selected, "nodes")
    try:
        raw = (
            tree_obj.get_tip_data(feature_names) if tips else tree_obj.get_node_data(feature_names)
        )
        table = _node_frame(raw, feature_names)
        if selected is not None:
            # one query per name: rows follow the order of the names, not toytree's own order
            wanted = [int(node.idx) for name in selected for node in tree_obj.get_nodes(name)]
            if not set(wanted) <= set(table["idx"].tolist()):
                raise ConfigurationError("nodes must be tips when tips_only is true.")
            table = table.set_index("idx", drop=False).loc[wanted].reset_index(drop=True)
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"toytree could not read the node data: {exc}") from exc
    columns = [str(column) for column in table.columns]
    rows: list[dict[str, object]] = []
    for record in table.to_dict("records"):
        row: dict[str, object] = {}
        for column in columns:
            value = record[column]
            if isinstance(value, str):
                row[column] = value
            elif column == "idx":
                row[column] = int(value)
            else:
                row[column] = _float_or_none(value)
        rows.append(row)
    written: int | None = None
    if target is not None:
        try:
            table.to_csv(target, sep=sep, index=False, na_rep="")
        except OSError as exc:
            raise BackendExecutionError(f"could not write {target.name}: {exc}") from exc
        written = target.stat().st_size
    return {
        "status": "completed",
        "operation": "export_tree_node_data",
        "result": {
            "tips_only": tips,
            "columns": columns,
            "row_count": len(rows),
            "rows": rows,
            "output_path": None if target is None else str(target),
            "output_bytes": written,
            "separator": sep,
        },
        "provenance": {
            "package": _package(),
            "convention": (
                "Rows follow toytree's node data (one per node index, or tips only); "
                "a node name selects its "
                "node index; NaN values are null in the rows and empty in the file."
            ),
        },
    }


def prune_tree_tips(
    tree: object,
    names: Sequence[str],
    *,
    preserve_dists: bool = True,
    require_root: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Keep the smallest subtree that connects the named nodes and return it as Newick (toytree
    prune).

    tree: Newick text or a path to a one-tree Newick file. names: at least two node names (tips or
    named internal nodes) to retain; every other tip is removed and the nodes that become unary are
    dropped. preserve_dists: keep the path lengths between retained nodes (toytree default); when
    false the retained edges keep their own lengths. require_root: keep the root of the original
    tree as a node of the result (toytree's option). The retained names and the Newick string of the
    pruned tree are returned.
    """

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    query = _names(names, "names", minimum=2)
    preserve = _flag(preserve_dists, "preserve_dists")
    keep_root = _flag(require_root, "require_root")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    _check_names_in_tree(source, query, "names")
    try:
        pruned = mod.prune(source, *query, preserve_dists=preserve, require_root=keep_root)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not prune the tree: {exc}") from exc
    return _newick_result(
        "prune_tree_tips",
        pruned,
        {"retained_names": query, "preserve_dists": preserve, "require_root": keep_root},
    )


def root_tree_on_outgroup(
    tree: object,
    outgroup: Sequence[str],
    *,
    edge_features: Sequence[str] | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Root a tree on the edge above an outgroup given by node names (toytree root), returned as
    Newick.

    tree: Newick text or a path to a one-tree Newick file. outgroup: one or more node names (tips or
    named internal nodes). With one name the root is placed on that node's edge; with several it is
    placed on the edge above their most recent common ancestor, as toytree's root does.
    edge_features: node features that toytree splits between the two edges of the new root (null for
    none). The outgroup names and the rooted Newick string are returned.
    """

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    query = _names(outgroup, "outgroup")
    features = _feature_list(edge_features)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    _check_names_in_tree(source, query, "outgroup")
    try:
        rooted = mod.root(source, *query, edge_features=features)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not root the tree on the outgroup: {exc}") from exc
    return _newick_result(
        "root_tree_on_outgroup",
        rooted,
        {"outgroup": query, "edge_features": features},
    )


__all__ = [
    "draw_tree_cloud_figure",
    "draw_tree_figure",
    "draw_tree_grid_figure",
    "draw_tree_node_feature_figure",
    "export_tree_node_data",
    "map_node_feature_values",
    "prune_tree_tips",
    "root_tree_on_outgroup",
]

_register_public_verbose(globals())
