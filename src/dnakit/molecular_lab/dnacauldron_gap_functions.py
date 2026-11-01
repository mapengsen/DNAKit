"""DNA Cauldron wrappers for the simulation PDF report and the leveled assembly-plan graph.

dnacauldron and matplotlib are imported only when a function runs.
``dnacauldron_write_simulation_pdf_report`` needs dnacauldron's optional ``reports`` extra
(pdf_reports, which needs WeasyPrint); when it is missing the function raises
BackendUnavailableError. ``dnacauldron_plot_leveled_graph`` draws with matplotlib on an Agg figure
and writes a new PNG. Output files are never overwritten.
"""

from __future__ import annotations

import importlib
import importlib.util
import math
import re
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-dnacauldron"
_PLOT_EXTRA = "bio-plots"
_PLOT_DPI = 150
_MAX_LEVELS = 20
_MAX_NODES = 500
_MAX_ROWS = 50
_NODE = re.compile(r"^[^\x00\n\r]{1,64}$")
_OUTCOME_COLUMN = "Outcome"
_COUNT_COLUMN = "Number of assemblies"


def _text(value: object, name: str, *, maximum: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum or "\x00" in value:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    return value


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not minimum <= float(value) <= maximum:
        raise ConfigurationError(f"{name} must be in [{minimum}, {maximum}].")
    return float(value)


def _new_output(value: object, suffix: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("output_path must be a path string.")
    target = Path(value).expanduser().resolve()
    if target.suffix.lower() != suffix:
        raise ConfigurationError(f"output_path must end with {suffix}.")
    if target.exists():
        raise ConfigurationError(
            f"Output file already exists and will not be overwritten: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _provenance(convention: str) -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    try:
        package = f"dnacauldron {version('dnacauldron')}"
    except PackageNotFoundError:
        package = "dnacauldron"
    return {"package": package, "convention": convention}


def _simulation_rows(value: object) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not 1 <= len(value) <= _MAX_ROWS:
        raise ConfigurationError(
            f"simulation_info must be a list of 1 to {_MAX_ROWS} outcome rows."
        )
    rows: list[dict[str, Any]] = []
    for item in value:
        if not isinstance(item, dict) or set(item) != {_OUTCOME_COLUMN, _COUNT_COLUMN}:
            raise ConfigurationError(
                "each simulation_info row must have exactly the keys "
                f"{_OUTCOME_COLUMN!r} and {_COUNT_COLUMN!r}."
            )
        outcome = _text(item[_OUTCOME_COLUMN], "Outcome", maximum=64)
        count = item[_COUNT_COLUMN]
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ConfigurationError("Number of assemblies must be a non-negative integer.")
        rows.append({_OUTCOME_COLUMN: outcome, _COUNT_COLUMN: count})
    return rows


def _pdf_reports_module() -> Any:
    """Import pdf_reports, which dnacauldron's report writer needs; raise a clear error when
    absent."""

    try:
        return importlib.import_module("pdf_reports")
    except ImportError as exc:
        raise BackendUnavailableError(
            f"pdf_reports could not be imported: {exc}",
            hint=(
                "Install the reports extra of dnacauldron: "
                'python -m pip install "dnacauldron[reports]" (pdf_reports needs WeasyPrint).'
            ),
        ) from exc


def dnacauldron_write_simulation_pdf_report(
    output_path: str,
    simulation_info: list[dict[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the DNA Cauldron simulation summary as a PDF report.

    ``simulation_info`` lists one row per assembly outcome, with the keys ``Outcome`` (for example
    ``Valid``, ``Cancelled`` or ``Errored``) and ``Number of assemblies``. These are the columns of
    the table dnacauldron's AssemblyPlanSimulation builds. Rendering needs dnacauldron's reports
    extra, which in turn needs WeasyPrint.
    """

    target = _new_output(output_path, ".pdf")
    rows = _simulation_rows(simulation_info)
    _pdf_reports_module()
    pandas = require("pandas", _EXTRA)
    reports = require("dnacauldron.reports", _EXTRA)
    frame = pandas.DataFrame(rows, columns=[_OUTCOME_COLUMN, _COUNT_COLUMN])
    reports.write_simulation_pdf_report(str(target), frame)
    return {
        "path": str(target),
        "data": {
            "rows": rows,
            "total_assemblies": sum(row[_COUNT_COLUMN] for row in rows),
            "provenance": _provenance(
                "PDF report written by dnacauldron.reports.write_simulation_pdf_report"
            ),
        },
    }


def _graph_inputs(
    levels: object, edges: object
) -> tuple[list[list[str]], list[list[str]], set[str]]:
    if not isinstance(levels, list) or not 1 <= len(levels) <= _MAX_LEVELS:
        raise ConfigurationError(f"levels must be a list of 1 to {_MAX_LEVELS} levels.")
    seen: set[str] = set()
    cleaned: list[list[str]] = []
    for level in levels:
        if not isinstance(level, list) or not level:
            raise ConfigurationError("each level must be a non-empty list of node names.")
        names = []
        for node in level:
            if not isinstance(node, str) or not _NODE.match(node):
                raise ConfigurationError(
                    "node names must be strings of 1 to 64 printable characters."
                )
            if node in seen:
                raise ConfigurationError(f"node {node!r} appears in more than one place in levels.")
            seen.add(node)
            names.append(node)
        cleaned.append(names)
    if len(seen) > _MAX_NODES:
        raise ConfigurationError(f"a graph can have at most {_MAX_NODES} nodes.")
    if not isinstance(edges, list):
        raise ConfigurationError("edges must be a list of [source, target] pairs.")
    pairs: list[list[str]] = []
    for edge in edges:
        if not isinstance(edge, list) or len(edge) != 2:
            raise ConfigurationError("each edge must be a [source, target] pair.")
        for node in edge:
            if node not in seen:
                raise ConfigurationError(f"edge endpoint {node!r} is not a node of levels.")
        pairs.append([str(edge[0]), str(edge[1])])
    return cleaned, pairs, seen


def dnacauldron_plot_leveled_graph(
    output_path: str,
    levels: list[list[str]],
    edges: list[list[str]],
    *,
    width_factor: float = 2.5,
    height_factor: float = 2.0,
    scale: float = 1.0,
    interlevel_shift: float = 0.0,
    margin: float | None = None,
    font_size: float = 9.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a leveled assembly graph (nodes in columns by level, edges between them) into a PNG.

    ``levels`` lists the nodes of each column and ``edges`` the (source, target) pairs. Each node is
    drawn as a labelled box at the position dnacauldron's layout gives it. Positions are returned
    in the data as fractions of the plot area.
    """

    target = _new_output(output_path, ".png")
    cleaned, pairs, _ = _graph_inputs(levels, edges)
    width = _number(width_factor, "width_factor", minimum=0.1, maximum=50.0)
    height = _number(height_factor, "height_factor", minimum=0.1, maximum=50.0)
    factor = _number(scale, "scale", minimum=0.01, maximum=20.0)
    shift = _number(interlevel_shift, "interlevel_shift", minimum=-1.0, maximum=1.0)
    padding = None if margin is None else _number(margin, "margin", minimum=0.0, maximum=1.0)
    size = _number(font_size, "font_size", minimum=1.0, maximum=72.0)
    graph = require("dnacauldron.AssemblyPlan.plot_leveled_graph", _EXTRA)
    require("matplotlib", _PLOT_EXTRA).use("Agg")
    plt = require("matplotlib.pyplot", _PLOT_EXTRA)

    def draw_node(x: float, y: float, node: str, ax: Any, **_: Any) -> None:
        ax.text(
            x,
            y,
            node,
            ha="center",
            va="center",
            fontsize=size,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "black"},
        )

    positions, axis = graph.plot_leveled_graph(
        levels=cleaned,
        edges=[tuple(pair) for pair in pairs],
        draw_node=draw_node,
        ax=None,
        width_factor=width,
        height_factor=height,
        scale=factor,
        interlevel_shift=shift,
        margin=padding,
    )
    figure = axis.figure
    try:
        figure.savefig(target, format="png", dpi=_PLOT_DPI, bbox_inches="tight")
    finally:
        plt.close(figure)
    return {
        "path": str(target),
        "data": {
            "level_count": len(cleaned),
            "node_count": sum(len(level) for level in cleaned),
            "edge_count": len(pairs),
            "levels": cleaned,
            "edges": pairs,
            "positions": {str(node): [float(x), float(y)] for node, (x, y) in positions.items()},
            "provenance": _provenance(
                "layout and drawing by dnacauldron.AssemblyPlan.plot_leveled_graph"
            ),
        },
    }


__all__ = [
    "dnacauldron_plot_leveled_graph",
    "dnacauldron_write_simulation_pdf_report",
]

_register_public_verbose(globals())
