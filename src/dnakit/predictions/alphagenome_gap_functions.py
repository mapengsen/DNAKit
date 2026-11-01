"""AlphaGenome wrappers: Atlas and DNA model clients, splice junction and track data, and PNG plots.

``create_alphagenome_atlas_client`` and ``create_alphagenome_dna_client`` open a secure gRPC channel
to
the AlphaGenome service with an API key. They are the only network calls in this module. Every other
function works offline: it builds alphagenome data objects (intervals, splice junctions,
transcripts,
track data) from JSON-friendly values, calls alphagenome's own plotting or data function, and either
writes a PNG to an explicit path or returns JSON data.

alphagenome, numpy, pandas and matplotlib are imported only when a function runs. The draw
primitives
that only work on an axis at internal layout coordinates (draw_transcript, draw_strand_arrows and
draw_interval) are not wrapped; plot_alphagenome_transcripts draws through them.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-alphagenome"
_PLOT_EXTRA = "bio-plots"
_NUMPY_EXTRA = "bio-analysis"
_PLOT_MODULE = "alphagenome.visualization.plot"
_COMPONENTS_MODULE = "alphagenome.visualization.plot_components"
_TRANSCRIPTS_MODULE = "alphagenome.visualization.plot_transcripts"
_GENOME_MODULE = "alphagenome.data.genome"
_PLOT_DPI = 150
_DEFAULT_ADDRESS = "dns:///gdmscience.googleapis.com:443"
_ATLAS_SERVICE = "alphagenome.atlas.atlas"
_DNA_SERVICE = "alphagenome.models.dna_client"
_DNA_MODEL_MODULE = "alphagenome.models.v1.dna_model"
_MODEL_VERSIONS = ("all_folds", "fold_0", "fold_1", "fold_2", "fold_3")
_LETTER_COLORS = ("default", "factorbook")
_YLIM_MODES = ("auto", "same")
_SPLICE_STRANDS = ("+", "-")
_STRAND_OPTIONS = ("+", "-", ".")
_COMPONENT_KEYS = {
    "tracks": frozenset({"type", "names", "values", "strands", "track_height", "filled"}),
    "seq_logo": frozenset({"type", "scores", "alphabet", "fig_height", "ylabel"}),
    "empty": frozenset({"type", "fig_height"}),
}
_TRANSCRIPT_KEYS = frozenset({"transcript_id", "strand", "exons", "cds", "info"})
_OPERATION_ERRORS = (ValueError, TypeError, IndexError, KeyError)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _PLOT_EXTRA)


def _secret(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("api_key must be a non-empty string.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _optional_text(value: object, name: str) -> str | None:
    return None if value is None else _text(value, name)


def _choice(value: object, options: tuple[str, ...], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of: {', '.join(options)}.")
    return value


def _integer(value: object, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    return value


def _number(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not math.isfinite(value) or (positive and value <= 0):
        raise ConfigurationError(
            f"{name} must be a finite number{' above zero' if positive else ''}."
        )
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _timeout(value: object) -> float | None:
    return None if value is None else _number(value, "timeout", positive=True)


def _array(np: Any, value: object, name: str, *, ndims: tuple[int, ...]) -> Any:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular list of numbers.") from exc
    if array.ndim not in ndims or array.size == 0:
        shapes = " or ".join(f"{ndim}D" for ndim in ndims)
        raise ConfigurationError(f"{name} must be a non-empty {shapes} list of numbers.")
    if not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must contain only finite numbers.")
    return array


def _matrix(np: Any, value: object, name: str) -> Any:
    array = _array(np, value, name, ndims=(2,))
    if array.shape[1] < 1:
        raise ConfigurationError(f"{name} must have at least one column.")
    return array


def _vector(np: Any, value: object, name: str, *, length: int | None = None) -> Any:
    array = _array(np, value, name, ndims=(1,))
    if length is not None and array.shape[0] != length:
        raise ConfigurationError(f"{name} must have {length} values.")
    return array


def _png_target(output_path: object) -> Path:
    if not isinstance(output_path, str) or not output_path.strip():
        raise ConfigurationError("output_path must be a non-empty string.")
    target = Path(output_path).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _render(output_path: object, operation: str, draw: Callable[[Any], Any]) -> Path:
    """Run ``draw`` (which creates a figure and returns it) and save that figure as a PNG."""

    target = _png_target(output_path)
    plt = _pyplot()
    before = set(plt.get_fignums())
    try:
        try:
            figure = draw(plt)
        except _OPERATION_ERRORS as exc:
            raise ConfigurationError(f"alphagenome could not draw {operation}: {exc}") from exc
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    finally:
        for number in plt.get_fignums():
            if number not in before:
                plt.close(number)
    return target


def _plot_result(target: Path, data: dict[str, Any]) -> dict[str, Any]:
    return {"path": str(target), "data": data}


def _pairs(value: object, name: str) -> list[tuple[int, int]]:
    if not isinstance(value, list) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of [start, end] pairs.")
    pairs: list[tuple[int, int]] = []
    for index, pair in enumerate(value):
        if not isinstance(pair, list) or len(pair) != 2:
            raise ConfigurationError(f"{name}[{index}] must be a [start, end] pair.")
        start = _integer(pair[0], f"{name}[{index}] start")
        end = _integer(pair[1], f"{name}[{index}] end")
        if end <= start:
            raise ConfigurationError(f"{name}[{index}] end must be greater than its start.")
        pairs.append((start, end))
    return pairs


def _span(genome: Any, chromosome: object, start: object, end: object) -> Any:
    name = _text(chromosome, "chromosome")
    low = _integer(start, "start")
    high = _integer(end, "end")
    if high <= low:
        raise ConfigurationError("end must be greater than start.")
    return genome.Interval(name, low, high)


def _span_dict(genome: Any, value: object, name: str) -> Any:
    if not isinstance(value, dict):
        raise ConfigurationError(f"{name} must be an object with chromosome, start and end.")
    return _span(genome, value.get("chromosome"), value.get("start"), value.get("end"))


def _junction_items(value: object, *, with_k: bool) -> list[dict[str, Any]]:
    if not isinstance(value, list):
        raise ConfigurationError("junctions must be a list of objects.")
    items: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        label = f"junctions[{index}]"
        if not isinstance(item, dict):
            raise ConfigurationError(f"{label} must be an object.")
        start = _integer(item.get("start"), f"{label}.start")
        end = _integer(item.get("end"), f"{label}.end")
        if end <= start:
            raise ConfigurationError(f"{label}.end must be greater than start.")
        record: dict[str, Any] = {
            "chromosome": _text(item.get("chromosome"), f"{label}.chromosome"),
            "start": start,
            "end": end,
            "strand": _choice(item.get("strand"), _SPLICE_STRANDS, f"{label}.strand"),
        }
        if with_k:
            count = _number(item.get("k"), f"{label}.k")
            if count < 0:
                raise ConfigurationError(f"{label}.k must not be negative.")
            record["k"] = count
        items.append(record)
    return items


def _junction_object(genome: Any, record: dict[str, Any], k: float | None = None) -> Any:
    if k is None:
        return genome.Junction(
            record["chromosome"], record["start"], record["end"], record["strand"]
        )
    return genome.Junction(
        record["chromosome"], record["start"], record["end"], record["strand"], k=k
    )


def _names(value: object, name: str, count: int) -> list[str]:
    if not isinstance(value, list) or len(value) != count:
        raise ConfigurationError(f"{name} must list {count} names.")
    names = [_text(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if len(set(names)) != len(names):
        raise ConfigurationError(f"{name} must not repeat a name.")
    return names


def _info(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        raise ConfigurationError(f"{name} must be an object with string keys.")
    for key, item in value.items():
        if item is not None and (
            isinstance(item, (list, dict)) or not isinstance(item, (str, int, float, bool))
        ):
            raise ConfigurationError(f"{name}.{key} must be a string, number, boolean or null.")
    return dict(value)


def create_alphagenome_atlas_client(
    api_key: str,
    *,
    timeout: float | None = None,
    address: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Connect to the AlphaGenome Atlas service with an API key and describe the connection.

    This is a network call: alphagenome opens a secure gRPC channel to the Atlas server (default
    ``gdmscience.googleapis.com:443``) and waits until it is ready, for up to ``timeout`` seconds
    (no limit when None). The API key is sent only to that server and is not part of the result.
    The client keeps its channel open, so the client object itself is not returned; query Atlas
    through alphagenome's client directly.
    """

    atlas = require(_ATLAS_SERVICE, _EXTRA)
    key = _secret(api_key)
    wait = _timeout(timeout)
    requested = _optional_text(address, "address")
    server = requested or _DEFAULT_ADDRESS
    client = atlas.create(key, timeout=wait, address=requested)
    return {
        "service": "alphagenome.atlas",
        "source": server,
        "result": {
            "client_type": type(client).__name__,
            "address": server,
            "timeout_seconds": wait,
        },
    }


def create_alphagenome_dna_client(
    api_key: str,
    *,
    model_version: str | None = None,
    timeout: float | None = None,
    address: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Connect to the hosted AlphaGenome DNA model service with an API key and describe the client.

    This is a network call that opens a secure gRPC channel and waits for it to be ready (up to
    ``timeout`` seconds, no limit when None). ``model_version`` is one of ``all_folds`` (the default
    when None), ``fold_0`` to ``fold_3``. The API key is sent only to the server and is not part of
    the result. The client keeps its channel open, so the client object is not returned; use
    alphagenome's client directly for predictions.
    """

    dna_client = require(_DNA_SERVICE, _EXTRA)
    dna_model = require(_DNA_MODEL_MODULE, _EXTRA)
    key = _secret(api_key)
    version = (
        None if model_version is None else _choice(model_version, _MODEL_VERSIONS, "model_version")
    )
    wait = _timeout(timeout)
    requested = _optional_text(address, "address")
    server = requested or _DEFAULT_ADDRESS
    enum_value = None if version is None else dna_model.ModelVersion[version.upper()]
    client = dna_client.create(key, model_version=enum_value, timeout=wait, address=requested)
    return {
        "service": "alphagenome.models.dna_client",
        "source": server,
        "result": {
            "client_type": type(client).__name__,
            "address": server,
            "model_version": version or "all_folds",
            "model_version_source": "requested" if version else "upstream default",
            "timeout_seconds": wait,
        },
    }


def get_alphagenome_junctions_to_plot(
    junctions: list[dict[str, Any]],
    values: list[list[float]],
    track_names: list[str],
    name: str,
    strand: str,
    *,
    k_threshold: float | None = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the splice junctions of one predicted track and strand, using alphagenome's rule.

    ``junctions`` lists ``{"chromosome", "start", "end", "strand"}`` objects, ``values`` holds one
    row
    per junction and one column per track, and ``track_names`` names the columns. Junctions of track
    ``name`` on ``strand`` with k at or above ``k_threshold`` are returned as JSON objects with k
    rounded to two decimals, as alphagenome does. When ``k_threshold`` is None, alphagenome uses 5 %
    of the maximum value. This is data selection only; use ``plot_alphagenome_sashimi`` to draw.
    """

    np = require("numpy", _NUMPY_EXTRA)
    pd = require("pandas", _EXTRA)
    genome = require(_GENOME_MODULE, _EXTRA)
    junction_data = require("alphagenome.data.junction_data", _EXTRA)
    records = _junction_items(junctions, with_k=False)
    names = _names(track_names, "track_names", len(track_names))
    track = _text(name, "name")
    chosen_strand = _choice(strand, _SPLICE_STRANDS, "strand")
    threshold = None if k_threshold is None else _number(k_threshold, "k_threshold")
    try:
        matrix = np.asarray(values, dtype=float).reshape(len(records), len(names))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "values must hold one row per junction and one column per track name."
        ) from exc
    if not np.isfinite(matrix).all():
        raise ConfigurationError("values must contain only finite numbers.")
    junction_array = np.empty(len(records), dtype=object)
    for index, record in enumerate(records):
        junction_array[index] = _junction_object(genome, record)
    data = junction_data.JunctionData(
        junctions=junction_array,
        values=matrix,
        metadata=pd.DataFrame({"name": names}),
    )
    try:
        selected = junction_data.get_junctions_to_plot(
            predictions=data, name=track, strand=chosen_strand, k_threshold=threshold
        )
    except ValueError as exc:
        raise ConfigurationError(f"junction selection failed: {exc}") from exc
    return {
        "track": track,
        "strand": chosen_strand,
        "k_threshold": threshold,
        "count": len(selected),
        "junctions": [
            {
                "chromosome": junction.chromosome,
                "start": int(junction.start),
                "end": int(junction.end),
                "strand": junction.strand,
                "k": float(junction.k),
            }
            for junction in selected
        ],
    }


def plot_alphagenome_seqlogo(
    letter_heights: list[list[float]],
    output_path: str,
    *,
    alphabet: str = "ACGT",
    one_based: bool = True,
    start: int = 0,
    letter_colors: str = "default",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a sequence logo from per-position letter heights to a PNG, using alphagenome's
    ``seqlogo``.

    ``letter_heights`` has one row per position and one column per letter of ``alphabet``. Positive
    heights stack above the axis and negative heights below it. ``letter_colors`` is ``default`` or
    ``factorbook``. The result holds the PNG path and the heights that were drawn.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    heights = _matrix(np, letter_heights, "letter_heights")
    letters = _text(alphabet, "alphabet")
    if heights.shape[1] != len(letters):
        raise ConfigurationError("letter_heights must have one column per alphabet letter.")
    origin = _integer(start, "start")
    one_based_flag = _flag(one_based, "one_based")
    colors = _choice(letter_colors, _LETTER_COLORS, "letter_colors")

    def draw(plt: Any) -> Any:
        figure, ax = plt.subplots(figsize=(max(4.0, 0.25 * len(heights)), 2.5))
        plot.seqlogo(
            heights,
            alphabet=letters,
            one_based=one_based_flag,
            start=origin,
            ax=ax,
            letter_colors=colors,
        )
        return figure

    target = _render(output_path, "seqlogo", draw)
    return _plot_result(
        target,
        {
            "letter_heights": heights.tolist(),
            "alphabet": letters,
            "one_based": one_based_flag,
            "start": origin,
            "letter_colors": colors,
            "positions": int(heights.shape[0]),
        },
    )


def plot_alphagenome_contact_map(
    contact_map: list[list[float]],
    output_path: str,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    square: bool = True,
    cbar_shrink: float = 0.4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a square contact map on a logarithmic colour scale to a PNG, using alphagenome's
    heatmap.

    The colour scale runs from ``vmin`` to ``vmax``; both must be positive with ``vmin`` below
    ``vmax``. When they are None they are taken from the matrix (smallest positive value and largest
    value), as alphagenome does. The result reports the scale that was used.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    matrix = _matrix(np, contact_map, "contact_map")
    if matrix.shape[0] != matrix.shape[1]:
        raise ConfigurationError("contact_map must be square.")
    positives = matrix[matrix > 0]
    low = (
        _number(vmin, "vmin", positive=True)
        if vmin is not None
        else (float(positives.min()) if positives.size else None)
    )
    high = _number(vmax, "vmax") if vmax is not None else float(matrix.max())
    if low is None:
        raise ConfigurationError(
            "contact_map needs at least one positive value, or an explicit vmin."
        )
    if low >= high:
        raise ConfigurationError("the colour scale needs vmin below vmax, both above zero.")
    flag = _flag(square, "square")
    shrink = _number(cbar_shrink, "cbar_shrink", positive=True)

    def draw(plt: Any) -> Any:
        figure, ax = plt.subplots(figsize=(6.0, 5.0))
        plot.plot_contact_map(
            matrix,
            vmin=low,
            vmax=high,
            square=flag,
            cbar_shrink=shrink,
            ax=ax,
        )
        return figure

    target = _render(output_path, "contact map", draw)
    return _plot_result(
        target,
        {
            "contact_map": matrix.tolist(),
            "vmin": low,
            "vmax": high,
            "square": flag,
            "cbar_shrink": shrink,
        },
    )


def plot_alphagenome_track(
    track: list[Any],
    output_path: str,
    *,
    x: list[float] | None = None,
    legend: bool = False,
    ylim: list[float] | None = None,
    color: str | list[str] | None = None,
    filled: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw one track to a PNG, choosing the rendering from its shape as alphagenome's
    ``plot_track`` does.

    A 1D list, or a 2D list with one column, is drawn as a line (filled when ``filled``). A 2D list
    with four columns is drawn as a sequence logo. A 2D list with two columns is drawn as two lines
    (positive and negative strand) with a legend when ``legend``. Other shapes are rejected.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    values = _array(np, track, "track", ndims=(1, 2))
    if values.ndim == 2 and values.shape[1] not in (1, 2, 4):
        raise ConfigurationError("track must have 1, 2 or 4 columns when it is 2D.")
    length = int(values.shape[0])
    positions = _vector(np, x, "x", length=length) if x is not None else None
    if values.ndim == 1 or values.shape[1] == 1:
        rendering = "line"
    elif values.shape[1] == 4:
        rendering = "sequence_logo"
    else:
        rendering = "two_strand_lines"
    limits = _limits(ylim)
    colour = _colour(color)
    legend_flag = _flag(legend, "legend")
    filled_flag = _flag(filled, "filled")

    def draw(plt: Any) -> Any:
        figure, ax = plt.subplots(figsize=(10.0, 3.0))
        plot.plot_track(
            values,
            ax,
            x=positions,
            legend=legend_flag,
            ylim=limits,
            color=colour,
            filled=filled_flag,
        )
        return figure

    target = _render(output_path, "track", draw)
    return _plot_result(
        target,
        {
            "track": values.tolist(),
            "length": length,
            "rendering": rendering,
            "x": None if positions is None else positions.tolist(),
            "ylim": None if limits is None else list(limits),
            "filled": filled_flag,
        },
    )


def _limits(value: object) -> tuple[float, float] | None:
    if value is None:
        return None
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ConfigurationError("ylim must be a [low, high] pair or null.")
    low = _number(value[0], "ylim[0]")
    high = _number(value[1], "ylim[1]")
    if high <= low:
        raise ConfigurationError("ylim high must be above low.")
    return (low, high)


def _colour(value: object) -> str | list[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return _text(value, "color")
    if (
        isinstance(value, list)
        and value
        and all(isinstance(item, str) and item.strip() for item in value)
    ):
        return [item.strip() for item in value]
    raise ConfigurationError("color must be a colour string, a list of colour strings, or null.")


def plot_alphagenome_tracks(
    tracks: dict[str, list[Any]],
    output_path: str,
    *,
    x: list[float] | None = None,
    title: str | None = None,
    legend: bool = False,
    fig_width: float = 20.0,
    fig_track_height: float | dict[str, float] = 1.5,
    ylim: str | list[float] | None = "auto",
    yticks_min_max_only: bool = False,
    ylab: bool = True,
    color: str | dict[str, str] | None = None,
    horizontal_ylab: bool = True,
    filled_tracks: list[str] | None = None,
    despine: bool = True,
    despine_keep_bottom: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw several tracks as stacked panels in one PNG, using alphagenome's ``plot_tracks``.

    ``tracks`` maps each track name to a list of values (1D, or 2D with one column per channel); all
    tracks must have the same length. ``ylim`` is ``auto`` (one scale per panel), ``same`` (one
    shared
    scale), a ``[low, high]`` pair or null. ``fig_track_height`` is the total height in inches, or a
    per-track dict of heights in the track order. ``color`` is one colour or a per-track dict.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    if not isinstance(tracks, dict) or not tracks:
        raise ConfigurationError("tracks must be a non-empty object mapping names to value lists.")
    series: dict[str, Any] = {}
    for name, values in tracks.items():
        label = _text(name, "track name")
        series[label] = _array(np, values, f"tracks[{label}]", ndims=(1, 2))
    if len({len(values) for values in series.values()}) != 1:
        raise ConfigurationError("all tracks must have the same length.")
    length = len(next(iter(series.values())))
    positions = _vector(np, x, "x", length=length) if x is not None else None
    if ylim is None or isinstance(ylim, (list, tuple)):
        limits: str | tuple[float, float] | None = _limits(ylim)
    elif isinstance(ylim, str):
        limits = _choice(ylim, _YLIM_MODES, "ylim")
    else:
        raise ConfigurationError("ylim must be 'auto', 'same', a [low, high] pair or null.")
    heights = _heights(fig_track_height, list(series))
    colours = _track_colours(color, list(series))
    filled = _names_subset(filled_tracks, "filled_tracks")
    width = _number(fig_width, "fig_width", positive=True)
    flags = {
        "legend": _flag(legend, "legend"),
        "ylab": _flag(ylab, "ylab"),
        "horizontal_ylab": _flag(horizontal_ylab, "horizontal_ylab"),
        "yticks_min_max_only": _flag(yticks_min_max_only, "yticks_min_max_only"),
        "despine": _flag(despine, "despine"),
        "despine_keep_bottom": _flag(despine_keep_bottom, "despine_keep_bottom"),
    }
    chart_title = _optional_text(title, "title")

    def draw(plt: Any) -> Any:
        return plot.plot_tracks(
            series,
            x=positions,
            title=chart_title,
            legend=flags["legend"],
            fig_width=width,
            fig_track_height=heights,
            ylim=limits,
            yticks_min_max_only=flags["yticks_min_max_only"],
            ylab=flags["ylab"],
            color=colours,
            horizontal_ylab=flags["horizontal_ylab"],
            filled_tracks=filled,
            despine=flags["despine"],
            despine_keep_bottom=flags["despine_keep_bottom"],
        )

    target = _render(output_path, "tracks", draw)
    return _plot_result(
        target,
        {
            "tracks": {name: values.tolist() for name, values in series.items()},
            "track_names": list(series),
            "length": length,
            "x": None if positions is None else positions.tolist(),
            "ylim": ylim,
            "title": chart_title,
            "filled_tracks": filled or [],
        },
    )


def _heights(value: object, names: list[str]) -> float | dict[str, float]:
    if isinstance(value, dict):
        if set(value) != set(names):
            raise ConfigurationError("fig_track_height must give one height per track name.")
        return {
            name: _number(value[name], f"fig_track_height[{name}]", positive=True) for name in names
        }
    return _number(value, "fig_track_height", positive=True)


def _track_colours(value: object, names: list[str]) -> str | dict[str, str] | None:
    if value is None:
        return None
    if isinstance(value, dict):
        if set(value) - set(names):
            raise ConfigurationError("color names a track that is not in tracks.")
        return {name: _text(value[name], f"color[{name}]") for name in value if name in names}
    return _text(value, "color")


def _names_subset(value: object, name: str) -> list[str] | None:
    if value is None:
        return None
    if not isinstance(value, list):
        raise ConfigurationError(f"{name} must be a list of track names.")
    return [_text(item, f"{name}[{index}]") for index, item in enumerate(value)]


def plot_alphagenome_sashimi(
    junctions: list[dict[str, Any]],
    output_path: str,
    *,
    interval: dict[str, Any] | None = None,
    filter_threshold: float = 0.01,
    annotate_counts: bool = True,
    seed: int | None = 0,
    color: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw splice junctions as a Sashimi plot to a PNG, using alphagenome's ``sashimi_plot``.

    ``junctions`` lists ``{"chromosome", "start", "end", "strand", "k"}`` objects, where k is the
    read
    count. Junctions with k below ``filter_threshold`` are not drawn. ``interval`` (an object with
    chromosome, start and end) limits the count labels to the region. ``seed`` fixes the random arc
    jitter so that the PNG is reproducible; None uses an unseeded generator. At least one k must be
    above zero, because alphagenome scales arc widths by the total read count.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    genome = require(_GENOME_MODULE, _EXTRA)
    records = _junction_items(junctions, with_k=True)
    threshold = _number(filter_threshold, "filter_threshold")
    if threshold < 0:
        raise ConfigurationError("filter_threshold must not be negative.")
    total = float(sum(record["k"] for record in records))
    if records and total <= 0:
        raise ConfigurationError("at least one junction needs k above zero.")
    window = None if interval is None else _span_dict(genome, interval, "interval")
    seed_value = None if seed is None else _integer(seed, "seed")
    annotate = _flag(annotate_counts, "annotate_counts")
    colour = _optional_text(color, "color")
    drawn = [record for record in records if record["k"] >= threshold]

    def draw(plt: Any) -> Any:
        figure, ax = plt.subplots(figsize=(8.0, 3.0))
        plot.sashimi_plot(
            [_junction_object(genome, record, k=record["k"]) for record in records],
            ax,
            interval=window,
            filter_threshold=threshold,
            annotate_counts=annotate,
            rng=np.random.default_rng(seed_value),
            color=colour,
        )
        return figure

    target = _render(output_path, "Sashimi plot", draw)
    return _plot_result(
        target,
        {
            "junctions": records,
            "drawn": drawn,
            "drawn_count": len(drawn),
            "total_k": total,
            "filter_threshold": threshold,
            "annotate_counts": annotate,
            "seed": seed_value,
            "ylim": [0.0, 0.15],
        },
    )


def pad_alphagenome_track(
    track: list[list[float]],
    new_len: int,
    *,
    value: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Centre a (positions x channels) track inside a longer window, filled with ``value``.

    This pads a model output so that it can be aligned with a wider genomic window. Alphagenome's
    ``pad_track`` puts the same padding on both sides when the difference is even; with an odd
    difference the extra position goes at the end. ``offset`` is the index of the first original
    row.
    """

    np = require("numpy", _NUMPY_EXTRA)
    plot = require(_PLOT_MODULE, _EXTRA)
    array = _matrix(np, track, "track")
    target_length = _integer(new_len, "new_len", minimum=1)
    if target_length < array.shape[0]:
        raise ConfigurationError("new_len must be at least the track length.")
    fill = _number(value, "value")
    padded = plot.pad_track(array, target_length, fill)
    return {
        "length": int(array.shape[0]),
        "new_length": target_length,
        "offset": (target_length - int(array.shape[0])) // 2,
        "value": fill,
        "padded": padded.tolist(),
    }


def plot_alphagenome_components(
    components: list[dict[str, Any]],
    output_path: str,
    *,
    chromosome: str,
    start: int,
    end: int,
    fig_width: float = 20.0,
    fig_height_scale: float = 1.0,
    title: str | None = None,
    despine: bool = True,
    despine_keep_bottom: bool = False,
    hspace: float = 0.3,
    xlabel: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw model outputs as stacked components over one interval to a PNG, using alphagenome's
    ``plot``.

    Each component is an object with a ``type``:

    - ``tracks``: ``names``, ``values`` (one row per base, one column per track), optional
      ``strands`` (default ``+``), ``track_height`` and ``filled``;
    - ``seq_logo``: ``scores`` (one row per base, one column per letter), optional ``alphabet``
      (default ``ACGT``), ``fig_height`` and ``ylabel``;
    - ``empty``: optional ``fig_height``, a blank spacer.

    The number of rows must equal ``end - start``. Transcript and contact-map components need
    alphagenome annotation objects and are drawn by ``plot_alphagenome_transcripts`` and
    ``plot_alphagenome_contact_map``; annotations are not supported here.
    """

    np = require("numpy", _NUMPY_EXTRA)
    pd = require("pandas", _EXTRA)
    genome = require(_GENOME_MODULE, _EXTRA)
    track_data = require("alphagenome.data.track_data", _EXTRA)
    components_module = require(_COMPONENTS_MODULE, _EXTRA)
    if not isinstance(components, list) or not components:
        raise ConfigurationError("components must be a non-empty list of objects.")
    span = _span(genome, chromosome, start, end)
    width = span.end - span.start
    objects: list[Any] = []
    described: list[dict[str, Any]] = []
    for index, item in enumerate(components):
        component, record = _component(
            np, pd, components_module, track_data, span, width, item, index
        )
        objects.append(component)
        described.append(record)
    chart_title = _optional_text(title, "title")
    label = _optional_text(xlabel, "xlabel")
    width_inches = _number(fig_width, "fig_width", positive=True)
    height_scale = _number(fig_height_scale, "fig_height_scale", positive=True)
    spacing = _number(hspace, "hspace")
    flags = (_flag(despine, "despine"), _flag(despine_keep_bottom, "despine_keep_bottom"))

    def draw(plt: Any) -> Any:
        return components_module.plot(
            objects,
            span,
            fig_width=width_inches,
            fig_height_scale=height_scale,
            title=chart_title,
            despine=flags[0],
            despine_keep_bottom=flags[1],
            hspace=spacing,
            xlabel=label,
        )

    target = _render(output_path, "components", draw)
    return _plot_result(
        target,
        {
            "interval": {"chromosome": span.chromosome, "start": span.start, "end": span.end},
            "components": described,
            "title": chart_title,
        },
    )


def _component(
    np: Any,
    pd: Any,
    components_module: Any,
    track_data: Any,
    span: Any,
    width: int,
    item: object,
    index: int,
) -> tuple[Any, dict[str, Any]]:
    label = f"components[{index}]"
    if not isinstance(item, dict):
        raise ConfigurationError(f"{label} must be an object.")
    kind = _choice(item.get("type"), tuple(_COMPONENT_KEYS), f"{label}.type")
    extra = set(item) - _COMPONENT_KEYS[kind]
    if extra:
        raise ConfigurationError(f"{label} has unsupported keys: {', '.join(sorted(extra))}.")
    if kind == "empty":
        height = _number(item.get("fig_height", 1.0), f"{label}.fig_height", positive=True)
        return components_module.EmptyComponent(fig_height=height), {
            "type": "empty",
            "fig_height": height,
        }
    if kind == "tracks":
        values = _matrix(np, item.get("values"), f"{label}.values")
        if values.shape[0] != width:
            raise ConfigurationError(f"{label}.values must have one row per base ({width}).")
        names = _names(item.get("names"), f"{label}.names", values.shape[1])
        strands = item.get("strands", ["+"] * len(names))
        if not isinstance(strands, list) or len(strands) != len(names):
            raise ConfigurationError(f"{label}.strands must give one strand per track.")
        strands = [_choice(strand, _STRAND_OPTIONS, f"{label}.strands") for strand in strands]
        if len(set(zip(names, strands, strict=True))) != len(names):
            raise ConfigurationError(f"{label} repeats a (name, strand) pair.")
        height = _number(item.get("track_height", 1.0), f"{label}.track_height", positive=True)
        filled = _flag(item.get("filled", False), f"{label}.filled")
        tdata = track_data.TrackData(
            values=values.astype(np.float32),
            metadata=pd.DataFrame({"name": names, "strand": strands}),
            resolution=1,
            interval=span,
        )
        component = components_module.Tracks(tdata, track_height=height, filled=filled)
        return component, {
            "type": "tracks",
            "names": names,
            "strands": strands,
            "values": values.tolist(),
            "track_height": height,
            "filled": filled,
        }
    scores = _matrix(np, item.get("scores"), f"{label}.scores")
    if scores.shape[0] != width:
        raise ConfigurationError(f"{label}.scores must have one row per base ({width}).")
    alphabet = _text(item.get("alphabet", "ACGT"), f"{label}.alphabet")
    if scores.shape[1] != len(alphabet):
        raise ConfigurationError(f"{label}.scores must have one column per alphabet letter.")
    height = _number(item.get("fig_height", 1.0), f"{label}.fig_height", positive=True)
    ylabel = item.get("ylabel", "")
    if not isinstance(ylabel, str):
        raise ConfigurationError(f"{label}.ylabel must be a string.")
    component = components_module.SeqLogo(
        scores=scores.astype(np.float32),
        scores_interval=span,
        fig_height=height,
        alphabet=alphabet,
        ylabel=ylabel,
    )
    return component, {
        "type": "seq_logo",
        "scores": scores.tolist(),
        "alphabet": alphabet,
        "fig_height": height,
        "ylabel": ylabel,
    }


def plot_alphagenome_transcripts(
    transcripts: list[dict[str, Any]],
    output_path: str,
    *,
    chromosome: str,
    start: int,
    end: int,
    zero_origin: bool = False,
    label_name: str | None = None,
    plot_labels_once: bool = False,
    fig_width: float = 12.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw transcript models over one interval to a PNG, using alphagenome's ``plot_transcripts``.

    Each transcript is an object with ``transcript_id``, ``strand`` (``+`` or ``-``), ``exons`` (a
    list
    of ``[start, end]`` pairs), optional ``cds`` (coding pairs; omit or null for non-coding) and an
    optional ``info`` object. ``label_name`` names the ``info`` key written next to each transcript,
    and every transcript must have it. ``zero_origin`` shifts the x axis to start at 0.
    """

    genome = require(_GENOME_MODULE, _EXTRA)
    transcript_module = require("alphagenome.data.transcript", _EXTRA)
    plot_transcripts = require(_TRANSCRIPTS_MODULE, _EXTRA)
    if not isinstance(transcripts, list) or not transcripts:
        raise ConfigurationError("transcripts must be a non-empty list of objects.")
    span = _span(genome, chromosome, start, end)
    label = _optional_text(label_name, "label_name")
    objects: list[Any] = []
    described: list[dict[str, Any]] = []
    for index, item in enumerate(transcripts):
        transcript, record = _transcript(genome, transcript_module, span.chromosome, item, index)
        if label is not None and label not in record["info"]:
            raise ConfigurationError(f"transcripts[{index}] has no info value named {label!r}.")
        objects.append(transcript)
        described.append(record)
    shift = _flag(zero_origin, "zero_origin")
    once = _flag(plot_labels_once, "plot_labels_once")
    width = _number(fig_width, "fig_width", positive=True)
    height = 1.0 + 0.6 * len(objects)

    def draw(plt: Any) -> Any:
        figure, ax = plt.subplots(figsize=(width, height))
        plot_transcripts.plot_transcripts(
            ax,
            objects,
            span,
            zero_origin=shift,
            label_name=label,
            plot_labels_once=once,
        )
        return figure

    target = _render(output_path, "transcripts", draw)
    return _plot_result(
        target,
        {
            "interval": {"chromosome": span.chromosome, "start": span.start, "end": span.end},
            "transcripts": described,
            "zero_origin": shift,
            "label_name": label,
            "plot_labels_once": once,
        },
    )


def _transcript(
    genome: Any, transcript_module: Any, chromosome: str, item: object, index: int
) -> tuple[Any, dict[str, Any]]:
    label = f"transcripts[{index}]"
    if not isinstance(item, dict):
        raise ConfigurationError(f"{label} must be an object.")
    extra = set(item) - _TRANSCRIPT_KEYS
    if extra:
        raise ConfigurationError(f"{label} has unsupported keys: {', '.join(sorted(extra))}.")
    transcript_id = _text(item.get("transcript_id"), f"{label}.transcript_id")
    strand = _choice(item.get("strand"), _SPLICE_STRANDS, f"{label}.strand")
    exons = _pairs(item.get("exons"), f"{label}.exons")
    cds_value = item.get("cds")
    cds = None if cds_value is None else _pairs(cds_value, f"{label}.cds")
    info = _info(item.get("info", {}), f"{label}.info")
    transcript = transcript_module.Transcript(
        exons=[genome.Interval(chromosome, low, high, strand) for low, high in exons],
        cds=None
        if cds is None
        else [genome.Interval(chromosome, low, high, strand) for low, high in cds],
        transcript_id=transcript_id,
        info=dict(info),
    )
    return transcript, {
        "transcript_id": transcript_id,
        "strand": strand,
        "exons": [list(pair) for pair in exons],
        "cds": None if cds is None else [list(pair) for pair in cds],
        "info": info,
    }


__all__ = [
    "create_alphagenome_atlas_client",
    "create_alphagenome_dna_client",
    "get_alphagenome_junctions_to_plot",
    "pad_alphagenome_track",
    "plot_alphagenome_components",
    "plot_alphagenome_contact_map",
    "plot_alphagenome_sashimi",
    "plot_alphagenome_seqlogo",
    "plot_alphagenome_track",
    "plot_alphagenome_tracks",
    "plot_alphagenome_transcripts",
]

_register_public_verbose(globals(), tuple(__all__))
