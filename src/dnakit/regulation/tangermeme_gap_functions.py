"""tangermeme wrappers for sequence logos, PWM logos, attribution plots and categorical scatters.

Each plot is drawn with tangermeme's own function on the matplotlib Agg backend and written to a new
PNG file. The function returns that path together with the numbers that were drawn, never a figure.
Attribution rows are positions, with columns A, C, G, T; tangermeme's (4, length) layout is used
internally. Models are never loaded here: ``plot_sequence_attributions`` takes a torch model that is
already in eval mode, as the other tangermeme wrappers in this package do.

The interactive logo needs mpld3 for its hover tooltips, which the bio-tangermeme extra does not
install. Its wrapper therefore reports a missing mpld3 as BackendUnavailableError and writes a
static
PNG, so hover annotations are not captured.

torch, numpy, pandas, matplotlib and tangermeme are imported only when a function runs. No network
is used.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.regulation.backends import (
    _dna,
    _flag,
    _integer,
    _one_hot_batch,
    _tensor_provenance,
    _torch_model,
)

_TANGERMEME_EXTRA = "bio-tangermeme"
_BASES = ("A", "C", "G", "T")
_PLOT_DPI = 150
_MAX_POSITIONS = 10_000
_MAX_SEQUENCES = 32
_MAX_LENGTH = 10_000
_MAX_ANNOTATIONS = 1_000
_MAX_SHUFFLES = 200
_MAX_PATH = 4_096
_MPLD3_HINT = (
    "Install mpld3 with: python -m pip install mpld3 (it is not part of dnakit[bio-tangermeme])"
)
_UPSTREAM_ERRORS = (ValueError, TypeError, IndexError, KeyError, RuntimeError)


def _numpy() -> Any:
    return require("numpy", _TANGERMEME_EXTRA)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _TANGERMEME_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _TANGERMEME_EXTRA)


def _new_png(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path string.")
    if len(str(value)) > _MAX_PATH or "\x00" in str(value) or not str(value).strip():
        raise ConfigurationError(f"{name} is not a usable path.")
    target = Path(str(value).strip()).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError(f"{name} must end with .png.")
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _position_rows(values: object, name: str, *, limit: int = _MAX_POSITIONS) -> list[list[float]]:
    """A non-empty list of positions, each with four finite numbers for A, C, G and T."""

    if isinstance(values, (str, bytes, Mapping)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of rows.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} may contain at most {limit} positions.")
    rows: list[list[float]] = []
    for index, row in enumerate(values):
        label = f"{name}[{index}]"
        if isinstance(row, (str, bytes, Mapping)) or not isinstance(row, Sequence) or len(row) != 4:
            raise ConfigurationError(f"{label} must hold four numbers for A, C, G and T.")
        numbers: list[float] = []
        for item in row:
            if (
                isinstance(item, bool)
                or not isinstance(item, (int, float))
                or not math.isfinite(item)
            ):
                raise ConfigurationError(f"{label} must hold finite numbers.")
            numbers.append(float(item))
        rows.append(numbers)
    return rows


def _window(start: object, end: object, length: int) -> tuple[int | None, int | None]:
    first = None if start is None else _integer(start, "start", 0, length)
    last = None if end is None else _integer(end, "end", 1, length)
    if first is not None and last is not None and first >= last:
        raise ConfigurationError("start must be smaller than end.")
    return first, last


def _symmetric_limit(value: object) -> float | None:
    """A bound B for the logo's y axis, which runs from -B to B."""

    if value is None:
        return None
    bound = _finite(value, "ylim")
    if bound <= 0:
        raise ConfigurationError("ylim must be a positive bound; the axis runs from -ylim to ylim.")
    return bound


def _ylim_pair(value: object) -> tuple[float, float] | None:
    """A (low, high) pair for the interactive logo, which applies it directly."""

    if value is None:
        return None
    if (
        isinstance(value, (str, bytes, Mapping))
        or not isinstance(value, Sequence)
        or len(value) != 2
    ):
        raise ConfigurationError("ylim must be None or a pair of numbers.")
    low, high = (_finite(item, "ylim") for item in value)
    if low >= high:
        raise ConfigurationError("ylim must be (low, high) with low < high.")
    return low, high


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _annotation_frame(annotations: object, length: int) -> Any:
    """Records with name, start, end and score, in the column order plot_logo reads (name first)."""

    pandas = require("pandas", _TANGERMEME_EXTRA)
    if annotations is None:
        return None
    if isinstance(annotations, (str, bytes, Mapping)) or not isinstance(annotations, Sequence):
        raise ConfigurationError("annotations must be None or a list of records.")
    if not annotations or len(annotations) > _MAX_ANNOTATIONS:
        raise ConfigurationError(f"annotations must hold 1 to {_MAX_ANNOTATIONS} records.")
    records: list[list[Any]] = []
    for index, item in enumerate(annotations):
        label = f"annotations[{index}]"
        if not isinstance(item, Mapping) or set(item) - {"name", "start", "end", "score"}:
            raise ConfigurationError(f"{label} must be a mapping with name, start, end and score.")
        if not all(key in item for key in ("name", "start", "end", "score")):
            raise ConfigurationError(f"{label} must have name, start, end and score.")
        name = item["name"]
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError(f"{label}.name must be a non-empty string.")
        start = _integer(item["start"], f"{label}.start", 0, length)
        end = _integer(item["end"], f"{label}.end", 1, length)
        if start >= end:
            raise ConfigurationError(f"{label} must have start smaller than end.")
        records.append([name.strip(), start, end, _finite(item["score"], f"{label}.score")])
    return pandas.DataFrame(records, columns=["motif_name", "start", "end", "score"])


def _logo_size(length: int, start: int | None, end: int | None) -> tuple[float, float]:
    span = (end or length) - (start or 0)
    return (min(max(4.0, 0.25 * span), 60.0), 2.5)


def _rows_to_tensor(torch: Any, rows: list[list[float]]) -> Any:
    return torch.tensor(rows, dtype=torch.float32).transpose(0, 1).contiguous()


def plot_attribution_logo(
    attributions: Sequence[Sequence[float]],
    output_path: str,
    *,
    start: int | None = None,
    end: int | None = None,
    ylim: float | None = None,
    annotations: Sequence[Mapping[str, Any]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a sequence logo of an (L, 4) attribution matrix (A, C, G, T) to a PNG with tangermeme.

    Characters are scaled by their attribution. ``start`` and ``end`` choose the window, and
    ``ylim`` is a symmetric bound: the axis runs from -ylim to ylim. ``annotations`` is a list of
    records with ``name``, ``start``, ``end`` and ``score``; they are drawn under the logo.
    tangermeme keeps only annotations that start after the window start and end before its end, so
    an annotation starting at position 0 or ending at the last position is not drawn.
    """

    rows = _position_rows(attributions, "attributions")
    first, last = _window(start, end, len(rows))
    bound = _symmetric_limit(ylim)
    frame = _annotation_frame(annotations, len(rows))
    target = _new_png(output_path, "output_path")
    torch = require("torch", _TANGERMEME_EXTRA)
    plot = require("tangermeme.plot", _TANGERMEME_EXTRA)
    plt = _pyplot()
    figure = plt.figure(figsize=_logo_size(len(rows), first, last))
    try:
        ax = figure.add_subplot(111)
        # tangermeme takes a scalar bound when annotations are drawn (a pair fails there) and
        # ignores ylim without annotations, so the limit is applied to the axis after drawing.
        plot.plot_logo(
            _rows_to_tensor(torch, rows),
            ax=ax,
            start=first,
            end=last,
            ylim=bound,
            annotations=frame,
        )
        if bound is not None:
            ax.set_ylim(-bound, bound)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"tangermeme could not draw the logo: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_attribution_logo",
        "path": str(target),
        "image_format": "png",
        "data": {
            "position_count": len(rows),
            "start": first,
            "end": last,
            "ylim": None if bound is None else [-bound, bound],
            "alphabet": list(_BASES),
            "attributions": rows,
            "annotations": None if annotations is None else [dict(item) for item in annotations],
        },
        "provenance": _tensor_provenance(torch, plot="tangermeme.plot.plot_logo"),
    }


def plot_pwm_logo(
    pwm: Sequence[Sequence[float]],
    output_path: str,
    *,
    name: str | None = None,
    eps: float = 1e-7,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a position weight matrix as an information-content logo PNG (tangermeme's plot_pwm).

    ``pwm`` has one row per position with four probabilities for A, C, G, T; every row must be
    non-negative and sum to 1 (to 1e-4). ``name`` is an optional title and ``eps`` guards the
    logarithm.
    """

    rows = _position_rows(pwm, "pwm")
    for index, row in enumerate(rows):
        if any(value < 0 for value in row) or abs(sum(row) - 1.0) > 1e-4:
            raise ConfigurationError(f"pwm[{index}] must be non-negative and sum to 1.")
    if name is not None and (not isinstance(name, str) or len(name) > 200):
        raise ConfigurationError("name must be None or a string of at most 200 characters.")
    if isinstance(eps, bool) or not isinstance(eps, (int, float)) or not 0.0 < eps < 1.0:
        raise ConfigurationError("eps must be a number in (0, 1).")
    target = _new_png(output_path, "output_path")
    numpy = _numpy()
    plot = require("tangermeme.plot", _TANGERMEME_EXTRA)
    plt = _pyplot()
    figure = plt.figure(figsize=(min(max(4.0, 0.3 * len(rows)), 60.0), 2.5))
    try:
        ax = figure.add_subplot(111)
        plot.plot_pwm(numpy.asarray(rows, dtype=float).T, ax=ax, name=name, eps=float(eps))
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"tangermeme could not draw the PWM: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_pwm_logo",
        "path": str(target),
        "image_format": "png",
        "data": {
            "position_count": len(rows),
            "alphabet": list(_BASES),
            "name": name,
            "pwm": rows,
        },
        "provenance": _tensor_provenance(
            require("torch", _TANGERMEME_EXTRA), plot="tangermeme.plot.plot_pwm"
        ),
    }


def plot_sequence_attributions(
    model: Any,
    sequences: Sequence[str],
    output_path: str,
    *,
    target: int = 0,
    n_shuffles: int = 20,
    random_state: int = 0,
    hypothetical: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute DeepLIFT-SHAP attributions for each sequence and draw one logo per sequence.

    ``model`` is a torch module already in eval mode; ``sequences`` are DNA strings of one length
    (up to 32). Attributions are tangermeme's ``deep_lift_shap`` against dinucleotide-shuffled
    references seeded with ``random_state``. ``target`` picks the model output. The PNG holds one
    logo per sequence, stacked. The result gives the attributions as (sequence, position, base)
    lists.
    """

    torch = _torch_model(model)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    shuffle = require("tangermeme.deep_lift_shap", _TANGERMEME_EXTRA)
    plot = require("tangermeme.plot", _TANGERMEME_EXTRA)
    if isinstance(sequences, (str, bytes, Mapping)) or not isinstance(sequences, Sequence):
        raise ConfigurationError("sequences must be a list of DNA strings.")
    if not 1 <= len(sequences) <= _MAX_SEQUENCES:
        raise ConfigurationError(f"sequences must hold 1 to {_MAX_SEQUENCES} DNA strings.")
    texts = [_dna(item, f"sequences[{index}]", _MAX_LENGTH) for index, item in enumerate(sequences)]
    length = len(texts[0])
    if any(len(text) != length for text in texts):
        raise ConfigurationError("sequences must all have the same length.")
    if isinstance(target, bool) or not isinstance(target, int) or target < 0:
        raise ConfigurationError("target must be a non-negative integer.")
    shuffles = _integer(n_shuffles, "n_shuffles", 1, _MAX_SHUFFLES)
    seed = _integer(random_state, "random_state", 0, 2**31 - 1)
    flag = _flag(hypothetical, "hypothetical")
    target_file = _new_png(output_path, "output_path")

    batch = _one_hot_batch(torch, utils, [(str(index), text) for index, text in enumerate(texts)])
    attribute_kwargs = {
        "target": target,
        "n_shuffles": shuffles,
        "hypothetical": flag,
        "references": partial(shuffle.dinucleotide_shuffle, random_state=seed),
    }
    plt = _pyplot()
    figure = plt.figure(figsize=(min(max(4.0, 0.25 * length), 60.0), 1.8 * len(texts)))
    try:
        try:
            _axes, attributions = plot.plot_attributions(
                model,
                batch,
                attribute_kwargs=attribute_kwargs,
                layout=None,
            )
        except _UPSTREAM_ERRORS as exc:
            raise BackendExecutionError(f"tangermeme attribution failed: {exc}") from exc
        figure.savefig(target_file, format="png", dpi=_PLOT_DPI)
    finally:
        plt.close(figure)
    values = attributions.permute(0, 2, 1).tolist()
    return {
        "status": "completed",
        "operation": "plot_sequence_attributions",
        "path": str(target_file),
        "image_format": "png",
        "data": {
            "sequence_count": len(texts),
            "length": length,
            "target": target,
            "n_shuffles": shuffles,
            "random_state": seed,
            "hypothetical": flag,
            "sequences": texts,
            "attributions": values,
        },
        "provenance": _tensor_provenance(
            torch, attribution="tangermeme.deep_lift_shap", plot="tangermeme.plot.plot_attributions"
        ),
    }


def plot_categorical_attribution_scatter(
    attributions: Sequence[Sequence[float]],
    output_path: str,
    *,
    size: float = 20.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scatter the summed attribution at every position, coloured by the dominant base (tangermeme).

    Each position is coloured by the base with the largest absolute value there (A green, C blue,
    G orange, T red), at height equal to the sum over the four bases. ``size`` is the marker size.
    Use it for long sequences that a logo cannot show legibly. Colours are fixed: tangermeme's
    ``colors`` argument fails with an UnboundLocalError, so it is not exposed.
    """

    rows = _position_rows(attributions, "attributions")
    marker = _finite(size, "size")
    if marker <= 0:
        raise ConfigurationError("size must be positive.")
    torch = require("torch", _TANGERMEME_EXTRA)
    plot = require("tangermeme.plot", _TANGERMEME_EXTRA)
    target = _new_png(output_path, "output_path")
    tensor = _rows_to_tensor(torch, rows)
    dominant = torch.abs(tensor).argmax(dim=0).tolist()
    totals = tensor.sum(dim=0).tolist()
    plt = _pyplot()
    figure = plt.figure(figsize=(min(max(4.0, 0.25 * len(rows)), 60.0), 2.5))
    try:
        ax = figure.add_subplot(111)
        plot.plot_categorical_scatter(tensor, ax=ax, s=marker)
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"tangermeme could not draw the scatter: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "plot_categorical_attribution_scatter",
        "path": str(target),
        "image_format": "png",
        "data": {
            "position_count": len(rows),
            "alphabet": list(_BASES),
            "dominant_base_index_by_position": [int(index) for index in dominant],
            "summed_attribution_by_position": [float(value) for value in totals],
            "attributions": rows,
        },
        "provenance": _tensor_provenance(torch, plot="tangermeme.plot.plot_categorical_scatter"),
    }


def render_interactive_logo_png(
    attributions: Sequence[Sequence[float]],
    output_path: str,
    *,
    start: int | None = None,
    end: int | None = None,
    ylim: Sequence[float] | None = None,
    annotations: Sequence[Mapping[str, Any]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw tangermeme's interactive logo to a static PNG (hover tooltips are not captured).

    The same inputs as ``plot_attribution_logo``. tangermeme's ``interactive_logo`` needs mpld3,
    which is not part of the bio-tangermeme extra; without it this function raises
    BackendUnavailableError and writes nothing.
    """

    rows = _position_rows(attributions, "attributions")
    first, last = _window(start, end, len(rows))
    limits = _ylim_pair(ylim)
    frame = _annotation_frame(annotations, len(rows))
    try:
        require("mpld3", _TANGERMEME_EXTRA)
    except BackendUnavailableError as exc:
        raise BackendUnavailableError(exc.message, hint=_MPLD3_HINT) from exc
    target = _new_png(output_path, "output_path")
    torch = require("torch", _TANGERMEME_EXTRA)
    plot = require("tangermeme.plot", _TANGERMEME_EXTRA)
    plt = _pyplot()
    figure = plt.figure(figsize=_logo_size(len(rows), first, last))
    try:
        ax = figure.add_subplot(111)
        plot.interactive_logo(
            _rows_to_tensor(torch, rows),
            ax=ax,
            start=first,
            end=last,
            ylim=limits,
            annotations=frame,
        )
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"tangermeme could not draw the interactive logo: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "render_interactive_logo_png",
        "path": str(target),
        "image_format": "png",
        "data": {
            "position_count": len(rows),
            "start": first,
            "end": last,
            "ylim": None if limits is None else list(limits),
            "alphabet": list(_BASES),
            "attributions": rows,
        },
        "provenance": _tensor_provenance(torch, plot="tangermeme.plot.interactive_logo"),
    }


__all__ = [
    "plot_attribution_logo",
    "plot_categorical_attribution_scatter",
    "plot_pwm_logo",
    "plot_sequence_attributions",
    "render_interactive_logo_png",
]

_register_public_verbose(globals())
