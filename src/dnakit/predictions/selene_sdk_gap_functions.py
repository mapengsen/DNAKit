"""selene-sdk wrappers for score heatmaps, feature ROC and precision-recall curves, YAML pipeline
runs and whole-model predictions on local files.

Plots are drawn with selene's own functions on the matplotlib Agg backend. Score heatmaps are
written to an explicit PNG path. selene's ROC and precision-recall functions write
``roc_curves.svg`` and ``precision_recall_curves.svg`` into the output directory they are given,
so those wrappers return the SVG path they wrote. Three plot functions are not wrapped because
the installed packages break them (see the skipped rows in the manifest): the interpret variant
scatter plot (plotly rejects its layout) and the interpret sequence logo (matplotlib rejects its
path effect).

Pipeline wrappers run selene's configuration-driven operations on a YAML file or a parsed dict,
with an explicit output directory. Prediction loads a whole torch model from a local file on
the CPU; that file is unpickled, so only trusted files should be used.

selene_sdk, numpy, matplotlib and torch are imported only when a function runs. No network is used.
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _clean,
    _dna_batch,
    _integer,
    _provenance,
    _text,
)
from dnakit.predictions.selene_sdk_functions import _is_sequence, _number, _rows

_SELENE_EXTRA = "bio-selene"
_BASES = ("A", "C", "G", "T")
_RUN_OPERATIONS = ("train", "evaluate", "analyze")
_PLOT_DPI = 150
_MAX_POSITIONS = 5_000
_MAX_SAMPLES = 1_000_000
_MAX_FEATURES = 1_000
_MAX_PREDICT_BATCH = 256
_MAX_PREDICT_LENGTH = 10_000
_MAX_OPERATIONS = 16
_MAX_PATH = 4_096
_MAX_COUNT = 10**9
_ROC_STYLE = "seaborn-v0_8-colorblind"
_UPSTREAM_ERRORS = (KeyError, TypeError, ValueError)


def _numpy() -> Any:
    return require("numpy", _SELENE_EXTRA)


def _pyplot() -> Any:
    matplotlib = require("matplotlib", _SELENE_EXTRA)
    matplotlib.use("Agg")
    return require("matplotlib.pyplot", _SELENE_EXTRA)


def _genome() -> Any:
    return require("selene_sdk.sequences", _SELENE_EXTRA).Genome


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _new_output(value: object, name: str, suffix: str) -> Path:
    """A new file with the given suffix; its parent directory is created."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path string.")
    target = Path(_text(str(value), name, maximum=_MAX_PATH)).expanduser().resolve()
    if target.suffix.lower() != suffix:
        raise ConfigurationError(f"{name} must end with {suffix}.")
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _directory(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a directory path string.")
    target = Path(_text(str(value), name, maximum=_MAX_PATH)).expanduser().resolve()
    if target.exists() and not target.is_dir():
        raise ConfigurationError(f"{name} exists and is not a directory: {target}")
    target.mkdir(parents=True, exist_ok=True)
    return target


def _subdirectories(path: Path) -> set[str]:
    return {child.name for child in path.iterdir() if child.is_dir()}


def _feature_table(values: object, name: str) -> list[list[float]]:
    """Sample rows of per-feature numbers; at least two samples and one feature."""

    if not _is_sequence(values) or len(values) < 2:
        raise ConfigurationError(f"{name} must be a list of at least two sample rows.")
    if len(values) > _MAX_SAMPLES:
        raise ConfigurationError(f"{name} may contain at most {_MAX_SAMPLES} rows.")
    width: int | None = None
    table: list[list[float]] = []
    for index, row in enumerate(values):
        label = f"{name}[{index}]"
        if not _is_sequence(row) or not row:
            raise ConfigurationError(f"{label} must be a non-empty list of numbers.")
        if width is None:
            width = len(row)
            if width > _MAX_FEATURES:
                raise ConfigurationError(f"{name} may have at most {_MAX_FEATURES} features.")
        if len(row) != width:
            raise ConfigurationError(f"{name} rows must all have {width} features.")
        table.append([_number(item, label) for item in row])
    return table


def _operations(value: object) -> list[str]:
    if not _is_sequence(value) or len(value) > _MAX_OPERATIONS:
        raise ConfigurationError(f"operations must be a list of up to {_MAX_OPERATIONS} names.")
    return [_choice(item, "operations", _RUN_OPERATIONS) for item in value]


def _config_copy(configs: object) -> dict[str, Any]:
    if not isinstance(configs, Mapping) or not all(isinstance(key, str) for key in configs):
        raise ConfigurationError("configs must be a mapping with string keys.")
    return copy.deepcopy(dict(configs))


def _learning_rate(value: object) -> float | None:
    if value is None:
        return None
    rate = _number(value, "lr")
    if rate <= 0:
        raise ConfigurationError("lr must be positive.")
    return rate


def _check_run_keys(configs: Mapping[str, Any], operations: list[str], lr: float | None) -> None:
    if "train" in operations and "lr" not in configs and lr is None:
        raise ConfigurationError("A learning rate is required when operations include train.")
    model = configs.get("model")
    if not isinstance(model, Mapping) or not isinstance(model.get("path"), str):
        raise ConfigurationError("configs must name the model file as model: {path: ...}.")
    model_path = Path(model["path"]).expanduser()
    if not model_path.exists():
        raise ConfigurationError(f"Model file or directory does not exist: {model_path}")


def _curve_metric(kind: str, truth: Sequence[float], scores: Sequence[float]) -> float:
    metrics = require("sklearn.metrics", _SELENE_EXTRA)
    if kind == "roc":
        return float(metrics.roc_auc_score(truth, scores))
    return float(metrics.average_precision_score(truth, scores))


def _draw_curves(
    kind: str,
    prediction: list[list[float]],
    truth: list[list[float]],
    directory: Path,
    *,
    minimum: int,
    style: str,
    fig_title: str,
) -> None:
    """Call selene's curve writer; restore the backend and close the figures it leaves open."""

    numpy = _numpy()
    module = require("selene_sdk.utils.performance_metrics", _SELENE_EXTRA)
    matplotlib = require("matplotlib", _SELENE_EXTRA)
    plt = _pyplot()
    existing = set(plt.get_fignums())
    upstream = (
        module.visualize_roc_curves if kind == "roc" else module.visualize_precision_recall_curves
    )
    try:
        with matplotlib.rc_context():
            upstream(
                numpy.asarray(prediction, dtype=float),
                numpy.asarray(truth, dtype=float),
                str(directory),
                report_gt_feature_n_positives=minimum,
                style=style,
                fig_title=fig_title,
                dpi=_PLOT_DPI,
            )
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"selene could not draw the {kind} curves: {exc}") from exc
    finally:
        _pyplot()
        for number in plt.get_fignums():
            if number not in existing:
                plt.close(number)


def _curves(
    kind: str,
    prediction: object,
    target: object,
    output_dir: object,
    *,
    report_gt_feature_n_positives: int,
    style: str,
    fig_title: str,
    operation: str,
) -> dict[str, Any]:
    pred = _feature_table(prediction, "prediction")
    truth = _feature_table(target, "target")
    if len(pred) != len(truth) or len(pred[0]) != len(truth[0]):
        raise ConfigurationError(
            "prediction and target must have the same number of rows and features."
        )
    if any(value not in (0.0, 1.0) for row in truth for value in row):
        raise ConfigurationError("target must hold only 0 and 1 labels.")
    minimum = _integer(
        report_gt_feature_n_positives,
        "report_gt_feature_n_positives",
        minimum=0,
        maximum=_MAX_COUNT,
    )
    style_name = _text(style, "style", maximum=128)
    if not isinstance(fig_title, str) or len(fig_title) > 200 or "\x00" in fig_title:
        raise ConfigurationError("fig_title must be a string of at most 200 characters.")
    directory = _directory(output_dir, "output_dir")
    filename = "roc_curves.svg" if kind == "roc" else "precision_recall_curves.svg"
    written = directory / filename
    if written.exists():
        raise ConfigurationError(f"Output file already exists: {written}")

    feature_count = len(pred[0])
    positives = [sum(1 for row in truth if row[index] == 1.0) for index in range(feature_count)]
    # Same rule as selene: both labels present and more positives than the threshold.
    plotted = [
        index
        for index in range(feature_count)
        if 0 < positives[index] < len(truth) and positives[index] > minimum
    ]
    values = [
        _curve_metric(kind, [row[index] for row in truth], [row[index] for row in pred])
        for index in plotted
    ]
    _draw_curves(
        kind, pred, truth, directory, minimum=minimum, style=style_name, fig_title=fig_title
    )
    if not written.is_file():
        raise BackendExecutionError(f"selene did not write {written}")
    return {
        "status": "completed",
        "operation": operation,
        "path": str(written),
        "image_format": "svg",
        "data": {
            "sample_count": len(pred),
            "feature_count": feature_count,
            "report_gt_feature_n_positives": minimum,
            "positives_per_feature": positives,
            "plotted_features": plotted,
            "metric": "roc_auc" if kind == "roc" else "average_precision",
            "metric_per_plotted_feature": values,
            "style": style_name,
        },
        "provenance": _provenance(
            "selene-sdk", figure="SVG written by selene_sdk.utils.performance_metrics"
        ),
    }


def selene_score_heatmap(
    score_matrix: Sequence[Sequence[float]],
    output_path: str,
    *,
    mask: Sequence[Sequence[bool]] | None = None,
    cmap: str = "Blues",
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw an (L, 4) in silico mutagenesis score matrix as a heatmap PNG.

    Rows are positions and columns are A, C, G, T; the figure shows the bases from top to bottom
    as T, G, C, A, as selene does. ``mask`` is an optional (L, 4) grid of True cells to hide.
    """

    rows = _rows(score_matrix, "score_matrix", limit=_MAX_POSITIONS)
    mask_rows: list[list[bool]] | None = None
    if mask is not None:
        if not _is_sequence(mask) or len(mask) != len(rows):
            raise ConfigurationError("mask must have the same number of rows as score_matrix.")
        mask_rows = []
        for index, row in enumerate(mask):
            if not _is_sequence(row) or len(row) != len(_BASES):
                raise ConfigurationError(f"mask[{index}] must hold four booleans.")
            if any(not isinstance(item, bool) for item in row):
                raise ConfigurationError(f"mask[{index}] must hold only True or False.")
            mask_rows.append(list(row))
    colour_map = _text(cmap, "cmap", maximum=64)
    target = _new_output(output_path, "output_path", ".png")
    numpy = _numpy()
    plt = _pyplot()
    vis = require("selene_sdk.interpret.vis", _SELENE_EXTRA)
    mask_array = None if mask_rows is None else numpy.asarray(mask_rows, dtype=bool)
    figure, ax = plt.subplots(figsize=(min(max(4.0, 0.4 * len(rows)), 40.0), 2.8))
    try:
        vis.heatmap(
            numpy.asarray(rows, dtype=float),
            mask=mask_array,
            sequence_type=_genome(),
            ax=ax,
            cmap=colour_map,
        )
        figure.savefig(target, format="png", dpi=_PLOT_DPI)
    except (ValueError, TypeError, KeyError, IndexError) as exc:
        raise ConfigurationError(f"selene could not draw the heatmap: {exc}") from exc
    finally:
        plt.close(figure)
    return {
        "status": "completed",
        "operation": "selene_score_heatmap",
        "path": str(target),
        "image_format": "png",
        "data": {
            "position_count": len(rows),
            "cmap": colour_map,
            "row_labels_top_to_bottom": list(reversed(_BASES)),
            "score_matrix": rows,
            "mask": mask_rows,
        },
        "provenance": _provenance(
            "selene-sdk", convention="rows are positions; columns A, C, G, T"
        ),
    }


def selene_roc_curves(
    prediction: Sequence[Sequence[float]],
    target: Sequence[Sequence[int]],
    output_dir: str,
    *,
    report_gt_feature_n_positives: int = 50,
    style: str = _ROC_STYLE,
    fig_title: str = "Feature ROC curves",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one ROC curve per model feature, as selene draws them, to ``output_dir``.

    ``prediction`` and ``target`` are (samples, features) tables; ``target`` holds 0 and 1.
    A feature is drawn only when both labels occur and it has more than
    ``report_gt_feature_n_positives`` positives. The file is ``roc_curves.svg``. ``style``
    defaults to ``seaborn-v0_8-colorblind`` because selene's default style name was removed
    from matplotlib. The result gives the plotted features and their ROC AUC.
    """

    return _curves(
        "roc",
        prediction,
        target,
        output_dir,
        report_gt_feature_n_positives=report_gt_feature_n_positives,
        style=style,
        fig_title=fig_title,
        operation="selene_roc_curves",
    )


def selene_precision_recall_curves(
    prediction: Sequence[Sequence[float]],
    target: Sequence[Sequence[int]],
    output_dir: str,
    *,
    report_gt_feature_n_positives: int = 50,
    style: str = _ROC_STYLE,
    fig_title: str = "Feature precision-recall curves",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one precision-recall curve per model feature, as selene draws them, to ``output_dir``.

    The rules match ``selene_roc_curves``; the file is ``precision_recall_curves.svg``. The
    result gives the plotted features and their average precision.
    """

    return _curves(
        "pr",
        prediction,
        target,
        output_dir,
        report_gt_feature_n_positives=report_gt_feature_n_positives,
        style=style,
        fig_title=fig_title,
        operation="selene_precision_recall_curves",
    )


def run_selene_config_file(
    config_path: str,
    output_dir: str,
    *,
    create_subdirectory: bool = True,
    lr: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the selene pipeline described by a YAML file, writing its outputs under ``output_dir``.

    The file is loaded with selene's ``load_path``; its ``output_dir`` key, if any, is replaced by
    ``output_dir``. The run then goes through selene's ``parse_configs_and_run``, which is what
    selene's ``load_and_parse_configs_and_run`` calls after it reads the file, except that it
    writes a re-dumped ``configuration.yaml`` instead of a copy of the original text.
    ``create_subdirectory`` puts the run in a new timestamped directory. ``lr`` is needed only
    when the pipeline trains. The result names the run directory and the files it holds.
    """

    path = input_path(config_path)
    flag = _flag(create_subdirectory, "create_subdirectory")
    rate = _learning_rate(lr)
    loader = require("selene_sdk.utils.config_utils", _SELENE_EXTRA)
    try:
        loaded = loader.load_path(str(path), instantiate=False)
    except (OSError, ValueError, TypeError) as exc:
        raise ConfigurationError(f"selene could not read the configuration {path}: {exc}") from exc
    if not isinstance(loaded, Mapping):
        raise ConfigurationError("the configuration file must hold a mapping at the top level.")
    configs = _config_copy(loaded)
    return _run_parsed(
        configs,
        output_dir,
        create_subdirectory=flag,
        lr=rate,
        source=str(path),
        operation="run_selene_config_file",
    )


def run_selene_configs(
    configs: Mapping[str, Any],
    output_dir: str,
    *,
    create_subdirectory: bool = True,
    lr: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a selene pipeline from a parsed configuration dict, writing under ``output_dir``.

    This is selene's ``parse_configs_and_run``. The dict is copied first, so the caller's object
    is not changed. Its ``ops`` list names the operations, and ``model: {path: ...}`` names the
    model file that selene copies into the run directory. The result names the run directory.
    """

    flag = _flag(create_subdirectory, "create_subdirectory")
    rate = _learning_rate(lr)
    settings = _config_copy(configs)
    return _run_parsed(
        settings,
        output_dir,
        create_subdirectory=flag,
        lr=rate,
        source="configs dict",
        operation="run_selene_configs",
    )


def _run_parsed(
    configs: dict[str, Any],
    output_dir: object,
    *,
    create_subdirectory: bool,
    lr: float | None,
    source: str,
    operation: str,
) -> dict[str, Any]:
    output = _directory(output_dir, "output_dir")
    operations = _operations(configs.get("ops"))
    _check_run_keys(configs, operations, lr)
    configs["output_dir"] = str(output)
    module = require("selene_sdk.utils.config_utils", _SELENE_EXTRA)
    before = _subdirectories(output)
    try:
        module.parse_configs_and_run(configs, create_subdirectory=create_subdirectory, lr=lr)
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"selene could not run the configuration: {exc}") from exc
    created = sorted(_subdirectories(output) - before)
    run_directory = output / created[0] if create_subdirectory and len(created) == 1 else output
    return {
        "status": "completed",
        "operation": operation,
        "source": source,
        "operations": operations,
        "output_dir": str(output),
        "run_directory": str(run_directory),
        "files": sorted(child.name for child in run_directory.iterdir() if child.is_file()),
        "provenance": _provenance("selene-sdk", pipeline="configuration file or dict, CPU run"),
    }


def run_selene_operations(
    operations: Sequence[str],
    configs: Mapping[str, Any],
    output_dir: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run selene's configuration-defined operations (train, evaluate, analyze) in order.

    This is selene's ``execute``: ``configs`` holds the sections the operations read, and
    ``output_dir`` receives their outputs. ``train`` needs ``lr`` in ``configs``. Operations
    are run as selene defines them, so each needs the sections it reads.
    """

    names = _operations(operations)
    settings = _config_copy(configs)
    if "train" in names and "lr" not in settings:
        raise ConfigurationError("configs must hold lr when operations include train.")
    output = _directory(output_dir, "output_dir")
    module = require("selene_sdk.utils.config_utils", _SELENE_EXTRA)
    try:
        module.execute(names, settings, str(output))
    except _UPSTREAM_ERRORS as exc:
        raise ConfigurationError(f"selene could not run the operations: {exc}") from exc
    return {
        "status": "completed",
        "operation": "run_selene_operations",
        "operations": names,
        "output_dir": str(output),
        "provenance": _provenance("selene-sdk", pipeline="configuration-defined operations"),
    }


def predict_selene_sequences(
    model_path: str,
    sequences: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a whole trained selene-style torch model on DNA sequences on the CPU.

    ``model_path`` must be a file written with ``torch.save(model)`` holding a whole
    ``torch.nn.Module``; state dicts are not accepted. The file is unpickled, so only trusted
    files should be used. Sequences must share one length; they are one-hot encoded by selene
    (N is a uniform 0.25 row) and passed to selene's ``predict`` as (batch, length, 4). The
    result holds the (batch, features) predictions.
    """

    path = input_path(model_path)
    items = _dna_batch(sequences, "sequences")
    if len(items) > _MAX_PREDICT_BATCH:
        raise ConfigurationError(f"sequences may contain at most {_MAX_PREDICT_BATCH} items.")
    width = len(items[0])
    if any(len(item) != width for item in items):
        raise ConfigurationError("sequences must all have the same length.")
    if width > _MAX_PREDICT_LENGTH:
        raise ConfigurationError(f"sequences may be at most {_MAX_PREDICT_LENGTH} bases long.")
    torch = require("torch", _SELENE_EXTRA)
    try:
        model = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:  # torch raises several unrelated types for unreadable files
        raise ConfigurationError(f"model file could not be loaded: {exc}") from exc
    if not isinstance(model, torch.nn.Module):
        raise ConfigurationError(
            "model file must hold a whole torch.nn.Module (torch.save(model)), not a state dict."
        )
    model.eval()
    numpy = _numpy()
    genome = _genome()
    batch = numpy.stack(
        [numpy.asarray(genome.sequence_to_encoding(item), dtype=numpy.float32) for item in items]
    )
    predict_module = require("selene_sdk.predict.model_predict", _SELENE_EXTRA)
    try:
        output = predict_module.predict(model, batch, use_cuda=False)
    except (RuntimeError, ValueError, TypeError) as exc:
        raise BackendExecutionError(f"selene prediction failed: {exc}") from exc
    array = numpy.asarray(output)
    return {
        "status": "completed",
        "operation": "predict_selene_sequences",
        "model_path": str(path),
        "sequence_count": len(items),
        "length": width,
        "shape": list(array.shape),
        "predictions": _clean(array),
        "provenance": _provenance("selene-sdk", device="cpu", mode="eval"),
    }


__all__ = [
    "predict_selene_sequences",
    "run_selene_config_file",
    "run_selene_configs",
    "run_selene_operations",
    "selene_precision_recall_curves",
    "selene_roc_curves",
    "selene_score_heatmap",
]

_register_public_verbose(globals())
