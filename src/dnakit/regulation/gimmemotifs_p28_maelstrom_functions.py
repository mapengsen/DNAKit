"""Motif-activity prediction, rank aggregation and heatmap functions of maelstrom (pack P28).

The functions work on precomputed tables: a region-by-motif score or count table, a label or value
table, or the output directory of an earlier maelstrom run. No genome is scanned and nothing is
downloaded. The end-to-end ``run_maelstrom`` pipeline is not wrapped (see the disposition file).
"""

from __future__ import annotations

import contextlib
import math
import random
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import _choose, _flag, _gimme_provenance, _integer

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_METHODS = {
    "bayesianridge": "BayesianRidge regression on motif scores (coefficients)",
    "xgboost": "XGBoost regression on motif scores (feature importances)",
    "mwu": "Mann-Whitney U test per cluster (-log10 BH-corrected p-value)",
    "rf": "random forest classifier (signed feature importances)",
    "multitasklasso": "multi-task lasso regression (coefficients)",
    "svr": "linear support vector regression (weights)",
}
_CLASSIFIERS = ("mwu", "rf")
_PREDICTORS = ("multitasklasso", "svr")
_RANK_METHODS = ("int_stouffer", "stuart")
_MAX_CELLS = 20_000_000


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except BaseException as exc:  # gimmemotifs, pandas and sklearn raise many exception types
        raise BackendExecutionError(f"gimmemotifs failed: {type(exc).__name__}: {exc}") from exc


@contextlib.contextmanager
def _seeded(seed: int | None) -> Iterator[None]:
    """Seed numpy's and Python's global generators for one call and restore them afterwards."""

    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    numpy_state = np.random.get_state()
    python_state = random.getstate()
    try:
        if seed is not None:
            np.random.seed(seed)
            random.seed(seed)
        yield
    finally:
        np.random.set_state(numpy_state)
        random.setstate(python_state)


def _seed(seed: object) -> int | None:
    if seed is None:
        return None
    return _integer(seed, "seed", 0, 2**32 - 1)


def _finite_rows(frame: Any) -> list[list[float | None]]:
    return [
        [(float(v) if math.isfinite(float(v)) else None) for v in row]
        for row in frame.to_numpy().tolist()
    ]


def _check_size(frame: Any, name: str) -> None:
    if frame.shape[0] == 0 or frame.shape[1] == 0 or frame.shape[0] * frame.shape[1] > _MAX_CELLS:
        raise ConfigurationError(f"{name} must be a non-empty table of at most {_MAX_CELLS} cells.")


def gimmemotifs_predict_motif_activity(
    labels_file: str | Path,
    motif_table: str | Path,
    *,
    method: str = "bayesianridge",
    predict: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict motif activity from a region-by-motif table with one maelstrom method.

    Calls ``gimmemotifs.maelstrom.moap``.

    ``method`` selects the predictor: ``bayesianridge`` (Bayesian ridge regression),
    ``xgboost`` (XGBoost regression; needs the ``xgboost`` package), ``mwu`` (Mann-Whitney U test),
    ``rf`` (random forest), ``multitasklasso`` (multi-task lasso) or ``svr`` (support vector
    regression). ``labels_file`` is a tab-separated table with a header, indexed by region name:
    ``mwu`` and ``rf`` need one column with the cluster label, the regression methods numeric
    columns. ``motif_table`` has the same index, one column per motif and the motif scores (or
    counts). Scores are scaled by the methods that do so upstream. The result has one row per motif
    and one column per label class or value column. With ``predict=True`` (``multitasklasso`` and
    ``svr`` only) the values the fitted model predicts for the regions are returned too. ``seed``
    makes the stochastic methods reproducible.
    """

    chosen = _choose(method, tuple(_METHODS), "method")
    wants_predictions = _flag(predict, "predict")
    if wants_predictions and chosen not in _PREDICTORS:
        raise ConfigurationError(f"predict=True is available for {', '.join(_PREDICTORS)} only.")
    if chosen == "xgboost":
        require("xgboost", "xgboost")
    moap_module = require("gimmemotifs.maelstrom.moap", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    pd = require("pandas", _GIMMEMOTIFS_EXTRA)
    labels = input_path(labels_file)
    table = input_path(motif_table)
    chosen_seed = _seed(seed)
    state = None if chosen_seed is None else np.random.RandomState(chosen_seed)
    with _seeded(chosen_seed):
        frame = _call(
            moap_module.moap,
            str(labels),
            method=chosen,
            scoring="score",
            outfile=None,
            motiffile=str(table),
            ncpus=1,
            zscore=False,
            gc=False,
            random_state=state,
            progress=False,
        )
    _check_size(frame, "activity table")
    columns = [str(name) for name in frame.columns]
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "predict_motif_activity",
        "method": chosen,
        "description": _METHODS[chosen],
        "seed": chosen_seed,
        "classes": columns,
        "motifs": [str(name) for name in frame.index.tolist()],
        "values": _finite_rows(frame),
        "provenance": _gimme_provenance(function="gimmemotifs.maelstrom.moap.moap", method=chosen),
    }
    if wants_predictions:
        scores = _call(pd.read_table, str(table), index_col=0, comment="#")
        values = _call(pd.read_table, str(labels), index_col=0, comment="#")
        values = values.loc[scores.index]
        model = _call(
            moap_module.Moap.create,
            chosen,
            ncpus=1,
            random_state=None if chosen_seed is None else np.random.RandomState(chosen_seed),
            disable_tqdm=True,
        )
        with _seeded(chosen_seed):
            _call(model.fit, scores, values)
            predicted = _call(model.predict, scores)
        result["predicted_regions"] = [str(name) for name in predicted.index.tolist()]
        result["predicted_columns"] = [str(name) for name in predicted.columns]
        result["predicted"] = _finite_rows(predicted)
    return result


def gimmemotifs_aggregate_motif_ranks(
    scores_table: str | Path | Mapping[str, Mapping[str, float]],
    *,
    method: str = "int_stouffer",
    include_reverse: bool = True,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Aggregate several motif-activity columns into one score per motif.

    Calls ``gimmemotifs.maelstrom.rank.rankagg``.

    ``scores_table`` is a tab-separated table with a header (rows are motifs, columns are the
    activities of different methods) or a mapping ``{motif: {column: value}}``. ``method`` is
    ``int_stouffer`` (inverse-normal transform of each column's ranks, combined with Stouffer's
    method into a z-score per motif) or ``stuart`` (robust rank aggregation, Kolde et al. 2012;
    ``include_reverse`` also aggregates the lowest ranks and subtracts them). Ties are broken at
    random; ``seed`` makes the result reproducible. Motifs are returned best first.
    """

    chosen = _choose(method, _RANK_METHODS, "method")
    reverse = _flag(include_reverse, "include_reverse")
    rank_module = require("gimmemotifs.maelstrom.rank", _GIMMEMOTIFS_EXTRA)
    pd = require("pandas", _GIMMEMOTIFS_EXTRA)
    if isinstance(scores_table, Mapping):
        frame = _call(
            pd.DataFrame.from_dict,
            {str(k): dict(v) for k, v in scores_table.items()},
            orient="index",
        )
    else:
        frame = _call(pd.read_table, str(input_path(scores_table)), index_col=0, comment="#")
    _check_size(frame, "scores_table")
    if frame.shape[1] < 1 or frame.isna().any().any():
        raise ConfigurationError(
            "scores_table must be complete: no missing values, at least one column."
        )
    chosen_seed = _seed(seed)
    with _seeded(chosen_seed):
        aggregated = _call(rank_module.rankagg, frame, method=chosen, include_reverse=reverse)
    series = aggregated.iloc[:, 0].sort_values(ascending=False, kind="stable")
    return {
        "status": "completed",
        "operation": "aggregate_motif_ranks",
        "method": chosen,
        "seed": chosen_seed,
        "value_name": str(aggregated.columns[0]),
        "n_motifs": int(frame.shape[0]),
        "n_columns": int(frame.shape[1]),
        "ranking": [{"motif": str(name), "value": float(value)} for name, value in series.items()],
        "provenance": _gimme_provenance(
            function="gimmemotifs.maelstrom.rank.rankagg", method=chosen
        ),
    }


def gimmemotifs_plot_maelstrom_heatmap(
    maelstrom_dir: str | Path,
    output_image: str | Path,
    *,
    kind: str = "final",
    min_freq: float = 0.01,
    threshold: float = 2.0,
    name: bool = True,
    indirect: bool = True,
    max_number_factors: int = 5,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the clustered heatmap of motif activity of a maelstrom run.

    Calls ``MaelstromResult.plot_heatmap``.

    ``maelstrom_dir`` is a maelstrom output directory: ``final.out.txt`` (rank-aggregated
    activity), ``motif.score.txt.gz``, ``motif.count.txt.gz`` and a PFM file, optionally
    ``motif.freq.txt`` and ``activity.*.txt``. ``kind`` is ``final`` (the aggregated result),
    ``freq`` (motif frequencies per cluster) or the name of one method's activity such as
    ``rf.score``. Motifs need an absolute activity of at least ``threshold`` and a frequency of at
    least ``min_freq``. With ``name`` the rows carry factor names instead of motif IDs (at most
    ``max_number_factors``; ``indirect`` adds predicted factors). Rows are ordered by Ward
    clustering with correlation distance. The PNG is written to ``output_image`` (must not exist);
    the plotted rows, columns and values are returned.
    """

    maelstrom = require("gimmemotifs.maelstrom", _GIMMEMOTIFS_EXTRA)
    matplotlib = require("matplotlib", _GIMMEMOTIFS_EXTRA)
    matplotlib.use("Agg")
    plt = require("matplotlib.pyplot", _GIMMEMOTIFS_EXTRA)
    directory = Path(maelstrom_dir).expanduser()
    if not directory.is_dir():
        raise ConfigurationError(f"maelstrom_dir is not a directory: {directory}")
    destination = Path(output_image).expanduser()
    if destination.exists():
        raise ConfigurationError(f"output_image already exists: {destination}")
    frequency = float(min_freq)
    cutoff = float(threshold)
    if not math.isfinite(frequency) or not math.isfinite(cutoff):
        raise ConfigurationError("min_freq and threshold must be finite numbers.")
    shown = _integer(max_number_factors, "max_number_factors", 1, 100)
    resolution = _integer(dpi, "dpi", 50, 600)
    result = _call(maelstrom.MaelstromResult, str(directory))
    allowed = ["final", "freq", *sorted(result.activity)]
    chosen = _choose(kind, allowed, "kind")
    grid = _call(
        result.plot_heatmap,
        kind=chosen,
        min_freq=frequency,
        threshold=cutoff,
        name=_flag(name, "name"),
        indirect=_flag(indirect, "indirect"),
        max_number_factors=shown,
    )
    if grid is None:
        plt.close("all")
        raise ConfigurationError("No motif passes the threshold and min_freq filters; lower them.")
    try:
        mesh = grid.collections[0]
        array = mesh.get_array()
        rows = [tick.get_text() for tick in grid.get_yticklabels()]
        columns = [tick.get_text() for tick in grid.get_xticklabels()]
        values = [
            [float(v) for v in row] for row in array.reshape(len(rows), len(columns)).tolist()
        ]
        destination.parent.mkdir(parents=True, exist_ok=True)
        grid.get_figure().savefig(str(destination), dpi=resolution)
    finally:
        plt.close("all")
    return {
        "status": "completed",
        "operation": "plot_maelstrom_heatmap",
        "kind": chosen,
        "rows": rows,
        "columns": columns,
        "values": values,
        "outputs": {"image": str(destination)},
        "provenance": _gimme_provenance(
            function="gimmemotifs.maelstrom.MaelstromResult.plot_heatmap"
        ),
    }


__all__ = [
    "gimmemotifs_aggregate_motif_ranks",
    "gimmemotifs_plot_maelstrom_heatmap",
    "gimmemotifs_predict_motif_activity",
]

_register_public_verbose(globals(), tuple(__all__))
