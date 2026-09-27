"""CSV-driven comparison plots for summarized and per-fold metrics."""

from __future__ import annotations

import csv
import math
import statistics
from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError

from .figure import save_figure


@dataclass(frozen=True, slots=True)
class _MetricPoint:
    source: str
    model: str
    mean: float
    se: float | None


def _read_csv(path: str | Path, required: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    source = Path(path)
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        headers = reader.fieldnames
        if not headers or len(set(headers)) != len(headers) or not required <= set(headers):
            raise InputFormatError(
                f"{source} must have unique headers including {', '.join(sorted(required))}."
            )
        rows = []
        try:
            for index, row in enumerate(reader, start=1):
                if index > 1000:
                    raise ConfigurationError("Metric comparison CSV exceeds 1000 rows.")
                if None in row or any(value is None for value in row.values()):
                    raise InputFormatError(f"Malformed CSV row {index + 1} in {source}.")
                rows.append(row)
        except csv.Error as exc:
            raise InputFormatError(f"Malformed CSV in {source}: {exc}.") from exc
    return headers, rows


def _number(value: str, name: str, source: str | Path) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise InputFormatError(f"{name} in {source} must be a finite number.") from exc
    if not math.isfinite(number):
        raise InputFormatError(f"{name} in {source} must be a finite number.")
    return number


def _label(value: str, name: str, source: str | Path) -> str:
    label = value.strip()
    if not label or len(label) > 256:
        raise InputFormatError(f"{name} in {source} must contain 1 to 256 characters.")
    return label


def plot_metric_comparison_csv(
    source: str | Path,
    *,
    output: str | Path,
    fold_csvs: Mapping[str, str | Path] | None = None,
    test_csv: str | Path | None = None,
    test_model_names: Mapping[str, str] | None = None,
    dataset: str | None = None,
    metric_column: str = "mcc",
    metric_label: str | None = None,
    fold_group: str = "Cross-validation",
    test_group: str = "Test",
    title: str | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> Path:
    """Read metric CSVs, compare their rows by source, and save a figure.

    The main CSV needs ``model,mean`` and optionally ``se,source,dataset``.
    With only this CSV, one call is enough to draw and save a comparison plot.
    Optional fold CSVs map displayed model names to files containing one
    ``metric_column`` value per fold; their sample standard error is computed.
    Optional test CSV rows need ``method,metric_column`` and, when filtering,
    ``dataset``. ``test_model_names`` maps method names to displayed names.
    Each source is drawn in a separate panel with a shared metric axis.
    """
    try:
        from matplotlib.figure import Figure
    except ImportError as exc:
        raise BackendUnavailableError(
            "Metric comparison plots require Matplotlib.",
            hint="Install DNAKit with the 'bio-plots' extra.",
        ) from exc

    if not isinstance(metric_column, str) or not metric_column.strip():
        raise ConfigurationError("metric_column must be a non-empty string.")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be a boolean.")
    headers, reference = _read_csv(source, {"model", "mean"})
    if dataset is not None and "dataset" not in headers:
        raise InputFormatError("dataset filtering requires a dataset column in the main CSV.")
    points: list[_MetricPoint] = []
    for row in reference:
        if dataset is not None and row["dataset"] != dataset:
            continue
        se = _number(row["se"], "se", source) if row.get("se", "").strip() else None
        if se is not None and se < 0:
            raise InputFormatError(f"se in {source} cannot be negative.")
        points.append(
            _MetricPoint(
                _label(row.get("source", "Reference"), "source", source),
                _label(row["model"], "model", source),
                _number(row["mean"], "mean", source),
                se,
            )
        )
    if dataset is not None and not points:
        raise InputFormatError(f"No reference metric rows for {dataset!r} in {source}.")

    for model, path in (fold_csvs or {}).items():
        _headers, rows = _read_csv(path, {metric_column})
        if not rows:
            raise InputFormatError(f"No fold metric rows in {path}.")
        values = [_number(row[metric_column], metric_column, path) for row in rows]
        se = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else None
        points.append(
            _MetricPoint(
                _label(fold_group, "fold_group", path),
                f"{_label(model, 'model', path)} (n={len(values)})",
                statistics.fmean(values),
                se,
            )
        )

    if test_csv is not None:
        required = {"method", metric_column}
        if dataset is not None:
            required.add("dataset")
        _headers, rows = _read_csv(test_csv, required)
        selected = [row for row in rows if dataset is None or row["dataset"] == dataset]
        if not selected:
            raise InputFormatError(f"No test metric rows for {dataset!r} in {test_csv}.")
        for row in selected:
            method = _label(row["method"], "method", test_csv)
            model = method if test_model_names is None else test_model_names.get(method, method)
            points.append(
                _MetricPoint(
                    _label(test_group, "test_group", test_csv),
                    _label(model, "model", test_csv),
                    _number(row[metric_column], metric_column, test_csv),
                    None,
                )
            )

    if not points:
        raise InputFormatError("No metric rows matched the requested dataset.")
    groups: OrderedDict[str, list[_MetricPoint]] = OrderedDict()
    for point in points:
        groups.setdefault(point.source, []).append(point)
    widths = [max(2, len(group)) for group in groups.values()]
    figure = Figure(figsize=(max(8, 2 + sum(widths)), 6))
    axes = figure.subplots(
        1, len(groups), sharey=True, squeeze=False, gridspec_kw={"width_ratios": widths}
    )[0]
    for axis, (group_name, group_points) in zip(axes, groups.items(), strict=True):
        for index, point in enumerate(group_points):
            axis.errorbar(index, point.mean, yerr=point.se, fmt="o", capsize=3)
        axis.set(
            xticks=range(len(group_points)),
            xticklabels=[point.model for point in group_points],
            title=group_name,
            ylabel=metric_label or metric_column.upper(),
        )
        axis.tick_params(axis="x", rotation=75)
    figure.suptitle(title or dataset or Path(source).stem)
    return save_figure(figure, output, overwrite=overwrite)


__all__ = ["plot_metric_comparison_csv"]

_register_public_verbose(globals())
