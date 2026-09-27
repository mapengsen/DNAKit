"""Learning curves from one epoch-level CSV file."""

from __future__ import annotations

import csv
import math
import statistics
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError

from .figure import save_figure

_METRICS = ("train_loss", "validation_loss", "train_mcc", "validation_mcc")
_OPTIONAL_GROUPS = ("dataset", "setting", "method", "fold")
_COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")


@dataclass(frozen=True, slots=True)
class _EpochRow:
    dataset: str
    setting: str
    method: str
    fold: str | None
    epoch: int
    values: dict[str, float | None]


def _number(value: str, *, column: str, line: int) -> float | None:
    if not value.strip():
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise InputFormatError(f"{column} in CSV row {line} must be a finite number.") from exc
    if not math.isfinite(result):
        raise InputFormatError(f"{column} in CSV row {line} must be a finite number.")
    return result


def _read_epochs(source: str | Path) -> OrderedDict[tuple[str, str], list[_EpochRow]]:
    path = Path(source)
    groups: OrderedDict[tuple[str, str], list[_EpochRow]] = OrderedDict()
    seen_folds: set[tuple[str, str, str, str, int]] = set()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        headers = reader.fieldnames
        if (
            not headers
            or len(headers) != len(set(headers))
            or not {"epoch", *_METRICS} <= set(headers)
        ):
            raise InputFormatError(
                "Learning-curve CSV needs unique headers including epoch, "
                "train_loss, validation_loss, train_mcc, and validation_mcc."
            )
        try:
            for line, record in enumerate(reader, start=2):
                if None in record or any(value is None for value in record.values()):
                    raise InputFormatError(f"Malformed learning-curve CSV row {line}.")
                for column in _OPTIONAL_GROUPS:
                    if column in headers and not record[column].strip():
                        raise InputFormatError(f"{column} in CSV row {line} cannot be blank.")
                dataset = record.get("dataset", path.stem).strip()
                setting = record.get("setting", "").strip()
                method = record.get("method", "model").strip()
                fold = record["fold"].strip() if "fold" in headers else None
                if any(len(value) > 128 for value in (dataset, setting, method)):
                    raise InputFormatError(f"Grouping label in CSV row {line} is too long.")
                try:
                    epoch = int(record["epoch"].strip())
                except ValueError as exc:
                    raise InputFormatError(f"epoch in CSV row {line} must be an integer.") from exc
                if epoch < 0:
                    raise InputFormatError(f"epoch in CSV row {line} cannot be negative.")
                values = {name: _number(record[name], column=name, line=line) for name in _METRICS}
                if all(value is None for value in values.values()):
                    raise InputFormatError(f"CSV row {line} has no metric values.")
                if fold is not None:
                    key = (dataset, setting, method, fold, epoch)
                    if key in seen_folds:
                        raise InputFormatError(f"Duplicate fold and epoch in CSV row {line}.")
                    seen_folds.add(key)
                groups.setdefault((dataset, setting), []).append(
                    _EpochRow(dataset, setting, method, fold, epoch, values)
                )
        except csv.Error as exc:
            raise InputFormatError(f"Malformed learning-curve CSV: {exc}.") from exc
    if not groups:
        raise InputFormatError("Learning-curve CSV has no epoch rows.")
    return groups


def _filename_label(value: str) -> str:
    cleaned = "".join(
        character if character.isalnum() or character in "-_." else "_" for character in value
    )
    return cleaned.strip("._")[:80] or "plot"


def plot_learning_curves_csv(
    source: str | Path,
    *,
    output: str | Path,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[tuple[str, str], Path]:
    """Plot mean training/validation loss and MCC by epoch from one CSV.

    Required columns are ``epoch``, ``train_loss``, ``validation_loss``,
    ``train_mcc``, and ``validation_mcc``. Optional ``dataset`` and ``setting``
    columns create a separate figure per pair. ``method`` separates curves;
    ``fold`` identifies replicates for mean and sample-SD bands. Metric cells
    may be blank, but at least one value is needed for each panel in a group.
    For one group, ``output`` may be an image filename or a directory. Multiple
    groups require a directory and produce one SVG per group.
    """
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be a boolean.")
    groups = _read_epochs(source)
    target = Path(output)
    if target.suffix and len(groups) != 1:
        raise ConfigurationError("Multiple learning-curve groups require an output directory.")
    if target.suffix:
        destinations = {next(iter(groups)): target}
    else:
        destinations = {
            key: target
            / (
                f"{_filename_label(key[0])}_{_filename_label(key[1])}.svg"
                if key[1]
                else f"{_filename_label(key[0])}.svg"
            )
            for key in groups
        }
        if len(set(destinations.values())) != len(destinations):
            raise ConfigurationError("Learning-curve group names produce duplicate filenames.")

    for key, records in groups.items():
        for metric, label in (("loss", "Loss"), ("mcc", "MCC")):
            if not any(
                row.values[f"{split}_{metric}"] is not None
                for row in records
                for split in ("train", "validation")
            ):
                raise InputFormatError(f"No {label} values for learning-curve group {key!r}.")

    try:
        from matplotlib import rc_context
        from matplotlib.figure import Figure
        from matplotlib.ticker import MaxNLocator
    except ImportError as exc:
        raise BackendUnavailableError(
            "Learning-curve plots require Matplotlib.",
            hint="Install DNAKit with the 'bio-plots' extra.",
        ) from exc

    with rc_context({"svg.fonttype": "none", "svg.hashsalt": "dnakit-learning-curves"}):
        for key, records in groups.items():
            methods = list(dict.fromkeys(row.method for row in records))
            figure = Figure(figsize=(12, 4.5))
            axes = figure.subplots(1, 2, squeeze=False)[0]
            for axis, metric, label in zip(axes, ("loss", "mcc"), ("Loss", "MCC"), strict=True):
                for method_index, method in enumerate(methods):
                    color = _COLORS[method_index % len(_COLORS)]
                    for split, style in (("train", "-"), ("validation", "--")):
                        column = f"{split}_{metric}"
                        by_epoch: dict[int, list[float]] = defaultdict(list)
                        for row in records:
                            value = row.values[column]
                            if row.method == method and value is not None:
                                by_epoch[row.epoch].append(value)
                        if not by_epoch:
                            continue
                        epochs = sorted(by_epoch)
                        means = [statistics.fmean(by_epoch[epoch]) for epoch in epochs]
                        deviations = [
                            statistics.stdev(by_epoch[epoch]) if len(by_epoch[epoch]) > 1 else 0.0
                            for epoch in epochs
                        ]
                        axis.plot(
                            epochs, means, label=f"{method}/{split}", color=color, linestyle=style
                        )
                        axis.fill_between(
                            epochs,
                            [mean - sd for mean, sd in zip(means, deviations, strict=True)],
                            [mean + sd for mean, sd in zip(means, deviations, strict=True)],
                            color=color,
                            alpha=0.12,
                        )
                axis.set(xlabel="Epoch", ylabel=label)
                axis.xaxis.set_major_locator(MaxNLocator(integer=True))
                axis.legend(fontsize=8)
            heading = "/".join(part for part in key if part)
            figure.suptitle(f"{heading} — mean ± SD")
            save_figure(figure, destinations[key], overwrite=overwrite)
    return destinations


__all__ = ["plot_learning_curves_csv"]

_register_public_verbose(globals())
