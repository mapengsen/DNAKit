"""Project DNA embeddings and highlight errors supplied by a CSV column."""

from __future__ import annotations

import csv
import math
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics.features import analyze_features, cluster_features
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError
from dnakit.io import load_arrays, save_arrays, write_rows

from .figure import save_figure

_METHOD_NAME = re.compile(r"^[A-Za-z0-9_-]+$")
_COLORS = ("#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1")


def _text(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _error_kind(value: str) -> str:
    normalized = value.strip().casefold().replace("-", "_").replace(" ", "_")
    if normalized in {"correct", "0", "false", "no"}:
        return "correct"
    if normalized in {"fp", "false_positive"}:
        return "FP"
    if normalized in {"fn", "false_negative"}:
        return "FN"
    if normalized in {"error", "incorrect", "1", "true", "yes"}:
        return "error"
    raise InputFormatError(f"Unknown error value {value!r}; use correct, FP, FN, or error.")


def _read_error_csv(source: str | Path) -> dict[str, dict[str, dict[str, str]]]:
    path = Path(source)
    methods: dict[str, dict[str, dict[str, str]]] = OrderedDict()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        headers = reader.fieldnames
        if not headers or len(set(headers)) != len(headers) or not {"id", "error"} <= set(headers):
            raise InputFormatError("Error CSV needs unique headers including id and error.")
        try:
            for index, row in enumerate(reader, start=1):
                if index > 1_000_000:
                    raise ConfigurationError("Error CSV exceeds one million rows.")
                if None in row or any(value is None for value in row.values()):
                    raise InputFormatError(f"Malformed error CSV row {index + 1}.")
                identifier = row["id"].strip()
                method = row.get("method", "errors").strip()
                if not identifier or not _METHOD_NAME.fullmatch(method):
                    raise InputFormatError(f"Invalid id or method in error CSV row {index + 1}.")
                row["error"] = _error_kind(row["error"])
                group = methods.setdefault(method, {})
                if identifier in group:
                    raise InputFormatError(f"Duplicate id {identifier!r} for method {method!r}.")
                group[identifier] = row
        except csv.Error as exc:
            raise InputFormatError(f"Malformed error CSV: {exc}.") from exc
    if not methods:
        raise InputFormatError("Error CSV has no samples.")
    return methods


def _validate_prediction(row: dict[str, str], label: str, identifier: str) -> None:
    if row.get("label", "").strip() and row["label"].strip() != label:
        raise InputFormatError(f"Label for id {identifier!r} differs from the embedding file.")
    error = row["error"]
    if (error == "FP" and label == "1") or (error == "FN" and label == "0"):
        raise InputFormatError(f"Error value for id {identifier!r} contradicts its true label.")
    prediction = row.get("prediction", "").strip()
    if not prediction or label not in {"0", "1"}:
        return
    if prediction not in {"0", "1"}:
        raise InputFormatError(f"Prediction for id {identifier!r} must be 0 or 1.")
    if label == "0" and prediction == "1":
        expected = "FP"
    elif label == "1" and prediction == "0":
        expected = "FN"
    else:
        expected = "correct"
    if error != expected and not (error == "error" and expected != "correct"):
        raise InputFormatError(f"Error value for id {identifier!r} contradicts its prediction.")


def plot_embedding_errors(
    embeddings_file: str | Path,
    errors_csv: str | Path,
    *,
    output_dir: str | Path | None = None,
    title: str | None = None,
    n_clusters: int | None = None,
    save_details: bool = False,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Path]:
    """Plot errors from one NPZ embedding file and one CSV with ``id,error``.

    The NPZ needs ``embeddings,ids`` and may include ``splits,labels``. When
    train/test splits exist, scaling and PCA fit on train and plot test rows;
    otherwise all rows are projected. A CSV ``method`` column creates one plot
    per method. Error values may be correct, FP, FN, or a generic error flag.
    Optional clustering uses standardized test embeddings, not PCA coordinates.
    By default, plots are written to a directory named after the CSV stem.
    ``save_details`` writes projection, point, and cluster-summary files.
    """
    try:
        import numpy as np
        from matplotlib.figure import Figure
    except ImportError as exc:
        raise BackendUnavailableError(
            "Embedding error plots require NumPy and Matplotlib.",
            hint="Install DNAKit with the 'bio-analysis' and 'bio-plots' extras.",
        ) from exc

    data = load_arrays(embeddings_file)
    if not {"embeddings", "ids"} <= data.keys():
        raise InputFormatError("Embedding NPZ needs embeddings and ids arrays.")
    embeddings = np.asarray(data["embeddings"], dtype=float)
    if embeddings.ndim != 2 or embeddings.shape[1] < 2 or not np.isfinite(embeddings).all():
        raise InputFormatError("embeddings must be a finite matrix with at least two columns.")
    raw_ids = np.asarray(data["ids"])
    if raw_ids.ndim != 1:
        raise InputFormatError("Embedding ids must be a one-dimensional array.")
    ids = [_text(value) for value in raw_ids.tolist()]
    if len(ids) != len(embeddings) or len(set(ids)) != len(ids):
        raise InputFormatError("Embedding ids must be unique and match the matrix rows.")
    raw_splits = np.asarray(data["splits"]) if "splits" in data else None
    if raw_splits is not None and raw_splits.ndim != 1:
        raise InputFormatError("Embedding splits must be a one-dimensional array.")
    splits = (
        [_text(value) for value in raw_splits.tolist()]
        if raw_splits is not None
        else ["test"] * len(ids)
    )
    if len(splits) != len(ids):
        raise InputFormatError("Embedding splits must match the matrix rows.")
    if set(splits) - {"train", "test"}:
        raise InputFormatError("Embedding splits may contain only train and test.")
    has_train = "train" in splits
    train = (
        np.array([split == "train" for split in splits]) if has_train else np.ones(len(ids), bool)
    )
    test = np.array([split == "test" for split in splits]) if has_train else np.ones(len(ids), bool)
    if train.sum() < 2 or not test.any():
        raise InputFormatError("At least two train and one test embedding are required.")
    raw_labels = np.asarray(data["labels"]) if "labels" in data else None
    if raw_labels is not None and raw_labels.ndim != 1:
        raise InputFormatError("Embedding labels must be a one-dimensional array.")
    labels = (
        [_text(value) for value in raw_labels.tolist()]
        if raw_labels is not None
        else [""] * len(ids)
    )
    if len(labels) != len(ids):
        raise InputFormatError("Embedding labels must match the matrix rows.")
    test_ids = [identifier for identifier, selected in zip(ids, test, strict=True) if selected]
    test_labels = [label for label, selected in zip(labels, test, strict=True) if selected]
    if len(set(test_ids)) != len(test_ids):
        raise InputFormatError("Test embedding ids must be unique.")

    methods = _read_error_csv(errors_csv)
    expected_ids = set(test_ids)
    for method, rows in methods.items():
        if set(rows) != expected_ids:
            missing = len(expected_ids - rows.keys())
            extra = len(rows.keys() - expected_ids)
            raise InputFormatError(
                f"Method {method!r} does not match test embeddings: "
                f"{missing} missing, {extra} extra ids."
            )

    folder = Path(output_dir) if output_dir is not None else Path(errors_csv).with_suffix("")
    paths = {method: folder / f"{method}.svg" for method in methods}
    detail_paths = (
        [folder / "projection.npz", *[folder / f"{method}_points.csv" for method in methods]]
        if save_details
        else []
    )
    if save_details and n_clusters is not None:
        detail_paths.append(folder / "cluster_errors.csv")
    if not overwrite:
        for path in [*paths.values(), *detail_paths]:
            if path.exists():
                raise FileExistsError(path)

    analysis = analyze_features(embeddings[train], parameters={"svd_solver": "full"})
    scaler = analysis.estimator.named_steps["scale"]
    pca = analysis.estimator.named_steps["analysis"]
    coordinates = analysis.transform(embeddings[test])
    clusters = None
    if n_clusters is not None:
        if (
            isinstance(n_clusters, bool)
            or not isinstance(n_clusters, int)
            or not 2 <= n_clusters <= len(test_ids)
        ):
            raise ConfigurationError("n_clusters must be between 2 and the test sample count.")
        clusters = cluster_features(
            scaler.transform(embeddings[test]),
            seed=42,
            parameters={"n_clusters": n_clusters, "n_init": 10},
        ).labels_

    dataset_name = folder.name if output_dir is not None else Path(embeddings_file).parent.name
    all_statistics = []
    points_by_method = {}
    for method, rows in methods.items():
        points = []
        for index, (identifier, label, xy) in enumerate(
            zip(test_ids, test_labels, coordinates, strict=True)
        ):
            row = rows[identifier]
            effective_label = label or row.get("label", "").strip()
            _validate_prediction(row, effective_label, identifier)
            probability = row.get("probability", "").strip()
            if probability:
                try:
                    numeric_probability = float(probability)
                except ValueError as exc:
                    raise InputFormatError(
                        f"Probability for id {identifier!r} must be a number."
                    ) from exc
                if not math.isfinite(numeric_probability) or not 0 <= numeric_probability <= 1:
                    raise InputFormatError(
                        f"Probability for id {identifier!r} must be between 0 and 1."
                    )
            points.append(
                dict(
                    id=identifier,
                    label=effective_label,
                    prediction=row.get("prediction", ""),
                    probability=probability,
                    error=row["error"],
                    cluster=int(clusters[index]) if clusters is not None else None,
                    pc1=float(xy[0]),
                    pc2=float(xy[1]),
                )
            )

        points_by_method[method] = points
        if clusters is not None:
            for cluster in range(n_clusters):
                group = [point for point in points if point["cluster"] == cluster]
                all_statistics.append(
                    dict(
                        dataset=dataset_name,
                        method=method,
                        cluster=cluster,
                        n=len(group),
                        positive=sum(point["label"] == "1" for point in group),
                        fp=sum(point["error"] == "FP" for point in group),
                        fn=sum(point["error"] == "FN" for point in group),
                        error_rate=(
                            sum(point["error"] != "correct" for point in group) / len(group)
                            if group
                            else None
                        ),
                    )
                )

    for method, points in points_by_method.items():
        figure = Figure(figsize=(6, 5))
        axis = figure.subplots()
        distinct_labels = list(dict.fromkeys(point["label"] for point in points))
        for index, label in enumerate(distinct_labels):
            group = [point for point in points if point["label"] == label]
            axis.scatter(
                [point["pc1"] for point in group],
                [point["pc2"] for point in group],
                c=_COLORS[index % len(_COLORS)],
                s=15,
                label=f"True {label}" if label else "Samples",
                alpha=0.6,
            )
        for error, marker in (("FP", "x"), ("FN", "+"), ("error", "x")):
            group = [point for point in points if point["error"] == error]
            if group:
                axis.scatter(
                    [point["pc1"] for point in group],
                    [point["pc2"] for point in group],
                    c="black",
                    marker=marker,
                    s=40,
                    label=error,
                )
        axis.set(
            title=f"{title or Path(embeddings_file).stem} / {method}",
            xlabel="PC1",
            ylabel="PC2",
        )
        axis.legend()
        save_figure(figure, paths[method], overwrite=overwrite)

        if save_details:
            write_rows(points, folder / f"{method}_points.csv", overwrite=overwrite)

    if save_details:
        arrays = dict(
            scaler_mean=scaler.mean_,
            scaler_scale=scaler.scale_,
            pca_mean=pca.mean_,
            components=pca.components_,
            explained_variance=pca.explained_variance_,
            train_ids=np.asarray(ids)[train],
            test_ids=np.asarray(test_ids),
            coordinates=coordinates,
        )
        if clusters is not None:
            arrays["clusters"] = clusters
        save_arrays(arrays, folder / "projection.npz", overwrite=overwrite)
        if clusters is not None:
            write_rows(all_statistics, folder / "cluster_errors.csv", overwrite=overwrite)
    return paths


__all__ = ["plot_embedding_errors"]

_register_public_verbose(globals())
