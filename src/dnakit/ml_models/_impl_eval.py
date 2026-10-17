"""Evaluation, calibration, uncertainty and table-export operations of ``ml_models``."""

from __future__ import annotations

import json
import math
import random
from collections.abc import Mapping, Sequence

from . import _impl_metrics as M
from ._impl_common import (
    bad_input,
    bad_option,
    float_list,
    get_choice,
    get_float,
    get_int,
    get_str,
    load_rows,
    output_format,
    render_rows,
    to_float,
    write_or_return,
)

_TRUE_NAMES = ("y_true", "label", "labels", "target")
_PRED_NAMES = ("y_pred", "prediction", "predictions")
_SCORE_NAMES = ("y_score", "y_prob", "score", "probability")


def _pick(
    columns: Mapping[str, Sequence[object]], explicit: object, names: Sequence[str]
) -> list[object] | None:
    keys = [str(explicit)] if explicit is not None else list(names)
    for key in keys:
        if key in columns:
            return list(columns[key])
    if explicit is not None:
        raise bad_input(f"column '{explicit}' was not found.", code="ML_COLUMN_NOT_FOUND")
    return None


def prediction_columns(
    value: object, options: Mapping[str, object]
) -> dict[str, list[object] | None]:
    """Extract ``y_true`` / ``y_pred`` / ``y_score`` columns from a mapping, records or a table file.

    Column names come from ``true_column`` / ``pred_column`` / ``score_column`` (defaults:
    ``y_true|label|target``, ``y_pred|prediction``, ``y_score|y_prob|score|probability``).
    """

    columns: Mapping[str, Sequence[object]]
    if (
        isinstance(value, Mapping)
        and any(isinstance(v, Sequence) and not isinstance(v, (str, bytes)) for v in value.values())
        and not any(k in value for k in ("records", "rows", "data"))
    ):
        columns = {
            str(k): v
            for k, v in value.items()
            if isinstance(v, Sequence) and not isinstance(v, (str, bytes))
        }
    else:
        rows = load_rows(value, what="predictions")
        names = list(dict.fromkeys(k for row in rows for k in row))
        columns = {name: [row.get(name) for row in rows] for name in names}
    y_true = _pick(columns, options.get("true_column"), _TRUE_NAMES)
    if y_true is None:
        raise bad_input(
            "input needs a y_true column.",
            code="ML_COLUMN_NOT_FOUND",
            hint="Pass {'y_true': [...], 'y_pred': [...]} or set options['true_column'].",
        )
    y_pred = _pick(columns, options.get("pred_column"), _PRED_NAMES)
    y_score = _pick(columns, options.get("score_column"), _SCORE_NAMES)
    if not y_true:
        raise bad_input("y_true is empty.")
    for name, column in (("y_pred", y_pred), ("y_score", y_score)):
        if column is not None and len(column) != len(y_true):
            raise bad_input(
                f"{name} has {len(column)} values but y_true has {len(y_true)}.",
                code="ML_LENGTH_MISMATCH",
            )
    for name, column in (("y_true", y_true), ("y_pred", y_pred), ("y_score", y_score)):
        if column is not None and any(
            v is None or (isinstance(v, str) and not v.strip()) for v in column
        ):
            raise bad_input(f"{name} contains missing values.", code="ML_MISSING_VALUES")
    return {"y_true": y_true, "y_pred": y_pred, "y_score": y_score}


def _numeric(value: object) -> float | None:
    if isinstance(value, bool):
        return float(value)
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _task(options: Mapping[str, object], y_true: Sequence[object]) -> str:
    task = get_choice(options, "task", "auto", ("auto", "classification", "regression"))
    if task != "auto":
        return task
    numbers = [_numeric(v) for v in y_true]
    if all(n is not None for n in numbers) and any(n != int(n) for n in numbers if n is not None):
        return "regression"
    return "classification"


def _binary_setup(y_true: Sequence[str], options: Mapping[str, object]) -> tuple[list[str], str]:
    classes = M.sorted_labels(y_true)
    given = options.get("positive_label")
    if given is not None:
        positive = M.label_text(given)
        if positive not in classes:
            classes = M.sorted_labels([*classes, positive])
    elif len(classes) == 2:
        positive = classes[1]
    elif len(classes) == 1 and classes[0] in {"0", "1"}:
        positive = "1"
        classes = ["0", "1"]
    else:
        raise bad_input(
            "binary scores need exactly two classes or options['positive_label'].",
            code="ML_AMBIGUOUS_POSITIVE",
        )
    if len(classes) != 2:
        raise bad_input(
            "score-based evaluation supports two classes only.", code="ML_MULTICLASS_SCORES"
        )
    return classes, positive


def evaluate(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Compute classification or regression metrics from predictions."""

    cols = prediction_columns(value, options)
    y_true_raw, y_pred_raw, y_score_raw = cols["y_true"], cols["y_pred"], cols["y_score"]
    assert y_true_raw is not None
    task = _task(options, y_true_raw)
    wanted = options.get("metrics")
    if task == "regression":
        if y_pred_raw is None:
            raise bad_input(
                "regression evaluation needs a y_pred column.", code="ML_COLUMN_NOT_FOUND"
            )
        t = [to_float(v, what="y_true") for v in y_true_raw]
        p = [to_float(v, what="y_pred") for v in y_pred_raw]
        metrics = M.regression_metrics(t, p)
        extras: dict[str, object] = {}
    else:
        y_true = [M.label_text(v) for v in y_true_raw]
        scores = None
        positive = None
        if y_score_raw is not None:
            if any(isinstance(v, (list, tuple)) for v in y_score_raw):
                raise bad_input(
                    "evaluate accepts binary scores (one number per row) only.",
                    code="ML_MULTICLASS_SCORES",
                )
            scores = [to_float(v, what="y_score") for v in y_score_raw]
        if y_pred_raw is not None:
            y_pred = [M.label_text(v) for v in y_pred_raw]
            classes = options.get("classes")
            if classes is None:
                class_list = M.sorted_labels([*y_true, *y_pred])
            else:
                if not isinstance(classes, Sequence) or isinstance(classes, str):
                    raise bad_option("option 'classes' must be a list of labels.")
                class_list = [M.label_text(c) for c in classes]
                unknown = (set(y_true) | set(y_pred)) - set(class_list)
                if unknown:
                    raise bad_input(
                        f"labels not in options['classes']: {sorted(unknown)}.",
                        code="ML_UNKNOWN_CLASS",
                    )
            if scores is not None:
                class_list, positive = _binary_setup([*y_true, *y_pred], options)
        elif scores is not None:
            threshold = get_float(options, "threshold", 0.5)
            class_list, positive = _binary_setup(y_true, options)
            negative = class_list[0] if positive == class_list[1] else class_list[1]
            y_pred = [positive if s >= threshold else negative for s in scores]
        else:
            raise bad_input(
                "classification evaluation needs y_pred or y_score.", code="ML_COLUMN_NOT_FOUND"
            )
        metrics = M.classification_metrics(
            y_true, y_pred, class_list, scores=scores, positive_label=positive
        )
        extras = {
            "per_class": metrics.pop("per_class"),
            "confusion_matrix": metrics.pop("confusion_matrix"),
        }
    if wanted is not None:
        if not isinstance(wanted, Sequence) or isinstance(wanted, str):
            raise bad_option("option 'metrics' must be a list of metric names.")
        unknown = [m for m in wanted if m not in metrics]
        if unknown:
            raise bad_option(
                f"unknown metrics {unknown}; available: {sorted(metrics)}.",
                code="UNKNOWN_ML_METRIC",
            )
        metrics = {m: metrics[m] for m in wanted}
    return {"task": task, "n": len(y_true_raw), "metrics": metrics, **extras}


def check_predictions(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Check prediction columns for length, missing values, finiteness and probability range."""

    errors: list[dict[str, str]] = []
    warnings: list[str] = []
    try:
        cols = prediction_columns(value, options)
    except Exception as exc:
        if not hasattr(exc, "code") or exc.code not in {
            "ML_LENGTH_MISMATCH",
            "ML_MISSING_VALUES",
        }:
            raise
        return {
            "valid": False,
            "n": None,
            "errors": [{"code": exc.code, "message": getattr(exc, "message", str(exc))}],
            "warnings": warnings,
        }
    y_true = cols["y_true"]
    assert y_true is not None
    for name in ("y_pred", "y_score"):
        column = cols[name]
        if column is None:
            continue
        for index, v in enumerate(column):
            if name == "y_score" and isinstance(v, (list, tuple)):
                continue
            if name == "y_score" and _numeric(v) is None:
                errors.append(
                    {
                        "code": "NON_FINITE_SCORE",
                        "message": f"y_score[{index}]={v!r} is not a finite number",
                    }
                )
                break
    scores = cols["y_score"]
    if scores is not None and not errors:
        numbers = [_numeric(v) for v in scores if not isinstance(v, (list, tuple))]
        if numbers and any(n is not None and not 0.0 <= n <= 1.0 for n in numbers):
            warnings.append(
                "y_score has values outside [0, 1]; they are treated as raw scores, not probabilities."
            )
    task = _task(options, y_true)
    if task == "regression":
        for name in ("y_true", "y_pred"):
            column = cols[name]
            if column is not None and any(_numeric(v) is None for v in column):
                errors.append(
                    {"code": "NON_NUMERIC", "message": f"{name} must be numeric for regression"}
                )
        if cols["y_pred"] is None:
            errors.append({"code": "MISSING_COLUMN", "message": "regression needs y_pred"})
    else:
        if cols["y_pred"] is None and scores is None:
            errors.append(
                {"code": "MISSING_COLUMN", "message": "classification needs y_pred or y_score"}
            )
        if cols["y_pred"] is not None:
            unseen = {M.label_text(v) for v in cols["y_pred"]} - {M.label_text(v) for v in y_true}
            if unseen:
                warnings.append(f"y_pred contains labels never seen in y_true: {sorted(unseen)}.")
    return {
        "valid": not errors,
        "n": len(y_true),
        "task": task,
        "errors": errors,
        "warnings": warnings,
    }


def validate_dataset(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Validate a labelled sequence dataset (records or table file) for ML training."""

    rows = load_rows(value, what="dataset")
    seq_col = get_str(options, "sequence_column", "sequence") or "sequence"
    label_col = get_str(options, "label_column", "label") or "label"
    id_col = get_str(options, "id_column", "id") or "id"
    split_col = get_str(options, "split_column", None)
    alphabet_text = "".join(dict.fromkeys(str(options.get("alphabet", "ACGTN")).upper()))  # order as given
    alphabet = set(alphabet_text)
    min_len = get_int(options, "min_length", 1, low=0)
    max_len = get_int(options, "max_length", 10_000_000, low=1)
    errors: list[dict[str, object]] = []
    warnings: list[dict[str, object]] = []

    def add(bucket: list[dict[str, object]], row: int | None, code: str, message: str) -> None:
        if len(bucket) < 100:
            bucket.append({"row": row, "code": code, "message": message})

    if not rows:
        add(errors, None, "EMPTY_DATASET", "dataset has no records")
    ids: dict[str, int] = {}
    by_sequence: dict[str, list[tuple[int, str, str | None]]] = {}
    labels: dict[str, int] = {}
    lengths: list[int] = []
    for index, row in enumerate(rows):
        seq = row.get(seq_col)
        if not isinstance(seq, str) or not seq.strip():
            add(errors, index, "MISSING_SEQUENCE", f"column '{seq_col}' is empty")
            continue
        seq = "".join(seq.split()).upper()
        lengths.append(len(seq))
        bad = sorted(set(seq) - alphabet)
        if bad:
            add(
                errors,
                index,
                "INVALID_CHARACTERS",
                f"characters {bad} are outside the alphabet {alphabet_text}",
            )
        if not min_len <= len(seq) <= max_len:
            add(
                errors,
                index,
                "LENGTH_OUT_OF_RANGE",
                f"length {len(seq)} outside [{min_len}, {max_len}]",
            )
        label = row.get(label_col)
        if label is None or (isinstance(label, str) and not label.strip()):
            add(errors, index, "MISSING_LABEL", f"column '{label_col}' is empty")
            label_text = None
        else:
            label_text = M.label_text(label)
            labels[label_text] = labels.get(label_text, 0) + 1
        if id_col in row:
            key = str(row[id_col])
            if key in ids:
                add(errors, index, "DUPLICATE_ID", f"id '{key}' already used by row {ids[key]}")
            ids[key] = index
        split = None
        if split_col is not None:
            split = str(row.get(split_col, "")).strip() or None
            if split is None:
                add(errors, index, "MISSING_SPLIT", f"column '{split_col}' is empty")
        by_sequence.setdefault(seq, []).append((index, label_text or "", split))
    duplicate_rows = 0
    for members in by_sequence.values():
        if len(members) < 2:
            continue
        duplicate_rows += len(members) - 1
        if len({m[1] for m in members}) > 1:
            add(
                errors,
                members[0][0],
                "CONFLICTING_LABELS",
                f"identical sequence appears in rows {[m[0] for m in members][:5]} with different labels",
            )
        if split_col is not None and len({m[2] for m in members}) > 1:
            add(
                errors,
                members[0][0],
                "SPLIT_LEAKAGE",
                f"identical sequence appears in several splits (rows {[m[0] for m in members][:5]})",
            )
    if duplicate_rows:
        add(
            warnings,
            None,
            "DUPLICATE_SEQUENCES",
            f"{duplicate_rows} rows repeat an earlier sequence",
        )
    if len(labels) == 1:
        add(warnings, None, "SINGLE_CLASS", "only one distinct label is present")
    summary = {
        "n_records": len(rows),
        "label_counts": dict(sorted(labels.items())),
        "length_min": min(lengths) if lengths else None,
        "length_max": max(lengths) if lengths else None,
        "duplicate_sequence_rows": duplicate_rows,
    }
    return {"valid": not errors, "errors": errors, "warnings": warnings, "summary": summary}


# ---------------------------------------------------------------- calibration


def fit_calibrator(
    value: object, options: Mapping[str, object]
) -> tuple[dict[str, object], list[int], list[float]]:
    """Fit a Platt or isotonic calibrator on (y_true, y_score); returns model, binary labels, raw scores."""

    cols = prediction_columns(value, options)
    if cols["y_score"] is None:
        raise bad_input("calibration needs a y_score column.", code="ML_COLUMN_NOT_FOUND")
    y_true = [M.label_text(v) for v in cols["y_true"] or []]
    classes, positive = _binary_setup(y_true, options)
    scores = [to_float(v, what="y_score") for v in cols["y_score"]]
    labels = [1 if t == positive else 0 for t in y_true]
    method = get_choice(options, "method", "platt", ("platt", "isotonic"))
    model: dict[str, object] = {
        "type": "dnakit.calibration",
        "version": 1,
        "method": method,
        "classes": classes,
        "positive_label": positive,
    }
    if method == "platt":
        a, b = M.platt_fit(labels, scores)
        model["parameters"] = {"A": a, "B": b, "formula": "p = 1 / (1 + exp(A * score + B))"}
    else:
        if len(set(labels)) < 2:
            raise bad_input(
                "isotonic calibration needs both classes in y_true.", code="CALIBRATION_ONE_CLASS"
            )
        blocks = M.pav_isotonic(scores, [float(x) for x in labels])
        model["parameters"] = {
            "blocks": [{"x_first": f, "x_last": last, "value": v} for f, last, v in blocks]
        }
    return model, labels, scores


def apply_calibrator(model: Mapping[str, object], scores: Sequence[float]) -> list[float]:
    """Apply a calibrator produced by :func:`fit_calibrator`."""

    params = model["parameters"]
    assert isinstance(params, Mapping)
    if model["method"] == "platt":
        return [M.platt_predict(float(params["A"]), float(params["B"]), s) for s in scores]
    blocks = [
        (float(b["x_first"]), float(b["x_last"]), float(b["value"])) for b in params["blocks"]
    ]
    return [M.isotonic_predict(blocks, s) for s in scores]


def calibrate(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Fit a probability calibrator and apply it."""

    model, labels, scores = fit_calibrator(value, options)
    bins = get_int(options, "n_bins", 10, low=1, high=1000)
    target = scores
    if options.get("apply_to") is not None:
        target = float_list(options["apply_to"], what="apply_to")
    calibrated = apply_calibrator(model, target)
    in_range = all(0.0 <= s <= 1.0 for s in scores)
    fitted = apply_calibrator(model, scores)
    return {
        "method": model["method"],
        "model": model,
        "n_fit": len(scores),
        "calibrated": calibrated,
        "ece_before": M.expected_calibration_error(labels, scores, bins) if in_range else None,
        "ece_after_on_fit_data": M.expected_calibration_error(labels, fitted, bins),
        "note": "ECE is measured on the data used to fit the calibrator (resubstitution); use held-out data for an unbiased estimate.",
    }


def export_calibration(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Fit a calibrator and write it as JSON."""

    model, _, _ = fit_calibrator(value, options)
    info = write_or_return(json.dumps(model, indent=2) + "\n", options)
    return {"model": model, **info}


# ---------------------------------------------------------------- uncertainty


def _depth(value: object) -> int:
    depth = 0
    while isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and value:
        value = value[0]
        depth += 1
    return depth


def uncertainty(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Predictive uncertainty of classification probabilities or regression ensembles."""

    if not isinstance(value, Mapping):
        raise bad_input("input must be a mapping with 'probabilities' or 'predictions'.")
    base = get_float(options, "log_base", math.e, low=1.0, inclusive_low=False)
    ddof = get_int(options, "ddof", 0, low=0, high=1)
    if "probabilities" in value:
        probs = value["probabilities"]
        depth = _depth(probs)
        if depth == 2:
            rows = [
                M.check_distribution(r, where=f"probabilities[{i}]") for i, r in enumerate(probs)
            ]  # type: ignore[arg-type]
            out = []
            for r in rows:
                ordered = sorted(r, reverse=True)
                out.append(
                    {
                        "predicted_class": r.index(ordered[0]),
                        "confidence": ordered[0],
                        "entropy": M.entropy(r, base),
                        "least_confidence": 1.0 - ordered[0],
                        "margin": ordered[0] - (ordered[1] if len(ordered) > 1 else 0.0),
                    }
                )
            return {"kind": "classification", "n": len(out), "samples": out}
        if depth == 3:
            members = probs  # members x samples x classes
            n_members = len(members)  # type: ignore[arg-type]
            n = len(members[0])  # type: ignore[index]
            k = len(members[0][0])  # type: ignore[index]
            table = [
                [
                    M.check_distribution(members[m][i], where=f"probabilities[{m}][{i}]")
                    for i in range(n)
                ]  # type: ignore[index]
                for m in range(n_members)
            ]
            if any(
                len(table[m]) != n or any(len(r) != k for r in table[m]) for m in range(n_members)
            ):
                raise bad_input(
                    "all ensemble members must have the same number of samples and classes.",
                    code="ML_LENGTH_MISMATCH",
                )
            out = []
            for i in range(n):
                mean = [
                    math.fsum(table[m][i][c] for m in range(n_members)) / n_members
                    for c in range(k)
                ]
                total = M.entropy(mean, base)
                aleatoric = (
                    math.fsum(M.entropy(table[m][i], base) for m in range(n_members)) / n_members
                )
                out.append(
                    {
                        "mean_probabilities": mean,
                        "predicted_class": mean.index(max(mean)),
                        "predictive_entropy": total,
                        "expected_entropy": aleatoric,
                        "mutual_information": max(0.0, total - aleatoric),
                    }
                )
            return {
                "kind": "classification_ensemble",
                "n_members": n_members,
                "n": n,
                "samples": out,
            }
        raise bad_input(
            "'probabilities' must be samples x classes, or members x samples x classes."
        )
    if "predictions" in value:
        preds = value["predictions"]
        if _depth(preds) != 2:
            raise bad_input("'predictions' must be members x samples (at least 2 members).")
        table = [float_list(row, what="predictions") for row in preds]  # type: ignore[union-attr]
        if len(table) < 1 + ddof or len({len(r) for r in table}) != 1:
            raise bad_input(
                "ensemble members must have equal lengths and at least 2 members when ddof=1.",
                code="ML_LENGTH_MISMATCH",
            )
        m = len(table)
        out = []
        for i in range(len(table[0])):
            column = [table[j][i] for j in range(m)]
            mean = math.fsum(column) / m
            var = math.fsum((x - mean) ** 2 for x in column) / (m - ddof)
            out.append(
                {
                    "mean": mean,
                    "variance": var,
                    "std": math.sqrt(var),
                    "min": min(column),
                    "max": max(column),
                }
            )
        return {
            "kind": "regression_ensemble",
            "n_members": m,
            "n": len(out),
            "ddof": ddof,
            "samples": out,
        }
    raise bad_input(
        "input mapping needs 'probabilities' or 'predictions'.", code="ML_COLUMN_NOT_FOUND"
    )


def conformal_model(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Fit a split-conformal uncertainty model from a held-out calibration set."""

    if not isinstance(value, Mapping) or "y_true" not in value:
        raise bad_input(
            "input must be a mapping with y_true and y_pred (regression) or probabilities (classification)."
        )
    alpha = get_float(
        options, "alpha", 0.1, low=0.0, high=1.0, inclusive_low=False, inclusive_high=False
    )
    if "probabilities" in value:
        probs = value["probabilities"]
        labels = [M.label_text(v) for v in value["y_true"]]  # type: ignore[union-attr]
        classes = [M.label_text(c) for c in value.get("classes", M.sorted_labels(labels))]  # type: ignore[arg-type]
        if len(labels) != len(probs):  # type: ignore[arg-type]
            raise bad_input("y_true and probabilities differ in length.", code="ML_LENGTH_MISMATCH")
        scores = []
        for i, (label, row) in enumerate(zip(labels, probs, strict=False)):  # type: ignore[arg-type]
            dist = M.check_distribution(row, where=f"probabilities[{i}]")
            if len(dist) != len(classes) or label not in classes:
                raise bad_input(
                    "probabilities must have one column per class and y_true must use known classes.",
                    code="ML_UNKNOWN_CLASS",
                )
            scores.append(1.0 - dist[classes.index(label)])
        task = "classification"
    else:
        if "y_pred" not in value:
            raise bad_input(
                "regression conformal fitting needs y_pred.", code="ML_COLUMN_NOT_FOUND"
            )
        t = float_list(value["y_true"], what="y_true")
        p = float_list(value["y_pred"], what="y_pred")
        if len(t) != len(p):
            raise bad_input("y_true and y_pred differ in length.", code="ML_LENGTH_MISMATCH")
        scores = [abs(a - b) for a, b in zip(t, p, strict=False)]
        classes = []
        task = "regression"
    if not scores:
        raise bad_input("the calibration set is empty.")
    k = M.quantile_index(len(scores), alpha)
    ordered = sorted(scores)
    model: dict[str, object] = {
        "type": "dnakit.split_conformal",
        "version": 1,
        "task": task,
        "alpha": alpha,
        "n_calibration": len(scores),
        "rank": k,
        "quantile": ordered[k - 1] if k is not None else None,
        "nonconformity": "absolute residual |y - y_hat|"
        if task == "regression"
        else "1 - p(true class)",
        "coverage_guarantee": f">= {1 - alpha:.6g} marginal coverage under exchangeability"
        if k is not None
        else "calibration set too small for this alpha: the prediction set is unbounded",
    }
    if classes:
        model["classes"] = classes
    return model


def export_uncertainty(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Fit a conformal uncertainty model and write it as JSON."""

    model = conformal_model(value, options)
    info = write_or_return(json.dumps(model, indent=2) + "\n", options)
    return {"model": model, **info}


# ---------------------------------------------------------------- tables

_LOWER_IS_BETTER = ("loss", "error", "rmse", "mae", "mse", "perplexity", "nll", "ece")


def _higher_better(options: Mapping[str, object], metric: str) -> bool:
    if "higher_is_better" in options:
        value = options["higher_is_better"]
        if not isinstance(value, bool):
            raise bad_option("option 'higher_is_better' must be a boolean.")
        return value
    return not any(token in metric.lower() for token in _LOWER_IS_BETTER)


def competition_ranks(values: Sequence[float], higher_is_better: bool) -> list[int]:
    """Standard competition ranking (1, 2, 2, 4): ties share the best rank."""

    order = sorted(values, reverse=higher_is_better)
    first: dict[float, int] = {}
    for position, v in enumerate(order, 1):
        first.setdefault(v, position)
    return [first[v] for v in values]


def leaderboard(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Rank models by a metric and render the table."""

    rows = load_rows(value, what="leaderboard")
    if not rows:
        raise bad_input("leaderboard input has no rows.")
    metric = get_str(options, "metric", "score") or "score"
    model_col = get_str(options, "model_column", "model") or "model"
    if any(metric not in r or model_col not in r for r in rows):
        raise bad_input(
            f"every row needs '{model_col}' and '{metric}'.", code="ML_COLUMN_NOT_FOUND"
        )
    values = [to_float(r[metric], what=metric) for r in rows]
    higher = _higher_better(options, metric)
    ranks = competition_ranks(values, higher)
    order = sorted(range(len(rows)), key=lambda i: (ranks[i], i))
    ranked = [{"rank": ranks[i], **rows[i]} for i in order]
    fmt = output_format(options, "tsv", ("tsv", "csv", "json", "jsonl"))
    info = write_or_return(render_rows(ranked, fmt), options)
    return {
        "metric": metric,
        "higher_is_better": higher,
        "format": fmt,
        "n": len(ranked),
        "ranking": ranked,
        **info,
    }


def benchmark(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Aggregate model-by-dataset results into per-dataset ranks and a model summary."""

    rows = load_rows(value, what="benchmark")
    model_col = get_str(options, "model_column", "model") or "model"
    data_col = get_str(options, "dataset_column", "dataset") or "dataset"
    value_col = get_str(options, "value_column", "value") or "value"
    if not rows or any(c not in r for r in rows for c in (model_col, data_col, value_col)):
        raise bad_input(
            f"every row needs '{model_col}', '{data_col}' and '{value_col}'.",
            code="ML_COLUMN_NOT_FOUND",
        )
    higher = _higher_better(options, str(options.get("metric", value_col)))
    seen: set[tuple[str, str]] = set()
    by_dataset: dict[str, list[tuple[str, float]]] = {}
    for r in rows:
        key = (str(r[data_col]), str(r[model_col]))
        if key in seen:
            raise bad_input(
                f"duplicate result for dataset '{key[0]}' and model '{key[1]}'.",
                code="ML_DUPLICATE_RESULT",
            )
        seen.add(key)
        by_dataset.setdefault(key[0], []).append((key[1], to_float(r[value_col], what=value_col)))
    long_rows: list[dict[str, object]] = []
    per_model: dict[str, dict[str, list[float]]] = {}
    for dataset, items in by_dataset.items():
        ranks = competition_ranks([v for _, v in items], higher)
        for (model, v), rank in zip(items, ranks, strict=False):
            long_rows.append({"dataset": dataset, "model": model, "value": v, "rank": rank})
            entry = per_model.setdefault(model, {"values": [], "ranks": [], "wins": []})
            entry["values"].append(v)
            entry["ranks"].append(float(rank))
            entry["wins"].append(1.0 if rank == 1 else 0.0)
    summary = [
        {
            "model": m,
            "n_datasets": len(e["values"]),
            "mean_value": math.fsum(e["values"]) / len(e["values"]),
            "mean_rank": math.fsum(e["ranks"]) / len(e["ranks"]),
            "wins": int(sum(e["wins"])),
        }
        for m, e in per_model.items()
    ]
    summary.sort(key=lambda r: (r["mean_rank"], str(r["model"])))
    fmt = output_format(options, "json", ("json", "tsv", "csv"))
    if fmt == "json":
        text = json.dumps({"per_dataset": long_rows, "summary": summary}, indent=2) + "\n"
    else:
        text = render_rows(long_rows, fmt)
    info = write_or_return(text, options)
    return {
        "higher_is_better": higher,
        "format": fmt,
        "per_dataset": long_rows,
        "summary": summary,
        **info,
    }


# ---------------------------------------------------------------- dataset export and splits

_SPLIT_ALIASES = {
    "train": "train",
    "training": "train",
    "validation": "validation",
    "valid": "validation",
    "val": "validation",
    "dev": "validation",
    "test": "test",
    "testing": "test",
}
SPLITS = ("train", "validation", "test")


def split_counts(n: int, ratios: Sequence[float]) -> list[int]:
    """Largest-remainder allocation of ``n`` rows to the splits (ties go to the earlier split)."""

    exact = [n * r for r in ratios]
    counts = [math.floor(x + 1e-9) for x in exact]
    for index in sorted(range(len(ratios)), key=lambda i: (-(exact[i] - counts[i]), i))[
        : n - sum(counts)
    ]:
        counts[index] += 1
    return counts


def assign_splits(
    rows: Sequence[Mapping[str, object]], options: Mapping[str, object]
) -> dict[str, list[int]]:
    """Assign every row to train/validation/test.

    An existing ``split_column`` (default ``split``) wins when every row has a value;
    otherwise rows (or whole ``group_column`` groups) are shuffled with ``random.Random(seed)``.
    """

    split_col = get_str(options, "split_column", "split") or "split"
    present = [split_col in r and r[split_col] not in (None, "") for r in rows]
    result: dict[str, list[int]] = {name: [] for name in SPLITS}
    if rows and all(present):
        for index, r in enumerate(rows):
            name = _SPLIT_ALIASES.get(str(r[split_col]).strip().lower())
            if name is None:
                raise bad_input(
                    f"unknown split value {r[split_col]!r} in row {index}.",
                    code="ML_UNKNOWN_SPLIT",
                    hint="Use train, validation (valid/val/dev) or test.",
                )
            result[name].append(index)
        return result
    if any(present):
        raise bad_input(
            f"column '{split_col}' is filled for some rows only.", code="ML_PARTIAL_SPLIT"
        )
    ratios = [
        get_float(options, f"{name}_ratio", default, low=0.0)
        for name, default in zip(SPLITS, (0.8, 0.1, 0.1), strict=False)
    ]
    if abs(sum(ratios) - 1.0) > 1e-6:
        raise bad_option(
            "train_ratio + validation_ratio + test_ratio must equal 1.", code="INVALID_SPLIT_RATIOS"
        )
    seed = get_int(options, "seed", 42)
    group_col = get_str(options, "group_column", None)
    units: list[list[int]]
    if group_col is None:
        units = [[i] for i in range(len(rows))]
    else:
        groups: dict[str, list[int]] = {}
        for index, r in enumerate(rows):
            if group_col not in r:
                raise bad_input(
                    f"row {index} has no '{group_col}' column.", code="ML_COLUMN_NOT_FOUND"
                )
            groups.setdefault(str(r[group_col]), []).append(index)
        units = list(groups.values())
    random.Random(seed).shuffle(units)
    targets = split_counts(len(rows), ratios)
    filled = [0, 0, 0]
    for unit in units:
        choice = max(range(3), key=lambda s: (targets[s] - filled[s], -s))
        result[SPLITS[choice]].extend(unit)
        filled[choice] += len(unit)
    for name in SPLITS:
        result[name].sort()
    return result


def export_dataset(value: object, options: Mapping[str, object]) -> dict[str, object]:
    """Write a dataset of records as csv/tsv/json/jsonl/fasta."""

    rows = load_rows(value, what="dataset")
    if not rows:
        raise bad_input("dataset has no records.")
    fmt = output_format(options, "csv", ("csv", "tsv", "json", "jsonl", "fasta"))
    columns = options.get("columns")
    if columns is not None and (not isinstance(columns, Sequence) or isinstance(columns, str)):
        raise bad_option("option 'columns' must be a list of column names.")
    info = write_or_return(render_rows(rows, fmt, columns), options)  # type: ignore[arg-type]
    return {"format": fmt, "n_records": len(rows), **info}


def export_split(value: object, options: Mapping[str, object], which: str) -> dict[str, object]:
    """Write one split (train, validation or test) of a dataset."""

    rows = load_rows(value, what="dataset")
    if not rows:
        raise bad_input("dataset has no records.")
    assignment = assign_splits(rows, options)
    selected = [rows[i] for i in assignment[which]]
    fmt = output_format(options, "csv", ("csv", "tsv", "json", "jsonl", "fasta"))
    names = list(dict.fromkeys(k for row in rows for k in row))
    info = write_or_return(render_rows(selected, fmt, options.get("columns") or names), options)  # type: ignore[arg-type]
    result: dict[str, object] = {
        "split": which,
        "format": fmt,
        "n_records": len(selected),
        "counts": {name: len(idx) for name, idx in assignment.items()},
        **info,
    }
    if len(selected) <= 10_000:
        result["row_indices"] = assignment[which]
    return result
