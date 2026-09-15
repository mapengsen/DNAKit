"""Target encoding, decoding, and comprehensive supervised-learning metrics."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, cast

from dnakit._dependency_hints import dependency_hint
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError

from .schema import ClassificationMode, TargetSchema, TaskType


def _numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - provided by the training extra.
        raise BackendUnavailableError(
            "Custom property training requires NumPy.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("ml"),
        ) from exc
    return np


def _sklearn_metrics() -> Any:
    try:
        from sklearn import metrics  # type: ignore[import-not-found]
    except ImportError as exc:
        raise BackendUnavailableError(
            "Custom property metrics require scikit-learn.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint=dependency_hint("ml"),
        ) from exc
    return metrics


def _class_text(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise InputFormatError(
                "Classification labels must be finite scalar values.",
                code="INVALID_CUSTOM_LABEL",
            )
        if value.is_integer():
            return str(int(value))
    if isinstance(value, (str, int, float)):
        text = str(value).strip()
        if text:
            return text
    raise InputFormatError(
        "Classification labels must be non-empty scalar values.",
        code="INVALID_CUSTOM_LABEL",
        context={"type": type(value).__name__},
    )


def _sorted_classes(values: Sequence[str]) -> tuple[str, ...]:
    def order(value: str) -> tuple[int, float | str]:
        try:
            return 0, float(value)
        except ValueError:
            return 1, value

    return tuple(sorted(set(values), key=order))


def _binary_value(value: object) -> int | None:
    text = _class_text(value).lower()
    if text in {"0", "false", "no", "negative", "neg"}:
        return 0
    if text in {"1", "true", "yes", "positive", "pos"}:
        return 1
    return None


def _classification_mode(
    labels: Sequence[Sequence[object]],
    label_names: tuple[str, ...],
    requested: ClassificationMode,
) -> str:
    if requested != "auto":
        return requested
    if len(label_names) == 1:
        classes = {_class_text(row[0]) for row in labels}
        return "binary" if len(classes) == 2 else "multiclass"
    binary_rows = [[_binary_value(value) for value in row] for row in labels]
    if all(all(value is not None for value in row) for row in binary_rows):
        if all(sum(cast(int, value) for value in row) == 1 for row in binary_rows):
            return "multiclass"
        return "multilabel"
    return "multioutput"


def fit_target_schema(
    labels: Sequence[Sequence[object]],
    label_names: tuple[str, ...],
    *,
    task_type: TaskType,
    classification_mode: ClassificationMode,
    train_indices: Sequence[int],
    verbose: bool = True,
) -> TargetSchema:
    """Fit target normalization and class vocabularies using only training rows."""

    np = _numpy()
    if not label_names or any(len(row) != len(label_names) for row in labels):
        raise InputFormatError(
            "Target rows must align with at least one label field.",
            code="INVALID_CUSTOM_LABEL",
        )
    if not train_indices:
        raise ConfigurationError("The training split is empty.", code="INVALID_CUSTOM_SPLIT")
    if task_type == "regression":
        try:
            values = np.asarray(labels, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise InputFormatError(
                "Regression labels must all be numeric.", code="INVALID_CUSTOM_LABEL"
            ) from exc
        if not bool(np.isfinite(values).all()):
            raise InputFormatError(
                "Regression labels must all be finite.", code="INVALID_CUSTOM_LABEL"
            )
        training = values[np.asarray(train_indices, dtype=np.int64)]
        means = training.mean(axis=0)
        scales = training.std(axis=0)
        scales = np.where(scales > 1e-12, scales, 1.0)
        return TargetSchema(
            "regression",
            "regression",
            label_names,
            regression_means=tuple(float(value) for value in means),
            regression_scales=tuple(float(value) for value in scales),
        )

    mode = _classification_mode(labels, label_names, classification_mode)
    training_rows = [labels[index] for index in train_indices]
    if mode == "binary":
        if len(label_names) != 1:
            raise ConfigurationError(
                "binary classification requires exactly one label field.",
                code="INVALID_CUSTOM_CLASSIFICATION_MODE",
            )
        classes = _sorted_classes([_class_text(row[0]) for row in training_rows])
        if len(classes) != 2:
            raise InputFormatError(
                "Binary classification requires exactly two classes in the training split.",
                code="INVALID_CUSTOM_LABEL",
                context={"classes": classes},
            )
        schema = TargetSchema("classification", "binary", label_names, (classes,))
    elif mode == "multiclass":
        if len(label_names) == 1:
            classes = _sorted_classes([_class_text(row[0]) for row in training_rows])
            if len(classes) < 2:
                raise InputFormatError(
                    "Multiclass training requires at least two classes in the training split.",
                    code="INVALID_CUSTOM_LABEL",
                )
        else:
            parsed = [[_binary_value(value) for value in row] for row in labels]
            if not all(
                all(value is not None for value in row)
                and sum(cast(int, value) for value in row) == 1
                for row in parsed
            ):
                raise InputFormatError(
                    "Multiple fields in multiclass mode must be one-hot encoded.",
                    code="INVALID_CUSTOM_LABEL",
                )
            present = {
                next(index for index, value in enumerate(row) if value == 1)
                for row in (parsed[index] for index in train_indices)
            }
            if present != set(range(len(label_names))):
                raise InputFormatError(
                    "Every one-hot class must appear in the training split.",
                    code="INVALID_CUSTOM_LABEL",
                    context={"present_class_indices": tuple(sorted(present))},
                )
            classes = label_names
        schema = TargetSchema("classification", "multiclass", label_names, (classes,))
    elif mode == "multilabel":
        if len(label_names) < 2:
            raise ConfigurationError(
                "multilabel classification requires at least two label fields.",
                code="INVALID_CUSTOM_CLASSIFICATION_MODE",
            )
        if any(_binary_value(value) is None for row in labels for value in row):
            raise InputFormatError(
                "Multilabel targets must contain binary values.",
                code="INVALID_CUSTOM_LABEL",
            )
        schema = TargetSchema(
            "classification",
            "multilabel",
            label_names,
            tuple(("0", "1") for _ in label_names),
        )
    elif mode == "multioutput":
        if len(label_names) < 2:
            raise ConfigurationError(
                "multioutput classification requires at least two label fields.",
                code="INVALID_CUSTOM_CLASSIFICATION_MODE",
            )
        classes_by_target = tuple(
            _sorted_classes([_class_text(row[column]) for row in training_rows])
            for column in range(len(label_names))
        )
        if any(len(classes) < 2 for classes in classes_by_target):
            raise InputFormatError(
                "Every multioutput target needs at least two training classes.",
                code="INVALID_CUSTOM_LABEL",
            )
        schema = TargetSchema("classification", "multioutput", label_names, classes_by_target)
    else:  # pragma: no cover - configuration validates the mode.
        raise AssertionError(mode)
    encode_targets(labels, schema=schema)
    return schema


def encode_targets(
    labels: Sequence[Sequence[object]], *, schema: TargetSchema, verbose: bool = True
) -> Any:
    """Encode labels according to a fitted target schema."""

    np = _numpy()
    if schema.task_type == "regression":
        values = np.asarray(labels, dtype=np.float32)
        means = np.asarray(schema.regression_means, dtype=np.float32)
        scales = np.asarray(schema.regression_scales, dtype=np.float32)
        return (values - means) / scales
    if schema.mode in {"binary", "multiclass"}:
        classes = schema.class_names[0]
        lookup = {value: index for index, value in enumerate(classes)}
        if len(schema.label_names) == 1:
            text = [_class_text(row[0]) for row in labels]
            unknown = sorted(set(text) - set(lookup))
            if unknown:
                raise InputFormatError(
                    "Validation or test data contains classes absent from training.",
                    code="CUSTOM_DATA_UNSEEN_CLASS",
                    context={"classes": tuple(unknown)},
                )
            return np.asarray([lookup[value] for value in text], dtype=np.int64)
        parsed = [[_binary_value(value) for value in row] for row in labels]
        return np.asarray(
            [next(index for index, value in enumerate(row) if value == 1) for row in parsed],
            dtype=np.int64,
        )
    if schema.mode == "multilabel":
        return np.asarray(
            [[cast(int, _binary_value(value)) for value in row] for row in labels],
            dtype=np.float32,
        )
    encoded = np.empty((len(labels), len(schema.label_names)), dtype=np.int64)
    for column, classes in enumerate(schema.class_names):
        lookup = {value: index for index, value in enumerate(classes)}
        for row_index, row in enumerate(labels):
            value = _class_text(row[column])
            if value not in lookup:
                raise InputFormatError(
                    "Validation or test data contains classes absent from training.",
                    code="CUSTOM_DATA_UNSEEN_CLASS",
                    context={"target": schema.label_names[column], "class": value},
                )
            encoded[row_index, column] = lookup[value]
    return encoded


def _finite(value: object) -> float | int | None:
    try:
        converted = float(cast(Any, value))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(converted):
        return None
    if isinstance(value, int):
        return int(value)
    return converted


def _softmax(values: Any) -> Any:
    np = _numpy()
    shifted = values - values.max(axis=1, keepdims=True)
    exponent = np.exp(shifted)
    return exponent / exponent.sum(axis=1, keepdims=True)


def _sigmoid(values: Any) -> Any:
    np = _numpy()
    values = np.clip(values, -80, 80)
    return 1.0 / (1.0 + np.exp(-values))


def prediction_arrays(
    logits: Any, *, schema: TargetSchema, threshold: float, verbose: bool = True
) -> tuple[Any, Any]:
    """Convert raw head outputs into encoded predictions and probabilities/values."""

    np = _numpy()
    raw = np.asarray(logits, dtype=np.float64)
    if raw.ndim != 2 or raw.shape[1] != schema.output_dimension:
        raise ConfigurationError(
            "Prediction logits do not match the saved target schema.",
            code="INVALID_CUSTOM_PREDICTION_OUTPUT",
            context={
                "shape": tuple(int(value) for value in raw.shape),
                "expected_outputs": schema.output_dimension,
            },
        )
    if schema.task_type == "regression":
        means = np.asarray(schema.regression_means, dtype=np.float64)
        scales = np.asarray(schema.regression_scales, dtype=np.float64)
        values = raw * scales + means
        return values, values
    if schema.mode in {"binary", "multiclass"}:
        probabilities = _softmax(raw)
        return probabilities.argmax(axis=1), probabilities
    if schema.mode == "multilabel":
        probabilities = _sigmoid(raw)
        return (probabilities >= threshold).astype(np.int64), probabilities
    predictions = np.empty((raw.shape[0], len(schema.label_names)), dtype=np.int64)
    probability_blocks: list[Any] = []
    start = 0
    for column, classes in enumerate(schema.class_names):
        end = start + len(classes)
        block = _softmax(raw[:, start:end])
        predictions[:, column] = block.argmax(axis=1)
        probability_blocks.append(block)
        start = end
    return predictions, tuple(probability_blocks)


def _rank(values: Any) -> Any:
    np = _numpy()
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0
        start = end
    return ranks


def _correlation(first: Any, second: Any) -> float | None:
    np = _numpy()
    if len(first) < 2 or float(np.std(first)) <= 1e-15 or float(np.std(second)) <= 1e-15:
        return None
    return _finite(np.corrcoef(first, second)[0, 1])


def _regression_metrics(
    y_true: Any, y_pred: Any, label_names: tuple[str, ...]
) -> dict[str, object]:
    np = _numpy()
    metrics = _sklearn_metrics()
    true = np.asarray(y_true, dtype=np.float64)
    predicted = np.asarray(y_pred, dtype=np.float64)
    if true.ndim == 1:
        true = true[:, None]
    if predicted.ndim == 1:
        predicted = predicted[:, None]

    def one(column_true: Any, column_predicted: Any) -> dict[str, object]:
        errors = column_predicted - column_true
        nonzero = np.abs(column_true) > 1e-12
        denominator = np.abs(column_true) + np.abs(column_predicted)
        smape_mask = denominator > 1e-12
        return {
            "mae": _finite(metrics.mean_absolute_error(column_true, column_predicted)),
            "mse": _finite(metrics.mean_squared_error(column_true, column_predicted)),
            "rmse": _finite(math.sqrt(metrics.mean_squared_error(column_true, column_predicted))),
            "r2": (
                _finite(metrics.r2_score(column_true, column_predicted))
                if len(column_true) >= 2
                else None
            ),
            "explained_variance": (
                _finite(metrics.explained_variance_score(column_true, column_predicted))
                if len(column_true) >= 2
                else None
            ),
            "median_absolute_error": _finite(
                metrics.median_absolute_error(column_true, column_predicted)
            ),
            "max_error": _finite(np.max(np.abs(errors))),
            "mean_bias_error": _finite(np.mean(errors)),
            "mape": (
                _finite(np.mean(np.abs(errors[nonzero] / column_true[nonzero])))
                if bool(nonzero.any())
                else None
            ),
            "smape": (
                _finite(np.mean(2.0 * np.abs(errors[smape_mask]) / denominator[smape_mask]))
                if bool(smape_mask.any())
                else None
            ),
            "pearson_r": _correlation(column_true, column_predicted),
            "spearman_r": _correlation(_rank(column_true), _rank(column_predicted)),
            "support": len(column_true),
        }

    per_target = {
        name: one(true[:, index], predicted[:, index]) for index, name in enumerate(label_names)
    }
    aggregate: dict[str, object] = {}
    for key in next(iter(per_target.values())):
        if key == "support":
            continue
        values = [target[key] for target in per_target.values() if target[key] is not None]
        aggregate[key] = _finite(np.mean(values)) if values else None
    aggregate["support"] = int(true.shape[0])
    aggregate["target_count"] = int(true.shape[1])
    return {"aggregate": aggregate, "per_target": per_target}


def _safe_metric(function: Any, *args: Any, **kwargs: Any) -> float | None:
    try:
        return _finite(function(*args, **kwargs))
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def _single_classification_metrics(
    y_true: Any,
    y_pred: Any,
    probabilities: Any,
    classes: tuple[str, ...],
) -> dict[str, object]:
    np = _numpy()
    metrics = _sklearn_metrics()
    labels = np.arange(len(classes))
    result: dict[str, object] = {
        "accuracy": _safe_metric(metrics.accuracy_score, y_true, y_pred),
        "balanced_accuracy": _safe_metric(metrics.balanced_accuracy_score, y_true, y_pred),
        "matthews_correlation_coefficient": _safe_metric(metrics.matthews_corrcoef, y_true, y_pred),
        "cohen_kappa": _safe_metric(metrics.cohen_kappa_score, y_true, y_pred),
        "log_loss": _safe_metric(metrics.log_loss, y_true, probabilities, labels=labels),
        "support": len(y_true),
        "class_count": len(classes),
        "confusion_matrix": metrics.confusion_matrix(y_true, y_pred, labels=labels).tolist(),
    }
    for average in ("macro", "micro", "weighted"):
        result[f"precision_{average}"] = _safe_metric(
            metrics.precision_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"recall_{average}"] = _safe_metric(
            metrics.recall_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"f1_{average}"] = _safe_metric(
            metrics.f1_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"jaccard_{average}"] = _safe_metric(
            metrics.jaccard_score, y_true, y_pred, average=average, zero_division=0
        )
    if len(classes) == 2:
        result["auroc"] = _safe_metric(metrics.roc_auc_score, y_true, probabilities[:, 1])
        result["average_precision"] = _safe_metric(
            metrics.average_precision_score, y_true, probabilities[:, 1]
        )
        matrix = cast(list[list[int]], result["confusion_matrix"])
        tn, fp, fn, tp = matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]
        result["specificity"] = _finite(tn / (tn + fp)) if tn + fp else None
        result["sensitivity"] = _finite(tp / (tp + fn)) if tp + fn else None
    else:
        one_hot = np.eye(len(classes), dtype=np.int64)[np.asarray(y_true, dtype=np.int64)]
        result["auroc_macro_ovr"] = _safe_metric(
            metrics.roc_auc_score,
            one_hot,
            probabilities,
            average="macro",
            multi_class="ovr",
        )
        result["auroc_weighted_ovr"] = _safe_metric(
            metrics.roc_auc_score,
            one_hot,
            probabilities,
            average="weighted",
            multi_class="ovr",
        )
        result["average_precision_macro"] = _safe_metric(
            metrics.average_precision_score, one_hot, probabilities, average="macro"
        )
        if len(classes) > 2:
            result["top_2_accuracy"] = _safe_metric(
                metrics.top_k_accuracy_score,
                y_true,
                probabilities,
                k=min(2, len(classes)),
                labels=labels,
            )
    precision, recall, f1, support = metrics.precision_recall_fscore_support(
        y_true, y_pred, labels=labels, zero_division=0
    )
    result["per_class"] = {
        name: {
            "precision": _finite(precision[index]),
            "recall": _finite(recall[index]),
            "f1": _finite(f1[index]),
            "support": int(support[index]),
        }
        for index, name in enumerate(classes)
    }
    return result


def _multilabel_metrics(
    y_true: Any,
    y_pred: Any,
    probabilities: Any,
    label_names: tuple[str, ...],
) -> dict[str, object]:
    metrics = _sklearn_metrics()
    result: dict[str, object] = {
        "subset_accuracy": _safe_metric(metrics.accuracy_score, y_true, y_pred),
        "hamming_loss": _safe_metric(metrics.hamming_loss, y_true, y_pred),
        "label_ranking_average_precision": _safe_metric(
            metrics.label_ranking_average_precision_score, y_true, probabilities
        ),
        "coverage_error": _safe_metric(metrics.coverage_error, y_true, probabilities),
        "auroc_macro": _safe_metric(metrics.roc_auc_score, y_true, probabilities, average="macro"),
        "auroc_micro": _safe_metric(metrics.roc_auc_score, y_true, probabilities, average="micro"),
        "average_precision_macro": _safe_metric(
            metrics.average_precision_score, y_true, probabilities, average="macro"
        ),
        "average_precision_micro": _safe_metric(
            metrics.average_precision_score, y_true, probabilities, average="micro"
        ),
        "support": len(y_true),
        "label_count": len(label_names),
    }
    for average in ("macro", "micro", "weighted", "samples"):
        result[f"precision_{average}"] = _safe_metric(
            metrics.precision_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"recall_{average}"] = _safe_metric(
            metrics.recall_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"f1_{average}"] = _safe_metric(
            metrics.f1_score, y_true, y_pred, average=average, zero_division=0
        )
        result[f"jaccard_{average}"] = _safe_metric(
            metrics.jaccard_score, y_true, y_pred, average=average, zero_division=0
        )
    per_label: dict[str, object] = {}
    for index, name in enumerate(label_names):
        matrix = metrics.confusion_matrix(y_true[:, index], y_pred[:, index], labels=[0, 1])
        tn, fp, fn, tp = (int(value) for value in matrix.ravel())
        per_label[name] = {
            "precision": _safe_metric(
                metrics.precision_score,
                y_true[:, index],
                y_pred[:, index],
                zero_division=0,
            ),
            "recall": _safe_metric(
                metrics.recall_score,
                y_true[:, index],
                y_pred[:, index],
                zero_division=0,
            ),
            "f1": _safe_metric(
                metrics.f1_score, y_true[:, index], y_pred[:, index], zero_division=0
            ),
            "auroc": _safe_metric(metrics.roc_auc_score, y_true[:, index], probabilities[:, index]),
            "average_precision": _safe_metric(
                metrics.average_precision_score, y_true[:, index], probabilities[:, index]
            ),
            "specificity": _finite(tn / (tn + fp)) if tn + fp else None,
            "confusion_matrix": [[tn, fp], [fn, tp]],
            "support_positive": int(tp + fn),
        }
    result["per_label"] = per_label
    return result


def compute_metrics(
    encoded_targets: Any,
    logits: Any,
    *,
    schema: TargetSchema,
    threshold: float = 0.5,
    verbose: bool = True,
) -> dict[str, object]:
    """Compute task-appropriate aggregate and per-target metrics."""

    np = _numpy()
    predictions, scores = prediction_arrays(logits, schema=schema, threshold=threshold)
    targets = np.asarray(encoded_targets)
    if schema.task_type == "regression":
        means = np.asarray(schema.regression_means, dtype=np.float64)
        scales = np.asarray(schema.regression_scales, dtype=np.float64)
        original_targets = targets * scales + means
        return _regression_metrics(original_targets, predictions, schema.label_names)
    if schema.mode in {"binary", "multiclass"}:
        return _single_classification_metrics(
            targets,
            predictions,
            scores,
            schema.class_names[0],
        )
    if schema.mode == "multilabel":
        return _multilabel_metrics(targets, predictions, scores, schema.label_names)
    per_target: dict[str, object] = {}
    scalar_metrics: dict[str, list[float]] = {}
    for index, name in enumerate(schema.label_names):
        target_metrics = _single_classification_metrics(
            targets[:, index],
            predictions[:, index],
            scores[index],
            schema.class_names[index],
        )
        per_target[name] = target_metrics
        for key, value in target_metrics.items():
            if isinstance(value, (int, float)) and key not in {"support", "class_count"}:
                scalar_metrics.setdefault(key, []).append(float(value))
    aggregate: dict[str, object] = {
        key: _finite(np.mean(values)) for key, values in scalar_metrics.items() if values
    }
    aggregate["exact_match_accuracy"] = _finite(np.mean(np.all(targets == predictions, axis=1)))
    aggregate["support"] = int(targets.shape[0])
    aggregate["target_count"] = int(targets.shape[1])
    return {"aggregate": aggregate, "per_target": per_target}


def decoded_prediction_rows(
    record_ids: Sequence[str],
    logits: Any,
    *,
    schema: TargetSchema,
    threshold: float,
    verbose: bool = True,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Decode raw outputs into JSON-ready prediction and probability rows."""

    predictions, scores = prediction_arrays(logits, schema=schema, threshold=threshold)
    prediction_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    for row_index, record_id in enumerate(record_ids):
        predicted: dict[str, object] = {"record_id": record_id}
        probability: dict[str, object] = {"record_id": record_id}
        if schema.task_type == "regression":
            for column, name in enumerate(schema.label_names):
                predicted[name] = float(predictions[row_index, column])
        elif schema.mode in {"binary", "multiclass"}:
            classes = schema.class_names[0]
            selected = int(predictions[row_index])
            output_name = schema.label_names[0] if len(schema.label_names) == 1 else "class"
            predicted[output_name] = classes[selected]
            probability[output_name] = {
                name: float(scores[row_index, index]) for index, name in enumerate(classes)
            }
        elif schema.mode == "multilabel":
            for column, name in enumerate(schema.label_names):
                predicted[name] = int(predictions[row_index, column])
                probability[name] = float(scores[row_index, column])
        else:
            for column, name in enumerate(schema.label_names):
                classes = schema.class_names[column]
                predicted[name] = classes[int(predictions[row_index, column])]
                probability[name] = {
                    class_name: float(scores[column][row_index, class_index])
                    for class_index, class_name in enumerate(classes)
                }
        prediction_rows.append(predicted)
        if len(probability) > 1:
            probability_rows.append(probability)
    return prediction_rows, probability_rows


__all__ = [
    "compute_metrics",
    "decoded_prediction_rows",
    "encode_targets",
    "fit_target_schema",
    "prediction_arrays",
]

_register_public_verbose(globals())
