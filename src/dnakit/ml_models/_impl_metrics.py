"""Dependency-free evaluation metrics, calibration and uncertainty measures.

Every formula follows its textbook / scikit-learn definition. Undefined values (for
example a correlation with zero variance) are reported as ``None`` so results stay JSON-safe.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence

from ._impl_common import bad_input, bad_option, to_float


def label_text(value: object) -> str:
    """Normalise a class label so that ``1``, ``1.0`` and ``'1'`` compare equal."""

    if isinstance(value, bool):
        return str(int(value))
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def sorted_labels(labels: Sequence[str]) -> list[str]:
    """Sort labels numerically when all are numbers, otherwise lexicographically."""

    unique = set(labels)
    try:
        return sorted(unique, key=float)
    except ValueError:
        return sorted(unique)


def average_ranks(values: Sequence[float]) -> list[float]:
    """1-based ranks with ties replaced by their average rank."""

    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start
        while end + 1 < len(order) and values[order[end + 1]] == values[order[start]]:
            end += 1
        average = (start + end) / 2.0 + 1.0
        for position in range(start, end + 1):
            ranks[order[position]] = average
        start = end + 1
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    """Pearson correlation coefficient, ``None`` when either variable is constant."""

    n = len(x)
    if n < 2:
        return None
    mx, my = math.fsum(x) / n, math.fsum(y) / n
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    if sxx == 0.0 or syy == 0.0:
        return None
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y, strict=False))
    return max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy)))


def spearman(x: Sequence[float], y: Sequence[float]) -> float | None:
    """Spearman rank correlation: Pearson correlation of the average ranks."""

    return pearson(average_ranks(x), average_ranks(y))


def regression_metrics(y_true: Sequence[float], y_pred: Sequence[float]) -> dict[str, object]:
    """MAE, MSE, RMSE, R2, explained variance, Pearson and Spearman correlation."""

    n = len(y_true)
    errors = [p - t for t, p in zip(y_true, y_pred, strict=False)]
    mse = math.fsum(e * e for e in errors) / n
    mean_true = math.fsum(y_true) / n
    ss_tot = math.fsum((t - mean_true) ** 2 for t in y_true)
    ss_res = math.fsum(e * e for e in errors)
    # constant targets follow scikit-learn's convention: 1.0 for a perfect fit, else 0.0
    r2 = (1.0 if ss_res == 0.0 else 0.0) if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    mean_error = math.fsum(errors) / n
    var_error = math.fsum((e - mean_error) ** 2 for e in errors) / n
    var_true = ss_tot / n
    if var_true == 0.0:
        explained = 1.0 if var_error == 0.0 else 0.0
    else:
        explained = 1.0 - var_error / var_true
    return {
        "n": n,
        "mae": math.fsum(abs(e) for e in errors) / n,
        "mse": mse,
        "rmse": math.sqrt(mse),
        "r2": r2,
        "explained_variance": explained,
        "pearson": pearson(y_true, y_pred),
        "spearman": spearman(y_true, y_pred),
    }


def auroc(positive: Sequence[int], scores: Sequence[float]) -> float | None:
    """Area under the ROC curve from the Mann-Whitney U statistic (ties count one half)."""

    n_pos = sum(positive)
    n_neg = len(positive) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    ranks = average_ranks(scores)
    rank_sum = math.fsum(r for r, p in zip(ranks, positive, strict=False) if p)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def average_precision(positive: Sequence[int], scores: Sequence[float]) -> float | None:
    """Area under the precision-recall curve as step-wise average precision.

    ``AP = sum_k (R_k - R_{k-1}) * P_k`` over the distinct score thresholds in descending order.
    """

    n_pos = sum(positive)
    if n_pos == 0:
        return None
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = fp = 0
    previous_recall = 0.0
    total = 0.0
    index = 0
    while index < len(order):
        threshold = scores[order[index]]
        while index < len(order) and scores[order[index]] == threshold:
            if positive[order[index]]:
                tp += 1
            else:
                fp += 1
            index += 1
        recall = tp / n_pos
        total += (recall - previous_recall) * (tp / (tp + fp))
        previous_recall = recall
    return total


def binary_log_loss(positive: Sequence[int], probabilities: Sequence[float]) -> float:
    """Mean negative log-likelihood, probabilities clipped to [1e-15, 1 - 1e-15]."""

    eps = 1e-15
    total = 0.0
    for label, p in zip(positive, probabilities, strict=False):
        p = min(max(p, eps), 1.0 - eps)
        total -= math.log(p if label else 1.0 - p)
    return total / len(positive)


def expected_calibration_error(
    positive: Sequence[int], probabilities: Sequence[float], n_bins: int = 10
) -> float:
    """Expected calibration error with ``n_bins`` equal-width bins on [0, 1].

    ``ECE = sum_b (n_b / n) * |mean(label in b) - mean(prob in b)|``; probability 1.0 falls in the last bin.
    """

    bins: list[list[int]] = [[] for _ in range(n_bins)]
    for index, p in enumerate(probabilities):
        bins[min(int(p * n_bins), n_bins - 1)].append(index)
    n = len(positive)
    total = 0.0
    for members in bins:
        if members:
            accuracy = sum(positive[i] for i in members) / len(members)
            confidence = math.fsum(probabilities[i] for i in members) / len(members)
            total += len(members) / n * abs(accuracy - confidence)
    return total


def confusion_matrix(
    y_true: Sequence[str], y_pred: Sequence[str], classes: Sequence[str]
) -> list[list[int]]:
    """Rows are true classes, columns predicted classes (scikit-learn orientation)."""

    index = {name: i for i, name in enumerate(classes)}
    matrix = [[0] * len(classes) for _ in classes]
    for t, p in zip(y_true, y_pred, strict=False):
        matrix[index[t]][index[p]] += 1
    return matrix


def classification_metrics(
    y_true: Sequence[str],
    y_pred: Sequence[str],
    classes: Sequence[str],
    *,
    scores: Sequence[float] | None = None,
    positive_label: str | None = None,
) -> dict[str, object]:
    """Accuracy, balanced accuracy, per-class and averaged precision/recall/F1, MCC, kappa.

    For two classes (or when ``positive_label`` is set with scores) it also reports specificity,
    AUROC, average precision (AUPRC) and log loss from ``scores`` (probability of the positive class).
    The multiclass MCC is Gorodkin's R_K coefficient; zero denominators give 0.
    """

    n = len(y_true)
    matrix = confusion_matrix(y_true, y_pred, classes)
    k = len(classes)
    row_sums = [sum(row) for row in matrix]
    col_sums = [sum(matrix[i][j] for i in range(k)) for j in range(k)]
    correct = sum(matrix[i][i] for i in range(k))
    per_class: dict[str, dict[str, object]] = {}
    for i, name in enumerate(classes):
        tp = matrix[i][i]
        precision = tp / col_sums[i] if col_sums[i] else 0.0
        recall = tp / row_sums[i] if row_sums[i] else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_class[name] = {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "support": row_sums[i],
        }
    present = [name for i, name in enumerate(classes) if row_sums[i] > 0]
    result: dict[str, object] = {
        "n": n,
        "accuracy": correct / n,
        "balanced_accuracy": math.fsum(float(per_class[c]["recall"]) for c in present)
        / len(present),  # type: ignore[arg-type]
        "precision_macro": math.fsum(float(per_class[c]["precision"]) for c in classes) / k,  # type: ignore[arg-type]
        "recall_macro": math.fsum(float(per_class[c]["recall"]) for c in classes) / k,  # type: ignore[arg-type]
        "f1_macro": math.fsum(float(per_class[c]["f1"]) for c in classes) / k,  # type: ignore[arg-type]
        "f1_weighted": math.fsum(
            float(per_class[c]["f1"]) * int(per_class[c]["support"]) for c in classes
        )
        / n,  # type: ignore[arg-type]
    }
    denominator = math.sqrt(float(n * n - sum(c * c for c in col_sums))) * math.sqrt(
        float(n * n - sum(r * r for r in row_sums))
    )
    numerator = correct * n - sum(r * c for r, c in zip(row_sums, col_sums, strict=False))
    result["mcc"] = numerator / denominator if denominator else 0.0
    expected = sum(r * c for r, c in zip(row_sums, col_sums, strict=False)) / (n * n)
    observed = correct / n
    result["cohen_kappa"] = (observed - expected) / (1 - expected) if expected != 1 else 0.0
    if k == 2 or positive_label is not None:
        pos = positive_label if positive_label is not None else classes[1]
        p = list(classes).index(pos)
        tp = matrix[p][p]
        fn = row_sums[p] - tp
        fp = col_sums[p] - tp
        tn = n - tp - fn - fp
        result.update(
            positive_label=pos,
            precision=per_class[pos]["precision"],
            recall=per_class[pos]["recall"],
            sensitivity=per_class[pos]["recall"],
            f1=per_class[pos]["f1"],
            specificity=tn / (tn + fp) if tn + fp else None,
            tp=tp,
            fp=fp,
            tn=tn,
            fn=fn,
        )
        if scores is not None:
            positive = [1 if t == pos else 0 for t in y_true]
            result["auroc"] = auroc(positive, scores)
            result["auprc"] = average_precision(positive, scores)
            if all(0.0 <= s <= 1.0 for s in scores):
                result["log_loss"] = binary_log_loss(positive, scores)
    result["per_class"] = per_class
    result["confusion_matrix"] = {"classes": list(classes), "matrix": matrix}
    return result


def pav_isotonic(x: Sequence[float], y: Sequence[float]) -> list[tuple[float, float, float]]:
    """Isotonic (non-decreasing) regression by pool-adjacent-violators.

    Returns blocks ``(x_first, x_last, value)`` over ascending distinct x; tied x values are
    first averaged. The fitted function is constant inside a block and linear between blocks.
    """

    order = sorted(range(len(x)), key=lambda i: x[i])
    groups: list[list[float]] = []  # [x, weight, mean_y]
    for i in order:
        if groups and groups[-1][0] == x[i]:
            g = groups[-1]
            g[2] = (g[2] * g[1] + y[i]) / (g[1] + 1)
            g[1] += 1
        else:
            groups.append([x[i], 1.0, float(y[i])])
    blocks: list[list[float]] = []  # [x_first, x_last, weight, value]
    for gx, gw, gy in groups:
        blocks.append([gx, gx, gw, gy])
        while len(blocks) > 1 and blocks[-2][3] >= blocks[-1][3]:
            b = blocks.pop()
            a = blocks[-1]
            total = a[2] + b[2]
            a[3] = (a[3] * a[2] + b[3] * b[2]) / total
            a[2] = total
            a[1] = b[1]
    return [(b[0], b[1], b[3]) for b in blocks]


def isotonic_predict(blocks: Sequence[tuple[float, float, float]], value: float) -> float:
    """Evaluate a fitted isotonic function with clipping outside the fitted range."""

    if value <= blocks[0][1]:
        return blocks[0][2]
    if value >= blocks[-1][0]:
        return blocks[-1][2]
    for left, right in itertools.pairwise(blocks):
        if left[0] <= value <= left[1]:
            return left[2]
        if left[1] < value < right[0]:
            span = right[0] - left[1]
            return left[2] + (right[2] - left[2]) * (value - left[1]) / span
    return blocks[-1][2]


def platt_fit(positive: Sequence[int], scores: Sequence[float]) -> tuple[float, float]:
    """Fit Platt scaling ``p = 1 / (1 + exp(A * f + B))`` (Lin, Lin & Weng 2007, Newton + backtracking).

    Targets are smoothed to ``(N+ + 1) / (N+ + 2)`` and ``1 / (N- + 2)`` so separable data stay finite.
    """

    n_pos = sum(positive)
    n_neg = len(positive) - n_pos
    if n_pos == 0 or n_neg == 0:
        raise bad_input("Platt scaling needs both classes in y_true.", code="CALIBRATION_ONE_CLASS")
    hi, lo = (n_pos + 1.0) / (n_pos + 2.0), 1.0 / (n_neg + 2.0)
    targets = [hi if p else lo for p in positive]
    a, b = 0.0, math.log((n_neg + 1.0) / (n_pos + 1.0))

    def objective(a_: float, b_: float) -> float:
        total = 0.0
        for f, t in zip(scores, targets, strict=False):
            z = f * a_ + b_
            if z >= 0:
                total += t * z + math.log1p(math.exp(-z))
            else:
                total += (t - 1.0) * z + math.log1p(math.exp(z))
        return total

    value = objective(a, b)
    for _ in range(100):
        h11 = h22 = 1e-12
        h21 = g1 = g2 = 0.0
        for f, t in zip(scores, targets, strict=False):
            z = f * a + b
            if z >= 0:
                p, q = math.exp(-z) / (1.0 + math.exp(-z)), 1.0 / (1.0 + math.exp(-z))
            else:
                p, q = 1.0 / (1.0 + math.exp(z)), math.exp(z) / (1.0 + math.exp(z))
            d2 = p * q
            h11 += f * f * d2
            h22 += d2
            h21 += f * d2
            d1 = t - p
            g1 += f * d1
            g2 += d1
        if abs(g1) < 1e-9 and abs(g2) < 1e-9:
            break
        det = h11 * h22 - h21 * h21
        da, db = -(h22 * g1 - h21 * g2) / det, -(-h21 * g1 + h11 * g2) / det
        gd = g1 * da + g2 * db
        step = 1.0
        while step >= 1e-10:
            new_a, new_b = a + step * da, b + step * db
            new_value = objective(new_a, new_b)
            if new_value < value + 1e-4 * step * gd:
                a, b, value = new_a, new_b, new_value
                break
            step /= 2.0
        else:
            break
    return a, b


def platt_predict(a: float, b: float, score: float) -> float:
    """Apply Platt scaling stably."""

    z = a * score + b
    if z >= 0:
        return math.exp(-z) / (1.0 + math.exp(-z))
    return 1.0 / (1.0 + math.exp(z))


def entropy(probabilities: Sequence[float], base: float = math.e) -> float:
    """Shannon entropy of a discrete distribution (``0 log 0 = 0``)."""

    return -math.fsum(p * math.log(p) for p in probabilities if p > 0.0) / math.log(base)


def check_distribution(row: Sequence[float], *, where: str) -> list[float]:
    """Validate that a probability row is non-negative and sums to 1 (tolerance 1e-6)."""

    values = [to_float(v, what=where) for v in row]
    if any(v < 0.0 for v in values) or abs(math.fsum(values) - 1.0) > 1e-6:
        raise bad_input(f"{where} must be a probability distribution (non-negative, summing to 1).")
    return values


def quantile_index(n: int, alpha: float) -> int | None:
    """Rank ``k = ceil((n + 1) * (1 - alpha))`` of the split-conformal quantile; ``None`` if ``k > n``."""

    k = math.ceil((n + 1) * (1.0 - alpha) - 1e-12)
    return k if k <= n else None


def option_mapping(options: Mapping[str, object], key: str) -> Mapping[str, object]:
    """Return a mapping-valued option or an empty mapping."""

    value = options.get(key, {})
    if not isinstance(value, Mapping):
        raise bad_option(f"option '{key}' must be a mapping.")
    return value
