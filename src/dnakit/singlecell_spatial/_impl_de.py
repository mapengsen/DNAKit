"""Library-size normalisation and two-group differential expression shared by single-cell and bulk RNA functions."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from dnakit.exceptions import ConfigurationError, InputFormatError

from . import _impl_stats as stats
from ._impl_matrix import CountMatrix, opt_choice, opt_float, opt_int


def scale_columns(cm: CountMatrix, target_sum: float) -> CountMatrix:
    """Return counts scaled so that every cell (column) sums to ``target_sum`` (all-zero cells stay zero)."""

    totals = cm.column_totals()
    factors = [target_sum / total if total > 0 else 0.0 for total in totals]
    return CountMatrix(list(cm.genes), list(cm.cells), [[x * f for x, f in zip(row, factors)] for row in cm.values])


def log1p_columns(cm: CountMatrix, target_sum: float) -> CountMatrix:
    """Library-size normalise to ``target_sum`` and apply the natural ``log1p``."""

    scaled = scale_columns(cm, target_sum)
    return CountMatrix(scaled.genes, scaled.cells, [[math.log1p(x) for x in row] for row in scaled.values])


def two_groups(cm: CountMatrix, labels: Sequence[str], options: Mapping[str, object]) -> tuple[str, str, list[int], list[int]]:
    """Resolve ``group1`` / ``group2`` (default: the two distinct labels, sorted) into column indices."""

    present = sorted(set(labels))
    group1 = options.get("group1")
    group2 = options.get("group2")
    if group1 is None and group2 is None:
        if len(present) != 2:
            raise ConfigurationError(
                f"options['group1'] and options['group2'] are required when there are {len(present)} distinct labels.",
                code="MISSING_OPTION", hint="Name the two labels to compare.",
            )
        group1, group2 = present
    elif group1 is None or group2 is None:
        raise ConfigurationError("Give both options['group1'] and options['group2'].", code="MISSING_OPTION", hint="Name both groups.")
    group1, group2 = str(group1), str(group2)
    if group1 == group2:
        raise ConfigurationError("group1 and group2 must differ.", code="INVALID_OPTION", hint="Choose two different labels.")
    first = [i for i, label in enumerate(labels) if label == group1]
    second = [i for i, label in enumerate(labels) if label == group2]
    for name, members in ((group1, first), (group2, second)):
        if not members:
            raise InputFormatError(f"No cell carries the label {name!r}.", code="UNKNOWN_GROUP", hint=f"Labels present: {present}.")
    return group1, group2, first, second


def differential_expression(
    cm: CountMatrix,
    first: Sequence[int],
    second: Sequence[int],
    *,
    method: str,
    target_sum: float | None,
    pseudocount: float,
) -> list[dict[str, object]]:
    """Per-gene comparison of two column sets.

    Values are ``norm = count / library_size * target_sum`` (``target_sum=None`` keeps raw counts) and the test is run on
    ``log2(norm + pseudocount)``. ``log2fc = log2(mean(norm, group1) + pc) - log2(mean(norm, group2) + pc)``.
    ``method='welch'`` is Welch's t-test, ``'wilcoxon'`` the Mann-Whitney U / rank-sum test (normal approximation with
    tie and continuity correction). ``q_value`` is the Benjamini-Hochberg adjustment over all genes tested.
    """

    if method == "welch" and (len(first) < 2 or len(second) < 2):
        raise InputFormatError(
            "Welch's t-test needs at least two cells/samples in each group.", code="GROUP_TOO_SMALL",
            hint="Use method='wilcoxon' or add replicates.",
        )
    scaled = scale_columns(cm, target_sum) if target_sum is not None else cm
    rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for g_index, (gene, row) in enumerate(zip(scaled.genes, scaled.values)):
        raw = cm.values[g_index]
        x1 = [row[j] for j in first]
        x2 = [row[j] for j in second]
        mean1, mean2 = stats.mean(x1), stats.mean(x2)
        y1 = [math.log2(v + pseudocount) for v in x1]
        y2 = [math.log2(v + pseudocount) for v in x2]
        if method == "welch":
            statistic, _df, p_value = stats.welch_t_test(y1, y2)
        else:
            statistic, p_value = stats.mann_whitney_u(y1, y2)
        record = {
            "gene": gene,
            "mean_group1": mean1,
            "mean_group2": mean2,
            "log2fc": math.log2(mean1 + pseudocount) - math.log2(mean2 + pseudocount),
            "pct_group1": sum(1 for j in first if raw[j] > 0) / len(first),
            "pct_group2": sum(1 for j in second if raw[j] > 0) / len(second),
            "statistic": statistic if statistic is None or math.isfinite(statistic) else None,
            "p_value": p_value,
        }
        rows.append(record)
        p_values.append(p_value)
    for record, q_value in zip(rows, stats.benjamini_hochberg(p_values)):
        record["q_value"] = q_value
    return rows


def rank_genes(rows: list[dict[str, object]], *, top_n: int | None, min_pct: float = 0.0, only_positive: bool = False) -> list[dict[str, object]]:
    kept = [
        row for row in rows
        if max(row["pct_group1"], row["pct_group2"]) >= min_pct and (not only_positive or row["log2fc"] > 0)  # type: ignore[operator]
    ]
    kept.sort(key=lambda row: (row["p_value"], -abs(row["log2fc"]), row["gene"]))  # type: ignore[operator]
    return kept if top_n is None else kept[:top_n]


def de_settings(options: Mapping[str, object], default_target_sum: float) -> dict[str, object]:
    """Common differential-expression options (documented in the public docstrings)."""

    normalization = opt_choice(options, "normalization", "cpm", ("cpm", "none"))
    return {
        "method": opt_choice(options, "method", "welch", ("welch", "wilcoxon")),
        "target_sum": opt_float(options, "target_sum", default_target_sum, 0.0, math.inf, exclusive_min=True)
        if normalization == "cpm" else None,
        "pseudocount": opt_float(options, "pseudocount", 1.0, 0.0, math.inf, exclusive_min=True),
        "top_n": opt_int(options, "top_n", 100, 1, 10**7),
        "min_pct": opt_float(options, "min_pct", 0.0, 0.0, 1.0),
        "normalization": normalization,
    }
