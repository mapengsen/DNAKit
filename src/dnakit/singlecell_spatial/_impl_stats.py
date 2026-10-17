"""Standard-library statistics used by the single-cell, spatial and bulk RNA functions.

Every function follows the textbook definition and is checked against scipy in the unit tests.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


def mean(values: Sequence[float]) -> float:
    return math.fsum(values) / len(values)


def variance(values: Sequence[float], ddof: int = 1) -> float:
    centre = mean(values)
    return math.fsum((x - centre) ** 2 for x in values) / (len(values) - ddof)


def median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def quantile(values: Sequence[float], q: float) -> float:
    """Linear-interpolation quantile (numpy's default ``method='linear'``)."""

    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    low = math.floor(position)
    high = math.ceil(position)
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (modified Lentz)."""

    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1.0 / (tiny if abs(d) < tiny else d)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (tiny if abs(d) < tiny else d)
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (tiny if abs(d) < tiny else d)
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta function I_x(a, b)."""

    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    """Two-sided p-value of Student's t distribution with ``df`` degrees of freedom."""

    if math.isinf(t):
        return 0.0
    return betainc(df / 2.0, 0.5, df / (df + t * t))


def welch_t_test(x: Sequence[float], y: Sequence[float]) -> tuple[float | None, float | None, float]:
    """Welch's unequal-variance t-test: returns ``(t, df, two-sided p)``.

    Zero variance in both groups gives ``p = 1`` when the means are equal and ``p = 0`` (t = +-inf)
    otherwise. Both groups need at least two observations.
    """

    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        raise ValueError("Welch's t-test needs at least two observations per group")
    m1, m2 = mean(x), mean(y)
    v1, v2 = variance(x), variance(y)
    se2 = v1 / n1 + v2 / n2
    if se2 == 0.0:
        if m1 == m2:
            return 0.0, None, 1.0
        return math.copysign(math.inf, m1 - m2), None, 0.0
    t = (m1 - m2) / math.sqrt(se2)
    df = se2**2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    return t, df, t_two_sided_p(t, df)


def average_ranks(values: Sequence[float]) -> list[float]:
    """Ranks starting at 1 with ties given their average rank."""

    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        average = (start + stop) / 2.0 + 1.0
        for position in range(start, stop + 1):
            ranks[order[position]] = average
        start = stop + 1
    return ranks


def mann_whitney_u(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    """Mann-Whitney U (Wilcoxon rank-sum) test: ``(U of x, two-sided p)``.

    The p-value uses the normal approximation with tie correction and continuity correction (the
    same definition as ``scipy.stats.mannwhitneyu(..., method='asymptotic')``).
    """

    n1, n2 = len(x), len(y)
    if n1 < 1 or n2 < 1:
        raise ValueError("the rank-sum test needs at least one observation per group")
    combined = list(x) + list(y)
    ranks = average_ranks(combined)
    u1 = math.fsum(ranks[:n1]) - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    tie_term = 0.0
    counts: dict[float, int] = {}
    for value in combined:
        counts[value] = counts.get(value, 0) + 1
    for count in counts.values():
        tie_term += count**3 - count
    variance_u = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else 0.0
    if variance_u <= 0.0:
        return u1, 1.0
    centre = n1 * n2 / 2.0
    z = (abs(u1 - centre) - 0.5) / math.sqrt(variance_u)
    z = max(z, 0.0)
    return u1, math.erfc(z / math.sqrt(2.0))


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values (q-values), in the input order."""

    n = len(p_values)
    order = sorted(range(n), key=lambda index: p_values[index])
    adjusted = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * n / rank)
        adjusted[index] = running
    return adjusted


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value of the 2x2 table [[a, b], [c, d]] (sum of tables as or less likely)."""

    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2
    denominator = math.comb(total, col1)

    def probability(k: int) -> float:
        return math.comb(row1, k) * math.comb(row2, col1 - k) / denominator

    low = max(0, col1 - row2)
    high = min(row1, col1)
    observed = probability(a)
    tolerance = 1.0 + 1e-7
    return min(1.0, math.fsum(p for p in (probability(k) for k in range(low, high + 1)) if p <= observed * tolerance))


def pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    mx, my = mean(x), mean(y)
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    if sxx == 0.0 or syy == 0.0:
        return None
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy)


def spearman(x: Sequence[float], y: Sequence[float]) -> float | None:
    return pearson(average_ranks(x), average_ranks(y))
