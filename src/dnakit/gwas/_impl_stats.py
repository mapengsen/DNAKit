"""Pure-standard-library numerical helpers shared by the GWAS/population/breeding/design candidates.

Everything here is a textbook algorithm (Numerical Recipes style series and continued fractions
for the incomplete gamma/beta functions, Acklam's inverse normal with one Halley step, Gaussian
elimination with partial pivoting, IRLS for logistic regression).  Each helper is checked against
scipy in ``others/tests/unit/test_gwas_population_design_candidates_real.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

MEDIAN_CHI2_1DF = 0.4549364231195724  # qchisq(0.5, df=1), used by the genomic-control lambda

_EPS = 1e-15
_TINY = 1e-300


# --------------------------------------------------------------------------- normal distribution
def normal_cdf(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def normal_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def normal_ppf(p: float) -> float:
    """Inverse standard normal CDF (Acklam's rational approximation plus one Halley refinement)."""

    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    )
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00, 3.754408661907416e00)
    low = 0.02425
    if p < low:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    elif p <= 1.0 - low:
        q = p - 0.5
        r = q * q
        x = (
            (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
            * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
        )
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    err = normal_cdf(x) - p
    u = err * math.sqrt(2.0 * math.pi) * math.exp(x * x / 2.0)
    return x - u / (1.0 + x * u / 2.0)


def two_sided_p_from_z(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def neglog10_p_from_z(z: float) -> float:
    """-log10 of the two-sided normal p-value; uses the asymptotic tail when erfc underflows."""

    p = two_sided_p_from_z(z)
    if p > 1e-300:
        return -math.log10(p)
    x = abs(z)
    # erfc(x/sqrt2) ~ sqrt(2/pi)/x * exp(-x^2/2) * (1 - 1/x^2 + 3/x^4)  for large x
    series = 1.0 - 1.0 / (x * x) + 3.0 / x**4
    return (x * x / 2.0 + math.log(x * math.sqrt(math.pi / 2.0)) - math.log(series)) / math.log(
        10.0
    )


def z_from_two_sided_p(p: float) -> float:
    """|z| such that the two-sided normal p-value equals ``p`` (p in (0, 1])."""

    if p <= 0.0:
        return math.inf
    if p >= 1.0:
        return 0.0
    if p > 1e-300:
        return -normal_ppf(p / 2.0)
    lo, hi = 0.0, 40.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if neglog10_p_from_z(mid) < -math.log10(p):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


# -------------------------------------------------------------------------- incomplete gamma / beta
def gammainc_lower_regularized(a: float, x: float) -> float:
    """Regularised lower incomplete gamma function P(a, x)."""

    if a <= 0.0:
        raise ValueError("a must be positive")
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * _EPS:
                break
        return total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    return 1.0 - _gamma_cf_upper(a, x)


def gammainc_upper_regularized(a: float, x: float) -> float:
    """Regularised upper incomplete gamma function Q(a, x)."""

    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - gammainc_lower_regularized(a, x)
    return _gamma_cf_upper(a, x)


def _gamma_cf_upper(a: float, x: float) -> float:
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _TINY:
            d = _TINY
        c = b + an / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def chi2_sf(x: float, df: float) -> float:
    """Survival function of the central chi-square distribution."""

    if x <= 0.0:
        return 1.0
    return gammainc_upper_regularized(df / 2.0, x / 2.0)


def chi2_cdf(x: float, df: float) -> float:
    if x <= 0.0:
        return 0.0
    return gammainc_lower_regularized(df / 2.0, x / 2.0)


def chi2_ppf(p: float, df: float) -> float:
    """Quantile of the chi-square distribution by bracketing and bisection."""

    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    hi = max(1.0, df)
    while chi2_cdf(hi, df) < p:
        hi *= 2.0
    lo = 0.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if chi2_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def ncx2_sf(x: float, df: float, nc: float) -> float:
    """Survival function of the non-central chi-square distribution (Poisson mixture of central
    ones).
    """

    if x <= 0.0:
        return 1.0
    if nc <= 0.0:
        return chi2_sf(x, df)
    half = nc / 2.0
    total = 0.0
    start = int(half)
    # sum Poisson-weighted terms outwards from the mode for numerical stability
    for direction in (1, -1):
        j = start if direction == 1 else start - 1
        while j >= 0:
            weight = math.exp(-half + j * math.log(half) - math.lgamma(j + 1))
            term = weight * chi2_sf(x, df + 2 * j)
            total += term
            if weight < 1e-18 and (j > half):
                break
            if weight < 1e-18 and direction == -1:
                break
            j += direction
            if j > start + 100000:
                break
    return min(1.0, total)


def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _TINY:
        d = _TINY
    d = 1.0 / d
    h = d
    for m in range(1, 10000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = _TINY if abs(d) < _TINY else d
        c = 1.0 + aa / c
        c = _TINY if abs(c) < _TINY else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = _TINY if abs(d) < _TINY else d
        c = 1.0 + aa / c
        c = _TINY if abs(c) < _TINY else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return h


def betainc_regularized(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta function I_x(a, b)."""

    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    """Two-sided p-value of Student's t statistic."""

    if df <= 0:
        raise ValueError("df must be positive")
    if math.isinf(t):
        return 0.0
    return betainc_regularized(df / 2.0, 0.5, df / (df + t * t))


def t_sf(t: float, df: float) -> float:
    """One-sided upper tail of Student's t distribution."""

    half = 0.5 * t_two_sided_p(t, df)
    return half if t >= 0 else 1.0 - half


def t_ppf(p: float, df: float) -> float:
    """Quantile of Student's t distribution (bisection on the CDF)."""

    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    if p == 0.5:
        return 0.0
    if p < 0.5:
        return -t_ppf(1.0 - p, df)
    lo, hi = 0.0, 1.0
    while t_sf(hi, df) > 1.0 - p:
        hi *= 2.0
        if hi > 1e12:
            break
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if t_sf(mid, df) > 1.0 - p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def nct_power(
    t_crit: float, df: float, delta: float, *, two_sided: bool = True, intervals: int = 4000
) -> float:
    """P(T' > t_crit) (+ P(T' < -t_crit) when two-sided) for a non-central t variable with ``df``
    and ``delta``.

    T' = (Z + delta) / sqrt(V / df), V ~ chi-square(df).  The expectation over V is evaluated by
    Simpson's rule in u = sqrt(V), so the integrand is smooth for every df >= 1.
    """

    spread = 14.0 * math.sqrt(2.0 * df)
    lo = math.sqrt(max(0.0, df - spread))
    hi = math.sqrt(df + spread + 50.0)
    n = intervals + (intervals % 2)
    h = (hi - lo) / n
    log_norm = math.log(2.0) - (df / 2.0) * math.log(2.0) - math.lgamma(df / 2.0)
    total = 0.0
    for i in range(n + 1):
        u = lo + i * h
        if u <= 0.0:
            density = 0.0 if df > 1 else math.exp(log_norm)
        else:
            density = math.exp(log_norm + (df - 1.0) * math.log(u) - u * u / 2.0)
        s = u / math.sqrt(df)
        value = normal_sf(t_crit * s - delta)
        if two_sided:
            value += normal_cdf(-t_crit * s - delta)
        weight = 1 if i in (0, n) else (4 if i % 2 else 2)
        total += weight * density * value
    return min(1.0, max(0.0, total * h / 3.0))


# --------------------------------------------------------------------------- descriptive statistics
def mean(values: Sequence[float]) -> float:
    return math.fsum(values) / len(values)


def median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    return ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2.0


def variance(values: Sequence[float], ddof: int = 1) -> float:
    m = mean(values)
    return math.fsum((v - m) ** 2 for v in values) / (len(values) - ddof)


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    mx, my = mean(x), mean(y)
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y, strict=False))
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return math.nan
    return sxy / math.sqrt(sxx * syy)


def quantile(values: Sequence[float], q: float) -> float:
    """Linear-interpolation (type 7) quantile."""

    ordered = sorted(values)
    pos = q * (len(ordered) - 1)
    lo = math.floor(pos)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] + (pos - lo) * (ordered[hi] - ordered[lo])


# --------------------------------------------------------------------------- small linear algebra
Matrix = list[list[float]]


def solve_linear(a: Matrix, b: Sequence[float]) -> list[float]:
    """Solve A x = b by Gaussian elimination with partial pivoting."""

    n = len(a)
    m = [[*list(map(float, a[i])), float(b[i])] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ArithmeticError("singular matrix")
        m[col], m[pivot] = m[pivot], m[col]
        for row in range(col + 1, n):
            factor = m[row][col] / m[col][col]
            if factor:
                for k in range(col, n + 1):
                    m[row][k] -= factor * m[col][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][k] * x[k] for k in range(i + 1, n))) / m[i][i]
    return x


def invert(a: Matrix) -> Matrix:
    """Matrix inverse by Gauss-Jordan elimination with partial pivoting."""

    n = len(a)
    m = [list(map(float, a[i])) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ArithmeticError("singular matrix")
        m[col], m[pivot] = m[pivot], m[col]
        scale = m[col][col]
        m[col] = [v / scale for v in m[col]]
        for row in range(n):
            if row != col and m[row][col]:
                factor = m[row][col]
                m[row] = [v - factor * w for v, w in zip(m[row], m[col], strict=False)]
    return [row[n:] for row in m]


def matvec(a: Matrix, x: Sequence[float]) -> list[float]:
    return [math.fsum(v * w for v, w in zip(row, x, strict=False)) for row in a]


# --------------------------------------------------------------------------- regression
def simple_regression(
    x: Sequence[float], y: Sequence[float], df_extra: int = 0
) -> dict[str, float]:
    """OLS of y on x with intercept: slope, standard error, t statistic and two-sided p.

    ``df_extra`` removes additional degrees of freedom (number of covariates already
    projected out).
    """

    n = len(x)
    mx, my = mean(x), mean(y)
    sxx = math.fsum((a - mx) ** 2 for a in x)
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y, strict=False))
    syy = math.fsum((b - my) ** 2 for b in y)
    df = n - 2 - df_extra
    if sxx <= 0 or df < 1:
        return {"beta": math.nan, "se": math.nan, "t": math.nan, "p": math.nan, "df": float(df)}
    beta = sxy / sxx
    rss = max(0.0, syy - beta * sxy)
    se = math.sqrt(rss / df / sxx)
    if se == 0.0:
        return {"beta": beta, "se": 0.0, "t": math.inf, "p": 0.0, "df": float(df)}
    t = beta / se
    return {"beta": beta, "se": se, "t": t, "p": t_two_sided_p(t, df), "df": float(df)}


def logistic_regression(
    design: Sequence[Sequence[float]], y: Sequence[float], *, max_iter: int = 50, tol: float = 1e-10
) -> dict[str, object]:
    """Maximum-likelihood logistic regression by IRLS (Newton-Raphson).  ``design`` includes the
    intercept column.
    """

    n, k = len(design), len(design[0])
    beta = [0.0] * k
    loglik_old = -math.inf
    converged = False
    cov: Matrix = []
    loglik = -math.inf
    for _ in range(max_iter):
        xtwx = [[0.0] * k for _ in range(k)]
        xtz = [0.0] * k
        loglik = 0.0
        for row, yi in zip(design, y, strict=False):
            eta = sum(b * v for b, v in zip(beta, row, strict=False))
            if eta >= 0:
                e = math.exp(-eta)
                p = 1.0 / (1.0 + e)
                loglik += -math.log1p(e) if yi == 1 else -eta - math.log1p(e)
            else:
                e = math.exp(eta)
                p = e / (1.0 + e)
                loglik += eta - math.log1p(e) if yi == 1 else -math.log1p(e)
            w = p * (1.0 - p)
            for i in range(k):
                xi = row[i]
                xtz[i] += xi * (yi - p)
                wi = w * xi
                for j in range(i, k):
                    xtwx[i][j] += wi * row[j]
        for i in range(k):
            for j in range(i):
                xtwx[i][j] = xtwx[j][i]
        try:
            step = solve_linear(xtwx, xtz)
            cov = invert(xtwx)
        except ArithmeticError as exc:
            raise ArithmeticError(
                "logistic design matrix is singular or the data are separated"
            ) from exc
        beta = [b + s for b, s in zip(beta, step, strict=False)]
        if max(abs(s) for s in step) < tol or abs(loglik - loglik_old) < tol:
            converged = True
            break
        loglik_old = loglik
        if max(abs(b) for b in beta) > 30:
            raise ArithmeticError("logistic regression diverged (perfect separation)")
    # log-likelihood at the final beta
    loglik = 0.0
    for row, yi in zip(design, y, strict=False):
        eta = sum(b * v for b, v in zip(beta, row, strict=False))
        loglik += (eta if yi == 1 else 0.0) - (max(eta, 0.0) + math.log1p(math.exp(-abs(eta))))
    se = [math.sqrt(max(cov[i][i], 0.0)) for i in range(k)]
    return {"beta": beta, "se": se, "cov": cov, "loglik": loglik, "converged": converged, "n": n}


def ols_multiple(design: Sequence[Sequence[float]], y: Sequence[float]) -> dict[str, object]:
    """OLS with an explicit design matrix (intercept supplied by the caller)."""

    n, k = len(design), len(design[0])
    xtx = [
        [math.fsum(design[r][i] * design[r][j] for r in range(n)) for j in range(k)]
        for i in range(k)
    ]
    xty = [math.fsum(design[r][i] * y[r] for r in range(n)) for i in range(k)]
    inv = invert(xtx)
    beta = matvec(inv, xty)
    rss = math.fsum(
        (y[r] - sum(b * v for b, v in zip(beta, design[r], strict=False))) ** 2 for r in range(n)
    )
    df = n - k
    sigma2 = rss / df if df > 0 else math.nan
    se = [math.sqrt(sigma2 * inv[i][i]) for i in range(k)]
    return {"beta": beta, "se": se, "rss": rss, "df": df, "sigma2": sigma2}


# --------------------------------------------------------------------------- multiple testing
def bonferroni(p_values: Sequence[float]) -> list[float]:
    m = len(p_values)
    return [min(1.0, p * m) for p in p_values]


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values (q-values controlling the false discovery rate)."""

    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    adjusted = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p_values[i] * m / rank)
        adjusted[i] = running
    return adjusted
