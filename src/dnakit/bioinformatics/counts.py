"""Native count inference using explicit designs, contrasts and numerical diagnostics.

ANCOM-BC implements the original three-component sampling-bias model. RNA count
inference uses an unshrunk NB2 GLM; it is not a replacement implementation of
DESeq2/edgeR and does not estimate their empirical-Bayes dispersion trends.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from dnakit.bioinformatics._common import matrix, positive, require
from dnakit.exceptions import ConfigurationError


def adjust_pvalues(values: Any) -> Any:
    """Benjamini-Hochberg FDR correction of a finite one-dimensional p-value array."""
    np = require("numpy", "bio-analysis")
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ConfigurationError("p-values must be a nonempty finite vector in [0,1].")
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    output = np.empty_like(q)
    output[order] = np.minimum(q, 1)
    return output


def _design(values: Any, samples: int) -> Any:
    np = require("numpy", "bio-analysis")
    design = matrix(values, min_rows=3)
    if len(design) != samples or len(design) <= design.shape[1]:
        raise ConfigurationError(
            "Design rows must match samples and leave residual degrees of freedom."
        )
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ConfigurationError(
            "Design must have full column rank; include an intercept explicitly."
        )
    return design


@dataclass(frozen=True)
class AbundanceResult:
    """Feature-by-design-column natural-log effects, robust SE, Wald and BH p-values."""

    coefficients: Any
    standard_errors: Any
    statistics: Any
    pvalues: Any
    qvalues: Any
    sampling_bias: Any
    iterations: tuple[int, ...]


def _sampling_bias(
    beta: Any, variance: Any, tolerance: float, iterations: int
) -> tuple[float, int]:
    np = require("numpy", "bio-analysis")
    special = require("scipy.special", "bio-analysis")
    optimize = require("scipy.optimize", "bio-analysis")
    quantiles = np.quantile(beta, [0.125, 0.25, 0.75, 0.875])
    center = beta[(beta >= quantiles[1]) & (beta <= quantiles[2])]
    delta = float(center.mean())
    tails = (beta[beta < quantiles[0]], beta[beta > quantiles[3]])
    offsets = np.array(
        [
            0.0,
            *(
                float(t.mean()) if len(t) else float(beta.min() if i == 0 else beta.max())
                for i, t in enumerate(tails)
            ),
        ]
    )
    excess = np.array(
        [0.0, *(max(float(t.var(ddof=1)), 1e-12) if len(t) > 1 else 1.0 for t in tails)]
    )
    weights = np.array([0.75, 0.125, 0.125])
    for iteration in range(1, iterations + 1):
        previous = np.r_[weights, delta, offsets[1:], excess[1:]]
        variances = variance[None, :] + excess[:, None]
        means = delta + offsets
        logp = np.log(np.maximum(weights, 1e-300))[:, None] - 0.5 * (
            np.log(variances) + (beta[None, :] - means[:, None]) ** 2 / variances
        )
        responsibilities = np.exp(logp - special.logsumexp(logp, axis=0))
        precision = responsibilities / variances
        new_delta = float((precision * (beta - offsets[:, None])).sum() / precision.sum())
        new_offsets = np.sum(precision[1:] * (beta - delta), axis=1) / precision[1:].sum(axis=1)
        for component in (1, 2):

            def objective(
                value: Any,
                component: int = component,
                responsibilities: Any = responsibilities,
                means: Any = means,
            ) -> float:
                spread = variance + float(value[0])
                return float(
                    0.5
                    * np.sum(
                        responsibilities[component]
                        * (np.log(spread) + (beta - means[component]) ** 2 / spread)
                    )
                )

            fit = optimize.minimize(
                objective,
                [excess[component]],
                method="Nelder-Mead",
                bounds=[(0.0, None)],
                options={"xatol": 1e-8, "fatol": 1e-8},
            )
            if not fit.success:
                raise ConfigurationError("ANCOM-BC mixture variance optimization did not converge.")
            excess[component] = fit.x[0]
        weights = np.maximum(responsibilities.mean(axis=1), 1e-12)
        weights /= weights.sum()
        offsets[1:] = [min(float(new_offsets[0]), 0.0), max(float(new_offsets[1]), 0.0)]
        delta = new_delta
        if np.linalg.norm(np.r_[weights, delta, offsets[1:], excess[1:]] - previous) <= tolerance:
            return delta, iteration
    raise ConfigurationError("ANCOM-BC sampling-bias EM did not converge; increase max_iterations.")


def ancom_bc(
    counts: Any, design: Any, *, tolerance: float = 1e-5, max_iterations: int = 1000
) -> AbundanceResult:
    """Original ANCOM-BC log-linear model with three-component EM bias correction.

    Input is samples by taxa; zeros require explicit upstream handling. Design
    columns are supplied numerically, including the intercept. At least ten
    taxa and positive residual variances are required for bias estimation.
    BH correction is per design column. No structural-zero detection, repeated
    measures, global test, ANCOM-BC2 or variance-of-bias adjustment is implied.
    """
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    data = matrix(counts, min_rows=3)
    x = _design(design, len(data))
    positive(max_iterations, "max_iterations")
    if not 0 < tolerance < 1 or (data <= 0).any() or data.shape[1] < 10:
        raise ConfigurationError(
            "ANCOM-BC requires positive counts, at least ten taxa, and tolerance in (0,1)."
        )
    inverse = np.linalg.pinv(x)
    response = np.log(data)
    beta = inverse @ response
    residual = response - x @ beta
    residual -= residual.mean(axis=1, keepdims=True)
    # Heteroscedastic sandwich diagonal: sum_s (X^+_{js} e_{si})^2.
    variance = (inverse**2 @ residual**2).T
    if (variance <= 1e-14).any():
        raise ConfigurationError(
            "ANCOM-BC has zero residual coefficient variance; remove degenerate taxa/designs."
        )
    estimated = [
        _sampling_bias(beta[j], variance[:, j], tolerance, max_iterations)
        for j in range(x.shape[1])
    ]
    bias = np.array([item[0] for item in estimated])
    corrected = beta.T - bias
    se = np.sqrt(variance)
    wald = corrected / se
    p = 2 * stats.norm.sf(np.abs(wald))
    q = np.column_stack([adjust_pvalues(column) for column in p.T])
    return AbundanceResult(corrected, se, wald, p, q, bias, tuple(item[1] for item in estimated))


@dataclass(frozen=True)
class DifferentialExpressionResult:
    """Gene-order results for an explicit contrast, plus model diagnostics."""

    log2_fold_change: Any
    standard_errors: Any
    statistics: Any
    pvalues: Any
    qvalues: Any
    dispersions: Any
    size_factors: Any
    coefficients: Any
    converged: Any
    method: str = "native_nb2_unshrunk"


def differential_expression(
    counts: Any,
    design: Any,
    contrast: Any,
    *,
    size_factors: Any = None,
    dispersions: Any = None,
    max_iterations: int = 1000,
) -> DifferentialExpressionResult:
    """Unshrunk NB2 GLMs with library offsets, Wald contrasts and BH correction.

    Samples-by-genes nonnegative integer counts; include the design intercept
    explicitly. Supply biological replicates (at least two residual degrees of
    freedom). Default size factors use median ratios over genes positive in all
    samples. Default dispersion is a within-design residual method-of-moments
    estimate with floor 1e-8. For small studies prefer externally validated
    dispersions. All-zero genes and failed fits have NaN inference, never p=0.
    """
    np = require("numpy", "bio-analysis")
    optimize = require("scipy.optimize", "bio-analysis")
    special = require("scipy.special", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    y = matrix(counts, min_rows=4)
    x = _design(design, len(y))
    positive(max_iterations, "max_iterations")
    c = np.asarray(contrast, dtype=float)
    if len(y) - x.shape[1] < 2 or (y < 0).any() or (y != np.floor(y)).any():
        raise ConfigurationError("Use integer counts and at least two residual degrees of freedom.")
    if c.shape != (x.shape[1],) or not np.isfinite(c).all() or not np.any(c):
        raise ConfigurationError(
            "Contrast must be a finite nonzero vector matching design columns."
        )
    if size_factors is None:
        eligible = (y > 0).all(axis=0)
        if not eligible.any():
            raise ConfigurationError(
                "No gene is positive in every sample; provide explicit size_factors."
            )
        logs = np.log(y[:, eligible])
        factors = np.median(np.exp(logs - logs.mean(axis=0)), axis=1)
    else:
        factors = np.asarray(size_factors, dtype=float)
    if factors.shape != (len(y),) or not np.isfinite(factors).all() or (factors <= 0).any():
        raise ConfigurationError("size_factors must be positive and match sample rows.")
    factors = factors / np.exp(np.log(factors).mean())
    normalized = y / factors[:, None]
    if dispersions is None:
        fitted = x @ np.linalg.lstsq(x, normalized, rcond=None)[0]
        residual_variance = ((normalized - fitted) ** 2).sum(axis=0) / (len(y) - x.shape[1])
        mean = normalized.mean(axis=0)
        dispersion = np.maximum(
            (residual_variance - (normalized / factors[:, None]).mean(axis=0))
            / np.maximum(mean**2, 1e-30),
            1e-8,
        )
    else:
        dispersion = np.broadcast_to(np.asarray(dispersions, dtype=float), (y.shape[1],)).copy()
    if not np.isfinite(dispersion).all() or (dispersion <= 0).any():
        raise ConfigurationError("Dispersions must be positive and match genes.")
    coefficients = np.full((y.shape[1], x.shape[1]), np.nan)
    effect, se, statistic, p = (np.full(y.shape[1], np.nan) for _ in range(4))
    converged = np.zeros(y.shape[1], dtype=bool)
    offset = np.log(factors)
    for gene in range(y.shape[1]):
        observed = y[:, gene]
        if not observed.any():
            continue
        shape = 1 / dispersion[gene]
        log_shape = np.log(shape)

        def objective(
            beta: Any, observed: Any = observed, shape: float = shape, log_shape: float = log_shape
        ) -> tuple[float, Any]:
            eta = x @ beta + offset
            normalizer = np.logaddexp(log_shape, eta)
            loss = np.sum((observed + shape) * normalizer - observed * eta)
            gradient = x.T @ ((observed + shape) * special.expit(eta - log_shape) - observed)
            return float(loss), gradient

        initial = np.linalg.lstsq(x, np.log(normalized[:, gene] + 0.5), rcond=None)[0]
        fit = optimize.minimize(
            objective,
            initial,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": max_iterations, "ftol": 1e-12, "gtol": 1e-7},
        )
        if not fit.success or not np.isfinite(fit.x).all():
            continue
        eta = x @ fit.x + offset
        weights = shape * special.expit(eta - log_shape)
        information = x.T @ (weights[:, None] * x)
        if np.linalg.matrix_rank(information) != x.shape[1]:
            continue
        variance = float(c @ np.linalg.solve(information, c))
        if variance <= 0 or not np.isfinite(variance):
            continue
        converged[gene] = True
        coefficients[gene] = fit.x
        effect[gene] = c @ fit.x / np.log(2)
        se[gene] = np.sqrt(variance) / np.log(2)
        statistic[gene] = effect[gene] / se[gene]
        p[gene] = 2 * stats.norm.sf(abs(statistic[gene]))
    q = np.full_like(p, np.nan)
    if converged.any():
        q[converged] = adjust_pvalues(p[converged])
    return DifferentialExpressionResult(
        effect, se, statistic, p, q, dispersion, factors, coefficients, converged
    )
