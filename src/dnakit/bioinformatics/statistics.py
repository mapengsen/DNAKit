"""Native distance statistics, ordination, composition and ecological diversity."""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import matrix, positive, require
from dnakit.exceptions import ConfigurationError


def distance_matrix(values: Any, *, verbose: bool = True) -> Any:
    np = require("numpy", "bio-analysis")
    data = matrix(values, min_rows=2)
    if (
        data.shape[0] != data.shape[1]
        or (data < 0).any()
        or not np.allclose(data, data.T)
        or not np.allclose(np.diag(data), 0)
    ):
        raise ConfigurationError(
            "Distance matrix must be square, symmetric, nonnegative, with zero diagonal."
        )
    return data


@dataclass(frozen=True)
class PermutationResult:
    method: str
    statistic: float
    pvalue: float | None
    permutations: int
    seed: int


def _permutation(
    method: str,
    statistic: Callable[[Any], float],
    observed_labels: Any,
    permutations: int,
    seed: int,
    *,
    alternative: str = "greater",
) -> PermutationResult:
    np = require("numpy", "bio-analysis")
    if type(permutations) is not int or not 0 <= permutations <= 100_000:
        raise ConfigurationError("permutations must be in [0, 100000].")
    observed = statistic(observed_labels)
    if not np.isfinite(observed):
        raise ConfigurationError("Statistic is undefined for these data (e.g. zero variance).")
    if alternative not in {"greater", "less", "two-sided"}:
        raise ConfigurationError("Unknown alternative hypothesis.")
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(permutations):
        value = statistic(rng.permutation(observed_labels))
        if alternative == "two-sided":
            count += abs(value) >= abs(observed)
        else:
            count += value >= observed if alternative == "greater" else value <= observed
    return PermutationResult(
        method,
        observed,
        (count + 1) / (permutations + 1) if permutations else None,
        permutations,
        seed,
    )


def distance_test(
    distances: Any,
    groups: Iterable[Any],
    *,
    method: str = "permanova",
    permutations: int = 999,
    seed: int = 0,
    verbose: bool = True,
) -> PermutationResult:
    """One-factor PERMANOVA, ANOSIM, or centroid PERMDISP with label permutations.

    Group order must match rows. PERMDISP requires Euclidean distances and
    at least two samples per group. No stratified/repeated-measures inference.
    """
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    data = distance_matrix(distances)
    labels = np.asarray(tuple(groups))
    if labels.shape != (len(data),):
        raise ConfigurationError("One group label is required per distance row.")
    _, labels = np.unique(labels, return_inverse=True)
    unique = np.unique(labels)
    if not 1 < len(unique) < len(labels):
        raise ConfigurationError("Need between two and sample_count - 1 groups.")
    upper = np.triu_indices(len(data), 1)
    ranked = stats.rankdata(data[upper])
    if method == "permdisp":
        ordination = pcoa(data)
        if ordination.negative_eigenvalues:
            raise ConfigurationError("Centroid PERMDISP currently requires Euclidean distances.")
        if any(sum(labels == group) < 2 for group in unique):
            raise ConfigurationError("PERMDISP requires two samples per group.")
        coordinates = ordination.samples
    elif method not in {"permanova", "anosim"}:
        raise ConfigurationError("method must be permanova, anosim, or permdisp.")

    def statistic(permuted: Any) -> float:
        if method == "anosim":
            same = permuted[upper[0]] == permuted[upper[1]]
            return float(
                (ranked[~same].mean() - ranked[same].mean()) / (len(data) * (len(data) - 1) / 4)
            )
        if method == "permdisp":
            samples = [coordinates[permuted == group] for group in unique]
            deviations = [
                np.linalg.norm(sample - sample.mean(axis=0), axis=1) for sample in samples
            ]
            return float(stats.f_oneway(*deviations).statistic)
        total = float((data**2).sum() / (2 * len(data)))
        within = sum(
            float((data[np.ix_(permuted == group, permuted == group)] ** 2).sum())
            / (2 * sum(permuted == group))
            for group in unique
        )
        if within == 0:
            return float("inf")
        return float(((total - within) / (len(unique) - 1)) / (within / (len(data) - len(unique))))

    return _permutation(method, statistic, labels, permutations, seed)


def mantel(
    left: Any,
    right: Any,
    *,
    method: str = "pearson",
    permutations: int = 999,
    seed: int = 0,
    alternative: str = "two-sided",
    verbose: bool = True,
) -> PermutationResult:
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    a, b = distance_matrix(left), distance_matrix(right)
    if a.shape != b.shape or method not in {"pearson", "spearman"}:
        raise ConfigurationError("Distance shapes must match; method is pearson or spearman.")
    upper = np.triu_indices(len(a), 1)
    correlation = stats.pearsonr if method == "pearson" else stats.spearmanr

    def statistic(indices: Any) -> float:
        return float(correlation(a[upper], b[np.ix_(indices, indices)][upper]).statistic)

    return _permutation(
        "mantel", statistic, np.arange(len(a)), permutations, seed, alternative=alternative
    )


@dataclass(frozen=True)
class OrdinationResult:
    method: str
    samples: Any
    eigenvalues: Any
    proportion_explained: Any
    features: Any = None
    negative_eigenvalues: tuple[float, ...] = ()


def pcoa(
    distances: Any, *, dimensions: int | None = None, verbose: bool = True
) -> OrdinationResult:
    """Classical principal coordinates; retain positive axes and report negative eigenvalues."""
    np = require("numpy", "bio-analysis")
    data = distance_matrix(distances)
    squared = data**2
    gram = -0.5 * (squared - squared.mean(0)[None, :] - squared.mean(1)[:, None] + squared.mean())
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    tolerance = max(1.0, float(np.abs(values).max())) * 1e-10
    negative = tuple(float(value) for value in values if value < -tolerance)
    mask = values > tolerance
    values, vectors = values[mask], vectors[:, mask]
    total = values.sum()
    if dimensions is not None:
        positive(dimensions, "dimensions")
        if dimensions > len(values):
            raise ConfigurationError("Requested more dimensions than positive PCoA axes.")
        values, vectors = values[:dimensions], vectors[:, :dimensions]
    return OrdinationResult(
        "pcoa",
        vectors * np.sqrt(values),
        values,
        values / total if total else values,
        negative_eigenvalues=negative,
    )


def ordination(
    values: Any, *, method: str = "ca", constraints: Any = None, verbose: bool = True
) -> OrdinationResult:
    """CA (row principal coordinates), CCA constrained axes, or centered RDA.

    CCA uses chi-square standardized responses projected on weighted centered
    constraints. Returns constrained axes only; residual axes are not appended.
    """
    np = require("numpy", "bio-analysis")
    data = matrix(values, min_rows=2)
    if method == "pcoa":
        return pcoa(data)
    if method not in {"ca", "cca", "rda"}:
        raise ConfigurationError("Unknown ordination method.")
    weights = np.full(len(data), 1 / len(data))
    feature_weights = np.ones(data.shape[1])
    if method in {"ca", "cca"}:
        if (data < 0).any() or (data.sum(0) == 0).any() or (data.sum(1) == 0).any():
            raise ConfigurationError(
                "Correspondence analysis requires nonnegative data with positive margins."
            )
        frequencies = data / data.sum()
        weights, feature_weights = frequencies.sum(1), frequencies.sum(0)
        response = (frequencies - weights[:, None] * feature_weights) / np.sqrt(
            weights[:, None] * feature_weights
        )
    else:
        response = data - data.mean(axis=0)
    total_inertia = float((response**2).sum()) / (len(data) - 1 if method == "rda" else 1)
    if method in {"cca", "rda"}:
        predictors = matrix(constraints, min_rows=2)
        if len(predictors) != len(data):
            raise ConfigurationError("Constraint rows must match response rows.")
        predictors = predictors - np.average(predictors, axis=0, weights=weights)
        if method == "cca":
            predictors = predictors * np.sqrt(weights[:, None])
        response = predictors @ np.linalg.lstsq(predictors, response, rcond=None)[0]
    u, singular, vt = np.linalg.svd(response, full_matrices=False)
    keep = singular > max(1.0, float(singular.max())) * 1e-10
    u, singular, vt = u[:, keep], singular[keep], vt[keep]
    eigenvalues = singular**2 / (len(data) - 1 if method == "rda" else 1)
    samples, features = u * singular, vt.T
    if method != "rda":
        samples = samples / np.sqrt(weights[:, None])
        features = features / np.sqrt(feature_weights[:, None])
    return OrdinationResult(
        method,
        samples,
        eigenvalues,
        eigenvalues / total_inertia if total_inertia else eigenvalues,
        features,
    )


def bioenv(
    distances: Any, environmental: Any, *, max_columns: int = 12, verbose: bool = True
) -> tuple[dict[str, Any], ...]:
    """Best standardized environmental-variable subset for each subset size (Spearman)."""
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    spatial = require("scipy.spatial.distance", "bio-analysis")
    data, env = distance_matrix(distances), matrix(environmental, min_rows=2)
    if len(data) != len(env) or env.shape[1] > min(positive(max_columns, "max_columns"), 16):
        raise ConfigurationError("Environment rows must match, with bounded column count.")
    sd = env.std(axis=0)
    if (sd == 0).any():
        raise ConfigurationError("Remove constant environmental columns.")
    env = (env - env.mean(axis=0)) / sd
    target = data[np.triu_indices(len(data), 1)]
    output = []
    for size in range(1, env.shape[1] + 1):
        candidates: list[dict[str, Any]] = []
        for columns in itertools.combinations(range(env.shape[1]), size):
            statistic = float(stats.spearmanr(target, spatial.pdist(env[:, columns])).statistic)
            if np.isfinite(statistic):
                candidates.append({"columns": columns, "correlation": statistic})
        if not candidates:
            raise ConfigurationError("BIOENV correlations are undefined for these distances.")
        output.append(max(candidates, key=lambda item: item["correlation"]))
    return tuple(output)


def composition_transform(
    values: Any,
    *,
    method: str = "clr",
    reference: int = 0,
    inverse: bool = False,
    verbose: bool = True,
) -> Any:
    """CLR/ILR/ALR and inverses; ILR uses an orthonormal Helmert basis.

    Zeros require an explicit upstream replacement choice; no pseudocount is
    silently added. CLR inverse uses closure(exp(x)); ALR omits reference.
    """
    np = require("numpy", "bio-analysis")
    special = require("scipy.special", "bio-analysis")
    linalg = require("scipy.linalg", "bio-analysis")
    data = matrix(values)
    if method not in {"clr", "ilr", "alr"}:
        raise ConfigurationError("Unknown log-ratio transform.")
    width = data.shape[1] + (1 if inverse and method in {"ilr", "alr"} else 0)
    if width < 2 or type(reference) is not int or not 0 <= reference < width:
        raise ConfigurationError("Need at least two parts and a valid reference index.")
    basis = linalg.helmert(width, full=False)
    if inverse:
        log_values = data @ basis if method == "ilr" else data
        if method == "alr":
            log_values = np.insert(data, reference, 0, axis=1)
        return np.exp(log_values - special.logsumexp(log_values, axis=1, keepdims=True))
    if (data <= 0).any():
        raise ConfigurationError(
            "Log-ratio input must be strictly positive; replace zeros explicitly."
        )
    logs = np.log(data)
    centered = logs - logs.mean(axis=1, keepdims=True)
    if method == "clr":
        return centered
    if method == "ilr":
        return centered @ basis.T
    return np.delete(logs - logs[:, reference, None], reference, axis=1)


def ancom(
    values: Any,
    groups: Iterable[Any],
    *,
    alpha: float = 0.05,
    pseudocount: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Native pairwise log-ratio ANCOM W statistics with per-taxon Holm correction.

    Uses two-sided Mann-Whitney for two groups and Kruskal-Wallis otherwise.
    Returns W and adjusted pairwise p-values; does not invent a W rejection cutoff.
    """
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    data = matrix(values, min_rows=3)
    labels = np.asarray(tuple(groups))
    if (
        labels.shape != (len(data),)
        or not 0 < alpha < 1
        or not np.isfinite(pseudocount)
        or pseudocount < 0
    ):
        raise ConfigurationError("Invalid ANCOM groups, alpha or pseudocount.")
    if (data < 0).any() or (data + pseudocount <= 0).any() or data.shape[1] < 2:
        raise ConfigurationError(
            "ANCOM needs two positive parts; choose zero replacement explicitly."
        )
    unique = np.unique(labels)
    if len(unique) < 2 or any(sum(labels == group) < 2 for group in unique):
        raise ConfigurationError("ANCOM needs at least two samples in each of two or more groups.")
    logs = np.log(data + pseudocount)
    pvalues = np.ones((data.shape[1], data.shape[1]))
    for left, right in itertools.combinations(range(data.shape[1]), 2):
        ratio = logs[:, left] - logs[:, right]
        samples = [ratio[labels == group] for group in unique]
        if np.ptp(ratio) == 0:
            pvalue = 1.0
        elif len(samples) == 2:
            pvalue = float(stats.mannwhitneyu(*samples, alternative="two-sided").pvalue)
        else:
            pvalue = float(stats.kruskal(*samples).pvalue)
        pvalues[left, right] = pvalues[right, left] = pvalue
    adjusted = np.ones_like(pvalues)
    for index in range(data.shape[1]):
        other = np.delete(np.arange(data.shape[1]), index)
        order = np.argsort(pvalues[index, other])
        corrected = np.minimum(
            np.maximum.accumulate(pvalues[index, other[order]] * np.arange(len(other), 0, -1)), 1
        )
        adjusted[index, other[order]] = corrected
    return {"W": (adjusted < alpha).sum(axis=1), "pvalues": pvalues, "adjusted_pvalues": adjusted}


def alpha_diversity(counts: Any, *, metric: str = "shannon", verbose: bool = True) -> Any:
    """Per-sample observed richness, Shannon (natural logs), Simpson, or inverse Simpson."""
    np = require("numpy", "bio-analysis")
    values = matrix(counts)
    if (values < 0).any() or (values.sum(1) == 0).any():
        raise ConfigurationError("Abundances must be nonnegative with positive sample totals.")
    frequency = values / values.sum(1, keepdims=True)
    if metric == "observed":
        return (values > 0).sum(1)
    if metric == "shannon":
        return -(frequency * np.log(np.where(frequency > 0, frequency, 1))).sum(1)
    if metric == "simpson":
        return 1 - (frequency**2).sum(1)
    if metric == "inverse_simpson":
        return 1 / (frequency**2).sum(1)
    raise ConfigurationError("Unknown alpha diversity metric.")


def beta_diversity(counts: Any, *, metric: str = "braycurtis", verbose: bool = True) -> Any:
    spatial = require("scipy.spatial.distance", "bio-analysis")
    values = matrix(counts, min_rows=2)
    if (values < 0).any() or (values.sum(1) == 0).any():
        raise ConfigurationError("Abundances must be nonnegative with positive sample totals.")
    if metric not in {"braycurtis", "jaccard", "euclidean", "cityblock", "cosine"}:
        raise ConfigurationError("Unknown beta diversity metric.")
    if metric == "jaccard":
        values = values > 0
    return spatial.squareform(spatial.pdist(values, metric=metric))


_register_public_verbose(
    globals(),
    (
        "distance_matrix",
        "distance_test",
        "mantel",
        "pcoa",
        "ordination",
        "bioenv",
        "composition_transform",
        "ancom",
        "alpha_diversity",
        "beta_diversity",
    ),
)
