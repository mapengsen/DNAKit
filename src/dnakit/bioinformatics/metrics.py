"""Generative metrics on caller-supplied embeddings, properties, or sequences."""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from itertools import pairwise
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import matrix, positive, require
from dnakit.exceptions import ConfigurationError


def sequence_properties(
    sequences: Mapping[str, Any],
    predictors: Mapping[str, Callable[[str], float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run scalar Python predictors with stable sequence/property order.

    Callbacks can wrap existing DNAKit predictors or user models. Returned
    values feed property_hit_rate, pareto_ranks, hypervolume and property_kl.
    The caller controls model fitting; this function never trains on inputs.
    """
    from dnakit.bioinformatics._common import sequence_text

    if not sequences or not predictors or any(not callable(p) for p in predictors.values()):
        raise ConfigurationError("Provide named sequences and callable scalar predictors.")
    texts = tuple(sequence_text(value) for value in sequences.values())
    values = matrix([[predictor(text) for predictor in predictors.values()] for text in texts])
    return {"ids": tuple(sequences), "names": tuple(predictors), "values": values}


def kernel_distance(
    generated: Any,
    reference: Any,
    *,
    metric: str = "mmd",
    sigma: float = 1.0,
    unbiased: bool = True,
    max_pairs: int = 4_000_000,
    verbose: bool = True,
) -> float:
    """MMD squared (Gaussian kernel) or KID (cubic polynomial kernel).

    No rescaling, subsampling, square root, or clipping is performed. Unbiased
    estimates may be negative. KID uses (x@y / dimension + 1)**3.
    """
    np = require("numpy", "bio-analysis")
    distance = require("scipy.spatial.distance", "bio-analysis")
    x = matrix(generated, min_rows=2 if unbiased else 1)
    y = matrix(reference, min_rows=2 if unbiased else 1)
    positive(max_pairs, "max_pairs")
    if x.shape[1] != y.shape[1]:
        raise ConfigurationError("Embedding dimensions must match.")
    if len(x) ** 2 + len(y) ** 2 + len(x) * len(y) > max_pairs:
        raise ConfigurationError("Kernel evaluation exceeds max_pairs; subsample explicitly.")
    if metric not in {"mmd", "kid"}:
        raise ConfigurationError("metric must be mmd or kid.")
    if not math.isfinite(sigma) or sigma <= 0:
        raise ConfigurationError("sigma must be finite and positive.")

    def kernel(a: Any, b: Any) -> Any:
        if metric == "kid":
            return ((a @ b.T) / a.shape[1] + 1.0) ** 3
        squared = distance.cdist(a, b, metric="sqeuclidean")
        return np.exp(-0.5 * (squared / sigma) / sigma)

    xx, yy, xy = kernel(x, x), kernel(y, y), kernel(x, y)
    if unbiased:
        value = (
            (xx.sum() - np.trace(xx)) / (len(x) * (len(x) - 1))
            + (yy.sum() - np.trace(yy)) / (len(y) * (len(y) - 1))
            - 2 * xy.mean()
        )
    else:
        value = xx.mean() + yy.mean() - 2 * xy.mean()
    if not np.isfinite(value):
        raise ConfigurationError("Kernel result overflowed; rescale the input embeddings.")
    if not np.isfinite(value):
        raise ConfigurationError("Kernel computation overflowed; rescale embeddings/bandwidth.")
    return float(value)


def manifold_metrics(
    generated: Any,
    reference: Any,
    *,
    neighbors: int = 5,
    max_pairs: int = 4_000_000,
    verbose: bool = True,
) -> dict[str, float]:
    """Euclidean kNN precision/recall and PRDC density/coverage.

    Self-neighbors are excluded from manifold radii. These are the conventional
    PRDC definitions, not seqme's ClippedDensity/ClippedCoverage variants.
    """
    np = require("numpy", "bio-analysis")
    distance = require("scipy.spatial.distance", "bio-analysis")
    k = positive(neighbors, "neighbors")
    x, y = matrix(generated, min_rows=k + 1), matrix(reference, min_rows=k + 1)
    if x.shape[1] != y.shape[1]:
        raise ConfigurationError("Embedding dimensions must match.")
    if len(x) ** 2 + len(y) ** 2 + len(x) * len(y) > positive(max_pairs, "max_pairs"):
        raise ConfigurationError("Manifold evaluation exceeds max_pairs; subsample explicitly.")
    xx, yy, xy = distance.cdist(x, x), distance.cdist(y, y), distance.cdist(x, y)
    np.fill_diagonal(xx, np.inf)
    np.fill_diagonal(yy, np.inf)
    radius_x = np.partition(xx, k - 1, axis=1)[:, k - 1]
    radius_y = np.partition(yy, k - 1, axis=1)[:, k - 1]
    in_reference = xy < radius_y[None, :]
    in_generated = xy < radius_x[:, None]
    return {
        "precision": float(in_reference.any(axis=1).mean()),
        "recall": float(in_generated.any(axis=0).mean()),
        "density": float(in_reference.sum(axis=1).mean() / k),
        "coverage": float((xy.min(axis=0) < radius_y).mean()),
    }


def pareto_ranks(
    values: Any, *, maximize: Iterable[bool], max_pairs: int = 4_000_000, verbose: bool = True
) -> Any:
    """Zero-based nondominated front numbers; preserve input row order and ties."""
    np = require("numpy", "bio-analysis")
    points = matrix(values)
    directions = tuple(maximize)
    if len(directions) != points.shape[1] or any(type(value) is not bool for value in directions):
        raise ConfigurationError("maximize must contain one boolean per objective.")
    if len(points) ** 2 > positive(max_pairs, "max_pairs"):
        raise ConfigurationError("Nondominated sorting exceeds max_pairs.")
    points = points * np.asarray([-1.0 if value else 1.0 for value in directions])
    dominates = np.zeros((len(points), len(points)), dtype=bool)
    for index, point in enumerate(points):
        dominates[index] = np.all(point <= points, axis=1) & np.any(point < points, axis=1)
    counts = dominates.sum(axis=0)
    ranks = np.full(len(points), -1, dtype=int)
    front = np.flatnonzero(counts == 0)
    rank = 0
    while len(front):
        ranks[front] = rank
        counts -= dominates[front].sum(axis=0)
        front = np.flatnonzero((counts == 0) & (ranks < 0))
        rank += 1
    return ranks


def hypervolume(
    values: Any, reference_point: Any, *, maximize: Iterable[bool], verbose: bool = True
) -> float:
    """Exact dominated volume by recursive axis slicing (up to five objectives)."""
    np = require("numpy", "bio-analysis")
    points = matrix(values)
    directions = tuple(maximize)
    reference = np.asarray(reference_point, dtype=float)
    if (
        len(directions) != points.shape[1]
        or reference.shape != (points.shape[1],)
        or not np.isfinite(reference).all()
        or any(type(value) is not bool for value in directions)
    ):
        raise ConfigurationError("Reference and directions must match the objective dimension.")
    signs = np.asarray([-1.0 if value else 1.0 for value in directions])
    points, reference = points * signs, reference * signs
    if (points > reference).any():
        raise ConfigurationError("Every point must be no worse than the hypervolume reference.")
    if points.shape[1] > 5 or len(points) > 1000:
        raise ConfigurationError("Exact hypervolume is limited to five objectives and 1000 points.")
    work = 0

    def volume(rows: Any, boundary: Any) -> float:
        nonlocal work
        work += len(rows)
        if work > 1_000_000:
            raise ConfigurationError("Exact hypervolume exceeded its subdivision limit.")
        if not len(rows):
            return 0.0
        if rows.shape[1] == 1:
            return float(boundary[0] - rows[:, 0].min())
        coordinates = np.unique(np.append(rows[:, 0], boundary[0]))
        return sum(
            float(hi - lo) * volume(rows[rows[:, 0] <= lo, 1:], boundary[1:])
            for lo, hi in pairwise(coordinates)
        )

    return volume(np.unique(points, axis=0), reference)


def property_hit_rate(
    values: Any, *, lower: Any, upper: Any, verbose: bool = True
) -> dict[str, Any]:
    """Inclusive property bounds, per-objective fractions, and joint hit rate."""
    np = require("numpy", "bio-analysis")
    points = matrix(values)
    lo, hi = np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)
    if lo.shape != (points.shape[1],) or hi.shape != lo.shape:
        raise ConfigurationError("Bounds must contain one value per property.")
    if np.isnan(lo).any() or np.isnan(hi).any() or (lo > hi).any():
        raise ConfigurationError("Property bounds must be ordered and cannot contain NaN.")
    passed = (points >= lo) & (points <= hi)
    return {
        "joint": float(passed.all(axis=1).mean()),
        "per_property": passed.mean(axis=0),
        "mask": passed,
    }


def authenticity(generated: Any, training: Any, *, verbose: bool = True) -> float:
    """Fraction whose nearest training neighbor is closer to another training point.

    The training set must be the generator's training data, not a held-out
    reference set. Exact duplicates have zero authenticity contribution.
    """
    np = require("numpy", "bio-analysis")
    spatial = require("scipy.spatial", "bio-analysis")
    x, y = matrix(generated), matrix(training, min_rows=2)
    if x.shape[1] != y.shape[1]:
        raise ConfigurationError("Embedding dimensions must match.")
    tree = spatial.cKDTree(y)
    own, _ = tree.query(y, k=2)
    distance, nearest = tree.query(x, k=1)
    return float(np.mean(distance > own[nearest, 1]))


def clipped_manifold_metrics(
    generated: Any,
    reference: Any,
    *,
    neighbors: int = 5,
    max_pairs: int = 4_000_000,
    verbose: bool = True,
) -> dict[str, float]:
    """Native Salvy et al. (2025) clipped density and calibrated clipped coverage.

    Closed balls, median-clipped radii for density, leave-one-out density
    calibration, and beta-binomial coverage calibration. Equal sample counts
    are required. Calibration includes both endpoints and clips results to [0,1].
    """
    np = require("numpy", "bio-analysis")
    stats = require("scipy.stats", "bio-analysis")
    spatial = require("scipy.spatial.distance", "bio-analysis")
    k = positive(neighbors, "neighbors")
    x, y = matrix(generated), matrix(reference, min_rows=k + 1)
    if x.shape != y.shape:
        raise ConfigurationError("Clipped metrics require matching sample counts and dimensions.")
    if 2 * len(y) ** 2 > positive(max_pairs, "max_pairs"):
        raise ConfigurationError("Clipped evaluation exceeds max_pairs.")
    yy, xy = spatial.cdist(y, y), spatial.cdist(x, y)
    np.fill_diagonal(yy, np.inf)
    radii = np.partition(yy, k - 1, axis=1)[:, k - 1]
    clipped = np.minimum(radii, np.median(radii))
    calibration = float(np.minimum((yy <= clipped[None, :]).sum(1) / k, 1).mean())
    if calibration == 0:
        raise ConfigurationError("Clipped-density calibration is zero for these samples.")
    density = float(np.minimum((xy <= clipped[None, :]).sum(1) / k, 1).mean()) / calibration
    coverage = float(np.minimum((xy <= radii[None, :]).sum(0) / k, 1).mean())
    counts = np.arange(len(x) + 1)[:, None]
    low_counts = np.arange(k)[None, :]
    expected = 1 - np.sum(
        (1 - low_counts / k) * stats.betabinom.pmf(low_counts, counts, k, len(y) - k), axis=1
    )
    expected[0] = 0
    expected = np.maximum.accumulate(expected)
    index = min(int(np.searchsorted(expected, coverage, side="left")), len(x))
    return {"clipped_density": min(density, 1.0), "clipped_coverage": index / len(x)}


def property_kl(
    generated: Any,
    reference: Any,
    *,
    bins: int = 20,
    pseudocount: float = 1e-8,
    verbose: bool = True,
) -> Any:
    """Per-property KL(generated || reference), natural logs and shared histogram edges.

    Both tails are represented using edges spanning the combined samples.
    The explicit pseudocount is added to each bin count before normalization.
    """
    np = require("numpy", "bio-analysis")
    x, y = matrix(generated), matrix(reference)
    positive(bins, "bins")
    if x.shape[1] != y.shape[1] or not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ConfigurationError("Dimensions must match and pseudocount must be positive.")
    results = []
    for index in range(x.shape[1]):
        edges = np.histogram_bin_edges(np.concatenate((x[:, index], y[:, index])), bins=bins)
        p = np.histogram(x[:, index], bins=edges)[0] + pseudocount
        q = np.histogram(y[:, index], bins=edges)[0] + pseudocount
        p, q = p / p.sum(), q / q.sum()
        results.append(float(np.sum(p * np.log(p / q))))
    return np.asarray(results)


def conformity_score(
    generated: Any,
    reference: Any,
    *,
    folds: int = 5,
    bandwidth: float = 1.0,
    seed: int = 0,
    max_pairs: int = 4_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cross-fitted Gaussian-KDE conformity of property vectors.

    Each fold fits only on the other reference folds. Per generated sample,
    count validation densities no greater than its density and divide by
    validation_count + 1. Return fold mean, standard error, and fold scores.
    """
    np = require("numpy", "bio-analysis")
    distance = require("scipy.spatial.distance", "bio-analysis")
    special = require("scipy.special", "bio-analysis")
    x, y = matrix(generated), matrix(reference, min_rows=2)
    if type(folds) is not int or not 2 <= folds <= len(y):
        raise ConfigurationError("folds must be in [2, reference_count].")
    if x.shape[1] != y.shape[1] or not math.isfinite(bandwidth) or bandwidth <= 0:
        raise ConfigurationError("Dimensions must match and bandwidth must be positive.")
    if folds * len(x) * len(y) + len(y) ** 2 > positive(max_pairs, "max_pairs"):
        raise ConfigurationError("KDE evaluation exceeds max_pairs.")
    indices = np.random.default_rng(seed).permutation(len(y))
    scores = []
    for validation in np.array_split(indices, folds):
        train = y[np.setdiff1d(indices, validation)]

        def log_density(points: Any, train: Any = train) -> Any:
            squared = distance.cdist(points, train, metric="sqeuclidean")
            return special.logsumexp(-squared / (2 * bandwidth**2), axis=1) - math.log(len(train))

        val_density = np.sort(log_density(y[validation]))
        ranks = np.searchsorted(val_density, log_density(x), side="right")
        scores.append(float(np.mean(ranks / (len(validation) + 1))))
    return {
        "value": float(np.mean(scores)),
        "standard_error": float(np.std(scores) / math.sqrt(folds)),
        "fold_scores": tuple(scores),
    }


_register_public_verbose(
    globals(),
    (
        "sequence_properties",
        "kernel_distance",
        "manifold_metrics",
        "pareto_ranks",
        "hypervolume",
        "property_hit_rate",
        "authenticity",
        "clipped_manifold_metrics",
        "property_kl",
        "conformity_score",
    ),
)
