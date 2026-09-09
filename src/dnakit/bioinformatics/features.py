"""Traditional DNA encodings and fitted, reusable feature analysis."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from dnakit.bioinformatics._common import matrix, positive, require
from dnakit.bioinformatics.encoding import DNA_ENCODINGS, EncodingResult, encode_dna
from dnakit.exceptions import ConfigurationError

__all__ = [
    "DNA_ENCODINGS",
    "EncodingResult",
    "FeatureAnalysis",
    "analyze_features",
    "cluster_features",
    "encode_dna",
]


@dataclass(frozen=True)
class FeatureAnalysis:
    """Fitted backend plus transformed training values; call transform on held-out data."""

    values: Any
    estimator: Any
    input_names: tuple[str, ...]
    output_names: tuple[str, ...]
    method: str

    def transform(self, values: Any, *, feature_names: Sequence[str] | None = None) -> Any:
        data = matrix(values)
        if data.shape[1] != len(self.input_names):
            raise ConfigurationError("Feature width differs from fitted training data.")
        if feature_names is not None and tuple(feature_names) != self.input_names:
            raise ConfigurationError("Feature order differs from fitted training data.")
        if not hasattr(self.estimator, "transform"):
            raise ConfigurationError(
                f"{self.method} has no out-of-sample transform; refitting is required."
            )
        return self.estimator.transform(data)


def analyze_features(
    values: Any,
    *,
    method: str = "pca",
    labels: Any = None,
    feature_names: Sequence[str] | None = None,
    scale: str = "standard",
    components: int = 2,
    selected_features: int = 10,
    seed: int = 0,
    parameters: Mapping[str, Any] | None = None,
) -> FeatureAnalysis:
    """Fit scaling followed by reduction or selection on training rows only.

    Methods: pca, tsne, lda (supervised), latent_dirichlet, chi2,
    f_classif, f_regression, mutual_info_classif, mutual_info_regression, information_gain,
    pearson, variance, and none. chi2/latent_dirichlet require nonnegative
    input after scaling; use scale='none' or 'minmax' when appropriate.
    """
    from functools import partial

    data = matrix(values, min_rows=2)
    np = require("numpy", "bio-analysis")
    preprocessing = require("sklearn.preprocessing", "bio-analysis")
    selection = require("sklearn.feature_selection", "bio-analysis")
    pipeline = require("sklearn.pipeline", "bio-analysis")
    names = (
        tuple(feature_names)
        if feature_names is not None
        else tuple(f"feature_{index}" for index in range(data.shape[1]))
    )
    if len(names) != data.shape[1] or len(set(names)) != len(names):
        raise ConfigurationError("Feature names must be unique and match the matrix width.")
    scalers = {
        "standard": preprocessing.StandardScaler,
        "minmax": preprocessing.MinMaxScaler,
        "l1": partial(preprocessing.Normalizer, norm="l1"),
        "l2": partial(preprocessing.Normalizer, norm="l2"),
    }
    if scale not in {*scalers, "none"}:
        raise ConfigurationError("Unknown feature scaling method.")
    options = dict(parameters or {})
    if method in {"pca", "tsne", "lda", "latent_dirichlet"}:
        positive(components, "components")
        modules = {
            "pca": ("sklearn.decomposition", "PCA"),
            "tsne": ("sklearn.manifold", "TSNE"),
            "lda": ("sklearn.discriminant_analysis", "LinearDiscriminantAnalysis"),
            "latent_dirichlet": ("sklearn.decomposition", "LatentDirichletAllocation"),
        }
        module, name = modules[method]
        options.setdefault("n_components", components)
        if method != "lda":
            options.setdefault("random_state", seed)
        estimator = getattr(require(module, "bio-analysis"), name)(**options)
    elif method == "variance":
        estimator = selection.VarianceThreshold(**options)
    elif method == "none":
        estimator = preprocessing.FunctionTransformer(validate=True)
    else:

        def information_gain(x: Any, y: Any) -> Any:
            """Empirical label entropy reduction after training-only quantile binning."""
            target = np.asarray(y)
            if target.shape != (len(x),):
                raise ConfigurationError("Information gain requires one class label per row.")
            _, classes = np.unique(target, return_inverse=True)
            output = []
            for column in x.T:
                edges = np.unique(np.quantile(column, np.linspace(0, 1, bins + 1)[1:-1]))
                encoded = np.searchsorted(edges, column, side="right")
                joint = np.zeros((len(edges) + 1, int(classes.max()) + 1))
                np.add.at(joint, (encoded, classes), 1)
                joint /= len(x)
                independent = joint.sum(axis=1, keepdims=True) * joint.sum(axis=0)
                present = joint > 0
                output.append(
                    float(np.sum(joint[present] * np.log2(joint[present] / independent[present])))
                )
            return np.asarray(output)

        bins = positive(options.pop("bins", 10), "bins") if method == "information_gain" else 10
        scores = {
            "chi2": selection.chi2,
            "f_classif": selection.f_classif,
            "f_regression": selection.f_regression,
            "mutual_info_classif": partial(selection.mutual_info_classif, random_state=seed),
            "mutual_info_regression": partial(selection.mutual_info_regression, random_state=seed),
            "pearson": lambda x, y: np.abs(selection.r_regression(x, y)),
            "information_gain": information_gain,
        }
        if method not in scores:
            raise ConfigurationError(f"Unknown feature analysis method {method!r}.")
        positive(selected_features, "selected_features")
        if selected_features > data.shape[1]:
            raise ConfigurationError("selected_features exceeds the input width.")
        if options:
            raise ConfigurationError(
                "Selection methods use selected_features, not extra parameters."
            )
        estimator = selection.SelectKBest(scores[method], k=selected_features)
    if (
        method
        in {
            "lda",
            "chi2",
            "f_classif",
            "f_regression",
            "mutual_info_classif",
            "mutual_info_regression",
            "pearson",
            "information_gain",
        }
        and labels is None
    ):
        raise ConfigurationError(f"{method} requires training labels.")
    steps = [] if scale == "none" else [("scale", scalers[scale]())]
    fitted = pipeline.Pipeline([*steps, ("analysis", estimator)])
    try:
        transformed = fitted.fit_transform(data, labels)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"Feature analysis failed: {exc}") from exc
    if hasattr(estimator, "get_support"):
        output_names = tuple(names[index] for index in estimator.get_support(indices=True))
    elif method == "none":
        output_names = names
    else:
        output_names = tuple(f"{method}_{index}" for index in range(transformed.shape[1]))
    return FeatureAnalysis(transformed, fitted, names, output_names, method)


def cluster_features(
    values: Any,
    *,
    method: str = "kmeans",
    seed: int = 0,
    parameters: Mapping[str, Any] | None = None,
) -> Any:
    """Return a fitted sklearn clusterer, or MCL index groups for an affinity matrix.

    No automatic scaling is applied. sklearn results expose labels_ (GMM:
    predict). MCL expects a nonnegative symmetric affinity matrix.
    """
    data = matrix(values, min_rows=2)
    options = dict(parameters or {})
    if method == "mcl":
        np = require("numpy", "bio-analysis")
        if data.shape[0] != data.shape[1] or (data < 0).any() or not np.allclose(data, data.T):
            raise ConfigurationError("MCL requires a symmetric nonnegative affinity matrix.")
        if unknown := options.keys() - {
            "expansion",
            "inflation",
            "iterations",
            "tolerance",
            "pruning",
        }:
            raise ConfigurationError(f"Unknown MCL parameters: {sorted(unknown)}.")
        expansion = positive(options.get("expansion", 2), "expansion")
        inflation = float(options.get("inflation", 2.0))
        iterations = positive(options.get("iterations", 100), "iterations")
        tolerance = float(options.get("tolerance", 1e-6))
        pruning = float(options.get("pruning", 1e-6))
        if (
            not 2 <= expansion <= 10
            or not 1 < inflation <= 10
            or not 0 < tolerance < 1
            or not 0 <= pruning < 1
        ):
            raise ConfigurationError("Invalid MCL expansion/inflation/tolerance/pruning.")
        if len(data) > 2000:
            raise ConfigurationError("Dense native MCL is limited to 2000 nodes.")
        transition = data.copy()
        np.fill_diagonal(transition, 1.0)
        transition /= transition.sum(axis=0)
        for _ in range(iterations):
            updated = np.linalg.matrix_power(transition, expansion) ** inflation
            maxima = updated.argmax(axis=0)
            saved = updated[maxima, np.arange(len(data))].copy()
            updated[updated < pruning] = 0
            updated[maxima, np.arange(len(data))] = saved
            updated /= updated.sum(axis=0)
            converged = np.allclose(updated, transition, atol=tolerance, rtol=0)
            transition = updated
            if converged:
                groups = {
                    tuple(np.flatnonzero(row > tolerance).tolist())
                    for index, row in enumerate(transition)
                    if row[index] > tolerance
                }
                return tuple(sorted(group for group in groups if group))
        raise ConfigurationError("MCL did not converge within iterations; increase the budget.")
    methods = {
        "kmeans": "KMeans",
        "minibatch_kmeans": "MiniBatchKMeans",
        "dbscan": "DBSCAN",
        "spectral": "SpectralClustering",
        "affinity_propagation": "AffinityPropagation",
        "mean_shift": "MeanShift",
        "agglomerative": "AgglomerativeClustering",
        "birch": "Birch",
        "gmm": "GaussianMixture",
    }
    if method not in methods:
        raise ConfigurationError(f"Unknown clustering method {method!r}.")
    if method in {"kmeans", "minibatch_kmeans", "spectral", "affinity_propagation", "gmm"}:
        options.setdefault("random_state", seed)
    backend = require("sklearn.mixture" if method == "gmm" else "sklearn.cluster", "bio-analysis")
    try:
        return getattr(backend, methods[method])(**options).fit(data)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"Clustering failed: {exc}") from exc
