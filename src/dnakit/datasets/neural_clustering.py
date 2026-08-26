"""DNA foundation-model representations followed by k-means clustering."""

from __future__ import annotations

import importlib
from collections.abc import Iterable
from dataclasses import replace
from typing import Any

from dnakit.core import DNARecord, DNASet
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, SequenceError
from dnakit.representations import RepresentationBackend, extract_representations

from .config import NeuralClusteringConfig
from .results import Cluster, NeuralClusteringResult


def _materialize(
    records: Iterable[DNARecord],
    *,
    max_records: int,
) -> tuple[DNARecord, ...]:
    materialized: list[DNARecord] = []
    for index, record in enumerate(records):
        if index >= max_records:
            raise ConfigurationError(
                "Neural clustering exceeded max_records.",
                code="NEURAL_CLUSTER_RECORD_LIMIT",
                context={"max_records": max_records},
            )
        if not isinstance(record, DNARecord):
            raise SequenceError(
                "Neural clustering accepts only DNARecord objects.",
                code="INVALID_RECORD_SEQUENCE",
                context={"index": index, "type": type(record).__name__},
            )
        materialized.append(record)
    if not materialized:
        raise ConfigurationError(
            "At least one DNARecord is required for neural clustering.",
            code="EMPTY_NEURAL_CLUSTER_INPUT",
        )
    return tuple(materialized)


def _require_scientific_stack() -> tuple[Any, Any, Any, Any]:
    try:
        import numpy as np

        cluster_module = importlib.import_module("sklearn.cluster")
        decomposition_module = importlib.import_module("sklearn.decomposition")
        metrics_module = importlib.import_module("sklearn.metrics")
    except ImportError as exc:
        raise BackendUnavailableError(
            "Neural clustering requires NumPy and scikit-learn.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
        ) from exc
    return (
        np,
        cluster_module.KMeans,
        decomposition_module.PCA,
        metrics_module.silhouette_score,
    )


def _l2_normalize(matrix: Any, np: Any) -> Any:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    zero_rows = np.flatnonzero(norms[:, 0] == 0)
    if len(zero_rows):
        raise ConfigurationError(
            "Cannot L2-normalize zero-valued model representations.",
            code="ZERO_REPRESENTATION_VECTOR",
            context={"row_indices": tuple(int(value) for value in zero_rows[:20])},
        )
    return matrix / norms


def neural_cluster_sequences(
    records: Iterable[DNARecord],
    *,
    config: NeuralClusteringConfig | None = None,
    backend: RepresentationBackend | None = None,
) -> NeuralClusteringResult:
    """Extract model representations and group them with reproducible k-means."""

    resolved = NeuralClusteringConfig() if config is None else config
    if not isinstance(resolved, NeuralClusteringConfig):
        raise ConfigurationError(
            "config must be NeuralClusteringConfig or None.",
            code="INVALID_NEURAL_CLUSTER_CONFIG",
        )
    materialized = _materialize(records, max_records=resolved.max_records)
    if resolved.n_clusters > len(materialized):
        raise ConfigurationError(
            "n_clusters cannot exceed the number of input records.",
            code="INVALID_NEURAL_CLUSTER_COUNT",
            context={"n_clusters": resolved.n_clusters, "input_count": len(materialized)},
        )

    representation_config = replace(
        resolved.representation,
        max_records=min(resolved.representation.max_records, resolved.max_records),
    )
    representation_result = extract_representations(
        materialized,
        config=representation_config,
        backend=backend,
    )
    np, kmeans_class, pca_class, silhouette_function = _require_scientific_stack()
    matrix = np.asarray(representation_result.representations, dtype=np.float64)
    if resolved.normalize:
        matrix = _l2_normalize(matrix, np)

    explained_variance: tuple[float, ...] = ()
    pca_components = resolved.pca_components
    if pca_components is not None:
        if len(materialized) < 2:
            raise ConfigurationError(
                "PCA requires at least two input records.",
                code="INVALID_NEURAL_CLUSTER_PCA",
            )
        maximum = min(int(matrix.shape[0]), int(matrix.shape[1]))
        if pca_components > maximum:
            raise ConfigurationError(
                "pca_components cannot exceed min(input_count, embedding_dimension).",
                code="INVALID_NEURAL_CLUSTER_PCA",
                context={"pca_components": pca_components, "maximum": maximum},
            )
        pca = pca_class(n_components=pca_components, svd_solver="full")
        matrix = pca.fit_transform(matrix)
        explained_variance = tuple(float(value) for value in pca.explained_variance_ratio_)

    estimator = kmeans_class(
        n_clusters=resolved.n_clusters,
        init="k-means++",
        n_init=resolved.n_init,
        max_iter=resolved.max_iter,
        tol=resolved.tolerance,
        random_state=resolved.seed,
        algorithm="lloyd",
    )
    raw_labels = estimator.fit_predict(matrix)
    raw_groups = {
        cluster_index: tuple(int(index) for index in np.flatnonzero(raw_labels == cluster_index))
        for cluster_index in range(resolved.n_clusters)
    }
    if any(not members for members in raw_groups.values()):
        raise ConfigurationError(
            "k-means returned an empty cluster.",
            code="EMPTY_NEURAL_CLUSTER",
        )
    old_cluster_order = tuple(
        sorted(raw_groups, key=lambda cluster_index: raw_groups[cluster_index][0])
    )
    old_to_new = {old: new for new, old in enumerate(old_cluster_order)}
    labels = tuple(old_to_new[int(label)] for label in raw_labels)
    reordered_centers = estimator.cluster_centers_[list(old_cluster_order)]

    clusters: list[Cluster] = []
    representatives: list[DNARecord] = []
    for cluster_index, old_cluster_index in enumerate(old_cluster_order):
        members = raw_groups[old_cluster_index]
        center = reordered_centers[cluster_index]
        representative_index = min(
            members,
            key=lambda index: (float(np.linalg.norm(matrix[index] - center)), index),
        )
        representatives.append(materialized[representative_index])
        clusters.append(
            Cluster(
                cluster_index,
                members,
                tuple(materialized[index].id for index in members),
                representative_index,
                materialized[representative_index].id,
            )
        )

    score: float | None = None
    unique_labels = len(set(labels))
    if 1 < unique_labels < len(materialized):
        score_kwargs: dict[str, object] = {"metric": "euclidean"}
        if len(materialized) > resolved.silhouette_sample_size:
            score_kwargs.update(
                {
                    "sample_size": resolved.silhouette_sample_size,
                    "random_state": resolved.seed,
                }
            )
        score = float(silhouette_function(matrix, labels, **score_kwargs))

    centers = tuple(
        tuple(float(value) for value in center)
        for center in np.asarray(reordered_centers, dtype=np.float64)
    )
    return NeuralClusteringResult(
        tuple(clusters),
        labels,
        DNASet(representatives),
        centers,
        representation_result.model_name,
        representation_result.checkpoint_path,
        representation_result.pooling,
        representation_result.embedding_dimension,
        int(matrix.shape[1]),
        resolved.n_clusters,
        float(estimator.inertia_),
        score,
        resolved.normalize,
        pca_components,
        explained_variance,
        resolved.seed,
        resolved.n_init,
        resolved.max_iter,
        int(estimator.n_iter_),
        len(materialized),
    )


__all__ = ["neural_cluster_sequences"]
