"""Label-stratified diversity sampling from frozen DNA model representations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from typing import Any, Literal, cast

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core import DNARecord, DNASet
from dnakit.core._json import JSONScalar
from dnakit.exceptions import ConfigurationError
from dnakit.representations import RepresentationBackend, extract_representations

from .config import DiversitySamplingConfig
from .neural_clustering import _l2_normalize, _materialize, _require_scientific_stack
from .results import Cluster, DiversitySamplingResult, SamplingAllocation


def _strata(
    records: tuple[DNARecord, ...], label_key: str | None
) -> tuple[tuple[JSONScalar, tuple[int, ...]], ...]:
    groups: dict[tuple[type, JSONScalar], list[int]] = {}
    for index, record in enumerate(records):
        label: JSONScalar = None
        if label_key is not None:
            value = record.metadata.get(label_key)
            if value is None or not isinstance(value, (str, int, float, bool)):
                raise ConfigurationError(
                    "Each record must have a non-null scalar class label in metadata.",
                    code="INVALID_DIVERSITY_SAMPLING_LABEL",
                    context={"index": index, "record_id": record.id, "label_key": label_key},
                )
            label = cast(JSONScalar, value)
        # Keep e.g. numeric 1, boolean True, and text "1" as distinct classes.
        groups.setdefault((type(label), label), []).append(index)
    return tuple((key[1], tuple(indices)) for key, indices in groups.items())


def _allocate(counts: tuple[int, ...], n_samples: int) -> tuple[int, ...]:
    """Hamilton allocation with exact integer arithmetic and input-order tie breaking."""
    total = sum(counts)
    quotas = [n_samples * count // total for count in counts]
    order = sorted(range(len(counts)), key=lambda i: (-(n_samples * counts[i] % total), i))
    for index in order[: n_samples - sum(quotas)]:
        quotas[index] += 1
    return tuple(quotas)


def _compact_neighborhood(
    matrix: Any,
    centers: Any,
    indices: tuple[int, ...],
    n_samples: int,
    np: Any,
) -> tuple[int, ...]:
    """Choose the tightest fixed-size neighborhood among the fitted centers.

    Each center proposes n_samples nearest records from the entire class. Score
    each proposal by its mean Euclidean distance to its own mean embedding.
    This can span cluster boundaries when a cluster alone has too few records.
    """
    candidates: list[tuple[float, tuple[int, ...]]] = []
    for center in centers:
        distances = np.linalg.norm(matrix - center, axis=1)
        nearest = np.argsort(distances, kind="stable")[:n_samples]
        subset = matrix[nearest]
        dispersion = float(np.linalg.norm(subset - subset.mean(axis=0), axis=1).mean())
        selected = tuple(sorted(indices[int(i)] for i in nearest))
        candidates.append((dispersion, selected))
    return min(candidates)[1]


def diversity_sample(
    records: Iterable[DNARecord],
    *,
    n_samples: int | None = None,
    diversity: Literal["high", "low"] | None = None,
    low_n_clusters: int | None = None,
    config: DiversitySamplingConfig | None = None,
    backend: RepresentationBackend | None = None,
    verbose: bool = True,
) -> DiversitySamplingResult:
    """Select exactly n_samples records with high or low embedding diversity.

    Explicit n_samples, diversity, and low_n_clusters override their configuration
    fields without modifying config (defaults: 400, 'high', and 5 respectively).

    Defaults to frozen DNABERT2 mean embeddings and L2 normalization. Labels come
    from metadata[label_key]; set label_key=None to sample a single unlabeled pool.
    Both modes preserve class proportions using largest-remainder allocation;
    a small class can receive zero slots. Selected records retain input order and
    metadata, and distance ties use original input positions.

    'high' sets the total cluster count to n_samples and selects the nearest actual
    record to each cluster center. 'low' fits low_n_clusters per sampled class
    (capped by its distinct embedding count). Each center proposes that class's
    quota of nearest records; choose the proposal with the smallest mean distance
    to its own mean embedding. This is a compactness heuristic, not a global optimum.

    High mode rejects insufficient distinct embeddings; low mode allows identical
    embeddings from different input rows. Empty fitted clusters raise an error.
    Neither mode repeats input rows or deduplicates sequences across labels.
    """
    resolved = DiversitySamplingConfig() if config is None else config
    if not isinstance(resolved, DiversitySamplingConfig):
        raise ConfigurationError(
            "config must be DiversitySamplingConfig or None.",
            code="INVALID_DIVERSITY_SAMPLING_CONFIG",
        )
    if n_samples is not None:
        resolved = replace(resolved, n_samples=n_samples)
    if diversity is not None:
        resolved = replace(resolved, diversity=diversity)
    if low_n_clusters is not None:
        resolved = replace(resolved, low_n_clusters=low_n_clusters)
    materialized = _materialize(
        records, max_records=min(resolved.max_records, resolved.representation.max_records)
    )
    if resolved.n_samples > len(materialized):
        raise ConfigurationError(
            "n_samples cannot exceed the number of input records.",
            code="INVALID_DIVERSITY_SAMPLE_SIZE",
            context={"n_samples": resolved.n_samples, "input_count": len(materialized)},
        )
    strata = _strata(materialized, resolved.label_key)
    quotas = _allocate(tuple(len(indices) for _, indices in strata), resolved.n_samples)
    # Fail on absent clustering dependencies before loading a foundation model.
    np, kmeans_class, _, _ = _require_scientific_stack()
    representation = extract_representations(
        materialized, config=resolved.representation, backend=backend, verbose=False
    )
    matrix = np.asarray(representation.representations, dtype=np.float64)
    if resolved.normalize:
        matrix = _l2_normalize(matrix, np)

    clusters: list[Cluster] = []
    centers: list[tuple[float, ...]] = []
    cluster_labels: list[int | None] = [None] * len(materialized)
    allocations: list[SamplingAllocation] = []
    selected: list[int] = []
    for (label, indices), quota in zip(strata, quotas, strict=True):
        allocations.append(SamplingAllocation(label, len(indices), quota))
        if quota == 0:
            continue
        class_matrix = matrix[list(indices)]
        distinct_count = int(np.unique(class_matrix, axis=0).shape[0])
        n_clusters = (
            quota if resolved.diversity == "high" else min(resolved.low_n_clusters, distinct_count)
        )
        if distinct_count < n_clusters:
            raise ConfigurationError(
                "A class has fewer distinct embeddings than its allocated sample count.",
                code="INSUFFICIENT_DIVERSITY_EMBEDDINGS",
                context={"label": label, "n_samples": quota, "distinct_count": distinct_count},
                hint="Reduce n_samples or inspect duplicate sequences and embeddings.",
            )
        estimator = kmeans_class(
            n_clusters=n_clusters,
            init="k-means++",
            n_init=resolved.n_init,
            max_iter=resolved.max_iter,
            tol=resolved.tolerance,
            random_state=resolved.seed,
            algorithm="lloyd",
        )
        assignments = estimator.fit_predict(class_matrix)
        groups = [
            tuple(int(i) for i in np.flatnonzero(assignments == k)) for k in range(n_clusters)
        ]
        if any(not members for members in groups):
            raise ConfigurationError(
                "k-means returned fewer non-empty clusters than requested.",
                code="EMPTY_DIVERSITY_CLUSTER",
                context={"label": label, "n_samples": quota, "n_clusters": n_clusters},
            )
        for old_index in sorted(range(n_clusters), key=lambda k: groups[k][0]):
            members = tuple(indices[i] for i in groups[old_index])
            center = estimator.cluster_centers_[old_index]
            representative_index = min(
                members,
                key=lambda i: (float(np.linalg.norm(matrix[i] - center)), i),
            )
            cluster_index = len(clusters)
            for index in members:
                cluster_labels[index] = cluster_index
            clusters.append(
                Cluster(
                    cluster_index,
                    members,
                    tuple(materialized[i].id for i in members),
                    representative_index,
                    materialized[representative_index].id,
                )
            )
            centers.append(tuple(float(value) for value in center))
            if resolved.diversity == "high":
                selected.append(representative_index)
        if resolved.diversity == "low":
            selected.extend(
                _compact_neighborhood(class_matrix, estimator.cluster_centers_, indices, quota, np)
            )

    selected_indices = tuple(sorted(selected))
    return DiversitySamplingResult(
        records=DNASet(materialized[i] for i in selected_indices),
        selected_indices=selected_indices,
        allocations=tuple(allocations),
        clusters=tuple(clusters),
        cluster_labels=tuple(cluster_labels),
        centers=tuple(centers),
        input_count=len(materialized),
        n_samples=resolved.n_samples,
        label_key=resolved.label_key,
        model_name=representation.model_name,
        checkpoint_path=representation.checkpoint_path,
        pooling=representation.pooling,
        embedding_dimension=representation.embedding_dimension,
        normalize=resolved.normalize,
        seed=resolved.seed,
        n_init=resolved.n_init,
        max_iter=resolved.max_iter,
        tolerance=resolved.tolerance,
        diversity=resolved.diversity,
        low_n_clusters=resolved.low_n_clusters if resolved.diversity == "low" else None,
        sampling_strategy=(
            "one-per-cluster"
            if resolved.diversity == "high"
            else "minimum-dispersion-centroid-neighborhood"
        ),
    )


__all__ = ["diversity_sample"]

_register_public_verbose(globals())
