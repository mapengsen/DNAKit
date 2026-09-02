"""Collection uniqueness, diversity, and redundancy evaluation."""

from __future__ import annotations

from collections.abc import Iterable

from dnakit.core import DNARecord
from dnakit.datasets import (
    ClusterConfig,
    IUPACDeduplicationConfig,
    cluster_sequences,
    deduplicate,
    deduplicate_iupac,
)
from dnakit.similarity import edit_distance

from ._shared import (
    EvaluationInput,
    enforce_pair_limit,
    materialize_input,
    mean,
    pair_count,
    pair_similarity,
    record_for,
    report,
    require_pairwise_compatible,
)
from .config import DiversityEvaluationConfig, UniquenessEvaluationConfig
from .results import EvaluationReport


def _records(value: EvaluationInput, config: UniquenessEvaluationConfig) -> tuple[DNARecord, ...]:
    return tuple(record_for(item) for item in materialize_input(value, limits=config.limits))


def evaluate_uniqueness(
    value: EvaluationInput,
    *,
    config: UniquenessEvaluationConfig | None = None,
) -> EvaluationReport:
    """Compute non-duplicate fraction under one explicitly configured equivalence."""

    resolved = UniquenessEvaluationConfig() if config is None else config
    if not isinstance(resolved, UniquenessEvaluationConfig):
        raise TypeError("config must be UniquenessEvaluationConfig or None.")
    records = _records(value, resolved)
    if resolved.equivalence == "iupac":
        iupac_result = deduplicate_iupac(
            records,
            config=IUPACDeduplicationConfig(
                max_records=resolved.limits.max_records,
                max_pairwise_comparisons=resolved.limits.max_pairwise_comparisons,
            ),
        )
        groups = tuple(group.member_ids for group in iupac_result.groups)
        unique_count = iupac_result.output_count
        comparisons = iupac_result.pairwise_comparison_count
        strategy = iupac_result.grouping_strategy
    elif resolved.equivalence == "approximate":
        cluster_method = (
            "identity" if resolved.approximate_method == "exact" else resolved.approximate_method
        )
        cluster_result = cluster_sequences(
            records,
            config=ClusterConfig(
                method=cluster_method,
                threshold=1.0 if resolved.approximate_method == "exact" else resolved.threshold,
                k=resolved.k,
                canonical=resolved.canonical,
                max_records=resolved.limits.max_records,
                max_pairwise_comparisons=resolved.limits.max_pairwise_comparisons,
                max_alignment_cells=resolved.limits.max_alignment_cells,
            ),
        )
        groups = tuple(cluster.member_ids for cluster in cluster_result.clusters)
        unique_count = len(cluster_result.clusters)
        comparisons = cluster_result.pairwise_comparison_count
        strategy = cluster_result.clustering_strategy
    else:
        dedup_result = deduplicate(records, equivalence=resolved.equivalence)
        groups = tuple(group.member_ids for group in dedup_result.groups)
        unique_count = len(dedup_result.groups)
        comparisons = 0
        strategy = "hash-equivalence-groups"
    duplicate_groups = tuple(group for group in groups if len(group) > 1)
    total = len(records)
    score = unique_count / total if total else 1.0
    return report(
        name="uniqueness",
        method=f"{resolved.equivalence}-equivalence",
        version="eval-uniqueness-v1",
        parameters={
            "equivalence": resolved.equivalence,
            "approximate_method": resolved.approximate_method,
            "threshold": resolved.threshold,
            "k": resolved.k,
            "canonical": resolved.canonical,
            "grouping_strategy": strategy,
            "limits": resolved.limits,
        },
        metrics={
            "score": score,
            "uniqueness_score": score,
            "record_count": total,
            "unique_count": unique_count,
            "duplicate_record_count": total - unique_count,
            "duplicate_group_count": len(duplicate_groups),
            "groups": groups,
            "duplicate_groups": duplicate_groups,
            "pairwise_comparison_count": comparisons,
        },
    )


def _similarity_state(
    records: tuple[DNARecord, ...], config: DiversityEvaluationConfig
) -> tuple[tuple[tuple[float, ...], ...], int]:
    comparisons = pair_count(len(records))
    enforce_pair_limit(comparisons, config.limits)
    matrix = [
        [1.0 if left == right else 0.0 for right in range(len(records))]
        for left in range(len(records))
    ]
    for left, right in _pair_indices(
        len(records),
        show_progress=config.show_progress,
    ):
        similarity, _ = pair_similarity(
            records[left],
            records[right],
            method=config.method,
            k=config.k,
            canonical=config.canonical,
            max_alignment_cells=config.limits.max_alignment_cells,
        )
        matrix[left][right] = matrix[right][left] = similarity
    return tuple(tuple(row) for row in matrix), comparisons


def _component_count(matrix: tuple[tuple[float, ...], ...], threshold: float) -> int:
    parents = list(range(len(matrix)))

    def find(value: int) -> int:
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value

    for left in range(len(matrix)):
        for right in range(left + 1, len(matrix)):
            if matrix[left][right] >= threshold:
                left_root, right_root = find(left), find(right)
                if left_root != right_root:
                    parents[right_root] = left_root
    return len({find(index) for index in range(len(matrix))})


def _pair_indices(
    record_count: int,
    *,
    show_progress: bool,
) -> Iterable[tuple[int, int]]:
    pairs = (
        (left, right) for left in range(record_count) for right in range(left + 1, record_count)
    )
    if not show_progress:
        return pairs
    from rich.progress import track

    return track(
        pairs,
        description="Diversity sequence pairs",
        total=pair_count(record_count),
    )


def _evaluate_levenshtein_diversity(
    records: tuple[DNARecord, ...],
    config: DiversityEvaluationConfig,
) -> EvaluationReport:
    comparisons = pair_count(len(records))
    enforce_pair_limit(comparisons, config.limits)
    for record in records:
        require_pairwise_compatible(record.sequence, role="diversity")
    distances = tuple(
        float(
            edit_distance(
                records[left],
                records[right],
                max_cells=config.limits.max_alignment_cells,
            ).distance
        )
        for left, right in _pair_indices(
            len(records),
            show_progress=config.show_progress,
        )
    )
    mean_distance = mean(distances)
    return report(
        name="diversity",
        method="mean-pairwise-levenshtein-distance",
        version="eval-diversity-levenshtein-v1",
        parameters={
            "calculation": config.calculation,
            "formula": "sum_{i != j} Levenshtein(x_i,x_j) / (n*(n-1))",
            "distance_units": "edit operations",
            "normalization": "none",
            "sequence_preprocessing": "none; no padding or truncation",
            "ambiguity_policy": "literal IUPAC symbols",
            "topology_policy": "linearized at stored origin; no circular rotation",
            "undefined_policy": "score=None when fewer than two records",
            "show_progress": config.show_progress,
            "limits": config.limits,
            "reference": "https://doi.org/10.1016/j.compbiomed.2024.109440",
        },
        metrics={
            "score": mean_distance,
            "diversity": mean_distance,
            "record_count": len(records),
            "mean_pair_distance": mean_distance,
            "mean_pairwise_levenshtein_distance": mean_distance,
            "minimum_pairwise_levenshtein_distance": min(distances) if distances else None,
            "maximum_pairwise_levenshtein_distance": max(distances) if distances else None,
            "pairwise_comparison_count": comparisons,
        },
    )


def evaluate_diversity(
    value: EvaluationInput,
    *,
    config: DiversityEvaluationConfig | None = None,
) -> EvaluationReport:
    """Evaluate diversity with similarity summaries or mean pairwise Levenshtein distance."""

    resolved = DiversityEvaluationConfig() if config is None else config
    if not isinstance(resolved, DiversityEvaluationConfig):
        raise TypeError("config must be DiversityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    records = tuple(record_for(item) for item in items)
    if resolved.calculation == "levenshtein":
        return _evaluate_levenshtein_diversity(records, resolved)
    matrix, comparisons = _similarity_state(records, resolved)
    distances = tuple(
        1.0 - matrix[left][right]
        for left in range(len(records))
        for right in range(left + 1, len(records))
    )
    nearest_distances = (
        tuple(
            min(1.0 - matrix[index][other] for other in range(len(records)) if other != index)
            for index in range(len(records))
        )
        if len(records) > 1
        else ()
    )
    cluster_count = _component_count(matrix, resolved.cluster_threshold)
    cluster_coverage = cluster_count / len(records) if records else 1.0
    mean_distance = mean(distances)
    mean_nearest = mean(nearest_distances)
    score = mean_nearest if mean_nearest is not None else (1.0 if len(records) <= 1 else 0.0)
    return report(
        name="diversity",
        method="exhaustive-pair-distance-and-threshold-components",
        version="eval-diversity-v1",
        parameters={
            "calculation": resolved.calculation,
            "similarity_method": resolved.method,
            "distance_definition": "1 - pair similarity",
            "k": resolved.k,
            "canonical": resolved.canonical,
            "cluster_threshold": resolved.cluster_threshold,
            "cluster_strategy": "connected-components-of-threshold-graph",
            "short_k_policy": "literal equality when either sequence is shorter than k",
            "topology_policy": "linearized at stored origin; no circular rotation",
            "singleton_policy": "score=1; pair metrics undefined",
            "show_progress": resolved.show_progress,
            "limits": resolved.limits,
        },
        metrics={
            "score": score,
            "record_count": len(records),
            "mean_pair_distance": mean_distance,
            "mean_nearest_neighbor_distance": mean_nearest,
            "minimum_nearest_neighbor_distance": min(nearest_distances)
            if nearest_distances
            else None,
            "cluster_count": cluster_count,
            "cluster_coverage": cluster_coverage,
            "pairwise_comparison_count": comparisons,
        },
    )


def evaluate_redundancy(
    value: EvaluationInput,
    *,
    config: DiversityEvaluationConfig | None = None,
) -> EvaluationReport:
    """Report exact, similarity-threshold, and cluster-level redundancy."""

    resolved = DiversityEvaluationConfig() if config is None else config
    if not isinstance(resolved, DiversityEvaluationConfig):
        raise TypeError("config must be DiversityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    records = tuple(record_for(item) for item in items)
    matrix, comparisons = _similarity_state(records, resolved)
    exact_unique_count = len(deduplicate(records, equivalence="exact").groups)
    near_pair_count = sum(
        matrix[left][right] >= resolved.cluster_threshold
        for left in range(len(records))
        for right in range(left + 1, len(records))
    )
    total_pairs = pair_count(len(records))
    cluster_count = _component_count(matrix, resolved.cluster_threshold)
    exact_fraction = (len(records) - exact_unique_count) / len(records) if records else 0.0
    cluster_fraction = (len(records) - cluster_count) / len(records) if records else 0.0
    near_pair_fraction = near_pair_count / total_pairs if total_pairs else 0.0
    score = (exact_fraction + near_pair_fraction + cluster_fraction) / 3.0
    return report(
        name="redundancy",
        method="exact-near-pair-threshold-component-mean",
        version="eval-redundancy-v1",
        parameters={
            "similarity_method": resolved.method,
            "threshold": resolved.cluster_threshold,
            "k": resolved.k,
            "canonical": resolved.canonical,
            "score_formula": (
                "mean(exact_duplicate_fraction,near_pair_fraction,cluster_redundancy_fraction)"
            ),
            "short_k_policy": "literal equality when either sequence is shorter than k",
            "topology_policy": "linearized at stored origin; no circular rotation",
            "limits": resolved.limits,
        },
        metrics={
            "score": score,
            "record_count": len(records),
            "exact_unique_count": exact_unique_count,
            "exact_duplicate_fraction": exact_fraction,
            "near_pair_count": near_pair_count,
            "pair_count": total_pairs,
            "near_pair_fraction": near_pair_fraction,
            "cluster_count": cluster_count,
            "cluster_redundancy_fraction": cluster_fraction,
            "pairwise_comparison_count": comparisons,
        },
    )


__all__ = ["evaluate_diversity", "evaluate_redundancy", "evaluate_uniqueness"]
