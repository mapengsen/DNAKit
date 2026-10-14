"""Collection uniqueness, diversity, and redundancy evaluation."""

from __future__ import annotations

import multiprocessing
import threading
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from math import fsum

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core import DNARecord, DNASequence
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
    pair_count,
    pair_similarity,
    record_for,
    report,
    require_pairwise_compatible,
)
from .config import DiversityEvaluationConfig, EvaluationLimits, UniquenessEvaluationConfig
from .results import EvaluationReport


def _records(value: EvaluationInput, config: UniquenessEvaluationConfig) -> tuple[DNARecord, ...]:
    return tuple(record_for(item) for item in materialize_input(value, limits=config.limits))


def evaluate_uniqueness(
    value: EvaluationInput,
    *,
    config: UniquenessEvaluationConfig | None = None,
    plot: bool = False,
    verbose: bool = True,
) -> EvaluationReport:
    """Compute non-duplicate fraction; ``plot=True`` saves ``uniqueness_ratio.png`` in CWD."""

    if not isinstance(plot, bool):
        raise TypeError("plot must be a boolean.")
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
    result = report(
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
    if plot:
        from ._plots import save_uniqueness_chart

        save_uniqueness_chart(score, unique_count, total, equivalence=resolved.equivalence)
    return result


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
    from dnakit._progress import track

    return track(
        pairs,
        description="Diversity sequence pairs",
        total=pair_count(record_count),
    )


_WORKER_SEQUENCES: tuple[DNASequence, ...] = ()
_WORKER_MAX_CELLS = 0


def _init_diversity_worker(sequences: tuple[str, ...], max_cells: int) -> None:
    global _WORKER_SEQUENCES, _WORKER_MAX_CELLS
    _WORKER_SEQUENCES = tuple(DNASequence(sequence, alphabet="iupac") for sequence in sequences)
    _WORKER_MAX_CELLS = max_cells


def _diversity_worker_distance(pair: tuple[int, int]) -> float:
    left, right = pair
    return edit_distance(
        _WORKER_SEQUENCES[left],
        _WORKER_SEQUENCES[right],
        max_cells=_WORKER_MAX_CELLS,
        verbose=False,
    ).distance


def _diversity_pair_distances(
    records: tuple[DNARecord, ...], config: DiversityEvaluationConfig
) -> Iterable[float]:
    comparisons = pair_count(len(records))
    if config.workers == 1 or comparisons < 2:
        for left, right in _pair_indices(len(records), show_progress=config.show_progress):
            yield edit_distance(
                records[left], records[right], max_cells=config.limits.max_alignment_cells
            ).distance
        return

    longest = sorted(records, key=lambda record: record.sequence.symbol_length, reverse=True)[:2]
    if (longest[0].sequence.symbol_length + 1) * (
        longest[1].sequence.symbol_length + 1
    ) > config.limits.max_alignment_cells:
        edit_distance(
            longest[0], longest[1], max_cells=config.limits.max_alignment_cells, verbose=False
        )

    sequences = tuple(record.sequence.symbols for record in records)
    # Forking a multi-threaded process (e.g. the MCP server, which runs tools in a worker
    # thread) can deadlock the children on inherited locks, so spawn there. Single-threaded
    # callers keep fork, which also works from scripts without a __main__ guard.
    context = "spawn" if threading.active_count() > 1 else None
    with ProcessPoolExecutor(
        max_workers=min(config.workers, comparisons),
        mp_context=multiprocessing.get_context(context),
        initializer=_init_diversity_worker,
        initargs=(sequences, config.limits.max_alignment_cells),
    ) as executor:
        distances = executor.map(
            _diversity_worker_distance,
            _pair_indices(len(records), show_progress=False),
            chunksize=16,
        )
        if config.show_progress:
            from dnakit._progress import track

            distances = track(
                distances,
                description="Diversity sequence pairs",
                total=pair_count(len(records)),
            )
        yield from distances


def _evaluate_levenshtein_diversity(
    records: tuple[DNARecord, ...],
    config: DiversityEvaluationConfig,
    *,
    plot: bool,
) -> EvaluationReport:
    comparisons = pair_count(len(records))
    enforce_pair_limit(comparisons, config.limits)
    for record in records:
        require_pairwise_compatible(record.sequence, role="diversity")
    pair_distances: list[float] | None = [] if plot else None
    total_distance = 0.0
    for distance in _diversity_pair_distances(records, config):
        total_distance += distance
        if pair_distances is not None:
            pair_distances.append(distance)
    mean_distance = total_distance / comparisons if comparisons else None
    result = report(
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
            "workers": config.workers,
            "limits": config.limits,
            "reference": "https://doi.org/10.1016/j.compbiomed.2024.109440",
        },
        metrics={
            "score": mean_distance,
            "diversity": mean_distance,
            "record_count": len(records),
            "mean_pair_distance": mean_distance,
            "mean_pairwise_levenshtein_distance": mean_distance,
            "pairwise_comparison_count": comparisons,
        },
    )
    if pair_distances is not None:
        from ._plots import save_diversity_charts

        save_diversity_charts(mean_distance, pair_distances)
    return result


def evaluate_diversity(
    value: EvaluationInput,
    *,
    config: DiversityEvaluationConfig | None = None,
    plot: bool = False,
    verbose: bool = True,
) -> EvaluationReport:
    """Evaluate mean pairwise Levenshtein distance; ``plot=True`` saves two PNGs in CWD."""

    if not isinstance(plot, bool):
        raise TypeError("plot must be a boolean.")
    resolved = DiversityEvaluationConfig() if config is None else config
    if not isinstance(resolved, DiversityEvaluationConfig):
        raise TypeError("config must be DiversityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    records = tuple(record_for(item) for item in items)
    return _evaluate_levenshtein_diversity(records, resolved, plot=plot)


def evaluate_kmer_diversity(
    value: EvaluationInput,
    *,
    k: int = 5,
    canonical: bool = True,
    limits: EvaluationLimits | None = None,
    verbose: bool = True,
) -> EvaluationReport:
    """Report mean pairwise and nearest-neighbor canonical k-mer Jaccard distances."""
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer.")
    if not isinstance(canonical, bool):
        raise TypeError("canonical must be boolean.")
    resolved_limits = EvaluationLimits() if limits is None else limits
    if not isinstance(resolved_limits, EvaluationLimits):
        raise TypeError("limits must be EvaluationLimits or None.")
    records = tuple(record_for(item) for item in materialize_input(value, limits=resolved_limits))
    comparisons = pair_count(len(records))
    enforce_pair_limit(comparisons, resolved_limits)
    nearest = [1.0] * len(records)
    distances = []
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            similarity, _ = pair_similarity(
                records[left],
                records[right],
                method="kmer",
                k=k,
                canonical=canonical,
                max_alignment_cells=resolved_limits.max_alignment_cells,
            )
            distance = 1.0 - similarity
            distances.append(distance)
            nearest[left] = min(nearest[left], distance)
            nearest[right] = min(nearest[right], distance)
    mean_pair = fsum(distances) / comparisons if comparisons else None
    mean_nearest = fsum(nearest) / len(records) if comparisons else None
    return report(
        name="kmer_diversity",
        method="canonical-kmer-jaccard-distance" if canonical else "kmer-jaccard-distance",
        version="eval-kmer-diversity-v1",
        parameters={"k": k, "canonical": canonical, "limits": resolved_limits},
        metrics={
            "score": mean_pair,
            "record_count": len(records),
            "pairwise_comparison_count": comparisons,
            "mean_pair_distance": mean_pair,
            "mean_nearest_neighbor_distance": mean_nearest,
        },
    )


def evaluate_redundancy(
    value: EvaluationInput, *, config: DiversityEvaluationConfig | None = None, verbose: bool = True
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


__all__ = [
    "evaluate_diversity",
    "evaluate_kmer_diversity",
    "evaluate_redundancy",
    "evaluate_uniqueness",
]

_register_public_verbose(globals())
