"""DNA adaptations of the MOSES Frag and SNN generative-set metrics."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping
from typing import TypeVar, cast

from dnakit.core import DNARecord
from dnakit.core._json import FrozenDict
from dnakit.descriptors import kmer_statistics
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints import hashed_kmer_fingerprint

from ._shared import (
    EvaluationInput,
    enforce_pair_limit,
    materialize_input,
    pair_count,
    record_for,
    report,
    require_nonempty,
)
from .config import EvaluationLimits, FragmentSimilarityConfig, SNNConfig
from .results import EvaluationEntry, EvaluationReport

_T = TypeVar("_T")


def _tracked(
    items: tuple[_T, ...],
    *,
    enabled: bool,
    description: str,
) -> Iterable[_T]:
    if not enabled:
        return items
    from rich.progress import track

    return cast(Iterable[_T], track(items, description=description, total=len(items)))


def _materialize_collections(
    generated: EvaluationInput,
    reference: EvaluationInput,
    *,
    operation: str,
    limits: EvaluationLimits,
) -> tuple[tuple[DNARecord, ...], tuple[DNARecord, ...]]:
    generated_items = materialize_input(generated, limits=limits)
    reference_items = materialize_input(reference, limits=limits)
    require_nonempty(generated_items, operation)
    require_nonempty(reference_items, operation)
    total_symbols = sum(item.sequence.symbol_length for item in generated_items + reference_items)
    if total_symbols > limits.max_total_symbols:
        raise ConfigurationError(
            f"Combined {operation} inputs exceed max_total_symbols.",
            code="EVALUATION_SYMBOL_LIMIT",
            context={
                "total_symbols": total_symbols,
                "max_total_symbols": limits.max_total_symbols,
            },
        )
    return (
        tuple(record_for(item) for item in generated_items),
        tuple(record_for(item) for item in reference_items),
    )


def _fragment_counts(
    records: tuple[DNARecord, ...],
    *,
    config: FragmentSimilarityConfig,
    description: str,
    initial_observations: int,
) -> tuple[Counter[str], int, int]:
    counts: Counter[str] = Counter()
    observations = initial_observations
    ignored_ambiguities = 0
    for record in _tracked(
        records,
        enabled=config.show_progress,
        description=description,
    ):
        result = kmer_statistics(
            record,
            config.k,
            canonical=config.canonical,
            ambiguity_policy=config.ambiguity_policy,
            cross_gaps=False,
        )
        observations += result.denominator
        if observations > config.max_kmer_observations:
            raise ConfigurationError(
                "Fragment similarity exceeds max_kmer_observations.",
                code="FRAGMENT_KMER_LIMIT",
                context={
                    "observation_count": observations,
                    "max_kmer_observations": config.max_kmer_observations,
                },
            )
        counts.update(cast(Mapping[str, int], result.counts))
        ignored_ambiguities += result.ignored_ambiguity_count
    return counts, observations, ignored_ambiguities


def _cosine_counts(left: Counter[str], right: Counter[str]) -> tuple[float, float, float, float]:
    left_squared = math.fsum(value * value for value in left.values())
    right_squared = math.fsum(value * value for value in right.values())
    if left_squared == 0 or right_squared == 0:
        raise ConfigurationError(
            "Fragment similarity requires at least one eligible k-mer in each collection.",
            code="EMPTY_FRAGMENT_VECTOR",
            context={"generated_empty": left_squared == 0, "reference_empty": right_squared == 0},
        )
    shared = set(left) & set(right)
    dot = math.fsum(left[word] * right[word] for word in shared)
    denominator = math.sqrt(left_squared * right_squared)
    similarity = min(1.0, max(0.0, dot / denominator))
    return similarity, dot, left_squared, right_squared


def evaluate_fragment_similarity(
    generated: EvaluationInput,
    reference: EvaluationInput,
    *,
    config: FragmentSimilarityConfig | None = None,
) -> EvaluationReport:
    """Compare exact DNA k-mer count distributions with cosine similarity."""

    resolved = FragmentSimilarityConfig() if config is None else config
    if not isinstance(resolved, FragmentSimilarityConfig):
        raise ConfigurationError(
            "config must be FragmentSimilarityConfig or None.",
            code="INVALID_FRAGMENT_SIMILARITY_CONFIG",
        )
    generated_records, reference_records = _materialize_collections(
        generated,
        reference,
        operation="fragment similarity",
        limits=resolved.limits,
    )
    generated_counts, generated_observations, generated_ignored = _fragment_counts(
        generated_records,
        config=resolved,
        description="Frag generated",
        initial_observations=0,
    )
    reference_counts, total_observations, reference_ignored = _fragment_counts(
        reference_records,
        config=resolved,
        description="Frag reference",
        initial_observations=generated_observations,
    )
    similarity, dot, generated_squared, reference_squared = _cosine_counts(
        generated_counts,
        reference_counts,
    )
    shared_fragments = len(set(generated_counts) & set(reference_counts))
    union_fragments = len(set(generated_counts) | set(reference_counts))
    return report(
        name="fragment_similarity",
        method="cosine-similarity-of-exact-dna-kmer-counts",
        version="eval-frag-dna-v1",
        parameters={
            "formula": "sum_f(c_G(f)*c_R(f)) / sqrt(sum_f(c_G(f)^2)*sum_f(c_R(f)^2))",
            "molecular_analog": "MOSES Frag over BRICS fragment count distributions",
            "dna_fragment_definition": "overlapping fixed-length A/C/G/T k-mers",
            "k": resolved.k,
            "canonical": resolved.canonical,
            "overlapping": True,
            "ambiguity_policy": resolved.ambiguity_policy,
            "gap_policy": "k-mers never cross explicit Gap objects",
            "topology_policy": "linearized at stored origin; no circular-origin wrap",
            "higher_is_better": True,
            "bounds": (0.0, 1.0),
            "max_kmer_observations": resolved.max_kmer_observations,
            "show_progress": resolved.show_progress,
            "limits": resolved.limits,
            "inference": (
                "DNA k-mer adaptation; not BRICS Frag and not an experimentally "
                "validated DNA quality score"
            ),
        },
        metrics={
            "score": similarity,
            "frag": similarity,
            "fragment_similarity": similarity,
            "dot_product": dot,
            "generated_squared_norm": generated_squared,
            "reference_squared_norm": reference_squared,
            "generated_count": len(generated_records),
            "reference_count": len(reference_records),
            "generated_fragment_observations": generated_observations,
            "reference_fragment_observations": total_observations - generated_observations,
            "generated_unique_fragments": len(generated_counts),
            "reference_unique_fragments": len(reference_counts),
            "shared_unique_fragments": shared_fragments,
            "union_unique_fragments": union_fragments,
            "generated_ignored_ambiguity_count": generated_ignored,
            "reference_ignored_ambiguity_count": reference_ignored,
        },
    )


def _bit_fingerprints(
    records: tuple[DNARecord, ...],
    *,
    config: SNNConfig,
    description: str,
) -> tuple[tuple[frozenset[str], ...], tuple[int, ...], int]:
    fingerprints: list[frozenset[str]] = []
    observation_counts: list[int] = []
    ignored_ambiguities = 0
    for record in _tracked(
        records,
        enabled=config.show_progress,
        description=description,
    ):
        fingerprint = hashed_kmer_fingerprint(
            record,
            k=config.k,
            n_bits=config.n_bits,
            canonical=config.canonical,
            seed=config.seed,
            representation="sparse",
            ambiguity_policy=config.ambiguity_policy,
            cross_gaps=False,
        )
        fingerprints.append(frozenset(fingerprint.sparse_values()))
        observation_counts.append(fingerprint.observation_count)
        ignored_ambiguities += fingerprint.ignored_ambiguity_count
    return tuple(fingerprints), tuple(observation_counts), ignored_ambiguities


def _binary_tanimoto(left: frozenset[str], right: frozenset[str]) -> float:
    union_size = len(left | right)
    if union_size == 0:
        return 1.0
    return len(left & right) / union_size


def evaluate_snn(
    generated: EvaluationInput,
    reference: EvaluationInput,
    *,
    config: SNNConfig | None = None,
) -> EvaluationReport:
    """Average each generated sequence's nearest-reference Tanimoto similarity."""

    resolved = SNNConfig() if config is None else config
    if not isinstance(resolved, SNNConfig):
        raise ConfigurationError(
            "config must be SNNConfig or None.",
            code="INVALID_SNN_CONFIG",
        )
    generated_records, reference_records = _materialize_collections(
        generated,
        reference,
        operation="SNN",
        limits=resolved.limits,
    )
    comparisons = pair_count(len(generated_records), len(reference_records))
    enforce_pair_limit(comparisons, resolved.limits)
    fingerprint_elements = (len(generated_records) + len(reference_records)) * resolved.n_bits
    if fingerprint_elements > resolved.max_fingerprint_elements:
        raise ConfigurationError(
            "SNN exceeds max_fingerprint_elements.",
            code="SNN_FINGERPRINT_LIMIT",
            context={
                "fingerprint_elements": fingerprint_elements,
                "max_fingerprint_elements": resolved.max_fingerprint_elements,
            },
        )
    generated_fingerprints, generated_observations, generated_ignored = _bit_fingerprints(
        generated_records,
        config=resolved,
        description="SNN generated fingerprints",
    )
    reference_fingerprints, reference_observations, reference_ignored = _bit_fingerprints(
        reference_records,
        config=resolved,
        description="SNN reference fingerprints",
    )

    entries: list[EvaluationEntry] = []
    nearest_values: list[float] = []
    generated_pairs = tuple(zip(generated_records, generated_fingerprints, strict=True))
    for generated_index, (record, fingerprint) in enumerate(
        _tracked(
            generated_pairs,
            enabled=resolved.show_progress,
            description="SNN nearest neighbors",
        )
    ):
        best_index = 0
        best_similarity = -1.0
        for reference_index, reference_fingerprint in enumerate(reference_fingerprints):
            similarity = _binary_tanimoto(fingerprint, reference_fingerprint)
            if similarity > best_similarity:
                best_index = reference_index
                best_similarity = similarity
        nearest_values.append(best_similarity)
        nearest_fingerprint = reference_fingerprints[best_index]
        shared_bit_count = len(fingerprint & nearest_fingerprint)
        union_bit_count = len(fingerprint | nearest_fingerprint)
        entries.append(
            EvaluationEntry(
                record.id,
                generated_index,
                FrozenDict(
                    {
                        "nearest_reference_id": reference_records[best_index].id,
                        "nearest_reference_index": best_index,
                        "nearest_similarity": best_similarity,
                        "generated_set_bit_count": len(fingerprint),
                        "nearest_reference_set_bit_count": len(nearest_fingerprint),
                        "nearest_shared_bit_count": shared_bit_count,
                        "nearest_union_bit_count": union_bit_count,
                    }
                ),
            )
        )
    snn = math.fsum(nearest_values) / len(nearest_values)
    return report(
        name="similarity_to_nearest_neighbor",
        method="mean-nearest-binary-tanimoto-over-hashed-dna-kmer-fingerprints",
        version="eval-snn-dna-v1",
        parameters={
            "formula": "mean_g max_r Tanimoto(fp(g), fp(r))",
            "molecular_analog": (
                "MOSES SNN over radius-2 Morgan fingerprints (1024 bits in the paper)"
            ),
            "dna_fingerprint": "SHA-256 hashed binary k-mer fingerprint",
            "k": resolved.k,
            "n_bits": resolved.n_bits,
            "canonical": resolved.canonical,
            "overlapping": True,
            "seed": resolved.seed,
            "ambiguity_policy": resolved.ambiguity_policy,
            "gap_policy": "fingerprint k-mers never cross explicit Gap objects",
            "topology_policy": "linearized at stored origin; no circular-origin wrap",
            "zero_vector_policy": "two empty fingerprints are identical; one empty is dissimilar",
            "tie_break": "reference-input-index-ascending",
            "higher_is_better": True,
            "bounds": (0.0, 1.0),
            "max_fingerprint_elements": resolved.max_fingerprint_elements,
            "show_progress": resolved.show_progress,
            "limits": resolved.limits,
            "inference": (
                "DNA k-mer fingerprint adaptation; not Morgan-fingerprint SNN and not an "
                "experimentally validated DNA quality score"
            ),
        },
        metrics={
            "score": snn,
            "snn": snn,
            "mean_nearest_similarity": snn,
            "minimum_nearest_similarity": min(nearest_values),
            "maximum_nearest_similarity": max(nearest_values),
            "generated_count": len(generated_records),
            "reference_count": len(reference_records),
            "pairwise_comparison_count": comparisons,
            "fingerprint_elements": fingerprint_elements,
            "generated_zero_fingerprint_count": sum(
                not fingerprint for fingerprint in generated_fingerprints
            ),
            "reference_zero_fingerprint_count": sum(
                not fingerprint for fingerprint in reference_fingerprints
            ),
            "generated_kmer_observation_count": sum(generated_observations),
            "reference_kmer_observation_count": sum(reference_observations),
            "generated_ignored_ambiguity_count": generated_ignored,
            "reference_ignored_ambiguity_count": reference_ignored,
        },
        entries=entries,
    )


__all__ = ["evaluate_fragment_similarity", "evaluate_snn"]
