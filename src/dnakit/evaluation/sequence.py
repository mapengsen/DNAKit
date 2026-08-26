"""Single-sequence validity, ambiguity, complexity, and quality evaluation."""

from __future__ import annotations

import math
from collections import Counter

from dnakit.core import Gap, IssueSeverity
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.descriptors import (
    exact_repeat_fraction,
    homopolymer_runs,
    linguistic_complexity,
    shannon_entropy,
)
from dnakit.standardize import ValidationConfig, validate

from ._shared import (
    EvaluationInput,
    InputItem,
    aggregate_numeric,
    as_float,
    issue,
    materialize_input,
    report,
)
from .config import (
    AmbiguityEvaluationConfig,
    ComplexityEvaluationConfig,
    EvaluationLimits,
    QualityEvaluationConfig,
)
from .results import EvaluationEntry, EvaluationReport


def _validation_entry(item: InputItem, index: int, config: ValidationConfig) -> EvaluationEntry:
    result = validate(item.record or item.sequence, config=config)
    return EvaluationEntry(
        item.subject_id,
        index,
        FrozenDict(
            {
                "is_valid": result.is_valid,
                "symbol_length": result.symbol_length,
                "coordinate_span": result.coordinate_span,
                "ambiguity_count": result.ambiguity.total_count,
                "ambiguity_fraction": result.ambiguity.fraction,
                "quality": to_json_compatible(result.quality),
            }
        ),
        result.issues,
    )


def evaluate_validity(
    value: EvaluationInput,
    *,
    config: ValidationConfig | None = None,
    limits: EvaluationLimits | None = None,
) -> EvaluationReport:
    """Evaluate object-level alphabet, length, gap, and optional record constraints."""

    resolved = ValidationConfig() if config is None else config
    resolved_limits = EvaluationLimits() if limits is None else limits
    if not isinstance(resolved, ValidationConfig):
        raise TypeError("config must be ValidationConfig or None.")
    if not isinstance(resolved_limits, EvaluationLimits):
        raise TypeError("limits must be EvaluationLimits or None.")
    items = materialize_input(value, limits=resolved_limits)
    entries = tuple(_validation_entry(item, index, resolved) for index, item in enumerate(items))
    valid_count = sum(entry.metrics["is_valid"] is True for entry in entries)
    return report(
        name="sequence_validity",
        method="dnakit-standardize-validation",
        version="eval-validity-v1",
        parameters={"validation": resolved, "limits": resolved_limits},
        metrics={
            "record_count": len(entries),
            "valid_count": valid_count,
            "invalid_count": len(entries) - valid_count,
            "valid_fraction": valid_count / len(entries) if entries else None,
        },
        entries=entries,
    )


def _ambiguity_entry(
    item: InputItem,
    index: int,
    config: AmbiguityEvaluationConfig,
) -> EvaluationEntry:
    sequence = item.sequence
    positions = tuple(
        index for index, symbol in enumerate(sequence.symbols) if symbol not in "ACGT"
    )
    counts = Counter(sequence.symbols[position] for position in positions)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    if config.gap_denominator_policy == "error" and gaps:
        from dnakit.exceptions import UnsupportedGapOperationError

        raise UnsupportedGapOperationError(
            "Ambiguity evaluation is configured to reject Gap objects.",
            code="AMBIGUITY_GAP_REJECTED",
            context={"subject_id": item.subject_id},
        )
    denominator: int | None = sequence.symbol_length
    if config.gap_denominator_policy == "include_known":
        denominator = sequence.coordinate_span
    weighted_count = math.fsum(
        counts[symbol] * float(config.symbol_weights[symbol]) for symbol in sorted(counts)
    )
    fraction = len(positions) / denominator if denominator else (0.0 if denominator == 0 else None)
    weighted_fraction = (
        weighted_count / denominator if denominator else (0.0 if denominator == 0 else None)
    )
    entry_issues = []
    if denominator is None:
        entry_issues.append(
            issue(
                "EVAL_AMBIGUITY_DENOMINATOR_UNKNOWN",
                IssueSeverity.WARNING,
                "Ambiguity fraction is undefined because a gap length is unknown.",
                gap_denominator_policy=config.gap_denominator_policy,
            )
        )
    elif fraction is not None and fraction > config.max_fraction:
        entry_issues.append(
            issue(
                "EVAL_AMBIGUITY_HIGH",
                IssueSeverity.WARNING,
                "Ambiguity fraction exceeds the configured threshold.",
                fraction=fraction,
                maximum=config.max_fraction,
            )
        )
    return EvaluationEntry(
        item.subject_id,
        index,
        FrozenDict(
            {
                "ambiguity_count": len(positions),
                "ambiguity_fraction": fraction,
                "weighted_ambiguity_count": weighted_count,
                "weighted_ambiguity_fraction": weighted_fraction,
                "denominator": denominator,
                "positions": positions,
                "counts_by_symbol": {symbol: counts[symbol] for symbol in sorted(counts)},
                "gap_count": len(gaps),
                "unknown_gap_count": sum(gap.length is None for gap in gaps),
                "passes_threshold": fraction is not None and fraction <= config.max_fraction,
            }
        ),
        tuple(entry_issues),
    )


def evaluate_ambiguity(
    value: EvaluationInput,
    *,
    config: AmbiguityEvaluationConfig | None = None,
) -> EvaluationReport:
    """Report ambiguity symbols, positions, weighted impact, and explicit denominator."""

    resolved = AmbiguityEvaluationConfig() if config is None else config
    if not isinstance(resolved, AmbiguityEvaluationConfig):
        raise TypeError("config must be AmbiguityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    entries = tuple(_ambiguity_entry(item, index, resolved) for index, item in enumerate(items))
    return report(
        name="ambiguity",
        method="weighted-iupac-symbol-fraction",
        version="eval-ambiguity-v1",
        parameters={
            "max_fraction": resolved.max_fraction,
            "symbol_weights": resolved.symbol_weights,
            "gap_denominator_policy": resolved.gap_denominator_policy,
            "limits": resolved.limits,
            "coordinate_system": "0-based-symbol",
        },
        metrics={
            "record_count": len(entries),
            "ambiguity_fraction_summary": aggregate_numeric(entries, "ambiguity_fraction"),
            "weighted_ambiguity_fraction_summary": aggregate_numeric(
                entries, "weighted_ambiguity_fraction"
            ),
            "failing_count": sum(entry.metrics["passes_threshold"] is False for entry in entries),
        },
        entries=entries,
    )


def _complexity_entry(
    item: InputItem,
    index: int,
    config: ComplexityEvaluationConfig,
) -> EvaluationEntry:
    entropy = shannon_entropy(item.sequence, ambiguity_policy="ignore")
    linguistic = linguistic_complexity(
        item.sequence,
        max_word_size=config.max_word_size,
        ambiguity_policy="ignore",
        max_observations=config.max_observations_per_sequence,
    )
    repeat = exact_repeat_fraction(
        item.sequence,
        max_unit_length=config.max_repeat_unit,
        min_repeats=config.min_repeat_count,
        ambiguity_policy="ignore",
        max_comparisons=config.max_comparisons_per_sequence,
    )
    homopolymer = homopolymer_runs(item.sequence, ambiguity_policy="ignore")
    entropy_normalized = min(1.0, max(0.0, entropy.entropy / 2.0))
    homopolymer_cleanliness = (
        1.0
        if homopolymer.longest_length <= config.acceptable_homopolymer_length
        else config.acceptable_homopolymer_length / homopolymer.longest_length
    )
    components = {
        "entropy": entropy_normalized,
        "linguistic": linguistic.score,
        "repeat_cleanliness": 1.0 - repeat.repeat_fraction,
        "homopolymer_cleanliness": homopolymer_cleanliness,
    }
    denominator = math.fsum(float(config.weights[name]) for name in components)
    score = (
        math.fsum(components[name] * float(config.weights[name]) for name in components)
        / denominator
    )
    return EvaluationEntry(
        item.subject_id,
        index,
        FrozenDict(
            {
                "score": score,
                "components": components,
                "shannon_entropy_bits": entropy.entropy,
                "linguistic_complexity": linguistic.score,
                "repeat_fraction": repeat.repeat_fraction,
                "longest_homopolymer": homopolymer.longest_length,
                "ignored_ambiguity_count": item.sequence.ambiguity_count,
                "repeat_comparisons": repeat.comparisons,
                "complexity_observations": linguistic.observation_count,
            }
        ),
    )


def evaluate_complexity(
    value: EvaluationInput,
    *,
    config: ComplexityEvaluationConfig | None = None,
) -> EvaluationReport:
    """Combine normalized entropy, vocabulary, repeat, and homopolymer components."""

    resolved = ComplexityEvaluationConfig() if config is None else config
    if not isinstance(resolved, ComplexityEvaluationConfig):
        raise TypeError("config must be ComplexityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    entries = tuple(_complexity_entry(item, index, resolved) for index, item in enumerate(items))
    return report(
        name="complexity",
        method="weighted-normalized-transparent-components",
        version="eval-complexity-v1",
        parameters={
            "weights": resolved.weights,
            "entropy_normalization": "min(1,max(0,H_base2/2))",
            "repeat_cleanliness": "1-exact_tandem_repeat_union_fraction",
            "homopolymer_cleanliness": "1 if longest<=acceptable else acceptable/longest",
            "max_word_size": resolved.max_word_size,
            "max_repeat_unit": resolved.max_repeat_unit,
            "min_repeat_count": resolved.min_repeat_count,
            "acceptable_homopolymer_length": resolved.acceptable_homopolymer_length,
            "max_observations_per_sequence": resolved.max_observations_per_sequence,
            "max_comparisons_per_sequence": resolved.max_comparisons_per_sequence,
            "limits": resolved.limits,
        },
        metrics={
            "record_count": len(entries),
            "score_summary": aggregate_numeric(entries, "score"),
            "score": (
                math.fsum(as_float(entry.metrics, "score") for entry in entries) / len(entries)
                if entries
                else None
            ),
        },
        entries=entries,
    )


def evaluate_quality(
    value: EvaluationInput,
    *,
    config: QualityEvaluationConfig | None = None,
) -> EvaluationReport:
    """Combine validity, ambiguity, complexity, length, and record completeness."""

    resolved = QualityEvaluationConfig() if config is None else config
    if not isinstance(resolved, QualityEvaluationConfig):
        raise TypeError("config must be QualityEvaluationConfig or None.")
    items = materialize_input(value, limits=resolved.complexity.limits)
    entries: list[EvaluationEntry] = []
    validation_config = ValidationConfig(
        allow_empty=resolved.min_length == 0,
        min_length=resolved.min_length,
        max_length=resolved.max_length,
        max_ambiguity_fraction=None,
    )
    for index, item in enumerate(items):
        validity = _validation_entry(item, index, validation_config)
        ambiguity = _ambiguity_entry(item, index, resolved.ambiguity)
        complexity = _complexity_entry(item, index, resolved.complexity)
        length = item.sequence.symbol_length
        length_ok = length >= resolved.min_length and (
            resolved.max_length is None or length <= resolved.max_length
        )
        completeness = 1.0
        if item.record is not None and "phred_quality" in item.record.letter_annotations:
            quality_values = item.record.letter_annotations["phred_quality"]
            completeness = (
                sum(value >= 20 for value in quality_values) / len(quality_values)
                if quality_values
                else 1.0
            )
        ambiguity_value = ambiguity.metrics["weighted_ambiguity_fraction"]
        ambiguity_score = (
            0.0
            if ambiguity_value is None
            else 1.0 - min(1.0, as_float(ambiguity.metrics, "weighted_ambiguity_fraction"))
        )
        components = {
            "validity": float(validity.metrics["is_valid"] is True),
            "ambiguity": ambiguity_score,
            "complexity": as_float(complexity.metrics, "score"),
            "length": float(length_ok),
            "completeness": completeness,
        }
        weight_sum = math.fsum(float(resolved.weights[name]) for name in components)
        score = (
            math.fsum(value * float(resolved.weights[name]) for name, value in components.items())
            / weight_sum
        )
        status = (
            "pass"
            if score >= resolved.pass_score
            else ("warning" if score >= resolved.warning_score else "fail")
        )
        entries.append(
            EvaluationEntry(
                item.subject_id,
                index,
                FrozenDict(
                    {
                        "score": score,
                        "status": status,
                        "components": components,
                        "symbol_length": length,
                        "phred_q20_fraction": completeness if item.record is not None else None,
                    }
                ),
                validity.issues + ambiguity.issues,
            )
        )
    materialized = tuple(entries)
    return report(
        name="sequence_quality",
        method="weighted-validity-ambiguity-complexity-length-completeness",
        version="eval-quality-v1",
        parameters={
            "weights": resolved.weights,
            "min_length": resolved.min_length,
            "max_length": resolved.max_length,
            "warning_score": resolved.warning_score,
            "pass_score": resolved.pass_score,
            "phred_completeness_rule": "fraction of provided phred_quality values >=20",
            "missing_phred_policy": "neutral completeness=1",
            "undefined_ambiguity_policy": "ambiguity component is zero",
            "length_definition": "symbol_length excluding Gap spans",
            "ambiguity_config": resolved.ambiguity,
            "complexity_config": resolved.complexity,
        },
        metrics={
            "record_count": len(materialized),
            "score_summary": aggregate_numeric(materialized, "score"),
            "score": (
                math.fsum(as_float(entry.metrics, "score") for entry in materialized)
                / len(materialized)
                if materialized
                else None
            ),
            "pass_count": sum(entry.metrics["status"] == "pass" for entry in materialized),
            "warning_count": sum(entry.metrics["status"] == "warning" for entry in materialized),
            "fail_count": sum(entry.metrics["status"] == "fail" for entry in materialized),
        },
        entries=materialized,
    )


__all__ = [
    "evaluate_ambiguity",
    "evaluate_complexity",
    "evaluate_quality",
    "evaluate_validity",
]
