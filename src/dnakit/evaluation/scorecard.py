"""Configurable transparent aggregation of evaluation and metric results."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TypeAlias

from dnakit.core import MetricResult
from dnakit.exceptions import ConfigurationError

from ._shared import report
from .config import ScorecardConfig, ScoreRule
from .results import EvaluationReport

ScoreInput: TypeAlias = EvaluationReport | MetricResult[object] | int | float | None


@dataclass(frozen=True, slots=True)
class _ResolvedScore:
    name: str
    raw_value: float | None
    source: str
    metric: str | None
    source_method: str | None
    source_algorithm_version: str | None


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            "Scorecard values must be finite numbers.",
            code="INVALID_SCORECARD_VALUE",
            context={"component": name, "value": value},
        )
    return float(value)


def _resolve(name: str, value: object, rule: ScoreRule) -> _ResolvedScore:
    if value is None:
        return _ResolvedScore(name, None, "missing", rule.metric, None, None)
    if isinstance(value, EvaluationReport):
        metric_name = rule.metric or "score"
        if metric_name not in value.metrics:
            return _ResolvedScore(
                name,
                None,
                f"EvaluationReport:{value.name}",
                metric_name,
                value.method,
                value.algorithm_version,
            )
        return _ResolvedScore(
            name,
            _number(value.metrics[metric_name], name),
            f"EvaluationReport:{value.name}",
            metric_name,
            value.method,
            value.algorithm_version,
        )
    if isinstance(value, MetricResult):
        if rule.metric is not None:
            raise ConfigurationError(
                "ScoreRule.metric is only used to select EvaluationReport.metrics.",
                code="INVALID_SCORECARD_RULE",
                context={"component": name},
            )
        return _ResolvedScore(
            name,
            _number(value.value, name),
            f"MetricResult:{value.name}",
            None,
            value.method,
            value.algorithm_version,
        )
    return _ResolvedScore(name, _number(value, name), "numeric", None, None, None)


def _normalize(value: float, rule: ScoreRule) -> tuple[float, bool]:
    unclipped = (value - rule.minimum) / (rule.maximum - rule.minimum)
    if rule.direction == "lower_is_better":
        unclipped = 1.0 - unclipped
    return min(1.0, max(0.0, unclipped)), not 0.0 <= unclipped <= 1.0


def evaluate_scorecard(
    values: Mapping[str, object],
    *,
    config: ScorecardConfig,
) -> EvaluationReport:
    """Normalize and aggregate named inputs while retaining every contribution."""

    if not isinstance(values, Mapping):
        raise ConfigurationError("values must be a mapping.", code="INVALID_SCORECARD_INPUT")
    if not isinstance(config, ScorecardConfig):
        raise TypeError("config must be ScorecardConfig.")
    unknown = set(values) - set(config.rules)
    if unknown:
        raise ConfigurationError(
            "Scorecard values contain unknown components.",
            code="UNKNOWN_SCORECARD_COMPONENT",
            context={"components": sorted(unknown)},
        )
    resolved = tuple(_resolve(name, values.get(name), rule) for name, rule in config.rules.items())
    missing = tuple(item.name for item in resolved if item.raw_value is None)
    if missing and config.missing_policy == "error":
        raise ConfigurationError(
            "Scorecard is missing required component values.",
            code="MISSING_SCORECARD_COMPONENT",
            context={"components": missing},
        )
    contributions: list[dict[str, object]] = []
    weighted_sum = 0.0
    weight_sum = 0.0
    for item in resolved:
        rule = config.rules[item.name]
        if not isinstance(rule, ScoreRule):
            raise AssertionError("ScorecardConfig accepted a non-ScoreRule value.")
        if item.raw_value is None and config.missing_policy == "omit":
            contributions.append(
                {
                    "name": item.name,
                    "source": item.source,
                    "source_metric": item.metric,
                    "source_method": item.source_method,
                    "source_algorithm_version": item.source_algorithm_version,
                    "raw_value": None,
                    "normalized_value": None,
                    "weight": rule.weight,
                    "weighted_contribution": None,
                    "direction": rule.direction,
                    "minimum": rule.minimum,
                    "maximum": rule.maximum,
                    "status": "omitted",
                    "clipped": False,
                }
            )
            continue
        raw_value = 0.0 if item.raw_value is None else item.raw_value
        normalized, clipped = (
            (0.0, False) if item.raw_value is None else _normalize(raw_value, rule)
        )
        contribution = normalized * rule.weight
        weighted_sum += contribution
        weight_sum += rule.weight
        contributions.append(
            {
                "name": item.name,
                "source": item.source,
                "source_metric": item.metric,
                "source_method": item.source_method,
                "source_algorithm_version": item.source_algorithm_version,
                "raw_value": item.raw_value,
                "imputed_value": raw_value if item.raw_value is None else None,
                "normalized_value": normalized,
                "weight": rule.weight,
                "weighted_contribution": contribution,
                "direction": rule.direction,
                "minimum": rule.minimum,
                "maximum": rule.maximum,
                "status": "zero-imputed" if item.raw_value is None else "included",
                "clipped": clipped,
            }
        )
    if weight_sum == 0:
        raise ConfigurationError(
            "No scorecard component remains after applying missing_policy.",
            code="EMPTY_SCORECARD",
        )
    score = weighted_sum / weight_sum
    status = (
        "pass"
        if score >= config.pass_score
        else ("warning" if score >= config.warning_score else "fail")
    )
    return report(
        name="scorecard",
        method="weighted-clipped-min-max-normalization",
        version="eval-scorecard-v1",
        parameters={
            "rules": {
                name: {
                    "direction": rule.direction,
                    "weight": rule.weight,
                    "minimum": rule.minimum,
                    "maximum": rule.maximum,
                    "metric": rule.metric,
                }
                for name, rule in config.rules.items()
            },
            "formula": "sum(weight_i * clipped_normalized_i) / sum(included_weight_i)",
            "higher_normalization": "clip((x-min)/(max-min),0,1)",
            "lower_normalization": "clip(1-(x-min)/(max-min),0,1)",
            "missing_policy": config.missing_policy,
            "zero_missing_definition": "normalized contribution is zero before weighting",
            "warning_score": config.warning_score,
            "pass_score": config.pass_score,
        },
        metrics={
            "score": score,
            "status": status,
            "weighted_sum": weighted_sum,
            "included_weight_sum": weight_sum,
            "missing_components": missing,
            "contributions": tuple(contributions),
        },
    )


__all__ = ["ScoreInput", "evaluate_scorecard"]
