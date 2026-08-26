"""Generic immutable result values shared by DNAKit domains."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, Generic, TypeVar, cast

from dnakit.core._json import FrozenDict, JSONValue, freeze_json, freeze_mapping, to_json_compatible
from dnakit.core.issues import Issue
from dnakit.core.provenance import Provenance
from dnakit.exceptions import ConfigurationError

T = TypeVar("T")


@dataclass(frozen=True, init=False)
class Uncertainty:
    """Optional uncertainty summary attached to a numerical result."""

    confidence_interval: tuple[float, float] | None
    standard_error: float | None
    method: str | None

    def __init__(
        self,
        confidence_interval: tuple[float, float] | None = None,
        standard_error: float | None = None,
        method: str | None = None,
    ) -> None:
        resolved_interval: tuple[float, float] | None = None
        if confidence_interval is not None:
            copied_interval = tuple(confidence_interval)
            if len(copied_interval) != 2:
                raise ConfigurationError("confidence_interval must contain exactly two bounds.")
            lower, upper = copied_interval
            if (
                isinstance(lower, bool)
                or not isinstance(lower, (int, float))
                or isinstance(upper, bool)
                or not isinstance(upper, (int, float))
                or not math.isfinite(lower)
                or not math.isfinite(upper)
                or lower > upper
            ):
                raise ConfigurationError(
                    "confidence_interval bounds must be finite and ordered.",
                    context={"lower": lower, "upper": upper},
                )
            resolved_interval = (float(lower), float(upper))
        if standard_error is not None and (
            isinstance(standard_error, bool)
            or not isinstance(standard_error, (int, float))
            or not math.isfinite(standard_error)
            or standard_error < 0
        ):
            raise ConfigurationError("standard_error must be a finite non-negative number.")
        if method is not None and (not isinstance(method, str) or not method.strip()):
            raise ConfigurationError("Uncertainty method must be None or a non-empty string.")

        object.__setattr__(self, "confidence_interval", resolved_interval)
        object.__setattr__(
            self,
            "standard_error",
            None if standard_error is None else float(standard_error),
        )
        object.__setattr__(self, "method", method)


@dataclass(frozen=True, init=False)
class MetricResult(Generic[T]):
    """A metric value plus method, conditions, provenance, and diagnostics."""

    name: str
    value: T
    unit: str | None
    method: str
    algorithm_version: str
    parameters: FrozenDict
    conditions: FrozenDict
    uncertainty: Uncertainty | None
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __init__(
        self,
        name: str,
        value: T,
        *,
        method: str,
        algorithm_version: str,
        provenance: Provenance | None = None,
        unit: str | None = None,
        parameters: Mapping[str, object] | None = None,
        conditions: Mapping[str, object] | None = None,
        uncertainty: Uncertainty | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        for field_name, field_value in (
            ("name", name),
            ("method", method),
            ("algorithm_version", algorithm_version),
        ):
            if not isinstance(field_value, str) or not field_value.strip():
                raise ConfigurationError(f"MetricResult {field_name} must be non-empty.")
        if unit is not None and (not isinstance(unit, str) or not unit.strip()):
            raise ConfigurationError("MetricResult unit must be None or a non-empty string.")
        if uncertainty is not None and not isinstance(uncertainty, Uncertainty):
            raise ConfigurationError("MetricResult uncertainty must be Uncertainty or None.")
        resolved_provenance = provenance or Provenance()
        if not isinstance(resolved_provenance, Provenance):
            raise ConfigurationError("MetricResult provenance must be Provenance.")
        issue_tuple = tuple(issues)
        if any(not isinstance(issue, Issue) for issue in issue_tuple):
            raise ConfigurationError("MetricResult issues must all be Issue objects.")
        frozen_value: JSONValue = freeze_json(value)

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "value", cast(T, frozen_value))
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "algorithm_version", algorithm_version)
        object.__setattr__(self, "parameters", freeze_mapping(parameters))
        object.__setattr__(self, "conditions", freeze_mapping(conditions))
        object.__setattr__(self, "uncertainty", uncertainty)
        object.__setattr__(self, "provenance", resolved_provenance)
        object.__setattr__(self, "issues", issue_tuple)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["MetricResult", "Uncertainty"]
