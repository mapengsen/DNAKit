"""Native Python DAG execution and interval-result boolean queries."""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from graphlib import CycleError, TopologicalSorter
from typing import Any

from dnakit.bioinformatics.reads import Region
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True)
class Step:
    """A Python callable and input-parameter -> predecessor-name mapping."""

    name: str
    function: Callable[..., Any]
    dependencies: Mapping[str, str] = field(default_factory=dict)
    parameters: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class StepResult:
    name: str
    status: str
    value: Any
    elapsed_seconds: float
    error: str | None = None


def run_workflow(steps: Iterable[Step], *, fail_fast: bool = True) -> tuple[StepResult, ...]:
    """Run a checked DAG; failures block dependents, independent steps can continue.

    Functions are supplied as actual Python callables, never evaluated/imported
    from untrusted configuration strings. Results retain intermediate objects.
    """
    steps = tuple(steps)
    lookup = {step.name: step for step in steps}
    if len(lookup) != len(steps) or any(
        not step.name or not callable(step.function) for step in steps
    ):
        raise ConfigurationError("Workflow names must be unique and functions callable.")
    for step in steps:
        if set(step.dependencies) & set(step.parameters) or any(
            name not in lookup for name in step.dependencies.values()
        ):
            raise ConfigurationError("Invalid dependency or overlapping input parameter.")
    try:
        order = tuple(
            TopologicalSorter(
                {step.name: set(step.dependencies.values()) for step in steps}
            ).static_order()
        )
    except CycleError as exc:
        raise ConfigurationError("Workflow contains a dependency cycle.") from exc
    results: dict[str, StepResult] = {}
    for name in order:
        step = lookup[name]
        if any(results[parent].status != "completed" for parent in step.dependencies.values()):
            results[name] = StepResult(name, "blocked", None, 0.0, "A predecessor failed.")
            continue
        started = time.monotonic()
        try:
            inputs = {
                parameter: results[parent].value for parameter, parent in step.dependencies.items()
            }
            value = step.function(**inputs, **step.parameters)
        except Exception as exc:
            if fail_fast:
                raise
            results[name] = StepResult(
                name, "failed", None, time.monotonic() - started, f"{type(exc).__name__}: {exc}"
            )
        else:
            results[name] = StepResult(name, "completed", value, time.monotonic() - started)
    return tuple(results[name] for name in order)


def query_regions(
    left: Iterable[Region], right: Iterable[Region], *, operation: str = "intersection"
) -> tuple[Region, ...]:
    """Combine intervals by reference and strand, merging adjacent spans."""
    if operation not in {"union", "intersection", "difference"}:
        raise ConfigurationError("Unknown region query operation.")
    groups: dict[tuple[str, str], list[tuple[int, int, int]]] = {}
    for side, regions in enumerate((left, right)):
        for region in regions:
            groups.setdefault((region.reference, region.strand), []).extend(
                ((region.start, side, 1), (region.end, side, -1))
            )
    result: list[Region] = []
    for (reference, strand), events in sorted(groups.items()):
        positions = sorted({event[0] for event in events})
        deltas: dict[int, list[int]] = {position: [0, 0] for position in positions}
        for position, side, delta in events:
            deltas[position][side] += delta
        active = [0, 0]
        start: int | None = None
        for position in positions:
            active = [active[side] + deltas[position][side] for side in (0, 1)]
            keep = (
                any(active)
                if operation == "union"
                else (
                    all(active) if operation == "intersection" else active[0] > 0 and active[1] == 0
                )
            )
            if keep and start is None:
                start = position
            elif not keep and start is not None:
                result.append(Region(reference, start, position, strand))
                start = None
    return tuple(result)
