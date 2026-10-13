"""Dependency-free assembly contiguity statistics.

The functions in this module operate only on explicit contig lengths.  They
do not read an assembly graph, inspect a reference genome, or run an external
assembler.  Returning the coordinate-free summary with provenance makes the
metric definition reproducible in a methods section.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit.exceptions import ConfigurationError


def _envelope(operation: str, result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "result": dict(result),
        "provenance": {
            "implementation": "dnakit.native_assembly_metrics",
            "backend": "local_dependency_free",
            "scientific_boundary": "contiguity statistics over explicit contig lengths; no reference or assembler",
        },
    }


def _lengths(value: object) -> list[int]:
    """Extract positive contig lengths from common JSON representations."""

    if isinstance(value, Mapping):
        for key in ("lengths", "contig_lengths", "contigs", "records", "items", "data"):
            if key in value:
                return _lengths(value[key])
        raise ConfigurationError(
            "assembly input must contain lengths, contig_lengths, contigs, or records."
        )
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("assembly input must be a sequence of positive contig lengths.")
    result: list[int] = []
    for index, item in enumerate(value):
        if isinstance(item, bool):
            raise ConfigurationError(f"contig length {index} must be a positive integer.")
        if isinstance(item, Mapping):
            length = item.get("length", item.get("size", item.get("contig_length")))
        else:
            length = item
        if isinstance(length, bool) or not isinstance(length, (int, float)):
            raise ConfigurationError(f"contig length {index} must be numeric.")
        if isinstance(length, float) and not length.is_integer():
            raise ConfigurationError(f"contig length {index} must be an integer.")
        normalized = int(length)
        if normalized <= 0:
            raise ConfigurationError(f"contig length {index} must be positive.")
        result.append(normalized)
    if not result:
        raise ConfigurationError("assembly input must contain at least one contig length.")
    return result


def _target_percentile(options: Mapping[str, object]) -> float:
    value = options.get("percentile", options.get("x", 50))
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < float(value) <= 100:
        raise ConfigurationError("options.percentile/options.x must be in the interval (0, 100].")
    return float(value)


def _nx(lengths: Sequence[int], percentile: float, denominator: int | None = None) -> tuple[int, int]:
    total = int(sum(lengths) if denominator is None else denominator)
    if total <= 0:
        raise ConfigurationError("the N-statistic denominator must be positive.")
    target = total * percentile / 100.0
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= target:
            return length, rank
    return 0, 0


def calculate_nx(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate an Nx statistic for an explicit contig-length list."""

    settings = dict(options or {})
    percentile = _target_percentile(settings)
    lengths = _lengths(input)
    value, rank = _nx(lengths, percentile)
    return _envelope(
        "calculate_nx",
        {
            "percentile": percentile,
            "nx": value,
            "contigs_to_threshold": rank,
            "contig_count": len(lengths),
            "total_length": sum(lengths),
        },
    )


def calculate_n50(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate N50 from explicit contig lengths."""

    settings = dict(options or {})
    lengths = _lengths(input)
    value, rank = _nx(lengths, 50.0)
    return _envelope(
        "calculate_n50",
        {
            "n50": value,
            "l50": rank,
            "contig_count": len(lengths),
            "total_length": sum(lengths),
            "percentile": 50.0,
            "options_used": sorted(settings),
        },
    )


def calculate_l50(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate L50, the rank of the contig crossing half the assembly."""

    settings = dict(options or {})
    percentile = _target_percentile(settings)
    lengths = _lengths(input)
    value, rank = _nx(lengths, percentile)
    return _envelope(
        "calculate_l50",
        {
            "l50": rank,
            "n_at_percentile": value,
            "percentile": percentile,
            "contig_count": len(lengths),
            "total_length": sum(lengths),
        },
    )


def calculate_ng50(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Calculate NG50 using an explicit expected genome size denominator."""

    settings = dict(options or {})
    expected = settings.get("genome_size", settings.get("expected_genome_size"))
    if isinstance(expected, bool) or not isinstance(expected, (int, float)) or expected <= 0:
        raise ConfigurationError("options.genome_size must be a positive number for NG50.")
    if isinstance(expected, float) and not expected.is_integer():
        expected_size = int(round(expected))
    else:
        expected_size = int(expected)
    percentile = _target_percentile(settings)
    lengths = _lengths(input)
    value, rank = _nx(lengths, percentile, expected_size)
    return _envelope(
        "calculate_ng50",
        {
            "ng": value,
            "ng50": value if percentile == 50.0 else None,
            "contigs_to_threshold": rank,
            "percentile": percentile,
            "expected_genome_size": expected_size,
            "assembly_length": sum(lengths),
        },
    )


def calculate_stats(input: object, *, options: Mapping[str, object] | None = None) -> dict[str, Any]:
    """Return common assembly length, Nx and Lx statistics."""

    settings = dict(options or {})
    lengths = _lengths(input)
    ordered = sorted(lengths, reverse=True)
    n50, l50 = _nx(lengths, 50.0)
    n90, l90 = _nx(lengths, 90.0)
    result: dict[str, Any] = {
        "contig_count": len(lengths),
        "total_length": sum(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "l90": l90,
        "sorted_lengths": ordered,
        "options_used": sorted(settings),
    }
    if "genome_size" in settings or "expected_genome_size" in settings:
        ng50 = calculate_ng50(lengths, options={**settings, "percentile": 50})
        result["ng50"] = ng50["result"]["ng50"]
    return _envelope("calculate_stats", result)


__all__ = [
    "calculate_l50",
    "calculate_n50",
    "calculate_ng50",
    "calculate_nx",
    "calculate_stats",
]
