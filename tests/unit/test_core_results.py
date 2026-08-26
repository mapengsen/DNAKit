"""Tests for structured diagnostics and provenance-bearing results."""

import json

import pytest

from dnakit.core import (
    BackendInfo,
    ImplementationInfo,
    ImplementationLabel,
    Interval,
    Issue,
    IssueSeverity,
    MetricResult,
    Provenance,
    Uncertainty,
)
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError


def test_issue_and_metric_result_are_json_serializable_and_defensive() -> None:
    raw_parameters = {"window": [1, 2]}
    issue = Issue(
        "AMBIGUOUS_BASE",
        IssueSeverity.WARNING,
        "Sequence contains N.",
        Interval(3, 4),
        {"symbol": "N"},
    )
    result = MetricResult(
        "gc_fraction",
        0.5,
        method="count",
        algorithm_version="1",
        unit="fraction",
        parameters=raw_parameters,
        uncertainty=Uncertainty(confidence_interval=(0.4, 0.6)),
        issues=[issue],
    )
    raw_parameters["window"].append(3)

    payload = result.to_dict()
    assert payload["parameters"] == {"window": [1, 2]}
    assert payload["issues"][0]["severity"] == "warning"
    assert json.loads(json.dumps(payload))["value"] == 0.5


def test_provenance_has_one_authoritative_implementation_and_backend() -> None:
    backend = BackendInfo("primer3", version="2.6.1", capabilities=["tm"])
    implementation = ImplementationInfo(label=ImplementationLabel.ADAPTER)
    provenance = Provenance(implementation=implementation, backend=backend)
    result = MetricResult(
        "tm",
        60.0,
        method="primer3",
        algorithm_version="2.6.1",
        provenance=provenance,
        unit="degC",
    )

    assert result.provenance.implementation.label is ImplementationLabel.ADAPTER
    assert result.provenance.backend is backend
    assert not hasattr(result, "implementation")
    assert not hasattr(implementation, "backend")


def test_backend_capabilities_have_stable_json_order() -> None:
    backend = BackendInfo(
        "example",
        capabilities=frozenset({"gamma", "alpha", "delta", "beta"}),
    )

    assert backend.to_dict()["capabilities"] == ["alpha", "beta", "delta", "gamma"]


def test_result_rejects_non_json_values_and_invalid_uncertainty() -> None:
    with pytest.raises(ConfigurationError, match="JSON-compatible"):
        MetricResult("bad", object(), method="test", algorithm_version="1")
    with pytest.raises(ConfigurationError, match="ordered"):
        Uncertainty(confidence_interval=(2.0, 1.0))
    with pytest.raises(ConfigurationError, match="finite"):
        Uncertainty(standard_error=float("nan"))


def test_uncertainty_defensively_copies_confidence_interval() -> None:
    raw_interval = [0.1, 0.9]
    uncertainty = Uncertainty(confidence_interval=raw_interval)  # type: ignore[arg-type]
    raw_interval[0] = 0.5

    assert uncertainty.confidence_interval == (0.1, 0.9)


def test_exceptions_expose_stable_codes_context_and_hint() -> None:
    error = InvalidAlphabetError(
        "Invalid symbol.",
        context={"symbol": "Z"},
        hint="Use IUPAC symbols.",
    )

    assert error.code == error.error_code == "INVALID_ALPHABET"
    assert error.context["symbol"] == "Z"
    assert "Use IUPAC symbols" in str(error)
    with pytest.raises(TypeError):
        error.context["symbol"] = "N"  # type: ignore[index]
