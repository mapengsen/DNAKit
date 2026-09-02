"""Tests for capability-named optional scientific functions."""

from __future__ import annotations

from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError
from types import ModuleType

import pytest

import dnakit.annotation.variants as variant_annotation
import dnakit.backends.scientific as scientific_backend
import dnakit.comparative.selection as comparative_selection
import dnakit.molbio.golden_gate as golden_gate
from dnakit.annotation import (
    annotate_rsid_vep,
    annotate_variant_vep,
    get_clinvar_significance,
    get_clinvar_variant,
    get_dbsnp_frequencies,
    get_dbsnp_variant,
    get_gnomad_population_frequencies,
    get_gnomad_variant,
    recode_variant,
    search_clinvar_variants,
    search_gnomad_variants,
)
from dnakit.comparative import calculate_dn_ds
from dnakit.core import DNASequence, Provenance, ProviderResult
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.molbio import assemble_golden_gate, design_golden_gate


def _result(
    operation: str,
    arguments: Mapping[str, object],
    *,
    parameters: Mapping[str, object] | None = None,
) -> ProviderResult:
    return ProviderResult(
        operation,
        "test-provider",
        "test-backend",
        {"arguments": dict(arguments)},
        Provenance(),
        parameters=parameters,
    )


def test_backend_loads_only_requested_tool_and_unwraps_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[object] = []

    class FakeToolUniverse:
        def __init__(self, *, workspace: str) -> None:
            events.append(("init", workspace))

        def load_tools(self, *, include_tools: list[str], quiet: bool) -> None:
            events.append(("load", tuple(include_tools), quiet))

        def run(self, call: dict[str, object], *, verbose: bool) -> dict[str, object]:
            events.append(("run", call, verbose))
            return {
                "status": "success",
                "data": {"dN_dS": 0.5},
                "metadata": {"engine": "nei_gojobori_1986"},
            }

    module = ModuleType("tooluniverse")
    module.ToolUniverse = FakeToolUniverse  # type: ignore[attr-defined]
    module.__version__ = "1.4.1"  # type: ignore[attr-defined]
    module.__file__ = "/optional/tooluniverse/__init__.py"
    monkeypatch.setattr(scientific_backend, "_tooluniverse_module", lambda: module)

    def missing_distribution(_name: str) -> str:
        raise PackageNotFoundError

    monkeypatch.setattr(scientific_backend, "version", missing_distribution)

    result = scientific_backend._run_scientific_function(
        "calculate_dn_ds",
        {"seq1": "ATG", "seq2": "ATA"},
        parameters={"sequence_length": 3},
    )

    assert events[0][0] == "init"  # type: ignore[index]
    assert isinstance(events[0][1], str)  # type: ignore[index]
    assert events[1] == ("load", ("Sequence_dn_ds",), True)
    assert result.operation == "calculate_dn_ds"
    assert result.backend == "ToolUniverse"
    assert result.data == {"dN_dS": 0.5}
    assert result.parameters == {"sequence_length": 3}
    assert result.metadata["engine"] == "nei_gojobori_1986"
    assert result.metadata["backend_tool"] == "Sequence_dn_ds"
    assert result.provenance.backend is not None
    assert result.provenance.backend.version == "1.4.1"
    assert result.to_dict()["data"] == {"dN_dS": 0.5}


def test_backend_reports_missing_dependency_and_provider_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def missing_module(_name: str) -> ModuleType:
        raise ModuleNotFoundError("tooluniverse")

    monkeypatch.setattr(
        "dnakit.backends.scientific.importlib.import_module",
        missing_module,
    )
    with pytest.raises(BackendUnavailableError) as missing:
        scientific_backend._run_scientific_function("get_dbsnp_variant", {"rsid": "rs1"})
    assert missing.value.code == "TOOLUNIVERSE_UNAVAILABLE"
    assert "dnakit[external-tools]" in str(missing.value)

    class FailedToolUniverse:
        def __init__(self, *, load_workspace: bool) -> None:
            del load_workspace

        def load_tools(self, *, include_tools: list[str], quiet: bool) -> None:
            del include_tools, quiet

        def run(self, call: dict[str, object], *, verbose: bool) -> dict[str, object]:
            del call, verbose
            return {"status": "error", "error": "remote job failed"}

    module = ModuleType("tooluniverse")
    module.ToolUniverse = FailedToolUniverse  # type: ignore[attr-defined]
    module.__version__ = "1.4.1"  # type: ignore[attr-defined]
    monkeypatch.setattr(scientific_backend, "_tooluniverse_module", lambda: module)
    with pytest.raises(BackendExecutionError) as failed:
        scientific_backend._run_scientific_function("get_dbsnp_variant", {"rsid": "rs1"})
    assert failed.value.code == "SCIENTIFIC_PROVIDER_ERROR"
    assert "remote job failed" in str(failed.value)


def test_variant_annotation_functions_use_explicit_capability_names(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, dict[str, object], dict[str, object]]] = []

    def fake_run(
        operation: str,
        arguments: Mapping[str, object],
        *,
        parameters: Mapping[str, object] | None = None,
    ) -> ProviderResult:
        calls.append((operation, dict(arguments), dict(parameters or {})))
        return _result(operation, arguments, parameters=parameters)

    monkeypatch.setattr(variant_annotation, "_run_scientific_function", fake_run)

    assert annotate_variant_vep("NM_000546.6:c.215C>G").operation == "annotate_variant_vep"
    annotate_rsid_vep("rs28934578")
    recode_variant("rs28934578")
    search_clinvar_variants(gene="TP53", variant_name="c.215C>G", limit=5)
    get_clinvar_variant("12345")
    get_clinvar_significance("12345")
    get_dbsnp_variant("rs28934578")
    get_dbsnp_frequencies("rs28934578")
    search_gnomad_variants("rs28934578")
    get_gnomad_variant("17-7675088-C-G")
    get_gnomad_population_frequencies("17-7675088-C-G")

    assert [call[0] for call in calls] == [
        "annotate_variant_vep",
        "annotate_rsid_vep",
        "recode_variant",
        "search_clinvar_variants",
        "get_clinvar_variant",
        "get_clinvar_significance",
        "get_dbsnp_variant",
        "get_dbsnp_frequencies",
        "search_gnomad_variants",
        "get_gnomad_variant",
        "get_gnomad_population_frequencies",
    ]
    assert calls[3][1]["max_results"] == 5
    assert calls[-1][1]["dataset"] == "gnomad_r4"

    with pytest.raises(ConfigurationError, match="at least one"):
        search_clinvar_variants()
    with pytest.raises(ConfigurationError, match="Unsupported gnomAD"):
        get_gnomad_variant("17-1-A-G", dataset="unknown")


def test_dn_ds_and_golden_gate_validate_then_route_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, dict[str, object], dict[str, object]]] = []

    def fake_run(
        operation: str,
        arguments: Mapping[str, object],
        *,
        parameters: Mapping[str, object] | None = None,
    ) -> ProviderResult:
        calls.append((operation, dict(arguments), dict(parameters or {})))
        return _result(operation, arguments, parameters=parameters)

    monkeypatch.setattr(comparative_selection, "_run_scientific_function", fake_run)
    monkeypatch.setattr(golden_gate, "_run_scientific_function", fake_run)

    calculate_dn_ds(DNASequence("ATGGCT"), DNASequence("ATGGCC"))
    design_golden_gate((DNASequence("ATGGCT"), DNASequence("GCCGAA")), enzyme="BbsI")
    assemble_golden_gate(
        (DNASequence("GGTCTCATGGCT"), DNASequence("GGTCTCGCCGAA")),
        labels=("vector", "insert"),
    )

    assert calls[0][0] == "calculate_dn_ds"
    assert calls[0][2] == {"sequence_length": 6, "codon_count": 2}
    assert calls[1][0] == "design_golden_gate"
    assert calls[1][2]["part_lengths"] == (6, 6)
    assert calls[2][0] == "assemble_golden_gate"
    assert calls[2][1]["labels"] == ["vector", "insert"]

    with pytest.raises(ConfigurationError, match="equal lengths"):
        calculate_dn_ds(DNASequence("ATG"), DNASequence("ATGGCT"))
    with pytest.raises(ConfigurationError, match="complete aligned codons"):
        calculate_dn_ds(DNASequence("ATGG"), DNASequence("ATGG"))
    with pytest.raises(ConfigurationError, match="at least two"):
        design_golden_gate((DNASequence("ATGGCT"),))
    with pytest.raises(ConfigurationError, match="one non-empty label"):
        assemble_golden_gate(
            (DNASequence("GGTCTCATGGCT"), DNASequence("GGTCTCGCCGAA")),
            labels=("only-one",),
        )
