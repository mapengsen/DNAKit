"""Allowlisted optional scientific-function backend used by domain APIs."""

from __future__ import annotations

import asyncio
import importlib
import inspect
import json
import tempfile
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from types import ModuleType
from typing import Any

from dnakit.core import (
    BackendInfo,
    Citation,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    ProviderResult,
    ReferenceInfo,
)
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
)

_MAX_RESULT_BYTES = 20_000_000
_TOOLUNIVERSE_CITATION = Citation(
    "tooluniverse",
    title="ToolUniverse: Democratizing AI Scientists",
    url="https://arxiv.org/abs/2509.23426",
)


@dataclass(frozen=True, slots=True)
class _ScientificFunction:
    operation: str
    tool_name: str
    provider: str
    citation_url: str
    remote: bool = True


_FUNCTIONS = {
    item.operation: item
    for item in (
        _ScientificFunction(
            "annotate_variant_vep",
            "EnsemblVEP_annotate_hgvs",
            "Ensembl VEP",
            "https://rest.ensembl.org/documentation/info/vep_hgvs_get",
        ),
        _ScientificFunction(
            "annotate_rsid_vep",
            "EnsemblVEP_annotate_rsid",
            "Ensembl VEP",
            "https://rest.ensembl.org/documentation/info/vep_id_get",
        ),
        _ScientificFunction(
            "recode_variant",
            "EnsemblVEP_variant_recoder",
            "Ensembl Variant Recoder",
            "https://rest.ensembl.org/documentation/info/variant_recoder",
        ),
        _ScientificFunction(
            "search_clinvar_variants",
            "ClinVar_search_variants",
            "NCBI ClinVar",
            "https://www.ncbi.nlm.nih.gov/clinvar/",
        ),
        _ScientificFunction(
            "get_clinvar_variant",
            "ClinVar_get_variant_details",
            "NCBI ClinVar",
            "https://www.ncbi.nlm.nih.gov/clinvar/",
        ),
        _ScientificFunction(
            "get_clinvar_significance",
            "ClinVar_get_clinical_significance",
            "NCBI ClinVar",
            "https://www.ncbi.nlm.nih.gov/clinvar/",
        ),
        _ScientificFunction(
            "get_dbsnp_variant",
            "dbsnp_get_variant_by_rsid",
            "NCBI dbSNP",
            "https://www.ncbi.nlm.nih.gov/snp/",
        ),
        _ScientificFunction(
            "get_dbsnp_frequencies",
            "dbsnp_get_frequencies",
            "NCBI dbSNP",
            "https://www.ncbi.nlm.nih.gov/snp/",
        ),
        _ScientificFunction(
            "search_gnomad_variants",
            "gnomad_search_variants",
            "gnomAD",
            "https://gnomad.broadinstitute.org/help/api",
        ),
        _ScientificFunction(
            "get_gnomad_variant",
            "gnomad_get_variant",
            "gnomAD",
            "https://gnomad.broadinstitute.org/help/api",
        ),
        _ScientificFunction(
            "get_gnomad_population_frequencies",
            "gnomad_get_variant_populations",
            "gnomAD",
            "https://gnomad.broadinstitute.org/help/api",
        ),
        _ScientificFunction(
            "calculate_dn_ds",
            "Sequence_dn_ds",
            "Nei-Gojobori dN/dS",
            "https://doi.org/10.1093/oxfordjournals.molbev.a040410",
            remote=False,
        ),
        _ScientificFunction(
            "design_golden_gate",
            "DNA_golden_gate_design",
            "Golden Gate assembly",
            "https://doi.org/10.1371/journal.pone.0003647",
            remote=False,
        ),
        _ScientificFunction(
            "assemble_golden_gate",
            "DNA_golden_gate_assemble",
            "Golden Gate assembly",
            "https://doi.org/10.1371/journal.pone.0003647",
            remote=False,
        ),
    )
}


def _tooluniverse_module() -> ModuleType:
    try:
        return importlib.import_module("tooluniverse")
    except (ImportError, ModuleNotFoundError) as exc:
        raise BackendUnavailableError(
            "Optional scientific functions require ToolUniverse.",
            code="TOOLUNIVERSE_UNAVAILABLE",
            hint=(
                'Install the optional backend with: python -m pip install "dnakit[external-tools]"'
            ),
        ) from exc
    except Exception as exc:
        raise BackendUnavailableError(
            "ToolUniverse could not be imported in this environment.",
            code="TOOLUNIVERSE_IMPORT_FAILED",
            context={"error_type": type(exc).__name__},
            hint=(
                "Reinstall the optional backend with: "
                'python -m pip install "dnakit[external-tools]"'
            ),
        ) from exc


def _backend_version(module: ModuleType) -> str:
    try:
        return version("tooluniverse")
    except PackageNotFoundError:
        candidate = getattr(module, "__version__", "unknown")
        return candidate if isinstance(candidate, str) and candidate.strip() else "unknown"


def _invoke(universe: Any, call: Mapping[str, object]) -> object:
    def execute() -> object:
        return universe.run(dict(call), verbose=False)

    try:
        asyncio.get_running_loop()
    except RuntimeError:
        result = execute()
    else:
        # ToolUniverse changes ``run`` into an async call when it sees a running
        # event loop.  These DNAKit APIs are synchronous, so execute in an
        # isolated worker and preserve the same public return type.
        with ThreadPoolExecutor(max_workers=1, thread_name_prefix="dnakit-provider") as pool:
            result = pool.submit(execute).result()
    if inspect.isawaitable(result):
        close = getattr(result, "close", None)
        if callable(close):
            close()
        raise BackendExecutionError(
            "Scientific provider returned an unexpected asynchronous result.",
            code="PROVIDER_ASYNC_RESULT",
        )
    return result


def _unwrap_result(result: object, *, tool_name: str) -> tuple[object, dict[str, object]]:
    metadata: dict[str, object] = {}
    if isinstance(result, Mapping):
        status = result.get("status")
        normalized_status = status.casefold() if isinstance(status, str) else None
        error = result.get("error")
        if normalized_status in {"error", "failed", "failure"} or error:
            message = error if isinstance(error, str) and error.strip() else "Provider call failed."
            raise BackendExecutionError(
                message,
                code="SCIENTIFIC_PROVIDER_ERROR",
                context={"tool": tool_name},
            )
        if normalized_status in {"ok", "success"} and "data" in result:
            raw_metadata = result.get("metadata")
            if isinstance(raw_metadata, Mapping):
                metadata.update({str(key): value for key, value in raw_metadata.items()})
            return result["data"], metadata
        return dict(result), metadata
    if result is None:
        raise BackendExecutionError(
            "Scientific provider returned no result.",
            code="EMPTY_SCIENTIFIC_PROVIDER_RESULT",
            context={"tool": tool_name},
        )
    return result, metadata


def _validate_result_size(data: object, metadata: Mapping[str, object]) -> None:
    try:
        encoded = json.dumps(
            {"data": data, "metadata": dict(metadata)},
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise BackendExecutionError(
            "Scientific provider returned a non-JSON-compatible result.",
            code="INVALID_SCIENTIFIC_PROVIDER_RESULT",
        ) from exc
    if len(encoded) > _MAX_RESULT_BYTES:
        raise BackendExecutionError(
            "Scientific provider result exceeded the 20 MB safety limit.",
            code="SCIENTIFIC_PROVIDER_RESULT_LIMIT",
            context={"observed_bytes": len(encoded), "max_bytes": _MAX_RESULT_BYTES},
        )


def _run_scientific_function(
    operation: str,
    arguments: Mapping[str, object],
    *,
    parameters: Mapping[str, object] | None = None,
) -> ProviderResult:
    """Run one allowlisted backend tool and return a DNAKit result envelope."""

    spec = _FUNCTIONS.get(operation)
    if spec is None:
        raise ConfigurationError(
            "Unknown optional scientific function.",
            code="UNKNOWN_SCIENTIFIC_FUNCTION",
            context={"operation": operation},
        )
    module = _tooluniverse_module()
    tooluniverse_class = getattr(module, "ToolUniverse", None)
    if not callable(tooluniverse_class):
        raise BackendUnavailableError(
            "Installed ToolUniverse does not expose ToolUniverse().",
            code="TOOLUNIVERSE_API_UNAVAILABLE",
        )
    backend_version = _backend_version(module)
    temporary_workspace: tempfile.TemporaryDirectory[str] | None = None
    try:
        constructor_parameters: Mapping[str, inspect.Parameter]
        try:
            constructor_parameters = inspect.signature(tooluniverse_class).parameters
        except (TypeError, ValueError):
            constructor_parameters = {}
        constructor_arguments: dict[str, object] = {}
        if "log_level" in constructor_parameters:
            constructor_arguments["log_level"] = "ERROR"
        if "load_workspace" in constructor_parameters:
            constructor_arguments["load_workspace"] = False
        elif "workspace" in constructor_parameters:
            temporary_workspace = tempfile.TemporaryDirectory(
                prefix=".dnakit-tooluniverse-",
                ignore_cleanup_errors=True,
            )
            profile = {
                "name": "dnakit-isolated",
                "tools": {"include_tools": [spec.tool_name]},
                "cache": {"enabled": False, "persist": False},
                "sources": [],
            }
            Path(temporary_workspace.name, "profile.yaml").write_text(
                json.dumps(profile, separators=(",", ":")),
                encoding="utf-8",
            )
            constructor_arguments["workspace"] = temporary_workspace.name
        universe = tooluniverse_class(**constructor_arguments)
        universe.load_tools(include_tools=[spec.tool_name], quiet=True)
        raw_result = _invoke(
            universe,
            {"name": spec.tool_name, "arguments": dict(arguments)},
        )
        data, upstream_metadata = _unwrap_result(raw_result, tool_name=spec.tool_name)
        _validate_result_size(data, upstream_metadata)
    except (BackendExecutionError, BackendUnavailableError):
        raise
    except Exception as exc:
        raise BackendExecutionError(
            "Optional scientific function execution failed.",
            code="SCIENTIFIC_FUNCTION_EXECUTION_FAILED",
            context={"operation": operation, "error_type": type(exc).__name__},
        ) from exc
    finally:
        if temporary_workspace is not None:
            temporary_workspace.cleanup()

    audit_parameters = dict(parameters or {})
    package_location = getattr(module, "__file__", None)
    backend = BackendInfo(
        "ToolUniverse",
        version=backend_version,
        package_location=(
            package_location
            if isinstance(package_location, str) and package_location.strip()
            else None
        ),
        license_expression="Apache-2.0",
        capabilities=(operation,),
        metadata={"tool_name": spec.tool_name},
    )
    provenance = Provenance(
        dependency_versions={"tooluniverse": backend_version},
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=(ExecutionMode.EXTERNAL if spec.remote else ExecutionMode.INTERNAL),
            origin_class=OriginClass.INTEGRATION,
            license_expression="Apache-2.0",
            citations=(
                _TOOLUNIVERSE_CITATION,
                Citation(
                    spec.provider.casefold().replace(" ", "-"),
                    title=spec.provider,
                    url=spec.citation_url,
                ),
            ),
        ),
        backend=backend,
        reference=ReferenceInfo(spec.provider, filters=audit_parameters),
    )
    metadata = {
        "backend_tool": spec.tool_name,
        "backend_version": backend_version,
        **upstream_metadata,
    }
    try:
        return ProviderResult(
            operation,
            spec.provider,
            "ToolUniverse",
            data,
            provenance,
            parameters=audit_parameters,
            metadata=metadata,
        )
    except ConfigurationError as exc:
        raise BackendExecutionError(
            "Scientific provider result could not be normalized.",
            code="INVALID_SCIENTIFIC_PROVIDER_RESULT",
            context={"operation": operation},
        ) from exc


__all__: list[str] = []
