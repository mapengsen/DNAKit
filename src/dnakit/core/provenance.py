"""Implementation, reference, artifact, and run provenance value objects."""

from __future__ import annotations

import platform as platform_module
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from typing import Any, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.backend_info import BackendInfo
from dnakit.core.enums import ExecutionMode, ImplementationLabel, OriginClass
from dnakit.core.issues import Issue
from dnakit.exceptions import ConfigurationError

_SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")


def _installed_dnakit_version() -> str:
    try:
        return version("dnakit")
    except PackageNotFoundError:
        return "unknown"


def _non_empty_optional(value: str | None, name: str) -> None:
    if value is not None and (not isinstance(value, str) or not value.strip()):
        raise ConfigurationError(f"{name} must be None or a non-empty string.")


@dataclass(frozen=True)
class Citation:
    """A compact citation identifier with optional resolvable metadata."""

    key: str
    title: str | None = None
    doi: str | None = None
    url: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.key, str) or not self.key.strip():
            raise ConfigurationError("Citation key must be a non-empty string.")
        for name in ("title", "doi", "url"):
            _non_empty_optional(getattr(self, name), f"Citation {name}")


@dataclass(frozen=True, init=False)
class ImplementationInfo:
    """The single authoritative classification of an implementation."""

    label: ImplementationLabel
    execution_mode: ExecutionMode
    origin_class: OriginClass
    license_expression: str | None
    citations: tuple[Citation, ...]

    def __init__(
        self,
        label: ImplementationLabel | str = ImplementationLabel.NATIVE,
        execution_mode: ExecutionMode | str = ExecutionMode.INTERNAL,
        origin_class: OriginClass | str = OriginClass.DNAKIT,
        license_expression: str | None = None,
        citations: Iterable[Citation] = (),
    ) -> None:
        try:
            resolved_label = (
                label if isinstance(label, ImplementationLabel) else ImplementationLabel(label)
            )
            resolved_execution = (
                execution_mode
                if isinstance(execution_mode, ExecutionMode)
                else ExecutionMode(execution_mode)
            )
            resolved_origin = (
                origin_class if isinstance(origin_class, OriginClass) else OriginClass(origin_class)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Invalid implementation classification.",
                context={
                    "label": label,
                    "execution_mode": execution_mode,
                    "origin_class": origin_class,
                },
            ) from exc
        _non_empty_optional(license_expression, "Implementation license_expression")
        citation_tuple = tuple(citations)
        if any(not isinstance(citation, Citation) for citation in citation_tuple):
            raise ConfigurationError("Implementation citations must all be Citation objects.")

        object.__setattr__(self, "label", resolved_label)
        object.__setattr__(self, "execution_mode", resolved_execution)
        object.__setattr__(self, "origin_class", resolved_origin)
        object.__setattr__(self, "license_expression", license_expression)
        object.__setattr__(self, "citations", citation_tuple)


@dataclass(frozen=True, init=False)
class ReferenceInfo:
    """Versioned reference database metadata."""

    name: str
    version: str | None
    date: str | None
    checksum: str | None
    filters: FrozenDict

    def __init__(
        self,
        name: str,
        *,
        version: str | None = None,
        date: str | None = None,
        checksum: str | None = None,
        filters: Mapping[str, object] | None = None,
    ) -> None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("Reference name must be a non-empty string.")
        for field_name, value in (
            ("version", version),
            ("date", date),
            ("checksum", checksum),
        ):
            _non_empty_optional(value, f"Reference {field_name}")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "date", date)
        object.__setattr__(self, "checksum", checksum)
        object.__setattr__(self, "filters", freeze_mapping(filters))


@dataclass(frozen=True, init=False)
class Provenance:
    """Environment and implementation identity attached to a result."""

    dnakit_version: str
    python_version: str
    platform: str
    dependency_versions: FrozenDict
    implementation: ImplementationInfo
    backend: BackendInfo | None
    reference: ReferenceInfo | None

    def __init__(
        self,
        dnakit_version: str | None = None,
        python_version: str | None = None,
        platform: str | None = None,
        dependency_versions: Mapping[str, str] | None = None,
        implementation: ImplementationInfo | None = None,
        backend: BackendInfo | None = None,
        reference: ReferenceInfo | None = None,
    ) -> None:
        resolved_dnakit = _installed_dnakit_version() if dnakit_version is None else dnakit_version
        resolved_python = (
            platform_module.python_version() if python_version is None else python_version
        )
        resolved_platform = platform_module.platform() if platform is None else platform
        for field_name, value in (
            ("dnakit_version", resolved_dnakit),
            ("python_version", resolved_python),
            ("platform", resolved_platform),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Provenance {field_name} must be non-empty.")
        if dependency_versions is not None and any(
            not isinstance(name, str)
            or not name.strip()
            or not isinstance(value, str)
            or not value.strip()
            for name, value in dependency_versions.items()
        ):
            raise ConfigurationError("Dependency versions must map non-empty strings to strings.")
        resolved_implementation = implementation or ImplementationInfo()
        if not isinstance(resolved_implementation, ImplementationInfo):
            raise ConfigurationError("Provenance implementation must be ImplementationInfo.")
        if backend is not None and not isinstance(backend, BackendInfo):
            raise ConfigurationError("Provenance backend must be BackendInfo or None.")
        if reference is not None and not isinstance(reference, ReferenceInfo):
            raise ConfigurationError("Provenance reference must be ReferenceInfo or None.")

        object.__setattr__(self, "dnakit_version", resolved_dnakit)
        object.__setattr__(self, "python_version", resolved_python)
        object.__setattr__(self, "platform", resolved_platform)
        object.__setattr__(self, "dependency_versions", freeze_mapping(dependency_versions))
        object.__setattr__(self, "implementation", resolved_implementation)
        object.__setattr__(self, "backend", backend)
        object.__setattr__(self, "reference", reference)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class ArtifactRef:
    """Integrity-checked reference to a persisted result artifact."""

    relative_path: str
    media_type: str
    schema_version: str
    sha256: str
    byte_size: int
    created_at: str

    def __post_init__(self) -> None:
        for name in ("relative_path", "media_type", "schema_version", "created_at"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Artifact {name} must be a non-empty string.")
        if not isinstance(self.sha256, str) or not _SHA256_PATTERN.fullmatch(self.sha256):
            raise ConfigurationError("Artifact sha256 must contain exactly 64 hexadecimal digits.")
        if (
            isinstance(self.byte_size, bool)
            or not isinstance(self.byte_size, int)
            or self.byte_size < 0
        ):
            raise ConfigurationError("Artifact byte_size must be a non-negative integer.")


@dataclass(frozen=True, init=False)
class RunManifest:
    """Resolved parameters and artifacts for a reproducible execution."""

    run_id: str
    command: tuple[str, ...]
    resolved_config: FrozenDict
    seed: int | None
    seed_derivation: str | None
    inputs: tuple[ArtifactRef, ...]
    outputs: tuple[ArtifactRef, ...]
    provenance: Provenance
    started_at: str
    finished_at: str | None
    status: str
    issues: tuple[Issue, ...]

    def __init__(
        self,
        run_id: str,
        command: Iterable[str],
        resolved_config: Mapping[str, object],
        provenance: Provenance,
        started_at: str,
        status: str,
        *,
        seed: int | None = None,
        seed_derivation: str | None = None,
        inputs: Iterable[ArtifactRef] = (),
        outputs: Iterable[ArtifactRef] = (),
        finished_at: str | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        for field_name, value in (
            ("run_id", run_id),
            ("started_at", started_at),
            ("status", status),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"RunManifest {field_name} must be non-empty.")
        _non_empty_optional(seed_derivation, "RunManifest seed_derivation")
        _non_empty_optional(finished_at, "RunManifest finished_at")
        if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
            raise ConfigurationError("RunManifest seed must be an integer or None.")
        command_tuple = tuple(command)
        if not command_tuple or any(
            not isinstance(argument, str) or not argument for argument in command_tuple
        ):
            raise ConfigurationError("RunManifest command must contain non-empty strings.")
        input_tuple = tuple(inputs)
        output_tuple = tuple(outputs)
        issue_tuple = tuple(issues)
        if any(not isinstance(item, ArtifactRef) for item in (*input_tuple, *output_tuple)):
            raise ConfigurationError("RunManifest artifacts must be ArtifactRef objects.")
        if any(not isinstance(issue, Issue) for issue in issue_tuple):
            raise ConfigurationError("RunManifest issues must be Issue objects.")
        if not isinstance(provenance, Provenance):
            raise ConfigurationError("RunManifest provenance must be Provenance.")

        object.__setattr__(self, "run_id", run_id)
        object.__setattr__(self, "command", command_tuple)
        object.__setattr__(self, "resolved_config", freeze_mapping(resolved_config))
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "seed_derivation", seed_derivation)
        object.__setattr__(self, "inputs", input_tuple)
        object.__setattr__(self, "outputs", output_tuple)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "started_at", started_at)
        object.__setattr__(self, "finished_at", finished_at)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "issues", issue_tuple)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = [
    "ArtifactRef",
    "Citation",
    "ImplementationInfo",
    "Provenance",
    "ReferenceInfo",
    "RunManifest",
]
