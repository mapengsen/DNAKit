"""Immutable configuration and manifests for public data downloads."""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import Provenance
from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.exceptions import ConfigurationError

_MD5 = re.compile(r"^[0-9a-fA-F]{32}$")
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


@dataclass(frozen=True, slots=True)
class DownloadConfig:
    """Resource limits and credentials for bounded public downloads."""

    timeout: float = 60.0
    chunk_size: int = 1_048_576
    max_file_bytes: int = 20_000_000_000
    max_total_bytes: int = 50_000_000_000
    max_files: int = 100
    overwrite: bool = False
    api_key: str | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.timeout, bool)
            or not isinstance(self.timeout, (int, float))
            or not math.isfinite(self.timeout)
            or not 0 < self.timeout <= 3_600
        ):
            raise ConfigurationError(
                "DownloadConfig timeout must be in (0, 3600].",
                code="INVALID_DOWNLOAD_CONFIG",
            )
        for name, value, maximum in (
            ("chunk_size", self.chunk_size, 64 * 1024 * 1024),
            ("max_file_bytes", self.max_file_bytes, 10**12),
            ("max_total_bytes", self.max_total_bytes, 10**12),
            ("max_files", self.max_files, 10_000),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(
                    f"DownloadConfig {name} must be in [1, {maximum}].",
                    code="INVALID_DOWNLOAD_CONFIG",
                )
        if self.max_total_bytes < self.max_file_bytes:
            raise ConfigurationError(
                "DownloadConfig max_total_bytes cannot be smaller than max_file_bytes.",
                code="INVALID_DOWNLOAD_CONFIG",
            )
        if not isinstance(self.overwrite, bool):
            raise ConfigurationError(
                "DownloadConfig overwrite must be boolean.", code="INVALID_DOWNLOAD_CONFIG"
            )
        if self.api_key is not None and (
            not isinstance(self.api_key, str) or not self.api_key.strip()
        ):
            raise ConfigurationError(
                "DownloadConfig api_key must be None or non-empty text.",
                code="INVALID_DOWNLOAD_CONFIG",
            )


@dataclass(frozen=True, slots=True)
class DownloadProgress:
    """One progress event for a streamed public file download."""

    url: str
    target: str
    file_index: int
    file_count: int
    bytes_completed: int
    total_bytes: int | None = None


@dataclass(frozen=True, slots=True)
class RemoteFile:
    """One explicit HTTPS resource and optional integrity expectations."""

    url: str
    filename: str | None = None
    expected_md5: str | None = None
    expected_sha256: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.url, str) or not self.url.strip():
            raise ConfigurationError("RemoteFile url must be non-empty text.")
        if self.filename is not None and (
            not isinstance(self.filename, str) or not self.filename.strip()
        ):
            raise ConfigurationError("RemoteFile filename must be None or non-empty text.")
        if self.expected_md5 is not None and _MD5.fullmatch(self.expected_md5) is None:
            raise ConfigurationError("RemoteFile expected_md5 must contain 32 hexadecimal digits.")
        if self.expected_sha256 is not None and _SHA256.fullmatch(self.expected_sha256) is None:
            raise ConfigurationError(
                "RemoteFile expected_sha256 must contain 64 hexadecimal digits."
            )


@dataclass(frozen=True, slots=True)
class DownloadedFile:
    """Integrity metadata for one atomically installed file."""

    url: str
    path: str
    byte_size: int
    md5: str
    sha256: str
    checksum_verified: bool

    def __post_init__(self) -> None:
        if (
            isinstance(self.byte_size, bool)
            or not isinstance(self.byte_size, int)
            or self.byte_size < 0
        ):
            raise ConfigurationError("DownloadedFile byte_size is invalid.")
        if _MD5.fullmatch(self.md5) is None or _SHA256.fullmatch(self.sha256) is None:
            raise ConfigurationError("DownloadedFile checksums are invalid.")


@dataclass(frozen=True, init=False)
class DatasetDownloadResult:
    """Auditable files and manifest produced by one logical download."""

    kind: str
    source: str
    output_directory: str
    files: tuple[DownloadedFile, ...]
    manifest_path: str
    metadata: FrozenDict
    provenance: Provenance

    def __init__(
        self,
        kind: str,
        source: str,
        output_directory: str,
        files: Iterable[DownloadedFile],
        manifest_path: str,
        provenance: Provenance,
        *,
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        for name, value in (
            ("kind", kind),
            ("source", source),
            ("output_directory", output_directory),
            ("manifest_path", manifest_path),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"DatasetDownloadResult {name} must be non-empty text.")
        resolved_files = tuple(files)
        if any(not isinstance(item, DownloadedFile) for item in resolved_files):
            raise ConfigurationError("DatasetDownloadResult files must be DownloadedFile objects.")
        if not isinstance(provenance, Provenance):
            raise ConfigurationError("DatasetDownloadResult provenance must be Provenance.")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "output_directory", output_directory)
        object.__setattr__(self, "files", resolved_files)
        object.__setattr__(self, "manifest_path", manifest_path)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))
        object.__setattr__(self, "provenance", provenance)

    @property
    def downloaded_bytes(self) -> int:
        return sum(item.byte_size for item in self.files)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class IndexArtifact:
    """One locally generated sequence-index artifact."""

    path: str
    byte_size: int
    sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.path, str) or not self.path.strip():
            raise ConfigurationError("IndexArtifact path must be non-empty text.")
        if (
            isinstance(self.byte_size, bool)
            or not isinstance(self.byte_size, int)
            or self.byte_size < 0
        ):
            raise ConfigurationError("IndexArtifact byte_size is invalid.")
        if _SHA256.fullmatch(self.sha256) is None:
            raise ConfigurationError("IndexArtifact sha256 is invalid.")


@dataclass(frozen=True, slots=True)
class IndexBuildResult:
    """Auditable artifacts produced by an explicit external index tool."""

    tool: str
    fasta_path: str
    output_prefix: str
    artifacts: tuple[IndexArtifact, ...]
    command_output: str
    elapsed_seconds: float
    manifest_path: str
    provenance: Provenance

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = [
    "DatasetDownloadResult",
    "DownloadConfig",
    "DownloadProgress",
    "DownloadedFile",
    "IndexArtifact",
    "IndexBuildResult",
    "RemoteFile",
]
