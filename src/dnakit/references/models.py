"""Value objects returned by reference-genome download adapters."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Literal, cast

from dnakit.core._json import to_json_compatible
from dnakit.core.provenance import Provenance
from dnakit.exceptions import ConfigurationError

DownloadPhase = Literal["resolve", "download", "extract"]
_MD5 = re.compile(r"^[0-9a-f]{32}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$", flags=re.IGNORECASE)


def _require_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"{name} must be a non-empty string.", code="INVALID_DOWNLOAD_RESULT"
        )
    return value


@dataclass(frozen=True, slots=True)
class GenomeAssembly:
    """Resolved NCBI genome assembly identity."""

    query: str
    accession: str
    organism: str | None = None
    assembly_name: str | None = None
    source_database: str | None = None

    def __post_init__(self) -> None:
        _require_text(self.query, "GenomeAssembly query")
        accession = _require_text(self.accession, "GenomeAssembly accession")
        if _ASSEMBLY_ACCESSION.fullmatch(accession) is None:
            raise ConfigurationError(
                "GenomeAssembly accession must be a versioned GCA_ or GCF_ assembly accession.",
                code="INVALID_ASSEMBLY_ACCESSION",
            )
        for name in ("organism", "assembly_name", "source_database"):
            value = getattr(self, name)
            if value is not None:
                _require_text(value, f"GenomeAssembly {name}")
        object.__setattr__(self, "accession", accession.upper())


@dataclass(frozen=True, slots=True)
class DownloadProgress:
    """One progress event emitted by a reference-data download."""

    phase: DownloadPhase
    query: str
    accession: str
    bytes_completed: int
    total_bytes: int | None = None
    current_file: str | None = None

    def __post_init__(self) -> None:
        if self.phase not in {"resolve", "download", "extract"}:
            raise ConfigurationError(
                "DownloadProgress phase is invalid.", code="INVALID_DOWNLOAD_PROGRESS"
            )
        _require_text(self.query, "DownloadProgress query")
        _require_text(self.accession, "DownloadProgress accession")
        for name, value in (
            ("bytes_completed", self.bytes_completed),
            ("total_bytes", self.total_bytes),
        ):
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int) or value < 0
            ):
                raise ConfigurationError(
                    f"DownloadProgress {name} must be a non-negative integer or None.",
                    code="INVALID_DOWNLOAD_PROGRESS",
                )
        if self.total_bytes is not None and self.bytes_completed > self.total_bytes:
            raise ConfigurationError(
                "DownloadProgress bytes_completed cannot exceed total_bytes.",
                code="INVALID_DOWNLOAD_PROGRESS",
            )
        if self.current_file is not None:
            _require_text(self.current_file, "DownloadProgress current_file")


@dataclass(frozen=True, slots=True)
class GenomeDownloadResult:
    """Auditable files and checksums produced by a genome download."""

    query: str
    accession: str
    organism: str | None
    assembly_name: str | None
    output_directory: str
    fasta_path: str
    metadata_path: str
    checksum_path: str
    package_path: str | None
    download_url: str
    downloaded_bytes: int
    package_sha256: str
    fasta_sha256: str
    fasta_md5: str
    provenance: Provenance

    def __post_init__(self) -> None:
        for name in (
            "query",
            "accession",
            "output_directory",
            "fasta_path",
            "metadata_path",
            "checksum_path",
            "download_url",
        ):
            _require_text(getattr(self, name), f"GenomeDownloadResult {name}")
        for name in ("organism", "assembly_name", "package_path"):
            value = getattr(self, name)
            if value is not None:
                _require_text(value, f"GenomeDownloadResult {name}")
        for name in ("package_sha256", "fasta_sha256", "fasta_md5"):
            value = _require_text(getattr(self, name), f"GenomeDownloadResult {name}").lower()
            pattern = _MD5 if name == "fasta_md5" else _SHA256
            if pattern.fullmatch(value) is None:
                raise ConfigurationError(
                    f"GenomeDownloadResult {name} must be a lowercase hexadecimal hash.",
                    code="INVALID_DOWNLOAD_RESULT",
                )
            object.__setattr__(self, name, value)
        if (
            isinstance(self.downloaded_bytes, bool)
            or not isinstance(self.downloaded_bytes, int)
            or self.downloaded_bytes < 1
        ):
            raise ConfigurationError(
                "GenomeDownloadResult downloaded_bytes must be positive.",
                code="INVALID_DOWNLOAD_RESULT",
            )
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError(
                "GenomeDownloadResult provenance must be Provenance.",
                code="INVALID_DOWNLOAD_RESULT",
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible download manifest."""

        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["DownloadPhase", "DownloadProgress", "GenomeAssembly", "GenomeDownloadResult"]
