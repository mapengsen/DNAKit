"""Shared validation and atomic persistence for accession sequence downloads."""

from __future__ import annotations

import hashlib
import io
import re
import tempfile
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search._http import redact_url, require_text
from dnakit.search._shared import adapter_provenance

from .files import _install_staged_files, _write_manifest
from .models import DatasetDownloadResult, DownloadConfig, DownloadedFile


def sequence_accession(value: str) -> str:
    """Accept one accession, preserving an explicitly requested version."""
    accession = require_text(value, "accession", max_length=128)
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*[0-9](?:\.[0-9]+)?", accession) is None:
        raise ConfigurationError(
            "Provide one sequence accession, optionally versioned.", code="INVALID_ACCESSION"
        )
    if accession.upper().startswith(("GCA_", "GCF_")):
        raise ConfigurationError(
            "Assembly accessions require genome_package().", code="INVALID_ACCESSION"
        )
    return accession


def sequence_target(output_path: str | Path, config: DownloadConfig) -> tuple[Path, Path]:
    raw = Path(output_path).expanduser().absolute()
    manifest = raw.with_name(raw.name + ".manifest.json")
    if raw.is_symlink() or manifest.is_symlink():
        raise ConfigurationError("Sequence output and manifest must not be symlinks.")
    target = raw.resolve()
    manifest = target.with_name(target.name + ".manifest.json")
    for path in (target, manifest):
        if path.exists():
            if not path.is_file():
                raise IsADirectoryError(str(path))
            if not config.overwrite:
                raise FileExistsError(f"Refusing to overwrite {path}")
    return target, manifest


def format_converter(source: str, target: str) -> Any:
    """Check the optional scientific format library before requesting data."""
    if source == target:
        return None
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ConfigurationError(
            "Sequence format conversion requires Biopython: pip install 'dnakit[bio-formats]'.",
            code="MISSING_OPTIONAL_DEPENDENCY",
        ) from exc
    return SeqIO


def validate_sequence(text: str, format: str) -> None:
    """Reject empty responses, error documents, and incomplete flat files."""
    stripped = text.strip()
    valid = False
    if format == "fasta" and stripped.startswith(">"):
        has_sequence = False
        valid = True
        for index, line in enumerate(stripped.splitlines()):
            if line.startswith(">"):
                if (index and not has_sequence) or not line[1:].strip():
                    valid = False
                    break
                has_sequence = False
            elif line.strip():
                if re.fullmatch(r"[A-Za-z*?.-]+", line.strip()) is None:
                    valid = False
                    break
                has_sequence = True
        valid = valid and has_sequence
    elif format == "genbank":
        valid = stripped.startswith("LOCUS ") and stripped.endswith("//") and "\nORIGIN" in text
    elif format == "embl":
        valid = stripped.startswith("ID ") and stripped.endswith("//") and "\nSQ " in text
    if not valid:
        raise DownloadError(
            "Provider returned an empty, incomplete, or invalid sequence file.",
            code="SEQUENCE_RESPONSE_ERROR",
        )


def sequence_identity(text: str, format: str, accession: str) -> tuple[str, int]:
    """Require one complete record with the requested accession and version."""
    if format == "fasta":
        headers = [line for line in text.splitlines() if line.startswith(">")]
        candidates = (
            re.findall(r"[A-Za-z][A-Za-z0-9_]*\.[0-9]+", headers[0].split()[0])
            if len(headers) == 1
            else []
        )
        length = sum(len(line.strip()) for line in text.splitlines() if not line.startswith(">"))
        expected_length = None
    elif format == "genbank":
        candidates = re.findall(r"^VERSION\s+(\S+)", text, re.M)
        origin = text.split("\nORIGIN", 1)[1].split("//", 1)[0]
        length = len(re.sub(r"[^A-Za-z]", "", origin))
        size = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+(?:bp|aa)\b", text, re.M)
        expected_length = int(size[1]) if size else None
    else:
        identifiers = re.findall(r"^ID\s+(\S+);\s+SV\s+(\d+);", text, re.M)
        candidates = [f"{name}.{version}" for name, version in identifiers]
        symbols = text.split("\nSQ ", 1)[1].split("\n", 1)[1].split("//", 1)[0]
        length = len(re.sub(r"[^A-Za-z]", "", symbols))
        size = re.search(r"^SQ\s+Sequence\s+(\d+)\s+BP", text, re.M)
        expected_length = int(size[1]) if size else None
    matched = [
        value
        for value in candidates
        if (
            value.casefold() == accession.casefold()
            if "." in accession
            else value.split(".", 1)[0].casefold() == accession.casefold()
        )
    ]
    if not matched or len(set(candidates)) != 1:
        raise DownloadError(
            "Sequence response does not match the requested accession/version.",
            code="SEQUENCE_IDENTITY_MISMATCH",
        )
    if not length or (expected_length is not None and length != expected_length):
        raise DownloadError(
            "Sequence symbols do not match the annotated record length.",
            code="SEQUENCE_LENGTH_MISMATCH",
        )
    return matched[0], length


def save_sequence(
    text: str,
    *,
    accession: str,
    database: str,
    provider: str,
    url: str,
    source_format: str,
    format: str,
    converter: Any,
    target: Path,
    manifest: Path,
    config: DownloadConfig,
    metadata: dict[str, object] | None = None,
) -> DatasetDownloadResult:
    validate_sequence(text, source_format)
    resolved_accession, sequence_length = sequence_identity(text, source_format, accession)
    if converter is not None:
        output = io.StringIO()
        try:
            count = converter.convert(io.StringIO(text), source_format, output, format)
            if count != 1:
                raise ValueError("Expected one sequence record")
            text = output.getvalue()
        except (ValueError, TypeError) as exc:
            raise DownloadError(
                "Could not convert the annotated sequence format.", code="SEQUENCE_CONVERSION_ERROR"
            ) from exc
        validate_sequence(text, format)
        converted_accession, converted_length = sequence_identity(text, format, accession)
        if (converted_accession, converted_length) != (resolved_accession, sequence_length):
            raise DownloadError(
                "Format conversion changed sequence identity or length.",
                code="SEQUENCE_CONVERSION_ERROR",
            )
    payload = text.encode("utf-8")
    if len(payload) > min(config.max_file_bytes, config.max_total_bytes):
        raise DownloadError(
            "Sequence output exceeds configured byte limits.", code="DOWNLOAD_SIZE_LIMIT"
        )
    safe_url = redact_url(url)
    details = {
        "accession": accession,
        "resolved_accession": resolved_accession,
        "sequence_length": sequence_length,
        "database": database,
        "format": format,
        "source_format": source_format,
        "converted": source_format != format,
        **(metadata or {}),
    }
    provenance = adapter_provenance(
        provider,
        citation_url=(
            "https://www.ncbi.nlm.nih.gov/books/NBK25499/"
            if provider == "NCBI Entrez"
            else "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html"
        ),
        filters=details,
    )
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    sha256 = hashlib.sha256(payload).hexdigest()
    downloaded = DownloadedFile(safe_url, str(target), len(payload), md5, sha256, False)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".dnakit-accession-", dir=target.parent) as raw:
        stage = Path(raw)
        staged_file, staged_manifest = stage / target.name, stage / manifest.name
        staged_file.write_bytes(payload)
        _write_manifest(
            {
                "kind": "accession_sequence",
                "source": provider,
                "request_url": safe_url,
                "metadata": details,
                "provenance": provenance.to_dict(),
                "file": {
                    "path": str(target),
                    "byte_size": len(payload),
                    "md5": md5,
                    "sha256": sha256,
                },
            },
            staged_manifest,
            overwrite=False,
        )
        _install_staged_files(
            (staged_file, staged_manifest), (target, manifest), overwrite=config.overwrite
        )
    return DatasetDownloadResult(
        "accession_sequence",
        provider,
        str(target.parent),
        (downloaded,),
        str(manifest),
        provenance,
        metadata=details,
    )
