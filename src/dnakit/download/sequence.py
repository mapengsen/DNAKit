"""Persist coordinate sequences returned by Ensembl as auditable FASTA."""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
from collections.abc import Callable, Sequence
from pathlib import Path

from dnakit.exceptions import DownloadError
from dnakit.search import QueryProgress, SearchConfig
from dnakit.search.ensembl import DEFAULT_ENSEMBL_REST_URL
from dnakit.search.ensembl import sequence as search_sequence

from .files import _install_staged_files, _write_manifest, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadedFile


def _fasta_identifier(value: object, fallback: str) -> str:
    raw = value if isinstance(value, str) and value.strip() else fallback
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(raw))[:255]


def sequence(
    species: str,
    region: str | Sequence[str],
    output_path: str | os.PathLike[str],
    *,
    strand: int = 1,
    upstream: int = 0,
    downstream: int = 0,
    mask: str | None = None,
    config: DownloadConfig | None = None,
    query_progress: Callable[[QueryProgress], None] | None = None,
    ensembl_api_base_url: str = DEFAULT_ENSEMBL_REST_URL,
) -> DatasetDownloadResult:
    """Fetch one or more 0-based half-open regions and atomically write FASTA."""

    resolved = resolved_config(config)
    target = Path(output_path).expanduser().resolve()
    if target.exists() and not target.is_file():
        raise IsADirectoryError(f"output_path must be a file: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    manifest = target.with_name(f"{target.name}.manifest.json")
    if manifest.exists() and not manifest.is_file():
        raise IsADirectoryError(f"Sequence manifest target must be a file: {manifest}")
    if not resolved.overwrite and (target.exists() or manifest.exists()):
        raise FileExistsError("Refusing to overwrite existing FASTA or manifest.")
    result = search_sequence(
        species,
        region,
        strand=strand,  # type: ignore[arg-type]
        upstream=upstream,
        downstream=downstream,
        mask=mask,  # type: ignore[arg-type]
        progress=query_progress,
        config=SearchConfig(
            timeout=min(float(resolved.timeout), 300.0),
            max_response_bytes=min(resolved.max_file_bytes, 100_000_000),
            max_records=resolved.max_files,
        ),
        api_base_url=ensembl_api_base_url,
    )
    lines: list[str] = []
    for index, record in enumerate(result.records, start=1):
        symbols = record.get("seq")
        if not isinstance(symbols, str) or not symbols:
            raise DownloadError(
                "Ensembl sequence response did not contain sequence text.",
                code="QUERY_RESPONSE_ERROR",
            )
        identifier = _fasta_identifier(record.get("id"), f"region_{index}")
        requested = record.get("requested_region")
        description = f" requested={requested}" if isinstance(requested, str) else ""
        lines.append(f">{identifier}{description}")
        lines.extend(symbols[offset : offset + 80] for offset in range(0, len(symbols), 80))
    payload = ("\n".join(lines) + "\n").encode("ascii")
    if len(payload) > resolved.max_file_bytes:
        raise DownloadError("FASTA output exceeds max_file_bytes.", code="DOWNLOAD_SIZE_LIMIT")
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    sha256 = hashlib.sha256(payload).hexdigest()
    downloaded = DownloadedFile(
        result.request_url,
        str(target),
        len(payload),
        md5,
        sha256,
        False,
    )
    with tempfile.TemporaryDirectory(prefix=".dnakit-sequence-", dir=str(target.parent)) as raw:
        stage = Path(raw)
        staged_target = stage / target.name
        staged_target.write_bytes(payload)
        staged_manifest = stage / manifest.name
        _write_manifest(
            {
                "kind": "coordinate_sequence",
                "source": "Ensembl",
                "query": result.to_dict(),
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
            (staged_target, staged_manifest),
            (target, manifest),
            overwrite=resolved.overwrite,
        )
    return DatasetDownloadResult(
        "coordinate_sequence",
        "Ensembl",
        str(target.parent),
        (downloaded,),
        str(manifest),
        result.provenance,
        metadata={
            "species": species,
            "region_count": len(result.records),
            "input_coordinate_system": "0-based half-open",
            "strand": strand,
        },
    )


__all__ = ["sequence"]
