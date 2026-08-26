"""Bounded, checksummed, atomic HTTPS downloads."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import replace
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlsplit
from urllib.request import Request, urlopen

from dnakit._version import __version__
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search._http import redact_url
from dnakit.search._shared import adapter_provenance

from .models import (
    DatasetDownloadResult,
    DownloadConfig,
    DownloadedFile,
    DownloadProgress,
    RemoteFile,
)

ProgressCallback = Callable[[DownloadProgress], None]
_SAFE_FILENAME = re.compile(r"^[^/\\\x00]{1,255}$")


def resolved_config(config: DownloadConfig | None) -> DownloadConfig:
    if config is None:
        return DownloadConfig()
    if not isinstance(config, DownloadConfig):
        raise TypeError("config must be DownloadConfig or None.")
    return config


def _https_url(url: str) -> str:
    if not isinstance(url, str) or not url.strip():
        raise ConfigurationError("Download URL must be non-empty text.")
    value = url.strip()
    parsed = urlsplit(value)
    if (
        parsed.scheme != "https"
        or not parsed.netloc
        or parsed.username
        or parsed.password
        or parsed.fragment
    ):
        raise ConfigurationError(
            "Downloads require an HTTPS URL without credentials or fragments.",
            code="INVALID_DOWNLOAD_URL",
        )
    return value


def _filename(resource: RemoteFile) -> str:
    raw = resource.filename
    if raw is None:
        raw = unquote(PurePosixPath(urlsplit(resource.url).path).name)
    if not isinstance(raw, str) or _SAFE_FILENAME.fullmatch(raw) is None or raw in {".", ".."}:
        raise ConfigurationError("Download filename is unsafe.", code="INVALID_DOWNLOAD_FILENAME")
    return raw


def _content_length(response: Any) -> int | None:
    headers = getattr(response, "headers", None)
    raw = headers.get("Content-Length") if headers is not None else None
    if not isinstance(raw, (str, bytes, int)) or isinstance(raw, bool):
        return None
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return None
    return value if value >= 0 else None


def download_file(
    resource: RemoteFile,
    target: str | os.PathLike[str],
    *,
    config: DownloadConfig | None = None,
    progress: ProgressCallback | None = None,
    file_index: int = 1,
    file_count: int = 1,
    headers: Mapping[str, str] | None = None,
    max_bytes_override: int | None = None,
) -> DownloadedFile:
    """Stream one HTTPS resource to an atomic, checksummed local file."""

    if not isinstance(resource, RemoteFile):
        raise TypeError("resource must be RemoteFile.")
    resolved = resolved_config(config)
    url = _https_url(resource.url)
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.")
    if (
        isinstance(file_index, bool)
        or not isinstance(file_index, int)
        or isinstance(file_count, bool)
        or not isinstance(file_count, int)
        or not 1 <= file_index <= file_count
    ):
        raise ConfigurationError("file_index/file_count are invalid.")
    max_bytes = resolved.max_file_bytes
    if max_bytes_override is not None:
        if (
            isinstance(max_bytes_override, bool)
            or not isinstance(max_bytes_override, int)
            or max_bytes_override < 1
        ):
            raise ConfigurationError("max_bytes_override must be positive.")
        max_bytes = min(max_bytes, max_bytes_override)
    destination = Path(target).expanduser().resolve()
    if destination.exists() and not resolved.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing download: {destination}")
    if destination.exists() and not destination.is_file():
        raise IsADirectoryError(f"Download target must be a file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    request_headers = {
        "Accept": "*/*",
        "User-Agent": f"DNAKit/{__version__} public-data-downloader",
        **dict(headers or {}),
    }
    request = Request(url, headers=request_headers, method="GET")
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.dnakit-",
            suffix=".part",
            dir=str(destination.parent),
            delete=False,
        ) as output:
            temp_path = Path(output.name)
            try:
                response_context = urlopen(request, timeout=float(resolved.timeout))
                with response_context as response:
                    final_url = getattr(response, "geturl", lambda: url)()
                    if urlsplit(final_url).scheme != "https":
                        raise DownloadError(
                            "Download redirected to a non-HTTPS URL.",
                            code="INSECURE_DOWNLOAD_REDIRECT",
                        )
                    total = _content_length(response)
                    if total is not None and total > max_bytes:
                        raise DownloadError(
                            "Remote file exceeds the configured byte limit.",
                            code="DOWNLOAD_SIZE_LIMIT",
                            context={"content_length": total, "max_bytes": max_bytes},
                        )
                    md5 = hashlib.md5(usedforsecurity=False)
                    sha256 = hashlib.sha256()
                    completed = 0
                    while True:
                        chunk = response.read(resolved.chunk_size)
                        if not chunk:
                            break
                        completed += len(chunk)
                        if completed > max_bytes:
                            raise DownloadError(
                                "Remote file exceeds the configured byte limit.",
                                code="DOWNLOAD_SIZE_LIMIT",
                                context={"downloaded_bytes": completed, "max_bytes": max_bytes},
                            )
                        output.write(chunk)
                        md5.update(chunk)
                        sha256.update(chunk)
                        if progress is not None:
                            progress(
                                DownloadProgress(
                                    redact_url(url),
                                    str(destination),
                                    file_index,
                                    file_count,
                                    completed,
                                    total,
                                )
                            )
            except HTTPError as exc:
                raise DownloadError(
                    "Remote server returned an HTTP error.",
                    code="DOWNLOAD_HTTP_ERROR",
                    context={"status": exc.code, "url": redact_url(url)},
                ) from exc
            except DownloadError:
                raise
            except (OSError, TimeoutError, URLError) as exc:
                raise DownloadError(
                    "Could not download the remote file.",
                    code="DOWNLOAD_NETWORK_ERROR",
                    context={"url": redact_url(url)},
                ) from exc
        actual_md5 = md5.hexdigest()
        actual_sha256 = sha256.hexdigest()
        if resource.expected_md5 is not None and actual_md5 != resource.expected_md5.lower():
            raise DownloadError(
                "Downloaded file failed its MD5 check.",
                code="CHECKSUM_MISMATCH",
                context={
                    "algorithm": "md5",
                    "expected": resource.expected_md5.lower(),
                    "actual": actual_md5,
                },
            )
        if (
            resource.expected_sha256 is not None
            and actual_sha256 != resource.expected_sha256.lower()
        ):
            raise DownloadError(
                "Downloaded file failed its SHA-256 check.",
                code="CHECKSUM_MISMATCH",
                context={
                    "algorithm": "sha256",
                    "expected": resource.expected_sha256.lower(),
                    "actual": actual_sha256,
                },
            )
        if resolved.overwrite:
            os.replace(temp_path, destination)
        else:
            os.link(temp_path, destination)
            temp_path.unlink()
        temp_path = None
        return DownloadedFile(
            redact_url(url),
            str(destination),
            completed,
            actual_md5,
            actual_sha256,
            resource.expected_md5 is not None or resource.expected_sha256 is not None,
        )
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def _write_manifest(
    result_values: Mapping[str, object],
    path: Path,
    *,
    overwrite: bool,
) -> None:
    payload = json.dumps(result_values, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing manifest: {path}")
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.dnakit-",
        suffix=".part",
        dir=str(path.parent),
        delete=False,
    ) as output:
        output.write(payload)
        temp_path = Path(output.name)
    try:
        if overwrite:
            os.replace(temp_path, path)
        else:
            os.link(temp_path, path)
            temp_path.unlink()
    finally:
        temp_path.unlink(missing_ok=True)


def _install_staged_files(
    staged: Sequence[Path],
    targets: Sequence[Path],
    *,
    overwrite: bool,
) -> None:
    """Install same-filesystem staged files as one rollback-capable transaction."""

    sources = tuple(staged)
    destinations = tuple(targets)
    if (
        not sources
        or len(sources) != len(destinations)
        or len(set(destinations)) != len(destinations)
    ):
        raise ConfigurationError("Staged output paths collide.", code="OUTPUT_COLLISION")
    if any(not source.is_file() for source in sources):
        raise ConfigurationError("Every staged output must be an existing file.")
    if not overwrite and any(target.exists() for target in destinations):
        raise FileExistsError("Refusing to overwrite existing outputs.")
    if any(target.exists() and not target.is_file() for target in destinations):
        raise IsADirectoryError("Every existing output target must be a file.")

    backup_dir = sources[0].parent / ".dnakit-backups"
    backup_dir.mkdir()
    backups: list[tuple[Path, Path]] = []
    installed: list[Path] = []
    try:
        for index, (source, target) in enumerate(zip(sources, destinations, strict=True)):
            if target.exists():
                backup = backup_dir / str(index)
                os.replace(target, backup)
                backups.append((target, backup))
            os.replace(source, target)
            installed.append(target)
    except BaseException:
        for target in installed:
            target.unlink(missing_ok=True)
        for target, backup in reversed(backups):
            if backup.exists():
                os.replace(backup, target)
        raise
    finally:
        for _, backup in backups:
            backup.unlink(missing_ok=True)
        backup_dir.rmdir()


def dataset(
    resources: Iterable[RemoteFile],
    output_dir: str | os.PathLike[str],
    *,
    kind: str = "dataset",
    source: str = "explicit HTTPS resources",
    metadata: Mapping[str, object] | None = None,
    config: DownloadConfig | None = None,
    progress: ProgressCallback | None = None,
    headers: Mapping[str, str] | None = None,
) -> DatasetDownloadResult:
    """Download an explicit list of public files and write an integrity manifest."""

    resolved = resolved_config(config)
    values = tuple(resources)
    if not values or len(values) > resolved.max_files:
        raise ConfigurationError(
            "resources must contain between 1 and max_files entries.",
            code="DOWNLOAD_FILE_LIMIT",
        )
    if any(not isinstance(resource, RemoteFile) for resource in values):
        raise TypeError("resources must contain RemoteFile objects.")
    kind_value = kind.strip() if isinstance(kind, str) else ""
    source_value = source.strip() if isinstance(source, str) else ""
    if not kind_value or not source_value:
        raise ConfigurationError("kind and source must be non-empty text.")
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(f"output_dir must be a directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    names = tuple(_filename(resource) for resource in values)
    if len(set(names)) != len(names):
        raise ConfigurationError("Download filenames collide.", code="OUTPUT_COLLISION")
    manifest = output / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', kind_value)}_manifest.json"
    prospective = (*tuple(output / name for name in names), manifest)
    if not resolved.overwrite and any(path.exists() for path in prospective):
        raise FileExistsError("Refusing to overwrite existing dataset outputs.")
    downloaded: list[DownloadedFile] = []
    total_bytes = 0
    with tempfile.TemporaryDirectory(prefix=".dnakit-dataset-", dir=str(output)) as stage_raw:
        stage = Path(stage_raw)
        stage_config = replace(resolved, overwrite=False)
        staged_files: list[Path] = []
        final_files: list[Path] = []
        for index, (resource, name) in enumerate(zip(values, names, strict=True), start=1):
            remaining = resolved.max_total_bytes - total_bytes
            if remaining < 1:
                raise DownloadError(
                    "Dataset exceeds max_total_bytes.", code="DOWNLOAD_TOTAL_SIZE_LIMIT"
                )
            final_path = output / name

            def report(event: DownloadProgress, *, target: Path = final_path) -> None:
                if progress is not None:
                    progress(
                        DownloadProgress(
                            event.url,
                            str(target),
                            event.file_index,
                            event.file_count,
                            event.bytes_completed,
                            event.total_bytes,
                        )
                    )

            staged_path = stage / name
            item = download_file(
                resource,
                staged_path,
                config=stage_config,
                progress=report if progress is not None else None,
                file_index=index,
                file_count=len(values),
                headers=headers,
                max_bytes_override=remaining,
            )
            downloaded.append(
                DownloadedFile(
                    item.url,
                    str(final_path),
                    item.byte_size,
                    item.md5,
                    item.sha256,
                    item.checksum_verified,
                )
            )
            staged_files.append(staged_path)
            final_files.append(final_path)
            total_bytes += item.byte_size
        provenance = adapter_provenance(
            source_value,
            citation_url=redact_url(values[0].url),
            filters={"kind": kind_value, "resource_count": len(values)},
        )
        manifest_payload = {
            "kind": kind_value,
            "source": source_value,
            "output_directory": str(output),
            "files": [
                {
                    "url": item.url,
                    "path": item.path,
                    "byte_size": item.byte_size,
                    "md5": item.md5,
                    "sha256": item.sha256,
                    "checksum_verified": item.checksum_verified,
                }
                for item in downloaded
            ],
            "metadata": dict(metadata or {}),
            "provenance": provenance.to_dict(),
        }
        staged_manifest = stage / manifest.name
        _write_manifest(manifest_payload, staged_manifest, overwrite=False)
        _install_staged_files(
            (*staged_files, staged_manifest),
            (*final_files, manifest),
            overwrite=resolved.overwrite,
        )
    return DatasetDownloadResult(
        kind_value,
        source_value,
        str(output),
        downloaded,
        str(manifest),
        provenance,
        metadata=metadata,
    )


def tracks(
    resources: Sequence[RemoteFile],
    output_dir: str | os.PathLike[str],
    *,
    source: str = "public genomic tracks",
    config: DownloadConfig | None = None,
    progress: ProgressCallback | None = None,
) -> DatasetDownloadResult:
    """Download explicit BED/bigBed/bedGraph/bigWig/peak resources."""

    return dataset(
        resources,
        output_dir,
        kind="tracks",
        source=source,
        config=config,
        progress=progress,
    )


__all__ = ["ProgressCallback", "dataset", "download_file", "resolved_config", "tracks"]
