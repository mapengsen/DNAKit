"""Download files selected from the UCSC Genome Browser file catalog."""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from urllib.parse import quote, urlsplit

from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import QueryResult
from dnakit.search._http import require_https_base

from .files import dataset, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadProgress, RemoteFile

DEFAULT_UCSC_DOWNLOAD_URL = "https://hgdownload.soe.ucsc.edu"


def _file_url(value: str, base_url: str) -> tuple[str, str]:
    raw = value.strip()
    parsed = urlsplit(raw)
    if parsed.scheme:
        if parsed.scheme != "https" or not parsed.netloc or parsed.fragment:
            raise ConfigurationError("UCSC file URL must use HTTPS.", code="INVALID_DOWNLOAD_URL")
        path = parsed.path
        url = raw
    else:
        relative = raw.lstrip("/")
        path = "/" + relative
        url = f"{require_https_base(base_url)}/{quote(relative, safe='/._-')}"
    parts = PurePosixPath(path).parts
    if not parts or ".." in parts:
        raise ConfigurationError("UCSC file path is unsafe.", code="INVALID_DOWNLOAD_URL")
    basename = PurePosixPath(path).name
    if not basename:
        raise ConfigurationError("UCSC file path has no filename.", code="INVALID_DOWNLOAD_URL")
    digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:12]
    maximum_basename = 240 - len(digest)
    filename = f"{digest}_{basename[:maximum_basename]}"
    return url, filename


def files(
    query: QueryResult,
    output_dir: str | Path,
    *,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    download_base_url: str = DEFAULT_UCSC_DOWNLOAD_URL,
) -> DatasetDownloadResult:
    """Download a bounded ``dnakit.search.ucsc_files`` result with local checksums."""

    if (
        not isinstance(query, QueryResult)
        or query.provider != "UCSC"
        or query.query_type != "files"
    ):
        raise ConfigurationError("query must be a dnakit.search.ucsc_files() result.")
    resolved = resolved_config(config)
    resources: list[RemoteFile] = []
    expected_total = 0
    source_paths: dict[str, str] = {}
    for record in query.records:
        raw_url = record.get("url")
        if not isinstance(raw_url, str) or not raw_url.strip():
            raise DownloadError("UCSC file record has no URL.", code="QUERY_RESPONSE_ERROR")
        url, filename = _file_url(raw_url, download_base_url)
        raw_size = record.get("sizeBytes")
        if isinstance(raw_size, int) and not isinstance(raw_size, bool) and raw_size >= 0:
            if raw_size > resolved.max_file_bytes:
                raise DownloadError("UCSC file exceeds max_file_bytes.", code="DOWNLOAD_SIZE_LIMIT")
            expected_total += raw_size
        resources.append(RemoteFile(url, filename=filename))
        source_paths[filename] = raw_url
    if not resources:
        raise DownloadError(
            "UCSC query selected no downloadable files.", code="UCSC_FILES_NOT_FOUND"
        )
    if len(resources) > resolved.max_files:
        raise DownloadError("UCSC files exceed max_files.", code="DOWNLOAD_FILE_LIMIT")
    if expected_total > resolved.max_total_bytes:
        raise DownloadError("UCSC files exceed max_total_bytes.", code="DOWNLOAD_TOTAL_SIZE_LIMIT")
    return dataset(
        resources,
        output_dir,
        kind="ucsc_files",
        source="UCSC Genome Browser Downloads",
        metadata={
            "catalog_url": query.request_url,
            "expected_total_bytes": expected_total,
            "source_paths": source_paths,
        },
        config=resolved,
        progress=progress,
    )


__all__ = ["DEFAULT_UCSC_DOWNLOAD_URL", "files"]
