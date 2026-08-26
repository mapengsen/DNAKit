"""Download public files discovered through ENCODE Portal search results."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from urllib.parse import quote, urlsplit

from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import QueryResult
from dnakit.search._http import require_https_base

from .files import dataset, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadProgress, RemoteFile

DEFAULT_ENCODE_DOWNLOAD_URL = "https://www.encodeproject.org"
_MD5 = re.compile(r"^[0-9a-fA-F]{32}$")


def _url(value: str, base_url: str) -> str:
    raw = value.strip()
    parsed = urlsplit(raw)
    if parsed.scheme:
        if parsed.scheme != "https" or not parsed.netloc or parsed.fragment:
            raise ConfigurationError("ENCODE file URL must use HTTPS.")
        return raw
    relative = raw.lstrip("/")
    if not relative or ".." in PurePosixPath(relative).parts:
        raise ConfigurationError("ENCODE file path is unsafe.")
    return f"{require_https_base(base_url)}/{quote(relative, safe='/._-@')}"


def files(
    query: QueryResult,
    output_dir: str | Path,
    *,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    download_base_url: str = DEFAULT_ENCODE_DOWNLOAD_URL,
) -> DatasetDownloadResult:
    """Download public ENCODE file records, verifying advertised MD5 values when present."""

    if not isinstance(query, QueryResult) or query.provider != "ENCODE":
        raise ConfigurationError("query must be a dnakit.search.encode_search() result.")
    resolved = resolved_config(config)
    resources: list[RemoteFile] = []
    expected_total = 0
    selected: list[dict[str, object]] = []
    for record in query.records:
        href = record.get("href")
        if not isinstance(href, str) or not href.strip():
            continue
        url = _url(href, download_base_url)
        filename = PurePosixPath(urlsplit(url).path).name
        if not filename:
            raise DownloadError("ENCODE file URL has no filename.", code="QUERY_RESPONSE_ERROR")
        raw_md5 = record.get("md5sum")
        expected_md5 = raw_md5 if isinstance(raw_md5, str) and _MD5.fullmatch(raw_md5) else None
        raw_size = record.get("file_size")
        if isinstance(raw_size, int) and not isinstance(raw_size, bool) and raw_size >= 0:
            if raw_size > resolved.max_file_bytes:
                raise DownloadError(
                    "ENCODE file exceeds max_file_bytes.", code="DOWNLOAD_SIZE_LIMIT"
                )
            expected_total += raw_size
        resources.append(RemoteFile(url, filename=filename, expected_md5=expected_md5))
        selected.append(
            {
                "accession": record.get("accession"),
                "file_format": record.get("file_format"),
                "output_type": record.get("output_type"),
                "href": href,
                "provider_md5_available": expected_md5 is not None,
            }
        )
    if not resources:
        raise DownloadError(
            "ENCODE query contains no downloadable file records; search object_type='File'.",
            code="ENCODE_FILES_NOT_FOUND",
        )
    if len(resources) > resolved.max_files:
        raise DownloadError("ENCODE files exceed max_files.", code="DOWNLOAD_FILE_LIMIT")
    if expected_total > resolved.max_total_bytes:
        raise DownloadError(
            "ENCODE files exceed max_total_bytes.", code="DOWNLOAD_TOTAL_SIZE_LIMIT"
        )
    return dataset(
        resources,
        output_dir,
        kind="encode_files",
        source="ENCODE Portal",
        metadata={
            "query_url": query.request_url,
            "expected_total_bytes": expected_total,
            "selected_files": selected,
        },
        config=resolved,
        progress=progress,
    )


__all__ = ["DEFAULT_ENCODE_DOWNLOAD_URL", "files"]
