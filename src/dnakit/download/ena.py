"""Download public read files discovered through ENA metadata."""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import QueryResult, SearchConfig
from dnakit.search.ena import DEFAULT_ENA_PORTAL_URL, DEFAULT_READ_FIELDS
from dnakit.search.ena import reads as search_reads

from .files import dataset, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadProgress, RemoteFile

_RUN_ACCESSION = re.compile(r"^[SED]RR\d+$", flags=re.IGNORECASE)


def _ena_https(value: str) -> str:
    item = value.strip()
    if item.startswith("ftp://"):
        return "https://" + item[len("ftp://") :]
    if item.startswith("https://"):
        return item
    if "://" not in item:
        return "https://" + item.lstrip("/")
    raise ConfigurationError(
        "ENA file URL uses an unsupported scheme.", code="INVALID_DOWNLOAD_URL"
    )


def _parts(record: object, key: str) -> tuple[str, ...]:
    if not isinstance(record, Mapping):
        return ()
    raw = record.get(key)
    if not isinstance(raw, str) or not raw.strip():
        return ()
    return tuple(item.strip() for item in raw.split(";") if item.strip())


def reads(
    query: str | QueryResult,
    output_dir: str | Path,
    *,
    file_kind: str = "fastq",
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    ena_api_base_url: str = DEFAULT_ENA_PORTAL_URL,
) -> DatasetDownloadResult:
    """Download public ENA FASTQ, submitted, or SRA files with advertised MD5 checks."""

    resolved = resolved_config(config)
    if file_kind not in {"fastq", "submitted", "sra"}:
        raise ConfigurationError("file_kind must be fastq, submitted, or sra.")
    if isinstance(query, QueryResult):
        result = query
        if result.provider != "ENA" or result.query_type != "reads":
            raise ConfigurationError("QueryResult must come from dnakit.search.reads().")
    elif isinstance(query, str):
        expression = (
            f'run_accession="{query}"' if _RUN_ACCESSION.fullmatch(query.strip()) else query
        )
        fields = tuple(DEFAULT_READ_FIELDS)
        if file_kind == "sra":
            fields += ("sra_ftp", "sra_md5", "sra_bytes")
        result = search_reads(
            expression,
            fields=fields,
            limit=resolved.max_files,
            config=SearchConfig(
                timeout=min(float(resolved.timeout), 300.0),
                max_response_bytes=min(resolved.max_file_bytes, 100_000_000),
                max_records=resolved.max_files,
            ),
            api_base_url=ena_api_base_url,
        )
    else:
        raise TypeError("query must be ENA query text, run accession, or QueryResult.")
    url_key = f"{file_kind}_ftp"
    md5_key = f"{file_kind}_md5"
    bytes_key = f"{file_kind}_bytes"
    resources: list[RemoteFile] = []
    expected_sizes: dict[str, int | None] = {}
    for record in result.records:
        urls = _parts(record, url_key)
        md5s = _parts(record, md5_key)
        sizes = _parts(record, bytes_key)
        if md5s and len(md5s) != len(urls):
            raise DownloadError(
                "ENA returned mismatched file and MD5 lists.", code="ENA_METADATA_MISMATCH"
            )
        if sizes and len(sizes) != len(urls):
            raise DownloadError(
                "ENA returned mismatched file and size lists.", code="ENA_METADATA_MISMATCH"
            )
        for index, raw_url in enumerate(urls):
            url = _ena_https(raw_url)
            filename = PurePosixPath(urlsplit(url).path).name
            checksum = md5s[index] if md5s else None
            resources.append(RemoteFile(url, filename=filename, expected_md5=checksum))
            try:
                expected_sizes[filename] = int(sizes[index]) if sizes else None
            except ValueError as exc:
                raise DownloadError(
                    "ENA returned an invalid file size.", code="ENA_METADATA_MISMATCH"
                ) from exc
    if not resources:
        raise DownloadError(
            f"ENA query returned no public {file_kind} file URLs.",
            code="ENA_FILES_NOT_FOUND",
        )
    if len(resources) > resolved.max_files:
        raise DownloadError("ENA files exceed max_files.", code="DOWNLOAD_FILE_LIMIT")
    return dataset(
        resources,
        output_dir,
        kind=f"ena_{file_kind}",
        source="European Nucleotide Archive",
        metadata={
            "query_url": result.request_url,
            "file_kind": file_kind,
            "expected_sizes": expected_sizes,
        },
        config=resolved,
        progress=progress,
    )


__all__ = ["reads"]
