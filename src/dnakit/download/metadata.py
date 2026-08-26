"""Atomically export query metadata to common interchange formats."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Literal, cast

from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import QueryResult

from .files import _install_staged_files, _write_manifest, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadedFile

MetadataFormat = Literal["json", "jsonl", "csv", "tsv", "xml"]


def _format(value: str | None, path: Path) -> MetadataFormat:
    selected = path.suffix.lstrip(".").casefold() if value is None else value.casefold()
    if selected not in {"json", "jsonl", "csv", "tsv", "xml"}:
        raise ConfigurationError("Metadata format must be json, jsonl, csv, tsv, or xml.")
    return cast(MetadataFormat, selected)


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (str, int, float)):
        return str(value)
    return json.dumps(to_json_compatible(value), ensure_ascii=False, sort_keys=True)


def _tabular(result: QueryResult, *, delimiter: str) -> bytes:
    fields = sorted({str(key) for record in result.records for key in record})
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
    if fields:
        writer.writerow(fields)
        for record in result.records:
            writer.writerow(_cell(record.get(field)) for field in fields)
    return stream.getvalue().encode("utf-8")


def _xml(result: QueryResult) -> bytes:
    root = ET.Element(
        "dnakit-query",
        {
            "provider": result.provider,
            "query-type": result.query_type,
            "request-url": result.request_url,
        },
    )
    for index, record in enumerate(result.records):
        node = ET.SubElement(root, "record", {"index": str(index)})
        for key in sorted(record):
            field = ET.SubElement(node, "field", {"name": str(key)})
            field.text = _cell(record[key])
    return cast(bytes, ET.tostring(root, encoding="utf-8", xml_declaration=True)) + b"\n"


def _payload(result: QueryResult, output_format: MetadataFormat) -> bytes:
    if output_format == "json":
        return (
            json.dumps(result.to_dict(), ensure_ascii=False, sort_keys=True, indent=2) + "\n"
        ).encode("utf-8")
    if output_format == "jsonl":
        return (
            "".join(
                json.dumps(to_json_compatible(record), ensure_ascii=False, sort_keys=True) + "\n"
                for record in result.records
            )
        ).encode("utf-8")
    if output_format in {"csv", "tsv"}:
        return _tabular(result, delimiter="," if output_format == "csv" else "\t")
    return _xml(result)


def metadata(
    query: QueryResult,
    output_path: str | os.PathLike[str],
    *,
    format: MetadataFormat | None = None,
    config: DownloadConfig | None = None,
) -> DatasetDownloadResult:
    """Export a ``QueryResult`` as JSON/JSONL/CSV/TSV/XML plus an audit manifest."""

    if not isinstance(query, QueryResult):
        raise TypeError("query must be QueryResult.")
    resolved = resolved_config(config)
    target = Path(output_path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not target.is_file():
        raise IsADirectoryError(f"Metadata output must be a file: {target}")
    manifest = target.with_name(f"{target.name}.manifest.json")
    if manifest.exists() and not manifest.is_file():
        raise IsADirectoryError(f"Metadata manifest must be a file: {manifest}")
    if not resolved.overwrite and (target.exists() or manifest.exists()):
        raise FileExistsError("Refusing to overwrite existing metadata output or manifest.")
    output_format = _format(format, target)
    payload = _payload(query, output_format)
    if len(payload) > resolved.max_file_bytes:
        raise DownloadError("Metadata output exceeds max_file_bytes.", code="DOWNLOAD_SIZE_LIMIT")
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    sha256 = hashlib.sha256(payload).hexdigest()
    item = DownloadedFile(
        query.request_url,
        str(target),
        len(payload),
        md5,
        sha256,
        False,
    )
    with tempfile.TemporaryDirectory(prefix=".dnakit-metadata-", dir=str(target.parent)) as raw:
        stage = Path(raw)
        staged_target = stage / target.name
        staged_target.write_bytes(payload)
        staged_manifest = stage / manifest.name
        _write_manifest(
            {
                "kind": "query_metadata",
                "source": query.provider,
                "query_type": query.query_type,
                "request_url": query.request_url,
                "returned_records": len(query.records),
                "reported_total_count": query.total_count,
                "format": output_format,
                "file": {
                    "path": str(target),
                    "byte_size": len(payload),
                    "md5": md5,
                    "sha256": sha256,
                },
                "query_metadata": to_json_compatible(query.metadata),
                "provenance": query.provenance.to_dict(),
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
        "query_metadata",
        query.provider,
        str(target.parent),
        (item,),
        str(manifest),
        query.provenance,
        metadata={"format": output_format, "query_type": query.query_type},
    )


__all__ = ["MetadataFormat", "metadata"]
