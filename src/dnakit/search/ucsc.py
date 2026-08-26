"""Bounded UCSC Genome Browser REST API adapters."""

from __future__ import annotations

import fnmatch
import re
from collections.abc import Mapping
from typing import cast

from dnakit.core import Provenance
from dnakit.exceptions import ConfigurationError, QueryError

from ._http import build_url, limited_records, redact_url, request_json, require_text
from ._shared import adapter_provenance, resolved_config
from .models import QueryResult, SearchConfig

DEFAULT_UCSC_API_URL = "https://api.genome.ucsc.edu"
_REGION = re.compile(r"^(?P<chrom>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)$")


def _name(value: str, field: str) -> str:
    return require_text(value, field, max_length=256)


def _region(value: str, *, max_span: int) -> tuple[str, int, int]:
    text = require_text(value, "region", max_length=512)
    match = _REGION.fullmatch(text)
    if match is None:
        raise ConfigurationError(
            "region must use 'chromosome:start-end' with 0-based half-open coordinates.",
            code="INVALID_REGION_QUERY",
        )
    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end <= start or end - start > max_span:
        raise ConfigurationError(
            "region is empty or exceeds the endpoint span limit.",
            code="REGION_QUERY_LIMIT",
            context={"span": end - start, "max_span": max_span},
        )
    return chrom, start, end


def _provenance(filters: Mapping[str, object], *, version: str | None = None) -> Provenance:
    return adapter_provenance(
        "UCSC Genome Browser REST API",
        citation_url="https://genome.ucsc.edu/goldenPath/help/api.html",
        filters=filters,
        version=version,
    )


def genomes(
    *,
    contains: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """List UCSC browser assemblies, optionally filtering their text fields."""

    resolved = resolved_config(config)
    needle = None if contains is None else require_text(contains, "contains").casefold()
    url = build_url(api_base_url, "/list/ucscGenomes")
    payload = request_json(url, resolved, provider="UCSC")
    root = payload.get("ucscGenomes") if isinstance(payload, Mapping) else None
    if not isinstance(root, Mapping):
        raise QueryError("UCSC returned an unexpected genome list.", code="QUERY_RESPONSE_ERROR")
    records: list[dict[str, object]] = []
    for identifier, raw in root.items():
        if not isinstance(identifier, str) or not isinstance(raw, Mapping):
            continue
        record: dict[str, object] = {"genome": identifier, **dict(raw)}
        if needle is None or needle in " ".join(str(value) for value in record.values()).casefold():
            records.append(record)
    records.sort(key=lambda item: str(item["genome"]))
    bounded = limited_records(records, resolved)
    return QueryResult(
        "assembly",
        "UCSC",
        redact_url(url),
        bounded,
        _provenance({"contains": contains}),
        total_count=len(bounded),
        metadata={"filter_applied_client_side": needle is not None},
    )


def chromosomes(
    genome: str,
    *,
    track: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """List chromosome, contig, scaffold, and organelle names and lengths."""

    resolved = resolved_config(config)
    assembly = _name(genome, "genome")
    track_name = None if track is None else _name(track, "track")
    url = build_url(
        api_base_url,
        "/list/chromosomes",
        (("genome", assembly), ("track", track_name)),
    )
    payload = request_json(url, resolved, provider="UCSC")
    raw = payload.get("chromosomes") if isinstance(payload, Mapping) else None
    if not isinstance(raw, Mapping):
        raise QueryError(
            "UCSC returned an unexpected chromosome list.", code="QUERY_RESPONSE_ERROR"
        )
    records = [
        {"name": name, "length": length}
        for name, length in raw.items()
        if isinstance(name, str)
        and isinstance(length, int)
        and not isinstance(length, bool)
        and length >= 0
    ]
    records.sort(key=lambda item: str(item["name"]))
    bounded = limited_records(records, resolved)
    data_time = payload.get("dataTime") if isinstance(payload, Mapping) else None
    return QueryResult(
        "chromosome",
        "UCSC",
        redact_url(url),
        bounded,
        _provenance(
            {"genome": assembly, "track": track_name},
            version=data_time if isinstance(data_time, str) else None,
        ),
        total_count=len(bounded),
        metadata={"genome": assembly, "track": track_name},
    )


def sequence(
    genome: str,
    region: str,
    *,
    reverse_complement: bool = False,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """Return a bounded UCSC assembly sequence for a 0-based half-open region."""

    resolved = resolved_config(config)
    assembly = _name(genome, "genome")
    chrom, start, end = _region(region, max_span=10_000_000)
    if not isinstance(reverse_complement, bool):
        raise ConfigurationError("reverse_complement must be boolean.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        "/getData/sequence",
        (
            ("genome", assembly),
            ("chrom", chrom),
            ("start", start),
            ("end", end),
            ("revComp", reverse_complement),
        ),
    )
    payload = request_json(url, resolved, provider="UCSC")
    dna = payload.get("dna") if isinstance(payload, Mapping) else None
    if not isinstance(dna, str):
        raise QueryError("UCSC sequence response has no DNA text.", code="QUERY_RESPONSE_ERROR")
    data_time = payload.get("dataTime") if isinstance(payload, Mapping) else None
    return QueryResult(
        "sequence",
        "UCSC",
        redact_url(url),
        (
            {
                "genome": assembly,
                "chromosome": chrom,
                "start_0based": start,
                "end_0based": end,
                "reverse_complement": reverse_complement,
                "seq": dna,
            },
        ),
        _provenance(
            {
                "genome": assembly,
                "region": region,
                "reverse_complement": reverse_complement,
            },
            version=data_time if isinstance(data_time, str) else None,
        ),
        total_count=1,
        metadata={"input_coordinate_system": "0-based half-open"},
    )


def tracks(
    genome: str,
    *,
    contains: str | None = None,
    leaves_only: bool = True,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """List UCSC track definitions, with an optional client-side text filter."""

    resolved = resolved_config(config)
    assembly = _name(genome, "genome")
    needle = None if contains is None else require_text(contains, "contains").casefold()
    if not isinstance(leaves_only, bool):
        raise ConfigurationError("leaves_only must be boolean.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        "/list/tracks",
        (("genome", assembly), ("trackLeavesOnly", leaves_only)),
    )
    payload = request_json(url, resolved, provider="UCSC")
    raw = payload.get(assembly) if isinstance(payload, Mapping) else None
    if not isinstance(raw, Mapping):
        raise QueryError("UCSC returned an unexpected track list.", code="QUERY_RESPONSE_ERROR")
    records: list[dict[str, object]] = []
    for track_name, definition in raw.items():
        if not isinstance(track_name, str) or not isinstance(definition, Mapping):
            continue
        record: dict[str, object] = {"track": track_name, **dict(definition)}
        if needle is None or needle in " ".join(str(value) for value in record.values()).casefold():
            records.append(record)
    records.sort(key=lambda item: str(item["track"]))
    bounded = limited_records(records, resolved)
    data_time = payload.get("dataTime") if isinstance(payload, Mapping) else None
    return QueryResult(
        "tracks",
        "UCSC",
        redact_url(url),
        bounded,
        _provenance(
            {"genome": assembly, "contains": contains, "leaves_only": leaves_only},
            version=data_time if isinstance(data_time, str) else None,
        ),
        total_count=len(bounded),
        metadata={"filter_applied_client_side": needle is not None},
    )


def track_data(
    genome: str,
    track: str,
    region: str,
    *,
    max_items: int = 1_000,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """Query bounded annotations, repeats, regulation, or signal from one UCSC track."""

    resolved = resolved_config(config)
    assembly = _name(genome, "genome")
    track_name = _name(track, "track")
    chrom, start, end = _region(region, max_span=5_000_000)
    if (
        isinstance(max_items, bool)
        or not isinstance(max_items, int)
        or not 1 <= max_items <= min(1_000_000, resolved.max_records)
    ):
        raise ConfigurationError(
            "max_items exceeds SearchConfig.max_records.", code="QUERY_RECORD_LIMIT"
        )
    url = build_url(
        api_base_url,
        "/getData/track",
        (
            ("genome", assembly),
            ("track", track_name),
            ("chrom", chrom),
            ("start", start),
            ("end", end),
            ("maxItemsOutput", max_items),
        ),
    )
    payload = request_json(url, resolved, provider="UCSC")
    raw = payload.get(track_name) if isinstance(payload, Mapping) else None
    if not isinstance(raw, list):
        raise QueryError("UCSC returned unexpected track data.", code="QUERY_RESPONSE_ERROR")
    records: list[dict[str, object]] = []
    for item in raw:
        if not isinstance(item, Mapping):
            continue
        record = dict(item)
        item_start = item.get("chromStart")
        item_end = item.get("chromEnd")
        if isinstance(item_start, int) and isinstance(item_end, int):
            record["start_0based"] = item_start
            record["end_0based"] = item_end
        records.append(record)
    bounded = limited_records(records, resolved)
    data_time = payload.get("dataTime") if isinstance(payload, Mapping) else None
    server_count = payload.get("itemsReturned") if isinstance(payload, Mapping) else None
    total = (
        server_count
        if isinstance(server_count, int) and not isinstance(server_count, bool)
        else len(bounded)
    )
    return QueryResult(
        "annotation",
        "UCSC",
        redact_url(url),
        bounded,
        _provenance(
            {"genome": assembly, "track": track_name, "region": region},
            version=data_time if isinstance(data_time, str) else None,
        ),
        total_count=max(total, len(bounded)),
        metadata={
            "input_coordinate_system": "0-based half-open",
            "track": track_name,
            "track_type": payload.get("trackType") if isinstance(payload, Mapping) else None,
        },
    )


def files(
    genome: str,
    *,
    pattern: str | None = None,
    limit: int = 1_000,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_UCSC_API_URL,
) -> QueryResult:
    """List UCSC downloadable files, optionally filtering paths with a glob."""

    resolved = resolved_config(config)
    assembly = _name(genome, "genome")
    if (
        isinstance(limit, bool)
        or not isinstance(limit, int)
        or not 1 <= limit <= min(1_000_000, resolved.max_records)
    ):
        raise ConfigurationError(
            "limit exceeds SearchConfig.max_records.", code="QUERY_RECORD_LIMIT"
        )
    glob = None if pattern is None else require_text(pattern, "pattern", max_length=1_024)
    url = build_url(
        api_base_url,
        "/list/files",
        (("genome", assembly), ("maxItemsOutput", limit)),
    )
    payload = request_json(url, resolved, provider="UCSC")
    raw = payload.get("urlList") if isinstance(payload, Mapping) else None
    if not isinstance(raw, list):
        raise QueryError("UCSC returned an unexpected file list.", code="QUERY_RESPONSE_ERROR")
    records = [
        dict(item)
        for item in raw
        if isinstance(item, Mapping)
        and isinstance(item.get("url"), str)
        and (glob is None or fnmatch.fnmatch(cast(str, item["url"]), glob))
    ]
    bounded = limited_records(records, resolved)
    data_time = payload.get("dataTime") if isinstance(payload, Mapping) else None
    return QueryResult(
        "files",
        "UCSC",
        redact_url(url),
        bounded,
        _provenance(
            {"genome": assembly, "pattern": pattern, "server_limit": limit},
            version=data_time if isinstance(data_time, str) else None,
        ),
        total_count=len(bounded),
        metadata={
            "filter_applied_client_side": glob is not None,
            "server_items_returned": payload.get("itemsReturned")
            if isinstance(payload, Mapping)
            else None,
            "server_limit_reached": payload.get("maxItemsLimit")
            if isinstance(payload, Mapping)
            else None,
        },
    )


__all__ = [
    "DEFAULT_UCSC_API_URL",
    "chromosomes",
    "files",
    "genomes",
    "sequence",
    "track_data",
    "tracks",
]
