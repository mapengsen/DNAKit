"""ENCODE Portal REST search adapter."""

from __future__ import annotations

import re
from collections.abc import Mapping

from dnakit.exceptions import ConfigurationError, QueryError

from ._http import (
    build_url,
    limited_records,
    mapping_records,
    redact_url,
    request_json,
    require_text,
)
from ._shared import adapter_provenance, resolved_config
from .models import QueryResult, SearchConfig

DEFAULT_ENCODE_URL = "https://www.encodeproject.org"
_FILTER_KEY = re.compile(r"^[A-Za-z0-9_.@!-]+$")


def encode_search(
    *,
    object_type: str = "Experiment",
    search_term: str | None = None,
    filters: Mapping[str, object] | None = None,
    limit: int = 25,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENCODE_URL,
    query_type: str = "regulation",
) -> QueryResult:
    """Search released ENCODE experiments, files, biosamples, or annotations."""

    resolved = resolved_config(config)
    type_value = require_text(object_type, "object_type", max_length=128)
    if (
        isinstance(limit, bool)
        or not isinstance(limit, int)
        or not 1 <= limit <= resolved.max_records
    ):
        raise ConfigurationError(
            "limit exceeds SearchConfig max_records.", code="QUERY_RECORD_LIMIT"
        )
    params: list[tuple[str, object]] = [
        ("type", type_value),
        ("format", "json"),
        ("frame", "object"),
        ("limit", limit),
    ]
    if search_term is not None:
        params.append(("searchTerm", require_text(search_term, "search_term")))
    safe_filters: dict[str, object] = {}
    for key, value in dict(filters or {}).items():
        if not isinstance(key, str) or _FILTER_KEY.fullmatch(key) is None:
            raise ConfigurationError("ENCODE filter name is invalid.", code="INVALID_QUERY")
        if isinstance(value, (str, int, float, bool)) and not isinstance(value, complex):
            safe_filters[key] = value
        elif isinstance(value, (tuple, list)) and all(
            isinstance(item, (str, int, float, bool)) and not isinstance(item, complex)
            for item in value
        ):
            safe_filters[key] = tuple(value)
        else:
            raise ConfigurationError(
                "ENCODE filter values must be scalars or scalar sequences.", code="INVALID_QUERY"
            )
        params.append((key, safe_filters[key]))
    url = build_url(api_base_url, "/search/", params)
    payload = request_json(url, resolved, provider="ENCODE")
    if not isinstance(payload, Mapping):
        raise QueryError("ENCODE returned an unexpected response.", code="QUERY_RESPONSE_ERROR")
    records = limited_records(mapping_records(payload, key="@graph"), resolved)
    raw_total = payload.get("total")
    total = (
        raw_total
        if isinstance(raw_total, int) and not isinstance(raw_total, bool)
        else len(records)
    )
    return QueryResult(
        query_type,
        "ENCODE",
        redact_url(url),
        records,
        adapter_provenance(
            "ENCODE Portal",
            citation_url="https://www.encodeproject.org/help/rest-api/",
            filters={
                "object_type": type_value,
                "search_term": search_term,
                "filters": safe_filters,
                "limit": limit,
            },
        ),
        total_count=max(total, len(records)),
        metadata={"object_type": type_value, "returned_count": len(records)},
    )


__all__ = ["DEFAULT_ENCODE_URL", "encode_search"]
