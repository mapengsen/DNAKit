"""Bounded transport for the P42-P45 public database adapters (no implicit paging)."""

from __future__ import annotations

import json
from collections.abc import Mapping
from datetime import datetime, timezone
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from dnakit.exceptions import ConfigurationError, QueryError
from dnakit.search._http import _read_limited, redact_url, require_https_base


def _request(
    url: str,
    *,
    body: Mapping[str, Any] | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    json_response: bool = True,
) -> dict[str, Any]:
    """Return the original response and provenance; raise on HTTP/service errors."""
    require_https_base(url.split("?", 1)[0])
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not 0 < timeout <= 300:
        raise ConfigurationError("timeout must be in (0, 300].")
    if (
        isinstance(max_bytes, bool)
        or not isinstance(max_bytes, int)
        or not 1 <= max_bytes <= 100_000_000
    ):
        raise ConfigurationError("max_bytes must be in [1, 100000000].")
    headers = {
        "Accept": "application/json" if json_response else "*/*",
        "User-Agent": "DNAKit public-data-adapter",
    }
    data = None
    if body is not None:
        try:
            data = json.dumps(dict(body), allow_nan=False).encode()
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("body must contain finite JSON data.") from exc
        if len(data) > 1_000_000:
            raise ConfigurationError("request body exceeds 1 MB.")
        headers["Content-Type"] = "application/json"
    request = Request(url, data=data, headers=headers, method="POST" if body is not None else "GET")
    try:
        with urlopen(request, timeout=timeout) as response:
            if not response.geturl().startswith("https://"):
                raise QueryError("Remote server redirected away from HTTPS.")
            raw = _read_limited(response, max_bytes=max_bytes)
            version = response.headers.get("Last-Modified") or response.headers.get("Date")
            content_type = response.headers.get("Content-Type", "")
    except (HTTPError, URLError, OSError, TimeoutError) as exc:
        raise QueryError(
            f"Remote request failed: {redact_url(url)} ({type(exc).__name__})."
        ) from exc
    try:
        text = raw.decode("utf-8")
        value = (json.loads(text) if text.strip() else None) if json_response else text
    except (UnicodeDecodeError, ValueError) as exc:
        raise QueryError("Remote response has an invalid encoding or JSON document.") from exc
    if isinstance(value, dict) and (value.get("error") or value.get("ERROR")):
        raise QueryError(
            f"Remote service rejected the request: {value.get('error') or value.get('ERROR')}."
        )
    if not json_response and (
        "<ERROR>" in text or text.lstrip().lower().startswith("<!doctype html")
    ):
        raise QueryError("Remote response contains an error or HTML instead of requested records.")
    return {
        "status": "completed",
        "data": value,
        "request_url": redact_url(url),
        "provenance": {
            "service_version": version,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "content_type": content_type,
            "paging": "one response; follow service page tokens explicitly",
        },
    }


def _params(parameters: Mapping[str, Any] | None, allowed: set[str]) -> dict[str, Any]:
    if parameters is None:
        return {}
    if not isinstance(parameters, Mapping) or any(not isinstance(k, str) for k in parameters):
        raise ConfigurationError("parameters must be a string-keyed mapping.")
    unknown = set(parameters) - allowed
    if unknown:
        raise ConfigurationError(
            f"Unsupported parameters: {sorted(unknown)}; allowed: {sorted(allowed)}."
        )
    return dict(parameters)
