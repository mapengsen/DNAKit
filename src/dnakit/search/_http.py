"""Bounded HTTPS helpers shared by public database adapters."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from typing import Any, BinaryIO, cast
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit
from urllib.request import Request, urlopen

from dnakit._version import __version__
from dnakit.exceptions import ConfigurationError, QueryError

from .models import SearchConfig

_SENSITIVE_QUERY_KEYS = frozenset(
    {
        "access_token",
        "api-key",
        "api_key",
        "apikey",
        "auth",
        "authorization",
        "email",
        "key",
        "sig",
        "signature",
        "token",
        "x-amz-credential",
        "x-amz-security-token",
        "x-amz-signature",
    }
)


def require_text(value: object, name: str, *, max_length: int = 4_096) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be non-empty text.", code="INVALID_QUERY")
    resolved = value.strip()
    if len(resolved) > max_length or "\x00" in resolved:
        raise ConfigurationError(
            f"{name} exceeds its text limit or contains NUL.", code="INVALID_QUERY"
        )
    return resolved


def require_https_base(base_url: str) -> str:
    resolved = require_text(base_url, "base_url", max_length=2_048).rstrip("/")
    parsed = urlsplit(resolved)
    if parsed.scheme != "https" or not parsed.netloc or parsed.username or parsed.password:
        raise ConfigurationError(
            "Remote database base_url must be an HTTPS origin without credentials.",
            code="INVALID_REMOTE_BASE_URL",
        )
    if parsed.query or parsed.fragment:
        raise ConfigurationError(
            "Remote database base_url cannot contain a query or fragment.",
            code="INVALID_REMOTE_BASE_URL",
        )
    return resolved


def build_url(
    base_url: str,
    path: str,
    params: Iterable[tuple[str, object]] = (),
) -> str:
    base = require_https_base(base_url)
    if not isinstance(path, str) or not path.startswith("/") or "\x00" in path:
        raise ConfigurationError("Remote API path must start with '/'.", code="INVALID_QUERY")
    encoded: list[tuple[str, str]] = []
    for key, raw_value in params:
        require_text(key, "query parameter name", max_length=128)
        values = raw_value if isinstance(raw_value, (tuple, list)) else (raw_value,)
        for value in values:
            if value is None:
                continue
            if isinstance(value, bool):
                encoded.append((key, "true" if value else "false"))
            elif isinstance(value, (str, int, float)) and not isinstance(value, complex):
                encoded.append((key, str(value)))
            else:
                raise ConfigurationError(
                    "Remote query parameters must be scalar values or scalar sequences.",
                    code="INVALID_QUERY",
                )
    query = urlencode(encoded, doseq=True)
    url = f"{base}{path}" + (f"?{query}" if query else "")
    if len(url) > 16_384:
        raise ConfigurationError("Remote query URL exceeds 16 KiB.", code="QUERY_URL_LIMIT")
    return url


def redact_url(url: str) -> str:
    parsed = urlsplit(url)
    redacted = [
        (key, "REDACTED" if key.casefold() in _SENSITIVE_QUERY_KEYS else value)
        for key, value in parse_qsl(parsed.query, keep_blank_values=True)
    ]
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path, urlencode(redacted), ""))


def _read_limited(stream: BinaryIO, *, max_bytes: int) -> bytes:
    chunks: list[bytes] = []
    observed = 0
    while True:
        chunk = stream.read(min(65_536, max_bytes - observed + 1))
        if not chunk:
            break
        observed += len(chunk)
        if observed > max_bytes:
            raise QueryError(
                "Remote database response exceeded max_response_bytes.",
                code="QUERY_RESPONSE_SIZE_LIMIT",
                context={"max_response_bytes": max_bytes},
            )
        chunks.append(chunk)
    return b"".join(chunks)


def request_json(
    url: str,
    config: SearchConfig,
    *,
    provider: str,
    headers: Mapping[str, str] | None = None,
) -> dict[str, Any] | list[Any]:
    request_headers = {
        "Accept": "application/json",
        "User-Agent": f"DNAKit/{__version__} public-database-adapter",
        **dict(headers or {}),
    }
    request = Request(url, headers=request_headers, method="GET")
    safe_url = redact_url(url)
    try:
        with urlopen(request, timeout=float(config.timeout)) as response:
            final_url = getattr(response, "geturl", lambda: url)()
            if urlsplit(final_url).scheme != "https":
                raise QueryError(
                    "Remote database redirected to a non-HTTPS URL.",
                    code="INSECURE_QUERY_REDIRECT",
                    context={"provider": provider},
                )
            payload = _read_limited(response, max_bytes=config.max_response_bytes)
    except HTTPError as exc:
        raise QueryError(
            "Remote database returned an HTTP error.",
            code="QUERY_HTTP_ERROR",
            context={"provider": provider, "status": exc.code, "url": safe_url},
        ) from exc
    except QueryError:
        raise
    except (OSError, TimeoutError, URLError) as exc:
        raise QueryError(
            "Could not reach the remote database.",
            code="QUERY_NETWORK_ERROR",
            context={"provider": provider, "url": safe_url},
        ) from exc
    try:
        decoded = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QueryError(
            "Remote database returned invalid JSON.",
            code="QUERY_RESPONSE_ERROR",
            context={"provider": provider, "url": safe_url},
        ) from exc
    if not isinstance(decoded, (dict, list)):
        raise QueryError(
            "Remote database returned an unexpected JSON root.",
            code="QUERY_RESPONSE_ERROR",
            context={"provider": provider, "url": safe_url},
        )
    return cast(dict[str, Any] | list[Any], decoded)


def request_text(
    url: str,
    config: SearchConfig,
    *,
    provider: str,
    headers: Mapping[str, str] | None = None,
    method: str = "GET",
    data: bytes | None = None,
) -> str:
    """Return one bounded UTF-8 text response from an HTTPS endpoint."""

    if method not in {"GET", "POST"}:
        raise ConfigurationError("Remote request method must be GET or POST.", code="INVALID_QUERY")
    if data is not None and (not isinstance(data, bytes) or len(data) > 5_000_000):
        raise ConfigurationError(
            "Remote request body must be bytes no larger than 5 MB.",
            code="QUERY_REQUEST_SIZE_LIMIT",
        )
    request_headers = {
        "Accept": "text/plain",
        "User-Agent": f"DNAKit/{__version__} public-database-adapter",
        **dict(headers or {}),
    }
    request = Request(url, headers=request_headers, data=data, method=method)
    safe_url = redact_url(url)
    try:
        with urlopen(request, timeout=float(config.timeout)) as response:
            final_url = getattr(response, "geturl", lambda: url)()
            if urlsplit(final_url).scheme != "https":
                raise QueryError(
                    "Remote database redirected to a non-HTTPS URL.",
                    code="INSECURE_QUERY_REDIRECT",
                    context={"provider": provider},
                )
            payload = _read_limited(response, max_bytes=config.max_response_bytes)
    except HTTPError as exc:
        raise QueryError(
            "Remote database returned an HTTP error.",
            code="QUERY_HTTP_ERROR",
            context={"provider": provider, "status": exc.code, "url": safe_url},
        ) from exc
    except QueryError:
        raise
    except (OSError, TimeoutError, URLError) as exc:
        raise QueryError(
            "Could not reach the remote database.",
            code="QUERY_NETWORK_ERROR",
            context={"provider": provider, "url": safe_url},
        ) from exc
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise QueryError(
            "Remote database returned invalid UTF-8 text.",
            code="QUERY_RESPONSE_ERROR",
            context={"provider": provider, "url": safe_url},
        ) from exc


def mapping_records(value: object, *, key: str | None = None) -> tuple[dict[str, object], ...]:
    selected = value
    if key is not None:
        if not isinstance(value, Mapping):
            return ()
        selected = value.get(key)
    if isinstance(selected, Mapping):
        return (dict(selected),)
    if not isinstance(selected, list):
        return ()
    return tuple(dict(item) for item in selected if isinstance(item, Mapping))


def limited_records(
    records: Iterable[Mapping[str, object]], config: SearchConfig
) -> tuple[dict[str, object], ...]:
    materialized = tuple(dict(record) for record in records)
    if len(materialized) > config.max_records:
        raise QueryError(
            "Remote database returned more records than max_records.",
            code="QUERY_RECORD_LIMIT",
            context={"max_records": config.max_records, "observed": len(materialized)},
        )
    return materialized


__all__ = [
    "build_url",
    "limited_records",
    "mapping_records",
    "redact_url",
    "request_json",
    "request_text",
    "require_https_base",
    "require_text",
]
