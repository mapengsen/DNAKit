"""Immutable configuration and results for public database queries."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import Provenance
from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class SearchConfig:
    """Resource and contact settings shared by remote query adapters.

    ``api_key`` is sent only to NCBI services.  It is never included in a
    ``QueryResult.request_url`` or serialized result.
    """

    timeout: float = 30.0
    max_response_bytes: int = 20_000_000
    max_records: int = 1_000
    api_key: str | None = None
    email: str | None = None
    tool: str = "dnakit"

    def __post_init__(self) -> None:
        if (
            isinstance(self.timeout, bool)
            or not isinstance(self.timeout, (int, float))
            or not math.isfinite(self.timeout)
            or not 0 < self.timeout <= 300
        ):
            raise ConfigurationError(
                "SearchConfig timeout must be in (0, 300].", code="INVALID_SEARCH_CONFIG"
            )
        for name, value, maximum in (
            ("max_response_bytes", self.max_response_bytes, 100_000_000),
            ("max_records", self.max_records, 10_000),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(
                    f"SearchConfig {name} must be in [1, {maximum}].",
                    code="INVALID_SEARCH_CONFIG",
                )
        for name in ("api_key", "email"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(
                    f"SearchConfig {name} must be None or non-empty text.",
                    code="INVALID_SEARCH_CONFIG",
                )
        if not isinstance(self.tool, str) or not self.tool.strip():
            raise ConfigurationError(
                "SearchConfig tool must be non-empty text.", code="INVALID_SEARCH_CONFIG"
            )


@dataclass(frozen=True, slots=True)
class QueryProgress:
    """One optional progress event emitted by a batch query."""

    completed: int
    total: int
    item: str

    def __post_init__(self) -> None:
        if (
            isinstance(self.completed, bool)
            or not isinstance(self.completed, int)
            or isinstance(self.total, bool)
            or not isinstance(self.total, int)
            or not 0 <= self.completed <= self.total
            or self.total < 1
        ):
            raise ConfigurationError(
                "QueryProgress counts are invalid.", code="INVALID_QUERY_PROGRESS"
            )
        if not isinstance(self.item, str) or not self.item:
            raise ConfigurationError(
                "QueryProgress item must be non-empty text.", code="INVALID_QUERY_PROGRESS"
            )


@dataclass(frozen=True, init=False)
class QueryResult:
    """Auditable, provider-neutral records returned by one logical query."""

    query_type: str
    provider: str
    request_url: str
    records: tuple[FrozenDict, ...]
    total_count: int | None
    next_page_token: str | None
    metadata: FrozenDict
    provenance: Provenance

    def __init__(
        self,
        query_type: str,
        provider: str,
        request_url: str,
        records: Iterable[Mapping[str, object]],
        provenance: Provenance,
        *,
        total_count: int | None = None,
        next_page_token: str | None = None,
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        for name, value in (
            ("query_type", query_type),
            ("provider", provider),
            ("request_url", request_url),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(
                    f"QueryResult {name} must be non-empty text.",
                    code="INVALID_QUERY_RESULT",
                )
        frozen_records = tuple(freeze_mapping(record) for record in records)
        if total_count is not None and (
            isinstance(total_count, bool)
            or not isinstance(total_count, int)
            or total_count < len(frozen_records)
        ):
            raise ConfigurationError(
                "QueryResult total_count must be None or at least the returned record count.",
                code="INVALID_QUERY_RESULT",
            )
        if next_page_token is not None and (
            not isinstance(next_page_token, str) or not next_page_token.strip()
        ):
            raise ConfigurationError(
                "QueryResult next_page_token must be None or non-empty text.",
                code="INVALID_QUERY_RESULT",
            )
        if not isinstance(provenance, Provenance):
            raise ConfigurationError(
                "QueryResult provenance must be Provenance.", code="INVALID_QUERY_RESULT"
            )
        object.__setattr__(self, "query_type", query_type)
        object.__setattr__(self, "provider", provider)
        object.__setattr__(self, "request_url", request_url)
        object.__setattr__(self, "records", frozen_records)
        object.__setattr__(self, "total_count", total_count)
        object.__setattr__(self, "next_page_token", next_page_token)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))
        object.__setattr__(self, "provenance", provenance)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible result without credentials."""

        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["QueryProgress", "QueryResult", "SearchConfig"]
