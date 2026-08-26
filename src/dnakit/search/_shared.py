"""Shared validation and provenance helpers for query adapters."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from urllib.parse import quote

from dnakit.core import (
    Citation,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    ReferenceInfo,
)
from dnakit.exceptions import ConfigurationError

from ._http import require_text
from .models import SearchConfig


def resolved_config(config: SearchConfig | None) -> SearchConfig:
    if config is None:
        return SearchConfig()
    if not isinstance(config, SearchConfig):
        raise TypeError("config must be SearchConfig or None.")
    return config


def query_values(
    value: str | int | Sequence[str | int],
    *,
    name: str = "query",
    maximum: int = 100,
) -> tuple[str, ...]:
    raw: Iterable[str | int]
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        raw = (value,)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        raw = value
    else:
        raise ConfigurationError(
            f"{name} must be text, an integer, or a sequence of them.", code="INVALID_QUERY"
        )
    resolved: list[str] = []
    for item in raw:
        if isinstance(item, bool) or not isinstance(item, (str, int)):
            raise ConfigurationError(
                f"{name} values must be text or integers.", code="INVALID_QUERY"
            )
        resolved.append(require_text(str(item), name, max_length=1_024))
        if len(resolved) > maximum:
            raise ConfigurationError(
                f"{name} accepts at most {maximum} values.", code="QUERY_VALUE_LIMIT"
            )
    if not resolved:
        raise ConfigurationError(f"{name} must not be empty.", code="INVALID_QUERY")
    return tuple(resolved)


def quoted_path_values(values: Sequence[str]) -> str:
    return quote(",".join(values), safe=",")


def page_size(value: int, config: SearchConfig, *, maximum: int = 1_000) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError("page_size must be a positive integer.", code="INVALID_QUERY")
    if value > maximum or value > config.max_records:
        raise ConfigurationError(
            "page_size exceeds the provider or SearchConfig record limit.",
            code="QUERY_RECORD_LIMIT",
            context={
                "page_size": value,
                "provider_maximum": maximum,
                "max_records": config.max_records,
            },
        )
    return value


def adapter_provenance(
    provider: str,
    *,
    citation_url: str,
    filters: Mapping[str, object] | None = None,
    version: str | None = None,
) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.EXTERNAL,
            origin_class=OriginClass.INTEGRATION,
            citations=(
                Citation(provider.casefold().replace(" ", "-"), title=provider, url=citation_url),
            ),
        ),
        reference=ReferenceInfo(provider, version=version, filters=filters),
    )


__all__ = [
    "adapter_provenance",
    "page_size",
    "query_values",
    "quoted_path_values",
    "resolved_config",
]
