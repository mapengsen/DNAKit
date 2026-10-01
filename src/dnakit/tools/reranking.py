"""Explicit, bounded HTTP reranking; no chat model or neural dependencies."""

from __future__ import annotations

import json
import math
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from http.client import HTTPException

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.llm.client import _NoRedirects, _strict_json
from dnakit.llm.models import LLMConfig


@dataclass(frozen=True, slots=True)
class RerankerConfig:
    """Configure a Cohere-compatible POST /rerank endpoint and fixed recall budget.

    The caller supplies the model and API prefix (e.g. a provider's /v2 prefix).
    Requests contain model, query, documents and top_n; responses must contain
    results with index and relevance_score. No credentials are read at construction.
    """

    model: str
    base_url: str
    api_key: str | None = field(default=None, repr=False)
    api_key_env: str = "DNAKIT_RERANK_API_KEY"
    candidate_limit: int = 50
    timeout_seconds: float = 60.0
    max_request_bytes: int = 1_000_000
    max_response_bytes: int = 2_000_000
    _connection: LLMConfig = field(init=False, repr=False)

    def __post_init__(self) -> None:
        try:
            connection = LLMConfig(
                provider="custom",
                model=self.model,
                base_url=self.base_url,
                api_key=self.api_key,
                api_key_env=self.api_key_env,
                timeout_seconds=self.timeout_seconds,
                max_request_bytes=self.max_request_bytes,
                max_response_bytes=self.max_response_bytes,
            )
        except ConfigurationError as exc:
            raise ConfigurationError(exc.message, code="INVALID_RERANKER_CONFIG") from None
        if (
            isinstance(self.candidate_limit, bool)
            or not isinstance(self.candidate_limit, int)
            or not 1 <= self.candidate_limit <= 1000
        ):
            raise ConfigurationError(
                "candidate_limit must be an integer in [1, 1000].",
                code="INVALID_RERANKER_CONFIG",
            )
        object.__setattr__(self, "model", connection.model)
        object.__setattr__(self, "base_url", connection.base_url)
        object.__setattr__(self, "_connection", connection)

    @property
    def endpoint(self) -> str:
        return self.base_url if self.base_url.endswith("/rerank") else f"{self.base_url}/rerank"

    def resolved_api_key(self) -> str:
        """Resolve only the reranking credential, never chat or embedding keys."""
        try:
            return self._connection.resolved_api_key()
        except ConfigurationError:
            raise ConfigurationError(
                "No valid reranker API key is configured.",
                code="MISSING_RERANKER_API_KEY",
                hint=f"Set {self.api_key_env} or pass RerankerConfig(api_key=...).",
            ) from None


def _rerank(query: str, documents: list[tuple[str, str]], config: RerankerConfig) -> list[str]:
    if not documents:
        return []
    # Rank the entire fixed pool so changing output limit does not change the experiment.
    payload = {
        "model": config.model,
        "query": query,
        "documents": [text for _, text in documents],
        "top_n": len(documents),
    }
    try:
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
    except UnicodeError:
        raise InputFormatError(
            "Reranker input contains invalid Unicode.", code="INVALID_RERANK_INPUT"
        ) from None
    if len(body) > config.max_request_bytes:
        raise InputFormatError(
            "Reranker input exceeds max_request_bytes; the request was not sent.",
            code="RERANK_REQUEST_LIMIT",
            hint="Reduce candidate_limit or increase max_request_bytes explicitly.",
        )
    request = urllib.request.Request(
        config.endpoint,
        data=body,
        method="POST",
        headers={
            "Authorization": f"Bearer {config.resolved_api_key()}",
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "DNAKit-Tool-Reranking/1",
        },
    )
    try:
        opener = urllib.request.build_opener(_NoRedirects())
        with opener.open(request, timeout=float(config.timeout_seconds)) as response:
            raw = response.read(config.max_response_bytes + 1)
    except urllib.error.HTTPError as exc:
        status = int(exc.code)
        exc.close()
        raise BackendExecutionError(
            "The reranker API rejected the request; no automatic retry was made.",
            code="RERANK_API_ERROR",
            context={"status": status},
        ) from None
    except (urllib.error.URLError, TimeoutError, OSError, HTTPException, ValueError):
        raise BackendExecutionError(
            "The reranker API could not be reached or its response could not be read.",
            code="RERANK_API_UNAVAILABLE",
        ) from None
    if len(raw) > config.max_response_bytes:
        raise BackendExecutionError(
            "The reranker response exceeds max_response_bytes.", code="RERANK_RESPONSE_LIMIT"
        )
    try:
        envelope = _strict_json(raw)
        if not isinstance(envelope, dict) or envelope.get("error"):
            raise ValueError
        results = envelope.get("results")
        if not isinstance(results, list) or len(results) != len(documents):
            raise ValueError
        scores: dict[int, float] = {}
        for item in results:
            if not isinstance(item, dict):
                raise ValueError
            index, score = item.get("index"), item.get("relevance_score")
            if (
                isinstance(index, bool)
                or not isinstance(index, int)
                or not 0 <= index < len(documents)
                or index in scores
                or isinstance(score, bool)
                or not isinstance(score, (int, float))
                or not math.isfinite(score)
            ):
                raise ValueError
            scores[index] = float(score)
        # Equal scores retain the first-stage ordering, independent of response row order.
        return [documents[i][0] for i in sorted(scores, key=lambda i: (-scores[i], i))]
    except (ValueError, TypeError, OverflowError, UnicodeError, RecursionError):
        raise BackendExecutionError(
            "The reranker API returned invalid scores or candidate indexes.",
            code="INVALID_RERANK_RESPONSE",
        ) from None
