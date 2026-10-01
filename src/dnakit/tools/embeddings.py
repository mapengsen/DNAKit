"""Optional cloud text embeddings for tool discovery, without neural dependencies."""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
import threading
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from http.client import HTTPException
from pathlib import Path
from typing import Literal

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.llm.client import _NoRedirects, _strict_json
from dnakit.llm.models import LLMConfig

ToolSearchMode = Literal["keyword", "semantic", "hybrid", "bm25", "bm25_hybrid", "rerank"]


def _default_cache_dir() -> Path:
    return (
        Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "dnakit/tool_embeddings"
    )


@dataclass(frozen=True, slots=True)
class EmbeddingConfig:
    """Configure an OpenAI-compatible /embeddings API independently of the chat LLM.

    Construction does not read credentials, create files, or contact the provider.
    Tool vectors are cached on disk; query vectors are not persisted. Set cache_dir
    to None for an in-memory tool cache only. Models and endpoints are explicit.
    """

    model: str
    base_url: str
    api_key: str | None = field(default=None, repr=False)
    api_key_env: str = "DNAKIT_EMBEDDING_API_KEY"
    dimensions: int | None = None
    batch_size: int = 10
    cache_dir: str | os.PathLike[str] | None = field(default_factory=_default_cache_dir)
    timeout_seconds: float = 60.0
    max_request_bytes: int = 1_000_000
    max_response_bytes: int = 16_000_000
    _connection: LLMConfig = field(init=False, repr=False)

    def __post_init__(self) -> None:
        # Reuse the existing URL, credential, timeout and byte-limit validation.
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
            raise ConfigurationError(exc.message, code="INVALID_EMBEDDING_CONFIG") from None
        for name, value, maximum in (
            ("batch_size", self.batch_size, 2048),
            ("dimensions", self.dimensions, 65_536),
        ):
            if name == "dimensions" and value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
                raise ConfigurationError(
                    f"{name} must be an integer in [1, {maximum}].",
                    code="INVALID_EMBEDDING_CONFIG",
                )
        if self.cache_dir is not None:
            try:
                if not os.fspath(self.cache_dir) or "\0" in os.fspath(self.cache_dir):
                    raise ValueError
                cache_dir = Path(self.cache_dir).expanduser()
            except (TypeError, ValueError, RuntimeError):
                raise ConfigurationError(
                    "cache_dir must be a filesystem path or None.",
                    code="INVALID_EMBEDDING_CONFIG",
                ) from None
            object.__setattr__(self, "cache_dir", cache_dir)
        object.__setattr__(self, "model", connection.model)
        object.__setattr__(self, "base_url", connection.base_url)
        object.__setattr__(self, "_connection", connection)

    @property
    def endpoint(self) -> str:
        """Append /embeddings once, retaining the configured API prefix."""
        return (
            self.base_url
            if self.base_url.endswith("/embeddings")
            else f"{self.base_url}/embeddings"
        )

    def resolved_api_key(self) -> str:
        """Resolve only the embedding credential, never the Agent's chat credential."""
        try:
            return self._connection.resolved_api_key()
        except ConfigurationError:
            raise ConfigurationError(
                "No valid embedding API key is configured.",
                code="MISSING_EMBEDDING_API_KEY",
                hint=f"Set {self.api_key_env} or pass EmbeddingConfig(api_key=...).",
            ) from None


def _unit_vector(value: object, dimensions: int | None = None) -> tuple[float, ...]:
    if not isinstance(value, list) or not 1 <= len(value) <= 65_536:
        raise ValueError("Expected a non-empty embedding vector.")
    if dimensions is not None and len(value) != dimensions:
        raise ValueError("Embedding dimensions do not match.")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in value):
        raise ValueError("Embedding components must be numbers.")
    numbers = tuple(float(x) for x in value)
    if not all(math.isfinite(x) for x in numbers):
        raise ValueError("Embedding components must be finite.")
    scale = max(abs(x) for x in numbers)
    if scale == 0:
        raise ValueError("Zero embedding vectors are unsupported.")
    scaled = tuple(x / scale for x in numbers)
    norm = math.sqrt(math.fsum(x * x for x in scaled))
    return tuple(x / norm for x in scaled)


def _embed(texts: list[str], config: EmbeddingConfig) -> list[tuple[float, ...]]:
    payload: dict[str, object] = {"model": config.model, "input": texts}
    if config.dimensions is not None:
        payload["dimensions"] = config.dimensions
    try:
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
    except UnicodeError:
        raise InputFormatError(
            "Embedding input contains invalid Unicode.", code="INVALID_EMBEDDING_INPUT"
        ) from None
    if len(body) > config.max_request_bytes:
        raise InputFormatError(
            "Embedding input exceeds max_request_bytes; this batch was not sent.",
            code="EMBEDDING_REQUEST_LIMIT",
            hint="Reduce batch_size or increase max_request_bytes explicitly.",
        )
    request = urllib.request.Request(
        config.endpoint,
        data=body,
        method="POST",
        headers={
            "Authorization": f"Bearer {config.resolved_api_key()}",
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "DNAKit-Tool-Embeddings/1",
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
            "The embedding API rejected the request; no automatic retry was made.",
            code="EMBEDDING_API_ERROR",
            context={"status": status},
            hint="Check the embedding URL, model, key, quota and supported batch size.",
        ) from None
    except (urllib.error.URLError, TimeoutError, OSError, HTTPException, ValueError):
        raise BackendExecutionError(
            "The embedding API could not be reached or its response could not be read.",
            code="EMBEDDING_API_UNAVAILABLE",
            hint="Check the network and timeout; no automatic retry was made.",
        ) from None
    if len(raw) > config.max_response_bytes:
        raise BackendExecutionError(
            "The embedding response exceeds max_response_bytes.", code="EMBEDDING_RESPONSE_LIMIT"
        )
    try:
        envelope = _strict_json(raw)
        if not isinstance(envelope, dict) or envelope.get("error"):
            raise ValueError
        data = envelope.get("data")
        if not isinstance(data, list) or len(data) != len(texts):
            raise ValueError
        by_index: dict[int, tuple[float, ...]] = {}
        dimensions = config.dimensions
        for item in data:
            if not isinstance(item, dict):
                raise ValueError
            index = item.get("index")
            if (
                isinstance(index, bool)
                or not isinstance(index, int)
                or not 0 <= index < len(texts)
                or index in by_index
            ):
                raise ValueError
            vector = _unit_vector(item.get("embedding"), dimensions)
            dimensions = len(vector)
            by_index[index] = vector
        return [by_index[i] for i in range(len(texts))]
    except (ValueError, TypeError, OverflowError, UnicodeError, RecursionError):
        raise BackendExecutionError(
            "The embedding API returned invalid vectors or response indexes.",
            code="INVALID_EMBEDDING_RESPONSE",
        ) from None


class _ToolEmbeddingIndex:
    """Lazy, content-addressed tool vectors with atomic JSON cache writes."""

    def __init__(self, config: EmbeddingConfig) -> None:
        self.config = config
        self._vectors: dict[str, tuple[float, ...]] = {}
        self._lock = threading.Lock()

    def _key(self, text: str) -> str:
        payload = [1, self.config.endpoint, self.config.model, self.config.dimensions, text]
        return hashlib.sha256(json.dumps(payload, ensure_ascii=True).encode()).hexdigest()

    def _load(self, key: str, dimensions: int) -> tuple[float, ...] | None:
        vector = self._vectors.get(key)
        if vector is not None and len(vector) == dimensions:
            return vector
        if self.config.cache_dir is None:
            return None
        try:
            with (Path(self.config.cache_dir) / f"{key}.json").open("rb") as stream:
                raw = stream.read(self.config.max_response_bytes + 1)
            if len(raw) > self.config.max_response_bytes:
                return None
            entry = _strict_json(raw)
            if not isinstance(entry, dict) or entry.get("key") != key:
                return None
            vector = _unit_vector(entry.get("embedding"), dimensions)
        except (OSError, ValueError, TypeError, OverflowError, RecursionError):
            return None
        self._vectors[key] = vector
        return vector

    def _save(self, key: str, vector: tuple[float, ...]) -> None:
        self._vectors[key] = vector
        if self.config.cache_dir is None:
            return
        temporary: str | None = None
        try:
            directory = Path(self.config.cache_dir)
            directory.mkdir(parents=True, exist_ok=True, mode=0o700)
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=directory, delete=False
            ) as stream:
                temporary = stream.name
                json.dump({"key": key, "embedding": vector}, stream, allow_nan=False)
            os.replace(temporary, directory / f"{key}.json")
        except OSError:
            raise ConfigurationError(
                "The tool embedding cache could not be written.",
                code="EMBEDDING_CACHE_ERROR",
                hint="Use a writable cache_dir or set cache_dir=None for an in-memory cache.",
            ) from None
        finally:
            if temporary is not None:
                Path(temporary).unlink(missing_ok=True)

    def rank(self, query: str, documents: list[tuple[str, str]]) -> list[str]:
        with self._lock:
            # Validate credentials and query before paying to embed the tool catalog.
            query_vector = _embed([query], self.config)[0]
            dimensions = len(query_vector)
            keys = [self._key(text) for _, text in documents]
            vectors = {key: self._load(key, dimensions) for key in keys}
            missing = [
                (key, text)
                for key, (_, text) in zip(keys, documents, strict=True)
                if vectors[key] is None
            ]
            for start in range(0, len(missing), self.config.batch_size):
                batch = missing[start : start + self.config.batch_size]
                embedded = _embed([text for _, text in batch], self.config)
                if any(len(vector) != dimensions for vector in embedded):
                    raise BackendExecutionError(
                        "Query and tool embedding dimensions differ.",
                        code="INVALID_EMBEDDING_RESPONSE",
                    )
                for (key, _), vector in zip(batch, embedded, strict=True):
                    self._save(key, vector)
                    vectors[key] = vector
            scores = []
            for (name, _), key in zip(documents, keys, strict=True):
                cached_vector = vectors[key]
                assert cached_vector is not None
                score = math.fsum(x * y for x, y in zip(query_vector, cached_vector, strict=True))
                scores.append((-score, name))
            scores.sort()
            return [name for _, name in scores]


__all__ = ["EmbeddingConfig", "ToolSearchMode"]
