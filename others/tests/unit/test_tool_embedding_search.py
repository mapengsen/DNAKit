"""Offline cloud-protocol, retrieval, cache and entry-point regression tests."""

from __future__ import annotations

import io
import json
import threading
import urllib.error
from dataclasses import replace
from email.message import Message
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.tools import AgentToolRegistry, EmbeddingConfig, default_tool_registry
from dnakit.tools.embeddings import _embed

RC = "dnakit.ops.reverse_complement"
DOWNLOAD = "dnakit.references.download_genome"
WEIGHT = "dnakit.thermodynamics.molecular_weight"
QUERY = "Find the strand paired in the opposite direction"


@pytest.fixture
def config(tmp_path: Path) -> EmbeddingConfig:
    return EmbeddingConfig(
        model="offline-embedding",
        base_url="https://embedding.example/v1",
        api_key="offline-test-key",
        cache_dir=tmp_path / "vectors",
        batch_size=2,
    )


def _registry(config: EmbeddingConfig | None) -> AgentToolRegistry:
    catalog = default_tool_registry(verbose=False)
    return AgentToolRegistry(
        tuple(catalog.resolve(name) for name in (RC, DOWNLOAD, WEIGHT)),
        embedding_config=config,
    )


@pytest.fixture
def endpoint(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []

    def open_request(request: Any, *, timeout: float) -> io.BytesIO:
        body = json.loads(request.data)
        assert request.full_url == "https://embedding.example/v1/embeddings"
        assert request.get_header("Authorization") == "Bearer offline-test-key"
        assert timeout == 60.0
        calls.append(body)
        vectors = []
        for text in body["input"]:
            if text.startswith(DOWNLOAD):
                vector = [0, 1, 0]
            elif text.startswith(WEIGHT):
                vector = [0, 0, 1]
            else:
                vector = [1, 0, 0]
            vectors.append(vector)
        # Providers may return rows out of order; indexes control correspondence.
        data = [{"index": i, "embedding": v} for i, v in enumerate(vectors)]
        return io.BytesIO(json.dumps({"data": list(reversed(data))}).encode())

    monkeypatch.setattr(
        "dnakit.tools.embeddings.urllib.request.build_opener",
        lambda *handlers: type("Endpoint", (), {"open": staticmethod(open_request)})(),
    )
    return calls


def test_keyword_default_never_reads_keys_calls_api_or_creates_cache(
    config: EmbeddingConfig,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fail = Mock(side_effect=AssertionError("Keyword search must remain local"))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", fail)
    config = replace(config, api_key=None)
    registry = _registry(config)
    assert registry.search("reverse_complement")[0].name == RC
    assert registry.search("reverse complement") == _registry(None).search("reverse complement")
    assert config.cache_dir is not None
    assert not Path(config.cache_dir).exists()
    fail.assert_not_called()


def test_semantic_finds_synonyms_without_keyword_matches_and_never_executes(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    registry = _registry(config)
    execute = Mock(side_effect=AssertionError("Discovery must not execute candidates"))
    monkeypatch.setattr(registry, "execute", execute)
    query = "把双螺旋的另一条链写出来"
    assert registry.search(query) == ()
    assert registry.search(query, mode="semantic", limit=1)[0].name == RC
    assert [len(call["input"]) for call in endpoint] == [1, 2, 1]
    execute.assert_not_called()


def test_hybrid_combines_ranks_and_keeps_exact_name_first(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
) -> None:
    registry = _registry(config)
    # Semantic ranking favors RC even for this query; keyword evidence favors download.
    assert registry.search("download_genome", mode="semantic", limit=1)[0].name == RC
    assert registry.search("download_genome", mode="hybrid", limit=1)[0].name == DOWNLOAD
    hybrid = registry.search("download", mode="hybrid")
    assert hybrid[0].name == DOWNLOAD
    assert {tool.name for tool in hybrid} == {RC, DOWNLOAD, WEIGHT}


def test_filters_apply_before_embedding_and_blank_queries_are_local(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
) -> None:
    registry = _registry(config)
    assert len(registry.search(" ", mode="semantic")) == 3
    assert registry.search(QUERY, category="annotation", mode="semantic") == ()
    assert endpoint == []
    assert [tool.name for tool in registry.search(QUERY, category="ops", mode="semantic")] == [RC]
    assert len(endpoint) == 2
    assert all(DOWNLOAD not in text for call in endpoint for text in call["input"])


def test_persistent_cache_reuses_vectors_and_refreshes_only_changed_descriptions(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
) -> None:
    registry = _registry(config)
    registry.search(QUERY, mode="semantic")
    endpoint.clear()
    fresh = _registry(config)
    fresh.search(QUERY, mode="semantic")
    assert len(endpoint) == 1  # Only the query goes to the provider.
    endpoint.clear()
    tools = tuple(
        replace(t, description=t.description + " Updated.") if t.name == RC else t
        for t in fresh.tools
    )
    AgentToolRegistry(tools, embedding_config=config).search(QUERY, mode="semantic")
    assert len(endpoint) == 2 and endpoint[1]["input"][0].startswith(RC)
    assert config.cache_dir is not None
    cached = "".join(path.read_text() for path in Path(config.cache_dir).glob("*.json"))
    assert QUERY not in cached and "offline-test-key" not in cached


@pytest.mark.parametrize("change", [{"model": "different-model"}, {"dimensions": 3}])
def test_cache_isolated_by_model_and_dimensions(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
    change: dict[str, Any],
) -> None:
    _registry(config).search(QUERY, mode="semantic")
    endpoint.clear()
    _registry(replace(config, **change)).search(QUERY, mode="semantic")
    assert [len(call["input"]) for call in endpoint] == [1, 2, 1]


def test_corrupt_cache_is_refreshed_and_memory_only_mode_writes_nothing(
    config: EmbeddingConfig,
    endpoint: list[dict[str, Any]],
    tmp_path: Path,
) -> None:
    _registry(config).search(QUERY, mode="semantic")
    assert config.cache_dir is not None
    for path in Path(config.cache_dir).glob("*.json"):
        path.write_text("broken")
    endpoint.clear()
    _registry(config).search(QUERY, mode="semantic")
    assert len(endpoint) == 3
    memory = _registry(replace(config, cache_dir=None))
    memory.search(QUERY, mode="semantic")
    endpoint.clear()
    memory.search(QUERY, mode="semantic")
    assert len(endpoint) == 1
    assert sorted(path.name for path in tmp_path.iterdir()) == ["vectors"]


@pytest.mark.parametrize("mode", ["semantic", "hybrid"])
def test_cloud_config_required_only_when_requested(mode: Any) -> None:
    registry = _registry(None)
    assert registry.search("download")
    with pytest.raises(ConfigurationError, match="MISSING_EMBEDDING_CONFIG"):
        registry.search(QUERY, mode=mode)


@pytest.mark.parametrize(
    "options",
    [
        {"base_url": "http://remote.example/v1"},
        {"base_url": "https://user:secret@example.com"},
        {"base_url": "https://example.com?key=secret"},
        {"model": ""},
        {"batch_size": True},
        {"batch_size": 0},
        {"dimensions": 0},
        {"dimensions": True},
        {"cache_dir": 42},
        {"timeout_seconds": float("inf")},
        {"api_key_env": "INVALID-ENV"},
    ],
)
def test_invalid_configuration_fails_locally(
    config: EmbeddingConfig,
    options: dict[str, Any],
) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_EMBEDDING_CONFIG"):
        replace(config, **options)


def test_missing_key_does_not_borrow_deepseek_key(
    config: EmbeddingConfig,
    monkeypatch: pytest.MonkeyPatch,
    endpoint: list[dict[str, Any]],
) -> None:
    monkeypatch.delenv("DNAKIT_EMBEDDING_API_KEY", raising=False)
    monkeypatch.setenv("DEEPSEEK_API_KEY", "chat-key-must-not-be-used")
    with pytest.raises(ConfigurationError, match="MISSING_EMBEDDING_API_KEY"):
        _registry(replace(config, api_key=None)).search(QUERY, mode="semantic")
    assert endpoint == []


@pytest.mark.parametrize(
    "data",
    [
        [],
        [{"index": 1, "embedding": [1, 0]}],
        [{"index": True, "embedding": [1]}],
        [{"index": 0, "embedding": []}],
        [{"index": 0, "embedding": [0, 0]}],
        [{"index": 0, "embedding": [True]}],
        [{"index": 0, "embedding": ["1"]}],
        [{"index": 0, "embedding": [float("nan")]}],
        [{"index": 0, "embedding": [float("inf")]}],
    ],
)
def test_invalid_vectors_are_rejected(
    config: EmbeddingConfig,
    monkeypatch: pytest.MonkeyPatch,
    data: Any,
) -> None:
    raw = json.dumps({"data": data}).encode()
    monkeypatch.setattr(
        "dnakit.tools.embeddings.urllib.request.build_opener",
        lambda *a: Mock(open=lambda *a, **kw: io.BytesIO(raw)),
    )
    with pytest.raises(BackendExecutionError, match="INVALID_EMBEDDING_RESPONSE"):
        _embed([QUERY], config)


def test_response_indexes_and_batch_dimensions_are_validated(
    config: EmbeddingConfig,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for data in (
        [{"index": 0, "embedding": [1, 0]}, {"index": 0, "embedding": [1, 0]}],
        [{"index": 0, "embedding": [1, 0]}, {"index": 1, "embedding": [1, 0, 0]}],
    ):
        raw = json.dumps({"data": data}).encode()
        monkeypatch.setattr(
            "dnakit.tools.embeddings.urllib.request.build_opener",
            lambda *a, raw=raw: Mock(open=lambda *a, **kw: io.BytesIO(raw)),
        )
        with pytest.raises(BackendExecutionError, match="INVALID_EMBEDDING_RESPONSE"):
            _embed([QUERY, "second query"], config)


def test_byte_limits_and_cloud_errors_do_not_fall_back_to_keyword(
    config: EmbeddingConfig,
    monkeypatch: pytest.MonkeyPatch,
    endpoint: list[dict[str, Any]],
) -> None:
    with pytest.raises(InputFormatError, match="EMBEDDING_REQUEST_LIMIT"):
        _embed([QUERY], replace(config, max_request_bytes=1))
    assert endpoint == []
    with pytest.raises(BackendExecutionError, match="EMBEDDING_RESPONSE_LIMIT"):
        _embed([QUERY], replace(config, max_response_bytes=1))
    for error, code in (
        (TimeoutError("secret must not leak"), "EMBEDDING_API_UNAVAILABLE"),
        (
            urllib.error.HTTPError(
                config.endpoint, 401, "secret", Message(), io.BytesIO(b"secret")
            ),
            "EMBEDDING_API_ERROR",
        ),
    ):
        open_request = Mock(side_effect=error)
        monkeypatch.setattr(
            "dnakit.tools.embeddings.urllib.request.build_opener",
            lambda *a, opener=open_request: Mock(open=opener),
        )
        with pytest.raises(BackendExecutionError, match=code) as caught:
            _registry(config).search("reverse complement", mode="hybrid")
        assert "secret" not in str(caught.value)
        open_request.assert_called_once()


def test_embedding_config_representation_hides_credentials(config: EmbeddingConfig) -> None:
    assert "offline-test-key" not in repr(config)
    assert config.endpoint == replace(config, base_url=config.endpoint).endpoint


def test_real_http_transport_and_redirect_rejection(config: EmbeddingConfig) -> None:
    requests: list[tuple[str, str | None]] = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self) -> None:
            requests.append((self.path, self.headers.get("Authorization")))
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            if self.path == "/redirect/embeddings":
                self.send_response(307)
                self.send_header("Location", "/must-not-receive-key")
                self.end_headers()
                return
            data = [
                {"index": i, "embedding": [0, 1] if text.startswith(DOWNLOAD) else [1, 0]}
                for i, text in enumerate(body["input"])
            ]
            raw = json.dumps({"data": data}).encode()
            self.send_response(200)
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

        def log_message(self, format: str, *args: Any) -> None:
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        base = f"http://127.0.0.1:{server.server_port}"
        registry = _registry(replace(config, base_url=base + "/v1"))
        assert registry.search(QUERY, mode="semantic", limit=1)[0].name == RC
        assert len(requests) == 3
        assert all(key == "Bearer offline-test-key" for _, key in requests)
        with pytest.raises(BackendExecutionError, match="EMBEDDING_API_ERROR"):
            _embed([QUERY], replace(config, base_url=base + "/redirect"))
        assert requests[-1][0] == "/redirect/embeddings"
        assert all(path != "/must-not-receive-key" for path, _ in requests)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
