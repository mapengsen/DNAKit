"""BM25 reference scores, fusion, bounded reranking and transport contracts."""

from __future__ import annotations

import io
import json
import math
import threading
import urllib.error
from dataclasses import replace
from email.message import Message
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any
from unittest.mock import Mock

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.tools import (
    AgentToolRegistry,
    BM25Config,
    EmbeddingConfig,
    RerankerConfig,
    default_tool_registry,
)
from dnakit.tools.reranking import _rerank
from dnakit.tools.retrieval import _bm25_scores, _rrf


@pytest.fixture
def embedding() -> EmbeddingConfig:
    return EmbeddingConfig(
        model="offline",
        base_url="https://embedding.example/v1",
        api_key="embed-key",
        cache_dir=None,
    )


@pytest.fixture
def reranker() -> RerankerConfig:
    return RerankerConfig(
        model="offline-reranker", base_url="https://reranker.example/v2", api_key="rerank-key"
    )


def _registry(
    count: int = 3,
    *,
    embedding: EmbeddingConfig | None = None,
    reranker: RerankerConfig | None = None,
) -> AgentToolRegistry:
    template = default_tool_registry(verbose=False).resolve("dnakit.ops.reverse_complement")
    return AgentToolRegistry(
        tuple(
            replace(
                template,
                name=f"dnakit.ops.candidate_{i:03}",
                description="searchable capability",
                search_text=None,
            )
            for i in range(count)
        ),
        embedding_config=embedding,
        reranker_config=reranker,
    )


def test_bm25_matches_hand_calculated_reference_and_length_normalization() -> None:
    documents = [("a", "alpha alpha beta"), ("b", "alpha beta beta beta"), ("c", "gamma")]
    # N=3, df(alpha)=2, average length=8/3; Robertson positive IDF=log(1.6).
    scores = _bm25_scores("alpha", documents, BM25Config())
    assert scores == pytest.approx(
        {"a": math.log(1.6) * 4.4 / 3.3125, "b": math.log(1.6) * 2.2 / 2.65}
    )
    assert scores["a"] > scores["b"]
    assert _bm25_scores("alpha alpha", documents, BM25Config()) == scores
    assert _bm25_scores("unknown", documents, BM25Config()) == {}
    assert _bm25_scores("", documents, BM25Config()) == {}
    assert _bm25_scores("alpha", [], BM25Config()) == {}
    assert _bm25_scores("alpha", documents, BM25Config(b=0))["a"] == pytest.approx(
        math.log(1.6) * 4.4 / 3.2
    )


def test_bm25_identifiers_chinese_and_stable_ties() -> None:
    documents = [("a", "Reverse_Complement 反向互补"), ("b", "download_genome 参考基因组")]
    for query in ("reverse complement", "REVERSE_COMPLEMENT", "反向互补", "互补"):
        assert set(_bm25_scores(query, documents, BM25Config())) == {"a"}
    assert set(_bm25_scores("参考基因组", documents, BM25Config())) == {"b"}
    registry = _registry()
    assert registry.search("capability", mode="bm25") == registry.tools


def test_bm25_is_local_uses_search_text_and_filters_before_ranking(
    embedding: EmbeddingConfig, reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    fail = Mock(side_effect=AssertionError("Local search must not use model credentials or APIs"))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", fail)
    monkeypatch.setattr("dnakit.tools.registry._rerank", fail)
    monkeypatch.setattr(EmbeddingConfig, "resolved_api_key", fail)
    monkeypatch.setattr(RerankerConfig, "resolved_api_key", fail)
    tools = _registry().tools
    registry = AgentToolRegistry(
        (
            replace(tools[0], description="obsolete", search_text="updated capability"),
            replace(tools[1], agent_compatible=False),
            replace(tools[2], category="references"),
        ),
        embedding_config=embedding,
        reranker_config=reranker,
    )
    assert registry.search("updated", mode="bm25")[0].name == tools[0].name
    assert registry.search("obsolete", mode="bm25") == ()
    assert len(registry.search("capability", mode="bm25", category="ops")) == 1
    assert len(registry.search("capability", mode="bm25", include_python_only=True)) == 3
    fail.assert_not_called()


def test_bm25_hybrid_fuses_both_lists_without_exact_name_override(
    embedding: EmbeddingConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    registry = _registry(embedding=embedding)
    a, b, c = [tool.name for tool in registry.tools]
    assert registry._embedding_index is not None
    monkeypatch.setattr(registry._embedding_index, "rank", lambda *args: [b, c, a])
    assert _rrf(([a, b], [b, c, a])) == [b, a, c]
    # Query matches each tool equally except for a's exact identifier component.
    assert registry.search("candidate_000 capability", mode="bm25")[0].name == a
    assert registry.search("candidate_000 capability", mode="bm25_hybrid")[0].name == b
    assert registry.search("candidate_000", mode="hybrid")[0].name == a


def test_rerank_uses_fixed_pool_and_changes_top_result_without_execution(
    embedding: EmbeddingConfig, reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    registry = _registry(60, embedding=embedding, reranker=reranker)
    monkeypatch.setattr(
        "dnakit.tools.embeddings._embed", lambda texts, config: [(1.0,)] * len(texts)
    )
    execute = Mock(side_effect=AssertionError("Search must not execute candidates"))
    monkeypatch.setattr(registry, "execute", execute)
    calls = []

    def open_request(request: Any, *, timeout: float) -> io.BytesIO:
        assert request.full_url == "https://reranker.example/v2/rerank"
        assert request.get_header("Authorization") == "Bearer rerank-key"
        assert timeout == 60.0
        payload = json.loads(request.data)
        calls.append(payload)
        assert payload["model"] == "offline-reranker"
        assert payload["top_n"] == len(payload["documents"]) == 50
        # Return rows in original order, with the last recalled document scored highest.
        return io.BytesIO(
            json.dumps(
                {"results": [{"index": i, "relevance_score": i / 50} for i in range(50)]}
            ).encode()
        )

    monkeypatch.setattr(
        "dnakit.tools.reranking.urllib.request.build_opener",
        lambda *handlers: Mock(open=open_request),
    )
    first = registry.search("capability", mode="rerank", limit=1)
    ten = registry.search("capability", mode="rerank", limit=10)
    assert first == ten[:1] == (registry.tools[49],)
    assert calls[0] == calls[1]
    assert all("candidate_050" not in text for text in calls[0]["documents"])
    execute.assert_not_called()


@pytest.mark.parametrize("mode", ["bm25_hybrid", "rerank"])
def test_cloud_filters_and_blank_queries_make_no_requests(
    embedding: EmbeddingConfig,
    reranker: RerankerConfig,
    monkeypatch: pytest.MonkeyPatch,
    mode: Any,
) -> None:
    registry = _registry(embedding=embedding, reranker=reranker)
    fail = Mock(side_effect=AssertionError("No cloud work expected"))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", fail)
    monkeypatch.setattr("dnakit.tools.registry._rerank", fail)
    assert registry.search(" ", mode=mode) == registry.tools
    assert registry.search("capability", mode=mode, category="annotation") == ()
    fail.assert_not_called()


def test_rerank_filters_candidates_and_preserves_settings_on_refresh(
    embedding: EmbeddingConfig, reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    tools = _registry().tools
    bm25 = BM25Config(k1=1.5, b=0.5)
    registry = AgentToolRegistry(
        (
            tools[0],
            replace(tools[1], agent_compatible=False),
            replace(tools[2], category="references"),
        ),
        embedding_config=embedding,
        bm25_config=bm25,
        reranker_config=reranker,
    )._refreshed()
    assert registry.bm25_config is bm25
    assert registry.reranker_config is reranker
    assert registry.embedding_config is embedding
    embed = Mock(side_effect=lambda texts, config: [(1.0,)] * len(texts))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", embed)
    rerank = Mock(side_effect=lambda query, documents, config: [name for name, _ in documents])
    monkeypatch.setattr("dnakit.tools.registry._rerank", rerank)
    result = registry.search("capability", category="ops", mode="rerank")
    assert [tool.name for tool in result] == [tools[0].name]
    assert len(rerank.call_args.args[1]) == 1
    assert all(
        "candidate_001" not in str(call) and "candidate_002" not in str(call)
        for call in embed.call_args_list
    )
    rerank.side_effect = BackendExecutionError("Unavailable", code="RERANK_API_UNAVAILABLE")
    with pytest.raises(BackendExecutionError, match="RERANK_API_UNAVAILABLE"):
        registry.search("capability", category="ops", mode="rerank")


def test_missing_configuration_and_keys_fail_before_embedding(
    embedding: EmbeddingConfig, reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(ConfigurationError, match="MISSING_EMBEDDING_CONFIG"):
        _registry().search("capability", mode="bm25_hybrid")
    with pytest.raises(ConfigurationError, match="MISSING_RERANKER_CONFIG"):
        _registry(embedding=embedding).search("capability", mode="rerank")
    monkeypatch.delenv("DNAKIT_RERANK_API_KEY", raising=False)
    monkeypatch.setenv("DNAKIT_EMBEDDING_API_KEY", "not-a-reranker-key")
    monkeypatch.setenv("DNAKIT_LLM_API_KEY", "not-a-reranker-key")
    fail = Mock(side_effect=AssertionError("Do not pay for embedding before key validation"))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", fail)
    registry = _registry(embedding=embedding, reranker=replace(reranker, api_key=None))
    with pytest.raises(ConfigurationError, match="MISSING_RERANKER_API_KEY"):
        registry.search("capability", mode="rerank")
    with pytest.raises(ConfigurationError, match="AGENT_TOOL_LIMIT_INVALID"):
        registry.search("capability", mode="rerank", limit=51)
    fail.assert_not_called()


@pytest.mark.parametrize(
    "results",
    [
        [],
        [{"index": 1, "relevance_score": 0.5}],
        [{"index": -1, "relevance_score": 0.5}],
        [{"index": True, "relevance_score": 0.5}],
        [{"index": 0, "relevance_score": True}],
        [{"index": 0, "relevance_score": "0.5"}],
        [{"index": 0, "relevance_score": float("nan")}],
        [{"index": 0, "relevance_score": float("inf")}],
        [{"index": 0}],
    ],
)
def test_invalid_reranker_results_are_rejected(
    reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch, results: Any
) -> None:
    raw = json.dumps({"results": results}).encode()
    monkeypatch.setattr(
        "dnakit.tools.reranking.urllib.request.build_opener",
        lambda *handlers: Mock(open=lambda *args, **kwargs: io.BytesIO(raw)),
    )
    with pytest.raises(BackendExecutionError, match="INVALID_RERANK_RESPONSE"):
        _rerank("query", [("a", "document")], reranker)


def test_duplicate_indexes_rejected_and_ties_follow_recall_order(
    reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    documents = [("z", "first"), ("a", "second")]
    for indexes, valid in (([0, 0], False), ([1, 0], True)):
        raw = json.dumps(
            {"results": [{"index": i, "relevance_score": 0.5} for i in indexes]}
        ).encode()
        monkeypatch.setattr(
            "dnakit.tools.reranking.urllib.request.build_opener",
            lambda *handlers, raw=raw: Mock(open=lambda *args, **kwargs: io.BytesIO(raw)),
        )
        if valid:
            assert _rerank("query", documents, reranker) == ["z", "a"]
        else:
            with pytest.raises(BackendExecutionError, match="INVALID_RERANK_RESPONSE"):
                _rerank("query", documents, reranker)


def test_reranker_limits_errors_and_credentials_are_bounded(
    reranker: RerankerConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    documents = [("a", "document")]
    assert "rerank-key" not in repr(reranker)
    assert replace(reranker, base_url=reranker.endpoint).endpoint == reranker.endpoint
    opener = Mock(open=Mock(return_value=io.BytesIO(b"large response")))
    monkeypatch.setattr("dnakit.tools.reranking.urllib.request.build_opener", lambda *a: opener)
    with pytest.raises(InputFormatError, match="RERANK_REQUEST_LIMIT"):
        _rerank("query", documents, replace(reranker, max_request_bytes=1))
    opener.open.assert_not_called()
    with pytest.raises(BackendExecutionError, match="RERANK_RESPONSE_LIMIT"):
        _rerank("query", documents, replace(reranker, max_response_bytes=1))
    for error, code in (
        (TimeoutError("secret"), "RERANK_API_UNAVAILABLE"),
        (
            urllib.error.HTTPError(
                reranker.endpoint, 401, "secret", Message(), io.BytesIO(b"secret")
            ),
            "RERANK_API_ERROR",
        ),
    ):
        opener.open.reset_mock()
        opener.open.side_effect = error
        with pytest.raises(BackendExecutionError, match=code) as caught:
            _rerank("query", documents, reranker)
        assert "secret" not in str(caught.value)
        opener.open.assert_called_once()


@pytest.mark.parametrize(
    "options",
    [
        {"candidate_limit": 0},
        {"candidate_limit": True},
        {"candidate_limit": 1001},
        {"base_url": "http://remote.example"},
        {"base_url": "https://user:secret@example.com"},
        {"model": ""},
        {"timeout_seconds": float("inf")},
        {"api_key_env": "bad-env"},
    ],
)
def test_reranker_configuration_validation(
    reranker: RerankerConfig, options: dict[str, Any]
) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_RERANKER_CONFIG"):
        replace(reranker, **options)


@pytest.mark.parametrize(
    "options",
    [
        {"k1": 0},
        {"k1": True},
        {"k1": float("nan")},
        {"b": -1},
        {"b": 1.1},
        {"b": "0.5"},
    ],
)
def test_bm25_configuration_validation(options: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_BM25_CONFIG"):
        BM25Config(**options)


def test_real_reranker_http_and_redirect_rejection(reranker: RerankerConfig) -> None:
    paths: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self) -> None:
            paths.append(self.path)
            assert self.headers["Authorization"] == "Bearer rerank-key"
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            assert body["top_n"] == 2
            if self.path == "/redirect/rerank":
                self.send_response(307)
                self.send_header("Location", "/must-not-receive-key")
                self.end_headers()
                return
            raw = (
                b'{"results":[{"index":1,"relevance_score":0.9},{"index":0,"relevance_score":0.1}]}'
            )
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
        documents = [("a", "first"), ("b", "second")]
        assert _rerank("query", documents, replace(reranker, base_url=base + "/v2")) == ["b", "a"]
        with pytest.raises(BackendExecutionError, match="RERANK_API_ERROR"):
            _rerank("query", documents, replace(reranker, base_url=base + "/redirect"))
        assert paths == ["/v2/rerank", "/redirect/rerank"]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
