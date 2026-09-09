"""Offline tests for reusable prompt/input/output LLM calls and their privacy boundaries."""

from __future__ import annotations

import asyncio
import io
import json
import runpy
import threading
import traceback
import urllib.error
import urllib.request
from dataclasses import replace
from email.message import Message
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.llm import LLMConfig, LLMResult, acall_llm, call_llm


def _envelope(content: str = "Processed.", **overrides: object) -> dict[str, object]:
    return {
        "id": "completion-test",
        "model": "provider-model",
        "choices": [
            {"message": {"role": "assistant", "content": content}, "finish_reason": "stop"}
        ],
        "usage": {"prompt_tokens": 12, "completion_tokens": 4, "total_tokens": 16},
        **overrides,
    }


class _Response(io.BytesIO):
    pass


class _FakeOpener:
    def __init__(self, payload: object) -> None:
        self.payload = payload
        self.requests: list[urllib.request.Request] = []
        self.timeouts: list[float] = []

    def open(self, request: urllib.request.Request, *, timeout: float) -> _Response:
        self.requests.append(request)
        self.timeouts.append(timeout)
        if isinstance(self.payload, Exception):
            raise self.payload
        raw = (
            self.payload
            if isinstance(self.payload, bytes)
            else json.dumps(self.payload, ensure_ascii=False).encode("utf-8")
        )
        return _Response(raw)


@pytest.fixture
def config() -> LLMConfig:
    return LLMConfig(model="user-model", api_key="test-header-only-secret")


def _mock_api(monkeypatch: pytest.MonkeyPatch, payload: object) -> _FakeOpener:
    opener = _FakeOpener(payload)
    monkeypatch.setattr(urllib.request, "build_opener", lambda *args: opener)
    return opener


@pytest.mark.parametrize(
    ("base_url", "endpoint"),
    [
        ("https://api.deepseek.com", "https://api.deepseek.com/chat/completions"),
        ("https://api.example/v1/", "https://api.example/v1/chat/completions"),
        ("https://api.example/v1/chat/completions/", "https://api.example/v1/chat/completions"),
        ("http://127.0.0.1:8080/v1", "http://127.0.0.1:8080/v1/chat/completions"),
        ("http://[::1]:8080/v1", "http://[::1]:8080/v1/chat/completions"),
    ],
)
def test_llm_endpoint_normalization(base_url: str, endpoint: str) -> None:
    assert LLMConfig(model="m", base_url=base_url).endpoint == endpoint


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("model", ""),
        ("model", "m\nInjected"),
        ("base_url", "http://remote.example/v1"),
        ("base_url", "file:///private/file"),
        ("base_url", "https://user:password@api.example"),
        ("base_url", "https://api.example/v1?api_key=private"),
        ("base_url", "https://api.example/v1#private"),
        ("base_url", "https://api.example:bad/v1"),
        ("base_url", "https://api.example:0/v1"),
        ("base_url", "https://api.example\\@elsewhere/"),
        ("base_url", "https://api.example/ has-space"),
        ("api_key", "secret\nInjected"),
        ("api_key", "private token"),
        ("api_key_env", "INVALID=NAME"),
        ("temperature", float("nan")),
        ("temperature", True),
        ("temperature", 2.1),
        ("max_output_tokens", 0),
        ("max_request_bytes", -1),
        ("max_response_bytes", 1.5),
        ("timeout_seconds", 0),
        ("timeout_seconds", float("inf")),
        ("token_limit_parameter", "other"),
        ("instruction_role", "user"),
    ],
)
def test_invalid_llm_config(name: str, value: Any) -> None:
    kwargs = {"model": "test", name: value}
    with pytest.raises(ConfigurationError):
        LLMConfig(**kwargs)


def test_text_call_custom_prompt_data_and_secret_handling(
    monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    opener = _mock_api(monkeypatch, _envelope("Processed private input."))
    result = call_llm(
        "Summarize the input without inventing values.",
        input_data={"records": [{"score": 0.5}, {"score": 0.7}]},
        config=config,
        system_prompt="You summarize experiments.",
    )
    assert isinstance(result, LLMResult)
    assert result.text == "Processed private input."
    assert result.data is None
    assert result.model == "provider-model"
    assert result.usage["total_tokens"] == 16
    assert "private input" not in repr(result)
    assert "test-header-only-secret" not in repr(config)
    assert "test-header-only-secret" not in json.dumps(result.to_dict())
    request = opener.requests[0]
    assert request.get_method() == "POST"
    assert request.get_header("Authorization") == "Bearer test-header-only-secret"
    assert isinstance(request.data, bytes)
    assert b"test-header-only-secret" not in request.data
    payload = json.loads(request.data)
    assert payload["model"] == "user-model"
    assert payload["stream"] is False
    assert payload["max_tokens"] == 2048
    assert "temperature" not in payload
    assert "response_format" not in payload
    assert payload["messages"][0] == {"role": "system", "content": "You summarize experiments."}
    assert payload["messages"][1]["content"] == "Summarize the input without inventing values."
    assert '"score": 0.7' in payload["messages"][2]["content"]
    assert opener.timeouts == [60.0]


def test_json_output_and_openai_parameter_variant(
    monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    opener = _mock_api(monkeypatch, _envelope('{"labels": ["a", "b"], "count": 2}'))
    result = call_llm(
        "Return labels and count.",
        input_data="a, b",
        output_format="json",
        config=replace(
            config,
            base_url="https://api.openai.com/v1",
            instruction_role="developer",
            token_limit_parameter="max_completion_tokens",
            max_output_tokens=100,
            temperature=0.2,
        ),
    )
    assert result.data == {"labels": ["a", "b"], "count": 2}
    mutable = result.data
    assert mutable is not None
    mutable["labels"] = ["changed"]
    assert result.data == {"labels": ["a", "b"], "count": 2}
    assert isinstance(opener.requests[0].data, bytes)
    payload = json.loads(opener.requests[0].data)
    assert payload["response_format"] == {"type": "json_object"}
    assert payload["max_completion_tokens"] == 100
    assert "max_tokens" not in payload
    assert payload["temperature"] == 0.2
    assert payload["messages"][0]["role"] == "developer"
    assert "JSON object" in payload["messages"][0]["content"]


def test_key_is_resolved_at_call_time(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("DNAKIT_TEST_LLM_KEY", raising=False)
    configured = LLMConfig(model="m", api_key_env="DNAKIT_TEST_LLM_KEY")
    opener = _mock_api(monkeypatch, _envelope())
    with pytest.raises(ConfigurationError, match="MISSING_LLM_API_KEY"):
        call_llm("Task", config=configured)
    assert not opener.requests
    monkeypatch.setenv("DNAKIT_TEST_LLM_KEY", "first-key")
    call_llm("Task", config=configured)
    monkeypatch.setenv("DNAKIT_TEST_LLM_KEY", "rotated-key")
    call_llm("Task", config=configured)
    assert [request.get_header("Authorization") for request in opener.requests] == [
        "Bearer first-key",
        "Bearer rotated-key",
    ]


def test_input_paths_are_not_implicitly_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    path = tmp_path / "local.csv"
    path.write_text("private data must remain local", encoding="utf-8")
    opener = _mock_api(monkeypatch, _envelope())
    with pytest.raises(InputFormatError):
        call_llm("Task", input_data=path, config=config)
    assert not opener.requests
    call_llm("Task", input_data=str(path), config=config)
    assert isinstance(opener.requests[0].data, bytes)
    assert b"private data must remain local" not in opener.requests[0].data


@pytest.mark.parametrize(
    "input_data", [float("nan"), {1: "ambiguous key"}, {"unsupported"}, b"binary", object()]
)
def test_invalid_input_is_rejected_before_network(
    input_data: object, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    opener = _mock_api(monkeypatch, _envelope())
    with pytest.raises(InputFormatError):
        call_llm("Task", input_data=input_data, config=config)
    assert not opener.requests


def test_recursive_input_and_request_size_limit(
    monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    opener = _mock_api(monkeypatch, _envelope())
    recursive: list[object] = []
    recursive.append(recursive)
    with pytest.raises(InputFormatError):
        call_llm("Task", input_data=recursive, config=config)
    with pytest.raises(InputFormatError, match="LLM_REQUEST_LIMIT"):
        call_llm("Task", input_data="A" * 1000, config=replace(config, max_request_bytes=100))
    assert not opener.requests


@pytest.mark.parametrize(
    ("content", "output_format", "code"),
    [
        ("", "text", "EMPTY_LLM_OUTPUT"),
        ("```json\n{}\n```", "json", "INVALID_LLM_JSON_OUTPUT"),
        ("not JSON", "json", "INVALID_LLM_JSON_OUTPUT"),
        ("[]", "json", "INVALID_LLM_JSON_OUTPUT"),
        ('{"a":NaN}', "json", "INVALID_LLM_JSON_OUTPUT"),
        ('{"a":1,"a":2}', "json", "INVALID_LLM_JSON_OUTPUT"),
    ],
)
def test_invalid_output_never_becomes_success(
    content: str, output_format: Any, code: str, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    _mock_api(monkeypatch, _envelope(content))
    with pytest.raises(BackendExecutionError) as caught:
        call_llm("Task", config=config, output_format=output_format)
    assert caught.value.code == code


@pytest.mark.parametrize(
    ("choice", "code"),
    [
        ({"message": {"content": "partial"}, "finish_reason": "length"}, "LLM_OUTPUT_TRUNCATED"),
        ({"message": {"refusal": "private"}, "finish_reason": "stop"}, "LLM_REFUSED"),
        ({"message": {}, "finish_reason": "content_filter"}, "LLM_REFUSED"),
        (
            {"message": {"tool_calls": [{"function": "do_not_execute"}]}, "finish_reason": "stop"},
            "UNSUPPORTED_LLM_OUTPUT",
        ),
        ({"message": {"content": "text"}}, "UNSUPPORTED_LLM_OUTPUT"),
    ],
)
def test_truncation_refusals_and_tool_calls_are_explicit_errors(
    choice: dict[str, object], code: str, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    _mock_api(monkeypatch, _envelope(choices=[choice]))
    with pytest.raises(BackendExecutionError) as caught:
        call_llm("Task", config=config)
    assert caught.value.code == code


@pytest.mark.parametrize(
    "payload", [b"invalid json", b"\xff", [], {}, {"choices": []}, _envelope(usage=[])]
)
def test_malformed_envelopes(
    payload: object, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    _mock_api(monkeypatch, payload)
    with pytest.raises(BackendExecutionError, match="INVALID_LLM_RESPONSE"):
        call_llm("Task", config=config)


@pytest.mark.parametrize("status", [301, 401, 429, 500])
def test_http_errors_are_sanitized_and_not_retried(
    status: int, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    secret = "sensitive-server-echo"
    error = urllib.error.HTTPError(
        f"https://api.example/?key={secret}", status, secret, Message(), io.BytesIO(secret.encode())
    )
    opener = _mock_api(monkeypatch, error)
    with pytest.raises(BackendExecutionError) as caught:
        call_llm("private input", config=config)
    assert caught.value.context == {"status": status}
    assert secret not in "".join(traceback.format_exception(caught.value))
    assert "private input" not in str(caught.value)
    assert len(opener.requests) == 1


@pytest.mark.parametrize("error", [TimeoutError("secret"), urllib.error.URLError("secret")])
def test_network_errors_and_response_limit(
    error: Exception, monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    _mock_api(monkeypatch, error)
    with pytest.raises(BackendExecutionError, match="LLM_API_UNAVAILABLE") as caught:
        call_llm("Task", config=config)
    assert "secret" not in "".join(traceback.format_exception(caught.value))
    _mock_api(monkeypatch, b"A" * 200)
    with pytest.raises(BackendExecutionError, match="LLM_RESPONSE_LIMIT"):
        call_llm("Task", config=replace(config, max_response_bytes=100))


def test_async_api_leaves_event_loop_free(
    monkeypatch: pytest.MonkeyPatch, config: LLMConfig
) -> None:
    loop_thread_id = threading.get_ident()
    request_threads: list[int] = []

    class ThreadOpener(_FakeOpener):
        def open(self, request: urllib.request.Request, *, timeout: float) -> _Response:
            request_threads.append(threading.get_ident())
            return super().open(request, timeout=timeout)

    monkeypatch.setattr(urllib.request, "build_opener", lambda *args: ThreadOpener(_envelope()))

    async def run() -> list[LLMResult]:
        return list(
            await asyncio.gather(acall_llm("A", config=config), acall_llm("B", config=config))
        )

    results = asyncio.run(run())
    assert [result.text for result in results] == ["Processed.", "Processed."]
    assert len(request_threads) == 2
    assert all(thread_id != loop_thread_id for thread_id in request_threads)


@pytest.mark.parametrize("status", [301, 302, 303, 307, 308])
def test_real_http_transport_refuses_redirects(config: LLMConfig, status: int) -> None:
    received: list[str] = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self) -> None:
            received.append(self.path)
            self.rfile.read(int(self.headers.get("Content-Length", "0")))
            self.send_response(status)
            self.send_header("Location", "/must-not-receive-credentials")
            self.end_headers()

        def do_GET(self) -> None:
            received.append(self.path)
            self.send_response(200)
            self.end_headers()

        def log_message(self, format: str, *args: object) -> None:
            pass

    with ThreadingHTTPServer(("127.0.0.1", 0), Handler) as server:
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            local = replace(config, base_url=f"http://127.0.0.1:{server.server_port}")
            with pytest.raises(BackendExecutionError, match="LLM_API_ERROR"):
                call_llm("Local fixture only", config=local)
            assert received == ["/chat/completions"]
        finally:
            server.shutdown()
            thread.join(timeout=5)


def test_editable_example_uses_one_explicit_call(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    opener = _mock_api(
        monkeypatch,
        _envelope(
            '{"items": [{"id": "sample-a", "summary": "Passed"}, '
            '{"id": "sample-b", "summary": "Repeat"}]}'
        ),
    )
    monkeypatch.setenv("DNAKIT_LLM_API_KEY", "example-secret")
    monkeypatch.setenv("DNAKIT_LLM_MODEL", "example-model")
    example = Path(__file__).parents[2] / "examples" / "custom_llm.py"
    runpy.run_path(str(example), run_name="__main__")
    assert len(opener.requests) == 1
    output = capsys.readouterr().out
    assert "sample-a Passed" in output
    assert "sample-b Repeat" in output
    assert "example-secret" not in output
