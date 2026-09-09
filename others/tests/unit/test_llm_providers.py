"""Offline provider contracts for the general-purpose LLM Python API."""

from __future__ import annotations

import asyncio
import io
import json
import urllib.request
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.llm import (
    LLMConfig,
    acall_llm,
    available_llm_providers,
    call_llm,
    get_llm_provider,
)

_PROVIDERS = (
    ("deepseek", "https://api.deepseek.com", "DEEPSEEK_API_KEY", "max_tokens"),
    (
        "gemini",
        "https://generativelanguage.googleapis.com/v1beta/openai",
        "GEMINI_API_KEY",
        "max_tokens",
    ),
    (
        "qwen",
        "https://dashscope.aliyuncs.com/compatible-mode/v1",
        "DASHSCOPE_API_KEY",
        "max_tokens",
    ),
    ("kimi", "https://api.moonshot.cn/v1", "MOONSHOT_API_KEY", "max_tokens"),
    ("openrouter", "https://openrouter.ai/api/v1", "OPENROUTER_API_KEY", "max_tokens"),
    ("openai", "https://api.openai.com/v1", "OPENAI_API_KEY", "max_completion_tokens"),
)


def _mock_api(monkeypatch: pytest.MonkeyPatch, content: str) -> list[urllib.request.Request]:
    requests: list[urllib.request.Request] = []

    def capture(request: urllib.request.Request, *, timeout: float) -> io.BytesIO:
        assert timeout == 60.0
        requests.append(request)
        return io.BytesIO(
            json.dumps(
                {
                    "choices": [{"message": {"content": content}, "finish_reason": "stop"}],
                    "usage": {"total_tokens": 42},
                }
            ).encode()
        )

    monkeypatch.setattr(urllib.request, "build_opener", lambda *args: SimpleNamespace(open=capture))
    return requests


@pytest.mark.parametrize(("provider", "base_url", "key_env", "token_field"), _PROVIDERS)
@pytest.mark.parametrize("output_format", ["text", "json"])
def test_provider_request_contract(
    provider: str,
    base_url: str,
    key_env: str,
    token_field: str,
    output_format: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(key_env, "provider-specific-test-key")
    # An unrelated general key must not override a selected provider's key.
    monkeypatch.setenv("DNAKIT_LLM_API_KEY", "wrong-generic-key")
    requests = _mock_api(monkeypatch, '{"summary": "ok"}' if output_format == "json" else "ok")
    config = LLMConfig(provider=provider, model="user-selected-model")
    assert config.base_url == base_url
    assert config.api_key_env == key_env
    result = call_llm(
        "Summarize", input_data={"value": 1}, config=config, output_format=output_format
    )
    assert result.data == ({"summary": "ok"} if output_format == "json" else None)
    assert result.model == "user-selected-model"
    request = requests[0]
    assert request.full_url == base_url + "/chat/completions"
    assert request.get_header("Authorization") == "Bearer provider-specific-test-key"
    assert isinstance(request.data, bytes)
    assert b"provider-specific-test-key" not in request.data
    assert b"wrong-generic-key" not in request.data
    payload = json.loads(request.data)
    assert payload["model"] == "user-selected-model"
    assert payload[token_field] == 2048
    assert {"max_tokens", "max_completion_tokens"} & payload.keys() == {token_field}
    assert payload["stream"] is False
    assert "temperature" not in payload
    assert payload.get("response_format") == (
        {"type": "json_object"} if output_format == "json" else None
    )


def test_registry_and_legacy_defaults_without_network() -> None:
    assert {item.name for item in available_llm_providers()} == {
        "deepseek",
        "gemini",
        "qwen",
        "kimi",
        "openrouter",
        "openai",
        "custom",
    }
    legacy = LLMConfig(model="unchanged")
    assert legacy.base_url == "https://api.deepseek.com"
    assert legacy.api_key_env == "DNAKIT_LLM_API_KEY"
    assert legacy.resolved_token_limit_parameter == "max_tokens"
    assert LLMConfig(model="m", provider="deepseek").api_key_env == "DEEPSEEK_API_KEY"


@pytest.mark.parametrize(
    ("alias", "canonical"),
    [
        ("Google", "gemini"),
        ("google-gemini", "gemini"),
        ("dashscope", "qwen"),
        ("alibaba", "qwen"),
        ("moonshot", "kimi"),
        ("openai-compatible", "custom"),
    ],
)
def test_provider_aliases(alias: str, canonical: str) -> None:
    assert get_llm_provider(alias).name == canonical
    config = LLMConfig(model="m", provider=alias, base_url="https://proxy.example/v1")
    assert config.provider == canonical


@pytest.mark.parametrize("provider", ["", "unknown", "openai-typo", 5, []])
def test_unknown_providers_fail_without_falling_back(provider: Any) -> None:
    with pytest.raises(ConfigurationError):
        LLMConfig(provider=provider, model="m")


@pytest.mark.parametrize(
    ("provider", "base_url"),
    [
        ("qwen", "https://workspace.cn-beijing.maas.aliyuncs.com/compatible-mode/v1"),
        ("qwen", "https://workspace.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1"),
        ("kimi", "https://api.moonshot.ai/v1"),
        ("custom", "http://127.0.0.1:8000/v1/chat/completions/"),
    ],
)
def test_custom_region_and_key_overrides(
    provider: str, base_url: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests = _mock_api(monkeypatch, "ok")
    monkeypatch.setenv("PROJECT_MODEL_KEY", "project-key")
    config = LLMConfig(
        provider=provider,
        model="my-model",
        base_url=base_url,
        api_key_env="PROJECT_MODEL_KEY",
        token_limit_parameter="max_tokens",
    )
    call_llm("Task", config=config)
    normalized = base_url.rstrip("/")
    expected = (
        normalized if normalized.endswith("/chat/completions") else normalized + "/chat/completions"
    )
    assert requests[0].full_url == expected
    assert requests[0].get_header("Authorization") == "Bearer project-key"


def test_custom_provider_needs_url() -> None:
    with pytest.raises(ConfigurationError, match="base_url"):
        LLMConfig(provider="custom", model="m")


def test_openai_url_auto_detection_and_explicit_token_override() -> None:
    assert (
        LLMConfig(model="m", base_url="https://api.openai.com/v1").resolved_token_limit_parameter
        == "max_completion_tokens"
    )
    assert (
        LLMConfig(
            model="m", provider="openai", base_url="https://proxy.example/v1"
        ).resolved_token_limit_parameter
        == "max_completion_tokens"
    )
    assert (
        LLMConfig(
            model="m", provider="openai", token_limit_parameter="max_tokens"
        ).resolved_token_limit_parameter
        == "max_tokens"
    )


@pytest.mark.parametrize(
    ("provider", "extra"),
    [
        ("qwen", {"enable_thinking": False}),
        ("kimi", {"reasoning_effort": "low"}),
        ("gemini", {"extra_body": {"google": {"thinking_config": {"thinking_level": "low"}}}}),
        ("openrouter", {"provider": {"require_parameters": True}}),
        ("openai", {"reasoning_effort": "low"}),
    ],
)
def test_provider_extensions_are_preserved_and_frozen(
    provider: str, extra: dict[str, object], monkeypatch: pytest.MonkeyPatch
) -> None:
    requests = _mock_api(monkeypatch, "ok")
    config = LLMConfig(provider=provider, model="m", api_key="test-key", extra_body=extra)
    expected = json.loads(json.dumps(extra))
    extra["changed_after_configuration"] = True
    call_llm("Task", config=config)
    assert isinstance(requests[0].data, bytes)
    payload = json.loads(requests[0].data)
    for key, value in expected.items():
        assert payload[key] == value
    assert "changed_after_configuration" not in payload


@pytest.mark.parametrize(
    "field",
    [
        "model",
        "messages",
        "stream",
        "n",
        "models",
        "max_tokens",
        "max_completion_tokens",
        "response_format",
        "temperature",
        "Authorization",
        "api_key",
        "headers",
        "tools",
        "tool_choice",
        "functions",
        "modalities",
        "audio",
    ],
)
def test_extensions_cannot_replace_core_or_enable_unsupported_features(field: str) -> None:
    with pytest.raises(ConfigurationError):
        LLMConfig(model="m", extra_body={field: "override"})


@pytest.mark.parametrize("extra", [{"option": float("nan")}, {1: "non-string-key"}, []])
def test_invalid_extension_data(extra: Any) -> None:
    with pytest.raises(ConfigurationError):
        LLMConfig(model="m", extra_body=extra)


def test_extensions_are_included_in_request_size_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    requests = _mock_api(monkeypatch, "ok")
    config = LLMConfig(
        model="m",
        api_key="test-key",
        extra_body={"metadata": {"text": "A" * 1000}},
        max_request_bytes=500,
    )
    with pytest.raises(InputFormatError, match="LLM_REQUEST_LIMIT"):
        call_llm("Task", config=config)
    assert not requests


def test_async_calls_use_provider_configuration(monkeypatch: pytest.MonkeyPatch) -> None:
    requests = _mock_api(monkeypatch, "ok")
    config = LLMConfig(provider="gemini", model="selected-gemini-model", api_key="local-test-key")
    result = asyncio.run(acall_llm("Task", config=config))
    assert result.text == "ok"
    assert (
        requests[0].full_url
        == "https://generativelanguage.googleapis.com/v1beta/openai/chat/completions"
    )
