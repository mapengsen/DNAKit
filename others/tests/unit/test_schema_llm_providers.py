"""Offline contracts for provider presets, request compatibility, and data privacy."""

from __future__ import annotations

import io
import json
import traceback
import urllib.error
import urllib.request
from collections.abc import Mapping
from email.message import Message
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.predictions import (
    LLMSchemaConfig,
    OpenAICompatibleSchemaClient,
    available_schema_llm_providers,
    get_schema_llm_provider,
    inspect_property_data,
)
from dnakit.predictions.custom_training import llm as llm_module

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
_PROPOSAL = {
    "format": "csv",
    "sequence_field": "sequence",
    "label_fields": ["label"],
    "id_field": "id",
}


def _response(content: object = None, *, finish_reason: str = "stop") -> io.BytesIO:
    return io.BytesIO(
        json.dumps(
            {
                "choices": [
                    {
                        "finish_reason": finish_reason,
                        "message": {
                            "content": json.dumps(_PROPOSAL) if content is None else content
                        },
                    }
                ]
            }
        ).encode()
    )


@pytest.mark.parametrize(("provider", "base_url", "key_env", "token_field"), _PROVIDERS)
def test_provider_requests_and_profile_privacy(
    provider: str,
    base_url: str,
    key_env: str,
    token_field: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv(key_env, "unit-test-key")
    calls: list[dict[str, Any]] = []

    def capture(request: urllib.request.Request, *, timeout: float) -> io.BytesIO:
        assert request.full_url == base_url + "/chat/completions"
        assert request.get_header("Authorization") == "Bearer unit-test-key"
        assert request.get_method() == "POST"
        assert timeout == 60.0
        assert isinstance(request.data, bytes)
        assert b"unit-test-key" not in request.data
        calls.append(json.loads(request.data))
        return _response()

    monkeypatch.setattr(llm_module, "_open_request", capture)
    config = LLMSchemaConfig(provider=provider, model="account-selected-model")
    assert config.base_url == base_url
    assert config.api_key_env == key_env
    source = tmp_path / "data.csv"
    source.write_text(
        "id,sequence,label\nr1,ACGT,1\nr2,TGCA,private-last-label\n", encoding="utf-8"
    )
    inspection = inspect_property_data(source, llm=config)
    assert inspection.mapping_origin == "llm"
    assert inspection.shared_sample_records == 1
    assert inspection.llm_model == "account-selected-model"
    assert len(calls) == 1
    payload = calls[0]
    assert payload["model"] == "account-selected-model"
    assert payload["response_format"] == {"type": "json_object"}
    assert payload[token_field] == 2_000
    assert {"max_tokens", "max_completion_tokens"} & payload.keys() == {token_field}
    assert payload["stream"] is False
    assert "temperature" not in payload
    assert "private-last-label" not in json.dumps(payload)
    assert "unit-test-key" not in json.dumps(inspection.to_dict())


def test_provider_registry_and_backward_compatible_defaults() -> None:
    assert {provider.name for provider in available_schema_llm_providers()} == {
        "deepseek",
        "gemini",
        "qwen",
        "kimi",
        "openrouter",
        "openai",
        "custom",
    }
    config = LLMSchemaConfig()
    assert config.model == "deepseek-v4-flash"
    assert config.base_url == "https://api.deepseek.com"
    assert config.api_key_env == "DEEPSEEK_API_KEY"
    assert config.resolved_token_limit_parameter == "max_tokens"


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
    assert get_schema_llm_provider(alias).name == canonical


def test_custom_and_regional_overrides(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PROJECT_LLM_KEY", "environment-key")
    config = LLMSchemaConfig(
        provider="custom",
        model="local-model",
        base_url="http://127.0.0.1:8000/v1/chat/completions/",
        api_key_env="PROJECT_LLM_KEY",
        api_key="explicit-key",
        token_limit_parameter="max_completion_tokens",
    )
    assert config.endpoint == "http://127.0.0.1:8000/v1/chat/completions"
    assert config.resolved_api_key() == "explicit-key"
    assert "explicit-key" not in repr(config)
    assert config.resolved_token_limit_parameter == "max_completion_tokens"
    kimi = LLMSchemaConfig(
        provider="kimi", model="chosen-model", base_url="https://api.moonshot.ai/v1"
    )
    assert kimi.endpoint == "https://api.moonshot.ai/v1/chat/completions"
    qwen = LLMSchemaConfig(
        provider="qwen",
        model="chosen-model",
        base_url="https://workspace.cn-beijing.maas.aliyuncs.com/compatible-mode/v1",
    )
    assert qwen.endpoint.endswith("/compatible-mode/v1/chat/completions")
    legacy = LLMSchemaConfig(model="chosen-model", base_url="https://api.openai.com/v1")
    assert legacy.resolved_token_limit_parameter == "max_completion_tokens"


@pytest.mark.parametrize(
    "options",
    [
        {"provider": "unknown"},
        {"provider": None},
        {"provider": "gemini"},
        {"provider": "custom", "model": "chosen-model"},
        {"temperature": -1},
        {"temperature": 3},
        {"temperature": True},
        {"temperature": float("nan")},
        {"json_mode": "true"},
        {"token_limit_parameter": "unknown"},
        {"token_limit_parameter": []},
        {"base_url": "https://user:password@example.com"},
        {"base_url": "https://example.com/v1?api_key=secret"},
        {"api_key_env": "INVALID-NAME"},
        {"api_key": "key\nheader"},
        {"extra_body": {"model": "substitute"}},
        {"extra_body": {"messages": []}},
        {"extra_body": {"Authorization": "secret"}},
        {"extra_body": {"stream": True}},
        {"extra_body": {"max_tokens": 9999}},
        {"extra_body": {"temperature": 0}},
        {"extra_body": {"response_format": {}}},
        {"extra_body": {"api_key": "secret"}},
        {"extra_body": {"value": object()}},
        {"extra_body": {"value": float("nan")}},
        {"extra_body": []},
    ],
)
def test_invalid_provider_options(options: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        LLMSchemaConfig(**options)


def test_provider_keys_are_not_read_from_another_provider(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("DEEPSEEK_API_KEY", "not-a-gemini-key")
    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    with pytest.raises(ConfigurationError) as error:
        LLMSchemaConfig(provider="gemini", model="chosen-model").resolved_api_key()
    assert error.value.code == "MISSING_LLM_API_KEY"
    assert "not-a-gemini-key" not in str(error.value)


def test_extra_body_controls_and_optional_json_mode(monkeypatch: pytest.MonkeyPatch) -> None:
    extra: dict[str, Any] = {"thinking": {"type": "disabled"}}
    config = LLMSchemaConfig(
        provider="kimi",
        model="chosen-model",
        api_key="unit-test-key",
        extra_body=extra,
        temperature=0.6,
        json_mode=False,
    )
    extra["thinking"]["type"] = "enabled"
    captured: list[dict[str, Any]] = []

    def capture(request: urllib.request.Request, *, timeout: float) -> io.BytesIO:
        assert isinstance(request.data, bytes)
        captured.append(json.loads(request.data))
        return _response("```json\n" + json.dumps(_PROPOSAL) + "\n```")

    monkeypatch.setattr(llm_module, "_open_request", capture)
    assert OpenAICompatibleSchemaClient().analyze({}, config=config) == _PROPOSAL
    assert captured[0]["thinking"] == {"type": "disabled"}
    assert captured[0]["temperature"] == 0.6
    assert "extra_body" not in captured[0]
    assert "response_format" not in captured[0]


@pytest.mark.parametrize(
    ("content", "reason", "error_code"),
    [
        ("{}", "length", "LLM_SCHEMA_OUTPUT_TRUNCATED"),
        ("{}", "content_filter", "LLM_SCHEMA_REFUSED"),
        ("[]", "stop", "INVALID_LLM_SCHEMA_RESPONSE"),
        ("explanation before {}", "stop", "INVALID_LLM_SCHEMA_RESPONSE"),
        ("", "stop", "INVALID_LLM_SCHEMA_RESPONSE"),
        ({"not": "text"}, "stop", "INVALID_LLM_SCHEMA_RESPONSE"),
    ],
)
def test_incomplete_or_invalid_responses(
    content: object,
    reason: str,
    error_code: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_open(request: Any, *, timeout: float) -> io.BytesIO:
        return _response(content, finish_reason=reason)

    monkeypatch.setattr(llm_module, "_open_request", fake_open)
    with pytest.raises(BackendExecutionError) as error:
        OpenAICompatibleSchemaClient().analyze({}, config=LLMSchemaConfig(api_key="unit-test-key"))
    assert error.value.code == error_code


@pytest.mark.parametrize("status", [302, 400, 401, 429, 500])
def test_http_errors_do_not_retry_or_echo_secrets(
    status: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = 0
    config = LLMSchemaConfig(api_key="unit-test-key")

    def fake_open(request: Any, *, timeout: float) -> io.BytesIO:
        nonlocal calls
        calls += 1
        raise urllib.error.HTTPError(
            request.full_url, status, "server echoed unit-test-key", Message(), None
        )

    monkeypatch.setattr(llm_module, "_open_request", fake_open)
    with pytest.raises(BackendExecutionError) as error:
        OpenAICompatibleSchemaClient().analyze({}, config=config)
    assert calls == 1
    assert error.value.code == "LLM_SCHEMA_API_ERROR"
    assert "unit-test-key" not in "".join(traceback.format_exception(error.value))


def test_redirects_cannot_forward_api_key() -> None:
    request = urllib.request.Request("https://api.openai.com/v1/chat/completions")
    request.add_header("Authorization", "Bearer unit-test-key")
    handler = llm_module._NoRedirectHandler()
    assert (
        handler.redirect_request(
            request, None, 302, "Found", {}, "https://another-host.invalid/chat/completions"
        )
        is None
    )


def test_structure_only_mode_with_a_provider(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def capture(request: urllib.request.Request, *, timeout: float) -> io.BytesIO:
        assert isinstance(request.data, bytes)
        assert b"private-sequence-id" not in request.data
        assert b"private-label" not in request.data
        payload = json.loads(request.data)
        assert "JSON" in payload["messages"][0]["content"]
        return _response()

    monkeypatch.setattr(llm_module, "_open_request", capture)
    source = tmp_path / "source.csv"
    source.write_text(
        "id,sequence,label\nprivate-sequence-id,ACGT,private-label\n", encoding="utf-8"
    )
    result = inspect_property_data(
        source,
        llm=LLMSchemaConfig(
            provider="openrouter",
            model="vendor/chosen-model",
            api_key="unit-test-key",
            include_sample_values=False,
        ),
    )
    # The audit counts the type-only sample, whose original values were not sent.
    assert result.shared_sample_records == 1
    assert isinstance(result.mapping.to_dict(), Mapping)
