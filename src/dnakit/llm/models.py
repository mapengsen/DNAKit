"""Configuration and results for general-purpose, explicitly invoked LLM tasks."""

from __future__ import annotations

import math
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, cast
from urllib.parse import urlsplit

from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError

from .providers import get_llm_provider

OutputFormat = Literal["text", "json"]


def _single_line(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ConfigurationError(
            f"{name} must be non-empty single-line text.", code="INVALID_LLM_CONFIG"
        )
    return value.strip()


def _api_key(value: object) -> str:
    key = _single_line(value, "api_key")
    if not key.isascii() or any(character.isspace() for character in key):
        raise ConfigurationError(
            "api_key must be ASCII text without whitespace.", code="INVALID_LLM_CONFIG"
        )
    return key


@dataclass(frozen=True, slots=True)
class LLMConfig:
    """Configure an OpenAI-compatible Chat Completions endpoint without an SDK.

    ``model`` is the provider's model identifier, not a DNA embedding model.
    Credentials are resolved only when making a request. Use ``api_key_env`` to
    avoid storing a credential in code. Remote endpoints require HTTPS; HTTP is
    allowed only for localhost. Requests never follow redirects.

    Select ``provider`` to fill the API URL and provider-specific key environment
    variable. Omitting it retains the legacy DeepSeek URL + DNAKIT_LLM_API_KEY.
    Explicit URL/key-variable overrides take precedence. ``extra_body`` supplies
    model-specific JSON parameters without replacing core request fields.
    """

    model: str
    base_url: str | None = None
    api_key: str | None = field(default=None, repr=False)
    api_key_env: str | None = None
    temperature: float | None = None
    max_output_tokens: int = 2_048
    token_limit_parameter: Literal["auto", "max_tokens", "max_completion_tokens"] = "auto"
    instruction_role: Literal["system", "developer"] = "system"
    timeout_seconds: float = 60.0
    max_request_bytes: int = 100_000
    max_response_bytes: int = 2_000_000
    provider: str | None = None
    extra_body: Mapping[str, object] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "model", _single_line(self.model, "model"))
        preset = get_llm_provider(self.provider if self.provider is not None else "deepseek")
        key_env = self.api_key_env
        if key_env is None:
            key_env = preset.api_key_env if self.provider is not None else "DNAKIT_LLM_API_KEY"
        object.__setattr__(self, "api_key_env", key_env)
        object.__setattr__(self, "provider", preset.name)
        requested_url = self.base_url if self.base_url is not None else preset.base_url
        if requested_url is None:
            raise ConfigurationError(
                "The custom provider requires an explicit base_url.", code="INVALID_LLM_CONFIG"
            )
        base_url = _single_line(requested_url, "base_url").rstrip("/")
        try:
            parsed = urlsplit(base_url)
            port = parsed.port
            valid_url = (
                bool(parsed.hostname)
                and parsed.username is None
                and parsed.password is None
                and not parsed.query
                and not parsed.fragment
                and not any(character.isspace() for character in base_url)
                and "\\" not in base_url
                and (port is None or port > 0)
                and (
                    parsed.scheme == "https"
                    or (
                        parsed.scheme == "http"
                        and parsed.hostname in {"localhost", "127.0.0.1", "::1"}
                    )
                )
            )
        except ValueError:
            valid_url = False
        if not valid_url:
            raise ConfigurationError(
                "base_url must be an HTTPS API URL (or localhost HTTP) without credentials, "
                "query parameters, or fragments.",
                code="INVALID_LLM_CONFIG",
            )
        object.__setattr__(self, "base_url", base_url)
        if self.api_key is not None:
            object.__setattr__(self, "api_key", _api_key(self.api_key))
        if not isinstance(self.api_key_env, str) or not re.fullmatch(
            r"[A-Za-z_][A-Za-z0-9_]*", self.api_key_env
        ):
            raise ConfigurationError(
                "api_key_env must be an environment-variable name.", code="INVALID_LLM_CONFIG"
            )
        if self.token_limit_parameter not in ("auto", "max_tokens", "max_completion_tokens"):
            raise ConfigurationError(
                "token_limit_parameter must be auto, max_tokens, or max_completion_tokens.",
                code="INVALID_LLM_CONFIG",
            )
        if self.instruction_role not in ("system", "developer"):
            raise ConfigurationError(
                "instruction_role must be system or developer.", code="INVALID_LLM_CONFIG"
            )
        for name, value, lower, upper in (
            ("max_output_tokens", self.max_output_tokens, 1, 1_000_000),
            ("max_request_bytes", self.max_request_bytes, 1, 100_000_000),
            ("max_response_bytes", self.max_response_bytes, 1, 100_000_000),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
                raise ConfigurationError(
                    f"{name} must be an integer in [{lower}, {upper}].", code="INVALID_LLM_CONFIG"
                )
        if (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, (int, float))
            or not math.isfinite(self.timeout_seconds)
            or not 0 < self.timeout_seconds <= 600
        ):
            raise ConfigurationError(
                "timeout_seconds must be a finite number in (0, 600].", code="INVALID_LLM_CONFIG"
            )
        if self.temperature is not None and (
            isinstance(self.temperature, bool)
            or not isinstance(self.temperature, (int, float))
            or not math.isfinite(self.temperature)
            or not 0 <= self.temperature <= 2
        ):
            raise ConfigurationError(
                "temperature must be None or a finite number in [0, 2].", code="INVALID_LLM_CONFIG"
            )
        protected = {
            "model",
            "models",
            "messages",
            "stream",
            "stream_options",
            "n",
            "max_tokens",
            "max_completion_tokens",
            "temperature",
            "response_format",
            "tools",
            "tool_choice",
            "functions",
            "function_call",
            "modalities",
            "audio",
            "authorization",
            "headers",
            "api_key",
            "api_key_env",
            "base_url",
        }
        if not isinstance(self.extra_body, Mapping) or any(
            not isinstance(key, str) or key.lower() in protected for key in self.extra_body
        ):
            raise ConfigurationError(
                "extra_body must be a JSON mapping without core request/authentication fields "
                "or unsupported tool/multimodal options.",
                code="INVALID_LLM_CONFIG",
            )
        try:
            object.__setattr__(self, "extra_body", FrozenDict(self.extra_body))
        except ConfigurationError:
            raise ConfigurationError(
                "extra_body must contain finite, bounded JSON-compatible values.",
                code="INVALID_LLM_CONFIG",
            ) from None

    @property
    def endpoint(self) -> str:
        """Append /chat/completions once, preserving optional /v1 or proxy prefixes."""

        assert self.base_url is not None  # Resolved and validated during construction.
        if self.base_url.endswith("/chat/completions"):
            return self.base_url
        return f"{self.base_url}/chat/completions"

    @property
    def resolved_token_limit_parameter(self) -> Literal["max_tokens", "max_completion_tokens"]:
        """Use OpenAI's current token field, including manually configured OpenAI URLs."""

        if self.token_limit_parameter != "auto":
            return self.token_limit_parameter
        if urlsplit(self.endpoint).hostname == "api.openai.com":
            return "max_completion_tokens"
        assert self.provider is not None
        return get_llm_provider(self.provider).token_limit_parameter

    def resolved_api_key(self) -> str:
        """Read and validate the key without logging its value."""

        assert self.api_key_env is not None  # Resolved and validated during construction.
        value = self.api_key if self.api_key is not None else os.environ.get(self.api_key_env)
        if value is None or not value.strip():
            raise ConfigurationError(
                "No API key is configured for this LLM request.",
                code="MISSING_LLM_API_KEY",
                hint=f"Set {self.api_key_env} or pass api_key explicitly.",
            )
        return _api_key(value)


@dataclass(frozen=True, slots=True)
class LLMResult:
    """A completed response with text, optional parsed JSON, and provider token usage.

    Prompt/input data and configuration are not retained. Text and parsed output
    are omitted from repr, but may still contain sensitive data supplied by the
    caller. ``data`` returns a fresh ordinary dict for downstream Python code.
    """

    text: str = field(repr=False)
    model: str
    response_id: str | None
    finish_reason: str
    output_format: OutputFormat
    usage: FrozenDict = field(default_factory=FrozenDict)
    _json_data: FrozenDict | None = field(default=None, repr=False)

    @property
    def data(self) -> dict[str, object] | None:
        """Return parsed JSON in json mode, or None in text mode."""

        if self._json_data is None:
            return None
        return cast(dict[str, object], to_json_compatible(self._json_data))

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable result; this includes the model's output text."""

        return {
            "text": self.text,
            "data": self.data,
            "model": self.model,
            "response_id": self.response_id,
            "finish_reason": self.finish_reason,
            "output_format": self.output_format,
            "usage": to_json_compatible(self.usage),
        }
