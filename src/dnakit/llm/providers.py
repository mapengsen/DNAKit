"""Connection presets for general-purpose LLM calls; models remain user-selected."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class LLMProvider:
    """Compatible API defaults, not a guarantee about every hosted model's features."""

    name: str
    base_url: str | None
    api_key_env: str
    token_limit_parameter: Literal["max_tokens", "max_completion_tokens"] = "max_tokens"


_PROVIDERS = (
    LLMProvider("deepseek", "https://api.deepseek.com", "DEEPSEEK_API_KEY"),
    LLMProvider(
        "gemini", "https://generativelanguage.googleapis.com/v1beta/openai", "GEMINI_API_KEY"
    ),
    # Beijing's legacy endpoint is still supported. Users may override base_url
    # with their region/workspace-specific Model Studio endpoint.
    LLMProvider("qwen", "https://dashscope.aliyuncs.com/compatible-mode/v1", "DASHSCOPE_API_KEY"),
    LLMProvider("kimi", "https://api.moonshot.cn/v1", "MOONSHOT_API_KEY"),
    LLMProvider("openrouter", "https://openrouter.ai/api/v1", "OPENROUTER_API_KEY"),
    LLMProvider("openai", "https://api.openai.com/v1", "OPENAI_API_KEY", "max_completion_tokens"),
    LLMProvider("custom", None, "DNAKIT_LLM_API_KEY"),
)
_ALIASES = {
    "google": "gemini",
    "google_gemini": "gemini",
    "dashscope": "qwen",
    "alibaba": "qwen",
    "moonshot": "kimi",
    "openai_compatible": "custom",
}


def available_llm_providers() -> tuple[LLMProvider, ...]:
    """List presets locally, without reading API keys or making requests."""

    return _PROVIDERS


def get_llm_provider(name: str) -> LLMProvider:
    """Resolve a provider name/alias, rejecting unknown names before a request."""

    if isinstance(name, str):
        normalized = name.strip().lower().replace("-", "_")
        canonical = _ALIASES.get(normalized, normalized)
        for provider in _PROVIDERS:
            if provider.name == canonical:
                return provider
    raise ConfigurationError(
        "Unknown LLM provider.",
        code="INVALID_LLM_CONFIG",
        hint="Choose deepseek, gemini, qwen, kimi, openrouter, openai, or custom.",
    )
