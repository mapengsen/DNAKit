"""Provider presets for privacy-bounded, OpenAI-compatible file schema analysis."""

from __future__ import annotations

from dataclasses import dataclass

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class SchemaLLMProvider:
    """Connection defaults, not a guarantee that every hosted model supports JSON."""

    name: str
    base_url: str | None
    api_key_env: str


_PROVIDERS = (
    SchemaLLMProvider("deepseek", "https://api.deepseek.com", "DEEPSEEK_API_KEY"),
    SchemaLLMProvider(
        "gemini", "https://generativelanguage.googleapis.com/v1beta/openai", "GEMINI_API_KEY"
    ),
    SchemaLLMProvider(
        "qwen", "https://dashscope.aliyuncs.com/compatible-mode/v1", "DASHSCOPE_API_KEY"
    ),
    SchemaLLMProvider("kimi", "https://api.moonshot.cn/v1", "MOONSHOT_API_KEY"),
    SchemaLLMProvider("openrouter", "https://openrouter.ai/api/v1", "OPENROUTER_API_KEY"),
    SchemaLLMProvider("openai", "https://api.openai.com/v1", "OPENAI_API_KEY"),
    SchemaLLMProvider("custom", None, "LLM_API_KEY"),
)
_ALIASES = {
    "google": "gemini",
    "google_gemini": "gemini",
    "dashscope": "qwen",
    "alibaba": "qwen",
    "moonshot": "kimi",
    "openai_compatible": "custom",
}


def available_schema_llm_providers(*, verbose: bool = True) -> tuple[SchemaLLMProvider, ...]:
    """List connection presets without reading keys or making network requests."""

    return _PROVIDERS


def get_schema_llm_provider(name: str, *, verbose: bool = True) -> SchemaLLMProvider:
    """Resolve a provider name or alias; model names remain caller-controlled."""

    if isinstance(name, str):
        normalized = name.strip().lower().replace("-", "_")
        canonical = _ALIASES.get(normalized, normalized)
        for provider in _PROVIDERS:
            if provider.name == canonical:
                return provider
    raise ConfigurationError(
        "Unknown schema LLM provider.",
        code="INVALID_LLM_SCHEMA_CONFIG",
        hint="Choose deepseek, gemini, qwen, kimi, openrouter, openai, or custom.",
    )


__all__ = ["SchemaLLMProvider", "available_schema_llm_providers", "get_schema_llm_provider"]

_register_public_verbose(globals())
