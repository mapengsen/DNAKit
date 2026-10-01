"""Search, describe, and execute the Agent-facing DNAKit tool catalog."""

from __future__ import annotations

import importlib
import inspect
import json
import re
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError

from .adapters import coerce_arguments
from .catalog import CATEGORY_LABELS, AgentToolSpec, _build_spec, build_tool_catalog
from .embeddings import EmbeddingConfig, ToolSearchMode, _ToolEmbeddingIndex
from .reranking import RerankerConfig, _rerank
from .retrieval import BM25Config, _bm25_scores, _rrf
from .serialization import AgentSerializationLimits, serialize_for_agent

_TOKEN = re.compile(r"[\w]+", re.UNICODE)
_MAX_ARGUMENT_BYTES = 2_000_000


class AgentToolRegistry:
    """Immutable searchable registry over public DNAKit functions."""

    def __init__(
        self,
        tools: tuple[AgentToolSpec, ...] | None = None,
        *,
        serialization_limits: AgentSerializationLimits | None = None,
        embedding_config: EmbeddingConfig | None = None,
        bm25_config: BM25Config | None = None,
        reranker_config: RerankerConfig | None = None,
    ) -> None:
        if embedding_config is not None and not isinstance(embedding_config, EmbeddingConfig):
            raise ConfigurationError(
                "embedding_config must be EmbeddingConfig or None.",
                code="INVALID_EMBEDDING_CONFIG",
            )
        self._embedding_config = embedding_config
        self._embedding_index = (
            _ToolEmbeddingIndex(embedding_config) if embedding_config is not None else None
        )
        if bm25_config is not None and not isinstance(bm25_config, BM25Config):
            raise ConfigurationError(
                "bm25_config must be BM25Config or None.", code="INVALID_BM25_CONFIG"
            )
        if reranker_config is not None and not isinstance(reranker_config, RerankerConfig):
            raise ConfigurationError(
                "reranker_config must be RerankerConfig or None.", code="INVALID_RERANKER_CONFIG"
            )
        self._bm25_config = BM25Config() if bm25_config is None else bm25_config
        self._reranker_config = reranker_config
        discovered = build_tool_catalog() if tools is None else tools
        names = [tool.name for tool in discovered]
        if len(names) != len(set(names)):
            duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
            raise ConfigurationError(
                "Agent tool names must be unique.",
                code="AGENT_TOOL_DUPLICATE",
                context={"tools": duplicates},
            )
        self._tools = tuple(sorted(discovered, key=lambda tool: tool.name))
        self._by_name = {tool.name: tool for tool in self._tools}
        aliases: defaultdict[str, list[AgentToolSpec]] = defaultdict(list)
        for tool in self._tools:
            aliases[tool.name.removeprefix("dnakit.")].append(tool)
            aliases[tool.name.rsplit(".", 1)[-1]].append(tool)
        self._aliases = {name: tuple(values) for name, values in aliases.items()}
        if serialization_limits is not None and not isinstance(
            serialization_limits, AgentSerializationLimits
        ):
            raise ConfigurationError(
                "serialization_limits must be AgentSerializationLimits or None.",
                code="AGENT_SERIALIZATION_LIMIT_INVALID",
            )
        self._serialization_limits = (
            AgentSerializationLimits() if serialization_limits is None else serialization_limits
        )

    @property
    def tools(self) -> tuple[AgentToolSpec, ...]:
        return self._tools

    @property
    def embedding_config(self) -> EmbeddingConfig | None:
        return self._embedding_config

    @property
    def bm25_config(self) -> BM25Config:
        return self._bm25_config

    @property
    def reranker_config(self) -> RerankerConfig | None:
        return self._reranker_config

    def _refreshed(
        self, additions: tuple[str, ...] = (), *, replacements: tuple[str, ...] = ()
    ) -> AgentToolRegistry:
        """Refresh descriptions and generated exports without broadening a caller's filter."""
        from .descriptions import apply_native_description

        current = {tool.name: tool for tool in build_tool_catalog(verbose=False)}
        if any(name not in self._by_name for name in replacements):
            raise ConfigurationError(
                "Only a selected registered tool may be explicitly replaced after repair.",
                code="AGENT_DEVELOPMENT_INVALID",
            )
        for name in replacements:
            if name in current:
                continue
            old_tool = self._by_name[name]
            try:
                module_name, export_name = name.rsplit(".", 1)
                if not module_name.startswith("dnakit."):
                    raise ValueError("The custom repair target is not a DNAKit module.")
                module = importlib.import_module(module_name)
                exported = getattr(module, export_name)
                old = inspect.unwrap(old_tool.function)
                updated = inspect.unwrap(exported)
                if not inspect.isfunction(old) or not inspect.isfunction(updated):
                    raise ValueError("Both old and repaired exports must be Python functions.")
                old_source, new_source = inspect.getsourcefile(old), inspect.getsourcefile(updated)
                if (
                    exported is old_tool.function
                    or not old_source
                    or not new_source
                    or Path(old_source).resolve() != Path(new_source).resolve()
                    or old.__module__ != updated.__module__
                    or old.__qualname__ != updated.__qualname__
                    or old.__globals__ is not updated.__globals__
                ):
                    raise ValueError("The custom export does not match the repaired definition.")
                current[name] = _build_spec(module_name, export_name, exported)
            except (AttributeError, ImportError, OSError, TypeError, ValueError) as exc:
                raise ConfigurationError(
                    "The selected custom tool was not replaced by a matching repaired export.",
                    code="AGENT_DEVELOPMENT_INVALID",
                    context={"tool_name": name},
                ) from exc
        selected: dict[str, AgentToolSpec] = {}
        for tool in self.tools:
            fresh = current.get(tool.name)
            if fresh is not None and (tool.name in replacements or fresh.function is tool.function):
                selected[tool.name] = replace(
                    fresh,
                    effect=tool.effect,
                    open_world=tool.open_world,
                    requires_confirmation=tool.requires_confirmation,
                    agent_compatible=tool.agent_compatible and fresh.agent_compatible,
                    incompatibilities=tuple(
                        dict.fromkeys((*tool.incompatibilities, *fresh.incompatibilities))
                    ),
                )
            elif tool.category == "generated":
                if fresh is not None:
                    selected[tool.name] = fresh
            else:
                selected[tool.name] = apply_native_description(tool)
        for name in additions:
            if name not in current or current[name].category != "generated":
                raise ConfigurationError(
                    "Only a validated generated tool may be added at runtime.",
                    code="AGENT_DEVELOPMENT_INVALID",
                )
            selected[name] = current[name]
        return AgentToolRegistry(
            tuple(selected.values()),
            serialization_limits=self._serialization_limits,
            embedding_config=self._embedding_config,
            bm25_config=self._bm25_config,
            reranker_config=self._reranker_config,
        )

    def _validate_search_mode(self, mode: ToolSearchMode) -> None:
        if mode not in ("keyword", "semantic", "hybrid", "bm25", "bm25_hybrid", "rerank"):
            raise ConfigurationError(
                "mode must be keyword, semantic, hybrid, bm25, bm25_hybrid, or rerank.",
                code="AGENT_SEARCH_MODE_INVALID",
            )
        if mode not in ("keyword", "bm25") and self._embedding_index is None:
            raise ConfigurationError(
                "Semantic, hybrid, bm25_hybrid and rerank require an embedding configuration.",
                code="MISSING_EMBEDDING_CONFIG",
                hint="Use AgentToolRegistry(embedding_config=EmbeddingConfig(...)).",
            )
        if mode == "rerank" and self._reranker_config is None:
            raise ConfigurationError(
                "Rerank search requires a reranker configuration.",
                code="MISSING_RERANKER_CONFIG",
                hint="Use AgentToolRegistry(reranker_config=RerankerConfig(...), ...).",
            )

    def _validate_search_limit(self, limit: int, mode: ToolSearchMode) -> None:
        if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
            raise ConfigurationError(
                "Search limit must be an integer in [1, 100].", code="AGENT_TOOL_LIMIT_INVALID"
            )
        if mode == "rerank":
            assert self._reranker_config is not None
            if limit > self._reranker_config.candidate_limit:
                raise ConfigurationError(
                    "Search limit must not exceed RerankerConfig.candidate_limit.",
                    code="AGENT_TOOL_LIMIT_INVALID",
                )

    def _validate_search_credentials(self, mode: ToolSearchMode) -> None:
        if mode == "rerank":
            assert self._reranker_config is not None
            self._reranker_config.resolved_api_key()
        if mode not in ("keyword", "bm25"):
            assert self._embedding_config is not None
            self._embedding_config.resolved_api_key()

    def stats(self) -> dict[str, object]:
        """Return aggregate catalog counts without loading optional backends."""

        compatible = sum(tool.agent_compatible for tool in self._tools)
        return {
            "total": len(self._tools),
            "agent_compatible": compatible,
            "python_only": len(self._tools) - compatible,
            "categories": len({tool.category for tool in self._tools}),
        }

    def categories(self) -> list[dict[str, object]]:
        """Return stable category names and callable counts."""

        totals = Counter(tool.category for tool in self._tools)
        compatible = Counter(tool.category for tool in self._tools if tool.agent_compatible)
        result: list[dict[str, object]] = []
        for category in sorted(totals):
            labels = CATEGORY_LABELS.get(category, (category, category))
            result.append(
                {
                    "category": category,
                    "label_zh": labels[0],
                    "label_en": labels[1],
                    "total": totals[category],
                    "agent_compatible": compatible[category],
                }
            )
        return result

    def resolve(self, name: str) -> AgentToolSpec:
        """Resolve a full name or an unambiguous category/function alias."""

        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError(
                "Agent tool name must be a non-empty string.",
                code="AGENT_TOOL_NAME_INVALID",
            )
        normalized = name.strip()
        exact = self._by_name.get(normalized)
        if exact is not None:
            return exact
        matches = self._aliases.get(normalized, ())
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ConfigurationError(
                "Agent tool alias is ambiguous.",
                code="AGENT_TOOL_AMBIGUOUS",
                context={"alias": normalized, "matches": [tool.name for tool in matches]},
                hint="Use the complete dnakit.<category>.<function> name.",
            )
        suggestions = [
            tool.name for tool in self.search(normalized, limit=5, include_python_only=True)
        ]
        raise ConfigurationError(
            "Agent tool was not found.",
            code="AGENT_TOOL_NOT_FOUND",
            context={"tool": normalized, "suggestions": suggestions},
        )

    def describe(self, name: str) -> dict[str, object]:
        """Return one full Agent tool manifest."""

        return self.resolve(name).to_dict()

    def search(
        self,
        query: str = "",
        *,
        category: str | None = None,
        limit: int = 20,
        include_python_only: bool = False,
        mode: ToolSearchMode = "keyword",
    ) -> tuple[AgentToolSpec, ...]:
        """Return candidates without executing them; keyword search is the default.

        Semantic search uses the configured cloud embeddings and local cosine
        similarity. Hybrid search combines lexical and semantic rankings with
        reciprocal rank fusion (k=60), prioritizing exact function names. Filters
        apply before either search. Blank queries list tools without API calls.
        BM25 is local; bm25_hybrid fuses BM25 and semantic ranks without name boosts.
        Rerank sends the top candidate_limit bm25_hybrid results to a configured
        reranker, then truncates to limit (which must not exceed candidate_limit).
        Cloud failures propagate; the requested mode never silently changes.
        """

        self._validate_search_mode(mode)
        if not isinstance(query, str):
            raise ConfigurationError(
                "Agent tool search query must be a string.",
                code="AGENT_TOOL_QUERY_INVALID",
            )
        self._validate_search_limit(limit, mode)
        if category is not None and (
            not isinstance(category, str) or category not in CATEGORY_LABELS
        ):
            raise ConfigurationError(
                "Unknown Agent tool category.",
                code="AGENT_TOOL_CATEGORY_INVALID",
                context={"category": category},
            )
        normalized = query.strip().casefold()
        query_tokens = set(_TOKEN.findall(normalized))
        ranked: list[tuple[int, str, AgentToolSpec]] = []
        documents: list[tuple[str, str]] = []
        for tool in self._tools:
            if category is not None and tool.category != category:
                continue
            if not include_python_only and not tool.agent_compatible:
                continue
            labels = CATEGORY_LABELS.get(tool.category, (tool.category, tool.category))
            retrieval_text = tool.search_text if tool.search_text is not None else tool.description
            documents.append((tool.name, "\n".join((tool.name, retrieval_text, *labels))))
            short_name = tool.name.rsplit(".", 1)[-1]
            haystack = " ".join(
                (
                    tool.name,
                    short_name.replace("_", " "),
                    retrieval_text,
                    labels[0],
                    labels[1],
                )
            ).casefold()
            if not normalized:
                score = 1
            else:
                haystack_tokens = set(_TOKEN.findall(haystack))
                score = 0
                if normalized == tool.name.casefold():
                    score += 1_000
                if normalized in {short_name.casefold(), short_name.replace("_", " ").casefold()}:
                    score += 800
                if normalized in haystack:
                    score += 300
                score += 40 * len(query_tokens & haystack_tokens)
                score += sum(10 for token in query_tokens if token in haystack)
            if score > 0:
                ranked.append((score, tool.name, tool))
        ranked.sort(key=lambda item: (-item[0], item[1]))
        if mode == "keyword" or not normalized or not documents:
            return tuple(item[2] for item in ranked[:limit])
        bm25_names: list[str] = []
        if mode in ("bm25", "bm25_hybrid", "rerank"):
            bm25_scores = _bm25_scores(query, documents, self._bm25_config)
            bm25_names = sorted(bm25_scores, key=lambda name: (-bm25_scores[name], name))
            if mode == "bm25":
                return tuple(self._by_name[name] for name in bm25_names[:limit])
        if mode == "rerank":
            # Fail before paying for embeddings if the reranker credential is missing.
            self._validate_search_credentials(mode)
        assert self._embedding_index is not None
        semantic_names = self._embedding_index.rank(query.strip(), documents)
        if mode == "semantic":
            return tuple(self._by_name[name] for name in semantic_names[:limit])
        lexical_names = [item[1] for item in ranked] if mode == "hybrid" else bm25_names
        names = _rrf((lexical_names, semantic_names))
        if mode == "hybrid":
            # Stable sort preserves the original hybrid's exact-name precedence.
            names.sort(
                key=lambda name: (
                    normalized not in (name.casefold(), name.rsplit(".", 1)[-1].casefold())
                )
            )
        if mode == "rerank":
            assert self._reranker_config is not None
            texts = dict(documents)
            candidates = names[: self._reranker_config.candidate_limit]
            names = _rerank(
                query.strip(), [(name, texts[name]) for name in candidates], self._reranker_config
            )
        return tuple(self._by_name[name] for name in names[:limit])

    def execute(
        self,
        name: str,
        arguments: Mapping[str, object] | None = None,
        *,
        allow_side_effects: bool = False,
    ) -> object:
        """Execute one compatible public function and return bounded JSON."""

        tool = self.resolve(name)
        if not tool.agent_compatible:
            raise ConfigurationError(
                "This DNAKit function requires live Python objects and is not Agent-callable.",
                code="AGENT_TOOL_PYTHON_ONLY",
                context={
                    "tool": tool.name,
                    "incompatibilities": list(tool.incompatibilities),
                },
                hint="Call this function through the regular Python API.",
            )
        if tool.requires_confirmation and not allow_side_effects:
            raise ConfigurationError(
                "This DNAKit tool may write files, download models, or run an external program.",
                code="AGENT_TOOL_CONFIRMATION_REQUIRED",
                context={"tool": tool.name, "effect": tool.effect},
                hint="Review the tool and set allow_side_effects=true to authorize this call.",
            )
        resolved_arguments: Mapping[str, object] = {} if arguments is None else arguments
        if not isinstance(resolved_arguments, Mapping) or any(
            not isinstance(key, str) for key in resolved_arguments
        ):
            raise ConfigurationError(
                "Agent tool arguments must be an object with string keys.",
                code="AGENT_ARGUMENTS_INVALID",
            )
        try:
            encoded_arguments = json.dumps(
                resolved_arguments,
                ensure_ascii=False,
                separators=(",", ":"),
            ).encode("utf-8")
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Agent tool arguments must contain only JSON values.",
                code="AGENT_ARGUMENTS_INVALID",
            ) from exc
        if len(encoded_arguments) > _MAX_ARGUMENT_BYTES:
            raise ConfigurationError(
                "Agent tool arguments exceed the byte limit.",
                code="AGENT_ARGUMENT_BYTE_LIMIT",
                context={
                    "argument_bytes": len(encoded_arguments),
                    "max_bytes": _MAX_ARGUMENT_BYTES,
                },
            )
        positional, keywords = coerce_arguments(tool.call_plan, resolved_arguments)
        result = tool.function(*positional, **keywords)
        return serialize_for_agent(result, limits=self._serialization_limits)


@lru_cache(maxsize=1)
def default_tool_registry(*, verbose: bool = True) -> AgentToolRegistry:
    """Return the process-wide immutable DNAKit Agent registry."""

    return AgentToolRegistry()


__all__ = ["AgentToolRegistry", "default_tool_registry"]

_register_public_verbose(globals())
