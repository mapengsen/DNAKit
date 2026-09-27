"""Bounded local-document retrieval with atomic JSON indexes and optional embeddings.

UTF-8 TXT, Markdown, RST, JSON/JSONL, CSV/TSV and YAML files are treated as
plain text. PDF, Office documents, images and compressed files must first be
converted to supported text. Retrieved text is evidence, never executable code
or trusted instructions. This module neither executes documents nor calls a
chat model; the caller supplies retrieved passages to its answering model.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._progress import Progress, track
from dnakit.exceptions import ConfigurationError, InputFormatError

from .embeddings import EmbeddingConfig, ToolSearchMode, _ToolEmbeddingIndex

_SUFFIXES = frozenset(
    {".txt", ".md", ".markdown", ".rst", ".json", ".jsonl", ".csv", ".tsv", ".yaml", ".yml"}
)
_TOKEN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]+|[a-z0-9]+(?:[_-][a-z0-9]+)*", re.I)
_CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
_INDEX_VERSION = 1
_MAX_SCAN_ENTRIES = 100_000


def _path(value: object, name: str) -> Path:
    try:
        if not isinstance(value, (str, os.PathLike)):
            raise TypeError
        raw = os.fspath(value)
        if not isinstance(raw, str) or not raw or "\0" in raw:
            raise ValueError
        return Path(raw).expanduser().resolve()
    except (TypeError, ValueError, OSError, RuntimeError):
        raise ConfigurationError(
            f"{name} must be a valid filesystem path.", code="INVALID_RAG_CONFIG"
        ) from None


@dataclass(frozen=True, slots=True)
class RAGConfig:
    """Configure a user-selected local document directory, without a database.

    Keyword mode stays entirely local. Semantic/hybrid modes require an explicit
    EmbeddingConfig and send document chunks and the question to that embedding
    endpoint; chat credentials are never reused. EmbeddingConfig.cache_dir
    controls vector persistence independently of this module's JSON index.

    Hidden files/directories, symbolic links and unsupported formats are skipped.
    Every retrieval fingerprints current source bytes, so additions, deletions
    and edits invalidate the index even when timestamps are preserved. Limits
    raise explicit errors instead of silently indexing a partial collection.
    """

    directory: str | os.PathLike[str]
    index_dir: str | os.PathLike[str] | None = None
    mode: ToolSearchMode = "keyword"
    embedding_config: EmbeddingConfig | None = None
    top_k: int = 5
    chunk_size: int = 1200
    chunk_overlap: int = 150
    max_context_chars: int = 12_000
    max_files: int = 2000
    max_file_bytes: int = 8_000_000
    max_total_bytes: int = 64_000_000
    max_chunks: int = 20_000

    def __post_init__(self) -> None:
        root = _path(self.directory, "directory")
        object.__setattr__(self, "directory", root)
        if self.index_dir is None:
            cache = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
            index_dir = cache / "dnakit/knowledge" / hashlib.sha256(str(root).encode()).hexdigest()
        else:
            index_dir = _path(self.index_dir, "index_dir")
        index_dir = index_dir.resolve()
        if index_dir == root:
            raise ConfigurationError(
                "index_dir must differ from directory; use a dedicated index subdirectory.",
                code="INVALID_RAG_CONFIG",
            )
        object.__setattr__(self, "index_dir", index_dir)
        if self.mode not in ("keyword", "semantic", "hybrid"):
            raise ConfigurationError(
                "RAG mode must be keyword, semantic, or hybrid.", code="INVALID_RAG_CONFIG"
            )
        if self.embedding_config is not None and not isinstance(
            self.embedding_config, EmbeddingConfig
        ):
            raise ConfigurationError(
                "embedding_config must be EmbeddingConfig or None.", code="INVALID_RAG_CONFIG"
            )
        if (
            self.embedding_config is not None
            and self.embedding_config.cache_dir is not None
            and Path(self.embedding_config.cache_dir).expanduser().resolve() == root
        ):
            raise ConfigurationError(
                "The embedding cache must use a directory separate from the knowledge source root.",
                code="INVALID_RAG_CONFIG",
            )
        if self.mode != "keyword" and self.embedding_config is None:
            raise ConfigurationError(
                "Semantic/hybrid RAG requires an explicit embedding_config.",
                code="MISSING_EMBEDDING_CONFIG",
            )
        for name, minimum, maximum in (
            ("top_k", 1, 100),
            ("chunk_size", 64, 100_000),
            ("chunk_overlap", 0, 99_999),
            ("max_context_chars", 1, 1_000_000),
            ("max_files", 1, 100_000),
            ("max_file_bytes", 1, 100_000_000),
            ("max_total_bytes", 1, 1_000_000_000),
            ("max_chunks", 1, 100_000),
        ):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or not minimum <= value <= maximum
            ):
                raise ConfigurationError(
                    f"{name} must be an integer in [{minimum}, {maximum}].",
                    code="INVALID_RAG_CONFIG",
                )
        if self.chunk_overlap >= self.chunk_size:
            raise ConfigurationError(
                "chunk_overlap must be smaller than chunk_size.", code="INVALID_RAG_CONFIG"
            )


def _files(config: RAGConfig) -> Iterator[Path]:
    root = Path(config.directory)
    excluded = {Path(config.index_dir)} if config.index_dir is not None else set()
    if config.embedding_config is not None and config.embedding_config.cache_dir is not None:
        excluded.add(Path(config.embedding_config.cache_dir).expanduser().resolve())
    scanned = 0
    supported = 0

    def fail(error: OSError) -> None:
        raise InputFormatError(
            "The knowledge directory could not be scanned.", code="RAG_SOURCE_READ_ERROR"
        ) from error

    for folder, directories, filenames in os.walk(root, followlinks=False, onerror=fail):
        parent = Path(folder)
        scanned += len(directories) + len(filenames)
        if scanned > _MAX_SCAN_ENTRIES:
            raise InputFormatError(
                "The knowledge directory exceeds the 100000-entry scan limit.",
                code="RAG_SOURCE_LIMIT",
            )
        directories[:] = sorted(
            name
            for name in directories
            if not name.startswith(".")
            and not (parent / name).is_symlink()
            and (parent / name) not in excluded
        )
        for name in sorted(filenames):
            path = parent / name
            if name.startswith(".") or path.suffix.lower() not in _SUFFIXES or path.is_symlink():
                continue
            if not path.is_file():
                continue
            supported += 1
            if supported > config.max_files:
                raise InputFormatError(
                    "The knowledge directory exceeds max_files.", code="RAG_SOURCE_LIMIT"
                )
            yield path


def _read_sources(config: RAGConfig, verbose: bool) -> tuple[dict[str, str], dict[str, str]]:
    root = Path(config.directory)
    if not root.is_dir():
        raise ConfigurationError(
            "The configured knowledge directory does not exist or is not a directory.",
            code="RAG_DIRECTORY_NOT_FOUND",
        )
    manifest: dict[str, str] = {}
    documents: dict[str, str] = {}
    total = 0
    for path in track(_files(config), description="Reading knowledge files", disable=not verbose):
        try:
            with path.open("rb") as stream:
                raw = stream.read(min(config.max_file_bytes, config.max_total_bytes - total) + 1)
        except OSError:
            raise InputFormatError(
                "A knowledge file could not be read.",
                code="RAG_SOURCE_READ_ERROR",
                context={"path": str(path)},
            ) from None
        total += len(raw)
        if len(raw) > config.max_file_bytes or total > config.max_total_bytes:
            raise InputFormatError(
                "Knowledge files exceed max_file_bytes or max_total_bytes.",
                code="RAG_SOURCE_LIMIT",
                context={"path": str(path)},
            )
        try:
            content = raw.decode("utf-8-sig")
            if "\0" in content:
                raise ValueError
        except (UnicodeError, ValueError):
            raise InputFormatError(
                "Knowledge files must contain UTF-8 text without NUL characters.",
                code="RAG_SOURCE_ENCODING",
                context={"path": str(path)},
            ) from None
        source = path.relative_to(root).as_posix()
        manifest[source] = hashlib.sha256(raw).hexdigest()
        documents[source] = content
    return manifest, documents


def _split(text: str, config: RAGConfig) -> Iterator[tuple[int, str]]:
    start = 0
    while start < len(text):
        end = min(start + config.chunk_size, len(text))
        if end < len(text):
            boundary = text.rfind("\n", start + config.chunk_size // 2, end)
            if boundary == -1:
                boundary = text.rfind(" ", start + config.chunk_size // 2, end)
            if boundary != -1:
                end = boundary + 1
        chunk = text[start:end]
        if chunk.strip():
            yield start, chunk
        if end == len(text):
            return
        # A large configured overlap must never prevent forward movement.
        start = max(start + 1, end - config.chunk_overlap)


def _digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, ensure_ascii=True, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def _chunks(config: RAGConfig, verbose: bool) -> list[dict[str, Any]]:
    manifest, documents = _read_sources(config, verbose)
    root = Path(config.directory)
    assert config.index_dir is not None
    index_dir = Path(config.index_dir)
    index_path = index_dir / f"index-{hashlib.sha256(str(root).encode()).hexdigest()}.json"
    identity = {
        "version": _INDEX_VERSION,
        "directory": str(root),
        "manifest": manifest,
        "chunk_size": config.chunk_size,
        "chunk_overlap": config.chunk_overlap,
    }
    # Account for escaped characters, chunk overlap, paths and per-chunk metadata.
    max_index_bytes = min(512_000_000, config.max_total_bytes * 12 + config.max_chunks * 4096)
    try:
        with index_path.open("rb") as stream:
            raw = stream.read(max_index_bytes + 1)
        if len(raw) <= max_index_bytes:
            entry = json.loads(raw)
            if isinstance(entry, dict) and entry.get("identity") == identity:
                chunks = entry.get("chunks")
                if (
                    isinstance(chunks, list)
                    and len(chunks) <= config.max_chunks
                    and entry.get("chunks_sha256") == _digest(chunks)
                    and all(
                        isinstance(chunk, dict)
                        and isinstance(chunk.get("text"), str)
                        and isinstance(chunk.get("chunk_id"), str)
                        and isinstance(chunk.get("source"), str)
                        and chunk["source"] in manifest
                        and chunk.get("path") == str(root / chunk["source"])
                        for chunk in chunks
                    )
                ):
                    return chunks
    except (OSError, UnicodeError, ValueError, TypeError, RecursionError):
        pass
    chunks = []
    for source, content in track(
        documents.items(), description="Indexing knowledge files", disable=not verbose
    ):
        for offset, text in _split(content, config):
            if len(chunks) >= config.max_chunks:
                raise InputFormatError(
                    "Knowledge text exceeds max_chunks; increase chunk_size or max_chunks.",
                    code="RAG_SOURCE_LIMIT",
                )
            chunks.append(
                {
                    "source": source,
                    "path": str(root / source),
                    "chunk_id": f"{source}:{offset}",
                    "offset": offset,
                    "text": text,
                }
            )
    temporary: str | None = None
    try:
        index_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=index_dir, delete=False
        ) as stream:
            temporary = stream.name
            json.dump(
                {"identity": identity, "chunks": chunks, "chunks_sha256": _digest(chunks)},
                stream,
                ensure_ascii=False,
                allow_nan=False,
            )
        os.replace(temporary, index_path)
    except OSError:
        raise ConfigurationError(
            "The local knowledge index could not be written.", code="RAG_INDEX_WRITE_ERROR"
        ) from None
    finally:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)
    return chunks


def _tokens(text: str, *, query: bool = False) -> list[str]:
    result: list[str] = []
    for word in _TOKEN.findall(text.casefold()):
        if _CJK.match(word):
            if not query or len(word) == 1:
                result.extend(word)
            result.extend(word[index : index + 2] for index in range(len(word) - 1))
        else:
            result.append(word)
            if "_" in word or "-" in word:
                result.extend(re.split(r"[_-]", word))
    return result


def _keyword_scores(query: str, chunks: list[dict[str, Any]]) -> dict[str, float]:
    terms = set(_tokens(query, query=True))
    if not terms or not chunks:
        return {}
    frequencies = [Counter(_tokens(f"{chunk['source']}\n{chunk['text']}")) for chunk in chunks]
    document_frequency: Counter[str] = Counter()
    for counts in frequencies:
        document_frequency.update(terms.intersection(counts))
    lengths = [sum(counts.values()) for counts in frequencies]
    average = sum(lengths) / len(lengths) or 1.0
    scores: dict[str, float] = {}
    for chunk, counts, length in zip(chunks, frequencies, lengths, strict=True):
        score = 0.0
        for term in terms.intersection(counts):
            count = counts[term]
            idf = math.log(
                1
                + (len(chunks) - document_frequency[term] + 0.5) / (document_frequency[term] + 0.5)
            )
            score += idf * count * 2.2 / (count + 1.2 * (0.25 + 0.75 * length / average))
        if score > 0:
            scores[chunk["chunk_id"]] = score
    return scores


def retrieve_knowledge(
    query: str, *, config: RAGConfig, verbose: bool = True
) -> list[dict[str, Any]]:
    """Return ranked source passages before the caller asks its chat model.

    Each result includes source (relative filename), path (absolute filename),
    chunk_id, offset, text and score. Keyword scores use BM25; semantic scores
    are reciprocal ranks and hybrid scores use reciprocal-rank fusion (k=60),
    not calibrated probabilities. Text is bounded by max_context_chars overall.
    Empty questions and collections return no passages; no model is called.
    """
    if not isinstance(config, RAGConfig):
        raise ConfigurationError("config must be RAGConfig.", code="INVALID_RAG_CONFIG")
    if not isinstance(query, str):
        raise ConfigurationError("RAG query must be a string.", code="INVALID_RAG_QUERY")
    if not query.strip():
        return []
    chunks = _chunks(config, verbose)
    if not chunks:
        return []
    lexical = _keyword_scores(query, chunks) if config.mode != "semantic" else {}
    scores = lexical
    if config.mode != "keyword":
        assert config.embedding_config is not None
        with Progress(disable=not verbose) as progress:
            task = progress.add_task("Retrieving knowledge embeddings", total=None)
            ranking = _ToolEmbeddingIndex(config.embedding_config).rank(
                query,
                [(chunk["chunk_id"], f"{chunk['source']}\n{chunk['text']}") for chunk in chunks],
            )
            progress.update(task, total=len(chunks), completed=len(chunks))
        scores = {key: 1.0 / (60 + rank) for rank, key in enumerate(ranking, 1)}
        if config.mode == "hybrid":
            lexical_ranking = sorted(lexical, key=lambda key: (-lexical[key], key))
            for rank, key in enumerate(lexical_ranking, 1):
                scores[key] = scores.get(key, 0.0) + 1.0 / (60 + rank)
    lookup = {chunk["chunk_id"]: chunk for chunk in chunks}
    remaining = config.max_context_chars
    result: list[dict[str, Any]] = []
    for key in sorted(scores, key=lambda key: (-scores[key], key))[: config.top_k]:
        if remaining <= 0:
            break
        chunk = dict(lookup[key])
        chunk["text"] = chunk["text"][:remaining]
        chunk["score"] = scores[key]
        remaining -= len(chunk["text"])
        result.append(chunk)
    return result
