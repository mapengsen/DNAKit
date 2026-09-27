"""Offline tests for bounded, file-backed local-document RAG retrieval."""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.tools.embeddings import EmbeddingConfig
from dnakit.tools.knowledge import RAGConfig, retrieve_knowledge


@pytest.fixture
def config(tmp_path: Path) -> RAGConfig:
    documents = tmp_path / "documents"
    documents.mkdir()
    return RAGConfig(directory=documents, index_dir=tmp_path / "index")


def _write(config: RAGConfig, source: str, text: str) -> Path:
    path = Path(config.directory) / source
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _index(config: RAGConfig) -> Path:
    assert config.index_dir is not None
    return next(Path(config.index_dir).glob("index-*.json"))


def test_keyword_retrieval_ranks_relevant_sources_locally(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write(config, "sequencing/qc.md", "Sequencing quality control removes low-quality reads.")
    _write(config, "plants.txt", "Photosynthesis converts light into chemical energy.")
    fail = Mock(side_effect=AssertionError("Keyword RAG must never call an embedding API"))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", fail)
    matches = retrieve_knowledge("sequencing quality", config=config, verbose=False)
    assert [match["source"] for match in matches] == ["sequencing/qc.md"]
    assert matches[0]["path"] == str(Path(config.directory) / "sequencing/qc.md")
    assert matches[0]["score"] > 0
    assert "low-quality" in matches[0]["text"]
    assert _index(config).is_file()
    fail.assert_not_called()


def test_chinese_keywords_use_substrings_without_whitespace(config: RAGConfig) -> None:
    _write(config, "quality.md", "测序质控需要检查碱基质量、接头污染与重复率。")
    _write(config, "biology.md", "植物光合作用产生氧气。")
    matches = retrieve_knowledge("如何进行测序质控?", config=config, verbose=False)
    assert matches[0]["source"] == "quality.md"
    assert retrieve_knowledge("污染", config=config, verbose=False)[0]["source"] == "quality.md"


def test_content_changes_invalidate_index_even_with_same_size_and_mtime(config: RAGConfig) -> None:
    source = _write(config, "sample.txt", "alpha marker")
    assert retrieve_knowledge("alpha", config=config, verbose=False)
    before = source.stat()
    index_before = _index(config).read_text()
    source.write_text("omega marker", encoding="utf-8")
    os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert retrieve_knowledge("alpha", config=config, verbose=False) == []
    assert retrieve_knowledge("omega", config=config, verbose=False)[0]["text"] == "omega marker"
    assert _index(config).read_text() != index_before


def test_added_deleted_sources_and_unchanged_indexes(config: RAGConfig) -> None:
    first = _write(config, "first.txt", "variant evidence")
    retrieve_knowledge("variant", config=config, verbose=False)
    initial_mtime = _index(config).stat().st_mtime_ns
    retrieve_knowledge("variant", config=config, verbose=False)
    assert _index(config).stat().st_mtime_ns == initial_mtime
    second = _write(config, "second.txt", "variant annotation")
    first.unlink()
    matches = retrieve_knowledge("variant", config=config, verbose=False)
    assert [match["source"] for match in matches] == [second.name]


def test_internal_indexes_hidden_files_symlinks_and_unsupported_formats_are_excluded(
    config: RAGConfig, tmp_path: Path
) -> None:
    config = replace(config, index_dir=Path(config.directory) / "generated")
    _write(config, "guide.txt", "variant annotation")
    _write(config, ".hidden.txt", "variant hidden")
    _write(config, ".private/secret.md", "variant private")
    _write(config, "paper.pdf", "variant PDF unsupported")
    outside = tmp_path / "outside.txt"
    outside.write_text("variant outside", encoding="utf-8")
    (Path(config.directory) / "symlink.txt").symlink_to(outside)
    (Path(config.directory) / "symlink_dir").symlink_to(tmp_path, target_is_directory=True)
    for _ in range(2):
        matches = retrieve_knowledge("variant", config=config, verbose=False)
        assert [match["source"] for match in matches] == ["guide.txt"]
    entry = json.loads(_index(config).read_text())
    assert list(entry["identity"]["manifest"]) == ["guide.txt"]


def test_corrupt_or_modified_index_is_rebuilt(config: RAGConfig) -> None:
    _write(config, "guide.txt", "variant annotation")
    retrieve_knowledge("variant", config=config, verbose=False)
    path = _index(config)
    path.write_text("{broken", encoding="utf-8")
    assert retrieve_knowledge("variant", config=config, verbose=False)[0]["text"] == (
        "variant annotation"
    )
    entry = json.loads(path.read_text())
    entry["chunks"][0]["text"] = "altered evidence"
    path.write_text(json.dumps(entry), encoding="utf-8")
    assert retrieve_knowledge("variant", config=config, verbose=False)[0]["text"] == (
        "variant annotation"
    )


def test_chunk_overlap_and_context_budget(config: RAGConfig) -> None:
    content = "variant DNA annotation evidence " * 20
    _write(config, "guide.txt", content)
    config = replace(config, chunk_size=80, chunk_overlap=20, top_k=5, max_context_chars=110)
    matches = retrieve_knowledge("variant", config=config, verbose=False)
    assert len(matches) == 2
    assert sum(len(match["text"]) for match in matches) == 110
    chunks = json.loads(_index(config).read_text())["chunks"]
    assert chunks[0]["text"][-20:] == chunks[1]["text"][:20]
    assert all(content[item["offset"] :].startswith(item["text"]) for item in chunks)


def test_retrieved_instructions_are_returned_as_text_and_never_executed(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = 'variant instructions: os.system("touch /tmp/unwanted-rag-command")'
    _write(config, "instructions.md", payload)
    execute = Mock(side_effect=AssertionError("Knowledge text must not execute"))
    monkeypatch.setattr("os.system", execute)
    assert retrieve_knowledge("variant", config=config, verbose=False)[0]["text"] == payload
    execute.assert_not_called()


def test_semantic_uses_explicit_embeddings_and_reuses_persisted_vectors(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _write(config, "blood.md", "Hematology examines erythrocytes.")
    _write(config, "plants.md", "Botany examines flowers.")
    embedding = EmbeddingConfig(
        model="offline-model", base_url="https://example.invalid/v1", cache_dir=tmp_path / "vectors"
    )
    config = replace(config, mode="semantic", embedding_config=embedding, top_k=1)
    calls: list[list[str]] = []

    def fake_embed(texts: list[str], config: EmbeddingConfig) -> list[tuple[float, ...]]:
        assert config is embedding
        calls.append(texts)
        return [(0.0, 1.0) if "flowers" in text else (1.0, 0.0) for text in texts]

    monkeypatch.setattr("dnakit.tools.embeddings._embed", fake_embed)
    matches = retrieve_knowledge("red cells", config=config, verbose=False)
    assert matches[0]["source"] == "blood.md"
    assert [len(batch) for batch in calls] == [1, 2]
    calls.clear()
    retrieve_knowledge("red cells", config=config, verbose=False)
    assert calls == [["red cells"]]


def test_hybrid_combines_lexical_and_semantic_evidence(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write(config, "blood.md", "Hematology examines erythrocytes.")
    _write(config, "plants.md", "Botany examines flowers.")
    embedding = EmbeddingConfig(
        model="offline", base_url="https://example.invalid/v1", cache_dir=None
    )
    monkeypatch.setattr(
        "dnakit.tools.embeddings._embed",
        lambda texts, config: [(0.0, 1.0) if "flowers" in text else (1.0, 0.0) for text in texts],
    )
    semantic = replace(config, mode="semantic", embedding_config=embedding)
    assert retrieve_knowledge("Botany", config=semantic, verbose=False)[0]["source"] == "blood.md"
    hybrid = replace(semantic, mode="hybrid")
    assert retrieve_knowledge("Botany", config=hybrid, verbose=False)[0]["source"] == "plants.md"


@pytest.mark.parametrize(
    "change",
    [
        {"directory": ""},
        {"directory": 42},
        {"index_dir": "\0"},
        {"mode": "invalid"},
        {"embedding_config": {}},
        {"top_k": True},
        {"top_k": 0},
        {"chunk_size": 63},
        {"chunk_overlap": 1200},
        {"chunk_overlap": -1},
        {"max_context_chars": 0},
        {"max_files": 0},
        {"max_file_bytes": False},
        {"max_total_bytes": 0},
        {"max_chunks": 0},
    ],
)
def test_invalid_configurations_fail_before_io(config: RAGConfig, change: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_RAG_CONFIG"):
        replace(config, **change)


@pytest.mark.parametrize("mode", ["semantic", "hybrid"])
def test_embedding_configuration_is_required_for_cloud_retrieval(
    config: RAGConfig, mode: Any
) -> None:
    with pytest.raises(ConfigurationError, match="MISSING_EMBEDDING_CONFIG"):
        replace(config, mode=mode)


def test_source_root_cannot_double_as_index_or_vector_cache(config: RAGConfig) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_RAG_CONFIG"):
        replace(config, index_dir=config.directory)
    embedding = EmbeddingConfig(
        model="offline", base_url="https://example.invalid/v1", cache_dir=config.directory
    )
    with pytest.raises(ConfigurationError, match="INVALID_RAG_CONFIG"):
        replace(config, embedding_config=embedding)


def test_internal_vector_cache_is_excluded(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write(config, "guide.txt", "variant annotation")
    embedding = EmbeddingConfig(
        model="offline",
        base_url="https://example.invalid/v1",
        cache_dir=Path(config.directory) / "vectors",
    )
    config = replace(config, mode="semantic", embedding_config=embedding)
    monkeypatch.setattr(
        "dnakit.tools.embeddings._embed", lambda texts, config: [(1.0,)] * len(texts)
    )
    for _ in range(2):
        matches = retrieve_knowledge("variant", config=config, verbose=False)
        assert [match["source"] for match in matches] == ["guide.txt"]


def test_directory_scan_is_bounded_even_for_unsupported_files(
    config: RAGConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write(config, "guide.txt", "variant annotation")
    _write(config, "unsupported.pdf", "PDF content")
    monkeypatch.setattr("dnakit.tools.knowledge._MAX_SCAN_ENTRIES", 1)
    with pytest.raises(InputFormatError, match="RAG_SOURCE_LIMIT"):
        retrieve_knowledge("variant", config=config, verbose=False)


@pytest.mark.parametrize(
    "change", [{"max_files": 1}, {"max_file_bytes": 5}, {"max_total_bytes": 15}, {"max_chunks": 1}]
)
def test_collection_limits_fail_without_silent_partial_index(
    config: RAGConfig, change: dict[str, Any]
) -> None:
    _write(config, "first.txt", "variant annotation")
    _write(config, "second.txt", "variant evidence")
    with pytest.raises(InputFormatError, match="RAG_SOURCE_LIMIT"):
        retrieve_knowledge("variant", config=replace(config, **change), verbose=False)
    assert config.index_dir is not None
    assert not Path(config.index_dir).exists()


def test_missing_directory_invalid_encoding_and_empty_inputs(config: RAGConfig) -> None:
    assert retrieve_knowledge("variant", config=config, verbose=False) == []
    assert retrieve_knowledge(" ", config=config, verbose=False) == []
    with pytest.raises(ConfigurationError, match="RAG_DIRECTORY_NOT_FOUND"):
        retrieve_knowledge(
            "variant",
            config=replace(config, directory=Path(config.directory) / "absent"),
            verbose=False,
        )
    (Path(config.directory) / "binary.txt").write_bytes(b"\xff\xfe")
    with pytest.raises(InputFormatError, match="RAG_SOURCE_ENCODING"):
        retrieve_knowledge("variant", config=config, verbose=False)


def test_progress_can_be_silenced(config: RAGConfig, capsys: pytest.CaptureFixture[str]) -> None:
    _write(config, "guide.txt", "variant evidence")
    retrieve_knowledge("variant", config=config, verbose=False)
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""
