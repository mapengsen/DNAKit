"""Offline integration checks for file checkpoints and retrieval before answering."""

from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, cast
from unittest.mock import Mock

import pytest

import dnakit.tools.agent as agent
from dnakit.core._json import FrozenDict
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.llm import LLMConfig, LLMResult
from dnakit.tools import (
    AgentResult,
    AgentToolRegistry,
    PythonExecutionConfig,
    RAGConfig,
    default_tool_registry,
    run_agent,
)

CONFIG = LLMConfig(provider="deepseek", model="offline-test", api_key="offline-placeholder")
REVERSE = "dnakit.ops.reverse_complement"
WEIGHT = "dnakit.thermodynamics.molecular_weight"


def _model(
    monkeypatch: pytest.MonkeyPatch,
    actions: Sequence[dict[str, Any] | BaseException],
) -> list[dict[str, Any]]:
    pending = iter(actions)
    inputs: list[dict[str, Any]] = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        inputs.append(deepcopy(kwargs["input_data"]))
        action = next(pending)
        if isinstance(action, BaseException):
            raise action
        return LLMResult(
            text="Offline action",
            model="offline-test",
            response_id=None,
            finish_reason="stop",
            output_format="json",
            _json_data=FrozenDict(action),
        )

    monkeypatch.setattr(agent, "call_llm", reply)
    return inputs


def _call(name: str, **arguments: Any) -> dict[str, Any]:
    return {"action": "call", "tool_name": name, "arguments": arguments}


def _discover(name: str) -> list[dict[str, Any]]:
    return [
        {"action": "search", "query": name},
        {"action": "describe", "tool_name": name},
    ]


def _finish(status: str = "completed") -> dict[str, Any]:
    return {"action": "finish", "status": status, "text": "Observed results checked"}


def _observation(result: AgentResult, index: int = -1) -> dict[str, Any]:
    return cast(dict[str, Any], result.steps[index]["observation"])


def _python() -> dict[str, Any]:
    return {
        "action": "execute_python",
        "reason": "Prepare a small FASTA file and record actual execution count",
        "code": (
            "from pathlib import Path\n"
            "counter = Path('counter.txt')\n"
            "counter.write_text(counter.read_text() + 'x' if counter.exists() else 'x')\n"
            "Path('prepared.fna').write_text('>sample\\nAACG\\n')\n"
        ),
        "inputs": {},
        "expected_outputs": ["prepared.fna"],
    }


@pytest.fixture
def registry() -> AgentToolRegistry:
    catalog = default_tool_registry(verbose=False)
    return AgentToolRegistry((catalog.resolve(REVERSE), catalog.resolve(WEIGHT)))


def test_three_completed_calls_resume_with_original_references_and_new_step_budget(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    execute = Mock(wraps=registry.execute)
    monkeypatch.setattr(registry, "execute", execute)
    _model(
        monkeypatch,
        [
            *_discover(REVERSE),
            _call(REVERSE, value="AACG"),
            _call(REVERSE, value="TTTT"),
            _call(REVERSE, value="GGCA"),
            KeyboardInterrupt(),
        ],
    )
    with pytest.raises(KeyboardInterrupt):
        run_agent(
            "Reverse three sequences and reverse the first back",
            config=CONFIG,
            registry=registry,
            cache_dir=tmp_path,
            verbose=False,
        )
    assert execute.call_count == 3

    inputs = _model(
        monkeypatch,
        [
            _call(
                REVERSE, value={"$ref": {"step": 3, "path": ["records", 0, "sequence", "parts", 0]}}
            ),
            _finish(),
        ],
    )
    result = run_agent(
        "Reverse three sequences and reverse the first back",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path,
        reuse_results=True,
        max_steps=2,
        verbose=False,
    )
    assert result.status == "completed" and result.reused_steps == 3
    assert execute.call_count == 4
    assert len(inputs[0]["history"]) == 5
    assert [entry["step"] for entry in result.steps] == list(range(1, 7))
    assert _observation(result)["result"]["records"][0]["sequence"]["parts"] == ["AACG"]
    assert result.cache_dir and (Path(result.cache_dir) / "checkpoint.json").is_file()


def test_step_limit_can_resume_without_rerunning_completed_call(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    execute = Mock(wraps=registry.execute)
    monkeypatch.setattr(registry, "execute", execute)
    action = _call(WEIGHT, sequence="ACGT")
    _model(monkeypatch, [*_discover(WEIGHT), action])
    with pytest.raises(BackendExecutionError, match="AGENT_STEP_LIMIT"):
        run_agent(
            "Calculate mass",
            config=CONFIG,
            registry=registry,
            cache_dir=tmp_path,
            max_steps=3,
            verbose=False,
        )
    _model(monkeypatch, [action, _finish()])
    result = run_agent(
        "Calculate mass",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path,
        reuse_results=True,
        max_steps=2,
        verbose=False,
    )
    assert result.reused_steps == 1
    assert _observation(result)["reused"] is True
    execute.assert_called_once()


def test_failed_call_is_retried_in_new_invocation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    real_execute = registry.execute
    attempts = 0

    def flaky(*args: Any, **kwargs: Any) -> object:
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError("Temporary input failure")
        return real_execute(*args, **kwargs)

    monkeypatch.setattr(registry, "execute", flaky)
    action = _call(WEIGHT, sequence="ACGT")
    _model(monkeypatch, [*_discover(WEIGHT), action, _finish("blocked")])
    first = run_agent(
        "Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False
    )
    assert _observation(first)["ok"] is False
    inputs = _model(monkeypatch, [action, _finish()])
    resumed = run_agent(
        "Calculate mass",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path,
        reuse_results=True,
        verbose=False,
    )
    assert resumed.status == "completed" and resumed.reused_steps == 0
    assert attempts == 2
    assert inputs[0]["history"][-1]["observation"]["ok"] is False
    assert _observation(resumed)["result"]["value_dalton"] == pytest.approx(1173.84)


@pytest.mark.parametrize("reuse,different_request", [(False, False), (True, True)])
def test_fresh_run_and_different_request_do_not_reuse_history(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
    reuse: bool,
    different_request: bool,
) -> None:
    execute = Mock(wraps=registry.execute)
    monkeypatch.setattr(registry, "execute", execute)
    actions = [*_discover(WEIGHT), _call(WEIGHT, sequence="ACGT"), _finish()]
    _model(monkeypatch, actions)
    first = run_agent(
        "Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False
    )
    inputs = _model(monkeypatch, actions)
    second = run_agent(
        "Calculate mass again" if different_request else "Calculate mass",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path,
        reuse_results=reuse,
        verbose=False,
    )
    assert second.reused_steps == 0 and inputs[0]["history"] == []
    assert execute.call_count == 2
    assert first.cache_dir != second.cache_dir
    assert len(list(tmp_path.rglob("checkpoint.json"))) == 2


@pytest.mark.parametrize("invalidate", ["input_changed", "output_deleted"])
def test_changed_input_or_deleted_output_truncates_only_affected_suffix(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
    invalidate: str,
) -> None:
    source = tmp_path / "input.fna"
    source.write_text(">sample\nAACG\n", encoding="utf-8")
    output = tmp_path / "reverse.fna"
    cache = tmp_path / "cache"
    execute = Mock(wraps=registry.execute)
    monkeypatch.setattr(registry, "execute", execute)
    file_action = _call(
        REVERSE, value=str(source), output_path=str(output), write_config={"overwrite": True}
    )
    last_action = _call(WEIGHT, sequence="TTTT")
    _model(
        monkeypatch,
        [
            *_discover(WEIGHT),
            _call(WEIGHT, sequence="ACGT"),
            *_discover(REVERSE),
            file_action,
            last_action,
            _finish(),
        ],
    )
    run_agent(
        "Calculate masses and reverse FASTA",
        config=CONFIG,
        registry=registry,
        cache_dir=cache,
        allow_side_effects=True,
        verbose=False,
    )
    assert execute.call_count == 3 and output.is_file()
    if invalidate == "input_changed":
        source.write_text(">sample\nAAAA\n", encoding="utf-8")
    else:
        output.unlink()
    inputs = _model(monkeypatch, [file_action, last_action, _finish()])
    result = run_agent(
        "Calculate masses and reverse FASTA",
        config=CONFIG,
        registry=registry,
        cache_dir=cache,
        reuse_results=True,
        allow_side_effects=True,
        verbose=False,
    )
    assert result.reused_steps == 1
    assert len(inputs[0]["history"]) == 5
    assert execute.call_count == 5
    expected = "TTTT" if invalidate == "input_changed" else "CGTT"
    assert output.read_text(encoding="utf-8") == f">sample\n{expected}\n"


def test_successful_python_step_is_not_executed_again_after_interruption(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    config = PythonExecutionConfig(tmp_path / "work")
    cache = tmp_path / "cache"
    _model(monkeypatch, [_python(), KeyboardInterrupt()])
    with pytest.raises(KeyboardInterrupt):
        run_agent(
            "Prepare FASTA",
            config=CONFIG,
            registry=registry,
            cache_dir=cache,
            code_execution=config,
            allow_side_effects=True,
            verbose=False,
        )
    _model(monkeypatch, [_python(), _finish()])
    result = run_agent(
        "Prepare FASTA",
        config=CONFIG,
        registry=registry,
        cache_dir=cache,
        reuse_results=True,
        code_execution=config,
        allow_side_effects=True,
        verbose=False,
    )
    assert result.reused_steps == 1 and result.status == "completed"
    assert _observation(result)["reused"] is True
    assert (tmp_path / "work/counter.txt").read_text() == "x"
    assert result.execution_dir and not list(Path(result.execution_dir).glob("attempt_*"))


@pytest.mark.parametrize("allow_side_effects", [False, True])
def test_disabling_python_cannot_reuse_its_previous_success(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
    allow_side_effects: bool,
) -> None:
    cache = tmp_path / "cache"
    _model(monkeypatch, [_python(), _finish()])
    run_agent(
        "Prepare FASTA",
        config=CONFIG,
        registry=registry,
        cache_dir=cache,
        code_execution=PythonExecutionConfig(tmp_path / "work"),
        allow_side_effects=True,
        verbose=False,
    )
    inputs = _model(monkeypatch, [_finish(), _finish("blocked")])
    result = run_agent(
        "Prepare FASTA",
        config=CONFIG,
        registry=registry,
        cache_dir=cache,
        reuse_results=True,
        allow_side_effects=allow_side_effects,
        verbose=False,
    )
    assert result.reused_steps == 0 and result.status == "blocked"
    assert inputs[0]["history"] == []
    assert _observation(result)["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"
    assert (tmp_path / "work/counter.txt").read_text() == "x"


def test_corrupt_checkpoint_is_reported_before_llm_call(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    _model(monkeypatch, [*_discover(WEIGHT), _call(WEIGHT, sequence="ACGT"), _finish()])
    first = run_agent(
        "Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False
    )
    assert first.cache_dir
    (Path(first.cache_dir) / "checkpoint.json").write_text("{broken", encoding="utf-8")
    inputs = _model(monkeypatch, [])
    with pytest.raises(ConfigurationError, match="AGENT_CHECKPOINT_INVALID"):
        run_agent(
            "Calculate mass",
            config=CONFIG,
            registry=registry,
            cache_dir=tmp_path,
            reuse_results=True,
            verbose=False,
        )
    assert inputs == []


def test_cache_does_not_store_config_credentials_or_transport_error_text(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    secret = "offline-placeholder"
    _model(
        monkeypatch,
        [
            *_discover(WEIGHT),
            _call(WEIGHT, sequence="ACGT"),
            BackendExecutionError(f"Transport rejected token {secret}"),
        ],
    )
    with pytest.raises(BackendExecutionError):
        run_agent(
            "Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False
        )
    files = [path for path in tmp_path.rglob("*") if path.is_file()]
    assert files
    for path in files:
        assert secret.encode() not in path.read_bytes()


def test_pending_interrupted_call_is_neither_success_nor_automatically_repeated(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    execute = Mock(side_effect=KeyboardInterrupt())
    monkeypatch.setattr(registry, "execute", execute)
    action = _call(WEIGHT, sequence="ACGT")
    _model(monkeypatch, [*_discover(WEIGHT), action])
    with pytest.raises(KeyboardInterrupt):
        run_agent(
            "Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False
        )
    inputs = _model(monkeypatch, [action, _finish(), _finish("blocked")])
    result = run_agent(
        "Calculate mass",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path,
        reuse_results=True,
        verbose=False,
    )
    execute.assert_called_once()
    assert result.status == "blocked" and result.reused_steps == 0
    assert inputs[0]["resume"]["uncertain_actions"]
    assert _observation(result, -2)["error"]["code"] == "AGENT_PREVIOUS_ACTION_UNCERTAIN"
    assert _observation(result)["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"


def test_rag_retrieves_sources_before_first_llm_and_supports_direct_answer(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    documents = tmp_path / "documents"
    documents.mkdir()
    source = documents / "formats.md"
    source.write_text("FASTA records begin with a header line starting with >.\n", encoding="utf-8")
    index = tmp_path / "index"
    execute = Mock(wraps=registry.execute)
    monkeypatch.setattr(registry, "execute", execute)
    inputs = _model(
        monkeypatch,
        [
            {
                "action": "finish",
                "status": "completed",
                "text": "FASTA headers start with >. Source: formats.md",
            }
        ],
    )
    result = run_agent(
        "How do FASTA records begin?",
        config=CONFIG,
        registry=registry,
        rag=RAGConfig(documents, index_dir=index),
        verbose=False,
    )
    assert result.status == "completed" and result.steps == ()
    assert result.knowledge_context
    assert inputs[0]["knowledge_context"] == list(result.knowledge_context)
    passage = result.knowledge_context[0]
    assert str(passage["source"]).endswith("formats.md") and passage["chunk_id"]
    assert "FASTA" in str(passage)
    assert list(index.rglob("*.json"))
    execute.assert_not_called()


def test_empty_rag_match_does_not_support_invented_completion(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    registry: AgentToolRegistry,
) -> None:
    documents = tmp_path / "documents"
    documents.mkdir()
    (documents / "formats.txt").write_text("FASTA headers begin with >.\n", encoding="utf-8")
    inputs = _model(monkeypatch, [_finish(), _finish("blocked")])
    result = run_agent(
        "zzzxxyyqq",
        config=CONFIG,
        registry=registry,
        rag=RAGConfig(documents, index_dir=tmp_path / "index"),
        verbose=False,
    )
    assert result.status == "blocked" and result.knowledge_context == ()
    assert inputs[0]["knowledge_context"] == []
    assert _observation(result, 0)["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"


def test_changed_knowledge_does_not_restore_work_based_on_old_evidence(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, registry: AgentToolRegistry
) -> None:
    documents = tmp_path / "documents"
    documents.mkdir()
    source = documents / "mass.md"
    source.write_text("Mass analysis uses the original lab procedure.", encoding="utf-8")
    rag = RAGConfig(documents, index_dir=tmp_path / "index")
    _model(monkeypatch, [*_discover(WEIGHT), _call(WEIGHT, sequence="ACGT"), _finish()])
    run_agent(
        "Calculate mass according to the lab procedure",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path / "cache",
        rag=rag,
        verbose=False,
    )
    source.write_text("Mass analysis uses the revised lab procedure.", encoding="utf-8")
    inputs = _model(monkeypatch, [_finish("blocked")])
    result = run_agent(
        "Calculate mass according to the lab procedure",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path / "cache",
        reuse_results=True,
        rag=rag,
        verbose=False,
    )
    assert result.reused_steps == 0 and inputs[0]["history"] == []
    assert "revised" in str(inputs[0]["knowledge_context"])


def test_rag_evidence_cannot_override_a_failed_tool_execution(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, registry: AgentToolRegistry
) -> None:
    documents = tmp_path / "documents"
    documents.mkdir()
    (documents / "mass.txt").write_text("Calculate mass from a DNA sequence.", encoding="utf-8")
    _model(monkeypatch, [*_discover(WEIGHT), _call(WEIGHT), _finish(), _finish("blocked")])
    result = run_agent(
        "Calculate mass",
        config=CONFIG,
        registry=registry,
        rag=RAGConfig(documents, index_dir=tmp_path / "index"),
        verbose=False,
    )
    assert result.status == "blocked" and result.knowledge_context
    assert _observation(result)["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"


def test_changed_tool_implementation_invalidates_its_previous_results(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, registry: AgentToolRegistry
) -> None:
    _model(monkeypatch, [*_discover(WEIGHT), _call(WEIGHT, sequence="ACGT"), _finish()])
    run_agent("Calculate mass", config=CONFIG, registry=registry, cache_dir=tmp_path, verbose=False)

    def replacement(sequence: str) -> float:
        return 123.0

    changed = AgentToolRegistry(
        (registry.resolve(REVERSE), replace(registry.resolve(WEIGHT), function=replacement))
    )
    inputs = _model(monkeypatch, [_finish("blocked")])
    result = run_agent(
        "Calculate mass",
        config=CONFIG,
        registry=changed,
        cache_dir=tmp_path,
        reuse_results=True,
        verbose=False,
    )
    assert result.reused_steps == 0
    assert len(inputs[0]["history"]) == 1
    assert inputs[0]["resume"]["invalidated_from"] == 2


def test_declared_output_is_checked_when_tool_result_does_not_include_its_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, registry: AgentToolRegistry
) -> None:
    output = tmp_path / "output.fna"

    def write(*args: Any, **kwargs: Any) -> int:
        output.write_text(">sample\nACGT\n", encoding="utf-8")
        return 1

    execute = Mock(side_effect=write)
    monkeypatch.setattr(registry, "execute", execute)
    action = _call(REVERSE, value="ACGT", output_path=str(output))
    _model(monkeypatch, [*_discover(REVERSE), action, _finish()])
    run_agent(
        "Export sequence",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path / "cache",
        verbose=False,
    )
    output.unlink()
    inputs = _model(monkeypatch, [action, _finish()])
    result = run_agent(
        "Export sequence",
        config=CONFIG,
        registry=registry,
        cache_dir=tmp_path / "cache",
        reuse_results=True,
        verbose=False,
    )
    assert result.reused_steps == 0 and len(inputs[0]["history"]) == 2
    assert execute.call_count == 2 and output.is_file()
