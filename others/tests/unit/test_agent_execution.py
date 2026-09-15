"""Real subprocesses and file tools, with deterministic offline model decisions."""

from __future__ import annotations

import json
import os
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

import dnakit.tools.agent as agent
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.llm import LLMConfig, LLMResult
from dnakit.tools import PythonExecutionConfig, run_agent
from dnakit.tools.execution import _PythonExecutor

CONFIG = LLMConfig(provider="deepseek", model="offline-test", api_key="offline-placeholder")


def _model(monkeypatch: pytest.MonkeyPatch, actions: list[dict[str, Any]]) -> list[Any]:
    pending = iter(actions)
    observations = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        observations.append(
            {"input_data": deepcopy(kwargs["input_data"]), "system_prompt": kwargs["system_prompt"]}
        )
        return LLMResult(
            text="action",
            model="offline-test",
            response_id=None,
            finish_reason="stop",
            output_format="json",
            _json_data=FrozenDict(next(pending)),
        )

    monkeypatch.setattr(agent, "call_llm", reply)
    return observations


def _action(code: str, outputs: list[str] | None = None, inputs: Any = None) -> dict[str, Any]:
    return {
        "action": "execute_python",
        "reason": "Convert observed data for the next tool",
        "code": code,
        "inputs": inputs or {},
        "expected_outputs": outputs or [],
    }


def _finish(status: str = "completed") -> dict[str, str]:
    return {"action": "finish", "status": status, "text": "Actual outputs checked"}


def test_tool_result_code_failure_repair_and_next_file_tool(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    name = "dnakit.ops.reverse_complement"
    inputs = {"observed": {"$ref": {"step": 3, "path": []}}}
    load = "import json, sys\nfrom pathlib import Path\nd=json.load(open(sys.argv[1]))\n"
    broken = _action(load + "sequence=d['wrong_key']\n", ["converted.fna"], inputs)
    corrected = _action(
        load + "sequence=''.join(d['observed']['records'][0]['sequence']['parts'])\n"
        "assert sequence == 'CGTT'\nPath('converted.fna').write_text('>x\\n'+sequence+'\\n')\n",
        ["converted.fna"],
        inputs,
    )
    calls = _model(
        monkeypatch,
        [
            {"action": "search", "query": name},
            {"action": "describe", "tool_name": name},
            {"action": "call", "tool_name": name, "arguments": {"value": "AACG"}},
            broken,
            corrected,
            {
                "action": "call",
                "tool_name": name,
                "arguments": {
                    "value": {"$ref": {"step": 5, "path": ["outputs", 0, "path"]}},
                    "output_path": str(tmp_path / "restored.fna"),
                },
            },
            _finish(),
        ],
    )
    result = run_agent(
        "Convert and reverse back",
        config=CONFIG,
        allow_side_effects=True,
        code_execution=PythonExecutionConfig(tmp_path),
        verbose=False,
    )
    assert result.status == "completed"
    assert (tmp_path / "restored.fna").read_text() == ">x\nAACG\n"
    failed = result.steps[3]["observation"]
    assert failed["ok"] is False
    assert "KeyError" in failed["result"]["stderr"]
    assert calls[4]["input_data"]["history"][-1]["observation"] == failed
    assert "Python fallback is enabled" in calls[0]["system_prompt"]
    assert result.execution_dir
    directory = Path(result.execution_dir)
    assert (
        json.loads((directory / "attempt_002/inputs.json").read_text())["observed"]
        == (result.steps[2]["observation"]["result"])
    )
    assert json.loads((directory / "summary.json").read_text())["code_attempts"] == 2
    assert len((directory / "history.jsonl").read_text().splitlines()) == 6
    assert "offline-placeholder" not in (directory / "settings.json").read_text()


def test_model_cannot_enable_execution(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    marker = tmp_path / "unauthorized.txt"
    _model(monkeypatch, [_action(f"open({str(marker)!r}, 'w').write('x')"), _finish("blocked")])
    result = run_agent("Task", config=CONFIG, allow_side_effects=True, verbose=False)
    assert result.steps[0]["observation"]["error"]["code"] == "AGENT_CODE_EXECUTION_DISABLED"
    assert not marker.exists()


def test_execution_requires_caller_side_effect_permission(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls = _model(monkeypatch, [])
    with pytest.raises(ConfigurationError, match="AGENT_TOOL_CONFIRMATION_REQUIRED"):
        run_agent(
            "Task", config=CONFIG, code_execution=PythonExecutionConfig(tmp_path), verbose=False
        )
    assert not calls
    assert not (tmp_path / ".dnakit_agent").exists()


def test_failed_code_cannot_support_completion_and_duplicate_is_not_reexecuted(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    action = _action(
        "from pathlib import Path\np=Path('counter')\n"
        "p.write_text(p.read_text()+'x' if p.exists() else 'x')\n"
        "raise ValueError('actual error')"
    )
    _model(monkeypatch, [action, {**action, "reason": "Try again"}, _finish(), _finish("blocked")])
    result = run_agent(
        "Task",
        config=CONFIG,
        allow_side_effects=True,
        code_execution=PythonExecutionConfig(tmp_path),
        verbose=False,
    )
    assert (tmp_path / "counter").read_text() == "x"
    assert result.steps[1]["observation"]["reused"] is True
    assert result.steps[2]["observation"]["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"


def test_attempt_limit_is_enforced_across_different_scripts(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _model(
        monkeypatch,
        [_action("raise ValueError('bad')"), _action("print('fixed')"), _finish("blocked")],
    )
    result = run_agent(
        "Task",
        config=CONFIG,
        allow_side_effects=True,
        code_execution=PythonExecutionConfig(tmp_path, max_attempts=1),
        verbose=False,
    )
    assert result.steps[1]["observation"]["error"]["code"] == "AGENT_CODE_ATTEMPT_LIMIT"
    assert len(list(Path(result.execution_dir).glob("attempt_*"))) == 1


def test_missing_artifact_is_failure_even_when_process_exits_zero(tmp_path: Path) -> None:
    executor = _PythonExecutor(PythonExecutionConfig(tmp_path))
    result = executor.execute(
        code="print('done')",
        reason="Save file",
        inputs={},
        expected_outputs=["missing.csv"],
        step=1,
    )
    assert result["status"] == "missing_outputs"
    assert result["returncode"] == 0


def test_traceback_and_output_are_saved_but_model_output_is_bounded(tmp_path: Path) -> None:
    executor = _PythonExecutor(PythonExecutionConfig(tmp_path, max_output_bytes=100))
    result = executor.execute(
        code="print('x'*10000)\nraise RuntimeError('bad format')",
        reason="Inspect",
        inputs={},
        expected_outputs=[],
        step=1,
    )
    assert result["status"] == "failed"
    assert len(result["stdout"]) == 100 and result["stdout_truncated"] is True
    assert Path(result["stdout_path"]).stat().st_size > 10000
    assert "RuntimeError: bad format" in result["stderr"]


@pytest.mark.skipif(os.name != "posix", reason="POSIX process group cleanup")
def test_timeout_kills_child_process_and_preserves_logs(tmp_path: Path) -> None:
    executor = _PythonExecutor(PythonExecutionConfig(tmp_path, timeout_seconds=1))
    result = executor.execute(
        code="import subprocess, sys, time\nfrom pathlib import Path\n"
        "p=subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
        "Path('child.pid').write_text(str(p.pid))\nprint('started', flush=True)\n"
        "time.sleep(60)",
        reason="Timeout exercise",
        inputs={},
        expected_outputs=[],
        step=1,
    )
    assert result["status"] == "timeout" and result["elapsed_seconds"] < 10
    assert "started" in result["stdout"]
    pid = int((tmp_path / "child.pid").read_text())
    # A reparented zombie may await init's reap, but it must not be running.
    import time

    for _ in range(100):
        status = Path(f"/proc/{pid}/status")
        try:
            process_status = status.read_text()
        except FileNotFoundError:
            break  # The terminated process can be reaped between filesystem probes.
        if "State:\tZ" in process_status:
            break
        time.sleep(0.01)
    else:
        pytest.fail("Child survived timeout")


def test_llm_failure_preserves_previous_execution_and_summary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _model(monkeypatch, [_action("print('observed')")])
    with pytest.raises(StopIteration):
        run_agent(
            "Task",
            config=CONFIG,
            allow_side_effects=True,
            code_execution=PythonExecutionConfig(tmp_path),
            verbose=False,
        )
    directory = next((tmp_path / ".dnakit_agent").iterdir())
    assert (directory / "attempt_001/result.json").is_file()
    assert json.loads((directory / "summary.json").read_text())["status"] == "interrupted"


@pytest.mark.parametrize(
    "options",
    [
        {"timeout_seconds": 0},
        {"max_attempts": True},
        {"max_output_bytes": 100000},
    ],
)
def test_invalid_limits(options: dict[str, Any], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="AGENT_EXECUTION_CONFIG"):
        PythonExecutionConfig(tmp_path, **options)
