"""Offline action-loop tests with real discovery and deterministic DNAKit tools."""

from __future__ import annotations

import io
from copy import deepcopy
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
    EmbeddingConfig,
    default_tool_registry,
    run_agent,
)

CONFIG = LLMConfig(provider="deepseek", model="offline-test", api_key="offline-placeholder")
WEIGHT = "dnakit.thermodynamics.molecular_weight"
DOWNLOAD = "dnakit.references.download_genome"


def _model(monkeypatch: pytest.MonkeyPatch, actions: list[dict[str, Any]]) -> list[Any]:
    pending = iter(actions)
    inputs = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        inputs.append(deepcopy(kwargs["input_data"]))
        return LLMResult(
            text="JSON action",
            model="offline-test",
            response_id=None,
            finish_reason="stop",
            output_format="json",
            _json_data=FrozenDict(next(pending)),
        )

    monkeypatch.setattr(agent, "call_llm", reply)
    return inputs


def _steps(name: str, arguments: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"action": "search", "query": name.rsplit(".", 1)[-1]},
        {"action": "describe", "tool_name": name},
        {"action": "call", "tool_name": name, "arguments": arguments},
    ]


def _finish(status: str = "completed") -> dict[str, str]:
    return {"action": "finish", "status": status, "text": "Summary of actual observations"}


def _observation(result: AgentResult, index: int) -> dict[str, Any]:
    return cast(dict[str, Any], result.steps[index]["observation"])


def test_natural_language_task_uses_multiple_categories(monkeypatch: pytest.MonkeyPatch) -> None:
    actions = _steps(WEIGHT, {"sequence": "ACGT"})
    actions += _steps("dnakit.ops.reverse_complement", {"value": "AACG"})
    inputs = _model(monkeypatch, [*actions, _finish()])
    result = run_agent(
        "Calculate mass for ACGT and reverse complement AACG", config=CONFIG, verbose=False
    )
    assert result.status == "completed"
    assert len(result.steps) == 6
    assert _observation(result, 2)["result"]["value_dalton"] == pytest.approx(1173.84)
    assert _observation(result, 5)["result"]["records"][0]["sequence"]["parts"] == ["CGTT"]
    assert "description" in inputs[2]["history"][1]["observation"]["tool"]["input_schema"]
    assert "api_key" not in inputs[0] and "config" not in inputs[0]


def test_invalid_arguments_are_returned_for_correction(monkeypatch: pytest.MonkeyPatch) -> None:
    actions = _steps(WEIGHT, {"sequence": "ACGT", "strand": "invalid"})
    actions.append({"action": "call", "tool_name": WEIGHT, "arguments": {"sequence": "ACGT"}})
    inputs = _model(monkeypatch, [*actions, _finish()])
    result = run_agent("Calculate molecular weight", config=CONFIG, verbose=False)
    assert inputs[3]["history"][-1]["observation"]["ok"] is False
    assert result.status == "completed"
    assert _observation(result, -1)["ok"] is True


def test_model_cannot_grant_itself_write_permission(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    actions = _steps(DOWNLOAD, {"query": "hg38", "output_dir": str(tmp_path / "genome")})
    actions[-1]["allow_side_effects"] = True
    _model(monkeypatch, [*actions, _finish("blocked")])
    result = run_agent("Download hg38", config=CONFIG, verbose=False)
    assert result.status == "blocked"
    assert _observation(result, -1)["error"]["code"] == "AGENT_TOOL_CONFIRMATION_REQUIRED"
    assert not (tmp_path / "genome").exists()


@pytest.mark.parametrize("failed", [False, True])
def test_identical_calls_do_not_repeat_side_effects(
    monkeypatch: pytest.MonkeyPatch,
    failed: bool,
) -> None:
    registry = default_tool_registry(verbose=False)
    execute = Mock(
        return_value={"fasta_path": "mock.fna"},
        side_effect=OSError("simulated partial write") if failed else None,
    )
    monkeypatch.setattr(registry, "execute", execute)
    arguments = {"query": "hg38", "output_dir": "mock-directory"}
    actions = _steps(DOWNLOAD, arguments)
    _model(monkeypatch, [*actions, actions[-1], _finish("blocked" if failed else "completed")])
    result = run_agent(
        "Download hg38", config=CONFIG, registry=registry, allow_side_effects=True, verbose=False
    )
    execute.assert_called_once_with(DOWNLOAD, arguments, allow_side_effects=True)
    assert _observation(result, -1)["reused"] is True
    assert _observation(result, -1)["ok"] is not failed


@pytest.mark.parametrize("discover", [False, True])
def test_calls_require_discovery_and_schema(
    monkeypatch: pytest.MonkeyPatch,
    discover: bool,
) -> None:
    actions = _steps(WEIGHT, {"sequence": "ACGT"})
    inputs = _model(
        monkeypatch, ([actions[0]] if discover else []) + [actions[2], _finish("blocked")]
    )
    result = run_agent("Calculate mass", config=CONFIG, verbose=False)
    code = "AGENT_TOOL_NOT_DESCRIBED" if discover else "AGENT_TOOL_NOT_DISCOVERED"
    assert inputs[-1]["history"][-1]["observation"]["error"]["code"] == code
    assert result.status == "blocked"


def test_completed_without_execution_is_rejected(monkeypatch: pytest.MonkeyPatch) -> None:
    inputs = _model(monkeypatch, [_finish(), _finish("blocked")])
    result = run_agent("Download hg38", config=CONFIG, verbose=False)
    assert result.status == "blocked"
    assert inputs[1]["history"][0]["observation"]["error"]["code"] == "AGENT_COMPLETION_UNSUPPORTED"


def test_step_limit_stops_without_extra_llm_requests(monkeypatch: pytest.MonkeyPatch) -> None:
    inputs = _model(monkeypatch, [{"action": "search", "query": "weight"}] * 2)
    with pytest.raises(BackendExecutionError, match="AGENT_STEP_LIMIT"):
        run_agent("Calculate mass", config=CONFIG, max_steps=2, verbose=False)
    assert len(inputs) == 2


@pytest.mark.parametrize("limit", [0, -1, True, 1.5])
def test_invalid_step_limits_fail_before_llm(monkeypatch: pytest.MonkeyPatch, limit: Any) -> None:
    inputs = _model(monkeypatch, [])
    with pytest.raises(ConfigurationError, match="AGENT_STEP_LIMIT_INVALID"):
        run_agent("Calculate mass", config=CONFIG, max_steps=limit, verbose=False)
    assert inputs == []


@pytest.mark.parametrize("error", [KeyboardInterrupt(), BackendExecutionError("LLM failed")])
def test_llm_errors_and_interrupts_propagate(
    monkeypatch: pytest.MonkeyPatch, error: BaseException
) -> None:
    llm = Mock(side_effect=error)
    monkeypatch.setattr(agent, "call_llm", llm)
    with pytest.raises(type(error)):
        run_agent("Calculate mass", config=CONFIG, verbose=False)
    llm.assert_called_once()


@pytest.mark.parametrize("mode", ["keyword", "semantic", "hybrid"])
def test_agent_uses_configured_search_without_sharing_embedding_credentials(
    monkeypatch: pytest.MonkeyPatch,
    mode: Any,
) -> None:
    embedding = EmbeddingConfig(
        model="embedding-only-model",
        base_url="https://embedding.example/v1",
        api_key="embedding-only-key",
        cache_dir=None,
    )
    registry = AgentToolRegistry(
        (default_tool_registry(verbose=False).resolve(WEIGHT),),
        embedding_config=embedding,
    )
    embed = Mock(side_effect=lambda texts, config: [(1.0, 0.0)] * len(texts))
    monkeypatch.setattr("dnakit.tools.embeddings._embed", embed)
    inputs = _model(monkeypatch, [*_steps(WEIGHT, {"sequence": "ACGT"}), _finish()])
    result = run_agent(
        "Calculate molecular weight",
        config=CONFIG,
        registry=registry,
        search_mode=mode,
        search_limit=1,
        verbose=False,
    )
    assert result.status == "completed"
    assert len(_observation(result, 0)["tools"]) == 1
    assert inputs[0]["search"] == {"mode": mode, "limit": 1}
    assert "embedding-only-key" not in str(inputs)
    assert "embedding-only-model" not in str(inputs)
    assert embed.call_count == (0 if mode == "keyword" else 2)


@pytest.mark.parametrize(
    "options,code",
    [
        ({"search_mode": "semantic"}, "MISSING_EMBEDDING_CONFIG"),
        ({"search_mode": "bad-mode"}, "AGENT_SEARCH_MODE_INVALID"),
        ({"search_limit": 0}, "AGENT_TOOL_LIMIT_INVALID"),
        ({"search_limit": True}, "AGENT_TOOL_LIMIT_INVALID"),
        ({"search_limit": 101}, "AGENT_TOOL_LIMIT_INVALID"),
    ],
)
def test_invalid_search_configuration_fails_before_chat_request(
    monkeypatch: pytest.MonkeyPatch,
    options: dict[str, Any],
    code: str,
) -> None:
    inputs = _model(monkeypatch, [])
    with pytest.raises(ConfigurationError, match=code):
        run_agent("Calculate mass", config=CONFIG, verbose=False, **options)
    assert inputs == []


def test_name_only_task_passes_observed_urls_and_files_through_full_ml_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mock only LLM/network; discover, download, fit, predict and plot for real."""
    np = pytest.importorskip("numpy")
    pytest.importorskip("sklearn")
    from dnakit.download import files
    from dnakit.llm.client import _request_body
    from dnakit.search import web

    rng = np.random.default_rng(5)
    payloads = {}
    for split, count in (("train", 100), ("test", 12)):
        rows = ["id,sequence,label"]
        for index in range(count):
            sequence = "".join(rng.choice(list("ACGT"), size=50))
            rows.append(f"{split}-{index},{sequence},{index % 2}")
        payloads[f"https://data.example/files/{split}.csv"] = ("\n".join(rows) + "\n").encode()
    monkeypatch.setattr(
        web,
        "_request",
        lambda *args, **kwargs: (
            "<p>Example revised classification dataset</p>"
            + "".join(f'<a href="{url}">{url.rsplit("/", 1)[-1]}</a>' for url in payloads),
            "https://data.example/dataset",
            "text/html",
            200,
        ),
    )

    class DownloadResponse(io.BytesIO):
        def __init__(self, url: str) -> None:
            super().__init__(payloads[url])
            self.url = url
            self.headers = {"Content-Length": str(len(payloads[url]))}

        def geturl(self) -> str:
            return self.url

    downloaded_urls: list[str] = []

    def download(request: Any, **kwargs: Any) -> DownloadResponse:
        downloaded_urls.append(request.full_url)
        return DownloadResponse(request.full_url)

    monkeypatch.setattr(files, "urlopen", download)
    observed: dict[str, Any] = {}
    observed_steps: dict[str, int] = {}
    train_name = "dnakit.predictions.train_custom_property_model"
    predict_name = "dnakit.predictions.predict_custom_properties"
    plot_name = "dnakit.visualization.plot_csv_distributions"
    download_name = "dnakit.download.download_url"
    read_name = "dnakit.search.read_webpage"

    def decisions() -> Any:
        # Discover the downloader first; resolve its missing URL using a candidate page.
        yield from _steps(download_name, {})[:2]
        yield from _steps(read_name, {"url": "https://data.example/dataset"})
        links = observed[read_name]["records"]
        for link in links:
            yield {
                "action": "call",
                "tool_name": download_name,
                "arguments": {
                    "url": link["url"],
                    "output_dir": str(tmp_path),
                    "verbose": False,
                },
            }
        yield from _steps(
            train_name,
            {
                "source": str(tmp_path / "train.csv"),
                "config": {
                    "task_type": "classification",
                    "model": "svm",
                    "finetune_method": "machine_learning",
                    "cv_folds": 10,
                    "model_params": {"kernel": "linear", "probability": True},
                    "output_dir": str(tmp_path / "model"),
                    "workers": 1,
                    "show_progress": False,
                },
                "verbose": False,
            },
        )
        yield from _steps(
            predict_name,
            {
                "source": str(tmp_path / "test.csv"),
                "model_dir": observed[train_name]["output_dir"],
                "trust_model": True,
                "output_path": str(tmp_path / "predictions.csv"),
                "verbose": False,
            },
        )
        yield from _steps(
            plot_name,
            {
                "source": observed[predict_name]["target_path"],
                "columns": ["label"],
                "task": "classification",
                "config": {"categories": ["0", "1"]},
                "verbose": False,
            },
        )
        yield from _steps(
            "dnakit.visualization.save_image",
            {
                "artifact": {"$ref": {"step": observed_steps[plot_name], "path": ["label"]}},
                "target": str(tmp_path / "distribution.svg"),
                "verbose": False,
            },
        )
        yield _finish()

    actions = decisions()
    request_sizes = []

    def reply(*args: Any, **kwargs: Any) -> LLMResult:
        history = kwargs["input_data"]["history"]
        if history and history[-1]["action"]["action"] == "call":
            assert history[-1]["observation"]["ok"], history[-1]["observation"]
            observed[history[-1]["action"]["tool_name"]] = history[-1]["observation"]["result"]
            observed_steps[history[-1]["action"]["tool_name"]] = history[-1]["step"]
        request_sizes.append(
            len(
                _request_body(
                    args[0],
                    kwargs["input_data"],
                    config=CONFIG,
                    output_format="json",
                    system_prompt=kwargs["system_prompt"],
                )
            )
        )
        return LLMResult(
            text="action",
            model="offline-test",
            response_id=None,
            finish_reason="stop",
            output_format="json",
            _json_data=FrozenDict(next(actions)),
        )

    monkeypatch.setattr(agent, "call_llm", reply)
    result = run_agent(
        f"Download Example revised to {tmp_path}, run 10-fold SVM on train, "
        "predict test and plot label proportions.",
        config=CONFIG,
        allow_side_effects=True,
        verbose=False,
    )
    assert result.status == "completed"
    assert set(downloaded_urls) == set(payloads)
    assert observed[train_name]["split_counts"] == {"train": 100, "validation": 0, "test": 0}
    assert observed[predict_name]["row_count"] == 12
    plot = observed[plot_name]["label"]["metadata"]
    assert plot["labels"] == ["0", "1"]
    assert sum(plot["counts"]) == 12
    assert sum(plot["proportions"]) == pytest.approx(1)
    assert (tmp_path / "distribution.svg").is_file()
    assert max(request_sizes) < 100_000
    assert "predictions" not in observed[predict_name]


def test_result_references_reuse_observed_values_and_share_call_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    actions = _steps("dnakit.ops.reverse_complement", {"value": "AACG"})
    actions += _steps(
        WEIGHT,
        {
            "sequence": {
                "$ref": {
                    "step": 3,
                    "path": ["records", 0, "sequence", "parts", 0],
                }
            }
        },
    )
    actions += [{"action": "call", "tool_name": WEIGHT, "arguments": {"sequence": "CGTT"}}]
    _model(monkeypatch, [*actions, _finish()])
    result = run_agent("Reverse complement and calculate its mass", config=CONFIG, verbose=False)
    assert result.status == "completed"
    assert _observation(result, -1)["reused"] is True
    assert result.steps[2]["step"] == 3


@pytest.mark.parametrize(
    "reference",
    [
        {"step": 0, "path": []},
        {"step": 2, "path": []},
        {"step": True, "path": []},
        {"step": 1, "path": ["missing"]},
        {"step": 1, "path": ["items", -1]},
        {"step": 1, "path": ["items", True]},
        {"step": 1, "path": "items"},
        {"step": 1},
    ],
)
def test_invalid_result_references_fail_before_execution(reference: Any) -> None:
    history = [
        {"action": {"action": "call"}, "observation": {"ok": True, "result": {"items": ["x"]}}}
    ]
    with pytest.raises(ConfigurationError, match="AGENT_REFERENCE_INVALID"):
        agent._resolve_arguments({"input": {"$ref": reference}}, history)


def test_reference_resolution_does_not_interpret_data_as_new_references() -> None:
    nested = {"$ref": {"step": 999, "path": []}}
    history = [{"action": {"action": "call"}, "observation": {"ok": True, "result": nested}}]
    assert agent._resolve_arguments({"$ref": {"step": 1, "path": []}}, history) == nested
    literal = {"$ref": "#/schema/Sequence"}
    assert agent._resolve_arguments(literal, []) == literal
    history[0]["observation"]["ok"] = False
    with pytest.raises(ConfigurationError, match="AGENT_REFERENCE_INVALID"):
        agent._resolve_arguments({"$ref": {"step": 1, "path": []}}, history)
