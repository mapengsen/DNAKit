"""Offline regression tests for on-demand weights and prediction feedback."""

from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import DownloadError
from dnakit.predictions import (
    BiologicalSequence,
    PredictionOutput,
    PropertyPredictionConfig,
    available_enformer_benchmark_tasks,
    ensure_prediction_checkpoint,
    get_enformer_benchmark_task,
    predict_enformer_benchmark,
)
from dnakit.representations import ensure_model_checkpoint


@pytest.mark.parametrize("task", available_enformer_benchmark_tasks(verbose=False))
def test_download_only_selected_task_then_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any, task: str
) -> None:
    spec = get_enformer_benchmark_task(task, verbose=False)
    calls: list[str] = []

    def download(**kwargs: Any) -> Any:
        calls.append(kwargs["id"])
        with zipfile.ZipFile(kwargs["output"], "w") as archive:
            archive.writestr("checkpoint/data.pkl", b"test fixture")
        return kwargs["output"]

    monkeypatch.setitem(sys.modules, "gdown", SimpleNamespace(download=download))
    monkeypatch.chdir(tmp_path)
    config = PropertyPredictionConfig(model="enformer", task=task)
    first = ensure_prediction_checkpoint(config)
    path = tmp_path / "ckpt" / "enformer-benchmarks" / spec.checkpoint_filename
    assert first.paths == (str(path),)
    assert first.sources == (spec.checkpoint_url,)
    assert first.downloaded
    messages = capsys.readouterr().err
    assert "Starting checkpoint download" in messages
    assert "Checkpoint download complete" in messages
    assert str(path) in messages and spec.checkpoint_url in messages
    second = ensure_prediction_checkpoint(config)
    assert not second.downloaded
    assert calls == [spec.checkpoint_file_id]
    assert "Using local checkpoint" in capsys.readouterr().err
    assert list(path.parent.iterdir()) == [path]


@pytest.mark.parametrize("payload", [b"", b"<html>Access denied</html>", b"PK\x03\x04partial"])
def test_invalid_download_never_becomes_cached_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any, payload: bytes
) -> None:
    def download(**kwargs: Any) -> Any:
        kwargs["output"].write(payload)
        return kwargs["output"]

    monkeypatch.setitem(sys.modules, "gdown", SimpleNamespace(download=download))
    config = PropertyPredictionConfig(model="enformer", task="enhancers", checkpoint_dir=tmp_path)
    with pytest.raises(DownloadError) as error:
        ensure_prediction_checkpoint(config)
    assert error.value.code == "MODEL_CHECKPOINT_DOWNLOAD_FAILED"
    assert "15fkN_7ZftfdHzi0In9_NcskOm2ZvjZdF" in error.value.context["download_url"]
    assert not list(tmp_path.iterdir())
    assert "Checkpoint download complete" not in capsys.readouterr().err


@pytest.mark.parametrize("error", [OSError("connection lost"), KeyboardInterrupt()])
def test_interrupted_download_cleans_partial_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error: BaseException
) -> None:
    def download(**kwargs: Any) -> Any:
        kwargs["output"].write(b"partial")
        raise error

    monkeypatch.setitem(sys.modules, "gdown", SimpleNamespace(download=download))
    config = PropertyPredictionConfig(model="enformer", task="enhancers", checkpoint_dir=tmp_path)
    with pytest.raises(
        KeyboardInterrupt if isinstance(error, KeyboardInterrupt) else DownloadError
    ):
        ensure_prediction_checkpoint(config, verbose=False)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("subdir", ["", "enformer-benchmarks", "nt-revised"])
def test_custom_layout_and_quiet_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any, subdir: str
) -> None:
    path = tmp_path / subdir / "H2AFZ.ckpt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"existing")
    monkeypatch.setitem(sys.modules, "gdown", None)
    config = PropertyPredictionConfig(model="enformer", task="h2afz", checkpoint_dir=tmp_path)
    info = ensure_prediction_checkpoint(config, verbose=False)
    assert info.paths == (str(path),)
    assert capsys.readouterr().err == ""


def test_explicit_missing_file_downloads_to_requested_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any
) -> None:
    path = tmp_path / "custom" / "enhancers.ckpt"

    def download(**kwargs: Any) -> Any:
        assert kwargs["quiet"] is True
        with zipfile.ZipFile(kwargs["output"], "w") as archive:
            archive.writestr("checkpoint/data.pkl", b"fixture")
        return kwargs["output"]

    monkeypatch.setitem(sys.modules, "gdown", SimpleNamespace(download=download))
    config = PropertyPredictionConfig(model="enformer", task="enhancers", checkpoint_path=path)
    info = ensure_prediction_checkpoint(config, verbose=False)
    assert info.paths == (str(path),) and info.downloaded
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize(
    "model,task", [("segmentnt", "genomic_segmentation"), ("evo2", "exon_probability")]
)
def test_direct_heads_download_then_report_local_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any, model: str, task: str
) -> None:
    import dnakit.predictions.checkpoints as checkpoints

    calls: list[str] = []

    def snapshot(repo_id: str, target: Path, **kwargs: Any) -> None:
        calls.append(repo_id)
        target.mkdir(parents=True, exist_ok=True)
        (target / "config.json").write_text("{}")
        (target / "modeling_segment_nt.py").write_text("# fixture")
        (target / "model.safetensors").write_bytes(b"weights")
        (target / "evo2_7b_base.pt").write_bytes(b"weights")

    monkeypatch.setattr(checkpoints, "_snapshot", snapshot)
    config = PropertyPredictionConfig(model=model, task=task, checkpoint_dir=tmp_path)
    first = ensure_prediction_checkpoint(config)
    assert first.downloaded
    messages = capsys.readouterr().err
    assert "Starting checkpoint download" in messages
    assert "Checkpoint download complete" in messages
    second = ensure_prediction_checkpoint(config)
    assert not second.downloaded and first.paths == second.paths
    messages = capsys.readouterr().err
    assert "Starting checkpoint download" not in messages
    assert all(path in messages for path in second.paths)
    assert len(calls) == (1 if model == "segmentnt" else 2)


def test_sharded_backbone_downloads_missing_shard_and_reuses(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any
) -> None:
    import dnakit.representations.checkpoints as checkpoints

    target = tmp_path / "dnabert2-117m"
    target.mkdir()
    (target / "config.json").write_text("{}")
    (target / "model.safetensors.index.json").write_text(
        json.dumps({"weight_map": {"a": "one.safetensors", "b": "two.safetensors"}})
    )
    (target / "one.safetensors").write_bytes(b"weights")
    calls: list[str] = []

    def download(spec: Any, destination: Path, **kwargs: Any) -> None:
        calls.append(spec.name)
        (destination / "two.safetensors").write_bytes(b"weights")

    monkeypatch.setattr(checkpoints, "_download_huggingface", download)
    first = ensure_model_checkpoint("dnabert2", checkpoint_dir=tmp_path)
    assert first.downloaded
    messages = capsys.readouterr().err
    assert "Starting checkpoint download" in messages
    assert "Checkpoint download complete" in messages
    second = ensure_model_checkpoint("dnabert2", checkpoint_dir=tmp_path)
    assert not second.downloaded
    assert calls == ["dnabert2"]
    assert str(target) in capsys.readouterr().err


def test_prediction_reports_checkpoint_before_backend_load(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any
) -> None:
    import dnakit.predictions.api as api

    path = tmp_path / "enhancers.ckpt"
    path.write_bytes(b"existing")

    class Backend:
        def predict(self, inputs: Any, *, show_progress: bool) -> Any:
            return [PredictionOutput((0.2, 0.8), ("0", "1")) for _ in inputs]

    def create_backend(checkpoints: Any, config: Any) -> Backend:
        assert checkpoints.paths == (str(path),)
        messages = capsys.readouterr().err
        assert "Using local checkpoint" in messages
        assert str(path) in messages
        return Backend()

    monkeypatch.setattr(api, "create_prediction_backend", create_backend)
    result = predict_enformer_benchmark(
        [BiologicalSequence("sequence-1", "ACGT" * 50)],
        task="enhancers",
        checkpoint_dir=tmp_path,
    )
    assert len(result.records) == 1
