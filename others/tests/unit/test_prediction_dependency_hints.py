"""Missing prediction dependencies produce guidance before the traceback."""

from __future__ import annotations

import builtins
import importlib
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import BackendUnavailableError
from dnakit.predictions import backends, checkpoints, get_enformer_benchmark_task
from dnakit.predictions.custom_training import encoder, ml_models, splitting, targets, training
from dnakit.representations import backends as representation_backends
from dnakit.representations import get_embedding_model


def _block_module(monkeypatch: pytest.MonkeyPatch, missing: str) -> None:
    original = importlib.import_module

    def load(name: str, package: str | None = None) -> Any:
        if name == missing:
            raise ModuleNotFoundError(f"No module named {name!r}", name=name)
        return original(name, package)

    monkeypatch.setattr(importlib, "import_module", load)


@pytest.mark.parametrize(
    ("model", "module", "extra"),
    [
        ("Enformer benchmark", "torch", "neural,neural-enformer"),
        ("Enformer benchmark", "enformer_pytorch", "neural,neural-enformer"),
        ("Enformer", "enformer_pytorch", "neural,neural-enformer"),
        ("SegmentNT", "transformers", "neural"),
        ("GENERator", "torch", "neural"),
        ("Evo 2", "evo2", "neural,neural-evo2"),
        ("Evo 2 exon classifier", "torch", "neural,neural-evo2"),
    ],
)
def test_neural_dependency_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    model: str,
    module: str,
    extra: str,
) -> None:
    _block_module(monkeypatch, module)
    with pytest.raises(BackendUnavailableError) as error:
        backends._require_module(module, model)
    captured = capsys.readouterr()
    command = f'python -m pip install -e ".[{extra}]"'
    assert command in captured.err
    assert command in str(error.value)
    assert f'python -m pip install "dnakit[{extra}]"' in captured.err
    assert error.value.code == "MISSING_NEURAL_DEPENDENCY"
    assert isinstance(error.value.__cause__, ModuleNotFoundError)
    assert captured.out == ""


@pytest.mark.parametrize(
    "module", ["sklearn.ensemble", "numpy", "joblib", "xgboost", "lightgbm", "catboost"]
)
def test_ml_guidance_does_not_require_neural_packages(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    module: str,
) -> None:
    _block_module(monkeypatch, module)
    with pytest.raises(BackendUnavailableError) as error:
        ml_models.require_ml_module(module)
    extra = module if module in {"xgboost", "lightgbm", "catboost"} else "ml"
    message = capsys.readouterr().err
    assert f'python -m pip install -e ".[{extra}]"' in message
    assert "neural" not in message
    assert error.value.code == "MISSING_MACHINE_LEARNING_DEPENDENCY"


@pytest.mark.parametrize("model", ["enformer", "caduceus", "dnabert2"])
def test_ml_embedding_backend_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    model: str,
) -> None:
    _block_module(monkeypatch, "torch")
    spec = get_embedding_model(model)
    with pytest.raises(BackendUnavailableError):
        representation_backends._require_module("torch", spec)
    extra = "neural" if model == "dnabert2" else f"neural,neural-{model}"
    assert f'python -m pip install -e ".[{extra}]"' in capsys.readouterr().err


def test_checkpoint_download_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    _block_module(monkeypatch, "gdown")
    spec = get_enformer_benchmark_task("enhancers")
    with pytest.raises(BackendUnavailableError):
        checkpoints._download_benchmark_checkpoint(tmp_path / "enhancers.ckpt", spec)
    assert 'python -m pip install -e ".[neural,neural-enformer]"' in capsys.readouterr().err
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("model", ["AlphaGenome", "LucaOneTasks", "JanusDNA"])
def test_official_source_backends_keep_source_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    model: str,
) -> None:
    _block_module(monkeypatch, "torch")
    with pytest.raises(BackendUnavailableError):
        backends._require_module("torch", model)
    message = capsys.readouterr().err
    assert "https://github.com/" in message
    assert "model_source_path" in message
    assert "neural-enformer" not in message


@pytest.mark.parametrize("module", ["torch", "transformers"])
def test_finetuning_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    module: str,
) -> None:
    _block_module(monkeypatch, module)
    with pytest.raises(BackendUnavailableError):
        encoder._require(module)
    assert 'python -m pip install -e ".[training]"' in capsys.readouterr().err


@pytest.mark.parametrize(
    "reader", [targets._numpy, splitting._numpy, training._numpy, training._torch]
)
def test_direct_import_checks_print_guidance(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    reader: Any,
) -> None:
    original = builtins.__import__
    missing = "torch" if reader is training._torch else "numpy"

    def load(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == missing:
            raise ModuleNotFoundError(f"No module named {name!r}", name=name)
        return original(name, *args, **kwargs)

    with monkeypatch.context() as context:
        context.setattr(builtins, "__import__", load)
        with pytest.raises(BackendUnavailableError):
            reader()
    extra = "training" if reader in (training._numpy, training._torch) else "ml"
    assert f'python -m pip install -e ".[{extra}]"' in capsys.readouterr().err


def test_installed_module_is_reused_silently(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    existing = SimpleNamespace(__version__="2.7.1")
    monkeypatch.setitem(sys.modules, "torch", existing)
    assert backends._require_module("torch", "Enformer benchmark", "torch>=2.5") is existing
    assert encoder._require("torch") is existing
    assert capsys.readouterr() == ("", "")


def test_hint_precedes_real_traceback_even_in_quiet_mode() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            """
import importlib
from dnakit.predictions.backends import _require_module
from dnakit._verbosity import with_verbose
def missing(name):
    raise ModuleNotFoundError("No module named 'torch'", name="torch")
importlib.import_module = missing
@with_verbose
def run(*, verbose=True):
    _require_module("torch", "Enformer benchmark", "torch>=2.5")
run(verbose=False)
""",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    command = 'python -m pip install -e ".[neural,neural-enformer]"'
    assert result.stderr.index(command) < result.stderr.index("Traceback")
    assert "MISSING_NEURAL_DEPENDENCY" in result.stderr
    assert "--upgrade" not in result.stderr
    assert "--force-reinstall" not in result.stderr
