"""Offline tests for dnakit.predictions.kipoi_precondition_b_functions (kipoi 0.8.6).

Each test writes a toy kipoi model into pytest's tmp_path and uses source="dir", so no test reads the model zoo. Socket
connections are blocked for every test in this module. The Singularity runtime is not installed here: the tests that
run a container put a small shim named ``singularity`` on PATH that drops the ``exec -B <dirs> <image>`` prefix and runs
the kipoi command on the host Python. The shim checks the wiring of the command line only; it is not a container.
"""

from __future__ import annotations

import inspect
import socket
import stat
import sys
from pathlib import Path

import pytest

import dnakit.predictions.kipoi_precondition_b_functions as kb
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

kipoi_gap = pytest.importorskip("dnakit.predictions.kipoi_gap_functions")
pytest.importorskip("kipoi")

MODEL_YAML = """type: custom
args:
  file: model.py
  object: GcModel
default_dataloader: .
info:
  authors:
    - name: DNAKit tests
  doc: Toy model returning the GC fraction of one-hot encoded sequences.
  version: 0.1.0
  license: MIT
dependencies:
  conda:
    - python
  pip:
    - numpy
schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot DNA (A,C,G,T)
  targets:
    name: gc
    shape: (1,)
    doc: GC fraction
"""

MODEL_PY = """import numpy as np
from kipoi.model import BaseModel


class GcModel(BaseModel):
    def predict_on_batch(self, x):
        x = np.asarray(x)
        return x[:, :, 1:3].sum(axis=(1, 2)) / x.shape[1]
"""

DATALOADER_YAML = """type: Dataset
defined_as: dataloader.py::SeqDataset
args:
  fasta_file:
    doc: FASTA file
    example: example.fa
info:
  authors:
    - name: DNAKit tests
  doc: One-hot encodes each FASTA record
output_schema:
  inputs:
    name: seq
    shape: (None, 4)
    doc: one-hot
  metadata:
    name:
      type: str
      doc: record id
dependencies:
  pip:
    - numpy
"""

DATALOADER_PY = """import numpy as np
from kipoi.data import Dataset


class SeqDataset(Dataset):
    def __init__(self, fasta_file):
        self.names, self.seqs = [], []
        name = None
        for line in open(fasta_file):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                self.names.append(name)
                self.seqs.append("")
            elif line:
                self.seqs[-1] += line.upper()

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        s = self.seqs[i][:20].ljust(20, "A")
        x = np.zeros((20, 4), dtype="float32")
        for j, c in enumerate(s):
            if c in "ACGT":
                x[j, "ACGT".index(c)] = 1
        return {"inputs": x, "metadata": {"name": self.names[i]}}
"""

FASTA = ">a\nGGGGCCCCGGGGCCCCGGGG\n>b\nAAAAAAAAAATTTTTTTTTT\n"

PUBLIC = (
    "predict_kipoi_model_in_singularity",
    "test_kipoi_model_in_singularity",
    "test_kipoi_model_in_current_env",
)

SHIM = """#!/bin/sh
if [ "$1" = "--version" ]; then echo "singularity shim 0.0"; exit 0; fi
if [ "$2" = "-B" ]; then shift 4; else shift 2; fi
shift 1
exec {python} -m kipoi "$@"
"""

FAILING_SHIM = """#!/bin/sh
if [ "$1" = "--version" ]; then echo "singularity shim 0.0"; exit 0; fi
echo "simulated container failure" >&2
exit 3
"""


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def blocked(self: socket.socket, *args: object, **kwargs: object) -> None:
        raise OSError("network access is blocked in this test module")

    monkeypatch.setattr(socket.socket, "connect", blocked)


@pytest.fixture()
def toy(tmp_path: Path) -> dict[str, Path]:
    model = tmp_path / "toy_model"
    model.mkdir()
    (model / "model.yaml").write_text(MODEL_YAML, encoding="utf-8")
    (model / "model.py").write_text(MODEL_PY, encoding="utf-8")
    (model / "dataloader.yaml").write_text(DATALOADER_YAML, encoding="utf-8")
    (model / "dataloader.py").write_text(DATALOADER_PY, encoding="utf-8")
    fasta = model / "example.fa"
    fasta.write_text(FASTA, encoding="utf-8")
    image = tmp_path / "kipoi.sif"
    image.write_bytes(b"not a real image; the shim ignores it")
    return {"model": model, "fasta": fasta, "image": image, "root": tmp_path}


def _install_shim(bin_dir: Path, body: str) -> None:
    bin_dir.mkdir(exist_ok=True)
    script = bin_dir / "singularity"
    script.write_text(body, encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture()
def shim_on_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bin_dir = tmp_path / "bin"
    _install_shim(bin_dir, SHIM.format(python=sys.executable))
    monkeypatch.setenv("PATH", f"{bin_dir}")
    return bin_dir


@pytest.fixture()
def no_singularity(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    empty = tmp_path / "empty_path"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))


def test_public_functions_have_keyword_only_verbose() -> None:
    for name in PUBLIC:
        function = getattr(kb, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name
    assert set(kb.__all__) == set(PUBLIC)


def test_missing_singularity_raises_backend_unavailable(
    toy: dict[str, Path], no_singularity: None
) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        kb.predict_kipoi_model_in_singularity(
            str(toy["image"]), str(toy["model"]), str(toy["root"] / "out.tsv"), verbose=False
        )
    assert "singularity" in str(info.value).lower()
    with pytest.raises(BackendUnavailableError):
        kb.test_kipoi_model_in_singularity(str(toy["image"]), str(toy["model"]), verbose=False)


def test_dry_run_returns_the_command_without_singularity(
    toy: dict[str, Path], no_singularity: None
) -> None:
    out = toy["root"] / "dry.tsv"
    result = kb.predict_kipoi_model_in_singularity(
        str(toy["image"]),
        str(toy["model"]),
        str(out),
        dataloader_args={"fasta_file": str(toy["fasta"])},
        batch_size=3,
        dry_run=True,
        verbose=False,
    )
    assert result["status"] == "dry_run"
    assert result["outputs"] == []
    command = result["command"]
    assert command[:2] == ["singularity", "exec"]
    assert str(toy["image"]) in command
    assert command[command.index("kipoi") : command.index("kipoi") + 2] == ["kipoi", "predict"]
    assert command[command.index("--batch_size") + 1] == "3"
    assert str(toy["fasta"].parent) in result["bind_directories"]
    assert not out.exists()

    test_result = kb.test_kipoi_model_in_singularity(
        str(toy["image"]), str(toy["model"]), dry_run=True, verbose=False
    )
    assert test_result["status"] == "dry_run"
    inner = test_result["command"][test_result["command"].index("kipoi") :]
    assert inner[:2] == ["kipoi", "test"]
    assert inner[inner.index("--batch_size") + 1] == "32"
    assert inner[-1] == str(toy["model"])


def test_predict_in_singularity_matches_predict_kipoi_model(
    toy: dict[str, Path], shim_on_path: Path
) -> None:
    direct = toy["root"] / "direct.tsv"
    contained = toy["root"] / "contained.tsv"
    kipoi_gap.predict_kipoi_model(
        str(toy["model"]),
        str(direct),
        source="dir",
        dataloader_args={"fasta_file": str(toy["fasta"])},
        batch_size=2,
        verbose=False,
    )
    result = kb.predict_kipoi_model_in_singularity(
        str(toy["image"]),
        str(toy["model"]),
        str(contained),
        dataloader_args={"fasta_file": str(toy["fasta"])},
        batch_size=2,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["return_code"] == 0
    assert result["outputs"] == [str(contained)]
    assert result["provenance"] == {"backend": "singularity", "version": "singularity shim 0.0"}
    assert contained.read_text(encoding="utf-8") == direct.read_text(encoding="utf-8")
    assert "preds" in contained.read_text(encoding="utf-8").splitlines()[0]


def test_predict_in_singularity_rejects_bad_output_name(
    toy: dict[str, Path], shim_on_path: Path
) -> None:
    with pytest.raises(ConfigurationError):
        kb.predict_kipoi_model_in_singularity(
            str(toy["image"]), str(toy["model"]), str(toy["root"] / "out.h5"), verbose=False
        )
    existing = toy["root"] / "exists.tsv"
    existing.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        kb.predict_kipoi_model_in_singularity(
            str(toy["image"]), str(toy["model"]), str(existing), verbose=False
        )


def test_test_in_singularity_runs_kipoi_test_through_shim(
    toy: dict[str, Path], shim_on_path: Path
) -> None:
    result = kb.test_kipoi_model_in_singularity(
        str(toy["image"]), str(toy["model"]), batch_size=2, verbose=False
    )
    assert result["status"] == "completed"
    assert result["return_code"] == 0
    assert result["outputs"] == []
    assert str(toy["model"]) in result["command"]


def test_test_in_singularity_nonzero_exit_raises(
    toy: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bin_dir = tmp_path / "failing_bin"
    _install_shim(bin_dir, FAILING_SHIM)
    monkeypatch.setenv("PATH", str(bin_dir))
    with pytest.raises(BackendExecutionError) as info:
        kb.test_kipoi_model_in_singularity(str(toy["image"]), str(toy["model"]), verbose=False)
    assert "exited with code 3" in str(info.value)


def test_current_env_test_passes_on_toy_model(toy: dict[str, Path]) -> None:
    result = kb.test_kipoi_model_in_current_env(str(toy["model"]), batch_size=2, verbose=False)
    assert result["status"] == "completed"
    assert result["warnings"] == 0
    assert result["batch_size"] == 2
    assert result["provenance"]["backend"] == "kipoi"


def test_current_env_test_fails_on_logged_warning(
    toy: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    import logging

    import kipoi.cli.main as kmain

    def noisy(command: str, raw_args: list[str]) -> None:
        logging.getLogger("kipoi.test_probe").warning("probe warning")

    monkeypatch.setattr(kmain, "cli_test", noisy)
    with pytest.raises(BackendExecutionError) as info:
        kb.test_kipoi_model_in_current_env(str(toy["model"]), verbose=False)
    assert "1 warning" in str(info.value)
    assert "probe warning" in str(info.value)


def test_current_env_test_reports_kipoi_failures(
    toy: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    import kipoi.cli.main as kmain

    def broken(command: str, raw_args: list[str]) -> None:
        raise ValueError("expected predictions differ")

    monkeypatch.setattr(kmain, "cli_test", broken)
    with pytest.raises(BackendExecutionError, match="expected predictions differ"):
        kb.test_kipoi_model_in_current_env(str(toy["model"]), verbose=False)

    def exits(command: str, raw_args: list[str]) -> None:
        raise SystemExit(2)

    monkeypatch.setattr(kmain, "cli_test", exits)
    with pytest.raises(BackendExecutionError, match="exited with code 2"):
        kb.test_kipoi_model_in_current_env(str(toy["model"]), verbose=False)


def test_current_env_test_rejects_missing_model_dir(toy: dict[str, Path]) -> None:
    with pytest.raises(ConfigurationError):
        kb.test_kipoi_model_in_current_env(str(toy["root"] / "missing"), verbose=False)
