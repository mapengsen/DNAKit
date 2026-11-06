"""Offline tests for dnakit.predictions.kipoi_gap_functions (kipoi 0.8.6).

Each test writes a toy kipoi model into pytest's tmp_path and reads it with source="dir" (a local directory), so
no test touches the network. Socket access is blocked for every test in this module.
"""

from __future__ import annotations

import inspect
import json
import socket
import sys
from pathlib import Path

import numpy as np
import pytest

import dnakit.predictions.kipoi_gap_functions as kg
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

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

DATALOADER_PY = """import numpy as np
from kipoi.data import Dataset


class SeqDataset(Dataset):
    def __init__(self, fasta_file):
        self.names, self.seqs = [], []
        for line in open(fasta_file):
            line = line.strip()
            if line.startswith(">"):
                self.names.append(line[1:].split()[0])
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

DATALOADER_REMOTE_YAML = DATALOADER_YAML.replace(
    "    example: example.fa\n",
    "    example:\n"
    "      url: https://example.invalid/data/example.fa\n"
    "      md5: 0123456789abcdef0123456789abcdef\n"
    "      name: example_remote.fa\n",
)

EXAMPLE_FA = ">a\nACGTACGTACGTACGTACGT\n>b\nAAAAAAAAAAAAAAAAAAAA\n>c\nGGGGGGGGGGGGGGGGGGGG\n"
# GC fraction of each example record: 10 of 20 positions, 0 of 20, 20 of 20.
EXPECTED_PREDICTIONS = np.array([0.5, 0.0, 1.0])


class _NetworkAccessError(RuntimeError):
    pass


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def blocked(*args: object, **kwargs: object) -> None:
        raise _NetworkAccessError("network access is not allowed in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", blocked)
    monkeypatch.setattr(socket.socket, "connect", blocked)


@pytest.fixture
def kipoi() -> object:
    return pytest.importorskip("kipoi")


def _write_model(
    directory: Path,
    *,
    remote_example: bool = False,
    model_yaml: str = MODEL_YAML,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "model.yaml").write_text(model_yaml)
    (directory / "model.py").write_text(MODEL_PY)
    (directory / "dataloader.py").write_text(DATALOADER_PY)
    (directory / "dataloader.yaml").write_text(
        DATALOADER_REMOTE_YAML if remote_example else DATALOADER_YAML
    )
    (directory / "example.fa").write_text(EXAMPLE_FA)
    return directory


@pytest.fixture
def zoo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A local model zoo with one toy model; the working directory is tmp_path."""
    root = tmp_path / "zoo"
    _write_model(root / "toy")
    monkeypatch.chdir(tmp_path)
    return root


def _example_inputs(zoo: Path) -> np.ndarray:
    dataloader_cls = kg.load_kipoi_dataloader("zoo/toy/.", source="dir")
    batches = list(
        dataloader_cls(fasta_file=str(zoo / "toy" / "example.fa")).batch_iter(batch_size=3)
    )
    return np.concatenate([batch["inputs"] for batch in batches])


def _write_expected(path: Path, inputs: np.ndarray, preds: np.ndarray) -> None:
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "w") as handle:
        handle.create_dataset("inputs", data=inputs)
        handle.create_dataset("preds", data=preds)


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(kg.__all__) == 10
    for name in kg.__all__:
        function = getattr(kg, name)
        assert callable(function)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_missing_kipoi_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "kipoi", None)
    with pytest.raises(BackendUnavailableError) as info:
        kg.load_kipoi_model("zoo/toy", source="dir")
    assert info.value.hint == "Install kipoi with: python -m pip install kipoi"


def test_load_model_from_local_directory(kipoi: object, zoo: Path) -> None:
    model = kg.load_kipoi_model("zoo/toy", source="dir", verbose=False)
    assert model.type == "custom"
    batch = np.zeros((2, 20, 4), dtype="float32")
    batch[0, :, 0] = 1.0  # all A: GC fraction 0
    batch[1, :, 2] = 1.0  # all G: GC fraction 1
    assert model.predict_on_batch(batch) == pytest.approx([0.0, 1.0])


def test_load_model_accepts_absolute_directory(kipoi: object, zoo: Path) -> None:
    model = kg.load_kipoi_model(str((zoo / "toy").resolve()), source="dir")
    assert model.type == "custom"


def test_load_model_arguments_are_keyword_only(kipoi: object, zoo: Path) -> None:
    with pytest.raises(TypeError):
        kg.load_kipoi_model("zoo/toy", "dir")  # type: ignore[misc]


def test_unknown_local_model_raises_backend_error(kipoi: object, zoo: Path) -> None:
    with pytest.raises(BackendExecutionError, match="kipoi get_model failed"):
        kg.load_kipoi_model("zoo/missing", source="dir")


def test_load_dataloader_class_batches_fasta(kipoi: object, zoo: Path) -> None:
    dataloader_cls = kg.load_kipoi_dataloader("zoo/toy/.", source="dir")
    batches = list(
        dataloader_cls(fasta_file=str(zoo / "toy" / "example.fa")).batch_iter(batch_size=2)
    )
    assert [batch["inputs"].shape for batch in batches] == [(2, 20, 4), (1, 20, 4)]
    names = [name for batch in batches for name in batch["metadata"]["name"]]
    assert names == ["a", "b", "c"]


def test_describe_dataloader_reports_arguments_as_json(kipoi: object, zoo: Path) -> None:
    result = kg.describe_kipoi_dataloader("zoo/toy/.", source="dir")
    json.dumps(result)
    description = result["description"]
    assert description["defined_as"] == "dataloader.py::SeqDataset"
    assert description["type"] == "Dataset"
    assert description["inputs"] == {"name": "seq", "shape": [None, 4], "doc": "one-hot"}
    argument = description["arguments"]["fasta_file"]
    assert argument["doc"] == "FASTA file"
    assert argument["example"] == "example.fa"
    assert argument["optional"] is False


def test_describe_dataloader_reports_remote_example_without_download(
    kipoi: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_model(tmp_path / "remote", remote_example=True)
    monkeypatch.chdir(tmp_path)
    result = kg.describe_kipoi_dataloader("remote/.", source="dir")
    assert result["description"]["arguments"]["fasta_file"]["example"] == {
        "url": "https://example.invalid/data/example.fa",
        "md5": "0123456789abcdef0123456789abcdef",
        "name": "example_remote.fa",
    }


def test_describe_component_reads_model_and_dataloader_yaml(kipoi: object, zoo: Path) -> None:
    model = kg.describe_kipoi_component(str(zoo / "toy"), which="model")
    assert model["type"] == "custom"
    assert model["targets"]["name"] == "gc"
    dataloader = kg.describe_kipoi_component(str(zoo / "toy"), which="dataloader")
    assert dataloader["defined_as"] == "dataloader.py::SeqDataset"
    json.dumps([model, dataloader])


def test_describe_component_rejects_bad_arguments(kipoi: object, zoo: Path) -> None:
    with pytest.raises(ConfigurationError):
        kg.describe_kipoi_component(str(zoo / "toy"), which="source")
    with pytest.raises(ConfigurationError):
        kg.describe_kipoi_component(str(zoo / "missing"))


def test_list_local_zoo_models_and_dataloaders(kipoi: object, zoo: Path) -> None:
    _write_model(zoo / "remote", remote_example=True)
    models = kg.list_kipoi_directory_models(str(zoo))
    assert models["total_models"] == 2
    assert sorted(row["model"] for row in models["models"]) == ["remote", "toy"]
    dataloaders = kg.list_kipoi_directory_dataloaders(str(zoo))
    assert sorted(row["dataloader"] for row in dataloaders["dataloaders"]) == ["remote", "toy"]
    json.dumps([models, dataloaders])


def test_list_empty_directory_returns_no_models(kipoi: object, tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    result = kg.list_kipoi_directory_models(str(empty))
    assert result["total_models"] == 0
    assert result["models"] == []


def test_predict_writes_tsv_matching_model_output(kipoi: object, zoo: Path, tmp_path: Path) -> None:
    out = tmp_path / "predictions.tsv"
    result = kg.predict_kipoi_model(
        "zoo/toy",
        str(out),
        source="dir",
        dataloader_args={"fasta_file": str(zoo / "toy" / "example.fa")},
        batch_size=2,
    )
    json.dumps(result)
    assert result["status"] == "completed"
    assert result["samples"] == 3
    assert result["batches"] == 2
    assert result["outputs"] == [str(out.resolve())]
    lines = out.read_text().splitlines()
    assert lines[0] == "preds"
    assert [float(line) for line in lines[1:]] == pytest.approx(EXPECTED_PREDICTIONS.tolist())


def test_predict_refuses_hdf5_existing_and_missing_folder(
    kipoi: object, zoo: Path, tmp_path: Path
) -> None:
    arguments = {"fasta_file": str(zoo / "toy" / "example.fa")}
    with pytest.raises(ConfigurationError, match="TSV"):
        kg.predict_kipoi_model(
            "zoo/toy", str(tmp_path / "out.h5"), source="dir", dataloader_args=arguments
        )
    existing = tmp_path / "exists.tsv"
    existing.write_text("x")
    with pytest.raises(ConfigurationError, match="must not exist"):
        kg.predict_kipoi_model("zoo/toy", str(existing), source="dir", dataloader_args=arguments)
    with pytest.raises(ConfigurationError, match="does not exist"):
        kg.predict_kipoi_model(
            "zoo/toy",
            str(tmp_path / "missing" / "out.tsv"),
            source="dir",
            dataloader_args=arguments,
        )


def test_predict_rejects_bad_batch_size_and_arguments(
    kipoi: object, zoo: Path, tmp_path: Path
) -> None:
    out = str(tmp_path / "out.tsv")
    with pytest.raises(ConfigurationError, match="batch_size"):
        kg.predict_kipoi_model("zoo/toy", out, source="dir", batch_size=0)
    with pytest.raises(ConfigurationError, match="dataloader_args"):
        kg.predict_kipoi_model("zoo/toy", out, source="dir", dataloader_args=["fasta_file"])  # type: ignore[arg-type]


def test_validate_example_without_expected_file(kipoi: object, zoo: Path, tmp_path: Path) -> None:
    out = tmp_path / "example_predictions.tsv"
    result = kg.validate_kipoi_model_example("zoo/toy", output_file=str(out))
    json.dumps(result)
    assert result["status"] == "completed"
    assert result["samples"] == 3
    assert result["prediction_shape"] == [3]
    assert result["matched"] is None
    assert result["expect_file"] is None
    assert result["outputs"] == [str(out.resolve())]
    assert out.is_file()


def test_validate_matches_expected_predictions(kipoi: object, zoo: Path, tmp_path: Path) -> None:
    expect = tmp_path / "expect.h5"
    _write_expected(expect, _example_inputs(zoo), EXPECTED_PREDICTIONS.astype("float32"))
    result = kg.validate_kipoi_model_example("zoo/toy", expect_file=str(expect))
    assert result["status"] == "completed"
    assert result["matched"] is True
    assert result["compared_samples"] == 3


def test_validate_reports_mismatch(kipoi: object, zoo: Path, tmp_path: Path) -> None:
    expect = tmp_path / "expect.h5"
    wrong = np.array([0.5, 0.0, 0.9], dtype="float32")
    _write_expected(expect, _example_inputs(zoo), wrong)
    result = kg.validate_kipoi_model_example("zoo/toy", expect_file=str(expect))
    assert result["status"] == "failed"
    assert result["matched"] is False


def test_validate_reports_remote_expected_file_without_downloading(
    kipoi: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model_yaml = MODEL_YAML + (
        "test:\n"
        "  expect:\n"
        "    url: https://example.invalid/expect.h5\n"
        "    md5: 0123456789abcdef0123456789abcdef\n"
    )
    _write_model(tmp_path / "remote_test", model_yaml=model_yaml)
    monkeypatch.chdir(tmp_path)
    result = kg.validate_kipoi_model_example("remote_test")
    assert result["status"] == "completed"
    assert result["matched"] is None
    assert "not downloaded" in result["note"]


def test_example_kwargs_return_local_example_as_written(
    kipoi: object, zoo: Path, tmp_path: Path
) -> None:
    result = kg.get_kipoi_example_kwargs(
        "zoo/toy", source="dir", output_dir=str(tmp_path / "examples")
    )
    json.dumps(result)
    assert result["example_kwargs"] == {"fasta_file": "example.fa"}


def test_example_kwargs_remote_dry_run_names_target_without_download(
    kipoi: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_model(tmp_path / "remote", remote_example=True)
    monkeypatch.chdir(tmp_path)
    cache = tmp_path / "cache"
    result = kg.get_kipoi_example_kwargs(
        "remote", source="dir", output_dir=str(cache), dry_run=True
    )
    assert result["dry_run"] is True
    assert result["example_kwargs"]["fasta_file"] == str(cache.absolute() / "example_remote.fa")
    assert not (cache / "example_remote.fa").exists()


def test_compare_predictions_uses_decimal_places(kipoi: object) -> None:
    assert kg.compare_kipoi_predictions([1.0, 2.0], [1.0, 2.0])["matched"] is True
    assert kg.compare_kipoi_predictions([1.0, 2.0], [1.0, 2.0004], decimal=3)["matched"] is True
    different = kg.compare_kipoi_predictions([1.0, 2.0], [1.0, 2.5])
    assert different["matched"] is False
    assert different["mismatch"]
    shorter = kg.compare_kipoi_predictions([1.0, 2.0], [1.0, 2.0, 3.0])
    assert shorter["matched"] is False
    assert shorter["expected_shape"] == [3]
    with pytest.raises(ConfigurationError):
        kg.compare_kipoi_predictions([1.0], [1.0], decimal=-1)
