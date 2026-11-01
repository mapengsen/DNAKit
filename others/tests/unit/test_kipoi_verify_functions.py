"""Tests for dnakit.predictions.kipoi_verify_functions (kipoi 0.8.6, offline).

Every expected value is produced by kipoi's own function on the same input. The models and dataloaders are written
into pytest's temporary directory, so nothing is read from the network or from the repository.
"""

from __future__ import annotations

import importlib
import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("kipoi")

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError  # noqa: E402

kv = importlib.import_module("dnakit.predictions.kipoi_verify_functions")

_MODEL_YAML = """type: custom
args:
  file: model.py
  object: {object}
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

_DATALOADER_PY = '''import numpy as np
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
'''

_DATALOADER_YAML = """type: Dataset
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

_MODEL_PY = '''import numpy as np
from kipoi.model import BaseModel


class GcModel(BaseModel):
    def predict_on_batch(self, x):
        x = np.asarray(x)
        return x[:, :, 1:3].sum(axis=(1, 2)) / x.shape[1]
'''

_BROKEN_MODEL_PY = '''from kipoi.model import BaseModel


class BrokenModel(BaseModel):
    def predict_on_batch(self, x):
        raise RuntimeError("deliberately broken toy model")
'''

_FASTA = ">a\nGGGGCCCCGGGGCCCCGGGG\n>b\nAAAAAAAAAATTTTTTTTTT\n"


def _write_model(root: Path, *, broken: bool = False) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "model.yaml").write_text(_MODEL_YAML.format(object="BrokenModel" if broken else "GcModel"))
    (root / "model.py").write_text(_BROKEN_MODEL_PY if broken else _MODEL_PY)
    (root / "dataloader.yaml").write_text(_DATALOADER_YAML)
    (root / "dataloader.py").write_text(_DATALOADER_PY)
    (root / "example.fa").write_text(_FASTA)
    return root


def _write_dataloader(root: Path, args_yaml: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    head, tail = _DATALOADER_YAML.split("args:\n", 1)
    _, rest = tail.split("info:\n", 1)
    (root / "dataloader.yaml").write_text(head + "args:\n" + args_yaml + "info:\n" + rest)
    (root / "dataloader.py").write_text(_DATALOADER_PY)
    (root / "example.fa").write_text(_FASTA)
    return root


def _batch(n: int) -> dict:
    return {
        "inputs": np.arange(n * 20 * 4, dtype="float32").reshape(n, 20, 4),
        "targets": np.zeros((n, 1), dtype="float32"),
        "metadata": {"name": np.array([f"r{i}" for i in range(n)])},
    }


def _same(a, b) -> bool:
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a), np.asarray(b))
    return a == b


def test_public_functions_take_keyword_only_verbose() -> None:
    assert set(kv.__all__) == {
        "prepare_kipoi_batch",
        "concatenate_kipoi_batches",
        "kipoi_dataloader_example_kwargs",
        "kipoi_dataloader_default_args",
        "test_kipoi_source_models",
    }
    for name in kv.__all__:
        parameter = inspect.signature(getattr(kv, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_prepare_batch_matches_kipoi_and_keeps_input_unchanged() -> None:
    import kipoi.cli.main as main

    predictions = np.array([0.25, 0.5])
    cases = [
        ({"inputs": np.ones((2, 20, 4)), "targets": np.zeros(2), "metadata": {"name": ["a", "b"]}}, {}),
        ({"inputs": np.ones((2, 20, 4)), "targets": np.zeros(2), "metadata": {"name": ["a", "b"]}},
         {"keep_inputs": True, "keep_metadata": True}),
        ({"inputs": np.ones((2, 20, 4))}, {"keep_metadata": True}),  # edge: no targets or metadata entries
    ]
    for batch, flags in cases:
        original = dict(batch)
        got = kv.prepare_kipoi_batch(batch, predictions, **flags)
        expected = main.prepare_batch(dict(original), predictions, **{k: v for k, v in flags.items()})
        assert _same(got, expected)
        assert batch.keys() == original.keys()


def test_prepare_batch_rejects_non_mapping() -> None:
    with pytest.raises(ConfigurationError):
        kv.prepare_kipoi_batch([1, 2, 3], np.zeros(3))  # type: ignore[arg-type]


def test_concatenate_batches_matches_kipoi() -> None:
    import kipoi.cli.main as main

    cases = [
        [_batch(2), _batch(3)],  # two batches of the same structure
        [_batch(4), _batch(1), _batch(2)],  # three batches with uneven sizes
        [_batch(5)],  # edge: a single batch is returned with the same samples
    ]
    for batches in cases:
        got = kv.concatenate_kipoi_batches(batches)
        expected = main.numpy_collate_concat(list(batches))
        assert _same(got, expected)
    joined = kv.concatenate_kipoi_batches(cases[0])
    assert joined["inputs"].shape == (5, 20, 4)
    assert list(joined["metadata"]["name"]) == ["r0", "r1", "r0", "r1", "r2"]


def test_concatenate_batches_errors() -> None:
    with pytest.raises(ConfigurationError):
        kv.concatenate_kipoi_batches([])
    with pytest.raises(ConfigurationError):
        kv.concatenate_kipoi_batches("not a batch list")
    with pytest.raises(BackendExecutionError):
        kv.concatenate_kipoi_batches([{"inputs": np.zeros((2, 20, 4))}, {"inputs": np.zeros((2, 7, 4))}])


def test_example_kwargs_plain_values_match_kipoi(tmp_path: Path) -> None:
    import kipoi.sources as sources
    import kipoi.specs as specs

    toy = _write_model(tmp_path / "toy")
    descr = sources.load_component_descr(str(toy), which="dataloader")
    expected = specs.example_kwargs(descr.args)
    got = kv.kipoi_dataloader_example_kwargs(str(toy))
    assert got["example_kwargs"] == {"fasta_file": "example.fa"}
    assert got["example_kwargs"] == dict(expected)
    assert got["status"] == "completed" and got["dry_run"] is False


def test_example_kwargs_remote_example_dry_run_and_no_example(tmp_path: Path) -> None:
    import kipoi.sources as sources
    import kipoi.specs as specs

    remote_args = (
        "  fasta_file:\n"
        "    doc: FASTA file\n"
        "    example:\n"
        "      url: https://example.invalid/records.fa\n"
        "      md5: abcdef0123456789abcdef0123456789\n"
        "      name: records.fa\n"
        "  window:\n"
        "    doc: no example given\n"
    )
    dataloader = _write_dataloader(tmp_path / "remote", remote_args)
    cache = tmp_path / "cache"

    # Without a cache directory the remote example is reported, not downloaded.
    no_cache = kv.kipoi_dataloader_example_kwargs(str(dataloader))
    assert no_cache["example_kwargs"] == {
        "fasta_file": {"url": "https://example.invalid/records.fa",
                       "md5": "abcdef0123456789abcdef0123456789", "name": "records.fa"}
    }
    # A dry run with a cache directory returns the target path and writes nothing.
    dry = kv.kipoi_dataloader_example_kwargs(str(dataloader), cache_dir=str(cache), dry_run=True)
    target = str(cache.absolute() / "records.fa")
    assert dry["example_kwargs"] == {"fasta_file": target}
    assert not (cache / "records.fa").exists()
    expected = specs.example_kwargs(
        sources.load_component_descr(str(dataloader), which="dataloader").args,
        cache_path=str(cache.absolute()),
        absolute_path=True,
        dry_run=True,
    )
    assert dict(expected) == {"fasta_file": target}


def test_example_kwargs_rejects_missing_directory(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        kv.kipoi_dataloader_example_kwargs(str(tmp_path / "missing"))


def test_default_args_match_kipoi(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import kipoi.sources as sources
    import kipoi.specs as specs

    monkeypatch.chdir(tmp_path)
    (tmp_path / "example.fa").write_text(_FASTA)  # the default file must exist in the current directory
    defaults = (
        "  fasta_file:\n"
        "    doc: FASTA file\n"
        "    example: example.fa\n"
        "  width:\n"
        "    doc: window width\n"
        "    default: 987654\n"
        "  reverse:\n"
        "    doc: reverse complement\n"
        "    default: true\n"
        "  label:\n"
        "    doc: free text\n"
        "    default: GC\n"
        "  sequence_file:\n"
        "    doc: a local file used as default (checked in the current directory)\n"
        "    default: example.fa\n"
    )
    dataloader = _write_dataloader(tmp_path / "defaults", defaults)
    got = kv.kipoi_dataloader_default_args(str(dataloader), str(tmp_path / "downloads"))
    expected = specs.download_default_args(
        sources.load_component_descr(str(dataloader), which="dataloader").args, str(tmp_path / "downloads")
    )
    assert got["overrides"] == dict(expected)
    assert got["overrides"]["width"] == 987654 and got["overrides"]["reverse"] is True
    assert got["overrides"]["label"] == "GC"
    assert got["overrides"]["sequence_file"] == str((tmp_path / "example.fa").absolute())
    assert "fasta_file" not in got["overrides"]  # an example is not a default
    assert not (tmp_path / "downloads").exists()  # nothing remote, so nothing created


def test_default_args_without_defaults_and_bad_input(tmp_path: Path) -> None:
    toy = _write_model(tmp_path / "toy")
    got = kv.kipoi_dataloader_default_args(str(toy), str(tmp_path / "out"))
    assert got["overrides"] == {}
    with pytest.raises(ConfigurationError):
        kv.kipoi_dataloader_default_args(str(tmp_path / "missing"), str(tmp_path / "out"))
    with pytest.raises(ConfigurationError):
        kv.kipoi_dataloader_default_args(str(toy), "")


def test_source_models_pass_on_toy_model(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_model(source / "toy")
    result = kv.test_kipoi_source_models(str(source))
    assert result["status"] == "passed"
    assert result["tested"] == 1 and result["passed"] == 1 and result["failed"] == []
    assert result["results"][0]["model"] == "toy"
    explicit = kv.test_kipoi_source_models(str(source), models=["toy"], batch_size=1)
    assert explicit["status"] == "passed" and explicit["batch_size"] == 1


def test_source_models_collect_failures_and_stop_on_failure(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_model(source / "good")
    _write_model(source / "bad", broken=True)
    result = kv.test_kipoi_source_models(str(source))
    assert result["status"] == "failed"
    assert result["failed"] == ["bad"]
    # models run in sorted order: "bad" before "good"
    assert [(r["model"], r["status"]) for r in result["results"]] == [("bad", "failed"), ("good", "passed")]
    assert "deliberately broken toy model" in result["results"][0]["error"]
    stopped = kv.test_kipoi_source_models(str(source), stop_on_failure=True)
    assert stopped["tested"] == 1 and stopped["failed"] == ["bad"]


def test_source_models_rejects_unknown_model_and_empty_source(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_model(source / "toy")
    with pytest.raises(ConfigurationError):
        kv.test_kipoi_source_models(str(source), models=["nope"])
    with pytest.raises(ConfigurationError):
        kv.test_kipoi_source_models(str(source), models=["../toy"])
    with pytest.raises(ConfigurationError):
        kv.test_kipoi_source_models(str(tmp_path / "empty_source_dir_missing"))
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ConfigurationError):
        kv.test_kipoi_source_models(str(empty))


def test_missing_kipoi_raises_backend_unavailable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    toy = _write_model(tmp_path / "toy")
    monkeypatch.setitem(sys.modules, "kipoi", None)  # makes `import kipoi` raise ImportError
    with pytest.raises(BackendUnavailableError):
        kv.kipoi_dataloader_default_args(str(toy), str(tmp_path / "out"))


def test_default_args_survive_numeric_defaults_that_look_like_open_descriptors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """kipoi 0.8.6 calls os.path.exists(20) for ``default: 20``; an open descriptor 20 then breaks abspath."""
    import kipoi.specs as specs

    monkeypatch.chdir(tmp_path)
    defaults = "  width:\n    doc: window width\n    default: 20\n  ratio:\n    doc: a ratio\n    default: 0.5\n"
    dataloader = _write_dataloader(tmp_path / "defaults", defaults)
    original_exists = specs.os.path.exists
    monkeypatch.setattr(  # an int is read as a file descriptor; an open one makes exists() true
        specs.os.path, "exists", lambda value: True if isinstance(value, int) else original_exists(value)
    )
    got = kv.kipoi_dataloader_default_args(str(dataloader), str(tmp_path / "downloads"))
    assert got["overrides"] == {"width": 20, "ratio": 0.5}
