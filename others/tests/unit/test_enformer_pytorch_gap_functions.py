"""Tests for dnakit.predictions.enformer_pytorch_gap_functions.

A tiny randomly initialised Enformer stands in for the published weights; the download test patches
``Enformer.from_pretrained`` so that nothing is fetched from the Hugging Face hub.
"""

from __future__ import annotations

import inspect
import json

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

torch = pytest.importorskip("torch")
numpy = pytest.importorskip("numpy")
enformer_pytorch = pytest.importorskip("enformer_pytorch")

from enformer_pytorch import Enformer, EnformerConfig

from dnakit.predictions import enformer_pytorch_gap_functions as facades

TINY = {
    "dim": 96,
    "depth": 1,
    "heads": 2,
    "num_downsamples": 5,
    "dim_divisible_by": 8,
    "target_length": 8,
    "output_heads": {"human": 2},
}


@pytest.fixture()
def tiny_weights(tmp_path):
    torch.manual_seed(0)
    model = Enformer(EnformerConfig(**TINY))
    path = tmp_path / "tiny.pt"
    torch.save(model.state_dict(), path)
    return str(path)


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_download_saves_state_and_config_without_network(tmp_path, monkeypatch):
    torch.manual_seed(1)
    tiny_model = Enformer(EnformerConfig(**TINY))
    requested = {}

    def fake_from_pretrained(name, **kwargs):
        requested["name"] = name
        requested["kwargs"] = kwargs
        return tiny_model

    monkeypatch.setattr(Enformer, "from_pretrained", fake_from_pretrained)
    out = tmp_path / "weights.pt"
    result = facades.download_pretrained_enformer(
        "EleutherAI/enformer-official-rough", str(out), verbose=False
    )
    assert requested["name"] == "EleutherAI/enformer-official-rough"
    assert requested["kwargs"] == {}
    attention = [m for m in tiny_model.modules() if type(m).__name__ == "Attention"]
    assert attention and all(m.use_tf_gamma for m in attention)
    assert out.is_file()
    config = json.loads((tmp_path / "weights.pt.config.json").read_text(encoding="utf-8"))
    assert config["dim"] == TINY["dim"]
    assert result["parameter_count"] == sum(p.numel() for p in tiny_model.parameters())
    assert result["output_heads"] == ["human"]


def test_download_refuses_existing_outputs(tmp_path):
    existing = tmp_path / "weights.pt"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        facades.download_pretrained_enformer("EleutherAI/x", str(existing), verbose=False)


def test_embeddings_have_expected_shape_and_are_saved(tiny_weights, tmp_path):
    sequences = ["ACGT" * 64, "tgca" * 64]
    out = tmp_path / "emb.npy"
    result = facades.enformer_embeddings(
        tiny_weights, sequences, str(out), config=TINY, verbose=False
    )
    array = numpy.load(out)
    assert array.shape == (2, TINY["target_length"], TINY["dim"] * 2)
    assert result["shape"] == list(array.shape)
    assert result["sequence_length"] == 256


def test_frozen_embeddings_match_trainable_forward(tiny_weights, tmp_path):
    sequences = ["ACGT" * 64]
    trainable = tmp_path / "trainable.npy"
    frozen = tmp_path / "frozen.npy"
    facades.enformer_embeddings(tiny_weights, sequences, str(trainable), config=TINY, verbose=False)
    result = facades.enformer_embeddings(
        tiny_weights, sequences, str(frozen), config=TINY, freeze=True, verbose=False
    )
    assert result["freeze"] is True
    assert numpy.allclose(numpy.load(trainable), numpy.load(frozen))


def test_embeddings_reject_bad_inputs(tiny_weights, tmp_path):
    out = tmp_path / "bad.npy"
    with pytest.raises(ConfigurationError):
        facades.enformer_embeddings(
            tiny_weights, ["ACGT", "ACG"], str(out), config=TINY, verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.enformer_embeddings(tiny_weights, ["ACGX"], str(out), config=TINY, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.enformer_embeddings(
            tiny_weights,
            ["ACGT" * 64],
            str(out),
            config=TINY,
            freeze=True,
            train_last_n_layers_only=1,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        facades.enformer_embeddings(
            tiny_weights,
            ["ACGT" * 64],
            str(out),
            config={"dim": 5, "unknown_option": 1},
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        facades.enformer_embeddings(
            str(tmp_path / "missing.pt"), ["ACGT" * 64], str(out), config=TINY, verbose=False
        )


def test_missing_torch_gives_install_hint(monkeypatch, tiny_weights, tmp_path):
    import importlib

    real = importlib.import_module

    def blocked(name, package=None):
        if name == "torch":
            raise ImportError("blocked for test")
        return real(name, package)

    monkeypatch.setattr("dnakit.bioinformatics._common.importlib.import_module", blocked)
    with pytest.raises(BackendUnavailableError) as info:
        facades.enformer_embeddings(
            tiny_weights, ["ACGT" * 64], str(tmp_path / "x.npy"), config=TINY, verbose=False
        )
    assert "enformer" in (info.value.hint or "")
