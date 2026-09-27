"""Native fine-tuning adapters must update backbones and reload independently."""

from __future__ import annotations

import json
import pickle
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.predictions.custom_training import encoder_native as native
from dnakit.representations.models import MODEL_REGISTRY

torch = pytest.importorskip("torch")


def _write_untrusted_marker(path: str) -> None:
    Path(path).write_text("untrusted checkpoint executed", encoding="utf-8")


class _UntrustedCheckpointMetadata:
    def __init__(self, marker: Path) -> None:
        self.marker = marker

    def __reduce__(self) -> Any:
        return _write_untrusted_marker, (str(self.marker),)


class _JanusConfig:
    def __init__(self, **kwargs: Any) -> None:
        self.values = kwargs

    def to_dict(self) -> dict[str, Any]:
        return self.values


class _TinyJanus(torch.nn.Module):
    def __init__(self, config: Any) -> None:
        super().__init__()
        self.embedding = torch.nn.Embedding(16, 4)
        self.layer = torch.nn.Linear(4, 4)

    def forward(self, input_ids: Any, **kwargs: Any) -> Any:
        # Sequence-wide context exposes padding-induced drift in a mixed batch.
        hidden = self.layer(self.embedding(input_ids))
        return SimpleNamespace(last_hidden_state=hidden + hidden.mean(dim=1, keepdim=True))


class _TinyJanusLM(torch.nn.Module):
    def __init__(self, config: Any) -> None:
        super().__init__()
        self.model = _TinyJanus(config)
        self.lm_head = torch.nn.Linear(4, 16)


class _TinyEvo(torch.nn.Module):
    def __init__(self, config: Any) -> None:
        super().__init__()
        self.config = config
        self.embedding_layer = torch.nn.Embedding(128, 4)
        self.blocks = torch.nn.ModuleList([torch.nn.Linear(4, 4), torch.nn.Linear(4, 4)])
        self.norm = torch.nn.LayerNorm(4)
        self.block_idx_to_device = {0: "cpu", 1: "cpu"}
        self.inference_mode = True
        self.log_poles = torch.nn.Parameter(torch.ones(1))
        self.residues = torch.nn.Parameter(torch.ones(1))
        self.h = torch.ones(1)
        self.t = torch.arange(3)

    def stateless_forward(self, hidden: Any) -> tuple[Any, None]:
        assert self.h is None
        hidden = hidden * self.log_poles + self.residues
        for block in self.blocks:
            hidden = torch.tanh(block(hidden))
        return hidden, None

    def forward(self, *args: Any, **kwargs: Any) -> Any:
        raise AssertionError("The inference wrapper/unembedding must not be called")


@pytest.fixture
def native_modules(monkeypatch: pytest.MonkeyPatch) -> None:
    def original_load(model: Any, path: str) -> None:
        # Simulate the original loader's attention conversion exactly once.
        state = native._load_weights(torch, Path(path))
        state["blocks.0.weight"] = state["blocks.0.weight"].flip(0)
        model.load_state_dict(state)
        with torch.inference_mode():
            model.blocks[0].weight = torch.nn.Parameter(model.blocks[0].weight.clone())

    class Tokenizer:
        def __init__(self, size: int) -> None:
            pass

        def tokenize(self, sequence: str) -> list[int]:
            return list(sequence.encode("ascii"))

    modules = {
        "torch": torch,
        "janusdna.configuration_janusdna": SimpleNamespace(JanusDNAConfig=_JanusConfig),
        "janusdna.modeling_janusdna": SimpleNamespace(
            JanusDNAModel=_TinyJanus,
            JanusDNAForCausalLM=_TinyJanusLM,
        ),
        "vortex.model.model": SimpleNamespace(StripedHyena=_TinyEvo),
        "vortex.model.utils": SimpleNamespace(dotdict=dict, load_checkpoint=original_load),
        "vortex.model.tokenizer": SimpleNamespace(CharLevelTokenizer=Tokenizer),
        "evo2.utils": SimpleNamespace(CONFIG_MAP={"evo2_7b": "config.yml"}),
        "yaml": SimpleNamespace(safe_load=lambda raw: {"vocab_size": 128}),
    }
    monkeypatch.setattr(native, "_require", modules.__getitem__)
    monkeypatch.setattr(native.pkgutil, "get_data", lambda *args: b"vocab_size: 128")
    # Lightweight substitutes do not need Vortex's CUDA-only execution kernels.
    monkeypatch.setattr(native, "_require_evo2_cuda", lambda device: None)


def _source(directory: Path, name: str) -> Path:
    directory.mkdir()
    if name == "janusdna":
        (directory / "model_config.json").write_text('{"config": {}}', encoding="utf-8")
        (directory / "checkpoints").mkdir()
        model = _TinyJanusLM(None)
        state = {f"model.{key}": value for key, value in model.state_dict().items()}
        torch.save({"state_dict": state}, directory / "checkpoints" / "last.ckpt")
    else:
        torch.save(_TinyEvo({}).state_dict(), directory / "evo2_7b.pt")
    return directory


def _encoder(name: str, path: Path, pooling: str = "mean", dtype: str = "float32") -> Any:
    cls = native.JanusDNASequenceEncoder if name == "janusdna" else native.Evo2SequenceEncoder
    return cls(
        MODEL_REGISTRY[name],
        path,
        device="cpu",
        dtype=dtype,
        pooling=pooling,
        max_length=8,
        allow_remote_code=False,
    )


@pytest.mark.parametrize("name", ["janusdna", "evo2"])
@pytest.mark.parametrize("pooling", ["mean", "max", "last", "cls"])
def test_native_gradient_and_independent_reload(
    native_modules: None,
    tmp_path: Path,
    name: str,
    pooling: str,
) -> None:
    torch.manual_seed(7)
    source = _source(tmp_path / "source", name)
    encoder = _encoder(name, source, pooling)
    encoder.train()
    sequences = ["A", "ACGT", "TCGATCGA"]
    initial = encoder.encode(sequences)
    initial.square().mean().backward()
    gradients = [p.grad for p in encoder.parameters()]
    assert all(g is not None and torch.isfinite(g).all() for g in gradients)
    assert any(g.abs().sum() > 0 for g in gradients)
    torch.optim.SGD(encoder.parameters(), lr=0.05).step()
    encoder.eval()
    expected = encoder.encode(sequences)
    assert not torch.equal(initial, expected)
    torch.testing.assert_close(expected, torch.cat([encoder.encode([s]) for s in sequences]))
    saved = tmp_path / "saved"
    encoder.save_pretrained(saved)
    shutil.rmtree(source)
    loaded = _encoder(name, saved, pooling)
    loaded.eval()
    torch.testing.assert_close(expected, loaded.encode(sequences), rtol=0, atol=0)
    model_id = id(encoder.model)
    with torch.no_grad():
        for parameter in encoder.parameters():
            parameter.add_(0.1)
    encoder.load_pretrained(saved)
    assert id(encoder.model) == model_id
    torch.testing.assert_close(expected, encoder.encode(sequences), rtol=0, atol=0)


@pytest.mark.parametrize("name", ["janusdna", "evo2"])
def test_native_length_is_bases_and_empty_input_rejected(
    native_modules: None,
    tmp_path: Path,
    name: str,
) -> None:
    encoder = _encoder(name, _source(tmp_path / "source", name))
    encoder.eval()
    assert encoder.max_length == 8
    torch.testing.assert_close(encoder.encode(["ACGTACGTAAAA"]), encoder.encode(["ACGTACGT"]))
    torch.testing.assert_close(encoder.encode(["ARYKGT"]), encoder.encode(["ANNNGT"]))
    for sequences in ([], [""]):
        with pytest.raises(ConfigurationError, match="non-empty"):
            encoder.encode(sequences)


def test_evo_training_clears_long_filter_cache_and_uses_final_blocks(
    native_modules: None,
    tmp_path: Path,
) -> None:
    encoder = _encoder("evo2", _source(tmp_path / "source", "evo2"))
    encoder.model.h = torch.ones(3)
    encoder.model.t = torch.ones(3)
    encoder.model.inference_mode = True
    encoder.train()
    assert encoder.model.h is None
    assert encoder.model.t is None
    assert encoder.model.inference_mode is False
    encoder.encode(["ACGT"]).square().sum().backward()
    assert encoder.model.blocks[-1].weight.grad.abs().sum() > 0
    assert encoder.model.norm.weight.grad.abs().sum() > 0
    encoder.save_pretrained(tmp_path / "saved")
    config = json.loads((tmp_path / "saved" / "dnakit_evo2_config.json").read_text())
    assert config["use_fp8_input_projections"] is False
    assert config["inference_mode"] is False


def test_saved_native_precision_preserved_with_auto(
    native_modules: None,
    tmp_path: Path,
) -> None:
    encoder = _encoder("evo2", _source(tmp_path / "source", "evo2"))
    encoder.model.to(dtype=torch.float64)
    encoder.eval()
    expected = encoder.encode(["ACGT"])
    encoder.save_pretrained(tmp_path / "saved")
    loaded = _encoder("evo2", tmp_path / "saved", dtype="auto")
    loaded.eval()
    assert next(loaded.parameters()).dtype == torch.float64
    torch.testing.assert_close(expected, loaded.encode(["ACGT"]), rtol=0, atol=0)


@pytest.mark.parametrize("device", ["cpu", "mps"])
def test_evo2_official_runtime_requires_cuda(device: str) -> None:
    with pytest.raises(ConfigurationError, match="requires CUDA"):
        native._require_evo2_cuda(torch.device(device))


def test_janus_frozen_embeddings_do_not_depend_on_chunk_batch_size() -> None:
    from dnakit.representations.backends import _JanusDNABackend

    backend = _JanusDNABackend.__new__(_JanusDNABackend)
    backend.torch = torch
    backend.device = torch.device("cpu")
    backend.model = _TinyJanus(None).eval()
    backend.base_ids = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
    backend.chunk_length = 4
    values = []
    for batch_size in (1, 2, 3):
        backend.config = SimpleNamespace(batch_size=batch_size, pooling="mean")
        values.append(backend._encode_sequence("ACGTTGCAA"))
    torch.testing.assert_close(values[0], values[1])
    torch.testing.assert_close(values[0], values[2])


def test_official_janus_metadata_uses_scoped_fixed_allowlist(tmp_path: Path) -> None:
    from dnakit.representations.backends import _load_janusdna_checkpoint

    omegaconf = pytest.importorskip("omegaconf")
    checkpoint = tmp_path / "last.ckpt"
    torch.save(
        {
            "state_dict": {"model.weight": torch.arange(4, dtype=torch.float32)},
            "hyper_parameters": omegaconf.OmegaConf.create(
                {"layers": [1, 2], "options": {"hidden_size": 4}}
            ),
        },
        checkpoint,
    )
    with torch.serialization.safe_globals([int]):
        before = list(torch.serialization.get_safe_globals())
        state = _load_janusdna_checkpoint(torch, checkpoint, MODEL_REGISTRY["janusdna"])
        assert torch.serialization.get_safe_globals() == before
    torch.testing.assert_close(state["state_dict"]["model.weight"], torch.arange(4).float())
    assert list(state["hyper_parameters"].layers) == [1, 2]


def test_janus_checkpoint_rejects_unlisted_executable_globals(tmp_path: Path) -> None:
    from dnakit.representations.backends import _load_janusdna_checkpoint

    pytest.importorskip("omegaconf")
    marker = tmp_path / "executed.txt"
    checkpoint = tmp_path / "last.ckpt"
    torch.save({"metadata": _UntrustedCheckpointMetadata(marker)}, checkpoint)
    before = list(torch.serialization.get_safe_globals())
    with pytest.raises(pickle.UnpicklingError, match="Unsupported global"):
        _load_janusdna_checkpoint(torch, checkpoint, MODEL_REGISTRY["janusdna"])
    assert not marker.exists()
    assert torch.serialization.get_safe_globals() == before
