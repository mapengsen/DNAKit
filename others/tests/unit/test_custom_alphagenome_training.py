"""Exercise the actual JAX/Torch gradient bridge with a small Haiku backbone."""

from __future__ import annotations

import shutil
import sys
from contextlib import nullcontext
from enum import Enum
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

from dnakit.predictions.custom_training import encoder_alphagenome as module
from dnakit.representations import RepresentationConfig, get_embedding_model

torch = pytest.importorskip("torch")
jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")
hk = pytest.importorskip("haiku")
ocp = pytest.importorskip("orbax.checkpoint")


@pytest.fixture
def small_backbone(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Any]:
    """Mock only the expensive model architecture, retaining both real runtimes."""

    @hk.transform_with_state
    def forward(dna: Any, organism: Any, is_training: bool) -> Any:
        del organism
        weight = hk.get_parameter("weight", (4, 3), init=hk.initializers.Constant(0.2))
        bias = hk.get_parameter("bias", (3,), init=hk.initializers.Constant(0.1))
        counter = hk.get_state("counter", (), dtype=jnp.float32, init=jnp.zeros)
        if is_training:
            hk.set_state("counter", counter + 1)
        output = jnp.tanh(dna @ weight + bias)
        if is_training:
            output = hk.dropout(hk.next_rng_key(), 0.1, output)
        return output

    cpu = jax.devices("cpu")[0]
    with jax.default_device(cpu):
        params, state = forward.init(
            jax.random.PRNGKey(11), jnp.ones((1, 2048, 4)), jnp.zeros((1,), jnp.int32), False
        )
        params["alphagenome/head/test"] = {"weight": jnp.full((2, 3), 7.0)}
        checkpoint = tmp_path / "original"
        with ocp.StandardCheckpointer() as checkpointer:
            checkpointer.save(str(checkpoint), (params, state))
            checkpointer.wait_until_finished()
    monkeypatch.setattr(module, "_build_training_forward", lambda **kwargs: (forward.apply, 0))
    return checkpoint, forward.apply


def _encoder(checkpoint: Path, *, pooling: str = "mean", seed: int = 42) -> Any:
    return module.AlphaGenomeSequenceEncoder(
        get_embedding_model("alphagenome"),
        checkpoint,
        device="cpu",
        dtype="float32",
        pooling=pooling,
        max_length=8,
        allow_remote_code=False,
        seed=seed,
    )


@pytest.mark.parametrize("pooling", ["mean", "max", "cls", "last"])
def test_real_jax_gradients_and_optimizer_update(small_backbone: Any, pooling: str) -> None:
    checkpoint, apply = small_backbone
    encoder = _encoder(checkpoint, pooling=pooling)
    parameters = tuple(encoder.parameters())
    values = [encoder._to_jax(parameter) for parameter in parameters]
    params = jax.tree_util.tree_unflatten(encoder._tree, values)
    sequences = ["ACGTACGTACGT", "TTT"]

    def reference(current: Any) -> Any:
        from dnakit.representations.backends import _AlphaGenomeBackend

        rows = []
        for sequence in sequences:
            chunks = []
            for start in range(0, len(sequence), 8):
                chunk = sequence[start : start + 8]
                dna = jax.device_put(
                    _AlphaGenomeBackend._one_hot(chunk + "N" * (2048 - len(chunk))),
                    encoder.jax_device,
                )
                output, _ = apply(
                    current, encoder.state, encoder._rng, dna, encoder.organism, False
                )
                if pooling == "mean":
                    chunks.append(output.mean(axis=1))
                elif pooling == "max":
                    chunks.append(output.max(axis=1))
                elif pooling == "cls":
                    chunks.append(output[:, 0, :])
                else:
                    chunks.append(output[:, -1, :])
            rows.append(jnp.stack(chunks).mean(axis=0).squeeze(0))
        return jnp.stack(rows)

    output = encoder.encode(sequences)
    np.testing.assert_allclose(output.detach().numpy(), np.asarray(reference(params)), atol=1e-6)
    expected = jax.grad(lambda current: (reference(current) ** 2).sum())(params)
    output.square().sum().backward()
    for parameter, expected_grad in zip(
        parameters, jax.tree_util.tree_leaves(expected), strict=True
    ):
        assert parameter.requires_grad
        assert parameter.grad is not None
        np.testing.assert_allclose(parameter.grad.numpy(), np.asarray(expected_grad), atol=1e-6)
    before = [parameter.detach().clone() for parameter in parameters]
    optimizer = torch.optim.AdamW(parameters, lr=0.01, weight_decay=0.0)
    optimizer.step()
    assert any(not torch.equal(old, new) for old, new in zip(before, parameters, strict=True))
    assert not torch.equal(output, encoder.encode(sequences))


def test_save_reload_restores_weights_state_and_rng(small_backbone: Any, tmp_path: Path) -> None:
    checkpoint, _ = small_backbone
    encoder = _encoder(checkpoint, seed=7)
    encoder.train()
    before_rng = np.asarray(encoder._rng).copy()
    output = encoder.encode(["ACGT"])
    output.square().sum().backward()
    torch.optim.AdamW(encoder.parameters(), lr=0.01).step()
    np.testing.assert_array_equal(encoder._head_params["alphagenome/head/test"]["weight"], 7.0)
    assert not np.array_equal(before_rng, np.asarray(encoder._rng))
    assert float(jax.tree_util.tree_leaves(encoder.state)[0]) == 1

    saved = tmp_path / "saved"
    encoder.save_pretrained(saved)
    encoder.save_pretrained(saved)  # Replacing the best checkpoint must be supported.
    shutil.rmtree(checkpoint)
    reloaded = _encoder(saved, seed=999)
    encoder.eval()
    with torch.inference_mode():
        torch.testing.assert_close(encoder.encode(["ACGT"]), reloaded.encode(["ACGT"]))
    assert reloaded.seed == 7
    np.testing.assert_array_equal(reloaded._head_params["alphagenome/head/test"]["weight"], 7.0)
    np.testing.assert_array_equal(encoder._rng, reloaded._rng)
    for before, after in zip(
        jax.tree_util.tree_leaves(encoder.state),
        jax.tree_util.tree_leaves(reloaded.state),
        strict=True,
    ):
        np.testing.assert_array_equal(before, after)
    encoder.train()
    reloaded.train()
    torch.testing.assert_close(encoder.encode(["ACGT"]), reloaded.encode(["ACGT"]))


def test_inplace_reload_preserves_parameters(small_backbone: Any, tmp_path: Path) -> None:
    checkpoint, _ = small_backbone
    encoder = _encoder(checkpoint)
    saved = tmp_path / "best"
    encoder.save_pretrained(saved)
    parameters = tuple(encoder.parameters())
    expected = encoder.encode(["AAAA"]).detach().clone()
    expected_rng = np.asarray(encoder._rng).copy()
    with torch.no_grad():
        for parameter in parameters:
            parameter.add_(1)
    encoder.train()
    encoder.encode(["AAAA"])
    encoder.load_pretrained(saved)
    encoder.eval()
    assert all(old is new for old, new in zip(parameters, encoder.parameters(), strict=True))
    with torch.no_grad():
        torch.testing.assert_close(encoder.encode(["AAAA"]), expected)
    np.testing.assert_array_equal(encoder._rng, expected_rng)
    assert float(jax.tree_util.tree_leaves(encoder.state)[0]) == 0


def test_training_is_reproducible_and_eval_has_no_state_updates(small_backbone: Any) -> None:
    checkpoint, _ = small_backbone
    first, second = _encoder(checkpoint, seed=123), _encoder(checkpoint, seed=123)
    first.train()
    second.train()
    torch.testing.assert_close(first.encode(["ACGT"]), second.encode(["ACGT"]))
    first.eval()
    state = first.state
    rng = np.asarray(first._rng).copy()
    with torch.no_grad():
        batched = first.encode(["ACGT", "T"])
        individual = torch.cat([first.encode(["ACGT"]), first.encode(["T"])])
    torch.testing.assert_close(batched, individual)
    assert first.state is state
    np.testing.assert_array_equal(first._rng, rng)


@pytest.mark.parametrize("sequence", ["ARYSWKMBDHVN", "aryswkmbdhvn"])
def test_full_encoder_masks_iupac_symbols_like_n(small_backbone: Any, sequence: str) -> None:
    checkpoint, _ = small_backbone
    encoder = _encoder(checkpoint)
    ambiguous = encoder.encode([sequence])
    masked = encoder.encode(["A" + "N" * 11])
    torch.testing.assert_close(ambiguous, masked, rtol=0, atol=0)
    ambiguous.square().sum().backward()
    assert all(parameter.grad is not None for parameter in encoder.parameters())
    assert any(parameter.grad.abs().sum() > 0 for parameter in encoder.parameters())


def test_alpha_linear_probe_matches_full_encoder(small_backbone: Any, monkeypatch: Any) -> None:
    from dnakit.representations import backends

    checkpoint, forward_apply = small_backbone

    class Organism(Enum):
        HOMO_SAPIENS = 10
        MUS_MUSCULUS = 20

    api = ModuleType("alphagenome.models.dna_model")
    api.Organism = Organism
    research = ModuleType("alphagenome_research.model.dna_model")
    research.convert_to_organism_index = lambda organism: 0 if organism.value == 10 else 1
    metadata = ModuleType("alphagenome_research.model.metadata.metadata")
    metadata.load = lambda organism: {}
    for mock in (api, research, metadata):
        monkeypatch.setitem(sys.modules, mock.__name__, mock)

    def apply(params: Any, state: Any, dna: Any, organism: Any) -> Any:
        embeddings, _ = forward_apply(params, state, None, dna, organism, False)
        return {"embeddings_1bp": embeddings}

    monkeypatch.setattr(backends, "_build_alphagenome_apply_fn", lambda *a, **kw: apply)
    config = RepresentationConfig(model="alphagenome", device="cpu", dtype="float32", max_length=8)
    frozen = backends._AlphaGenomeBackend(get_embedding_model("alphagenome"), checkpoint, config)
    assert frozen.organism_index == 0
    assert len(frozen._pad("ACGT")) == 2048
    encoder = _encoder(checkpoint)
    sequences = ["ACGTACGTACGT", "G"]
    with torch.no_grad():
        expected = encoder.encode(sequences).numpy()
    np.testing.assert_allclose(frozen.extract(sequences, show_progress=False), expected, atol=1e-6)


@pytest.mark.parametrize("count", [3, 5])
def test_official_apply_return_signature_compatibility(count: int) -> None:
    from dnakit.representations.backends import _build_alphagenome_apply_fn

    apply = object()
    result = [None, apply, *([None] * (count - 2))]
    official = SimpleNamespace(create_model=lambda metadata: result)
    assert _build_alphagenome_apply_fn(official, {}, force_float32=False) is apply


def test_alpha_full_training_and_saved_inference(small_backbone: Any, tmp_path: Path) -> None:
    pytest.importorskip("sklearn")
    from dnakit.predictions import (
        CustomPropertyTrainingConfig,
        predict_custom_properties,
        train_custom_property_model,
    )

    checkpoint, _ = small_backbone
    source = tmp_path / "training.csv"
    source.write_text(
        "id,sequence,label\n"
        + "".join(f"s{i},{'G' * (i % 5 + 1)}{'A' * (6 - i % 5)},{i % 2}\n" for i in range(20)),
        encoding="utf-8",
    )
    output = tmp_path / "trained"
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="classification",
            model="alphagenome",
            finetune_method="full",
            checkpoint_path=checkpoint,
            output_dir=output,
            epochs=1,
            batch_size=4,
            workers=0,
            device="cpu",
            dtype="float32",
            max_length=8,
            show_progress=False,
        ),
    )
    assert result.finetune_method == "full"
    assert (output / "backbone" / "alphagenome_encoder.json").is_file()
    shutil.rmtree(checkpoint)
    prediction = predict_custom_properties(
        source, model_dir=output, batch_size=4, device="cpu", show_progress=False
    )
    assert prediction.record_ids == tuple(f"s{i}" for i in range(20))
    assert len(prediction.predictions) == 20


@pytest.mark.parametrize("modern", [False, True], ids=["legacy-call", "modern-trunk"])
def test_official_api_versions_train_and_reload(
    modern: bool, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise the actual builder, including the old inference-only norm API."""

    class LegacyAlphaGenome(hk.Module):
        def __init__(self, metadata: Any) -> None:
            super().__init__(name="alphagenome")
            assert metadata == {}  # Original output heads must not be constructed.

        def __call__(self, dna: Any, organism: Any) -> Any:
            del organism
            w = hk.get_parameter("w", (4, 3), init=hk.initializers.Constant(0.2))
            scale = hk.get_parameter("scale", (3,), init=jnp.ones)
            variance = hk.get_state("variance", (), init=jnp.ones)
            embeddings = SimpleNamespace(embeddings_1bp=(dna @ w + 0.1) * scale / variance)
            return {}, embeddings

    class ModernAlphaGenome(LegacyAlphaGenome):
        @hk.name_like("__call__")
        def forward_trunk(self, dna: Any, organism: Any, *, is_training: bool) -> Any:
            _, embeddings = self(dna, organism)
            updates = hk.get_state("updates", (), init=jnp.zeros)
            if is_training:
                hk.set_state("updates", updates + 1)
            return embeddings

    model_type = ModernAlphaGenome if modern else LegacyAlphaGenome
    modules = {
        "alphagenome_research.model.model": SimpleNamespace(AlphaGenome=model_type),
        "alphagenome_research.model.dna_model": SimpleNamespace(
            convert_to_organism_index=lambda organism: 0
        ),
        "alphagenome.models.dna_model": SimpleNamespace(
            Organism=SimpleNamespace(HOMO_SAPIENS=SimpleNamespace(value=10))
        ),
    }
    require = module._require
    monkeypatch.setattr(
        module, "_require", lambda name: modules[name] if name in modules else require(name)
    )

    @hk.transform_with_state
    def initialize(dna: Any, organism: Any) -> Any:
        model = model_type({})
        if modern:
            return model.forward_trunk(dna, organism, is_training=False).embeddings_1bp
        return model(dna, organism)[1].embeddings_1bp

    checkpoint = tmp_path / "native"
    with jax.default_device(jax.devices("cpu")[0]):
        params, state = initialize.init(
            jax.random.PRNGKey(1), jnp.ones((1, 2048, 4)), jnp.zeros((1,), jnp.int32)
        )
        with ocp.StandardCheckpointer() as checkpointer:
            checkpointer.save(str(checkpoint), (params, state))
            checkpointer.wait_until_finished()
    with nullcontext() if modern else pytest.warns(UserWarning, match="running statistics"):
        encoder = _encoder(checkpoint)
    before = [parameter.detach().clone() for parameter in encoder.parameters()]
    initial_state = encoder.state
    encoder.train()
    encoder.encode(["ACGT"]).square().sum().backward()
    assert all(parameter.grad.abs().sum() > 0 for parameter in encoder.parameters())
    torch.optim.AdamW(encoder.parameters(), lr=0.01).step()
    assert all(
        not torch.equal(old, new) for old, new in zip(before, encoder.parameters(), strict=True)
    )
    policy = "updated_during_training" if modern else "fixed_checkpoint_statistics"
    assert encoder.normalization_state_policy == policy
    if not modern:
        for old, new in zip(
            jax.tree_util.tree_leaves(initial_state),
            jax.tree_util.tree_leaves(encoder.state),
            strict=True,
        ):
            np.testing.assert_array_equal(old, new)
    saved = tmp_path / "saved"
    encoder.save_pretrained(saved)
    shutil.rmtree(checkpoint)
    with nullcontext() if modern else pytest.warns(UserWarning, match="running statistics"):
        reloaded = _encoder(saved)
    encoder.eval()
    with torch.no_grad():
        torch.testing.assert_close(encoder.encode(["ACGT"]), reloaded.encode(["ACGT"]))
