"""Gradient, padding and artifact checks for Transformers training adapters."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from dnakit.predictions.custom_training.encoder import TransformersSequenceEncoder
from dnakit.representations import get_embedding_model
from dnakit.representations.backends import _forward_transformers_batch


def _tiny_checkpoint(directory: Path, *, causal: bool = False) -> None:
    transformers = pytest.importorskip("transformers")
    directory.mkdir()
    vocabulary = directory / "vocab.txt"
    vocabulary.write_text("[PAD]\n[UNK]\n[CLS]\n[SEP]\n[MASK]\nA\nC\nG\nT\nN\n")
    transformers.BertTokenizer(vocab_file=str(vocabulary), do_lower_case=False).save_pretrained(
        directory
    )
    if causal:
        model = transformers.LlamaModel(
            transformers.LlamaConfig(
                vocab_size=10,
                hidden_size=8,
                intermediate_size=16,
                num_hidden_layers=1,
                num_attention_heads=2,
                num_key_value_heads=2,
                max_position_embeddings=32,
            )
        )
    else:
        model = transformers.BertModel(
            transformers.BertConfig(
                vocab_size=10,
                hidden_size=8,
                intermediate_size=16,
                num_hidden_layers=1,
                num_attention_heads=2,
                max_position_embeddings=32,
            )
        )
    model.save_pretrained(directory)


@pytest.mark.parametrize("causal", [False, True])
def test_transformers_gradient_export_and_in_place_restore(tmp_path: Path, causal: bool) -> None:
    torch = pytest.importorskip("torch")
    checkpoint = tmp_path / "checkpoint"
    _tiny_checkpoint(checkpoint, causal=causal)
    encoder = TransformersSequenceEncoder(
        get_embedding_model("generator" if causal else "grover"),
        checkpoint,
        device="cpu",
        dtype="float32",
        pooling="mean",
        max_length=64,
        allow_remote_code=True,
    )
    assert encoder.max_length == 32
    encoder.train()
    optimizer = torch.optim.AdamW(encoder.parameters(), lr=0.01)
    before = encoder.model.get_input_embeddings().weight.detach().clone()
    encoder.encode(["A C G T", "T G"]).square().sum().backward()
    optimizer.step()
    assert not torch.equal(before, encoder.model.get_input_embeddings().weight)
    encoder.eval()
    expected = encoder.encode(["A C G T", "T G"]).detach()
    saved = tmp_path / "saved"
    encoder.save_pretrained(saved)
    with torch.no_grad():
        for parameter in encoder.parameters():
            parameter.add_(1)
    model_identity = id(encoder.model)
    encoder.load_pretrained(saved)
    assert id(encoder.model) == model_identity
    torch.testing.assert_close(encoder.encode(["A C G T", "T G"]), expected)
    reloaded = TransformersSequenceEncoder(
        encoder.spec,
        saved,
        device="cpu",
        dtype="float32",
        pooling="mean",
        max_length=64,
        allow_remote_code=True,
    )
    reloaded.eval()
    torch.testing.assert_close(reloaded.encode(["A C G T", "T G"]), expected)


def test_maskless_bidirectional_forward_excludes_padding_and_retains_gradients() -> None:
    torch = pytest.importorskip("torch")

    class ReverseEncoder(torch.nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.scale = torch.nn.Parameter(torch.tensor(0.5))

        def forward(self, input_ids: Any) -> Any:
            return input_ids.flip(1).cumsum(1).flip(1).unsqueeze(-1).float() * self.scale

    model = ReverseEncoder()
    inputs = {
        "input_ids": torch.tensor([[1, 2, 99], [3, 4, 5]]),
        "attention_mask": torch.tensor([[1, 1, 0], [1, 1, 1]]),
        "token_type_ids": torch.zeros(2, 3, dtype=torch.long),
    }
    pooled = _forward_transformers_batch(
        model, inputs, torch=torch, pooling="mean", special_tokens_mask=None
    )
    torch.testing.assert_close(pooled, torch.tensor([[1.25], [13 / 3]]))
    pooled.sum().backward()
    assert model.scale.grad is not None and model.scale.grad.item() > 0


_LOCAL_CHECKPOINTS = {
    "dnabert2": "DNABERT_2/weights/DNABERT-2-117M",
    "ntv2": "nucleotide-transformer/weights/nucleotide-transformer-v2-500m-multi-species",
    "hyenadna": "hyena-dna/weights/hyenadna-medium-450k-seqlen-hf",
    "caduceus": "caduceus/weights/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",
    "grover": "GROVER/model",
    "lucaone": "LucaOne/weights/LucaOne-gene-step36.8M",
    "generator": "GENERator/weights/GENERator-v2-eukaryote-1.2b-base",
}


@pytest.mark.parametrize("name", tuple(_LOCAL_CHECKPOINTS))
def test_local_official_architecture_gradient_reload_and_frozen_features(
    tmp_path: Path, name: str
) -> None:
    """Use installed official source/tokenizers, never download pretrained weights."""

    torch = pytest.importorskip("torch")
    transformers = pytest.importorskip("transformers")
    source = (
        Path(__file__).resolve().parents[3] / "compared_models/models" / _LOCAL_CHECKPOINTS[name]
    )
    if not source.is_dir():
        pytest.skip("Optional local official checkpoint source is unavailable")
    if name == "caduceus":
        pytest.importorskip("mamba_ssm")
        if not torch.cuda.is_available():
            pytest.skip("Caduceus selective-scan kernels require CUDA")
    checkpoint = tmp_path / name
    checkpoint.mkdir()
    for path in source.iterdir():
        if path.suffix in {".json", ".py", ".txt", ".model"} and ".index." not in path.name:
            shutil.copy2(path, checkpoint / path.name)
    config = transformers.AutoConfig.from_pretrained(
        checkpoint, trust_remote_code=True, local_files_only=True
    )
    sizes = {
        "hidden_size": 16,
        "intermediate_size": 32,
        "num_hidden_layers": 1,
        "num_attention_heads": 2,
        "num_key_value_heads": 2,
        "ffn_dim": 32,
        "d_model": 16,
        "d_inner": 32,
        "n_layer": 1,
        "filter_order": 8,
        "max_position_embeddings": 64,
        "max_seq_len": 64,
        "emb_dim": 3,
    }
    for field, value in sizes.items():
        if hasattr(config, field):
            setattr(config, field, value)
    if name == "caduceus":
        config.fused_add_norm = False
        config.rms_norm = False
        config.ssm_cfg["use_fast_path"] = False
    config.save_pretrained(checkpoint)
    model = transformers.AutoModel.from_config(config, trust_remote_code=True)
    # Caduceus shares bidirectional weights, which must be serialized explicitly.
    model.save_pretrained(
        checkpoint,
        state_dict={key: value.clone() for key, value in model.state_dict().items()},
    )
    del model
    encoder = TransformersSequenceEncoder(
        get_embedding_model(name),
        checkpoint,
        device="cuda:0" if name == "caduceus" else "cpu",
        dtype="float32",
        pooling="mean",
        max_length=32,
        allow_remote_code=True,
    )
    sequences = ["ACGTACGTACGTAC", "TTACGTA"]
    from dnakit.representations.backends import _TransformersBackend
    from dnakit.representations.models import RepresentationConfig

    encoder.eval()
    untrained_backend = _TransformersBackend(
        encoder.spec,
        checkpoint,
        RepresentationConfig(model=name, device=str(encoder.device), max_length=32),
    )
    with torch.no_grad():
        untrained_full = encoder.encode(sequences)
    torch.testing.assert_close(
        torch.from_numpy(untrained_backend.extract(sequences, show_progress=False)).to(
            encoder.device
        ),
        untrained_full,
        atol=1e-5,
        rtol=1e-4,
    )
    del untrained_backend
    encoder.train()
    optimizer = torch.optim.AdamW(encoder.parameters(), lr=0.001)
    output = encoder.encode(sequences)
    assert output.shape == (2, 16)
    output[:, 0].sum().backward()
    assert any(p.grad is not None and bool(p.grad.abs().sum() > 0) for p in encoder.parameters())
    optimizer.step()
    encoder.eval()
    with torch.no_grad():
        expected = encoder.encode(sequences)
        individually = torch.cat([encoder.encode([sequence]) for sequence in sequences])
    torch.testing.assert_close(expected, individually, atol=1e-5, rtol=1e-4)
    saved = tmp_path / "saved"
    encoder.save_pretrained(saved)
    reloaded = TransformersSequenceEncoder(
        encoder.spec,
        saved,
        device=str(encoder.device),
        dtype="float32",
        pooling="mean",
        max_length=32,
        allow_remote_code=True,
    )
    reloaded.eval()
    with torch.no_grad():
        torch.testing.assert_close(reloaded.encode(sequences), expected, atol=1e-5, rtol=1e-4)
    if name == "caduceus":
        assert (saved / "tokenization_caduceus.py").is_file()
    if name == "lucaone":
        assert encoder.tokenizer.encode("ACGT", seq_type="gene", add_special_tokens=False) == [
            5,
            7,
            8,
            6,
        ]
    if name == "generator":
        assert json.loads((saved / "config.json").read_text())["use_cache"] is False
    backend = _TransformersBackend(
        encoder.spec,
        saved,
        RepresentationConfig(
            model=name,
            device=str(encoder.device),
            max_length=32,
            allow_remote_code=True,
            show_progress=False,
        ),
    )
    frozen_before = [parameter.detach().clone() for parameter in backend.model.parameters()]
    vectors = backend.extract(sequences, show_progress=False)
    assert vectors.shape == (2, 16)
    torch.testing.assert_close(
        torch.from_numpy(vectors).to(expected.device), expected, atol=1e-5, rtol=1e-4
    )
    for before, parameter in zip(frozen_before, backend.model.parameters(), strict=True):
        torch.testing.assert_close(before, parameter)
        assert parameter.grad is None


def test_fast_tokenizer_windows_preserve_special_tokens_and_tail(tmp_path: Path) -> None:
    torch = pytest.importorskip("torch")
    transformers = pytest.importorskip("transformers")
    tokenizers = pytest.importorskip("tokenizers")
    from dnakit.representations.backends import _TransformersBackend
    from dnakit.representations.models import RepresentationConfig

    torch.manual_seed(19)
    checkpoint = tmp_path / "fast-tokenizer"
    checkpoint.mkdir()
    vocabulary = {
        token: index
        for index, token in enumerate(["[PAD]", "[UNK]", "[CLS]", "[SEP]", "A", "C", "G", "T"])
    }
    native = tokenizers.Tokenizer(tokenizers.models.WordLevel(vocabulary, unk_token="[UNK]"))
    native.pre_tokenizer = tokenizers.pre_tokenizers.Split("", behavior="isolated")
    native.post_processor = tokenizers.processors.TemplateProcessing(
        single="[CLS] $A [SEP]", special_tokens=[("[CLS]", 2), ("[SEP]", 3)]
    )
    transformers.PreTrainedTokenizerFast(
        tokenizer_object=native,
        unk_token="[UNK]",
        pad_token="[PAD]",
        cls_token="[CLS]",
        sep_token="[SEP]",
        model_max_length=32,
    ).save_pretrained(checkpoint)
    transformers.BertModel(
        transformers.BertConfig(
            vocab_size=8,
            hidden_size=8,
            intermediate_size=16,
            num_hidden_layers=1,
            num_attention_heads=2,
            max_position_embeddings=32,
        )
    ).save_pretrained(checkpoint)
    spec = get_embedding_model("grover")
    encoder = TransformersSequenceEncoder(
        spec,
        checkpoint,
        device="cpu",
        dtype="float32",
        pooling="mean",
        max_length=6,
        allow_remote_code=False,
    )
    encoder.eval()
    backend = _TransformersBackend(
        spec,
        checkpoint,
        RepresentationConfig(model="grover", device="cpu", max_length=6, batch_size=8),
    )
    chunks, kind = backend._chunks("ACGTAC")
    batch = backend._tokenize(chunks, kind)
    assert batch["input_ids"].tolist() == [[2, 4, 5, 6, 7, 3], [2, 4, 5, 3, 0, 0]]
    assert batch["special_tokens_mask"].tolist() == [[1, 0, 0, 0, 0, 1], [1, 0, 0, 1, 1, 1]]
    with torch.no_grad():
        expected = encoder.encode(["ACGT", "AC"]).mean(dim=0)
        single = encoder.encode(["ACGT"])[0]
    actual = backend.extract(["ACGTAC", "ACGT"], show_progress=False)
    torch.testing.assert_close(torch.from_numpy(actual[0]), expected)
    torch.testing.assert_close(torch.from_numpy(actual[1]), single)
    # A one-window batch must cover and aggregate exactly the same tail bases.
    backend.config = RepresentationConfig(model="grover", device="cpu", max_length=6, batch_size=1)
    torch.testing.assert_close(
        torch.from_numpy(backend.extract(["ACGTAC"], show_progress=False)[0]), expected
    )
