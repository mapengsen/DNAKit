"""Real, reduced Enformer architecture tests for both training modes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    predict_custom_properties,
    train_custom_property_model,
)
from dnakit.predictions.custom_training.encoder_enformer import EnformerSequenceEncoder
from dnakit.representations import RepresentationConfig, get_embedding_model
from dnakit.representations.backends import create_representation_backend


@pytest.fixture
def checkpoint(tmp_path: Path) -> Path:
    torch = pytest.importorskip("torch")
    library = pytest.importorskip("enformer_pytorch")
    safetensors = pytest.importorskip("safetensors.torch")
    torch.manual_seed(11)
    model = library.Enformer.from_hparams(
        dim=48,
        depth=1,
        heads=2,
        attn_dim_key=8,
        num_downsamples=3,
        dim_divisible_by=8,
        output_heads={"human": 2},
        target_length=2,
        use_tf_gamma=False,
        dropout_rate=0.0,
        attn_dropout=0.0,
        pos_dropout=0.0,
    )
    directory = tmp_path / "original"
    model.config.save_pretrained(directory)
    # Upstream Enformer has duplicate trunk registrations. Clone this tiny
    # fixture's aliases to write a complete safe HF checkpoint for its loader.
    safetensors.save_file(
        {name: tensor.clone() for name, tensor in model.state_dict().items()},
        directory / "model.safetensors",
        metadata={"format": "pt"},
    )
    return directory


def _encoder(checkpoint: Path, pooling: str = "mean") -> EnformerSequenceEncoder:
    return EnformerSequenceEncoder(
        get_embedding_model("enformer"),
        checkpoint,
        device="cpu",
        dtype="float32",
        pooling=pooling,
        max_length=32,
        allow_remote_code=False,
    )


@pytest.mark.parametrize("pooling", ["mean", "max", "cls", "last"])
def test_enformer_gradients_tail_chunks_and_portable_reload(
    checkpoint: Path, tmp_path: Path, pooling: str
) -> None:
    import torch

    encoder = _encoder(checkpoint, pooling)
    assert not encoder.model.heads
    sequences = ["ACGT" * 9, "G", "TTGCA" * 4]
    encoder.eval()
    batched = encoder.encode(sequences)
    single = torch.cat([encoder.encode([sequence]) for sequence in sequences])
    torch.testing.assert_close(batched, single, atol=1e-6, rtol=1e-5)
    before = encoder.model.stem[0].weight.detach().clone()
    optimizer = torch.optim.AdamW(encoder.parameters(), lr=0.001)
    encoder.train()
    encoder.encode(sequences).square().mean().backward()
    assert encoder.model.stem[0].weight.grad is not None
    assert torch.isfinite(encoder.model.stem[0].weight.grad).all()
    optimizer.step()
    assert not torch.equal(before, encoder.model.stem[0].weight)
    destination = tmp_path / "saved"
    encoder.save_pretrained(destination)
    encoder.eval()
    expected = encoder.encode(sequences).detach()
    with torch.no_grad():
        encoder.model.stem[0].weight.add_(10)
    encoder.load_pretrained(destination)
    torch.testing.assert_close(expected, encoder.encode(sequences))
    reloaded = _encoder(destination)
    reloaded.eval()
    torch.testing.assert_close(expected, reloaded.encode(sequences))


def _data(tmp_path: Path) -> Path:
    path = tmp_path / "sequences.csv"
    rows = ["id,sequence,label"]
    for index in range(24):
        sequence = ("ACGT" if index % 2 else "TTGC") * (3 + index % 4)
        rows.append(f"r{index},{sequence},{index % 2}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


@pytest.mark.parametrize("method", ["full", "linear_probe"])
def test_enformer_training_and_prediction_roundtrip(
    checkpoint: Path, tmp_path: Path, method: str
) -> None:
    torch = pytest.importorskip("torch")
    pytest.importorskip("sklearn")
    source = _data(tmp_path)
    output = tmp_path / "trained"
    backend = None
    frozen = {}
    if method == "linear_probe":
        backend = create_representation_backend(
            get_embedding_model("enformer"),
            checkpoint,
            RepresentationConfig(model="enformer", device="cpu", max_length=32),
        )
        frozen = {key: value.clone() for key, value in backend.model.state_dict().items()}
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="classification",
            model="enformer",
            finetune_method=method,
            checkpoint_path=checkpoint,
            output_dir=output,
            epochs=2,
            batch_size=5,
            workers=0,
            monitor="mcc",
            lr_scheduler="reduce_on_plateau",
            scheduler_patience=1,
            scheduler_factor=0.5,
            early_stopping_patience=5,
            gradient_accumulation_steps=1,
            max_length=32,
            device="cpu",
            show_progress=False,
        ),
        representation_backend=backend,
        verbose=False,
    )
    assert result.finetune_method == method
    if backend is not None:
        for key, value in backend.model.state_dict().items():
            assert torch.equal(frozen[key], value)
        assert all(parameter.grad is None for parameter in backend.model.parameters())
    else:
        assert (output / "backbone" / "enformer_encoder.json").is_file()
        state = torch.load(output / "backbone" / "encoder.pt", weights_only=True)
        original = _encoder(checkpoint).model.state_dict()
        assert not torch.equal(state["stem.0.weight"], original["stem.0.weight"])
    prediction = predict_custom_properties(
        source, model_dir=output, device="cpu", batch_size=3, show_progress=False, verbose=False
    )
    repeated = predict_custom_properties(
        source, model_dir=output, device="cpu", batch_size=7, show_progress=False, verbose=False
    )
    assert prediction.record_ids == repeated.record_ids
    assert prediction.predictions == repeated.predictions
    for left, right in zip(prediction.probabilities, repeated.probabilities, strict=True):
        assert left.keys() == right.keys()
        assert left["record_id"] == right["record_id"]
        for name in left["label"]:
            np.testing.assert_allclose(left["label"][name], right["label"][name], atol=1e-6)
    manifest = json.loads((output / "model.json").read_text(encoding="utf-8"))
    assert manifest["dtype"] == "auto"
    settings = manifest["training"]
    assert settings["monitor"] == "mcc"
    assert settings["monitor_mode"] == "max"
    assert settings["lr_scheduler"] == "reduce_on_plateau"
    assert settings["scheduler_patience"] == 1
    assert settings["scheduler_factor"] == 0.5
    assert settings["early_stopping_patience"] == 5
    assert settings["gradient_accumulation_steps"] == 1
    history = json.loads((output / "history.json").read_text(encoding="utf-8"))
    assert all(row["primary_metric_name"] == "mcc" for row in history)
    assert all(row["learning_rate"] == 1e-4 for row in history)


@pytest.mark.parametrize("method", ["full", "linear_probe"])
@pytest.mark.parametrize("validation", [True, False])
def test_enformer_cv_and_final_refit_without_test(
    checkpoint: Path, tmp_path: Path, method: str, validation: bool
) -> None:
    """Fold validation and final refit never require or fabricate a test partition."""
    from dnakit.predictions import DataMapping

    source = tmp_path / "fold.csv"
    source.write_text(
        "id,sequence,label,split\n"
        + "\n".join(
            f"r{i},{'ACGT' * (3 + i % 3)},{i % 2},"
            + ("valid" if validation and i >= 16 else "train")
            for i in range(20)
        )
        + "\n"
    )
    output = tmp_path / "model"
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="classification",
            classification_mode="binary",
            model="enformer",
            finetune_method=method,
            checkpoint_path=checkpoint,
            output_dir=output,
            mapping=DataMapping(
                sequence_field="sequence",
                id_field="id",
                label_fields=("label",),
                split_field="split",
            ),
            epochs=2,
            batch_size=8,
            max_length=32,
            device="cpu",
            workers=0,
            early_stopping_patience=5 if validation else 0,
            lr_scheduler="reduce_on_plateau" if validation else "none",
            monitor="mcc",
            log_train_metrics=True,
            show_progress=False,
        ),
        verbose=False,
    )
    history = json.loads((output / "history.json").read_text())
    assert len(history) == 2
    assert all(isinstance(row["train_mcc"], float) for row in history)
    assert "test" not in result.metrics
    assert not (output / "test_predictions.jsonl").exists()
    if validation:
        predictions = [
            json.loads(line)
            for line in (output / "validation_predictions.jsonl").read_text().splitlines()
        ]
        assert [r["record_id"] for r in predictions] == [f"r{i}" for i in range(16, 20)]
        assert all(isinstance(row["validation_mcc"], float) for row in history)
        assert "validation" in result.metrics
    else:
        assert "validation" not in result.metrics
        assert not (output / "validation_predictions.jsonl").exists()
        assert all(
            row["validation_loss"] is None
            and row["validation_mcc"] is None
            and row["monitor"] == "train_loss"
            and row["monitor_mode"] == "min"
            for row in history
        )
    prediction = predict_custom_properties(
        source, model_dir=output, device="cpu", show_progress=False, verbose=False
    )
    assert len(prediction.record_ids) == 20


def test_zero_holdout_ratios_preserve_every_training_record(tmp_path: Path) -> None:
    from dnakit.predictions.custom_training.data import load_property_data
    from dnakit.predictions.custom_training.splitting import split_property_data

    source = _data(tmp_path)
    data = load_property_data(source)
    for valid_ratio, test_ratio in [(0.0, 0.0), (0.2, 0.0), (0.0, 0.2)]:
        config = CustomPropertyTrainingConfig(
            task_type="classification",
            train_ratio=1 - valid_ratio - test_ratio,
            validation_ratio=valid_ratio,
            test_ratio=test_ratio,
            early_stopping_patience=0,
        )
        split = split_property_data(
            data,
            train_ratio=config.train_ratio,
            validation_ratio=valid_ratio,
            test_ratio=test_ratio,
            seed=42,
            task_type="classification",
            stratify=True,
        )
        assert len(split.train) + len(split.validation) + len(split.test) == 24
        assert bool(split.validation) == bool(valid_ratio)
        assert bool(split.test) == bool(test_ratio)
