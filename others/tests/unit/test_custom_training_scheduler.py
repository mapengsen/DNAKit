"""Validation-driven scheduling, early stopping, and best-checkpoint restoration."""

from __future__ import annotations

import json
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.predictions import CustomPropertyTrainingConfig
from dnakit.predictions.custom_training import training
from dnakit.predictions.custom_training._loading import _PropertyBatchLoader
from dnakit.predictions.custom_training.schema import TargetSchema
from dnakit.predictions.custom_training.splitting import SplitAssignment


@pytest.mark.parametrize(
    ("task", "monitor", "direction"),
    [
        ("classification", "validation_loss", "min"),
        ("classification", "mcc", "max"),
        ("classification", "f1_macro", "max"),
        ("classification", "accuracy", "max"),
        ("regression", "rmse", "min"),
        ("regression", "mae", "min"),
    ],
)
def test_monitor_direction(task: Any, monitor: Any, direction: str) -> None:
    config = CustomPropertyTrainingConfig(task_type=task, monitor=monitor)
    assert config.resolved_monitor_mode == direction
    assert (
        CustomPropertyTrainingConfig(
            task_type=task, monitor=monitor, monitor_mode="max"
        ).resolved_monitor_mode
        == "max"
    )


@pytest.mark.parametrize(
    "options",
    [
        {"lr_scheduler": "cosine"},
        {"scheduler_patience": -1},
        {"scheduler_patience": True},
        {"scheduler_patience": 1.5},
        {"scheduler_factor": 0},
        {"scheduler_factor": 1},
        {"scheduler_factor": float("nan")},
        {"scheduler_factor": True},
        {"monitor": "test_mcc"},
        {"monitor": "rmse"},
        {"monitor_mode": "up"},
        {"monitor_min_delta": -1},
        {"monitor_min_delta": float("inf")},
    ],
)
def test_invalid_scheduler_settings(options: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        CustomPropertyTrainingConfig(task_type="classification", **options)


@pytest.mark.parametrize("full", [False, True])
def test_mcc_controls_lr_early_stopping_and_best_weights(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, full: bool
) -> None:
    torch = pytest.importorskip("torch")

    class TinyEncoder:
        def __init__(self) -> None:
            self.device = torch.device("cpu")
            self.model = torch.nn.Linear(3, 3)

        def parameters(self) -> Any:
            return self.model.parameters()

        def train(self) -> None:
            self.model.train()

        def eval(self) -> None:
            self.model.eval()

        def encode(self, sequences: Any) -> Any:
            return self.model(
                torch.tensor(
                    [[len(s), s.count("G"), s.count("C")] for s in sequences], dtype=torch.float32
                )
                / 4
            )

        def save_pretrained(self, path: Path) -> None:
            path.mkdir(parents=True)
            torch.save(self.model.state_dict(), path / "encoder.pt")

        def load_pretrained(self, path: Path) -> None:
            self.model.load_state_dict(torch.load(path / "encoder.pt", weights_only=True))

    optimizers = []

    class RecordingAdamW(torch.optim.AdamW):  # type: ignore[name-defined, misc]
        def __init__(self, parameters: Any, **kwargs: Any) -> None:
            super().__init__(parameters, **kwargs)
            self.step_count = 0
            optimizers.append(self)

        def step(self, closure: Any = None) -> Any:
            self.step_count += 1
            return super().step(closure)

    monkeypatch.setattr(torch.optim, "AdamW", RecordingAdamW)
    encoder = TinyEncoder() if full else None
    snapshots: list[tuple[Any, Any]] = []
    before = encoder.model.weight.detach().clone() if encoder is not None else None

    def evaluate(*args: Any, **kwargs: Any) -> tuple[float, Any, Any]:
        snapshots.append(
            (
                {k: v.detach().clone() for k, v in kwargs["head"].state_dict().items()},
                {k: v.detach().clone() for k, v in encoder.model.state_dict().items()}
                if encoder is not None
                else None,
            )
        )
        # Loss keeps improving, whereas MCC peaks in epoch 2 and then plateaus.
        return 1.0 / len(snapshots), np.zeros(4), np.zeros((4, 2))

    scores = iter([-0.5, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25])
    monkeypatch.setattr(training, "_evaluate", evaluate)
    monkeypatch.setattr(
        training,
        "compute_metrics",
        lambda *a, **kw: {
            "matthews_correlation_coefficient": next(scores),
            "f1_macro": 0.9,
        },
    )
    schema = TargetSchema("classification", "binary", ("label",), (("0", "1"),))
    config = CustomPropertyTrainingConfig(
        task_type="classification",
        model="enformer",
        finetune_method="full" if full else "linear_probe",
        epochs=20,
        batch_size=64,
        workers=0,
        device="cpu",
        show_progress=False,
        learning_rate=1e-4,
        beta1=0.9,
        beta2=0.95,
        weight_decay=0.1,
        gradient_accumulation_steps=1,
        early_stopping_patience=5,
        monitor="mcc",
        lr_scheduler="reduce_on_plateau",
        scheduler_patience=1,
        scheduler_factor=0.5,
    )
    split = SplitAssignment(tuple(range(130)), (130, 131, 132, 133), (134,), "test")
    with closing(
        _PropertyBatchLoader(
            torch,
            ["ACGT"] * 135,
            np.arange(135) % 2,
            representations=None if full else np.ones((135, 3), dtype=np.float32),
            batch_size=64,
            workers=0,
            seed=42,
        )
    ) as loader:
        head, best_encoder, history = training._train_head_and_encoder(
            split,
            schema,
            config,
            tmp_path,
            batch_loader=loader,
            encoder=encoder,
            input_dimension=3,
        )
    assert len(history) == 7  # Five bad epochs after the epoch-2 maximum MCC.
    assert [row["learning_rate"] for row in history] == pytest.approx(
        [1e-4, 1e-4, 1e-4, 1e-4, 5e-5, 5e-5, 2.5e-5]
    )
    assert [row["next_learning_rate"] for row in history] == pytest.approx(
        [1e-4, 1e-4, 1e-4, 5e-5, 5e-5, 2.5e-5, 2.5e-5]
    )
    assert all(row["monitor_mode"] == "max" for row in history)
    assert history[1]["primary_metric"] == -0.25
    assert history[1]["primary_metric_name"] == "mcc"
    assert optimizers[0].step_count == 7 * 3  # 64, 64, 2; one update per batch.
    assert optimizers[0].defaults["betas"] == (0.9, 0.95)
    assert optimizers[0].defaults["weight_decay"] == 0.1
    for name, value in head.state_dict().items():
        torch.testing.assert_close(value, snapshots[1][0][name])
    assert not torch.equal(head.weight, snapshots[-1][0]["weight"])
    if full:
        assert best_encoder is encoder
        assert encoder is not None
        assert not torch.equal(before, encoder.model.weight)
        for name, value in encoder.model.state_dict().items():
            torch.testing.assert_close(value, snapshots[1][1][name])


@pytest.mark.parametrize("invalid_score", [None, float("nan"), float("inf")])
def test_unavailable_monitor_is_not_silently_replaced_by_loss(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, invalid_score: Any
) -> None:
    torch = pytest.importorskip("torch")
    monkeypatch.setattr(training, "_evaluate", lambda *a, **kw: (0.1, None, None))
    monkeypatch.setattr(
        training,
        "compute_metrics",
        lambda *a, **kw: {
            "matthews_correlation_coefficient": invalid_score,
        },
    )
    schema = TargetSchema("classification", "binary", ("label",), (("0", "1"),))
    config = CustomPropertyTrainingConfig(
        task_type="classification",
        finetune_method="linear_probe",
        monitor="mcc",
        workers=0,
        device="cpu",
        show_progress=False,
    )
    with (
        closing(
            _PropertyBatchLoader(
                torch,
                ["ACGT"] * 6,
                np.arange(6) % 2,
                representations=np.ones((6, 3), dtype=np.float32),
                batch_size=32,
                workers=0,
                seed=42,
            )
        ) as loader,
        pytest.raises(BackendExecutionError, match="monitored validation metric"),
    ):
        training._train_head_and_encoder(
            SplitAssignment((0, 1), (2, 3), (4, 5), "test"),
            schema,
            config,
            tmp_path,
            batch_loader=loader,
            encoder=None,
            input_dimension=3,
        )


@pytest.mark.parametrize("scheduler", ["none", "reduce_on_plateau"])
def test_regression_min_mode_and_optional_scheduler(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scheduler: Any
) -> None:
    torch = pytest.importorskip("torch")
    validation_scores = iter([0.5, 0.2, 0.3, 0.4])
    losses = iter([4.0, 3.0, 2.0, 1.0])
    monkeypatch.setattr(training, "_evaluate", lambda *a, **kw: (next(losses), None, None))
    monkeypatch.setattr(
        training,
        "compute_metrics",
        lambda *a, **kw: {"aggregate": {"rmse": next(validation_scores)}},
    )
    config = CustomPropertyTrainingConfig(
        task_type="regression",
        finetune_method="linear_probe",
        monitor="rmse",
        lr_scheduler=scheduler,
        scheduler_patience=0,
        scheduler_factor=0.5,
        early_stopping_patience=2,
        epochs=10,
        workers=0,
        device="cpu",
        show_progress=False,
    )
    schema = TargetSchema(
        "regression",
        "regression",
        ("label",),
        regression_means=(0.0,),
        regression_scales=(1.0,),
    )
    with closing(
        _PropertyBatchLoader(
            torch,
            ["ACGT"] * 6,
            np.ones((6, 1), dtype=np.float32),
            representations=np.ones((6, 3), dtype=np.float32),
            batch_size=32,
            workers=0,
            seed=42,
        )
    ) as loader:
        _, _, history = training._train_head_and_encoder(
            SplitAssignment((0, 1), (2, 3), (4, 5), "test"),
            schema,
            config,
            tmp_path,
            batch_loader=loader,
            encoder=None,
            input_dimension=3,
        )
    assert len(history) == 4
    assert all(row["monitor_mode"] == "min" for row in history)
    assert history[1]["primary_metric"] == 0.2
    expected = [1e-4] * 4 if scheduler == "none" else [1e-4, 1e-4, 5e-5, 2.5e-5]
    assert [row["next_learning_rate"] for row in history] == pytest.approx(expected)


def test_cli_passes_requested_finetuning_settings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from typer.testing import CliRunner

    from dnakit.cli.app import app

    captured = []

    def train(source: Path, *, config: CustomPropertyTrainingConfig) -> Any:
        captured.append(config)
        return SimpleNamespace(to_dict=lambda: {"monitor": config.monitor})

    monkeypatch.setattr("dnakit.predictions.train_custom_property_model", train)
    source = tmp_path / "labels.csv"
    source.write_text("sequence,label\nACGT,1\n")
    result = CliRunner().invoke(
        app,
        [
            "train-property",
            str(source),
            str(tmp_path / "model"),
            "--task-type",
            "classification",
            "--model",
            "enformer",
            "--finetune-method",
            "full",
            "--batch-size",
            "64",
            "--learning-rate",
            "1e-4",
            "--beta1",
            "0.9",
            "--beta2",
            "0.95",
            "--weight-decay",
            "0.1",
            "--gradient-accumulation-steps",
            "1",
            "--early-stopping-patience",
            "5",
            "--lr-scheduler",
            "reduce_on_plateau",
            "--scheduler-patience",
            "1",
            "--scheduler-factor",
            "0.5",
            "--monitor",
            "mcc",
            "--monitor-mode",
            "auto",
        ],
    )
    assert result.exit_code == 0, result.output
    config = captured[0]
    assert (config.batch_size, config.gradient_accumulation_steps) == (64, 1)
    assert config.resolved_monitor_mode == "max"
    assert config.lr_scheduler == "reduce_on_plateau"
    assert (config.scheduler_patience, config.scheduler_factor) == (1, 0.5)
    assert config.early_stopping_patience == 5
    assert json.loads(result.output)["monitor"] == "mcc"
