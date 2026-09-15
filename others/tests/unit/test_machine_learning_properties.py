"""Exercise every ML family, saved inference, target modes and feature isolation."""

from __future__ import annotations

import builtins
import csv
import importlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError
from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    DataMapping,
    available_machine_learning_models,
    predict_custom_properties,
    train_custom_property_model,
)
from dnakit.predictions.custom_training import ml_training
from dnakit.predictions.custom_training.ml_models import create_estimator
from dnakit.tools import serialize_for_agent


def _source(tmp_path: Path, task: str, mode: str = "single") -> tuple[Path, DataMapping]:
    rng = np.random.default_rng(42)
    x = rng.uniform(0.1, 3, (90, 4))
    rows: list[dict[str, object]] = []
    label_names: tuple[str, ...] = (
        ("target",) if mode in {"single", "binary"} else ("target", "second")
    )
    if mode == "onehot":
        label_names = ("a", "b", "c")
    for index, values in enumerate(x):
        row: dict[str, object] = {
            "id": f"row-{index}",
            "split": "train" if index < 60 else "validation" if index < 75 else "test",
            "sequence": "A" * (index % 9 + 1) + "CGT" * (index % 7 + 1),
            **{f"x{i}": float(value) for i, value in enumerate(values)},
        }
        if task == "regression":
            row.update(target=float(2 * values[0] - values[1]), second=float(values[2] * 4))
        elif mode == "multilabel":
            row.update(target=int(values[0] > 1.5), second=1)
        elif mode == "onehot":
            row.update({name: int(index % 3 == i) for i, name in enumerate(label_names)})
        elif mode == "binary":
            row.update(target=str(index % 2), second=str(index % 2))
        else:
            row.update(target=str(index % 3), second=str(index % 2))
        rows.append(row)
    path = tmp_path / "data.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path, DataMapping(
        sequence_field=None,
        feature_fields=tuple(f"x{i}" for i in range(4)),
        label_fields=label_names,
        id_field="id",
        split_field="split",
    )


def _params(name: str) -> dict[str, object]:
    if name in {
        "random_forest",
        "extra_trees",
        "gradient_boosting",
        "adaboost",
        "bagging",
        "xgboost",
        "lightgbm",
    }:
        return {"n_estimators": 4}
    if name == "catboost":
        return {"iterations": 4, "depth": 2}
    if name == "hist_gradient_boosting":
        return {"max_iter": 4}
    if name == "mlp":
        return {"hidden_layer_sizes": [8], "max_iter": 20}
    return {}


_CASES = [(spec.name, task) for spec in available_machine_learning_models() for task in spec.tasks]


@pytest.mark.parametrize("task", ["classification", "regression"])
def test_complete_cv_is_isolated_and_refits_all_training_rows(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    task: Any,
) -> None:
    source, mapping = _source(tmp_path, task, mode="binary")
    mapping = replace(mapping, split_field=None)
    observed: list[tuple[Any, Any]] = []
    original_fit = ml_training._fit_estimators

    def fit(values: Any, targets: Any, config: Any) -> Any:
        observed.append((values.copy(), targets.copy()))
        return original_fit(values, targets, config)

    monkeypatch.setattr(ml_training, "_fit_estimators", fit)
    trained = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type=task,
            model="svm" if task == "classification" else "ridge",
            finetune_method="machine_learning",
            feature_source="numeric",
            mapping=mapping,
            model_params={"kernel": "linear", "probability": True}
            if task == "classification"
            else {},
            cv_folds=10,
            train_ratio=0.9,
            validation_ratio=0.1,
            test_ratio=0.0,
            workers=1,
            output_dir=tmp_path / "cv-model",
            show_progress=False,
        ),
    )
    assert len(observed) == 11
    assert [len(values) for values, _ in observed] == [81] * 10 + [90]
    assert trained.split_counts == {"train": 90, "validation": 0, "test": 0}
    with (tmp_path / "cv-model" / "splits.csv").open() as stream:
        assignments = list(csv.DictReader(stream))
    ids = [int(row["validation_fold"]) for row in assignments]
    assert sorted(int(row["row_index"]) for row in assignments) == list(range(90))
    assert [ids.count(fold) for fold in range(10)] == [9] * 10
    all_values = observed[-1][0]
    for fold in range(10):
        expected = all_values[[index for index, assigned in enumerate(ids) if assigned != fold]]
        np.testing.assert_allclose(observed[fold][0], expected)
        if task == "regression":
            np.testing.assert_allclose(observed[fold][1].mean(axis=0), 0, atol=1e-6)
    metrics = trained.metrics["cross_validation"]
    assert metrics["n_splits"] == 10
    assert len(metrics["folds"]) == 10
    assert metrics["refit_record_count"] == 90
    assert "test" not in trained.metrics
    oof = (tmp_path / "cv-model" / "oof_predictions.jsonl").read_text().splitlines()
    assert len(oof) == 90
    assert len({json.loads(line)["record_id"] for line in oof}) == 90
    assert not (tmp_path / "cv-model" / "test_predictions.jsonl").exists()

    # Prediction writes through the existing table API, keeping its Agent result small.
    prediction_path = tmp_path / "predictions.csv"
    predicted = predict_custom_properties(
        source,
        model_dir=trained.output_dir,
        trust_model=True,
        output_path=prediction_path,
    )
    summary = serialize_for_agent(predicted)
    assert summary["target_path"] == str(prediction_path)
    assert summary["row_count"] == 90
    assert "predictions" not in summary
    assert len(json.dumps(summary)) < 2_000
    with prediction_path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(predicted.predictions) == 90
    assert rows[0]["record_id"] == predicted.predictions[0]["record_id"]
    if task == "classification":
        assert {row["target"] for row in rows} <= {"0", "1"}
    saved = prediction_path.read_bytes()
    with pytest.raises(FileExistsError):
        predict_custom_properties(
            source, model_dir=trained.output_dir, trust_model=True, output_path=prediction_path
        )
    assert prediction_path.read_bytes() == saved
    with pytest.raises(ConfigurationError, match="UNSAFE_CUSTOM_OUTPUT_PATH"):
        predict_custom_properties(
            source,
            model_dir=trained.output_dir,
            trust_model=True,
            output_path=source,
            overwrite=True,
        )


@pytest.mark.parametrize("folds", [1, True, 1.5, 101])
def test_cv_config_rejects_invalid_fold_counts(folds: Any) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_CUSTOM_CV"):
        CustomPropertyTrainingConfig(
            task_type="classification",
            model="svm",
            finetune_method="machine_learning",
            cv_folds=folds,
        )


def test_cv_rejects_mapped_splits_without_consuming_test_rows(tmp_path: Path) -> None:
    source, mapping = _source(tmp_path, "classification", mode="binary")
    with pytest.raises(InputFormatError, match="INVALID_CUSTOM_CV"):
        train_custom_property_model(
            source,
            config=CustomPropertyTrainingConfig(
                task_type="classification",
                model="svm",
                finetune_method="machine_learning",
                feature_source="numeric",
                mapping=mapping,
                cv_folds=10,
                output_dir=tmp_path / "invalid-model",
            ),
        )
    assert not (tmp_path / "invalid-model").exists()


@pytest.mark.parametrize(("model", "task"), _CASES)
def test_every_model_trains_and_reloads(tmp_path: Path, model: str, task: Any) -> None:
    if model in {"xgboost", "lightgbm", "catboost"}:
        pytest.importorskip(model)
    source, mapping = _source(tmp_path, task)
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type=task,
            model=model,
            finetune_method="machine_learning",
            feature_source="numeric",
            mapping=mapping,
            model_params=_params(model),
            output_dir=tmp_path / "model",
            workers=1,
        ),
    )
    assert result.split_counts == {"train": 60, "validation": 15, "test": 15}
    assert result.history == ()
    predicted = predict_custom_properties(source, model_dir=result.output_dir, trust_model=True)
    stored = [
        json.loads(line)
        for line in (Path(result.output_dir) / "test_predictions.jsonl").read_text().splitlines()
    ]
    for saved, actual in zip(stored, predicted.predictions[75:], strict=True):
        if task == "classification":
            assert saved["target"] == actual["target"]
        else:
            assert saved["target"] == pytest.approx(actual["target"])
    assert len(predicted.predictions) == 90
    test_metrics = cast(Mapping[str, Any], result.metrics["test"])
    if task == "classification":
        for row in predicted.probabilities:
            probabilities = row["target"]
            assert isinstance(probabilities, Mapping)
            assert sum(cast(Mapping[str, float], probabilities).values()) == pytest.approx(1.0)
        assert test_metrics["accuracy"] is not None
    else:
        assert test_metrics["aggregate"]["rmse"] is not None


@pytest.mark.parametrize(
    "model", ["random_forest", "linear_svm", "xgboost", "lightgbm", "catboost"]
)
@pytest.mark.parametrize(
    ("task", "mode"),
    [
        ("classification", "binary"),
        ("classification", "multilabel"),
        ("classification", "multioutput"),
        ("classification", "onehot"),
        ("regression", "multioutput"),
    ],
)
def test_multiple_targets(tmp_path: Path, model: str, task: Any, mode: str) -> None:
    if model in {"xgboost", "lightgbm", "catboost"}:
        pytest.importorskip(model)
    source, mapping = _source(tmp_path, task, mode)
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type=task,
            model=model,
            finetune_method="machine_learning",
            feature_source="numeric",
            mapping=mapping,
            classification_mode="multilabel" if mode == "multilabel" else "auto",
            model_params=_params(model),
            output_dir=tmp_path / "model",
            workers=1,
        ),
    )
    predicted = predict_custom_properties(source, model_dir=result.output_dir, trust_model=True)
    assert len(predicted.predictions) == 90
    if mode == "multilabel":
        assert all(row["second"] == 1 for row in predicted.predictions)
        assert any("constant" in warning for warning in result.warnings)
    if mode == "multioutput" and task == "classification":
        test_metrics = cast(Mapping[str, Any], result.metrics["test"])
        assert test_metrics["aggregate"]["exact_match_accuracy"] is not None


class _Embedding:
    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        return np.asarray([[len(seq), seq.count("A")] for seq in sequences], dtype=np.float32)


@pytest.mark.parametrize("feature_source", ["sequence", "embedding", "numeric"])
def test_feature_sources_without_torch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, feature_source: Any
) -> None:
    source, numeric_mapping = _source(tmp_path, "regression")
    mapping = (
        numeric_mapping
        if feature_source == "numeric"
        else DataMapping(
            sequence_field="sequence",
            label_fields=("target",),
            split_field="split",
            id_field="id",
        )
    )
    original_import = builtins.__import__

    def guarded(name: str, *args: Any, **kwargs: Any) -> Any:
        if name.split(".")[0] in {"torch", "transformers", "huggingface_hub"}:
            raise AssertionError(f"Unnecessary deep-learning import: {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    backend = _Embedding() if feature_source == "embedding" else None
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="regression",
            model="ridge",
            finetune_method="machine_learning",
            feature_source=feature_source,
            embedding_model="dnabert2",
            mapping=mapping,
            output_dir=tmp_path / "model",
            workers=1,
        ),
        representation_backend=backend,
    )
    prediction = predict_custom_properties(
        source, model_dir=result.output_dir, trust_model=True, representation_backend=backend
    )
    assert len(prediction.predictions) == 90


def test_scaler_uses_only_training_rows_and_unlabelled_prediction(tmp_path: Path) -> None:
    joblib = pytest.importorskip("joblib")

    source, mapping = _source(tmp_path, "regression")
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="regression",
            model="linear_regression",
            finetune_method="machine_learning",
            feature_source="numeric",
            mapping=mapping,
            output_dir=tmp_path / "model",
            workers=1,
        ),
    )
    estimator = joblib.load(Path(result.output_dir) / "estimator.joblib")["estimators"][0]
    with source.open() as handle:
        rows = list(csv.DictReader(handle))
    training = np.asarray([[float(row[f"x{i}"]) for i in range(4)] for row in rows[:60]])
    np.testing.assert_allclose(estimator.named_steps["scale"].mean_, training.mean(axis=0))
    query = tmp_path / "query.csv"
    query.write_text("id,x0,x1,x2,x3\nnew,1000,2,3,4\n")
    predicted = predict_custom_properties(
        query,
        model_dir=result.output_dir,
        trust_model=True,
    )
    assert predicted.predictions[0]["target"] == pytest.approx(1998, rel=1e-5)
    with pytest.raises(ConfigurationError, match="trust_model"):
        predict_custom_properties(query, model_dir=result.output_dir)
    with pytest.raises(InputFormatError) as mismatch:
        predict_custom_properties(
            query,
            model_dir=result.output_dir,
            trust_model=True,
            mapping=DataMapping(
                sequence_field=None,
                feature_fields=tuple(reversed(mapping.feature_fields)),
                id_field="id",
            ),
        )
    assert mismatch.value.code == "CUSTOM_FEATURE_SCHEMA_MISMATCH"
    with pytest.raises(FileExistsError):
        train_custom_property_model(
            source,
            config=CustomPropertyTrainingConfig(
                task_type="regression",
                model="ridge",
                finetune_method="machine_learning",
                feature_source="numeric",
                mapping=mapping,
                output_dir=result.output_dir,
            ),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"model": "unknown"},
        {"model": "lasso"},
        {"feature_source": "unknown"},
        {"feature_source": "numeric"},
        {"kmer_k": 7},
        {"kmer_k": True},
        {"model_params": {"bad": float("nan")}},
        {"model_params": {"bad": object()}},
    ],
)
def test_invalid_configuration(kwargs: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        CustomPropertyTrainingConfig(
            **{
                "task_type": "classification",
                "model": "random_forest",
                "finetune_method": "machine_learning",
                **kwargs,
            }
        )


def test_missing_optional_dependency(monkeypatch: pytest.MonkeyPatch) -> None:
    original = importlib.import_module

    def missing(name: str) -> Any:
        if name == "catboost":
            raise ImportError("not installed")
        return original(name)

    monkeypatch.setattr(importlib, "import_module", missing)
    config = CustomPropertyTrainingConfig(
        task_type="regression", model="catboost", finetune_method="machine_learning"
    )
    with pytest.raises(BackendUnavailableError) as error:
        create_estimator(config)
    assert error.value.code == "MISSING_MACHINE_LEARNING_DEPENDENCY"


@pytest.mark.parametrize("value", ["", "NaN", "inf", "text"])
def test_invalid_numeric_features(tmp_path: Path, value: str) -> None:
    source = tmp_path / "data.csv"
    source.write_text(f"x,target\n{value},1\n")
    with pytest.raises(InputFormatError) as error:
        train_custom_property_model(
            source,
            config=CustomPropertyTrainingConfig(
                task_type="regression",
                model="ridge",
                finetune_method="machine_learning",
                feature_source="numeric",
                mapping=DataMapping(
                    sequence_field=None, feature_fields=("x",), label_fields=("target",)
                ),
                output_dir=tmp_path / "model",
            ),
        )
    assert error.value.code == "INVALID_CUSTOM_FEATURE"


def test_feature_label_overlap_is_rejected() -> None:
    with pytest.raises(ConfigurationError):
        DataMapping(feature_fields=("target",), label_fields=("target",))


def test_weight_parameters_survive_config_replacement(tmp_path: Path) -> None:
    source, mapping = _source(tmp_path, "classification", "binary")
    config = CustomPropertyTrainingConfig(
        task_type="classification",
        model="logistic_regression",
        finetune_method="machine_learning",
        feature_source="numeric",
        mapping=mapping,
        model_params={"class_weight": {0: 1.0, 1: 2.0}},
        output_dir=tmp_path / "model",
    )
    copied = replace(config, workers=1)
    assert create_estimator(copied).class_weight == {0: 1.0, 1: 2.0}
    result = train_custom_property_model(source, config=copied)
    assert result.target_schema.mode == "binary"


def test_failed_fit_preserves_previous_output(tmp_path: Path) -> None:
    source, mapping = _source(tmp_path, "classification")
    # One training feature is negative, violating multinomial NB's input contract.
    with source.open() as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["x0"] = "-1"
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = tmp_path / "model"
    output.mkdir()
    sentinel = output / "old.txt"
    sentinel.write_text("existing user result")
    with pytest.raises(InputFormatError) as error:
        train_custom_property_model(
            source,
            config=CustomPropertyTrainingConfig(
                task_type="classification",
                model="multinomial_nb",
                finetune_method="machine_learning",
                feature_source="numeric",
                mapping=mapping,
                output_dir=output,
                overwrite=True,
            ),
        )
    assert error.value.code == "NEGATIVE_MACHINE_LEARNING_FEATURE"
    assert sentinel.read_text() == "existing user result"
    assert list(output.iterdir()) == [sentinel]


def test_calibration_requires_enough_training_samples(tmp_path: Path) -> None:
    source = tmp_path / "tiny.csv"
    source.write_text("x,target,split\n1,a,train\n2,b,train\n3,a,validation\n4,b,test\n")
    with pytest.raises(InputFormatError) as error:
        train_custom_property_model(
            source,
            config=CustomPropertyTrainingConfig(
                task_type="classification",
                model="linear_svm",
                finetune_method="machine_learning",
                feature_source="numeric",
                mapping=DataMapping(
                    sequence_field=None,
                    feature_fields=("x",),
                    label_fields=("target",),
                    split_field="split",
                ),
                output_dir=tmp_path / "model",
            ),
        )
    assert error.value.code == "INSUFFICIENT_CALIBRATION_SAMPLES"


def test_sequence_features_handle_ambiguity(tmp_path: Path) -> None:
    from dnakit.predictions import load_property_data
    from dnakit.predictions.custom_training.ml_features import machine_learning_features

    source = tmp_path / "sequence.csv"
    source.write_text("sequence,target\nANN,1\nNNN,2\nACG,3\n")
    config = CustomPropertyTrainingConfig(
        task_type="regression",
        model="ridge",
        finetune_method="machine_learning",
        kmer_k=2,
        mapping=DataMapping(label_fields=("target",)),
    )
    values, schema = machine_learning_features(
        load_property_data(source, mapping=config.mapping), config, None
    )
    assert values.shape == (3, 19)
    np.testing.assert_allclose(values[0, :3], [3, 0, 2 / 3])
    assert values[1, 3:].sum() == 0
    names = schema["names"]
    assert isinstance(names, list)
    assert values[2, names.index("kmer_AC")] == 0.5
    assert values[2, names.index("kmer_CG")] == 0.5


def test_missing_saved_embedding_path_uses_on_demand_resolver(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: Any
) -> None:
    from dnakit.predictions.custom_training import ml_training

    source, _ = _source(tmp_path, "regression")
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="regression",
            model="ridge",
            finetune_method="machine_learning",
            feature_source="embedding",
            mapping=DataMapping(label_fields=("target",), split_field="split", id_field="id"),
            output_dir=tmp_path / "model",
        ),
        representation_backend=_Embedding(),
        verbose=False,
    )
    manifest_path = Path(result.output_dir) / "model.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["features"]["checkpoint_path"] = str(tmp_path / "old-computer" / "weights")
    manifest_path.write_text(json.dumps(manifest))
    original = ml_training.machine_learning_features
    observed = []

    def features(data: Any, config: Any, backend: Any) -> Any:
        observed.append(config.checkpoint_path)
        return original(data, config, backend)

    monkeypatch.setattr(ml_training, "machine_learning_features", features)
    predicted = predict_custom_properties(
        source,
        model_dir=result.output_dir,
        trust_model=True,
        representation_backend=_Embedding(),
    )
    assert len(predicted.predictions) == 90
    assert observed == [None]  # Allows the normal backbone resolver to download.
    assert str(Path(result.output_dir) / "estimator.joblib") in capsys.readouterr().err
