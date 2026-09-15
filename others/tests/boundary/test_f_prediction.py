"""F4: model contracts with real encoding and controlled, offline model doubles."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    DownloadError,
    SequenceError,
)
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PredictionOutput,
    PropertyPredictionConfig,
    predict_properties,
)
from dnakit.predictions.backends import _EnformerBenchmarkBackend, _resolve_torch_device
from dnakit.predictions.enformer_benchmarks import get_enformer_benchmark_task
from dnakit.representations import RepresentationConfig, extract_representations


class FixedBackend:
    def __init__(self, values: Any) -> None:
        self.values = values

    def predict(self, inputs: Any, *, show_progress: Any) -> Any:
        return self.values


@pytest.fixture
def benchmark() -> Any:
    torch = pytest.importorskip("torch")
    backend = object.__new__(_EnformerBenchmarkBackend)
    backend.torch = torch
    backend.device = torch.device("cpu")
    backend.config = PropertyPredictionConfig(
        model="enformer", task="enhancers", batch_size=32, show_progress=False
    )
    backend.spec = get_enformer_benchmark_task("enhancers")
    backend.labels = ("0", "1")
    backend.checkpoint_metadata = {
        "format": "test-double",
        "saved_at_utc": "",
        "fold_id_one_based": 1,
    }
    backend.model = lambda one_hot, lengths: torch.stack(
        (lengths.float() / 100, -lengths.float() / 100), dim=1
    )
    return backend


@pytest.mark.parametrize("length", [1, 199, 200, 201, 196607, 196608])
def test_enformer_length_boundaries_without_loading_weights(benchmark: Any, length: Any) -> None:
    result = predict_properties(
        [BiologicalSequence("sample", "A" * length)],
        config=benchmark.config,
        backend=benchmark,
        verbose=False,
    )
    output = result.records[0].output
    assert output.metadata["sequence_length"] == length
    assert output.metadata["valid_bins"] == (length + 127) // 128
    assert np.isfinite(output.values).all()
    assert output.values.sum() == pytest.approx(1)


def test_enformer_rejects_model_length_overflow(benchmark: Any) -> None:
    with pytest.raises(SequenceError, match="PREDICTION_SEQUENCE_TOO_LONG") as error:
        predict_properties(
            [BiologicalSequence("long", "A" * 196609)],
            config=benchmark.config,
            backend=benchmark,
            verbose=False,
        )
    assert error.value.context["record_id"] == "long"


def test_empty_sequence_reports_record_id() -> None:
    with pytest.raises(SequenceError, match="EMPTY_PREDICTION_INPUT") as error:
        BiologicalSequence("empty_record", "")
    assert error.value.context["record_id"] == "empty_record"


@pytest.mark.parametrize("symbols", ["ACXT", "AC1T", "AC T"])
def test_invalid_prediction_symbols_identify_record(symbols: Any) -> None:
    with pytest.raises(SequenceError, match="INVALID_PREDICTION_INPUT_SYMBOL") as error:
        BiologicalSequence("record_51", symbols)
    assert error.value.context["record_id"] == "record_51"


@pytest.mark.parametrize("symbols", ["N" * 200, "ARYN", "A" * 200])
def test_low_information_sequences_follow_declared_policy(benchmark: Any, symbols: Any) -> None:
    inputs = [BiologicalSequence("ambiguous", symbols)]
    output = predict_properties(inputs, config=benchmark.config, backend=benchmark, verbose=False)
    assert output.record_ids == ("ambiguous",)
    if set(symbols) - set("ACGT"):
        benchmark.config = PropertyPredictionConfig(
            model="enformer", task="enhancers", ambiguity_policy="error"
        )
        with pytest.raises(SequenceError, match="MODEL_AMBIGUOUS_INPUT"):
            predict_properties(inputs, config=benchmark.config, backend=benchmark, verbose=False)


def test_100_inputs_tail_batch_and_id_order(benchmark: Any) -> None:
    # Unequal lengths force the real Enformer adapter to sort then restore order.
    inputs = [BiologicalSequence(f"sample_{i}", "A" * (201 - i % 3)) for i in range(100)]
    batch_sizes: list[int] = []
    model = benchmark.model

    def recording_model(one_hot: Any, lengths: Any) -> Any:
        batch_sizes.append(len(lengths))
        return model(one_hot, lengths)

    benchmark.model = recording_model
    result = predict_properties(inputs, config=benchmark.config, backend=benchmark, verbose=False)
    assert batch_sizes == [32, 32, 32, 4]
    assert result.record_ids == tuple(item.id for item in inputs)
    assert [r.output.metadata["sequence_length"] for r in result.records] == [
        len(i.sequence) for i in inputs
    ]


def test_illegal_middle_record_prevents_inference(benchmark: Any) -> None:
    def inputs() -> Any:
        for i in range(100):
            yield BiologicalSequence(f"sample_{i}", "ACXT" if i == 50 else "A" * 200)

    calls: list[int] = []
    benchmark.model = lambda *args: calls.append(1)
    with pytest.raises(SequenceError) as error:
        predict_properties(inputs(), config=benchmark.config, backend=benchmark, verbose=False)
    assert error.value.context["record_id"] == "sample_50"
    assert calls == []


@pytest.mark.parametrize("values", [[np.nan], [np.inf], [-np.inf], []])
def test_invalid_numeric_model_outputs(values: Any) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_PREDICTION_OUTPUT"):
        PredictionOutput(values)


@pytest.mark.parametrize(
    "task,model",
    [
        ("enhancers", "enformer"),
        ("genomic_segmentation", "segmentnt"),
        ("exon_probability", "evo2"),
    ],
)
@pytest.mark.parametrize("value", [-0.1, 1.1])
def test_probability_range_is_checked_across_prediction_models(
    task: Any, model: Any, value: Any
) -> None:
    config = PropertyPredictionConfig(model=model, task=task)
    sequence = BiologicalSequence("sample", "ACGT")
    item = BiologicalSequencePair("pair", sequence, sequence) if model == "evo2" else sequence
    with pytest.raises(BackendExecutionError, match="INVALID_PREDICTION_OUTPUT"):
        predict_properties(
            [item], config=config, backend=FixedBackend([PredictionOutput([value])]), verbose=False
        )


def test_missing_checkpoint_fails_without_fake_predictions(tmp_path: Any, monkeypatch: Any) -> None:
    gdown = pytest.importorskip("gdown")

    def failed_download(*args: Any, **kwargs: Any) -> Any:
        raise ConnectionError("Simulated model server failure")

    monkeypatch.setattr(gdown, "download", failed_download)
    config = PropertyPredictionConfig(
        model="enformer", task="enhancers", checkpoint_path=tmp_path / "missing.ckpt"
    )
    with pytest.raises(DownloadError, match="MODEL_CHECKPOINT_DOWNLOAD_FAILED"):
        predict_properties([BiologicalSequence("sample", "ACGT")], config=config, verbose=False)


@pytest.mark.parametrize("payload", [b"", b"corrupt checkpoint bytes"])
def test_corrupt_checkpoint_is_rejected(benchmark: Any, tmp_path: Any, payload: Any) -> None:
    benchmark.checkpoint = tmp_path / "invalid.ckpt"
    benchmark.checkpoint.write_bytes(payload)
    with pytest.raises(BackendExecutionError, match=r"INVALID_MODEL_CHECKPOINT|MODEL_LOAD_FAILED"):
        benchmark._load_checkpoint()


def test_checkpoint_task_mismatch(benchmark: Any) -> None:
    from dnakit.predictions.backends import _ENFORMER_BENCHMARK_FORMAT

    benchmark.payload = {
        "format": _ENFORMER_BENCHMARK_FORMAT,
        "family": benchmark.spec.family,
        "task": "wrong-task",
    }
    with pytest.raises(BackendExecutionError, match="MODEL_CHECKPOINT_TASK_MISMATCH"):
        benchmark._validate_metadata()


def test_unavailable_cuda_is_rejected() -> None:
    torch = pytest.importorskip("torch")
    fake = SimpleNamespace(device=torch.device, cuda=SimpleNamespace(is_available=lambda: False))
    with pytest.raises(ConfigurationError, match="INVALID_PREDICTION_DEVICE"):
        _resolve_torch_device(fake, "cuda")


@pytest.mark.parametrize("stage", ["encode", "inference"])
def test_simulated_enformer_oom_is_actionable(benchmark: Any, stage: Any) -> None:
    def oom(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError("CUDA out of memory")

    if stage == "encode":
        benchmark._encode = oom
    else:
        benchmark.model = oom
    with pytest.raises(BackendExecutionError, match="MODEL_OUT_OF_MEMORY") as error:
        predict_properties(
            [BiologicalSequence("sample", "ACGT")],
            config=benchmark.config,
            backend=benchmark,
            verbose=False,
        )
    assert error.value.hint is not None and "batch_size" in error.value.hint


@pytest.mark.parametrize("kind", ["prediction", "embedding"])
def test_shared_neural_api_reports_memory_failure(kind: Any) -> None:
    def oom(*args: Any, **kwargs: Any) -> Any:
        raise MemoryError("model allocation")

    with pytest.raises(BackendExecutionError, match="MODEL_OUT_OF_MEMORY"):
        if kind == "prediction":
            predict_properties(
                [BiologicalSequence("sample", "ACGT")],
                config=PropertyPredictionConfig(model="enformer", task="enhancers"),
                backend=SimpleNamespace(predict=oom),
                verbose=False,
            )
        else:
            extract_representations(
                [DNARecord(DNASequence("ACGT"), "sample")],
                config=RepresentationConfig(),
                backend=SimpleNamespace(extract=oom),
                verbose=False,
            )


@pytest.mark.parametrize("count", [0, 2])
def test_wrong_output_count_is_rejected(count: Any) -> None:
    with pytest.raises(BackendExecutionError, match="INVALID_PREDICTION_OUTPUT"):
        predict_properties(
            [BiologicalSequence("sample", "ACGT")],
            config=PropertyPredictionConfig(model="enformer", task="enhancers"),
            backend=FixedBackend([PredictionOutput([0.5])] * count),
            verbose=False,
        )


@pytest.mark.parametrize("values", [[1 + 2j], [10**1000]])
def test_complex_and_unrepresentable_model_values_are_rejected(values: Any) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_PREDICTION_OUTPUT"):
        PredictionOutput(values)


def test_duplicate_prediction_names_cannot_overwrite_csv_columns() -> None:
    with pytest.raises(ConfigurationError, match="unique"):
        PredictionOutput([0.2, 0.8], ("same", "same"))


@pytest.mark.parametrize("values", [[0.9, 0.9], [0.5], [[0.2, 0.8]]])
def test_enformer_class_probabilities_match_output_contract(values: Any) -> None:
    with pytest.raises(BackendExecutionError, match="INVALID_PREDICTION_OUTPUT"):
        predict_properties(
            [BiologicalSequence("sample", "ACGT")],
            config=PropertyPredictionConfig(model="enformer", task="enhancers"),
            backend=FixedBackend([PredictionOutput(values)]),
            verbose=False,
        )


@pytest.mark.parametrize("kind", ["prediction", "embedding"])
def test_invalid_cuda_index_is_rejected_before_loading(kind: Any) -> None:
    from dnakit.representations.backends import _resolve_torch_device as embedding_device

    torch = pytest.importorskip("torch")
    fake = SimpleNamespace(
        device=torch.device, cuda=SimpleNamespace(is_available=lambda: True, device_count=lambda: 1)
    )
    resolver = _resolve_torch_device if kind == "prediction" else embedding_device
    with pytest.raises(ConfigurationError, match="device index"):
        resolver(fake, "cuda:9")


def test_custom_length_limit_cannot_override_architecture_limit(benchmark: Any) -> None:
    benchmark.config = PropertyPredictionConfig(
        model="enformer", task="enhancers", max_length=200000
    )
    with pytest.raises(SequenceError, match="PREDICTION_SEQUENCE_TOO_LONG"):
        predict_properties(
            [BiologicalSequence("long", "A" * 196609)],
            config=benchmark.config,
            backend=benchmark,
            verbose=False,
        )


@pytest.mark.parametrize("kind", ["prediction", "embedding"])
def test_wrapped_backend_oom_remains_actionable(kind: Any) -> None:
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        try:
            raise RuntimeError("CUDA out of memory")
        except RuntimeError as cause:
            raise BackendExecutionError(
                "Third-party model failed", code="MODEL_LOAD_FAILED"
            ) from cause

    with pytest.raises(BackendExecutionError, match="MODEL_OUT_OF_MEMORY"):
        if kind == "prediction":
            predict_properties(
                [BiologicalSequence("sample", "ACGT")],
                config=PropertyPredictionConfig(model="enformer", task="enhancers"),
                backend=SimpleNamespace(predict=wrapped),
                verbose=False,
            )
        else:
            extract_representations(
                [DNARecord(DNASequence("ACGT"), "sample")],
                config=RepresentationConfig(),
                backend=SimpleNamespace(extract=wrapped),
                verbose=False,
            )
