"""Tests for user-defined classification/regression training and schema ingestion."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from contextlib import closing
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from numpy.typing import NDArray

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    DataMapping,
    LLMSchemaConfig,
    OpenAICompatibleSchemaClient,
    available_custom_training_models,
    compute_metrics,
    inspect_property_data,
    load_property_data,
    predict_custom_properties,
    train_custom_property_model,
)
from dnakit.predictions.custom_training._loading import _PropertyBatchLoader
from dnakit.predictions.custom_training.llm import mapping_from_llm_proposal
from dnakit.predictions.custom_training.targets import encode_targets, fit_target_schema


class _FakeSchemaClient:
    def __init__(self, proposal: Mapping[str, object]) -> None:
        self.proposal = proposal
        self.profile: Mapping[str, object] = {}

    def analyze(
        self,
        profile: Mapping[str, object],
        *,
        config: LLMSchemaConfig,
    ) -> Mapping[str, object]:
        del config
        self.profile = profile
        return self.proposal


class _FakeRepresentationBackend:
    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> NDArray[np.float32]:
        del show_progress
        return np.asarray(
            [[len(sequence), sequence.count("G"), sequence.count("C")] for sequence in sequences],
            dtype=np.float32,
        )


class _FakeHTTPResponse:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload

    def __enter__(self) -> _FakeHTTPResponse:
        return self

    def __exit__(self, *args: object) -> None:
        del args

    def read(self, maximum: int = -1) -> bytes:
        return self.payload if maximum < 0 else self.payload[:maximum]


def test_custom_training_capabilities_and_secret_repr() -> None:
    capabilities = {item.name: item for item in available_custom_training_models()}
    assert capabilities["dnabert2"].full_finetune is True
    assert capabilities["alphagenome"].full_finetune is False
    assert all(item.linear_probe for item in capabilities.values())

    llm = LLMSchemaConfig(api_key="do-not-display")
    assert "do-not-display" not in repr(llm)
    assert llm.endpoint == "https://api.deepseek.com/chat/completions"
    with pytest.raises(ConfigurationError) as unsupported:
        CustomPropertyTrainingConfig(
            task_type="classification",
            model="alphagenome",
            finetune_method="full",
        )
    assert unsupported.value.code == "CUSTOM_FULL_FINETUNE_UNSUPPORTED"
    with pytest.raises(ConfigurationError):
        LLMSchemaConfig(api_key="unsafe\nheader")
    with pytest.raises(ConfigurationError):
        DataMapping(delimiter="\n")


@pytest.mark.parametrize("method", ["full", "linear_probe"])
def test_custom_training_optimizer_and_loader_defaults(method: Any) -> None:
    config = CustomPropertyTrainingConfig(task_type="regression", finetune_method=method)
    assert config.workers == 6
    assert config.batch_size == 32
    assert config.learning_rate == config.resolved_learning_rate == 1e-4
    assert (config.beta1, config.beta2) == (0.9, 0.95)
    assert config.weight_decay == 0.1
    assert config.seed == 42
    assert (config.train_ratio, config.validation_ratio, config.test_ratio) == (0.8, 0.1, 0.1)
    assert (
        CustomPropertyTrainingConfig(
            task_type="regression", finetune_method=method, learning_rate=None
        ).resolved_learning_rate
        == 1e-4
    )


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("workers", -1),
        ("workers", 257),
        ("workers", 1.5),
        ("workers", True),
        ("workers", None),
        ("batch_size", 0),
        ("batch_size", 1.5),
        ("beta1", -0.1),
        ("beta1", 1.0),
        ("beta1", True),
        ("beta1", "0.9"),
        ("beta2", float("nan")),
        ("beta2", float("inf")),
        ("beta2", None),
        ("learning_rate", 0.0),
        ("learning_rate", float("inf")),
        ("weight_decay", -0.1),
        ("weight_decay", float("nan")),
    ],
)
def test_invalid_training_hyperparameters_are_rejected(option: str, value: Any) -> None:
    with pytest.raises(ConfigurationError, match=option):
        CustomPropertyTrainingConfig(task_type="regression", **{option: value})


def test_zero_workers_betas_and_weight_decay_are_supported() -> None:
    config = CustomPropertyTrainingConfig(
        task_type="regression", workers=0, beta1=0.0, beta2=0.0, weight_decay=0.0
    )
    assert config.workers == config.beta1 == config.beta2 == config.weight_decay == 0


@pytest.mark.parametrize("workers", [0, 2])
def test_seeded_batch_loading_and_worker_cleanup(workers: int) -> None:
    torch = pytest.importorskip("torch")
    sequences = tuple("A" * (index + 1) for index in range(10))
    targets = np.arange(20, dtype=np.float32).reshape(10, 2)
    representations = targets + 10
    processes: list[Any] = []
    with (
        pytest.raises(RuntimeError, match="test cleanup"),
        closing(
            _PropertyBatchLoader(
                torch,
                sequences,
                targets,
                representations=representations,
                batch_size=4,
                workers=workers,
                seed=42,
            )
        ) as loader,
    ):
        assert loader.loader.num_workers == workers
        for seed in (43, 44, 43):
            batches = list(loader.batches(tuple(range(10)), shuffle=True, seed=seed))
            assert [len(batch["targets"]) for batch in batches] == [4, 4, 2]
            order = torch.cat([batch["row_index"] for batch in batches]).tolist()
            assert order == np.random.default_rng(seed).permutation(10).tolist()
            for batch in batches:
                indices = batch["row_index"].tolist()
                assert batch["sequence"] == [sequences[index] for index in indices]
                assert batch["targets"].device.type == "cpu"
                assert batch["representations"].device.type == "cpu"
                np.testing.assert_array_equal(batch["targets"].numpy(), targets[indices])
                np.testing.assert_array_equal(
                    batch["representations"].numpy(), representations[indices]
                )
        validation = list(loader.batches((9, 5, 2), shuffle=False, seed=0))
        assert validation[0]["row_index"].tolist() == [9, 5, 2]
        if workers:
            processes = list(loader.loader._iterator._workers)
            assert len(processes) == workers
            assert all(process.is_alive() and process.pid != os.getpid() for process in processes)
        raise RuntimeError("test cleanup")
    assert all(not process.is_alive() for process in processes)


def test_csv_auto_mapping_and_seeded_default_split(tmp_path: Path) -> None:
    source = tmp_path / "labels.csv"
    source.write_text(
        "index,seq,label\n"
        + "".join(
            f"sample-{index},ACGTACG{'A' if index % 2 == 0 else 'T'},{index % 2}\n"
            for index in range(20)
        ),
        encoding="utf-8",
    )
    inspection = inspect_property_data(source)
    assert inspection.detected_format == "csv"
    assert inspection.mapping.sequence_field == "seq"
    assert inspection.mapping.id_field == "index"
    assert inspection.mapping.label_fields == ("label",)
    assert inspection.estimated_records == 20

    data = load_property_data(source, inspection=inspection)
    assert data.record_count == 20
    assert data.record_ids[0] == "sample-0"
    assert data.label_names == ("label",)


def test_nested_json_and_fasta_header_mapping(tmp_path: Path) -> None:
    nested = tmp_path / "nested.json"
    nested.write_text(
        json.dumps(
            {
                "experiment": {
                    "samples": [
                        {"sample": "a", "input": {"dna": "ACGT"}, "targets": {"score": 1.5}},
                        {"sample": "b", "input": {"dna": "TGCA"}, "targets": {"score": 2.5}},
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    mapping = DataMapping(
        format="json",
        record_path="experiment.samples",
        id_field="sample",
        sequence_field="input.dna",
        label_fields=("targets.score",),
    )
    data = load_property_data(nested, mapping=mapping)
    assert data.sequences == ("ACGT", "TGCA")
    assert data.labels == ((1.5,), (2.5,))

    fasta = tmp_path / "labeled.fa"
    fasta.write_text(
        ">r1 class=promoter|split=train\nACGT\n>r2 class=enhancer|split=test\nTGCA\n",
        encoding="utf-8",
    )
    fasta_mapping = DataMapping(
        format="fasta",
        id_field="id",
        sequence_field="sequence",
        label_fields=("header.class",),
        split_field="header.split",
    )
    fasta_data = load_property_data(fasta, mapping=fasta_mapping)
    assert fasta_data.labels == (("promoter",), ("enhancer",))
    assert fasta_data.supplied_splits == ("train", "test")


def test_llm_receives_only_bounded_profile_and_mapping_is_validated(tmp_path: Path) -> None:
    source = tmp_path / "unusual.records"
    source.write_text(
        "sample|bases|response\n" + "\n".join(f"r{i}|ACGTACGT|{i}" for i in range(50)),
        encoding="utf-8",
    )
    client = _FakeSchemaClient(
        {
            "format": "delimited",
            "sequence_field": "bases",
            "label_fields": ["response"],
            "id_field": "sample",
            "split_field": None,
            "record_path": None,
            "delimiter": "|",
            "header_row": 0,
            "sheet_name": None,
            "table_name": None,
            "fasta_header_separator": "|",
            "recommendations": ["Verify response units."],
            "warnings": [],
            "confidence": 0.98,
        }
    )
    llm = LLMSchemaConfig(
        api_key="test-only",
        sample_rows=3,
        max_sample_characters=256,
    )
    inspection = inspect_property_data(source, llm=llm, llm_client=client)
    assert inspection.mapping_origin == "llm"
    assert inspection.mapping.delimiter == "|"
    assert inspection.shared_sample_records <= 3
    policy = client.profile["sample_policy"]
    assert isinstance(policy, Mapping)
    assert policy["whole_file_shared"] is False
    assert "test-only" not in json.dumps(client.profile)

    malicious = _FakeSchemaClient({"format": "python", "sequence_field": "exec(payload)"})
    with pytest.raises(ConfigurationError):
        inspect_property_data(source, llm=llm, llm_client=malicious)

    small_source = tmp_path / "small.csv"
    small_source.write_text(
        "sample,bases,response\nr1,ACGT,1\nr2,TGCA,withheld-value\n",
        encoding="utf-8",
    )
    small_client = _FakeSchemaClient(client.proposal)
    small_inspection = inspect_property_data(
        small_source,
        llm=LLMSchemaConfig(api_key="test-only", sample_rows=8),
        llm_client=small_client,
    )
    assert small_inspection.shared_sample_records == 1
    assert "withheld-value" not in json.dumps(small_client.profile)
    assert small_client.profile["field_types"] == {
        "sample": ["str"],
        "bases": ["str"],
        "response": ["str"],
    }


def test_openai_compatible_request_and_strict_mapping_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeHTTPResponse:
        captured["url"] = request.full_url
        captured["authorization"] = request.get_header("Authorization")
        captured["payload"] = json.loads(request.data.decode("utf-8"))
        captured["timeout"] = timeout
        proposal = {
            "format": "csv",
            "sequence_field": "bases",
            "label_fields": ["score"],
            "id_field": None,
            "split_field": None,
            "record_path": None,
            "delimiter": None,
            "header_row": 0,
            "sheet_name": None,
            "table_name": None,
            "fasta_header_separator": "|",
        }
        envelope = {"choices": [{"message": {"content": json.dumps(proposal)}}]}
        return _FakeHTTPResponse(json.dumps(envelope).encode())

    monkeypatch.setattr("dnakit.predictions.custom_training.llm._open_request", fake_urlopen)
    config = LLMSchemaConfig(api_key="request-only-secret", timeout_seconds=12)
    proposal = OpenAICompatibleSchemaClient().analyze(
        {"columns_or_paths": ["sample", "bases", "score"]}, config=config
    )
    mapping, _, _ = mapping_from_llm_proposal(
        proposal,
        fallback=DataMapping(id_field="sample", delimiter=","),
    )
    assert captured["url"] == "https://api.deepseek.com/chat/completions"
    assert captured["authorization"] == "Bearer request-only-secret"
    assert captured["timeout"] == 12.0
    request_payload = captured["payload"]
    assert isinstance(request_payload, Mapping)
    assert request_payload["response_format"] == {"type": "json_object"}
    assert "request-only-secret" not in json.dumps(request_payload)
    assert mapping.id_field is None
    assert mapping.delimiter is None

    with pytest.raises(ConfigurationError):
        mapping_from_llm_proposal(
            {"sequence_field": None},
            fallback=DataMapping(),
        )


def test_classification_and_regression_metrics_are_task_appropriate() -> None:
    labels = [("negative",), ("positive",), ("negative",), ("positive",)]
    classification = fit_target_schema(
        labels,
        ("activity",),
        task_type="classification",
        classification_mode="auto",
        train_indices=(0, 1, 2, 3),
    )
    encoded = encode_targets(labels, schema=classification)
    classification_metrics = compute_metrics(
        encoded,
        np.asarray([[3.0, 0.1], [0.2, 2.0], [2.0, 0.1], [0.1, 2.5]]),
        schema=classification,
    )
    assert classification_metrics["accuracy"] == pytest.approx(1.0)
    assert classification_metrics["f1_macro"] == pytest.approx(1.0)
    assert "auroc" in classification_metrics

    regression_labels = [(1.0, 10.0), (2.0, 20.0), (3.0, 30.0), (4.0, 40.0)]
    regression = fit_target_schema(
        regression_labels,
        ("yield", "temperature"),
        task_type="regression",
        classification_mode="auto",
        train_indices=(0, 1, 2, 3),
    )
    regression_targets = encode_targets(regression_labels, schema=regression)
    regression_metrics = compute_metrics(
        regression_targets,
        regression_targets,
        schema=regression,
    )
    aggregate = regression_metrics["aggregate"]
    assert isinstance(aggregate, Mapping)
    assert aggregate["rmse"] == pytest.approx(0.0)
    assert aggregate["r2"] == pytest.approx(1.0)
    assert "f1" not in aggregate


def test_one_hot_multiclass_multilabel_and_multioutput_schemas() -> None:
    one_hot = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0)]
    one_hot_schema = fit_target_schema(
        one_hot,
        ("promoter", "enhancer", "other"),
        task_type="classification",
        classification_mode="auto",
        train_indices=(0, 1, 2, 3),
    )
    assert one_hot_schema.mode == "multiclass"
    assert one_hot_schema.output_dimension == 3

    multilabel = [(1, 0), (1, 1), (0, 1), (0, 0)]
    multilabel_schema = fit_target_schema(
        multilabel,
        ("active", "conserved"),
        task_type="classification",
        classification_mode="auto",
        train_indices=(0, 1, 2, 3),
    )
    assert multilabel_schema.mode == "multilabel"
    encoded_multilabel = encode_targets(multilabel, schema=multilabel_schema)
    multilabel_metrics = compute_metrics(
        encoded_multilabel,
        np.where(encoded_multilabel == 1, 8.0, -8.0),
        schema=multilabel_schema,
    )
    assert multilabel_metrics["subset_accuracy"] == pytest.approx(1.0)
    assert multilabel_metrics["f1_micro"] == pytest.approx(1.0)

    multioutput = [("low", "A"), ("high", "B"), ("low", "B"), ("high", "A")]
    multioutput_schema = fit_target_schema(
        multioutput,
        ("activity", "family"),
        task_type="classification",
        classification_mode="auto",
        train_indices=(0, 1, 2, 3),
    )
    assert multioutput_schema.mode == "multioutput"
    assert multioutput_schema.output_dimension == 4


@pytest.mark.parametrize("workers", [0, 2])
def test_linear_probe_training_artifacts_and_saved_inference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, workers: int
) -> None:
    torch = pytest.importorskip("torch")
    pytest.importorskip("sklearn")
    observed_optimizers: list[dict[str, Any]] = []
    adamw = torch.optim.AdamW

    def record_adamw(parameters: Any, **kwargs: Any) -> Any:
        observed_optimizers.append(kwargs)
        return adamw(parameters, **kwargs)

    monkeypatch.setattr(torch.optim, "AdamW", record_adamw)
    source = tmp_path / "training.csv"
    source.write_text(
        "id,sequence,label\n"
        + "".join(
            f"r{index},{'G' * (index % 6 + 1)}{'A' * (7 - index % 6)},{index % 2}\n"
            for index in range(40)
        ),
        encoding="utf-8",
    )
    output = tmp_path / "trained-model"
    config = CustomPropertyTrainingConfig(
        task_type="classification",
        model="dnabert2",
        finetune_method="linear_probe",
        output_dir=output,
        epochs=2,
        batch_size=7,
        workers=workers,
        learning_rate=3e-4,
        beta1=0.85,
        beta2=0.97,
        weight_decay=0.03,
        device="cpu",
        show_progress=False,
    )
    backend = _FakeRepresentationBackend()
    result = train_custom_property_model(
        source,
        config=config,
        representation_backend=backend,
    )
    assert dict(result.split_counts) == {"train": 32, "validation": 4, "test": 4}
    assert result.target_schema.mode == "binary"
    assert (output / "head.npz").is_file()
    assert (output / "metrics.json").is_file()
    assert (output / "splits.csv").is_file()
    manifest_text = (output / "model.json").read_text(encoding="utf-8")
    assert "api_key" not in manifest_text
    assert observed_optimizers == [{"lr": 3e-4, "betas": (0.85, 0.97), "weight_decay": 0.03}]
    settings = json.loads(manifest_text)["training"]
    assert settings["workers"] == workers
    assert settings["batch_size"] == 7
    assert settings["batch_size_scope"] == "per_device"
    assert settings["optimizer"] == "AdamW"
    assert settings["learning_rate"] == 3e-4
    assert (settings["beta1"], settings["beta2"]) == (0.85, 0.97)
    assert settings["weight_decay"] == 0.03

    prediction = predict_custom_properties(
        source,
        model_dir=output,
        representation_backend=backend,
        device="cpu",
        show_progress=False,
    )
    assert prediction.record_ids == tuple(f"r{index}" for index in range(40))
    assert len(prediction.predictions) == 40


def test_tiny_transformers_full_finetune_and_reload(tmp_path: Path) -> None:
    pytest.importorskip("torch")
    transformers = pytest.importorskip("transformers")
    pytest.importorskip("sklearn")

    checkpoint = tmp_path / "tiny-bert"
    checkpoint.mkdir()
    vocabulary = checkpoint / "vocab.txt"
    vocabulary.write_text(
        "[PAD]\n[UNK]\n[CLS]\n[SEP]\n[MASK]\nA\nC\nG\nT\nN\n",
        encoding="utf-8",
    )
    tokenizer = transformers.BertTokenizer(vocab_file=str(vocabulary))
    tokenizer.save_pretrained(checkpoint)
    backbone = transformers.BertModel(
        transformers.BertConfig(
            vocab_size=10,
            hidden_size=8,
            num_hidden_layers=1,
            num_attention_heads=2,
            intermediate_size=16,
            max_position_embeddings=64,
        )
    )
    backbone.save_pretrained(checkpoint)

    source = tmp_path / "full-training.csv"
    source.write_text(
        "id,sequence,label\n"
        + "".join(
            f"r{index},{'ACGT' if index % 2 else 'TGCA'},{'active' if index % 2 else 'inactive'}\n"
            for index in range(30)
        ),
        encoding="utf-8",
    )
    output = tmp_path / "full-model"
    result = train_custom_property_model(
        source,
        config=CustomPropertyTrainingConfig(
            task_type="classification",
            model="grover",
            finetune_method="full",
            checkpoint_path=checkpoint,
            output_dir=output,
            epochs=1,
            batch_size=6,
            workers=0,
            device="cpu",
            show_progress=False,
        ),
    )
    assert result.finetune_method == "full"
    assert (output / "backbone" / "config.json").is_file()
    assert (output / "head.npz").is_file()
    settings = json.loads((output / "model.json").read_text(encoding="utf-8"))["training"]
    assert settings["learning_rate"] == 1e-4
    assert (settings["beta1"], settings["beta2"]) == (0.9, 0.95)
    assert settings["weight_decay"] == 0.1

    prediction = predict_custom_properties(
        source,
        model_dir=output,
        batch_size=6,
        device="cpu",
        show_progress=False,
    )
    assert len(prediction.predictions) == 30


def test_missing_or_invalid_mapped_fields_fail_before_training(tmp_path: Path) -> None:
    source = tmp_path / "bad.csv"
    source.write_text("id,text,label\na,not-dna,1\n", encoding="utf-8")
    with pytest.raises(InputFormatError) as error:
        load_property_data(
            source,
            mapping=DataMapping(
                format="csv",
                sequence_field="missing",
                label_fields=("label",),
            ),
        )
    assert error.value.code == "CUSTOM_DATA_FIELD_MISSING"


def test_xml_doctype_and_unsafe_training_output_are_rejected(tmp_path: Path) -> None:
    xml_source = tmp_path / "unsafe.xml"
    xml_source.write_text(
        '<!DOCTYPE records [<!ENTITY dna "ACGT">]>'
        "<records><row><sequence>&dna;</sequence><label>1</label></row></records>",
        encoding="utf-8",
    )
    with pytest.raises(InputFormatError) as xml_error:
        load_property_data(
            xml_source,
            mapping=DataMapping(
                format="xml",
                sequence_field="sequence",
                label_fields=("label",),
            ),
        )
    assert xml_error.value.code == "UNSAFE_CUSTOM_DATA_XML"

    csv_source = tmp_path / "source.csv"
    csv_source.write_text("sequence,label\nACGT,1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError) as output_error:
        train_custom_property_model(
            csv_source,
            config=CustomPropertyTrainingConfig(
                task_type="regression",
                finetune_method="linear_probe",
                output_dir=tmp_path,
            ),
        )
    assert output_error.value.code == "UNSAFE_CUSTOM_OUTPUT_DIR"
