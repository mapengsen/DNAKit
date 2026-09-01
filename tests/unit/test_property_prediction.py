"""Tests for the pretrained property-prediction registry and public contract."""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.exceptions import BackendExecutionError, ConfigurationError, SequenceError
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PredictionOutput,
    PropertyPredictionConfig,
    VariantContext,
    available_prediction_models,
    available_prediction_tasks,
    get_prediction_task,
    predict_properties,
)
from dnakit.predictions.backends import _LUCA_PRESETS, _LucaOneTasksBackend
from dnakit.predictions.models import PredictionInput


class _FakeBackend:
    def __init__(self, outputs: Sequence[PredictionOutput]) -> None:
        self.outputs = tuple(outputs)
        self.inputs: tuple[PredictionInput, ...] = ()
        self.show_progress: bool | None = None

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        self.inputs = tuple(inputs)
        self.show_progress = show_progress
        return self.outputs


def test_direct_prediction_registry_excludes_embedding_only_models() -> None:
    assert available_prediction_models() == (
        "alphagenome",
        "enformer",
        "evo2",
        "generator",
        "lucaone",
        "segmentnt",
    )
    assert len(available_prediction_tasks("alphagenome")) == 11
    assert len(available_prediction_tasks("lucaone")) == 10
    assert available_prediction_tasks("evo2") == ("exon_probability", "variant_effect")
    assert "segmentnt:genomic_segmentation" in available_prediction_tasks()
    assert get_prediction_task("segment-nt", "segmentation").name == "genomic_segmentation"

    for embedding_only in ("dnabert2", "ntv2", "hyenadna", "caduceus", "janusdna"):
        with pytest.raises(ConfigurationError) as error:
            available_prediction_tasks(embedding_only)
        assert error.value.code == "INVALID_PREDICTION_MODEL"


def test_prediction_inputs_normalize_and_validate_variant_context() -> None:
    sequence = BiologicalSequence(" gene-1 ", "acgu", "gene")
    assert sequence.id == "gene-1"
    assert sequence.sequence == "ACGU"

    variant = VariantContext("v1", "AACG", "AATG")
    assert variant.variant_index == 2
    assert variant.reference_base == "C"
    assert variant.alternate_base == "T"

    with pytest.raises(SequenceError) as multiple:
        VariantContext("bad", "AAAA", "ACCA")
    assert multiple.value.code == "INVALID_VARIANT_CONTEXT"

    with pytest.raises(SequenceError) as noncanonical:
        VariantContext("bad", "AANA", "AACA")
    assert noncanonical.value.code == "INVALID_VARIANT_CONTEXT"


def test_custom_backend_returns_read_only_outputs_for_dna_records() -> None:
    backend = _FakeBackend(
        [
            PredictionOutput(
                [[0.25, 0.75], [0.9, 0.1]],
                ("exon", "intron"),
                {"axes": ("position", "feature")},
            )
        ]
    )
    config = PropertyPredictionConfig(
        model="segmentnt",
        task="genomic_segmentation",
        show_progress=False,
    )
    result = predict_properties(
        [DNARecord(DNASequence("ACGT"), "record-1")],
        config=config,
        backend=backend,
    )

    assert backend.show_progress is False
    assert isinstance(backend.inputs[0], BiologicalSequence)
    assert result.record_ids == ("record-1",)
    assert result.metadata["fine_tuning_performed"] is False
    assert result.records[0].output.values.shape == (2, 2)
    assert not result.records[0].output.values.flags.writeable
    assert result.to_dict(include_values=False)["records"][0]["output"]["shape"] == (2, 2)  # type: ignore[index]


def test_variant_and_pair_task_input_contracts_are_explicit() -> None:
    variant = VariantContext("v1", "AACG", "AATG")
    variant_backend = _FakeBackend(
        [
            PredictionOutput(
                (-1.0, -1.2, -0.2),
                ("reference", "alternate", "delta"),
            )
        ]
    )
    result = predict_properties(
        [variant],
        config=PropertyPredictionConfig(model="evo2", task="variant_effect"),
        backend=variant_backend,
    )
    assert result.input_kind == "variant"
    assert result.records[0].output.values[2] == pytest.approx(-0.2)

    pair = BiologicalSequencePair(
        "position-1",
        BiologicalSequence("forward", "ACGT"),
        BiologicalSequence("reverse", "TGCA"),
    )
    pair_backend = _FakeBackend([PredictionOutput((0.8,), ("exon_probability",))])
    pair_result = predict_properties(
        [pair],
        config=PropertyPredictionConfig(model="evo2", task="exon_probability"),
        backend=pair_backend,
    )
    assert pair_result.input_kind == "pair"

    with pytest.raises(ConfigurationError) as mismatch:
        predict_properties(
            [BiologicalSequence("sequence", "ACGT")],
            config=PropertyPredictionConfig(model="evo2", task="variant_effect"),
            backend=_FakeBackend([PredictionOutput((0.0,), ("score",))]),
        )
    assert mismatch.value.code == "PREDICTION_INPUT_KIND_MISMATCH"


def test_remote_checkpoint_code_requires_opt_in_for_standard_backend() -> None:
    config = PropertyPredictionConfig(
        model="segmentnt",
        task="genomic_segmentation",
        allow_remote_code=False,
    )
    with pytest.raises(ConfigurationError) as error:
        predict_properties([BiologicalSequence("a", "ACGT")], config=config)
    assert error.value.code == "MODEL_REMOTE_CODE_NOT_ALLOWED"


def test_prediction_rejects_duplicate_ids_and_bad_backend_output() -> None:
    config = PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation")
    with pytest.raises(ConfigurationError) as duplicate:
        predict_properties(
            [BiologicalSequence("a", "ACGT"), BiologicalSequence("a", "TGCA")],
            config=config,
            backend=_FakeBackend(
                [PredictionOutput((0.1,), ("p",)), PredictionOutput((0.2,), ("p",))]
            ),
        )
    assert duplicate.value.code == "DUPLICATE_PREDICTION_ID"

    with pytest.raises(BackendExecutionError) as count_error:
        predict_properties(
            [BiologicalSequence("a", "ACGT")],
            config=config,
            backend=_FakeBackend([]),
        )
    assert getattr(count_error.value, "code", None) == "INVALID_PREDICTION_OUTPUT"

    array = np.asarray([[1.0, float("nan")]])
    with pytest.raises(ConfigurationError) as finite_error:
        PredictionOutput(array, ("a", "b"))
    assert finite_error.value.code == "INVALID_PREDICTION_OUTPUT"


def test_lucaone_csv_uses_official_prot_type_and_safe_ids(tmp_path: Path) -> None:
    backend = object.__new__(_LucaOneTasksBackend)
    backend.preset = _LUCA_PRESETS["protein_location"]
    output = tmp_path / "input.csv"

    backend._write_input(
        output,
        [BiologicalSequence("protein-1", "MKWV", "protein")],
    )
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows == [
        ["seq_id", "seq_type", "seq"],
        ["protein-1", "prot", "MKWV"],
    ]

    with pytest.raises(SequenceError) as unsafe:
        backend._write_input(
            output,
            [BiologicalSequence("../escape", "MKWV", "protein")],
        )
    assert unsafe.value.code == "INVALID_PREDICTION_INPUT"
