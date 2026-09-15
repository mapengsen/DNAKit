"""K-fold partitions match sklearn and remain fixed across validation rounds."""

from __future__ import annotations

import builtins
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import pytest
from typer.testing import CliRunner

from dnakit.core import DNARecord, DNASequence
from dnakit.datasets import FileSplitResult, SplitConfig, split
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.io import read_inputs

sklearn_selection = pytest.importorskip("sklearn.model_selection")


def _records(count: int) -> list[DNARecord]:
    return [DNARecord(DNASequence("ACGT"), f"r{index}") for index in range(count)]


@pytest.mark.parametrize("n_splits", [2, 5, 10])
@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("preserve_order", [False, True])
def test_all_folds_match_sklearn_and_cover_each_validation_sample_once(
    n_splits: int, shuffle: bool, preserve_order: bool
) -> None:
    records = _records(37)
    config = SplitConfig(
        method="kfold", n_splits=n_splits, shuffle=shuffle, preserve_order=preserve_order
    )
    reference = sklearn_selection.KFold(
        n_splits=n_splits, shuffle=shuffle, random_state=42 if shuffle else None
    )
    seen: set[str] = set()
    sizes = []
    for fold_id, (train_indices, valid_indices) in enumerate(reference.split(records)):
        result = split(records, config=replace(config, fold_id=fold_id), verbose=False)
        train, valid = result.get("train"), result.get("valid")
        assert train.ids == tuple(records[index].id for index in train_indices)
        assert valid.ids == tuple(records[index].id for index in valid_indices)
        assert not set(train.ids) & set(valid.ids)
        assert set(train.ids) | set(valid.ids) == {record.id for record in records}
        assert not seen & set(valid.ids)
        seen.update(valid.ids)
        sizes.append(len(valid))
        audit = json.loads(json.dumps(result.to_dict()))
        assert (audit["seed"], audit["n_splits"], audit["fold_id"]) == (42, n_splits, fold_id)
        assert audit["assignment_strategy"] == "sklearn_kfold"
    assert seen == {record.id for record in records}
    assert max(sizes) - min(sizes) <= 1
    # The selected fold does not depend on the order in which rounds are requested.
    again = split(records, config=replace(config, fold_id=n_splits - 1), verbose=False)
    assert again.assignments == result.assignments


def test_test_set_stays_outside_every_fold_and_seed_is_fixed() -> None:
    outer = split(
        _records(100), config=SplitConfig(ratios={"train": 0.8, "test": 0.2}), verbose=False
    )
    training_pool, test = outer.get("train"), outer.get("test")
    original_test = test.ids
    assert SplitConfig().seed == 42
    config = SplitConfig(method="kfold")
    assert (config.n_splits, config.fold_id, config.seed) == (10, 0, 42)
    for fold_id in range(10):
        result = split(training_pool, config=replace(config, fold_id=fold_id), verbose=False)
        assert result.counts == {"train": 72, "valid": 8}
        assert not set(test.ids) & (set(result.get("train").ids) | set(result.get("valid").ids))
        assert test.ids == original_test
    first = split(training_pool, config=config, verbose=False)
    changed = split(training_pool, config=replace(config, seed=43), verbose=False)
    assert first.get("valid").ids != changed.get("valid").ids


def test_unshuffled_folds_are_consecutive_and_ignore_seed() -> None:
    records = _records(11)
    config = SplitConfig(method="kfold", n_splits=3, fold_id=1, shuffle=False)
    first = split(records, config=config, verbose=False)
    second = split(records, config=replace(config, seed=19), verbose=False)
    assert first.get("valid").ids == ("r4", "r5", "r6", "r7")
    assert first.assignments == second.assignments


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_splits": 0},
        {"n_splits": 1},
        {"n_splits": True},
        {"n_splits": 3.5},
        {"n_splits": "10"},
        {"fold_id": -1},
        {"fold_id": 10},
        {"fold_id": True},
        {"fold_id": 0.5},
        {"seed": -1},
        {"seed": 2**32},
        {"seed": True},
        {"ratios": {"train": 0.8, "test": 0.2}},
        {"ratios": {"train": 0.5, "valid": 0.5}},
        {"shuffle": "yes"},
    ],
)
def test_invalid_kfold_configuration_is_rejected(kwargs: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        SplitConfig(method="kfold", **kwargs)


@pytest.mark.parametrize("method", ["random", "hash", "stratified", "group", "similarity"])
def test_fold_arguments_cannot_silently_fall_back_to_ratio_splitting(method: Any) -> None:
    with pytest.raises(ConfigurationError, match='require method="kfold"'):
        SplitConfig(method=method, n_splits=10, fold_id=0)


@pytest.mark.parametrize("count", [0, 1, 9])
def test_too_few_records_cannot_generate_empty_validation_folds(count: int) -> None:
    with pytest.raises(ConfigurationError) as error:
        split(_records(count), config=SplitConfig(method="kfold"), verbose=False)
    assert error.value.code == "KFOLD_SIZE_LIMIT"


def test_missing_dependency_has_installation_hint_and_ratio_split_still_works(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_import = builtins.__import__

    def missing_sklearn(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "sklearn.model_selection":
            raise ImportError("simulated missing sklearn")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", missing_sklearn)
    with pytest.raises(BackendUnavailableError) as error:
        split(_records(20), config=SplitConfig(method="kfold"), verbose=False)
    assert error.value.code == "MISSING_KFOLD_DEPENDENCY"
    assert "bio-analysis" in str(error.value.hint)
    assert split(_records(20), config=SplitConfig(), verbose=False).counts == {
        "train": 16,
        "valid": 2,
        "test": 2,
    }


@pytest.mark.parametrize("suffix", ["csv", "fasta"])
@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("preserve_order", [False, True])
def test_file_folds_match_memory_and_preserve_full_records(
    tmp_path: Path, suffix: str, shuffle: bool, preserve_order: bool
) -> None:
    source = tmp_path / f"input.{suffix}"
    content = (
        "id,sequence,label\n" + "".join(f"r{i},ACGT,{i % 2}\n" for i in range(23))
        if suffix == "csv"
        else "".join(f">r{i} sample {i}\nACGT\n" for i in range(23))
    )
    source.write_text(content)
    with read_inputs(source) as stream:
        original = list(stream)
    config = SplitConfig(
        method="kfold", n_splits=3, fold_id=2, shuffle=shuffle, preserve_order=preserve_order
    )
    expected = split(original, config=config, verbose=False)
    result = cast(Any, split)(
        source, config=config, output_path=tmp_path / "fold", workers=1, verbose=False
    )
    assert isinstance(result, FileSplitResult)
    assert result.to_dict() == expected.to_dict()
    for name in ("train", "valid"):
        assert tuple(result.get(name)) == tuple(expected.get(name))
        assert tuple(result.get(name)) == tuple(result.get(name))
    assert source.read_text() == content
    assert not list(tmp_path.glob(".dnakit-split-spool-*"))


def test_cli_accepts_fold_arguments_and_exports_the_selected_pair(tmp_path: Path) -> None:
    from dnakit.cli.app import app

    source = tmp_path / "train.fasta"
    source.write_text("".join(f">r{i}\nACGT\n" for i in range(25)))
    destination = tmp_path / "fold3"
    result = CliRunner().invoke(
        app,
        [
            "split",
            str(source),
            str(destination),
            "--method",
            "kfold",
            "--n-splits",
            "5",
            "--fold-id",
            "3",
            "--no-shuffle",
        ],
    )
    assert result.exit_code == 0, result.output
    manifest = json.loads((destination / "assignments.json").read_text())
    assert (manifest["seed"], manifest["n_splits"], manifest["fold_id"]) == (42, 5, 3)
    assert manifest["subsets"]["valid"] == [f"r{i}" for i in range(15, 20)]
    assert (destination / "train.fasta").is_file()
    assert (destination / "valid.fasta").is_file()
    assert not (destination / "test.fasta").exists()


def test_agent_discovers_and_executes_fold_configuration() -> None:
    from dnakit.tools import default_tool_registry

    registry = default_tool_registry(verbose=False)
    manifest = registry.describe("dnakit.datasets.split")
    assert "fold_id" in str(manifest["input_schema"])
    result = registry.execute(
        "dnakit.datasets.split",
        {
            "records": [{"sequence": "ACGT", "id": f"r{i}"} for i in range(20)],
            "config": {"method": "kfold", "n_splits": 5, "fold_id": 2},
        },
        allow_side_effects=True,
    )
    assert isinstance(result, dict)
    assert result["counts"] == {"train": 16, "valid": 4}
    assert (result["seed"], result["n_splits"], result["fold_id"]) == (42, 5, 2)


def test_workflow_passes_fold_arguments_and_global_seed(tmp_path: Path) -> None:
    from dnakit.workflows import run_workflow

    source = tmp_path / "input.fasta"
    source.write_text("".join(f">r{i}\nACGT\n" for i in range(20)))
    config_path = tmp_path / "workflow.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "dnakit-workflow-v1",
                "run_id": "kfold-test",
                "input": {"path": "input.fasta", "format": "fasta"},
                "output_dir": "workflow-output",
                "seed": 42,
                "steps": [
                    {
                        "id": "fold",
                        "operation": "split",
                        "input": "input",
                        "params": {"method": "kfold", "n_splits": 5, "fold_id": 2},
                    },
                    {
                        "id": "written",
                        "operation": "write",
                        "input": "fold",
                        "params": {"target": "valid.fasta", "subset": "valid"},
                    },
                ],
            }
        )
    )
    result = run_workflow(config_path, verbose=False)
    assert result.status == "succeeded"
    with read_inputs(tmp_path / "workflow-output" / "valid.fasta") as records:
        actual = tuple(record.id for record in records)
    expected = split(
        _records(20), config=SplitConfig(method="kfold", n_splits=5, fold_id=2), verbose=False
    )
    assert actual == expected.get("valid").ids
