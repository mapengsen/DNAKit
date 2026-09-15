"""Actual spawned workers, deterministic outputs, and disk-backed split contracts."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from dnakit import datasets, ops, standardize
from dnakit._parallel import ordered_map, pack, resolve_workers, unpack
from dnakit.core import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.exceptions import ConfigurationError, CoordinateError
from dnakit.io import WriteConfig, read_inputs


def _pid_job(value: int) -> tuple[int, int]:
    return value, os.getpid()


def _failure_job(value: int) -> int:
    if value == 70:
        raise ConfigurationError("worker failed", code="TEST_WORKER_FAILURE", context={"value": 70})
    return value


def test_process_execution_preserves_order() -> None:
    results = list(ordered_map(_pid_job, range(256), workers=2))
    assert [value for value, _ in results] == list(range(256))
    if resolve_workers(2) > 1:
        assert all(pid != os.getpid() for _, pid in results)
    assert len({pid for _, pid in results}) <= 2


def test_process_exception_retains_error_details() -> None:
    with pytest.raises(ConfigurationError) as error:
        list(ordered_map(_failure_job, range(256), workers=2))
    assert error.value.code == "TEST_WORKER_FAILURE"
    assert error.value.context["value"] == 70


@pytest.mark.parametrize("workers", [0, -2, True, 1.5, "2"])
def test_invalid_workers_fail_before_creating_output(tmp_path: Path, workers: object) -> None:
    with pytest.raises(ConfigurationError, match="workers"):
        ops.reverse_complement(tmp_path / "missing.fna", workers=workers)
    with pytest.raises(ConfigurationError, match="workers"):
        ops.reverse("ACGT", workers=workers)


@pytest.mark.parametrize(
    "operation", ["reverse_complement", "trim", "transcribe", "translate", "evolution_generate"]
)
def test_sequence_processes_match_serial_files(tmp_path: Path, operation: str) -> None:
    source = tmp_path / "input.csv"
    source.write_text(
        "id,sequence,label\n" + "".join(f"s{i},ATGGCC,positive\n" for i in range(160))
    )
    kwargs = {"left": 1, "right": 1} if operation == "trim" else {}
    if operation == "evolution_generate":
        kwargs = {"seed": 42}
    function = getattr(ops, operation)
    first, second = tmp_path / "serial.csv", tmp_path / "parallel.csv"
    function(source, output_path=first, workers=1, verbose=False, **kwargs)
    function(source, output_path=second, workers=2, verbose=False, **kwargs)
    assert first.read_bytes() == second.read_bytes()


def test_normalization_handles_invalid_raw_text_in_processes(tmp_path: Path) -> None:
    source = tmp_path / "input.txt"
    source.write_text("a cXgt\n" * 160)
    for workers in (1, 2):
        standardize.normalize(
            source, output_path=tmp_path / f"{workers}.txt", workers=workers, verbose=False
        )
    assert (tmp_path / "1.txt").read_bytes() == (tmp_path / "2.txt").read_bytes()
    assert (tmp_path / "2.txt").read_text().splitlines() == ["ACGT"] * 160


def test_immutable_annotation_transport() -> None:
    record = DNARecord(
        DNASequence("ACGT"),
        "s",
        features=[DNAFeature("gene", Interval(0, 3))],
        metadata={"labels": {"x": [1]}},
        letter_annotations={"quality": [1, 2, 3, 4]},
    )
    assert unpack(pack(record)) == record


@pytest.mark.parametrize(
    "operation", ["deduplicate", "deduplicate_iupac", "cluster_sequences", "hierarchical_cluster"]
)
def test_dataset_processes_preserve_native_results(operation: str) -> None:
    records = [DNARecord(DNASequence("AAGCTT" if i % 2 else "TTGGCC"), f"s{i}") for i in range(18)]
    if operation == "deduplicate":
        records *= 4  # enough key tasks to exercise a process batch
    function = getattr(datasets, operation)
    first = function(records, workers=1, verbose=False)
    second = function(records, workers=2, verbose=False)
    assert first.to_dict() == second.to_dict()


@pytest.mark.parametrize("preserve_order", [True, False])
@pytest.mark.parametrize("shuffle", [True, False])
def test_file_split_matches_native_and_reopens_subsets(
    tmp_path: Path, preserve_order: bool, shuffle: bool
) -> None:
    source = tmp_path / "input.csv"
    source.write_text("id,dna,label\n" + "".join(f"s{i},AAGC,{i % 2}\n" for i in range(30)))
    with read_inputs(source, sequence_column="dna") as records:
        original = list(records)
    config = datasets.SplitConfig(
        ratios={"train": 0.8, "valid": 0.1, "test": 0.1},
        seed=42,
        preserve_order=preserve_order,
        shuffle=shuffle,
    )
    native = datasets.split(original, config=config, workers=1, verbose=False)
    for workers in (1, 2):
        result = datasets.split(
            source,
            config=config,
            sequence_column="dna",
            workers=workers,
            output_path=tmp_path / str(workers),
            verbose=False,
        )
        assert isinstance(result, datasets.FileSplitResult)
        assert result.to_dict() == native.to_dict()
        for name in result.counts:
            assert tuple(result.get(name)) == tuple(native.get(name))
            assert tuple(result.get(name)) == tuple(result.get(name))
    assert not list(tmp_path.glob(".dnakit-split-spool-*"))


def test_failed_parallel_file_does_not_install_output(tmp_path: Path) -> None:
    source = tmp_path / "input.fna"
    source.write_text("".join(f">s{i}\n{'A' if i == 70 else 'AAGC'}\n" for i in range(160)))
    output = tmp_path / "output.fna"
    with pytest.raises(CoordinateError):
        ops.subsequence(source, start=0, end=3, output_path=output, workers=2, verbose=False)
    assert not output.exists()
    assert not list(tmp_path.glob(".output.fna.*"))


def test_large_direction_blocks_match_serial_with_partial_tail() -> None:
    from dnakit._parallel import transform_text, worker_scope

    text = "ACGTRYMKBDHVN" * 660_000 + "AAT"
    for reverse, complement in ((True, True), (True, False), (False, True)):
        with worker_scope(1):
            expected = transform_text(text, reverse=reverse, complement=complement)
        with worker_scope(2):
            actual = transform_text(text, reverse=reverse, complement=complement)
        assert actual == expected


def test_file_split_compressed_output_reopens(tmp_path: Path) -> None:
    source = tmp_path / "input.fna"
    source.write_text("".join(f">s{i}\n{'AAGC' * 100}\n" for i in range(10)))
    result = datasets.split(
        source,
        config=datasets.SplitConfig(ratios={"train": 0.8, "test": 0.2}, seed=42),
        output_path=tmp_path / "split",
        write_config=WriteConfig(compression="gzip"),
        workers=2,
        verbose=False,
    )
    assert result.get("train").path.suffix == ".gz"
    assert len(tuple(result.get("train"))) == 8


def test_cli_workers_validation_and_scope(tmp_path: Path) -> None:
    from typer.testing import CliRunner

    from dnakit._parallel import _workers
    from dnakit.cli.app import app

    runner = CliRunner()
    assert runner.invoke(app, ["--workers", "0", "info"]).exit_code != 0
    assert _workers.get() == -1
    source = tmp_path / "input.fna"
    source.write_text(">s\nAAGC\n")
    output = tmp_path / "output.fna"
    result = runner.invoke(app, ["--workers", "2", "deduplicate", str(source), str(output)])
    assert result.exit_code == 0, result.output
    assert _workers.get() == -1
