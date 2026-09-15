"""Public verbosity, nested calls, and lazy/concurrent operation lifecycles."""

from __future__ import annotations

import asyncio
import importlib
import inspect
import io
from pathlib import Path
from typing import Any, cast

import pytest

from dnakit._verbosity import with_verbose
from dnakit.batch import BatchConfig, run_batch
from dnakit.core import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.io import RecordSource, read
from dnakit.ops import reverse_complement


@pytest.mark.parametrize(
    "module_name",
    [
        "dnakit",
        "dnakit.alignment",
        "dnakit.annotation",
        "dnakit.batch",
        "dnakit.chunking",
        "dnakit.comparative",
        "dnakit.config",
        "dnakit.core",
        "dnakit.datasets",
        "dnakit.descriptors",
        "dnakit.download",
        "dnakit.databases",
        "dnakit.encode",
        "dnakit.evaluation",
        "dnakit.fingerprints",
        "dnakit.io",
        "dnakit.llm",
        "dnakit.molbio",
        "dnakit.ops",
        "dnakit.patterns",
        "dnakit.predictions",
        "dnakit.references",
        "dnakit.representations",
        "dnakit.search",
        "dnakit.secondary_structure",
        "dnakit.similarity",
        "dnakit.standardize",
        "dnakit.structure3d",
        "dnakit.thermodynamics",
        "dnakit.visualization",
        "dnakit.workflows",
        "dnakit.bioinformatics.encoding",
        "dnakit.bioinformatics.ngs",
        "dnakit.bioinformatics.cloning",
    ],
)
def test_public_function_signatures(module_name: str) -> None:
    module = importlib.import_module(module_name)
    names = getattr(module, "__all__", [name for name in vars(module) if not name.startswith("_")])
    checked = 0
    for name in names:
        function = getattr(module, name)
        if not inspect.isfunction(function) or not function.__module__.startswith("dnakit"):
            continue
        if not hasattr(module, "__all__") and function.__module__ != module_name:
            continue
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.default is True
        assert parameter.annotation is bool
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        native_parameter = inspect.signature(inspect.unwrap(function)).parameters["verbose"]
        assert native_parameter.default is True
        assert native_parameter.kind is inspect.Parameter.KEYWORD_ONLY
        checked += 1
    assert checked > 0


@pytest.mark.parametrize("value", ["AATCGC", DNASequence("AATCGC")])
def test_text_and_object_start_finish_and_quiet(value: Any, capsys: Any) -> None:
    loud = reverse_complement(value)
    messages = capsys.readouterr()
    assert messages.out == ""
    assert messages.err.count("Starting") == 1
    assert messages.err.count("complete") == 1
    quiet = reverse_complement(value, verbose=False)
    assert loud == quiet
    assert capsys.readouterr() == ("", "")


def test_split_counts_are_part_of_completion(capsys: Any) -> None:
    records = [DNARecord(DNASequence("ACGT"), f"r{index}") for index in range(20)]
    config = SplitConfig(ratios={"train": 0.8, "valid": 0.1, "test": 0.1}, seed=42)
    result = split(records, config=config)
    messages = capsys.readouterr()
    assert "train=16, valid=2, test=2" in messages.err
    assert messages.err.count("Starting") == 1
    silent = split(records, config=config, verbose=False)
    assert silent.assignments == result.assignments
    assert capsys.readouterr() == ("", "")


def test_file_progress_and_nested_io_are_quiet(tmp_path: Path, capsys: Any) -> None:
    source = tmp_path / "input.fna"
    source.write_text(">r1\nAATCGC\n")
    result = reverse_complement(source, verbose=False, show_progress=True)
    assert result.record_count == 1
    assert capsys.readouterr() == ("", "")
    reverse_complement(source, output_path=tmp_path / "loud.fna")
    messages = capsys.readouterr().err
    assert messages.count("Starting") == 1
    assert messages.count("complete") == 1
    assert "Starting read" not in messages
    assert "Starting write" not in messages


def test_read_finishes_on_exhaustion_and_retains_recordsource(capsys: Any) -> None:
    records = read(io.StringIO(">r1\nACGT\n>r2\nAAGC\n"), format="fasta")
    assert isinstance(records, RecordSource)
    initial = capsys.readouterr().err
    assert "Starting read" in initial
    assert "complete" not in initial
    next(records)
    assert "complete" not in capsys.readouterr().err
    assert len(records.collect()) == 1
    assert records.closed
    assert capsys.readouterr().err.count("Read complete") == 1


def test_partial_stream_close_is_not_reported_as_complete(capsys: Any) -> None:
    records = read(io.StringIO(">r1\nACGT\n>r2\nAAGC\n"), format="fasta")
    next(records)
    records.close()
    output = capsys.readouterr().err
    assert "Read closed" in output
    assert "Read complete" not in output


def test_quiet_stream_does_not_silence_unrelated_operations(capsys: Any) -> None:
    records = read(io.StringIO(">r1\nACGT\n>r2\nAAGC\n"), format="fasta", verbose=False)
    next(records)
    assert capsys.readouterr() == ("", "")
    reverse_complement("ACGT")
    assert "complete" in capsys.readouterr().err
    records.close()
    assert capsys.readouterr() == ("", "")


def test_exception_and_invalid_verbose(capsys: Any) -> None:
    with pytest.raises(InvalidAlphabetError):
        reverse_complement("ACGX", verbose=False)
    assert capsys.readouterr() == ("", "")
    with pytest.raises(InvalidAlphabetError):
        reverse_complement("ACGX")
    output = capsys.readouterr().err
    assert "failed" in output
    assert "complete" not in output
    with pytest.raises(ConfigurationError, match="verbose"):
        reverse_complement("ACGT", verbose=cast(Any, "false"))


def test_async_verbosity_is_isolated(capsys: Any) -> None:
    @with_verbose
    async def quiet_task() -> str:
        await asyncio.sleep(0)
        result = reverse_complement("AATCGC", verbose=True)
        assert isinstance(result, str)
        return result

    @with_verbose
    async def loud_task() -> str:
        await asyncio.sleep(0)
        result = reverse_complement("AATCGC")
        assert isinstance(result, str)
        return result

    async def execute() -> Any:
        return await asyncio.gather(quiet_task(verbose=False), loud_task())

    assert asyncio.run(execute()) == ["GCGATT", "GCGATT"]
    output = capsys.readouterr().err
    assert "quiet task" not in output
    assert "Starting loud task" in output
    assert "Loud task complete" in output
    assert "reverse complement" not in output


def test_threaded_batch_inherits_quiet_setting(capsys: Any) -> None:
    def operation(record: DNARecord, _context: Any) -> Any:
        return reverse_complement(record)

    records = [DNARecord(DNASequence("ACGT"), f"r{index}") for index in range(4)]
    run_batch(
        records,
        operation,
        name="test reverse complement",
        config=BatchConfig(execution_mode="thread", jobs=2),
        verbose=False,
    )
    assert capsys.readouterr() == ("", "")
