"""F2/F3 and sibling operations: alphabet, atomic output and dataset identity."""

from __future__ import annotations

import math
from typing import Any, cast

import pytest

from dnakit import datasets, ops
from dnakit.core import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split
from dnakit.exceptions import (
    ConfigurationError,
    CoordinateError,
    InputFormatError,
    InvalidAlphabetError,
)
from dnakit.io import ReadConfig, WriteConfig, read_set


@pytest.mark.parametrize(
    "source,expected",
    [("A", "T"), ("AAGC", "GCTT"), ("ARYN", "NRYT"), ("acgt", "ACGT"), ("aa gc", "GCTT"), ("", "")],
)
def test_reverse_complement_minimal_and_normalized_text(source: Any, expected: Any) -> None:
    assert ops.reverse_complement(source, verbose=False) == expected


@pytest.mark.parametrize("source", ["ACXT", "AC1T", "AC💡T"])
@pytest.mark.parametrize(
    "operation", [ops.reverse, ops.complement, ops.reverse_complement, ops.transcribe]
)
def test_sequence_operations_reject_invalid_text(source: Any, operation: Any) -> None:
    with pytest.raises(InvalidAlphabetError):
        operation(source, verbose=False)


@pytest.mark.parametrize("content,count", [("", 0), (">empty\n", 1)])
def test_empty_fasta_is_distinct_from_empty_record(tmp_path: Any, content: Any, count: Any) -> None:
    source, target = tmp_path / "input.fna", tmp_path / "output.fna"
    source.write_text(content)
    result = ops.reverse_complement(source, target, workers=1, verbose=False)
    assert result.record_count == count
    records = list(read_set(target))
    assert len(records) == count
    if records:
        assert records[0].id == "empty" and records[0].sequence.symbols == ""


OPERATIONS = [
    (ops.reverse, {}),
    (ops.complement, {}),
    (ops.reverse_complement, {}),
    (ops.transcribe, {}),
    (ops.translate, {}),
    (ops.subsequence, {"start": 0, "end": 2}),
    (ops.insert, {"position": 1, "fragment": "A"}),
    (ops.delete, {"start": 0, "end": 1}),
    (ops.mask, {"intervals": [(0, 1)]}),
]


@pytest.mark.parametrize("operation,kwargs", OPERATIONS, ids=[f.__name__ for f, _ in OPERATIONS])
@pytest.mark.parametrize("suffix", ["fna", "csv", "txt"])
@pytest.mark.parametrize("workers", [1, 2])
def test_shared_file_operations_reject_late_invalid_record_atomically(
    tmp_path: Any, operation: Any, kwargs: Any, suffix: Any, workers: Any
) -> None:
    source, target = tmp_path / f"input.{suffix}", tmp_path / f"output.{suffix}"
    source.write_text(
        {
            "fna": ">good\nACGTAC\n>bad\nAC1TAC\n",
            "csv": "id,sequence,label\ngood,ACGTAC,1\nbad,AC1TAC,0\n",
            "txt": "ACGTAC\nAC1TAC\n",
        }[suffix]
    )
    target.write_text("previous complete output")
    with pytest.raises(InputFormatError):
        operation(
            source,
            output_path=target,
            workers=workers,
            write_config=WriteConfig(overwrite=True),
            verbose=False,
            **kwargs,
        )
    assert target.read_text() == "previous complete output"
    assert set(tmp_path.iterdir()) == {source, target}


@pytest.mark.parametrize(
    "kwargs", [{"start": -1, "end": 2}, {"start": 3, "end": 2}, {"start": 0, "end": 5}]
)
def test_invalid_coordinates_are_rejected(kwargs: Any) -> None:
    with pytest.raises(CoordinateError):
        cast(Any, ops.subsequence)("AAGC", verbose=False, **kwargs)


@pytest.mark.parametrize(
    "size,counts",
    [
        (0, (0, 0, 0)),
        (1, (1, 0, 0)),
        (2, (2, 0, 0)),
        (9, (7, 1, 1)),
        (10, (8, 1, 1)),
        (11, (9, 1, 1)),
        (100, (80, 10, 10)),
    ],
)
@pytest.mark.parametrize("file_input", [False, True])
def test_small_splits_preserve_sequence_label_identity(
    tmp_path: Any, size: Any, counts: Any, file_input: Any
) -> None:
    records = [
        DNARecord(DNASequence("ACGT" + "A" * i), f"r{i}", metadata={"label": i})
        for i in range(size)
    ]
    if file_input:
        source = tmp_path / "input.csv"
        source.write_text(
            "id,sequence,label\n"
            + "".join(f"{r.id},{r.sequence.symbols},{i}\n" for i, r in enumerate(records))
        )
        result = cast(Any, split)(
            source,
            config=SplitConfig(seed=42),
            output_path=tmp_path / "split",
            workers=1,
            verbose=False,
        )
    else:
        result = cast(Any, split)(records, config=SplitConfig(seed=42), workers=1, verbose=False)
    assert tuple(result.counts[name] for name in ["train", "valid", "test"]) == counts
    observed: list[int] = []
    for subset in result.subsets:
        for record in subset.records:
            index = int(record.id[1:])
            assert record.sequence.symbols == records[index].sequence.symbols
            assert int(str(record.metadata["label"])) == index
            observed.append(index)
    assert sorted(observed) == list(range(size))
    assert sorted(a.input_index for a in result.assignments) == list(range(size))


@pytest.mark.parametrize(
    "ratios",
    [
        {"train": -0.1, "test": 1.1},
        {"train": 0, "test": 0},
        {"train": math.nan, "test": 0.2},
        {"train": math.inf, "test": 0.2},
        {"train": 8, "valid": 1, "test": 1},
        {"train": 0.8, "test": 0.3},
    ],
)
def test_invalid_split_ratios(ratios: Any) -> None:
    with pytest.raises(ConfigurationError):
        SplitConfig(ratios=ratios)


def test_duplicate_ids_remain_traceable_by_input_index() -> None:
    records = [
        DNARecord(DNASequence(sequence), "same", metadata={"label": i})
        for i, sequence in enumerate(["AAGC", "TTTT", "CCCC"])
    ]
    result = cast(Any, split)(records, config=SplitConfig(seed=42), verbose=False)
    assert sorted(a.input_index for a in result.assignments) == [0, 1, 2]
    assert sorted(
        (r.sequence.symbols, r.metadata["label"])
        for subset in result.subsets
        for r in subset.records
    ) == sorted((r.sequence.symbols, r.metadata["label"]) for r in records)


def test_group_split_keeps_duplicates_rc_and_overlaps_together() -> None:
    records = [
        DNARecord(DNASequence(sequence), f"r{group}_{index}", metadata={"source_group": str(group)})
        for group in range(10)
        for index, sequence in enumerate(["AAGCTTAC", "AAGCTTAC", "GTAAGCTT", "GCTTAC"])
    ]
    result = cast(Any, split)(
        records,
        config=SplitConfig(method="group", metadata_key="source_group", seed=42),
        verbose=False,
    )
    memberships: dict[Any, set[str]] = {}
    for subset in result.subsets:
        for record in subset.records:
            memberships.setdefault(record.metadata["source_group"], set()).add(subset.name)
    assert len(memberships) == 10
    assert all(len(names) == 1 for names in memberships.values())


@pytest.mark.parametrize(
    "operation,kwargs",
    [
        (datasets.split, {"config": SplitConfig(seed=42)}),
        (datasets.deduplicate, {}),
        (datasets.cluster_sequences, {}),
    ],
)
def test_collection_families_reject_invalid_files_without_replacing_output(
    tmp_path: Any, operation: Any, kwargs: Any
) -> None:
    source = tmp_path / "input.csv"
    source.write_text("id,sequence\ngood,ACGT\nbad,AC1T\n")
    target = tmp_path / "result"
    target.mkdir()
    sentinel = target / "existing.txt"
    sentinel.write_text("keep")
    with pytest.raises(InputFormatError):
        operation(
            source,
            output_path=target,
            write_config=WriteConfig(overwrite=True),
            workers=1,
            verbose=False,
            **kwargs,
        )
    assert sentinel.read_text() == "keep"
    assert list(target.iterdir()) == [sentinel]


def test_explicit_reader_limit_stays_effective(tmp_path: Any) -> None:
    source = tmp_path / "input.fna"
    source.write_text(">r\nACGTA\n")
    with pytest.raises(InputFormatError):
        ops.reverse_complement(
            source, read_config=ReadConfig(max_sequence_symbols=4), verbose=False
        )
    assert list(tmp_path.iterdir()) == [source]


@pytest.mark.parametrize(
    "suffix,content",
    [
        ("fastq", "@r\nACGT\n+\n!!!\n"),
        ("fna", "ACGT\n>r\nTT\n"),
        ("csv", "id,sequence\nr,ACGT,extra\n"),
        ("jsonl", '{"id":"r","sequence":"ACGT"\n'),
    ],
)
@pytest.mark.parametrize(
    "operation,kwargs",
    [(ops.reverse_complement, {}), (datasets.split, {"config": SplitConfig(seed=42)})],
)
def test_malformed_file_structure_does_not_create_output(
    tmp_path: Any, suffix: Any, content: Any, operation: Any, kwargs: Any
) -> None:
    source = tmp_path / f"input.{suffix}"
    source.write_text(content)
    target = tmp_path / ("split" if operation is datasets.split else "output.fna")
    with pytest.raises(InputFormatError):
        operation(source, output_path=target, workers=1, verbose=False, **kwargs)
    assert not target.exists()


def test_leakage_detector_finds_cross_split_duplicates_and_reverse_complements() -> None:
    from dnakit.core import DNASet
    from dnakit.datasets import LeakageConfig, detect_leakage

    report = detect_leakage(
        {
            "train": DNASet([DNARecord(DNASequence("AAGC"), "original")]),
            "test": DNASet(
                [DNARecord(DNASequence("AAGC"), "duplicate"), DNARecord(DNASequence("GCTT"), "rc")]
            ),
        },
        config=LeakageConfig(method="kmer", k=2, canonical=True, threshold=1),
        verbose=False,
    )
    assert report.has_leakage
    assert {event.right_id for event in report.events} == {"duplicate", "rc"}


@pytest.mark.parametrize("operation", ["read", "reverse_complement", "normalize", "split"])
def test_truncated_gzip_has_controlled_error_and_no_partial_output(
    tmp_path: Any, operation: Any
) -> None:
    import gzip

    from dnakit import standardize
    from dnakit.io import read

    source, target = tmp_path / "input.fna.gz", tmp_path / "output.fna"
    source.write_bytes(gzip.compress(b">first\nACGT\n>second\nTTTT\n")[:-5])
    with pytest.raises(InputFormatError, match="INPUT_READ_FAILED"):
        if operation == "read":
            with read(source, verbose=False) as records:
                list(records)
        elif operation == "split":
            cast(Any, split)(
                source,
                config=SplitConfig(seed=42),
                output_path=tmp_path / "split",
                workers=1,
                verbose=False,
            )
        else:
            function = standardize.normalize if operation == "normalize" else ops.reverse_complement
            cast(Any, function)(source, output_path=target, workers=1, verbose=False)
    assert list(tmp_path.iterdir()) == [source]
