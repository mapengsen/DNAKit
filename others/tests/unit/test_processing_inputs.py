"""Cross-format public processing contracts and metadata-preserving outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, get_type_hints

import pytest

from dnakit import datasets, ops, standardize
from dnakit.chunking import ChunkingConfig, chunk_sequences
from dnakit.core import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io import ReadConfig, WriteConfig, read_inputs, read_set, write


def test_processing_wrappers_preserve_resolvable_type_annotations() -> None:
    hints = get_type_hints(datasets.split)
    assert hints["config"] is datasets.SplitConfig
    assert hints["return"] is datasets.SplitResult
    assert hints["verbose"] is bool


def fasta(path: Path, *sequences: str) -> Path:
    path.write_text("".join(f">seq{i}\n{seq}\n" for i, seq in enumerate(sequences)))
    return path


def symbols(path: Path) -> list[str]:
    return [record.sequence.to_string() for record in read_set(path)]


@pytest.mark.parametrize(
    "name,kwargs,expected",
    [
        ("reverse", {}, "CGTAAT"),
        ("complement", {}, "ATTACG"),
        ("insert", {"position": 1, "fragment": "GG"}, "TGG AATGC".replace(" ", "")),
        ("delete", {"start": 1, "end": 3}, "TTGC"),
        ("substitute", {"start": 1, "end": 3, "fragment": "CC"}, "TCCTGC"),
        ("trim", {"left": 1, "right": 2}, "AAT"),
        ("subsequence", {"start": 1, "end": 4}, "AAT"),
        ("mask", {"intervals": [(1, 3)]}, "TNNTGC"),
        ("mutate", {"position": 0, "replacement": "A"}, "AAATGC"),
        ("rotate", {"offset": 2, "topology": "circular"}, "ATGCTA"),
        ("canonical_origin", {"topology": "circular"}, "AATGCT"),
    ],
)
def test_sequence_families_share_file_io(
    tmp_path: Path, name: str, kwargs: Any, expected: str
) -> None:
    source = fasta(tmp_path / "input.fna", "TAATGC")
    target = tmp_path / "output.fna"
    result = getattr(ops, name)(source, output_path=target, show_progress=False, **kwargs)
    assert result.record_count == 1
    assert symbols(target) == [expected]


@pytest.mark.parametrize("suffix", ["txt", "fna", "csv"])
def test_normalize_raw_file_before_dna_validation(tmp_path: Path, suffix: str) -> None:
    source = tmp_path / f"input.{suffix}"
    source.write_text(
        {"txt": "a cXgt\n", "fna": ">s\na cXgt\n", "csv": "id,sequence,label\ns,a cXgt,positive\n"}[
            suffix
        ]
    )
    target = tmp_path / f"output.{suffix}"
    result = standardize.normalize(source, output_path=target, show_progress=False)
    assert result.record_count == 1
    assert symbols(target) == ["ACGT"]
    if suffix == "csv":
        with target.open() as handle:
            assert next(csv.DictReader(handle))["label"] == "positive"


def test_csv_column_names_order_and_labels_survive_operations(tmp_path: Path) -> None:
    source = tmp_path / "input.csv"
    source.write_text(
        'label,id,dna,chromosome\n"positive, strong",s1,AAGC,chr1\nnegative,s2,TTAA,chr2\n'
    )
    target = tmp_path / "output.csv"
    ops.reverse_complement(source, target, sequence_column="dna", show_progress=False)
    with target.open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        assert reader.fieldnames == ["label", "id", "dna", "chromosome"]
    assert rows[0] == {"label": "positive, strong", "id": "s1", "dna": "GCTT", "chromosome": "chr1"}


@pytest.mark.parametrize("operation,expected", [("transcribe", "AUGGCC"), ("translate", "MA")])
@pytest.mark.parametrize("suffix", ["txt", "fna", "csv"])
def test_rna_and_protein_output_do_not_use_dna_alphabet(
    tmp_path: Path, operation: str, expected: str, suffix: str
) -> None:
    source = tmp_path / f"input.{suffix}"
    source.write_text(
        {"txt": "ATGGCC\n", "fna": ">s\nATGGCC\n", "csv": "id,sequence,label\ns,ATGGCC,1\n"}[suffix]
    )
    target = tmp_path / f"output.{suffix}"
    getattr(ops, operation)(source, output_path=target, show_progress=False)
    if suffix == "csv":
        with target.open() as handle:
            row = next(csv.DictReader(handle))
        assert row["sequence"] == expected and row["label"] == "1"
    else:
        assert target.read_text().strip().splitlines()[-1] == expected


def test_collection_split_and_dedup_preserve_labels(tmp_path: Path) -> None:
    source = tmp_path / "input.csv"
    source.write_text("id,sequence,label\n" + "".join(f"s{i},AAGC,{i % 2}\n" for i in range(10)))
    output = tmp_path / "split"
    result = datasets.split(
        source,
        config=datasets.SplitConfig(ratios={"train": 0.8, "valid": 0.1, "test": 0.1}, seed=42),
        output_path=output,
        show_progress=False,
    )
    assert dict(result.counts) == {"train": 8, "valid": 1, "test": 1}
    assert (output / "result.json").exists()
    ids = set()
    for name, count in result.counts.items():
        with (output / f"{name}.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == count
        for row in rows:
            assert row["label"] == str(int(row["id"][1:]) % 2)
            ids.add(row["id"])
    assert len(ids) == 10
    dedup = datasets.deduplicate(source, show_progress=False)
    assert dedup.output_count == 1
    assert (tmp_path / "input_deduplicate" / "deduplicate.csv").exists()


def test_cluster_filter_and_representatives(tmp_path: Path) -> None:
    source = tmp_path / "input.csv"
    source.write_text("id,sequence,chromosome\na,AAGC,chr1\nb,AAGC,chr2\nc,TTTT,chr3\n")
    result = datasets.cluster_sequences(source, show_progress=False)
    assert len(result.labels) == 3
    assert (tmp_path / "input_cluster_sequences" / "result.json").exists()
    kept = datasets.exclude_chromosomes(source, ["chr2"], show_progress=False)
    assert kept.ids == ("a", "c")
    representatives = datasets.select_representatives(source, [0, 0, 1], show_progress=False)
    assert representatives.representative_ids == ("a", "c")


def test_multi_file_inputs_and_concatenation(tmp_path: Path) -> None:
    first = fasta(tmp_path / "first.fna", "AAGC")
    second = fasta(tmp_path / "second.fna", "TTAA")
    with read_inputs([first, second, "ACGT"]) as records:
        assert [r.sequence.symbols for r in records] == ["AAGC", "TTAA", "ACGT"]
    target = tmp_path / "concat.fna"
    ops.concat([first, second], output_path=target, show_progress=False)
    assert symbols(target) == ["AAGCTTAA"]


def test_crossover_pairs_records_and_rejects_count_mismatch(tmp_path: Path) -> None:
    first = fasta(tmp_path / "first.fna", "AAAA", "CCCC")
    second = fasta(tmp_path / "second.fna", "TTTT", "GGGG")
    target = tmp_path / "children.fna"
    ops.crossover(first, second, position=2, output_path=target, show_progress=False)
    assert symbols(target) == ["AATT", "CCGG"]
    second.write_text(">s\nTTTT\n")
    target.unlink()
    with pytest.raises(ValueError):
        ops.crossover(first, second, position=2, output_path=target, show_progress=False)
    assert not target.exists()


def test_chunks_and_validation_reports(tmp_path: Path) -> None:
    source = fasta(tmp_path / "input.fna", "AAGCTTAA")
    target = tmp_path / "chunks.fna"
    result = chunk_sequences(
        source, config=ChunkingConfig(length=4), output_path=target, show_progress=False
    )
    assert result.record_count == 2
    assert symbols(target) == ["AAGC", "TTAA"]
    report = tmp_path / "validation.jsonl"
    standardize.validate(source, output_path=report, show_progress=False)
    assert json.loads(report.read_text())["id"] == "seq0"


def test_native_objects_keep_native_results() -> None:
    assert isinstance(ops.reverse(DNASequence("AAGC")), DNASequence)
    assert standardize.normalize("a cXg.t").sequence.symbols == "ACGT"
    assert ops.translate("ATGGCC") == "MA"


def test_output_protection_and_explicit_limits(tmp_path: Path) -> None:
    source = fasta(tmp_path / "input.fna", "AAGC")
    with pytest.raises(ConfigurationError):
        ops.reverse(
            source,
            output_path=source,
            write_config=WriteConfig(overwrite=True),
            show_progress=False,
        )
    with pytest.raises(InputFormatError):
        ops.reverse(source, read_config=ReadConfig(max_sequence_symbols=3), show_progress=False)
    assert not (tmp_path / "input_reverse.fna").exists()


@pytest.mark.parametrize(
    "operation,expected", [("reverse_complement", "GGCCAT"), ("transcribe", "AUGGCC")]
)
def test_custom_csv_schema_keeps_aliases_and_extra_canonical_names(
    tmp_path: Path, operation: str, expected: str
) -> None:
    source = tmp_path / "input.csv"
    with source.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "dna", "desc", "meta", "sequence", "label"])
        writer.writerow(["s", "ATGGCC", "example", '{"batch": "b1"}', "extra", "1"])
    target = tmp_path / "result.csv"
    config = ReadConfig(
        id_column="name", sequence_column="dna", description_column="desc", metadata_column="meta"
    )
    getattr(ops, operation)(source, output_path=target, read_config=config, show_progress=False)
    with target.open() as handle:
        reader = csv.DictReader(handle)
        row = next(reader)
        assert reader.fieldnames == ["name", "dna", "desc", "meta", "sequence", "label"]
    assert row["dna"] == expected
    assert row["desc"] == "example" and row["sequence"] == "extra" and row["label"] == "1"
    assert json.loads(row["meta"])["batch"] == "b1"


def test_transcribe_dnakit_schema_and_chunk_annotation_mapping(tmp_path: Path) -> None:
    source = tmp_path / "input.csv"
    write([DNARecord(DNASequence("ATGGCC"), "s", metadata={"nested": {"x": [1]}})], source)
    ops.transcribe(source, show_progress=False)
    with (tmp_path / "input_transcribe.csv").open() as handle:
        reader = csv.DictReader(handle)
        row = next(reader)
        assert "parts" not in reader.fieldnames and "alphabet" not in reader.fieldnames
    assert row["sequence"] == "AUGGCC"
    assert json.loads(row["metadata"])["nested"] == {"x": [1]}

    annotated = tmp_path / "annotated.jsonl"
    write(
        [
            DNARecord(
                DNASequence("AACCGGTT"),
                "s",
                features=[DNAFeature("gene", Interval(1, 7))],
                letter_annotations={"quality": list(range(8))},
                metadata={"label": "positive"},
            )
        ],
        annotated,
    )
    chunk_sequences(annotated, config=ChunkingConfig(length=4), show_progress=False)
    records = read_set(tmp_path / "annotated_chunk_sequences.jsonl")
    assert [record.features[0].location for record in records] == [Interval(1, 4), Interval(0, 3)]
    assert tuple(records[1].letter_annotations["quality"]) == (4, 5, 6, 7)
    assert records[1].metadata["chunking"]["source_start"] == 4
    assert records[1].metadata["label"] == "positive"


def test_normalize_multiple_files_and_record_text(tmp_path: Path) -> None:
    first = fasta(tmp_path / "first.fna", "a cXgt")
    second = fasta(tmp_path / "second.fna", "tYta")
    target = tmp_path / "normalized.fna"
    standardize.normalize([first, second], output_path=target, show_progress=False)
    assert symbols(target) == ["ACGT", "TYTA"]
    assert ops.trim_record("AAGCTT", left=1, right=1).record.sequence.symbols == "AGCT"
