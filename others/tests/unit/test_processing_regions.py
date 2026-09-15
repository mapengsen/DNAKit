"""Region processing retains native records, reports, boundaries and atomic writes."""

from pathlib import Path

import pytest

from dnakit import datasets, ops, standardize
from dnakit._processing.inputs import input_config
from dnakit.io import WriteConfig, read_inputs


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize(
    "operation",
    ["reverse", "complement", "trim", "transcribe", "translate", "normalize", "validate"],
)
def test_regions_match_codec(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, workers: int, operation: str
) -> None:
    from dnakit._processing import regions

    monkeypatch.setattr(regions, "_MIN_BYTES", 0)
    monkeypatch.setattr(regions, "_SHARD_BYTES", 100)
    source = tmp_path / "input.fna"
    source.write_bytes(
        b"".join(f">s{i} description\r\n".encode() + b"aCGTryN" * 30 + b"\r\n" for i in range(8))
    )
    suffix = ".jsonl" if operation == "validate" else ".fna"
    function = getattr(standardize if operation in {"normalize", "validate"} else ops, operation)
    kwargs = {"left": 2, "right": 3} if operation == "trim" else {}
    first, second = tmp_path / ("native" + suffix), tmp_path / ("regions" + suffix)
    function(
        source,
        output_path=first,
        read_config=input_config(source, None),
        workers=1,
        verbose=False,
        **kwargs,
    )
    function(source, output_path=second, workers=workers, verbose=False, **kwargs)
    assert first.read_bytes() == second.read_bytes()
    assert not list(tmp_path.glob(".*regions-*"))


def test_region_failure_keeps_previous_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from dnakit._processing import regions

    monkeypatch.setattr(regions, "_MIN_BYTES", 0)
    monkeypatch.setattr(regions, "_SHARD_BYTES", 100)
    source = tmp_path / "bad.fna"
    source.write_text(">a\n" + "ACGT" * 50 + "\n>b\n" + "ACGT" * 50 + "X\n")
    target = tmp_path / "out.fna"
    target.write_text("previous")
    with pytest.raises(Exception, match="symbols incompatible"):
        ops.reverse(
            source,
            output_path=target,
            write_config=WriteConfig(overwrite=True),
            workers=2,
            verbose=False,
        )
    assert target.read_text() == "previous"


@pytest.mark.parametrize("preserve_order", [True, False])
@pytest.mark.parametrize("workers", [1, 2])
def test_indexed_split_matches_native(
    tmp_path: Path, workers: int, preserve_order: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    from dnakit._processing import split_exports

    monkeypatch.setattr(split_exports, "_MIN_PARALLEL_BYTES", 0)
    paths = [tmp_path / "a.fna", tmp_path / "b.fna"]
    for prefix, path in enumerate(paths):
        path.write_text(
            "".join(f">s{prefix}_{i} description\n" + "acgTNY" * 30 + "\n" for i in range(13))
        )
    with read_inputs(paths, verbose=False) as stream:
        records = list(stream)
    config = datasets.SplitConfig(
        method="random",
        ratios={"train": 0.8, "valid": 0.1, "test": 0.1},
        seed=42,
        preserve_order=preserve_order,
    )
    expected = datasets.split(records, config=config, workers=1, verbose=False)
    result = datasets.split(
        paths, output_path=tmp_path / "split", config=config, workers=workers, verbose=False
    )
    assert expected.to_dict() == result.to_dict()
    for name in result.counts:
        assert tuple(result.get(name)) == tuple(expected.get(name))


def test_default_normalization_file_matches_audited_api(tmp_path: Path) -> None:
    import random

    rng = random.Random(42)
    # Include characters whose Unicode uppercase expands or becomes valid DNA.
    texts = ["".join(rng.choices("ACGTacgtryN Uu-.,ßſａｔ\u200b", k=200)) for _ in range(100)]
    source = tmp_path / "raw.txt"
    source.write_text("\n".join(texts) + "\n")
    target = tmp_path / "normalized.txt"
    standardize.normalize(source, output_path=target, workers=2, verbose=False)
    expected = [standardize.normalize(text, verbose=False).sequence.symbols for text in texts]
    assert target.read_text().splitlines() == expected


def test_constructor_retains_first_invalid_offset() -> None:
    from dnakit.core import DNASequence, Gap
    from dnakit.exceptions import InvalidAlphabetError

    with pytest.raises(InvalidAlphabetError) as error:
        DNASequence(["ACGT", Gap(2), "ACGT" * 1000 + "éX"], alphabet="iupac")
    assert error.value.context["part_index"] == 2
    assert error.value.context["part_offset"] == 4000
    assert error.value.context["symbol"] == "é"


def test_large_validation_report_keeps_all_positions(tmp_path: Path) -> None:
    import json

    source = tmp_path / "ambiguous.fna"
    source.write_text(">s\n" + "N" * 150_000 + "\n")
    target = tmp_path / "report.jsonl"
    standardize.validate(source, output_path=target, workers=1, verbose=False)
    report = json.loads(target.read_text())["result"]["ambiguity"]
    assert report["total_count"] == 150_000
    assert report["by_symbol"][0]["positions"] == list(range(150_000))
    assert len(report["occurrences"]) == 150_000
    assert report["occurrences"][-1]["symbol_offset"] == 149_999


@pytest.mark.parametrize("canonical", [False, True])
@pytest.mark.parametrize("method", ["kmer", "fingerprint"])
def test_cached_kmer_matrix_matches_pair_api(canonical: bool, method: str) -> None:
    from dnakit.core import DNARecord, DNASequence
    from dnakit.datasets._advanced_shared import pair_similarity, similarity_matrix

    texts = (
        ["ACGTACGT", "TGCANTGA", "NNNNNNNN", "ACGTACGT"]
        if method == "kmer"
        else ["ACGTACGT", "TGCATTGA", "AAAAAAAA", "ACGTACGT"]
    )
    records = tuple(
        DNARecord(DNASequence(text, alphabet="iupac"), f"s{i}") for i, text in enumerate(texts)
    )
    options = dict(method=method, k=3, canonical=canonical, max_alignment_cells=10000)
    result = similarity_matrix(records, **options)
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            assert result[left][right] == pair_similarity(records[left], records[right], **options)


@pytest.mark.parametrize("operation", ["reverse", "complement"])
def test_direction_regions_match_native(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    from dnakit.ops import _fasta_direction

    monkeypatch.setattr(_fasta_direction, "_MIN_BYTES", 0)
    monkeypatch.setattr(_fasta_direction, "_BLOCK_BYTES", 67)
    source = tmp_path / "input.fna"
    source.write_text(">first\n" + "acgtryN" * 100 + "\n>second\nACGTT\n")
    native, parallel = tmp_path / "native.fna", tmp_path / "parallel.fna"
    function = getattr(ops, operation)
    function(
        source, output_path=native, read_config=input_config(source, None), workers=1, verbose=False
    )
    function(source, output_path=parallel, workers=2, verbose=False)
    assert native.read_bytes() == parallel.read_bytes()


def test_file_rng_state_is_preserved(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import random

    from dnakit._processing import regions

    monkeypatch.setattr(regions, "_MIN_BYTES", 0)
    monkeypatch.setattr(regions, "_SHARD_BYTES", 100)
    source = tmp_path / "canonical.fna"
    source.write_text("".join(f">s{i}\n" + "ACGT" * 100 + "\n" for i in range(80)))
    first, second = random.Random(42), random.Random(42)
    ops.evolution_generate(
        source, rng=first, workers=1, verbose=False, output_path=tmp_path / "one.fna"
    )
    ops.evolution_generate(
        source, rng=second, workers=2, verbose=False, output_path=tmp_path / "two.fna"
    )
    assert (tmp_path / "one.fna").read_bytes() == (tmp_path / "two.fna").read_bytes()
    assert first.getstate() == second.getstate()


@pytest.mark.parametrize("encoding", ["utf-8", "utf-16"])
@pytest.mark.parametrize("body", ["acgt\nAC\n", "ACGT\nA\nCGT\n", "\u017fAC\n"])
def test_positioned_split_matches_codec_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, encoding: str, body: str
) -> None:
    from dnakit._processing import split_exports

    monkeypatch.setattr(split_exports, "_MIN_PARALLEL_BYTES", 0)
    source = tmp_path / "input.fna"
    source.write_text('>\u03b1 quoted "description"\n' + body + ">empty\n>third\nTGCA\n")
    config = datasets.SplitConfig(
        method="random", ratios={"train": 0.8, "valid": 0.1, "test": 0.1}, seed=42
    )
    writer = WriteConfig(encoding=encoding)
    for name, reader, workers in [
        ("codec", input_config(source, None), 1),
        ("positioned", None, 2),
    ]:
        datasets.split(
            source,
            config=config,
            read_config=reader,
            workers=workers,
            write_config=writer,
            output_path=tmp_path / name,
            verbose=False,
        )
    for path in (tmp_path / "codec").iterdir():
        assert path.read_bytes() == (tmp_path / "positioned" / path.name).read_bytes()


@pytest.mark.parametrize("operation", ["reverse", "split"])
def test_region_errors_keep_file_line_numbers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    from dnakit._processing import regions, split_exports
    from dnakit.exceptions import InputFormatError

    monkeypatch.setattr(regions, "_MIN_BYTES", 0)
    monkeypatch.setattr(regions, "_SHARD_BYTES", 100)
    monkeypatch.setattr(split_exports, "_MIN_PARALLEL_BYTES", 0)
    source = tmp_path / "bad.fna"
    source.write_text(">first\n" + "ACGT" * 100 + "\n>second\nACGTX\n")
    function = ops.reverse if operation == "reverse" else datasets.split
    kwargs = (
        {} if operation == "reverse" else {"config": datasets.SplitConfig(method="random", seed=42)}
    )
    contexts = []
    for name, reader, workers in [("codec", input_config(source, None), 1), ("regions", None, 2)]:
        with pytest.raises(InputFormatError) as error:
            function(
                source,
                output_path=tmp_path / (name + ".fna" if operation == "reverse" else name),
                read_config=reader,
                workers=workers,
                verbose=False,
                **kwargs,
            )
        contexts.append(dict(error.value.context))
    assert contexts[0] == contexts[1]
    assert contexts[1]["line_number"] == 3
