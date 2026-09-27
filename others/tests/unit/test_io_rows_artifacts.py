"""Streaming experiment artifacts retain values and do not leave partial files."""

import gzip
import hashlib
import json
from collections.abc import Iterator
from pathlib import Path

import pytest

from dnakit import read
from dnakit.io import (
    RowWriter,
    file_digest,
    iter_rows,
    load_arrays,
    load_json,
    save_arrays,
    save_json,
    write_rows,
)


def test_single_pass_multiple_formats_and_gzip(tmp_path: Path) -> None:
    paths = [
        tmp_path / f"records.{suffix}" for suffix in ("csv", "tsv.gz", "json", "jsonl", "fasta")
    ]
    consumed = []

    def records() -> Iterator[dict[str, object]]:
        for i, sequence in enumerate(["ACGT", "NRY"]):
            consumed.append(i)
            yield {"id": f"s{i}", "sequence": sequence, "label": i}

    assert write_rows(records(), paths) == 2
    assert consumed == [0, 1]
    assert list(iter_rows(paths[0])) == list(iter_rows(paths[1]))
    assert (
        list(iter_rows(paths[2]))
        == list(iter_rows(paths[3]))
        == [
            {"id": "s0", "sequence": "ACGT", "label": 0},
            {"id": "s1", "sequence": "NRY", "label": 1},
        ]
    )
    with read(paths[4], verbose=False) as sequences:
        assert [(r.id, r.sequence.to_string()) for r in sequences] == [
            ("s0", "ACGT"),
            ("s1", "NRY"),
        ]


@pytest.mark.parametrize("suffix", ["csv", "tsv", "json", "jsonl", "fasta", "csv.gz"])
def test_empty_files(tmp_path: Path, suffix: str) -> None:
    target = tmp_path / f"empty.{suffix}"
    assert write_rows([], target, columns=["id", "sequence"]) == 0
    if suffix == "fasta":
        assert target.read_text() == ""
    else:
        assert list(iter_rows(target)) == []


def test_bad_generator_preserves_all_destinations(tmp_path: Path) -> None:
    paths = [tmp_path / "a.csv", tmp_path / "a.json"]
    for path in paths:
        path.write_text("original")

    def records() -> Iterator[dict[str, object]]:
        yield {"id": "a"}
        raise RuntimeError("source failed")

    with pytest.raises(RuntimeError, match="source failed"):
        write_rows(records(), paths, overwrite=True)
    assert all(p.read_text() == "original" for p in paths)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.csv", "a.json"]


def test_writer_schema_and_overwrite(tmp_path: Path) -> None:
    path = tmp_path / "a.csv"
    with RowWriter(path, columns=["id", "score"]) as writer:
        writer.write({"score": None, "id": "a"})
    assert list(iter_rows(path)) == [{"id": "a", "score": ""}]
    with pytest.raises(FileExistsError):
        write_rows([{"id": "b"}], path)
    with pytest.raises(ValueError, match="keys"):
        write_rows([{"id": "b"}, {"id": "c", "unexpected": 1}], path, overwrite=True)
    assert list(iter_rows(path)) == [{"id": "a", "score": ""}]
    with pytest.raises(ValueError, match="columns"):
        write_rows([], tmp_path / "unknown.csv")
    with pytest.raises(ValueError, match="unique"):
        write_rows([{"id": "a"}], [path, path])


def test_parquet_batches_roundtrip(tmp_path: Path) -> None:
    pytest.importorskip("pyarrow")
    path = tmp_path / "rows.parquet"
    rows = [{"id": str(i), "label": i, "value": float(i)} for i in range(9)]
    assert write_rows(iter(rows), path, batch_size=2) == 9
    assert list(iter_rows(path, batch_size=3)) == rows
    empty = tmp_path / "empty.parquet"
    assert write_rows([], empty, columns=["id", "sequence"]) == 0
    assert list(iter_rows(empty)) == []


def test_json_results_and_digest(tmp_path: Path) -> None:
    class Result:
        def to_dict(self) -> dict[str, object]:
            return {"title": "序列", "folds": [{"mcc": 0.7}, {"mcc": None}]}

    path = tmp_path / "nested" / "result.json"
    save_json(Result(), path)
    assert load_json(path) == Result().to_dict()
    assert file_digest(path) == hashlib.sha256(path.read_bytes()).hexdigest()
    assert file_digest(path, "md5") == hashlib.md5(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="JSON compliant"):
        save_json({"value": float("nan")}, path, overwrite=True)
    assert load_json(path) == Result().to_dict()


def test_npz_without_pickle(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    path = tmp_path / "arrays.npz"
    save_arrays({"ids": ["a", "b"], "embeddings": np.eye(2)}, path)
    result = load_arrays(path)
    np.testing.assert_array_equal(result["embeddings"], np.eye(2))
    assert result["ids"].tolist() == ["a", "b"]
    with pytest.raises(ValueError, match="Object"):
        save_arrays({"objects": np.array([{}], dtype=object)}, path, overwrite=True)
    np.testing.assert_array_equal(load_arrays(path)["embeddings"], np.eye(2))


def test_jsonl_nested_and_gzip(tmp_path: Path) -> None:
    path = tmp_path / "predictions.jsonl.gz"
    rows = [{"id": "a", "probabilities": {"label": {"1": 0.75}}}]
    write_rows(rows, path)
    assert list(iter_rows(path)) == rows
    with gzip.open(path, "rt") as source:
        assert json.loads(source.readline()) == rows[0]


def test_save_matplotlib_formats(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    from dnakit.visualization import save_figure

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0.2, 0.7])
    for suffix in ("svg", "png", "pdf"):
        path = tmp_path / "figures" / f"curve.{suffix}"
        save_figure(fig, path, close=False)
        assert path.stat().st_size > 100
    plt.close(fig)
