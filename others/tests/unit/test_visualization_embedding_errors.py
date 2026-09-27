"""Prediction-error plots from an embedding NPZ and a CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from dnakit.exceptions import InputFormatError
from dnakit.visualization import plot_embedding_errors


def _inputs(folder: Path) -> tuple[Path, Path]:
    np = pytest.importorskip("numpy")
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")
    embedding_file = folder / "embeddings.npz"
    np.savez(
        embedding_file,
        embeddings=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [10.0, 10.0, 0.0],
                [0.2, 0.8, 0.0],
                [0.9, 0.2, 0.0],
            ]
        ),
        ids=np.array(["a", "b", "c", "d", "t1", "t2", "t3"]),
        splits=np.array(["train"] * 4 + ["test"] * 3),
        labels=np.array([0, 1, 0, 1, 0, 1, 0]),
    )
    errors_csv = folder / "predictions.csv"
    errors_csv.write_text(
        "id,method,error,label,prediction,probability\n"
        "t3,full,FP,0,1,0.8\n"
        "t1,full,correct,0,0,0.2\n"
        "t2,full,FN,1,0,0.4\n"
        "t2,linear_probe,correct,1,1,0.7\n"
        "t3,linear_probe,correct,0,0,0.1\n"
        "t1,linear_probe,FP,0,1,0.6\n",
        encoding="utf-8",
    )
    return embedding_file, errors_csv


def test_embedding_error_plot_aligns_ids_and_fits_train_only(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    embedding_file, errors_csv = _inputs(tmp_path)
    paths = plot_embedding_errors(embedding_file, errors_csv, n_clusters=2, save_details=True)
    assert set(paths) == {"full", "linear_probe"}
    assert all(
        path.read_text(encoding="utf-8").lstrip().startswith("<?xml") for path in paths.values()
    )
    folder = tmp_path / "predictions"
    with np.load(folder / "projection.npz") as projection:
        assert projection["test_ids"].tolist() == ["t1", "t2", "t3"]
        assert projection["scaler_mean"].tolist() == pytest.approx([0.5, 0.5, 0.0])
        assert projection["coordinates"].shape == (3, 2)
    with (folder / "full_points.csv").open(newline="", encoding="utf-8") as stream:
        points = list(csv.DictReader(stream))
    assert [(row["id"], row["error"]) for row in points] == [
        ("t1", "correct"),
        ("t2", "FN"),
        ("t3", "FP"),
    ]
    with (folder / "cluster_errors.csv").open(newline="", encoding="utf-8") as stream:
        summary = list(csv.DictReader(stream))
    assert len(summary) == 4
    assert sum(int(row["n"]) for row in summary if row["method"] == "full") == 3
    assert sum(int(row["fp"]) for row in summary if row["method"] == "full") == 1
    assert sum(int(row["fn"]) for row in summary if row["method"] == "full") == 1
    with pytest.raises(FileExistsError):
        plot_embedding_errors(embedding_file, errors_csv)


def test_embedding_error_plot_accepts_flag_column_without_splits(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    pytest.importorskip("sklearn")
    pytest.importorskip("matplotlib")
    embedding_file = tmp_path / "values.npz"
    np.savez(
        embedding_file,
        embeddings=np.array([[0.0, 1.0], [1.0, 0.0], [0.3, 0.7]]),
        ids=np.array(["a", "b", "c"]),
    )
    errors_csv = tmp_path / "errors.csv"
    errors_csv.write_text("id,error\nc,1\na,0\nb,0\n", encoding="utf-8")
    paths = plot_embedding_errors(embedding_file, errors_csv)
    assert paths == {"errors": tmp_path / "errors" / "errors.svg"}
    assert paths["errors"].is_file()


@pytest.mark.parametrize(
    "replacement",
    [
        "t3,full,FP,0,1,0.8\n",  # one method omits an embedding ID
        "t3,full,correct,0,1,0.8\n",  # error contradicts the prediction
        "t3,full,FP,0,1,not-a-number\n",  # invalid optional probability
    ],
)
def test_embedding_error_plot_rejects_bad_csv_before_saving(
    tmp_path: Path, replacement: str
) -> None:
    embedding_file, errors_csv = _inputs(tmp_path)
    content = errors_csv.read_text(encoding="utf-8")
    if replacement == "t3,full,FP,0,1,0.8\n":
        content = content.replace("t2,full,FN,1,0,0.4\n", "")
    else:
        content = content.replace("t3,full,FP,0,1,0.8\n", replacement)
    errors_csv.write_text(content, encoding="utf-8")
    with pytest.raises(InputFormatError):
        plot_embedding_errors(embedding_file, errors_csv)
    assert not (tmp_path / "predictions").exists()
