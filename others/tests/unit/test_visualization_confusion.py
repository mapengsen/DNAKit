"""Confusion matrices from one prediction CSV."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from dnakit.exceptions import InputFormatError
from dnakit.visualization import plot_confusion_matrix_csv


def test_confusion_matrix_counts_proportions_and_save(tmp_path: Path) -> None:
    source = tmp_path / "predictions.csv"
    source.write_text(
        "id,label,prediction\na,0,0\nb,0,0\nc,0,1\nd,1,0\ne,1,1\n",
        encoding="utf-8",
    )
    target = tmp_path / "figures" / "confusion.svg"
    artifact = plot_confusion_matrix_csv(
        source, labels=("0", "1"), title="Test predictions", output=target
    )

    assert artifact.kind == "confusion_matrix"
    assert artifact.metadata["labels"] == ("0", "1")
    assert artifact.metadata["counts"] == ((2, 1), (1, 1))
    assert artifact.metadata["row_totals"] == (3, 2)
    assert artifact.metadata["row_normalized"] == (
        (pytest.approx(2 / 3), pytest.approx(1 / 3)),
        (pytest.approx(0.5), pytest.approx(0.5)),
    )
    assert artifact.metadata["sample_count"] == 5
    assert target.read_text(encoding="utf-8") == artifact.svg

    root = ET.fromstring(artifact.svg)
    cells = root.findall(".//{http://www.w3.org/2000/svg}rect[@class='confusion-cell']")
    assert len(cells) == 8
    assert [cell.attrib["data-value"] for cell in cells[:4]] == ["2", "1", "1", "1"]
    assert {cell.attrib["data-panel"] for cell in cells} == {"counts", "row-normalized"}
    with pytest.raises(FileExistsError):
        plot_confusion_matrix_csv(source, output=target)
    assert (
        plot_confusion_matrix_csv(
            source, labels=("0", "1"), title="Test predictions", output=target, overwrite=True
        ).svg
        == artifact.svg
    )


def test_confusion_matrix_preserves_unobserved_class(tmp_path: Path) -> None:
    source = tmp_path / "one_class.csv"
    source.write_text("label,prediction\n0,0\n0,1\n", encoding="utf-8")

    artifact = plot_confusion_matrix_csv(source, labels=("0", "1"))

    assert artifact.metadata["counts"] == ((1, 1), (0, 0))
    assert artifact.metadata["row_normalized"] == ((0.5, 0.5), (0.0, 0.0))


def test_confusion_matrix_accepts_named_columns_and_multiple_classes(tmp_path: Path) -> None:
    source = tmp_path / "predictions.csv"
    source.write_text("truth,guess\nB,A\nA,C\nC,C\n", encoding="utf-8")

    artifact = plot_confusion_matrix_csv(
        source, true_column="truth", predicted_column="guess", labels=("A", "B", "C")
    )

    assert artifact.metadata["counts"] == ((0, 0, 1), (1, 0, 0), (0, 0, 1))
    assert artifact.metadata["row_normalized"] == (
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
    )


@pytest.mark.parametrize(
    "content",
    [
        "label,other\n0,1\n",
        "label,label,prediction\n0,1,1\n",
        "label,prediction\n",
        "label,prediction\n0,\n",
        "label,prediction\n0,1,extra\n",
    ],
)
def test_confusion_matrix_rejects_invalid_csv(tmp_path: Path, content: str) -> None:
    source = tmp_path / "invalid.csv"
    source.write_text(content, encoding="utf-8")

    with pytest.raises(InputFormatError):
        plot_confusion_matrix_csv(source, output=tmp_path / "plot.svg")
    assert not (tmp_path / "plot.svg").exists()


def test_confusion_matrix_rejects_class_outside_explicit_labels(tmp_path: Path) -> None:
    source = tmp_path / "predictions.csv"
    source.write_text("label,prediction\n0,2\n", encoding="utf-8")

    with pytest.raises(InputFormatError, match="observed class"):
        plot_confusion_matrix_csv(source, labels=("0", "1"))
