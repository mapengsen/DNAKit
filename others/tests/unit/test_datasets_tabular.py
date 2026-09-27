"""File-based grouped and K-fold exports agree with record-based splits."""

from __future__ import annotations

import hashlib
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, export_grouped_split, export_kfold_table, split
from dnakit.io import iter_rows, write_rows


def test_grouped_export_keeps_pairs_together_and_samples_deterministically(tmp_path: Path) -> None:
    source = tmp_path / "chunks.csv"
    rows = [
        {
            "id": f"pair-{index}:{strand}",
            "pair_id": f"pair-{index}",
            "strand": strand,
            "sequence": "NNNN" if index == 3 else "ACGT",
        }
        for index in range(40)
        for strand in ("+", "-")
    ]
    write_rows(rows, source, verbose=False)
    config = SplitConfig(seed=7)
    options = {
        "group_field": "pair_id",
        "config": config,
        "sample_sizes": {"train": 2, "valid": 2, "test": 2},
        "sample_filter": lambda row: row["strand"] == "+" and set(row["sequence"]) <= set("ACGT"),
        "sample_digest_field": "sequence",
    }
    summary = export_grouped_split(source, tmp_path / "first", **options, verbose=False)
    export_grouped_split(source, tmp_path / "second", **options, verbose=False)

    manifest = list(iter_rows(tmp_path / "first/assignment.csv", verbose=False))
    assignments = {row["pair_id"]: row["split"] for row in manifest}
    representatives = [DNARecord(DNASequence("A"), f"pair-{index}") for index in range(40)]
    expected = split(representatives, config=config, verbose=False)
    assert assignments == {item.record_id: item.split for item in expected.assignments}
    assert summary == {
        "seed": 7,
        "group_counts": {"train": 32, "valid": 4, "test": 4},
        "analysis_counts": {"train": 2, "valid": 2, "test": 2},
    }
    for name in ("train", "valid", "test"):
        target = tmp_path / "first" / f"{name}.csv"
        assert list(iter_rows(target, verbose=False)) == [
            row for row in rows if assignments[row["pair_id"]] == name
        ]
        sample = list(iter_rows(tmp_path / "first/analysis" / f"{name}.csv", verbose=False))
        assert len(sample) == 2
        assert all(row["strand"] == "+" and row["sequence"] == "ACGT" for row in sample)
        assert all(
            row["sequence_sha256"] == hashlib.sha256(row["sequence"].encode()).hexdigest()
            for row in sample
        )
        assert (tmp_path / "first/analysis" / f"{name}.csv").read_bytes() == (
            tmp_path / "second/analysis" / f"{name}.csv"
        ).read_bytes()


def test_kfold_table_matches_record_splits_and_counts(tmp_path: Path) -> None:
    pytest.importorskip("sklearn.model_selection")
    source = tmp_path / "train.csv"
    rows = [{"id": f"r{index}", "sequence": "ACGT", "label": str(index % 3)} for index in range(23)]
    write_rows(rows, source, verbose=False)
    config = SplitConfig(method="kfold", n_splits=5, seed=7)
    counts = export_kfold_table(
        source, tmp_path / "folds.csv", config=config, label_field="label", verbose=False
    )
    annotated = list(iter_rows(tmp_path / "folds.csv", verbose=False))
    assert [row["id"] for row in annotated] == [row["id"] for row in rows]

    records = [DNARecord(DNASequence(row["sequence"]), row["id"]) for row in rows]
    for fold in range(5):
        expected = split(records, config=replace(config, fold_id=fold), verbose=False)
        valid = {item.record_id for item in expected.assignments if item.split == "valid"}
        actual = {row["id"] for row in annotated if int(row["fold"]) == fold}
        assert actual == valid
        assert counts[fold] == {
            "fold": fold,
            "train": dict(Counter(row["label"] for row in rows if row["id"] not in valid)),
            "valid": dict(Counter(row["label"] for row in rows if row["id"] in valid)),
        }
