from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.model_selection import train_test_split


def available_genomic_benchmark_tasks(raw_root: str | Path) -> list[str]:
    root = Path(raw_root)
    if not root.exists():
        raise FileNotFoundError(f"genomic_benchmarks 原始数据目录不存在: {root}")
    task_names: list[str] = []
    for task_dir in sorted(root.iterdir()):
        if not task_dir.is_dir():
            continue
        if (task_dir / "train").exists() and (task_dir / "test").exists():
            task_names.append(task_dir.name)
    return task_names


def _load_split(raw_root: Path, task_name: str, split_name: str) -> pd.DataFrame:
    split_dir = raw_root / task_name / split_name
    if not split_dir.exists():
        raise FileNotFoundError(f"缺少 genomic_benchmarks split 目录: {split_dir}")

    records: list[dict[str, Any]] = []
    for label_name, label_dir in enumerate(sorted(path for path in split_dir.iterdir() if path.is_dir())):
        for sequence_path in sorted(label_dir.glob("*.txt")):
            sequence = sequence_path.read_text(encoding="utf-8").strip().upper()
            if not sequence:
                continue
            records.append(
                {
                    "sequence": sequence,
                    "label": int(label_name),
                    "label_original": label_dir.name,
                    "name": sequence_path.stem,
                }
            )

    if not records:
        raise RuntimeError(f"{split_dir} 下未找到任何可用序列。")
    frame = pd.DataFrame.from_records(records)
    frame["raw_length"] = frame["sequence"].str.len().astype(int)
    return frame


def _split_train_valid(
    train_df: pd.DataFrame,
    *,
    valid_ratio: float,
    split_seed: int,
    stratified: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if valid_ratio <= 0:
        return train_df.reset_index(drop=True), train_df.iloc[0:0].copy()
    if not 0 < float(valid_ratio) < 1:
        raise ValueError(f"valid_ratio 必须位于 (0, 1) 内，收到 {valid_ratio!r}")

    stratify = None
    if stratified:
        label_counts = train_df["label"].value_counts()
        if (label_counts >= 2).all():
            stratify = train_df["label"]

    train_part, valid_part = train_test_split(
        train_df,
        test_size=float(valid_ratio),
        random_state=int(split_seed),
        shuffle=True,
        stratify=stratify,
    )
    return train_part.reset_index(drop=True), valid_part.reset_index(drop=True)


def _subset_split(df: pd.DataFrame, max_samples: int | None) -> pd.DataFrame:
    if max_samples is None or max_samples <= 0:
        return df
    return df.iloc[: int(max_samples)].reset_index(drop=True)


def _canonicalize_split(
    df: pd.DataFrame,
    *,
    task_name: str,
    split_name: str,
    num_labels: int,
) -> pd.DataFrame:
    canonical = df.copy().reset_index(drop=True)
    canonical.insert(0, "sample_id", [f"{task_name}-{split_name}-{idx:08d}" for idx in range(len(canonical))])
    canonical["split"] = split_name
    canonical["task_name"] = task_name
    canonical["task_type"] = "classification"
    canonical["input_schema"] = "single_sequence"
    canonical["num_labels"] = int(num_labels)
    keep_columns = [
        "sample_id",
        "task_name",
        "task_type",
        "input_schema",
        "split",
        "sequence",
        "raw_length",
        "label",
        "label_original",
        "num_labels",
        "name",
    ]
    return canonical.loc[:, keep_columns]


def _length_summary(df: pd.DataFrame) -> dict[str, float | int]:
    if df.empty:
        return {"count": 0, "min": 0, "max": 0, "mean": 0.0, "median": 0.0}
    lengths = df["raw_length"].astype(int)
    return {
        "count": int(len(lengths)),
        "min": int(lengths.min()),
        "max": int(lengths.max()),
        "mean": float(lengths.mean()),
        "median": float(lengths.median()),
    }


def prepare_genomic_benchmark_task(
    *,
    task_name: str,
    raw_root: str | Path,
    output_root: str | Path,
    valid_ratio: float = 0.1,
    split_seed: int = 42,
    stratified: bool = True,
    overwrite: bool = False,
    max_samples_per_split: int | None = None,
) -> dict[str, Any]:
    raw_root_path = Path(raw_root)
    output_root_path = Path(output_root)
    task_output_dir = output_root_path / task_name
    meta_path = task_output_dir / "meta.json"
    if meta_path.exists() and not overwrite:
        return json.loads(meta_path.read_text(encoding="utf-8"))

    train_df = _load_split(raw_root_path, task_name, "train")
    test_df = _load_split(raw_root_path, task_name, "test")
    train_df, valid_df = _split_train_valid(
        train_df,
        valid_ratio=valid_ratio,
        split_seed=split_seed,
        stratified=stratified,
    )
    train_df = _subset_split(train_df, max_samples_per_split)
    valid_df = _subset_split(valid_df, max_samples_per_split)
    test_df = _subset_split(test_df, max_samples_per_split)

    label_mapping = {
        label_name: int(label_id)
        for label_id, label_name in enumerate(sorted(pd.concat([train_df["label_original"], test_df["label_original"]]).drop_duplicates()))
    }
    num_labels = len(label_mapping)
    train_canonical = _canonicalize_split(train_df, task_name=task_name, split_name="train", num_labels=num_labels)
    valid_canonical = _canonicalize_split(valid_df, task_name=task_name, split_name="valid", num_labels=num_labels)
    test_canonical = _canonicalize_split(test_df, task_name=task_name, split_name="test", num_labels=num_labels)

    meta = {
        "task_name": task_name,
        "task_family": "genomic_benchmarks",
        "task_type": "classification",
        "input_schema": "single_sequence",
        "num_labels": int(num_labels),
        "label_mapping": label_mapping,
        "inverse_label_mapping": {str(v): k for k, v in label_mapping.items()},
        "raw_root": str(raw_root_path),
        "source_task_dir": str(raw_root_path / task_name),
        "valid_ratio": float(valid_ratio),
        "split_seed": int(split_seed),
        "stratified_split": bool(stratified),
        "max_samples_per_split": None if max_samples_per_split is None else int(max_samples_per_split),
        "splits": {
            "train": _length_summary(train_canonical),
            "valid": _length_summary(valid_canonical),
            "test": _length_summary(test_canonical),
        },
    }

    task_output_dir.mkdir(parents=True, exist_ok=True)
    train_canonical.to_parquet(task_output_dir / "train.parquet", index=False)
    valid_canonical.to_parquet(task_output_dir / "valid.parquet", index=False)
    test_canonical.to_parquet(task_output_dir / "test.parquet", index=False)
    train_canonical.to_csv(task_output_dir / "train.csv", index=False)
    valid_canonical.to_csv(task_output_dir / "valid.csv", index=False)
    test_canonical.to_csv(task_output_dir / "test.csv", index=False)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return meta
