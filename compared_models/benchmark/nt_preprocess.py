from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from sklearn.model_selection import train_test_split


def available_nt_tasks(raw_root: str | Path) -> list[str]:
    root = Path(raw_root)
    if not root.exists():
        raise FileNotFoundError(f"NT 原始数据目录不存在: {root}")

    task_names: list[str] = []
    for task_dir in sorted(root.iterdir()):
        if not task_dir.is_dir():
            continue
        if (task_dir / "train.parquet").exists() and (task_dir / "test.parquet").exists():
            task_names.append(task_dir.name)
    return task_names


def _load_nt_split(raw_root: Path, task_name: str, split: str) -> pd.DataFrame:
    split_path = raw_root / task_name / f"{split}.parquet"
    if not split_path.exists():
        raise FileNotFoundError(f"缺少 NT {task_name} 的 {split} 文件: {split_path}")
    df = pd.read_parquet(split_path)
    required = {"sequence", "label"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{split_path} 缺少必要列: {missing}")
    return df.copy()


def _normalize_sequences(df: pd.DataFrame) -> pd.DataFrame:
    normalized = df.copy()
    normalized["sequence"] = normalized["sequence"].astype(str).str.upper()
    normalized["raw_length"] = normalized["sequence"].str.len().astype(int)
    if "name" not in normalized.columns:
        normalized["name"] = ""
    normalized["name"] = normalized["name"].fillna("").astype(str)
    return normalized


def _build_label_mapping(train_df: pd.DataFrame, test_df: pd.DataFrame) -> dict[Any, int]:
    labels = pd.concat([train_df["label"], test_df["label"]], axis=0)
    unique_labels = sorted(labels.drop_duplicates().tolist())
    return {label: index for index, label in enumerate(unique_labels)}


def _apply_label_mapping(df: pd.DataFrame, label_mapping: dict[Any, int]) -> pd.DataFrame:
    mapped = df.copy()
    mapped["label_original"] = mapped["label"]
    mapped["label"] = mapped["label"].map(label_mapping).astype(int)
    return mapped


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
    return (
        train_part.reset_index(drop=True),
        valid_part.reset_index(drop=True),
    )


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
    canonical["num_labels"] = int(num_labels)
    keep_columns = [
        "sample_id",
        "task_name",
        "task_type",
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
        return {
            "count": 0,
            "min": 0,
            "max": 0,
            "mean": 0.0,
            "median": 0.0,
        }
    lengths = df["raw_length"].astype(int)
    return {
        "count": int(len(lengths)),
        "min": int(lengths.min()),
        "max": int(lengths.max()),
        "mean": float(lengths.mean()),
        "median": float(lengths.median()),
    }


def _build_meta(
    *,
    task_name: str,
    train_df: pd.DataFrame,
    valid_df: pd.DataFrame,
    test_df: pd.DataFrame,
    label_mapping: dict[Any, int],
    raw_root: Path,
    valid_ratio: float,
    split_seed: int,
    stratified: bool,
) -> dict[str, Any]:
    inverse_mapping = {str(index): str(label) for label, index in label_mapping.items()}
    return {
        "task_name": task_name,
        "task_family": "NT",
        "task_type": "classification",
        "num_labels": int(len(label_mapping)),
        "label_mapping": {str(label): int(index) for label, index in label_mapping.items()},
        "inverse_label_mapping": inverse_mapping,
        "raw_root": str(raw_root),
        "source_train": str(raw_root / task_name / "train.parquet"),
        "source_test": str(raw_root / task_name / "test.parquet"),
        "valid_ratio": float(valid_ratio),
        "split_seed": int(split_seed),
        "stratified_split": bool(stratified),
        "splits": {
            "train": _length_summary(train_df),
            "valid": _length_summary(valid_df),
            "test": _length_summary(test_df),
        },
    }


def prepare_nt_task(
    *,
    task_name: str,
    raw_root: str | Path,
    output_root: str | Path,
    valid_ratio: float = 0.1,
    split_seed: int = 42,
    stratified: bool = True,
    overwrite: bool = False,
) -> dict[str, Any]:
    raw_root_path = Path(raw_root)
    output_root_path = Path(output_root)
    task_output_dir = output_root_path / task_name

    meta_path = task_output_dir / "meta.json"
    if meta_path.exists() and not overwrite:
        return json.loads(meta_path.read_text(encoding="utf-8"))

    train_df = _normalize_sequences(_load_nt_split(raw_root_path, task_name, "train"))
    test_df = _normalize_sequences(_load_nt_split(raw_root_path, task_name, "test"))
    label_mapping = _build_label_mapping(train_df, test_df)
    train_df = _apply_label_mapping(train_df, label_mapping)
    test_df = _apply_label_mapping(test_df, label_mapping)

    train_df, valid_df = _split_train_valid(
        train_df,
        valid_ratio=valid_ratio,
        split_seed=split_seed,
        stratified=stratified,
    )

    num_labels = len(label_mapping)
    train_canonical = _canonicalize_split(
        train_df,
        task_name=task_name,
        split_name="train",
        num_labels=num_labels,
    )
    valid_canonical = _canonicalize_split(
        valid_df,
        task_name=task_name,
        split_name="valid",
        num_labels=num_labels,
    )
    test_canonical = _canonicalize_split(
        test_df,
        task_name=task_name,
        split_name="test",
        num_labels=num_labels,
    )

    meta = _build_meta(
        task_name=task_name,
        train_df=train_canonical,
        valid_df=valid_canonical,
        test_df=test_canonical,
        label_mapping=label_mapping,
        raw_root=raw_root_path,
        valid_ratio=valid_ratio,
        split_seed=split_seed,
        stratified=stratified,
    )

    task_output_dir.mkdir(parents=True, exist_ok=True)
    train_canonical.to_parquet(task_output_dir / "train.parquet", index=False)
    valid_canonical.to_parquet(task_output_dir / "valid.parquet", index=False)
    test_canonical.to_parquet(task_output_dir / "test.parquet", index=False)
    train_canonical.to_csv(task_output_dir / "train.csv", index=False)
    valid_canonical.to_csv(task_output_dir / "valid.csv", index=False)
    test_canonical.to_csv(task_output_dir / "test.csv", index=False)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return meta
