from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from torch.utils.data.distributed import DistributedSampler


def resolve_processed_task_names(processed_root: str | Path) -> list[str]:
    root = Path(processed_root)
    if not root.exists():
        raise FileNotFoundError(f"处理后的 benchmark 数据目录不存在: {root}")
    task_names = [
        task_dir.name
        for task_dir in sorted(root.iterdir())
        if task_dir.is_dir()
        and ((task_dir / "meta.json").exists() or (task_dir / "raw" / "meta.json").exists())
    ]
    return task_names


def load_processed_task_meta(processed_root: str | Path, task_name: str) -> dict[str, Any]:
    task_dir = Path(processed_root) / task_name
    meta_path = task_dir / "meta.json"
    if not meta_path.exists():
        meta_path = task_dir / "raw" / "meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"缺少 meta.json: {meta_path}")
    return json.loads(meta_path.read_text(encoding="utf-8"))


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"不支持的数据文件格式: {path}")


def _candidate_raw_table_paths(
    task_dir: Path,
    task_name: str,
    meta: dict[str, Any] | None,
) -> list[Path]:
    candidates: list[Path] = []
    if meta is not None:
        raw_files = meta.get("raw_files")
        if isinstance(raw_files, dict):
            for key in ("parquet", "csv"):
                value = raw_files.get(key)
                if value:
                    candidates.append(Path(value))

    raw_dir = task_dir / "raw"
    candidates.extend(
        [
            raw_dir / f"{task_name}.parquet",
            raw_dir / f"{task_name}.csv",
            raw_dir / f"{task_dir.name}.parquet",
            raw_dir / f"{task_dir.name}.csv",
        ]
    )
    return candidates


def _normalize_materialized_raw_dataframe(
    dataframe: pd.DataFrame,
    *,
    split_name: str,
    task_name: str,
) -> pd.DataFrame:
    if "split" not in dataframe.columns:
        raise ValueError("materialized raw 单表缺少 split 列，无法拆分 train/valid/test。")

    split_frame = dataframe[dataframe["split"].astype(str) == str(split_name)].copy()
    if split_frame.empty:
        raise FileNotFoundError(f"{task_name} 的 materialized raw 表中没有 split={split_name!r}。")

    rename_map = {
        "seq": "sequence",
        "sequence_length": "raw_length",
    }
    split_frame = split_frame.rename(
        columns={
            old_name: new_name
            for old_name, new_name in rename_map.items()
            if old_name in split_frame.columns and new_name not in split_frame.columns
        }
    )
    if "sample_id" not in split_frame.columns and "name" in split_frame.columns:
        split_frame["sample_id"] = split_frame["name"].astype(str)
    if "task_name" not in split_frame.columns:
        split_frame["task_name"] = str(task_name)
    return split_frame.reset_index(drop=True)


def load_processed_task_dataframe(
    processed_root: str | Path,
    task_name: str,
    split_name: str,
) -> pd.DataFrame:
    task_dir = Path(processed_root) / task_name
    split_path = task_dir / f"{split_name}.parquet"
    csv_path = split_path.with_suffix(".csv")
    if split_path.exists():
        try:
            return pd.read_parquet(split_path)
        except Exception:
            if csv_path.exists():
                return pd.read_csv(csv_path)
            raise
    if csv_path.exists():
        return pd.read_csv(csv_path)

    meta = None
    try:
        meta = load_processed_task_meta(processed_root, task_name)
    except FileNotFoundError:
        meta = None

    read_errors: list[str] = []
    for raw_table_path in _candidate_raw_table_paths(task_dir, task_name, meta):
        if raw_table_path.exists():
            try:
                raw_frame = _read_table(raw_table_path)
            except Exception as exc:
                read_errors.append(f"{raw_table_path}: {type(exc).__name__}: {exc}")
                continue
            return _normalize_materialized_raw_dataframe(
                raw_frame,
                split_name=split_name,
                task_name=task_name,
            )
    if read_errors:
        raise RuntimeError(
            f"{task_name} 的 materialized raw 表均读取失败；"
            f"已尝试: {' | '.join(read_errors)}"
        )
    raise FileNotFoundError(f"缺少 split 文件: {split_path}")


class SequenceBenchmarkDataset(Dataset):
    def __init__(self, dataframe: pd.DataFrame) -> None:
        frame = dataframe.reset_index(drop=True).copy()
        self.is_pair = {"sequence_ref", "sequence_alt"}.issubset(frame.columns)
        if self.is_pair:
            required = {"sample_id", "sequence_ref", "sequence_alt", "label"}
        else:
            required = {"sample_id", "sequence", "label"}
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"统一 benchmark 数据缺少必要列: {missing}")
        self.dataframe = frame

    def __len__(self) -> int:
        return int(len(self.dataframe))

    def __getitem__(self, index: int) -> dict[str, Any]:
        row = self.dataframe.iloc[int(index)]
        item = {
            "sample_id": str(row["sample_id"]),
            "label": int(row["label"]),
            "task_name": str(row.get("task_name", "")),
            "split": str(row.get("split", "")),
            "name": str(row.get("name", "")),
        }
        if self.is_pair:
            item["sequence_ref"] = str(row["sequence_ref"])
            item["sequence_alt"] = str(row["sequence_alt"])
            item["raw_length_ref"] = int(row.get("raw_length_ref", len(item["sequence_ref"])))
            item["raw_length_alt"] = int(row.get("raw_length_alt", len(item["sequence_alt"])))
        else:
            item["sequence"] = str(row["sequence"])
            item["raw_length"] = int(row.get("raw_length", len(item["sequence"])))
        return item


class SequenceBenchmarkCollator:
    def __call__(self, batch: list[dict[str, Any]]) -> dict[str, Any]:
        labels = torch.tensor([int(item["label"]) for item in batch], dtype=torch.long)
        sample_ids = [item["sample_id"] for item in batch]
        payload = {
            "labels": labels,
            "sample_ids": sample_ids,
            "names": [item["name"] for item in batch],
        }
        if "sequence_ref" in batch[0]:
            sequences_ref = [item["sequence_ref"] for item in batch]
            sequences_alt = [item["sequence_alt"] for item in batch]
            payload.update(
                {
                    "images": sequences_ref,
                    "sequences_ref": sequences_ref,
                    "sequences_alt": sequences_alt,
                    "raw_lengths_ref": torch.tensor(
                        [int(item["raw_length_ref"]) for item in batch],
                        dtype=torch.long,
                    ),
                    "raw_lengths_alt": torch.tensor(
                        [int(item["raw_length_alt"]) for item in batch],
                        dtype=torch.long,
                    ),
                }
            )
            return payload

        sequences = [item["sequence"] for item in batch]
        payload.update(
            {
                "images": sequences,
                "sequences": sequences,
                "raw_lengths": torch.tensor(
                    [int(item["raw_length"]) for item in batch],
                    dtype=torch.long,
                ),
            }
        )
        return payload


def build_sequence_dataloader(
    dataset: Dataset,
    *,
    batch_size: int,
    num_workers: int,
    shuffle: bool,
    is_distributed: bool,
    pin_memory: bool,
) -> DataLoader:
    sampler = None
    if is_distributed:
        sampler = DistributedSampler(
            dataset,
            shuffle=shuffle,
            drop_last=False,
        )

    return DataLoader(
        dataset,
        batch_size=int(batch_size),
        shuffle=bool(shuffle) and sampler is None,
        sampler=sampler,
        num_workers=int(num_workers),
        pin_memory=bool(pin_memory),
        persistent_workers=bool(num_workers > 0),
        collate_fn=SequenceBenchmarkCollator(),
        drop_last=False,
    )
