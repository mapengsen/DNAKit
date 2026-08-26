from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import torch
import pandas as pd

try:
    from tqdm.auto import tqdm
except ModuleNotFoundError:  # pragma: no cover
    def tqdm(iterable, **_: Any):
        return iterable


PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import (
    SequenceBenchmarkDataset,
    build_linear_probe_model,
    build_sequence_dataloader,
    get_model_spec,
    load_processed_task_dataframe,
    load_processed_task_meta,
    validate_finetune_method,
)
from visualdna.utils.seed_utils import set_seed


SPLIT_NAMES = ("train", "valid", "test")
SPLIT_TO_ID = {name: index for index, name in enumerate(SPLIT_NAMES)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="为 frozen linear probe 抽取 NT 对比模型的 backbone embedding cache。",
    )
    parser.add_argument("--task-name", required=True, help="处理后的 benchmark 任务名。")
    parser.add_argument("--processed-download", required=True, help="处理后的 NT benchmark 根目录。")
    parser.add_argument("--cache-path", required=True, help="输出 cache .pt 路径。")
    parser.add_argument("--model-name", required=True, help="对比模型名。")
    parser.add_argument("--raw-download", default=None, help="原始数据目录，仅写入元信息。")
    parser.add_argument("--checkpoint-override", default=None, help="覆盖模型默认 checkpoint。")
    parser.add_argument(
        "--token-readout",
        default="auto",
        choices=["auto", "mean", "last"],
        help="token 级读出方式；auto 使用模型默认策略。",
    )
    parser.add_argument(
        "--finetune-method",
        default="frozen_linear_probe",
        help="仅支持 frozen_linear_probe；保留参数是为了与原训练入口兼容。",
    )
    parser.add_argument(
        "--chunk-forward-batch-size",
        type=int,
        default=8,
        help="backbone 内部 chunk 微批大小。",
    )
    parser.add_argument("--batch-size", type=int, default=32, help="抽取 embedding 时的 DataLoader batch size。")
    parser.add_argument("--num-workers", type=int, default=4, help="DataLoader worker 数。")
    parser.add_argument("--device", default="cuda", help="抽取 embedding 使用的设备，例如 cuda / cpu。")
    parser.add_argument("--seed", type=int, default=42, help="随机种子。")
    parser.add_argument("--cache-dir", default=None, help="Hugging Face 缓存目录。")
    parser.add_argument(
        "--feature-dtype",
        default="float32",
        choices=["float32", "float16"],
        help="落盘 feature dtype；默认 float32 保持数值最稳。",
    )
    parser.add_argument(
        "--allow-remote-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="本地 checkpoint 加载失败时是否允许回退到 Hugging Face model id。",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="覆盖已有 cache。",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="如果 cache 已存在则直接退出。",
    )
    parser.add_argument(
        "--lock-timeout-seconds",
        type=int,
        default=86400,
        help="等待其他进程生成同一个 cache 的最长秒数。",
    )
    parser.add_argument(
        "--lock-poll-seconds",
        type=int,
        default=10,
        help="等待 cache lock 时的轮询间隔秒数。",
    )
    parser.add_argument("--max-train-samples", type=int, default=None, help="可选：仅抽取前 N 个 train 样本。")
    parser.add_argument("--max-valid-samples", type=int, default=None, help="可选：仅抽取前 N 个 valid 样本。")
    parser.add_argument("--max-test-samples", type=int, default=None, help="可选：仅抽取前 N 个 test 样本。")
    return parser


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def subset_dataframe(dataframe: pd.DataFrame, max_samples: int | None) -> pd.DataFrame:
    if max_samples is None or int(max_samples) <= 0:
        return dataframe
    return dataframe.iloc[: int(max_samples)].reset_index(drop=True)


def dataframe_fingerprint(dataframe: pd.DataFrame, *, pair_mode: bool) -> str:
    hasher = hashlib.sha256()
    columns = ["sample_id", "label"]
    if pair_mode:
        columns.extend(["sequence_ref", "sequence_alt"])
    else:
        columns.append("sequence")

    hasher.update(json.dumps(columns, ensure_ascii=False).encode("utf-8"))
    for row in dataframe[columns].itertuples(index=False, name=None):
        for value in row:
            hasher.update(str(value).encode("utf-8"))
            hasher.update(b"\0")
        hasher.update(b"\n")
    return hasher.hexdigest()


class CacheLock:
    def __init__(
        self,
        cache_path: Path,
        *,
        timeout_seconds: int,
        poll_seconds: int,
    ) -> None:
        self.cache_path = cache_path
        self.lock_path = cache_path.with_suffix(cache_path.suffix + ".lock")
        self.timeout_seconds = int(timeout_seconds)
        self.poll_seconds = max(int(poll_seconds), 1)
        self.acquired = False

    def __enter__(self):
        started_at = time.monotonic()
        last_report_at = 0.0
        self.lock_path.parent.mkdir(parents=True, exist_ok=True)
        while True:
            try:
                os.mkdir(self.lock_path)
            except FileExistsError:
                if self.cache_path.exists():
                    return self
                now = time.monotonic()
                if now - started_at > self.timeout_seconds:
                    raise TimeoutError(f"等待 embedding cache lock 超时: {self.lock_path}")
                if last_report_at == 0.0 or now - last_report_at >= max(self.poll_seconds, 30):
                    owner_text = ""
                    owner_path = self.lock_path / "owner.json"
                    try:
                        owner = json.loads(owner_path.read_text(encoding="utf-8"))
                        owner_text = f", owner_pid={owner.get('pid')}"
                    except Exception:
                        owner_text = ""
                    print(
                        f"[cache] waiting for lock: {self.lock_path} "
                        f"({int(now - started_at)}s elapsed{owner_text})",
                        flush=True,
                    )
                    last_report_at = now
                time.sleep(self.poll_seconds)
                continue
            self.acquired = True
            write_json(
                self.lock_path / "owner.json",
                {
                    "pid": os.getpid(),
                    "created_at": datetime.now(timezone.utc).isoformat(),
                    "cache_path": str(self.cache_path),
                },
            )
            return self

    def __exit__(self, exc_type, exc, tb):
        if not self.acquired:
            return False
        owner_path = self.lock_path / "owner.json"
        try:
            owner_path.unlink(missing_ok=True)
            self.lock_path.rmdir()
        except OSError:
            pass
        return False


def load_split_dataframes(args: argparse.Namespace) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    max_samples = {
        "train": args.max_train_samples,
        "valid": args.max_valid_samples,
        "test": args.max_test_samples,
    }
    dataframes: dict[str, pd.DataFrame] = {}
    fingerprints: dict[str, str] = {}
    meta = load_processed_task_meta(args.processed_download, args.task_name)
    pair_mode = meta.get("input_schema") == "pair_sequence"

    for split_name in SPLIT_NAMES:
        dataframe = load_processed_task_dataframe(args.processed_download, args.task_name, split_name)
        dataframe = subset_dataframe(dataframe, max_samples[split_name])
        dataframes[split_name] = dataframe
        fingerprints[split_name] = dataframe_fingerprint(dataframe, pair_mode=pair_mode)
    return dataframes, fingerprints


def encode_batch_features(model, batch: dict[str, Any]) -> torch.Tensor:
    token_readout = model.token_readout
    if "sequences_ref" in batch and "sequences_alt" in batch:
        features_ref = model.backbone.encode(batch["sequences_ref"], token_readout=token_readout)
        features_alt = model.backbone.encode(batch["sequences_alt"], token_readout=token_readout)
        return torch.cat([features_ref, features_alt, features_alt - features_ref], dim=-1)
    return model.backbone.encode(batch["sequences"], token_readout=token_readout)


def collect_split_features(
    *,
    model,
    dataframe: pd.DataFrame,
    split_name: str,
    args: argparse.Namespace,
    feature_dtype: torch.dtype,
) -> dict[str, Any]:
    dataset = SequenceBenchmarkDataset(dataframe)
    loader = build_sequence_dataloader(
        dataset,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        shuffle=False,
        is_distributed=False,
        pin_memory=torch.cuda.is_available() and str(args.device).startswith("cuda"),
    )

    feature_chunks: list[torch.Tensor] = []
    label_chunks: list[torch.Tensor] = []
    sample_ids: list[str] = []
    names: list[str] = []
    row_indices: list[int] = []

    progress = tqdm(loader, desc=f"extract {args.model_name}/{args.task_name}/{split_name}", unit="batch")
    cursor = 0
    print(
        f"[cache] extract split {split_name}: {len(dataset)} samples, "
        f"batch_size={args.batch_size}",
        flush=True,
    )
    for batch in progress:
        with torch.inference_mode():
            features = encode_batch_features(model, batch)
        labels = batch["labels"].view(-1).long().cpu()
        feature_chunks.append(features.detach().cpu().to(dtype=feature_dtype))
        label_chunks.append(labels)
        sample_ids.extend(str(value) for value in batch["sample_ids"])
        names.extend(str(value) for value in batch.get("names", [""] * labels.numel()))
        row_indices.extend(range(cursor, cursor + int(labels.numel())))
        cursor += int(labels.numel())

    feature_dim = int(feature_chunks[0].shape[1]) if feature_chunks else 0
    print(f"[cache] finished split {split_name}: {cursor} samples, feature_dim={feature_dim}", flush=True)
    return {
        "features": torch.cat(feature_chunks, dim=0) if feature_chunks else torch.empty((0, 0), dtype=feature_dtype),
        "labels": torch.cat(label_chunks, dim=0) if label_chunks else torch.empty((0,), dtype=torch.long),
        "sample_ids": sample_ids,
        "names": names,
        "row_indices": row_indices,
    }


def build_cache_payload(args: argparse.Namespace) -> dict[str, Any]:
    set_seed(args.seed)
    resolved_finetune_method = validate_finetune_method(get_model_spec(args.model_name), args.finetune_method)
    if resolved_finetune_method != "frozen_linear_probe":
        raise ValueError("embedding cache 仅支持 frozen_linear_probe。")

    print(f"[cache] load data: {args.model_name}/{args.task_name}", flush=True)
    task_meta = load_processed_task_meta(args.processed_download, args.task_name)
    if task_meta.get("task_type") != "classification":
        raise ValueError(f"当前 cache 脚本仅支持分类任务，收到 task_type={task_meta.get('task_type')!r}")
    pair_mode = task_meta.get("input_schema") == "pair_sequence"
    dataframes, fingerprints = load_split_dataframes(args)
    print(
        "[cache] split sizes: "
        + ", ".join(f"{name}={len(dataframes[name])}" for name in SPLIT_NAMES),
        flush=True,
    )
    spec = get_model_spec(args.model_name)
    print(f"[cache] build model: {args.model_name}", flush=True)
    model = build_linear_probe_model(
        model_name=args.model_name,
        num_labels=int(task_meta["num_labels"]),
        token_readout=args.token_readout,
        finetune_method=resolved_finetune_method,
        cache_dir=args.cache_dir,
        allow_remote_fallback=args.allow_remote_fallback,
        chunk_forward_batch_size=args.chunk_forward_batch_size,
        checkpoint_override=args.checkpoint_override,
        pair_mode=pair_mode,
    )
    model = model.to(args.device)
    model.eval()
    print(f"[cache] model ready on {args.device}", flush=True)

    feature_dtype = torch.float16 if args.feature_dtype == "float16" else torch.float32
    split_payloads = {
        split_name: collect_split_features(
            model=model,
            dataframe=dataframes[split_name],
            split_name=split_name,
            args=args,
            feature_dtype=feature_dtype,
        )
        for split_name in SPLIT_NAMES
    }

    features = torch.cat([split_payloads[name]["features"] for name in SPLIT_NAMES], dim=0)
    labels = torch.cat([split_payloads[name]["labels"] for name in SPLIT_NAMES], dim=0)
    split_ids = torch.cat(
        [
            torch.full(
                (int(split_payloads[name]["labels"].numel()),),
                SPLIT_TO_ID[name],
                dtype=torch.long,
            )
            for name in SPLIT_NAMES
        ],
        dim=0,
    )
    sample_ids: list[str] = []
    names: list[str] = []
    row_indices: list[int] = []
    split_names: list[str] = []
    for split_name in SPLIT_NAMES:
        count = int(split_payloads[split_name]["labels"].numel())
        sample_ids.extend(split_payloads[split_name]["sample_ids"])
        names.extend(split_payloads[split_name]["names"])
        row_indices.extend(split_payloads[split_name]["row_indices"])
        split_names.extend([split_name] * count)

    cache_meta = {
        "cache_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "task_name": args.task_name,
        "model_name": args.model_name,
        "model_display_name": spec.display_name,
        "model_spec": {
            "loader_kind": spec.loader_kind,
            "chunk_length": spec.chunk_length,
            "chunk_unit": spec.chunk_unit,
            "default_token_readout": spec.default_token_readout,
            "local_checkpoint": spec.local_checkpoint,
            "hf_model_id": spec.hf_model_id,
            "allowed_finetune_methods": list(spec.allowed_finetune_methods),
        },
        "processed_download": args.processed_download,
        "raw_download": args.raw_download,
        "checkpoint_override": args.checkpoint_override,
        "token_readout": model.token_readout,
        "requested_token_readout": args.token_readout,
        "finetune_method": resolved_finetune_method,
        "chunk_forward_batch_size": int(args.chunk_forward_batch_size),
        "feature_dtype": str(features.dtype).replace("torch.", ""),
        "feature_dim": int(features.shape[1]) if features.ndim == 2 else 0,
        "num_labels": int(task_meta["num_labels"]),
        "input_schema": task_meta.get("input_schema", "single_sequence"),
        "split_counts": {name: int(split_payloads[name]["labels"].numel()) for name in SPLIT_NAMES},
        "split_fingerprints": fingerprints,
        "seed": int(args.seed),
    }
    return {
        "features": features.contiguous(),
        "labels": labels.contiguous(),
        "split_ids": split_ids.contiguous(),
        "split_names": split_names,
        "row_indices": torch.tensor(row_indices, dtype=torch.long),
        "sample_ids": sample_ids,
        "names": names,
        "meta": cache_meta,
    }


def save_cache(cache_path: Path, payload: dict[str, Any]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
    torch.save(payload, tmp_path)
    tmp_path.replace(cache_path)
    write_json(cache_path.with_suffix(".json"), payload["meta"])


def main() -> None:
    args = build_parser().parse_args()
    cache_path = Path(args.cache_path)
    if cache_path.exists() and args.skip_existing and not args.overwrite:
        print(f"[cache] exists, skip: {cache_path}")
        return
    if cache_path.exists() and not args.overwrite:
        print(f"[cache] exists, reuse: {cache_path}")
        return

    with CacheLock(
        cache_path,
        timeout_seconds=args.lock_timeout_seconds,
        poll_seconds=args.lock_poll_seconds,
    ):
        if cache_path.exists() and not args.overwrite:
            print(f"[cache] created by another process, reuse: {cache_path}")
            return
        print(f"[cache] creating: {cache_path}", flush=True)
        payload = build_cache_payload(args)
        save_cache(cache_path, payload)
        print(
            "[cache] wrote "
            f"{cache_path} with {int(payload['labels'].numel())} samples, "
            f"feature_dim={payload['meta']['feature_dim']}"
        )


if __name__ == "__main__":
    main()
