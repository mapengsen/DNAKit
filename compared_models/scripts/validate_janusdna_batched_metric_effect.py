from __future__ import annotations

import argparse
import json
import math
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from tqdm.auto import tqdm

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import (  # noqa: E402
    SequenceBenchmarkDataset,
    build_linear_probe_model,
    build_sequence_dataloader,
    get_model_spec,
    load_processed_task_dataframe,
    load_processed_task_meta,
    resolve_processed_task_names,
    validate_finetune_method,
)
from compared_models.scripts.train_nt_cached_linear_probe import (  # noqa: E402
    CachedFeatureClassificationTrainer,
    CachedFeatureDataset,
)
from visualdna.metrics import classification_report  # noqa: E402
from visualdna.utils.seed_utils import set_seed  # noqa: E402
from visualdna.utils.training import build_optimizer_and_scheduler  # noqa: E402


SPLIT_NAMES = ("train", "valid", "test")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "JanusDNA 专用验证：比较 serial encode 与 batched encode 的 embedding 差异，"
            "并训练相同线性探针判断差异是否影响最终分类指标。"
        ),
    )
    parser.add_argument(
        "--task-name",
        required=True,
        help="NT task 名；支持 all 或逗号分隔，例如 promoter_all,enhancers。",
    )
    parser.add_argument("--processed-download", required=True, help="处理后的 NT benchmark 根目录。")
    parser.add_argument("--log-root", required=True, help="验证输出目录。")
    parser.add_argument("--raw-download", default=None, help="原始数据目录，仅写入元信息。")
    parser.add_argument("--model-name", default="janusdna", choices=["janusdna"], help="固定为 janusdna。")
    parser.add_argument("--checkpoint-override", default=None, help="覆盖默认 JanusDNA checkpoint。")
    parser.add_argument("--cache-dir", default=None, help="Hugging Face 缓存目录。")
    parser.add_argument(
        "--allow-remote-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="本地 checkpoint 加载失败时是否允许回退到 Hugging Face。",
    )
    parser.add_argument(
        "--token-readout",
        default="auto",
        choices=["auto", "mean", "last"],
        help="token 级读出方式。",
    )
    parser.add_argument("--device", default="cuda", help="运行设备，例如 cuda、cuda:2 或 cpu。")
    parser.add_argument("--num-workers", type=int, default=0, help="DataLoader worker 数。")
    parser.add_argument("--extract-batch-size", type=int, default=32, help="抽取 embedding 的 DataLoader batch size。")
    parser.add_argument(
        "--chunk-forward-batch-size",
        type=int,
        default=32,
        help="batched encode 时一次 forward 的 chunk 数上限。",
    )
    parser.add_argument(
        "--max-train-samples",
        type=int,
        default=4096,
        help="每个 task 抽取的 train 样本数；<=0 表示全量。",
    )
    parser.add_argument(
        "--max-valid-samples",
        type=int,
        default=1024,
        help="每个 task 抽取的 valid 样本数；<=0 表示全量。",
    )
    parser.add_argument(
        "--max-test-samples",
        type=int,
        default=1024,
        help="每个 task 抽取的 test 样本数；<=0 表示全量。",
    )
    parser.add_argument(
        "--sample-strategy",
        default="stratified",
        choices=["head", "stratified"],
        help="抽样方式；stratified 会按 label 分层抽样。",
    )
    parser.add_argument("--epochs", type=int, default=20, help="线性探针最大训练 epoch。")
    parser.add_argument("--batch-size", type=int, default=64, help="线性探针训练 batch size。")
    parser.add_argument("--lr", type=float, default=1e-3, help="线性探针学习率。")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="线性探针权重衰减。")
    parser.add_argument("--optimizer-name", default="adamw", choices=["adam", "adamw"], help="优化器。")
    parser.add_argument("--warmup-ratio", type=float, default=0.1, help="warmup 比例。")
    parser.add_argument("--min-lr-ratio", type=float, default=0.1, help="cosine annealing 最终 lr 比例。")
    parser.add_argument("--early-stopping-patience", type=int, default=5, help="早停 patience。")
    parser.add_argument(
        "--monitor",
        default="accuracy",
        choices=["accuracy", "acc", "f1", "mcc", "auroc", "loss"],
        help="早停监控指标。",
    )
    parser.add_argument("--seed", type=int, default=42, help="抽样和训练随机种子。")
    parser.add_argument(
        "--feature-train-dtype",
        default="float32",
        choices=["float32", "float16"],
        help="线性探针训练时的 feature dtype。",
    )
    parser.add_argument(
        "--metric-delta-threshold",
        type=float,
        default=0.005,
        help="判断指标差异是否可忽略的绝对阈值。",
    )
    parser.add_argument(
        "--prediction-flip-threshold",
        type=float,
        default=0.01,
        help="判断 test 预测翻转率是否可忽略的阈值。",
    )
    parser.add_argument(
        "--tensorboard",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="验证线性探针是否写 TensorBoard。",
    )
    parser.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否显示进度条。",
    )
    return parser


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_task_names(task_name: str, processed_root: str | Path) -> list[str]:
    if task_name == "all":
        task_names = resolve_processed_task_names(processed_root)
    else:
        task_names = [item.strip() for item in task_name.split(",") if item.strip()]
    if not task_names:
        raise ValueError("没有解析出任何 task。")
    return task_names


def _sample_group_counts(group_sizes: dict[Any, int], total: int) -> dict[Any, int]:
    available = sum(group_sizes.values())
    target = min(int(total), int(available))
    if target <= 0:
        return {key: 0 for key in group_sizes}

    raw_counts = {
        key: target * (count / available)
        for key, count in group_sizes.items()
    }
    counts = {
        key: min(group_sizes[key], int(math.floor(value)))
        for key, value in raw_counts.items()
    }
    remaining = target - sum(counts.values())
    fractional_order = sorted(
        group_sizes,
        key=lambda key: (raw_counts[key] - math.floor(raw_counts[key]), group_sizes[key]),
        reverse=True,
    )
    while remaining > 0:
        progressed = False
        for key in fractional_order:
            if counts[key] >= group_sizes[key]:
                continue
            counts[key] += 1
            remaining -= 1
            progressed = True
            if remaining <= 0:
                break
        if not progressed:
            break
    return counts


def sample_dataframe(
    dataframe: pd.DataFrame,
    *,
    max_samples: int | None,
    seed: int,
    strategy: str,
) -> pd.DataFrame:
    if max_samples is None or int(max_samples) <= 0 or len(dataframe) <= int(max_samples):
        return dataframe.reset_index(drop=True)
    sample_count = int(max_samples)
    if strategy == "head" or "label" not in dataframe.columns:
        return dataframe.iloc[:sample_count].reset_index(drop=True)

    pieces: list[pd.DataFrame] = []
    group_sizes = {
        label: len(group)
        for label, group in dataframe.groupby("label", sort=True)
    }
    counts = _sample_group_counts(group_sizes, sample_count)
    for label, group in dataframe.groupby("label", sort=True):
        count = counts.get(label, 0)
        if count <= 0:
            continue
        pieces.append(group.sample(n=count, random_state=int(seed)))
    sampled = pd.concat(pieces, axis=0).sort_index()
    return sampled.reset_index(drop=True)


def load_sampled_dataframes(args: argparse.Namespace, task_name: str) -> dict[str, pd.DataFrame]:
    max_samples = {
        "train": args.max_train_samples,
        "valid": args.max_valid_samples,
        "test": args.max_test_samples,
    }
    return {
        split_name: sample_dataframe(
            load_processed_task_dataframe(args.processed_download, task_name, split_name),
            max_samples=max_samples[split_name],
            seed=int(args.seed) + SPLIT_NAMES.index(split_name),
            strategy=args.sample_strategy,
        )
        for split_name in SPLIT_NAMES
    }


def build_extraction_loader(
    dataframe: pd.DataFrame,
    *,
    batch_size: int,
    num_workers: int,
    device: str,
) -> DataLoader:
    return build_sequence_dataloader(
        SequenceBenchmarkDataset(dataframe),
        batch_size=int(batch_size),
        num_workers=int(num_workers),
        shuffle=False,
        is_distributed=False,
        pin_memory=torch.cuda.is_available() and str(device).startswith("cuda"),
    )


def encode_batch(model, batch: dict[str, Any], *, mode: str, pair_mode: bool) -> torch.Tensor:
    encode_fn = model.backbone.encode_serial if mode == "serial" else model.backbone.encode_batched
    token_readout = model.token_readout
    if pair_mode:
        ref_features = encode_fn(batch["sequences_ref"], token_readout=token_readout)
        alt_features = encode_fn(batch["sequences_alt"], token_readout=token_readout)
        return torch.cat([ref_features, alt_features, alt_features - ref_features], dim=-1)
    return encode_fn(batch["sequences"], token_readout=token_readout)


def collect_features(
    *,
    model,
    dataframes: dict[str, pd.DataFrame],
    mode: str,
    task_name: str,
    pair_mode: bool,
    args: argparse.Namespace,
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    for split_name in SPLIT_NAMES:
        loader = build_extraction_loader(
            dataframes[split_name],
            batch_size=args.extract_batch_size,
            num_workers=args.num_workers,
            device=args.device,
        )
        feature_chunks: list[torch.Tensor] = []
        label_chunks: list[torch.Tensor] = []
        sample_ids: list[str] = []
        iterator = tqdm(
            loader,
            desc=f"{task_name}/{mode}/{split_name}",
            unit="batch",
            disable=not args.progress,
        )
        for batch in iterator:
            with torch.inference_mode():
                features = encode_batch(model, batch, mode=mode, pair_mode=pair_mode)
            feature_chunks.append(features.detach().cpu().float())
            label_chunks.append(batch["labels"].view(-1).detach().cpu().long())
            sample_ids.extend(str(value) for value in batch["sample_ids"])
        results[split_name] = {
            "features": torch.cat(feature_chunks, dim=0),
            "labels": torch.cat(label_chunks, dim=0),
            "sample_ids": sample_ids,
        }
    return results


def compare_features(
    serial: dict[str, dict[str, Any]],
    batched: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for split_name in SPLIT_NAMES:
        left = serial[split_name]["features"].float()
        right = batched[split_name]["features"].float()
        diff = (left - right).abs()
        per_sample_max = diff.max(dim=1).values if diff.numel() else torch.empty(0)
        cosine = torch.nn.functional.cosine_similarity(left, right, dim=1) if left.numel() else torch.empty(0)
        output[split_name] = {
            "shape": list(left.shape),
            "max_abs_diff": float(diff.max().item()) if diff.numel() else 0.0,
            "mean_abs_diff": float(diff.mean().item()) if diff.numel() else 0.0,
            "p95_abs_diff": float(torch.quantile(diff.flatten(), 0.95).item()) if diff.numel() else 0.0,
            "p99_abs_diff": float(torch.quantile(diff.flatten(), 0.99).item()) if diff.numel() else 0.0,
            "per_sample_max_abs_mean": float(per_sample_max.mean().item()) if per_sample_max.numel() else 0.0,
            "per_sample_max_abs_p95": (
                float(torch.quantile(per_sample_max, 0.95).item())
                if per_sample_max.numel()
                else 0.0
            ),
            "cosine_similarity_mean": float(cosine.mean().item()) if cosine.numel() else 0.0,
            "cosine_similarity_min": float(cosine.min().item()) if cosine.numel() else 0.0,
            "fraction_abs_diff_gt_1e_3": float((diff > 1e-3).float().mean().item()) if diff.numel() else 0.0,
            "fraction_abs_diff_gt_1e_4": float((diff > 1e-4).float().mean().item()) if diff.numel() else 0.0,
        }
    return output


def build_feature_loader(
    features: torch.Tensor,
    labels: torch.Tensor,
    *,
    batch_size: int,
    num_workers: int,
    device: str,
    shuffle: bool,
    feature_dtype: torch.dtype,
) -> DataLoader:
    dataset = CachedFeatureDataset(features.to(dtype=feature_dtype), labels.long())
    return DataLoader(
        dataset,
        batch_size=int(batch_size),
        shuffle=bool(shuffle),
        num_workers=int(num_workers),
        pin_memory=torch.cuda.is_available() and str(device).startswith("cuda"),
        persistent_workers=bool(num_workers > 0),
        drop_last=False,
    )


def build_feature_loaders(
    features_by_split: dict[str, dict[str, Any]],
    *,
    args: argparse.Namespace,
) -> dict[str, DataLoader]:
    feature_dtype = torch.float16 if args.feature_train_dtype == "float16" else torch.float32
    return {
        split_name: build_feature_loader(
            features_by_split[split_name]["features"],
            features_by_split[split_name]["labels"],
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            device=args.device,
            shuffle=(split_name == "train"),
            feature_dtype=feature_dtype,
        )
        for split_name in SPLIT_NAMES
    }


def load_best_state(model: nn.Module, checkpoint_path: Path) -> None:
    try:
        payload = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    except TypeError:
        payload = torch.load(checkpoint_path, map_location="cpu")
    state = payload.get("model_state_dict") if isinstance(payload, dict) else None
    if not isinstance(state, dict):
        raise ValueError(f"checkpoint 缺少 model_state_dict: {checkpoint_path}")
    model.load_state_dict(state)


def collect_predictions(
    model: nn.Module,
    loader: DataLoader,
    *,
    device: str,
) -> dict[str, Any]:
    model.eval()
    labels: list[torch.Tensor] = []
    predictions: list[torch.Tensor] = []
    probabilities: list[torch.Tensor] = []
    model_dtype = next(model.parameters()).dtype
    with torch.inference_mode():
        for batch in loader:
            features = batch["features"].to(device, dtype=model_dtype, non_blocking=True)
            batch_labels = batch["labels"].view(-1).long().to(device, non_blocking=True)
            logits = model(features)
            probs = torch.softmax(logits, dim=-1)
            labels.append(batch_labels.cpu())
            predictions.append(logits.argmax(dim=-1).cpu())
            probabilities.append(probs.cpu())
    label_array = torch.cat(labels, dim=0).numpy()
    prediction_array = torch.cat(predictions, dim=0).numpy()
    probability_array = torch.cat(probabilities, dim=0).numpy()
    return {
        "labels": label_array,
        "predictions": prediction_array,
        "probabilities": probability_array,
        "metrics": classification_report(label_array, prediction_array, probability_array),
    }


def train_linear_probe(
    *,
    mode: str,
    features_by_split: dict[str, dict[str, Any]],
    num_labels: int,
    task_name: str,
    log_dir: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    set_seed(int(args.seed))
    loaders = build_feature_loaders(features_by_split, args=args)
    feature_dim = int(features_by_split["train"]["features"].shape[1])
    model = nn.Linear(feature_dim, int(num_labels))
    optimizer, scheduler = build_optimizer_and_scheduler(
        optimizer_name=args.optimizer_name,
        params=model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
        train_loader=loaders["train"],
        epochs=args.epochs,
        warmup_ratio=args.warmup_ratio,
        min_lr_ratio=args.min_lr_ratio,
    )
    trainer = CachedFeatureClassificationTrainer(
        model=model,
        optimizer=optimizer,
        log_dir=log_dir,
        device=args.device,
        criterion=nn.CrossEntropyLoss(),
        scheduler=scheduler,
        use_tensorboard=args.tensorboard,
        scheduler_step_per_batch=True,
        enable_progress_bar=args.progress,
    )
    monitor_mode = "min" if args.monitor == "loss" else "max"
    fit_result = trainer.fit(
        train_loader=loaders["train"],
        valid_loader=loaders["valid"],
        test_loader=loaders["test"],
        epochs=args.epochs,
        checkpoint_path=log_dir / "checkpoints" / "best.pt",
        monitor=args.monitor,
        monitor_mode=monitor_mode,
        early_stopping_patience=args.early_stopping_patience,
    )
    load_best_state(model, log_dir / "checkpoints" / "best.pt")
    model.to(args.device)
    predictions = {
        split_name: collect_predictions(model, loaders[split_name], device=args.device)
        for split_name in ("valid", "test")
    }
    write_json(
        log_dir / "metrics.json",
        {
            "status": "ok",
            "mode": mode,
            "task_name": task_name,
            "model_name": args.model_name,
            "feature_dim": feature_dim,
            "num_labels": int(num_labels),
            "fit_result": fit_result,
            "best_predictions": {
                split_name: predictions[split_name]["metrics"]
                for split_name in predictions
            },
        },
    )
    return {
        "mode": mode,
        "fit_result": fit_result,
        "best_valid": predictions["valid"]["metrics"],
        "best_test": predictions["test"]["metrics"],
        "test_predictions": predictions["test"]["predictions"],
        "test_probabilities": predictions["test"]["probabilities"],
    }


def metric_delta(
    serial_metrics: dict[str, float],
    batched_metrics: dict[str, float],
) -> dict[str, float | None]:
    keys = sorted(set(serial_metrics) | set(batched_metrics))
    deltas: dict[str, float | None] = {}
    for key in keys:
        left = serial_metrics.get(key)
        right = batched_metrics.get(key)
        if left is None or right is None or math.isnan(float(left)) or math.isnan(float(right)):
            deltas[key] = None
        else:
            deltas[key] = float(right) - float(left)
    return deltas


def summarize_impact(
    serial_result: dict[str, Any],
    batched_result: dict[str, Any],
    *,
    metric_delta_threshold: float,
    prediction_flip_threshold: float,
) -> dict[str, Any]:
    serial_pred = np.asarray(serial_result["test_predictions"])
    batched_pred = np.asarray(batched_result["test_predictions"])
    if serial_pred.shape != batched_pred.shape:
        raise ValueError("serial 与 batched 的 test prediction shape 不一致。")
    flip_rate = float(np.mean(serial_pred != batched_pred)) if serial_pred.size else 0.0
    deltas = metric_delta(serial_result["best_test"], batched_result["best_test"])
    finite_abs_deltas = [
        abs(float(value))
        for key, value in deltas.items()
        if key in {"accuracy", "acc", "f1", "mcc", "auroc"} and value is not None
    ]
    max_abs_metric_delta = max(finite_abs_deltas) if finite_abs_deltas else 0.0
    status = (
        "negligible"
        if max_abs_metric_delta <= float(metric_delta_threshold)
        and flip_rate <= float(prediction_flip_threshold)
        else "non_negligible"
    )
    return {
        "status": status,
        "test_prediction_flip_rate": flip_rate,
        "max_abs_test_metric_delta": max_abs_metric_delta,
        "test_metric_delta_batched_minus_serial": deltas,
        "metric_delta_threshold": float(metric_delta_threshold),
        "prediction_flip_threshold": float(prediction_flip_threshold),
    }


def validate_task(args: argparse.Namespace, task_name: str) -> dict[str, Any]:
    started_at = time.perf_counter()
    task_log_dir = Path(args.log_root) / task_name
    task_log_dir.mkdir(parents=True, exist_ok=True)

    task_meta = load_processed_task_meta(args.processed_download, task_name)
    if task_meta.get("task_type") != "classification":
        raise ValueError(f"{task_name} 不是分类任务: {task_meta.get('task_type')!r}")
    pair_mode = task_meta.get("input_schema") == "pair_sequence"
    dataframes = load_sampled_dataframes(args, task_name)

    spec = get_model_spec(args.model_name)
    resolved_finetune_method = validate_finetune_method(spec, "frozen_linear_probe")
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
    if not hasattr(model.backbone, "encode_serial") or not hasattr(model.backbone, "encode_batched"):
        raise AttributeError("JanusDNA backbone 必须同时支持 encode_serial 和 encode_batched。")
    model = model.to(args.device).eval()

    serial_features = collect_features(
        model=model,
        dataframes=dataframes,
        mode="serial",
        task_name=task_name,
        pair_mode=pair_mode,
        args=args,
    )
    batched_features = collect_features(
        model=model,
        dataframes=dataframes,
        mode="batched",
        task_name=task_name,
        pair_mode=pair_mode,
        args=args,
    )
    del model
    if str(args.device).startswith("cuda") and torch.cuda.is_available():
        torch.cuda.empty_cache()

    feature_diff = compare_features(serial_features, batched_features)
    serial_train_result = train_linear_probe(
        mode="serial",
        features_by_split=serial_features,
        num_labels=int(task_meta["num_labels"]),
        task_name=task_name,
        log_dir=task_log_dir / "linear_probe_serial",
        args=args,
    )
    batched_train_result = train_linear_probe(
        mode="batched",
        features_by_split=batched_features,
        num_labels=int(task_meta["num_labels"]),
        task_name=task_name,
        log_dir=task_log_dir / "linear_probe_batched",
        args=args,
    )
    impact = summarize_impact(
        serial_train_result,
        batched_train_result,
        metric_delta_threshold=args.metric_delta_threshold,
        prediction_flip_threshold=args.prediction_flip_threshold,
    )
    payload = {
        "status": "ok",
        "task_name": task_name,
        "model_name": args.model_name,
        "raw_download": args.raw_download,
        "processed_download": args.processed_download,
        "sample_strategy": args.sample_strategy,
        "sample_counts": {
            split_name: int(len(dataframes[split_name]))
            for split_name in SPLIT_NAMES
        },
        "config": {
            "extract_batch_size": int(args.extract_batch_size),
            "chunk_forward_batch_size": int(args.chunk_forward_batch_size),
            "epochs": int(args.epochs),
            "batch_size": int(args.batch_size),
            "lr": float(args.lr),
            "weight_decay": float(args.weight_decay),
            "optimizer_name": args.optimizer_name,
            "warmup_ratio": float(args.warmup_ratio),
            "min_lr_ratio": float(args.min_lr_ratio),
            "early_stopping_patience": int(args.early_stopping_patience),
            "seed": int(args.seed),
            "token_readout": args.token_readout,
            "resolved_finetune_method": resolved_finetune_method,
            "feature_train_dtype": args.feature_train_dtype,
        },
        "feature_diff": feature_diff,
        "serial_best_valid": serial_train_result["best_valid"],
        "serial_best_test": serial_train_result["best_test"],
        "batched_best_valid": batched_train_result["best_valid"],
        "batched_best_test": batched_train_result["best_test"],
        "impact": impact,
        "elapsed_seconds": round(time.perf_counter() - started_at, 6),
    }
    write_json(task_log_dir / "summary.json", payload)
    return payload


def main() -> None:
    args = build_parser().parse_args()
    started_at = time.perf_counter()
    log_root = Path(args.log_root)
    log_root.mkdir(parents=True, exist_ok=True)
    task_names = parse_task_names(args.task_name, args.processed_download)
    results: list[dict[str, Any]] = []

    for task_name in task_names:
        try:
            result = validate_task(args, task_name)
        except Exception as exc:
            result = {
                "status": "failed",
                "task_name": task_name,
                "model_name": args.model_name,
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
            }
            write_json(log_root / task_name / "summary.json", result)
            raise
        finally:
            if str(args.device).startswith("cuda") and torch.cuda.is_available():
                torch.cuda.empty_cache()
        results.append(result)

    aggregate = {
        "status": "ok",
        "model_name": args.model_name,
        "tasks": task_names,
        "results": results,
        "elapsed_seconds": round(time.perf_counter() - started_at, 6),
    }
    write_json(log_root / "janusdna_batched_metric_effect_summary.json", aggregate)
    print(json.dumps(aggregate, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
