from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader, Dataset

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import get_model_spec, load_processed_task_meta
from visualdna.metrics import classification_report
from visualdna.trainer.base_trainer import BaseTrainer
from visualdna.utils.parsing import print_runtime_locations
from visualdna.utils.seed_utils import set_seed
from visualdna.utils.training import build_optimizer_and_scheduler


SPLIT_NAMES = ("train", "valid", "test")
SPLIT_TO_ID = {name: index for index, name in enumerate(SPLIT_NAMES)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从 frozen backbone embedding cache 训练 NT 对比模型线性探针。",
    )
    parser.add_argument("--task-name", required=True, help="处理后的 benchmark 任务名。")
    parser.add_argument("--processed-download", required=True, help="处理后的 NT benchmark 根目录。")
    parser.add_argument("--raw-download", default=None, help="原始或物化数据目录，仅用于日志记录。")
    parser.add_argument("--log-root", required=True, help="训练日志目录。")
    parser.add_argument("--cache-path", required=True, help="extract_nt_embedding_cache.py 生成的 .pt cache。")
    parser.add_argument("--model-name", required=True, help="对比模型名。")
    parser.add_argument("--device", default="cuda", help="训练设备，例如 cuda / cpu。")
    parser.add_argument("--epochs", type=int, default=50, help="最大训练轮数。")
    parser.add_argument("--batch-size", type=int, default=32, help="线性头训练 batch size。")
    parser.add_argument("--num-workers", type=int, default=4, help="DataLoader worker 数。")
    parser.add_argument("--lr", type=float, default=1e-3, help="学习率。")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="权重衰减。")
    parser.add_argument("--optimizer-name", default="adamw", choices=["adam", "adamw"], help="优化器名称。")
    parser.add_argument("--warmup-ratio", type=float, default=0.1, help="warmup 比例。")
    parser.add_argument("--min-lr-ratio", type=float, default=0.1, help="cosine annealing 最终 lr 比例。")
    parser.add_argument("--early-stopping-patience", type=int, default=15, help="早停 patience。")
    parser.add_argument("--seed", type=int, default=42, help="随机种子。")
    parser.add_argument(
        "--monitor",
        default="accuracy",
        choices=["accuracy", "acc", "f1", "mcc", "auroc", "loss"],
        help="早停监控指标。",
    )
    parser.add_argument(
        "--tensorboard",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否写入 TensorBoard 日志。",
    )
    parser.add_argument(
        "--feature-train-dtype",
        default="float32",
        choices=["float32", "float16"],
        help="送入线性头训练时的 feature dtype；默认 float32。",
    )
    parser.add_argument(
        "--eval-only",
        action="store_true",
        help="只加载 --checkpoint-path 中的 clean 线性探针并评估当前 cache 的 test split。",
    )
    parser.add_argument(
        "--checkpoint-path",
        default=None,
        help="eval-only 模式下要加载的 clean linear probe best.pt。",
    )
    parser.add_argument(
        "--checkpoint-strict",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="加载 eval-only checkpoint 时是否严格匹配线性头权重键。",
    )
    parser.add_argument(
        "--eval-output-name",
        default="metrics.json",
        help="eval-only 模式下写入 log-root 的指标文件名。",
    )
    return parser


def torch_load(path: Path) -> dict[str, Any]:
    try:
        payload = torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        payload = torch.load(path, map_location="cpu")
    if not isinstance(payload, dict):
        raise TypeError(f"cache payload 不是 dict: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def load_linear_probe_checkpoint(
    model: nn.Module,
    checkpoint_path: Path,
    *,
    strict: bool,
) -> dict[str, Any]:
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"找不到 checkpoint: {checkpoint_path}")
    payload = torch_load(checkpoint_path)
    state_dict = payload.get("model_state_dict")
    if not isinstance(state_dict, dict):
        raise ValueError(f"checkpoint 中缺少 model_state_dict: {checkpoint_path}")
    load_result = model.load_state_dict(state_dict, strict=bool(strict))
    return {
        "checkpoint_path": str(checkpoint_path),
        "checkpoint_epoch": payload.get("epoch"),
        "checkpoint_monitor": payload.get("monitor"),
        "checkpoint_score": payload.get("score"),
        "missing_keys": list(getattr(load_result, "missing_keys", [])),
        "unexpected_keys": list(getattr(load_result, "unexpected_keys", [])),
        "strict": bool(strict),
    }


def build_eval_only_result(test_metrics: dict[str, float]) -> dict[str, Any]:
    epoch_result = {
        "epoch": 0,
        "train": None,
        "valid": None,
        "test": test_metrics,
    }
    return {
        "best": epoch_result,
        "history": [epoch_result],
        "epochs_ran": 0,
        "eval_only": True,
    }


def cuda_device_index(device: str) -> int | None:
    if not str(device).startswith("cuda") or not torch.cuda.is_available() or torch.cuda.device_count() <= 0:
        return None
    index = torch.device(device).index
    if index is None:
        index = torch.cuda.current_device()
    if index < 0 or index >= torch.cuda.device_count():
        return None
    return int(index)


def safe_reset_peak_memory_stats(device: str) -> str | None:
    device_index = cuda_device_index(device)
    if device_index is None:
        return None
    try:
        torch.cuda.reset_peak_memory_stats(device_index)
    except RuntimeError as exc:
        return f"{type(exc).__name__}: {exc}"
    return None


def runtime_payload(*, start_time: float, device: str) -> dict[str, Any]:
    payload = {
        "wall_time_seconds": round(time.perf_counter() - start_time, 6),
        "device": device,
        "distributed": False,
        "rank": 0,
        "world_size": 1,
        "cached_embedding": True,
    }
    device_index = cuda_device_index(device)
    if device_index is None:
        return payload
    try:
        max_allocated = torch.cuda.max_memory_allocated(device_index)
        max_reserved = torch.cuda.max_memory_reserved(device_index)
        properties = torch.cuda.get_device_properties(device_index)
    except RuntimeError as exc:
        payload["cuda_memory_stats_error"] = f"{type(exc).__name__}: {exc}"
        return payload
    payload.update(
        {
            "max_cuda_memory_allocated_bytes": int(max_allocated),
            "max_cuda_memory_allocated_mib": round(float(max_allocated) / (1024 ** 2), 3),
            "max_cuda_memory_reserved_bytes": int(max_reserved),
            "max_cuda_memory_reserved_mib": round(float(max_reserved) / (1024 ** 2), 3),
            "cuda_device_name": properties.name,
            "cuda_device_total_memory_bytes": int(properties.total_memory),
            "cuda_device_total_memory_mib": round(float(properties.total_memory) / (1024 ** 2), 3),
        }
    )
    return payload


class CachedFeatureDataset(Dataset):
    def __init__(self, features: torch.Tensor, labels: torch.Tensor) -> None:
        if features.ndim != 2:
            raise ValueError(f"features 必须是二维张量，收到 shape={tuple(features.shape)}")
        if labels.ndim != 1:
            raise ValueError(f"labels 必须是一维张量，收到 shape={tuple(labels.shape)}")
        if int(features.shape[0]) != int(labels.shape[0]):
            raise ValueError("features 与 labels 样本数不一致。")
        self.features = features.contiguous()
        self.labels = labels.long().contiguous()

    def __len__(self) -> int:
        return int(self.labels.numel())

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        index = int(index)
        return {
            "features": self.features[index],
            "labels": self.labels[index],
        }


class CachedFeatureClassificationTrainer(BaseTrainer):
    def _init_epoch_state(self) -> dict[str, Any]:
        state = super()._init_epoch_state()
        state.update(
            {
                "total_correct": 0,
                "labels": [],
                "predictions": [],
                "probabilities": [],
            }
        )
        return state

    def _compute_step(self, batch: dict[str, Any]) -> dict[str, Any]:
        features = batch["features"]
        labels = batch["labels"]
        if not isinstance(features, torch.Tensor) or not isinstance(labels, torch.Tensor):
            raise TypeError("cached feature trainer 需要 Tensor 类型的 features / labels。")
        model_dtype = next(self.model.parameters()).dtype
        features = features.to(self.device, dtype=model_dtype, non_blocking=True)
        labels = labels.view(-1).long().to(self.device, non_blocking=True)
        logits = self.model(features)
        loss = self.criterion(logits, labels)
        predictions = logits.argmax(dim=-1)
        probabilities = torch.softmax(logits.detach(), dim=-1)
        return {
            "loss": loss,
            "batch_size": labels.numel(),
            "correct": int((predictions == labels).sum().item()),
            "labels": labels.detach().cpu(),
            "predictions": predictions.detach().cpu(),
            "probabilities": probabilities.detach().cpu(),
        }

    def _update_epoch_state(self, state: dict[str, Any], step_result: dict[str, Any]) -> None:
        super()._update_epoch_state(state, step_result)
        state["total_correct"] += int(step_result["correct"])
        state["labels"].append(step_result["labels"])
        state["predictions"].append(step_result["predictions"])
        state["probabilities"].append(step_result["probabilities"])

    def _progress_postfix(self, state: dict[str, Any]) -> dict[str, str]:
        total_samples = max(int(state["total_samples"]), 1)
        return {
            "loss": f"{state['total_loss'] / total_samples:.4f}",
            "acc": f"{state['total_correct'] / total_samples:.4f}",
        }

    def _finalize_epoch_metrics(self, state: dict[str, Any]) -> dict[str, float]:
        total_loss, total_samples, total_correct = self._all_reduce_scalars(
            float(state["total_loss"]),
            float(state["total_samples"]),
            float(state["total_correct"]),
        )
        labels = torch.cat(state["labels"], dim=0).cpu().numpy()
        predictions = torch.cat(state["predictions"], dim=0).cpu().numpy()
        probabilities = torch.cat(state["probabilities"], dim=0).cpu().numpy()
        gathered_labels = self._all_gather_object(labels)
        gathered_predictions = self._all_gather_object(predictions)
        gathered_probabilities = self._all_gather_object(probabilities)
        labels = np.concatenate(gathered_labels, axis=0)
        predictions = np.concatenate(gathered_predictions, axis=0)
        probabilities = np.concatenate(gathered_probabilities, axis=0)
        metrics = classification_report(labels, predictions, probabilities)
        metrics["loss"] = float(total_loss) / max(int(total_samples), 1)
        metrics["accuracy"] = float(total_correct) / max(int(total_samples), 1)
        metrics["acc"] = metrics["accuracy"]
        return metrics


def split_tensors(
    payload: dict[str, Any],
    *,
    split_name: str,
    feature_dtype: torch.dtype,
) -> tuple[torch.Tensor, torch.Tensor]:
    features = payload["features"]
    labels = payload["labels"]
    split_ids = payload["split_ids"]
    if not isinstance(features, torch.Tensor) or not isinstance(labels, torch.Tensor) or not isinstance(split_ids, torch.Tensor):
        raise TypeError("cache 必须包含 Tensor 类型的 features / labels / split_ids。")
    mask = split_ids.long() == SPLIT_TO_ID[split_name]
    return features[mask].to(dtype=feature_dtype).contiguous(), labels[mask].long().contiguous()


def build_loader(
    payload: dict[str, Any],
    *,
    split_name: str,
    batch_size: int,
    num_workers: int,
    feature_dtype: torch.dtype,
    pin_memory: bool,
    shuffle: bool,
) -> DataLoader:
    features, labels = split_tensors(payload, split_name=split_name, feature_dtype=feature_dtype)
    dataset = CachedFeatureDataset(features, labels)
    return DataLoader(
        dataset,
        batch_size=int(batch_size),
        shuffle=bool(shuffle),
        num_workers=int(num_workers),
        pin_memory=bool(pin_memory),
        persistent_workers=bool(num_workers > 0),
        drop_last=False,
    )


def build_model_profile(model: nn.Module, *, feature_dim: int, num_labels: int, cache_meta: dict[str, Any]) -> dict[str, Any]:
    total_params = sum(parameter.numel() for parameter in model.parameters())
    trainable_params = sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)
    return {
        "status": "cached_linear_probe",
        "cached_embedding": True,
        "feature_dim": int(feature_dim),
        "num_labels": int(num_labels),
        "total_params": int(total_params),
        "total_params_m": round(float(total_params) / 1_000_000, 6),
        "trainable_params": int(trainable_params),
        "trainable_params_m": round(float(trainable_params) / 1_000_000, 6),
        "backbone_model_name": cache_meta.get("model_name"),
    }


def validate_cache(args: argparse.Namespace, payload: dict[str, Any]) -> dict[str, Any]:
    meta = payload.get("meta")
    if not isinstance(meta, dict):
        raise ValueError("cache 中缺少 meta。")
    if meta.get("task_name") != args.task_name:
        raise ValueError(f"cache task_name 不匹配: {meta.get('task_name')} != {args.task_name}")
    if meta.get("model_name") != args.model_name:
        raise ValueError(f"cache model_name 不匹配: {meta.get('model_name')} != {args.model_name}")
    if meta.get("finetune_method") != "frozen_linear_probe":
        raise ValueError("cache 不是 frozen_linear_probe embedding。")
    return meta


def write_run_config(
    args: argparse.Namespace,
    *,
    cache_meta: dict[str, Any],
    dataset_meta: dict[str, Any],
    log_dir: Path,
) -> None:
    spec = get_model_spec(args.model_name)
    run_config = vars(args).copy()
    run_config["resolved_finetune_method"] = "frozen_linear_probe"
    run_config["token_readout"] = cache_meta.get("token_readout")
    run_config["cached_embedding"] = True
    run_config["cache_meta"] = cache_meta
    run_config["model_display_name"] = spec.display_name
    run_config["model_spec"] = {
        "loader_kind": spec.loader_kind,
        "chunk_length": spec.chunk_length,
        "chunk_unit": spec.chunk_unit,
        "default_token_readout": spec.default_token_readout,
        "local_checkpoint": spec.local_checkpoint,
        "hf_model_id": spec.hf_model_id,
        "allowed_finetune_methods": list(spec.allowed_finetune_methods),
        "ia3_target_modules": list(spec.ia3_target_modules),
        "ia3_feedforward_modules": list(spec.ia3_feedforward_modules),
    }
    run_config["num_labels"] = dataset_meta["num_labels"]
    run_config["input_schema"] = dataset_meta.get("input_schema", "single_sequence")
    write_json(log_dir / "run_config.json", run_config)
    write_json(log_dir / "dataset_meta.json", dataset_meta)


def run(args: argparse.Namespace) -> dict[str, Any]:
    started_at = time.perf_counter()
    cuda_memory_stats_error = safe_reset_peak_memory_stats(args.device)
    set_seed(args.seed)

    cache_path = Path(args.cache_path)
    payload = torch_load(cache_path)
    cache_meta = validate_cache(args, payload)
    dataset_meta = load_processed_task_meta(args.processed_download, args.task_name)
    log_dir = Path(args.log_root)
    log_dir.mkdir(parents=True, exist_ok=True)
    write_run_config(args, cache_meta=cache_meta, dataset_meta=dataset_meta, log_dir=log_dir)

    print_runtime_locations(
        task_name=args.task_name,
        raw_download=args.raw_download,
        processed_download=args.processed_download,
        log_dir=log_dir,
        hf_endpoint=None,
        hf_cache_dir=None,
        decoder_name="CachedLinearClassifier",
        decoder_hf_model_id=cache_meta.get("model_spec", {}).get("hf_model_id"),
    )

    feature_dtype = torch.float16 if args.feature_train_dtype == "float16" else torch.float32
    pin_memory = torch.cuda.is_available() and str(args.device).startswith("cuda")
    loaders = {
        "train": build_loader(
            payload,
            split_name="train",
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            feature_dtype=feature_dtype,
            pin_memory=pin_memory,
            shuffle=True,
        ),
        "valid": build_loader(
            payload,
            split_name="valid",
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            feature_dtype=feature_dtype,
            pin_memory=pin_memory,
            shuffle=False,
        ),
        "test": build_loader(
            payload,
            split_name="test",
            batch_size=args.batch_size,
            num_workers=args.num_workers,
            feature_dtype=feature_dtype,
            pin_memory=pin_memory,
            shuffle=False,
        ),
    }

    feature_dim = int(payload["features"].shape[1])
    num_labels = int(cache_meta["num_labels"])
    model = nn.Linear(feature_dim, num_labels)
    profile = build_model_profile(
        model,
        feature_dim=feature_dim,
        num_labels=num_labels,
        cache_meta=cache_meta,
    )
    write_json(log_dir / "model_profile.json", profile)

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
        log_dir=str(log_dir),
        device=args.device,
        criterion=nn.CrossEntropyLoss(),
        scheduler=scheduler,
        use_tensorboard=args.tensorboard,
        scheduler_step_per_batch=True,
    )
    monitor_mode = "min" if args.monitor == "loss" else "max"
    spec = get_model_spec(args.model_name)

    if args.eval_only:
        if args.checkpoint_path is None:
            raise ValueError("启用 --eval-only 时必须提供 --checkpoint-path。")
        checkpoint_path = Path(args.checkpoint_path)
        try:
            checkpoint_meta = load_linear_probe_checkpoint(
                model,
                checkpoint_path,
                strict=args.checkpoint_strict,
            )
            test_metrics = trainer.evaluate(loaders["test"], "test")
            fit_result = build_eval_only_result(test_metrics)
            fit_result["monitor"] = args.monitor
            fit_result["monitor_mode"] = monitor_mode
            runtime = runtime_payload(start_time=started_at, device=args.device)
            if cuda_memory_stats_error is not None:
                runtime["cuda_memory_stats_reset_error"] = cuda_memory_stats_error
            write_json(log_dir / "runtime.json", runtime)
            return trainer.report_run_summary(
                fit_result,
                output_name=args.eval_output_name,
                status="ok",
                task_name=args.task_name,
                model_name=args.model_name,
                model_display_name=spec.display_name,
                log_dir=str(log_dir),
                num_labels=num_labels,
                device=args.device,
                distributed=False,
                rank=0,
                world_size=1,
                token_readout=cache_meta.get("token_readout"),
                finetune_method="frozen_linear_probe",
                cached_embedding=True,
                cache_path=str(cache_path),
                checkpoint_path=str(checkpoint_path),
                checkpoint_meta=checkpoint_meta,
                model_profile=profile,
                runtime=runtime,
                monitor=args.monitor,
                chunk_length=cache_meta.get("model_spec", {}).get("chunk_length"),
                chunk_unit=cache_meta.get("model_spec", {}).get("chunk_unit"),
                input_schema=cache_meta.get("input_schema", "single_sequence"),
            )
        finally:
            trainer.close()

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
    runtime = runtime_payload(start_time=started_at, device=args.device)
    if cuda_memory_stats_error is not None:
        runtime["cuda_memory_stats_reset_error"] = cuda_memory_stats_error
    write_json(log_dir / "runtime.json", runtime)

    return trainer.report_run_summary(
        fit_result,
        status="ok",
        task_name=args.task_name,
        model_name=args.model_name,
        model_display_name=spec.display_name,
        log_dir=str(log_dir),
        num_labels=num_labels,
        device=args.device,
        distributed=False,
        rank=0,
        world_size=1,
        token_readout=cache_meta.get("token_readout"),
        finetune_method="frozen_linear_probe",
        cached_embedding=True,
        cache_path=str(cache_path),
        model_profile=profile,
        runtime=runtime,
        monitor=args.monitor,
        chunk_length=cache_meta.get("model_spec", {}).get("chunk_length"),
        chunk_unit=cache_meta.get("model_spec", {}).get("chunk_unit"),
        input_schema=cache_meta.get("input_schema", "single_sequence"),
    )


def main() -> None:
    args = build_parser().parse_args()
    try:
        run(args)
    except Exception as exc:
        log_dir = Path(args.log_root)
        log_dir.mkdir(parents=True, exist_ok=True)
        write_json(
            log_dir / "metrics.json",
            {
                "task_name": args.task_name,
                "model_name": args.model_name,
                "status": "failed",
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
                "cached_embedding": True,
                "cache_path": args.cache_path,
            },
        )
        raise


if __name__ == "__main__":
    main()
