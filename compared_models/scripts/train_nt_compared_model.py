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
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import (
    available_compared_models,
    build_linear_probe_model,
    build_sequence_dataloader,
    get_model_spec,
    load_processed_task_dataframe,
    load_processed_task_meta,
    resolve_processed_task_names,
    validate_finetune_method,
)
from visualdna.trainer import barrier, cleanup_distributed, maybe_wrap_model, setup_distributed_context
from visualdna.trainer.sequence_classification_trainer import SequenceClassificationTrainer
from visualdna.utils.parsing import print_runtime_locations
from visualdna.utils.seed_utils import set_seed
from visualdna.utils.training import build_optimizer_and_scheduler


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="使用统一 benchmark 协议验证 compared_models 中的 DNA foundation model 线性探针。",
    )
    parser.add_argument(
        "--task-name",
        required=True,
        help="处理后的 benchmark 任务名，例如 enhancers；或使用 all 一次跑完整个 processed root。",
    )
    parser.add_argument(
        "--processed-download",
        required=True,
        help="统一 benchmark 数据目录，例如 /zengxiangxiang/mps/ood_imageDNA/data/low_similarity_sequence_csv_original_parquet/nt。",
    )
    parser.add_argument(
        "--raw-download",
        default=None,
        help="原始或物化数据目录，仅用于日志记录，可为空。",
    )
    parser.add_argument(
        "--log-root",
        required=True,
        help="训练日志根目录；单任务会直接写这里，多任务/多模型时会自动分层建目录。",
    )
    parser.add_argument(
        "--model-name",
        required=True,
        help="对比模型名，例如 grover / ntv2；或使用 all 一次跑完全部已登记模型。",
    )
    parser.add_argument(
        "--checkpoint-override",
        default=None,
        help="可选：覆盖模型默认 checkpoint 路径或 Hugging Face model id。",
    )
    parser.add_argument(
        "--token-readout",
        default="auto",
        choices=["auto", "mean", "last"],
        help="token 级读出方式；auto 使用模型默认策略。",
    )
    parser.add_argument(
        "--finetune-method",
        default="frozen_linear_probe",
        help=(
            "下游微调方式：frozen_linear_probe / full / ia3；"
            "具体以当前模型允许方式为准，Evo2 和 AlphaGenome 仅允许 frozen_linear_probe。"
        ),
    )
    parser.add_argument(
        "--ia3-target-modules",
        default=None,
        help="可选：覆盖 IA3 target_modules，多个模块名用英文逗号分隔。",
    )
    parser.add_argument(
        "--ia3-feedforward-modules",
        default=None,
        help="可选：覆盖 IA3 feedforward_modules，多个模块名用英文逗号分隔。",
    )
    parser.add_argument(
        "--chunk-forward-batch-size",
        type=int,
        default=8,
        help="单个样本内部 chunk 微批大小，控制冻结 backbone 前向时的显存占用。",
    )
    parser.add_argument(
        "--allow-remote-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="本地 checkpoint 加载失败时，是否允许回退到 Hugging Face model id。",
    )
    parser.add_argument("--device", default="auto", help="训练设备，例如 auto / cpu / cuda。")
    parser.add_argument("--dist-backend", default=None, help="分布式后端，默认自动选择。")
    parser.add_argument("--epochs", type=int, default=100, help="最大训练轮数。")
    parser.add_argument("--batch-size", type=int, default=16, help="batch size。")
    parser.add_argument("--num-workers", type=int, default=4, help="DataLoader worker 数。")
    parser.add_argument(
        "--tensorboard",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否写入 TensorBoard 日志。",
    )
    parser.add_argument("--lr", type=float, default=1e-3, help="学习率。")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="权重衰减。")
    parser.add_argument(
        "--optimizer-name",
        default="adamw",
        choices=["adam", "adamw"],
        help="优化器名称。",
    )
    parser.add_argument("--warmup-ratio", type=float, default=0.1, help="warmup 比例。")
    parser.add_argument(
        "--min-lr-ratio",
        type=float,
        default=0.1,
        help="cosine annealing 最终 lr 比例。",
    )
    parser.add_argument(
        "--early-stopping-patience",
        type=int,
        default=30,
        help="验证指标连续多少个 epoch 未提升后提前停止。",
    )
    parser.add_argument("--seed", type=int, default=42, help="随机种子。")
    parser.add_argument(
        "--monitor",
        default="accuracy",
        choices=["accuracy", "acc", "f1", "mcc", "auroc", "loss"],
        help="早停监控指标。",
    )
    parser.add_argument(
        "--max-train-samples",
        type=int,
        default=None,
        help="可选：仅截取前 N 个 train 样本，便于冒烟测试。",
    )
    parser.add_argument(
        "--max-valid-samples",
        type=int,
        default=None,
        help="可选：仅截取前 N 个 valid 样本。",
    )
    parser.add_argument(
        "--max-test-samples",
        type=int,
        default=None,
        help="可选：仅截取前 N 个 test 样本。",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Hugging Face 缓存目录。",
    )
    parser.add_argument(
        "--continue-on-error",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="当某个模型或任务失败时，是否记录错误后继续执行剩余组合。",
    )
    return parser


def trainable_parameters(model: nn.Module) -> list[nn.Parameter]:
    return [parameter for parameter in model.parameters() if parameter.requires_grad]


def resolve_task_names(task_name: str, processed_root: str | Path) -> list[str]:
    if task_name != "all":
        return [task_name]
    task_names = resolve_processed_task_names(processed_root)
    if not task_names:
        raise RuntimeError(f"在 {processed_root} 下未找到任何处理后的 benchmark 任务。")
    return task_names


def resolve_model_names(model_name: str) -> list[str]:
    if model_name != "all":
        return [model_name]
    return available_compared_models(include_unsupported=True)


def subset_dataframe(df: pd.DataFrame, max_samples: int | None) -> pd.DataFrame:
    if max_samples is None or max_samples <= 0:
        return df
    return df.iloc[:int(max_samples)].reset_index(drop=True)


def build_datasets(
    *,
    processed_root: str | Path,
    task_name: str,
    max_train_samples: int | None,
    max_valid_samples: int | None,
    max_test_samples: int | None,
):
    from compared_models.benchmark.data import SequenceBenchmarkDataset

    train_df = subset_dataframe(
        load_processed_task_dataframe(processed_root, task_name, "train"),
        max_train_samples,
    )
    valid_df = subset_dataframe(
        load_processed_task_dataframe(processed_root, task_name, "valid"),
        max_valid_samples,
    )
    test_df = subset_dataframe(
        load_processed_task_dataframe(processed_root, task_name, "test"),
        max_test_samples,
    )
    return {
        "train": SequenceBenchmarkDataset(train_df),
        "valid": SequenceBenchmarkDataset(valid_df),
        "test": SequenceBenchmarkDataset(test_df),
    }


def build_loaders(
    *,
    datasets: dict[str, Any],
    batch_size: int,
    num_workers: int,
    dist_context,
) -> dict[str, Any]:
    pin_memory = torch.cuda.is_available()
    return {
        "train": build_sequence_dataloader(
            datasets["train"],
            batch_size=batch_size,
            num_workers=num_workers,
            shuffle=True,
            is_distributed=dist_context.is_distributed,
            pin_memory=pin_memory,
        ),
        "valid": build_sequence_dataloader(
            datasets["valid"],
            batch_size=batch_size,
            num_workers=num_workers,
            shuffle=False,
            is_distributed=dist_context.is_distributed,
            pin_memory=pin_memory,
        ),
        "test": build_sequence_dataloader(
            datasets["test"],
            batch_size=batch_size,
            num_workers=num_workers,
            shuffle=False,
            is_distributed=dist_context.is_distributed,
            pin_memory=pin_memory,
        ),
    }


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


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


def clear_profile_hooks(model: nn.Module) -> None:
    for module in model.modules():
        module._forward_hooks.clear()
        module._forward_pre_hooks.clear()
        for attr_name in ("total_ops", "total_params"):
            if hasattr(module, attr_name):
                delattr(module, attr_name)


def safe_profile_model_from_loader(
    model: nn.Module,
    train_loader,
    *,
    device: str,
    output_path: Path,
    pair_mode: bool,
) -> dict[str, Any]:
    try:
        from visualdna.utils.model_profile import profile_model_from_loader
    except ModuleNotFoundError as exc:
        payload = {
            "status": "skipped",
            "reason": f"missing_optional_dependency: {exc}",
        }
        write_json(output_path, payload)
        return payload

    try:
        return profile_model_from_loader(
            model,
            train_loader,
            device=device,
            pair_mode=pair_mode,
            output_path=output_path,
        )
    except Exception as exc:
        payload = {
            "status": "partial",
            "reason": f"profile_error: {type(exc).__name__}: {exc}",
        }
        write_json(output_path, payload)
        return payload
    finally:
        clear_profile_hooks(model)


def write_run_config(
    args: argparse.Namespace,
    *,
    task_name: str,
    model_name: str,
    spec,
    meta: dict[str, Any],
    log_dir: Path,
    dist_context,
) -> None:
    if not dist_context.is_main_process:
        return
    run_config = vars(args).copy()
    run_config["task_name"] = task_name
    run_config["model_name"] = model_name
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
    run_config["resolved_finetune_method"] = validate_finetune_method(spec, args.finetune_method)
    run_config["dataset_meta_path"] = str(log_dir / "dataset_meta.json")
    run_config["distributed"] = dist_context.is_distributed
    run_config["rank"] = dist_context.rank
    run_config["world_size"] = dist_context.world_size
    run_config["num_labels"] = meta["num_labels"]
    run_config["input_schema"] = meta.get("input_schema", "single_sequence")
    write_json(log_dir / "run_config.json", run_config)
    write_json(log_dir / "dataset_meta.json", meta)


def runtime_payload(
    *,
    start_time: float,
    dist_context,
) -> dict[str, Any]:
    payload = {
        "wall_time_seconds": round(time.perf_counter() - start_time, 6),
        "device": dist_context.device,
        "distributed": dist_context.is_distributed,
        "rank": dist_context.rank,
        "world_size": dist_context.world_size,
    }
    device_index = cuda_device_index(dist_context.device)
    if device_index is not None:
        try:
            max_allocated = torch.cuda.max_memory_allocated(device_index)
            max_reserved = torch.cuda.max_memory_reserved(device_index)
            properties = torch.cuda.get_device_properties(device_index)
        except RuntimeError as exc:
            payload["cuda_memory_stats_error"] = f"{type(exc).__name__}: {exc}"
        else:
            payload["max_cuda_memory_allocated_bytes"] = int(max_allocated)
            payload["max_cuda_memory_allocated_mib"] = round(float(max_allocated) / (1024 ** 2), 3)
            payload["max_cuda_memory_reserved_bytes"] = int(max_reserved)
            payload["max_cuda_memory_reserved_mib"] = round(float(max_reserved) / (1024 ** 2), 3)
            payload["cuda_device_name"] = properties.name
            payload["cuda_device_total_memory_bytes"] = int(properties.total_memory)
            payload["cuda_device_total_memory_mib"] = round(float(properties.total_memory) / (1024 ** 2), 3)
    return payload


def run_single(
    args: argparse.Namespace,
    *,
    task_name: str,
    model_name: str,
    log_dir: Path,
    dist_context,
) -> dict[str, Any]:
    start_time = time.perf_counter()
    cuda_memory_stats_error = safe_reset_peak_memory_stats(dist_context.device)

    set_seed(args.seed)
    spec = get_model_spec(model_name)
    resolved_finetune_method = validate_finetune_method(spec, args.finetune_method)
    meta = load_processed_task_meta(args.processed_download, task_name)
    if meta.get("task_type") != "classification":
        raise ValueError(f"当前脚本仅支持分类任务，收到 {task_name} 的 task_type={meta.get('task_type')!r}")
    pair_mode = meta.get("input_schema") == "pair_sequence"
    write_run_config(
        args,
        task_name=task_name,
        model_name=model_name,
        spec=spec,
        meta=meta,
        log_dir=log_dir,
        dist_context=dist_context,
    )
    if dist_context.is_main_process:
        print_runtime_locations(
            task_name=task_name,
            raw_download=args.raw_download,
            processed_download=args.processed_download,
            log_dir=log_dir,
            hf_endpoint=None,
            hf_cache_dir=args.cache_dir,
            decoder_name="LinearClassifier",
            decoder_hf_model_id=spec.hf_model_id,
        )

    datasets = build_datasets(
        processed_root=args.processed_download,
        task_name=task_name,
        max_train_samples=args.max_train_samples,
        max_valid_samples=args.max_valid_samples,
        max_test_samples=args.max_test_samples,
    )
    loaders = build_loaders(
        datasets=datasets,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        dist_context=dist_context,
    )
    model = build_linear_probe_model(
        model_name=model_name,
        num_labels=int(meta["num_labels"]),
        token_readout=args.token_readout,
        finetune_method=resolved_finetune_method,
        cache_dir=args.cache_dir,
        allow_remote_fallback=args.allow_remote_fallback,
        chunk_forward_batch_size=args.chunk_forward_batch_size,
        checkpoint_override=args.checkpoint_override,
        pair_mode=pair_mode,
        ia3_target_modules=args.ia3_target_modules,
        ia3_feedforward_modules=args.ia3_feedforward_modules,
    )

    model_profile = None
    if dist_context.is_main_process:
        model_profile = safe_profile_model_from_loader(
            model,
            loaders["train"],
            device=dist_context.device,
            output_path=log_dir / "model_profile.json",
            pair_mode=pair_mode,
        )
    if dist_context.is_distributed:
        barrier()
    model = maybe_wrap_model(model, dist_context)
    optimizer_params = trainable_parameters(model)
    if not optimizer_params:
        raise RuntimeError(f"{model_name} 在 {resolved_finetune_method} 模式下没有可训练参数。")
    optimizer, scheduler = build_optimizer_and_scheduler(
        optimizer_name=args.optimizer_name,
        params=optimizer_params,
        lr=args.lr,
        weight_decay=args.weight_decay,
        train_loader=loaders["train"],
        epochs=args.epochs,
        warmup_ratio=args.warmup_ratio,
        min_lr_ratio=args.min_lr_ratio,
    )
    trainer = SequenceClassificationTrainer(
        model=model,
        optimizer=optimizer,
        log_dir=str(log_dir),
        device=dist_context.device,
        criterion=nn.CrossEntropyLoss(),
        scheduler=scheduler,
        is_main_process=dist_context.is_main_process,
        use_tensorboard=args.tensorboard,
        scheduler_step_per_batch=True,
    )
    monitor = args.monitor
    monitor_mode = "min" if monitor == "loss" else "max"
    fit_result = trainer.fit(
        train_loader=loaders["train"],
        valid_loader=loaders["valid"],
        test_loader=loaders["test"],
        epochs=args.epochs,
        checkpoint_path=log_dir / "checkpoints" / "best.pt",
        monitor=monitor,
        monitor_mode=monitor_mode,
        early_stopping_patience=args.early_stopping_patience,
    )
    runtime = runtime_payload(start_time=start_time, dist_context=dist_context)
    if cuda_memory_stats_error is not None:
        runtime["cuda_memory_stats_reset_error"] = cuda_memory_stats_error
    if dist_context.is_main_process:
        write_json(log_dir / "runtime.json", runtime)

    summary = trainer.report_run_summary(
        fit_result,
        status="ok",
        task_name=task_name,
        model_name=model_name,
        model_display_name=spec.display_name,
        log_dir=str(log_dir),
        num_labels=meta["num_labels"],
        device=dist_context.device,
        distributed=dist_context.is_distributed,
        rank=dist_context.rank,
        world_size=dist_context.world_size,
        token_readout=args.token_readout,
        finetune_method=resolved_finetune_method,
        model_profile=model_profile,
        runtime=runtime,
        monitor=monitor,
        chunk_length=spec.chunk_length,
        chunk_unit=spec.chunk_unit,
        input_schema=meta.get("input_schema", "single_sequence"),
    )
    barrier()
    return summary


def skipped_summary(
    *,
    task_name: str,
    model_name: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "task_name": task_name,
        "model_name": model_name,
        "status": "skipped",
        "reason": reason,
    }


def failed_summary(
    *,
    task_name: str,
    model_name: str,
    exc: Exception,
) -> dict[str, Any]:
    return {
        "task_name": task_name,
        "model_name": model_name,
        "status": "failed",
        "error_type": type(exc).__name__,
        "error_message": str(exc),
        "traceback": traceback.format_exc(),
    }


def main() -> None:
    args = build_parser().parse_args()
    dist_context = setup_distributed_context(
        args.device,
        backend=args.dist_backend,
    )
    task_names = resolve_task_names(args.task_name, args.processed_download)
    model_names = resolve_model_names(args.model_name)
    base_log_root = Path(args.log_root)
    all_results: list[dict[str, Any]] = []

    try:
        for task_name in task_names:
            for model_name in model_names:
                spec = get_model_spec(model_name)
                task_log_dir = base_log_root
                if len(task_names) > 1 or len(model_names) > 1:
                    task_log_dir = base_log_root / task_name / model_name

                if spec.loader_kind == "unsupported":
                    summary = skipped_summary(
                        task_name=task_name,
                        model_name=model_name,
                        reason=spec.unsupported_reason or "当前脚本未接入该模型。",
                    )
                    all_results.append(summary)
                    if dist_context.is_main_process:
                        task_log_dir.mkdir(parents=True, exist_ok=True)
                        write_json(task_log_dir / "metrics.json", summary)
                    continue

                try:
                    summary = run_single(
                        args,
                        task_name=task_name,
                        model_name=model_name,
                        log_dir=task_log_dir,
                        dist_context=dist_context,
                    )
                    summary["status"] = "ok"
                    all_results.append(summary)
                except Exception as exc:
                    summary = failed_summary(
                        task_name=task_name,
                        model_name=model_name,
                        exc=exc,
                    )
                    all_results.append(summary)
                    if dist_context.is_main_process:
                        task_log_dir.mkdir(parents=True, exist_ok=True)
                        write_json(task_log_dir / "metrics.json", summary)
                    if not args.continue_on_error:
                        raise

        if dist_context.is_main_process and (len(task_names) > 1 or len(model_names) > 1):
            base_log_root.mkdir(parents=True, exist_ok=True)
            write_json(
                base_log_root / "all_tasks_metrics.json",
                {
                    "mode": "benchmark",
                    "tasks": task_names,
                    "models": model_names,
                    "results": all_results,
                },
            )
    finally:
        cleanup_distributed()


if __name__ == "__main__":
    main()
