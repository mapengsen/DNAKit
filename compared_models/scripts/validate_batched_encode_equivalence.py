from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import torch
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
    load_processed_task_dataframe,
    load_processed_task_meta,
    validate_finetune_method,
    get_model_spec,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="比较 compared model 的逐条 encode 与真 batch encode 是否数值一致。",
    )
    parser.add_argument("--task-name", required=True, help="处理后的 benchmark 任务名。")
    parser.add_argument(
        "--processed-download",
        required=True,
        help="统一 benchmark 数据根目录，例如 /zengxiangxiang/mps/ood_imageDNA/data/low_similarity_sequence_csv_original_parquet/nt。",
    )
    parser.add_argument("--log-root", required=True, help="验证结果输出目录。")
    parser.add_argument("--model-name", required=True, help="模型名。")
    parser.add_argument(
        "--split-name",
        default="valid",
        choices=["train", "valid", "test"],
        help="从哪个 split 抽取 batch。",
    )
    parser.add_argument("--batch-size", type=int, default=8, help="用于一致性验证的序列数。")
    parser.add_argument("--num-workers", type=int, default=0, help="DataLoader worker 数。")
    parser.add_argument(
        "--token-readout",
        default="auto",
        choices=["auto", "mean", "last"],
        help="token 级读出方式。",
    )
    parser.add_argument(
        "--finetune-method",
        default="frozen_linear_probe",
        help="下游微调方式。",
    )
    parser.add_argument(
        "--chunk-forward-batch-size",
        type=int,
        default=8,
        help="真 batch encode 时一次 forward 的 chunk 数上限。",
    )
    parser.add_argument(
        "--embedding-atol",
        type=float,
        default=1e-3,
        help="embedding 最大绝对误差阈值。",
    )
    parser.add_argument(
        "--logit-atol",
        type=float,
        default=1e-3,
        help="logit 最大绝对误差阈值。",
    )
    parser.add_argument(
        "--allow-remote-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="本地 checkpoint 不完整时，是否回退到 Hugging Face。",
    )
    parser.add_argument("--checkpoint-override", default=None, help="覆盖默认 checkpoint。")
    parser.add_argument("--cache-dir", default=None, help="Hugging Face 缓存目录。")
    parser.add_argument("--device", default="cuda", help="运行设备。")
    parser.add_argument(
        "--fail-on-mismatch",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="误差超过阈值时是否返回非零退出码。",
    )
    parser.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否显示 tqdm 阶段进度条。",
    )
    return parser


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _max_memory_payload(device: str) -> dict[str, Any]:
    if not str(device).startswith("cuda") or not torch.cuda.is_available():
        return {}
    device_index = torch.device(device).index
    if device_index is None:
        device_index = torch.cuda.current_device()
    return {
        "max_cuda_memory_allocated_mib": round(
            torch.cuda.max_memory_allocated(device_index) / (1024 ** 2),
            3,
        ),
        "max_cuda_memory_reserved_mib": round(
            torch.cuda.max_memory_reserved(device_index) / (1024 ** 2),
            3,
        ),
    }


def _move_features_to_head(model: torch.nn.Module, features: torch.Tensor) -> torch.Tensor:
    head = model.head
    if features.device != head.weight.device or features.dtype != head.weight.dtype:
        features = features.to(device=head.weight.device, dtype=head.weight.dtype)
    return features


def _logits_from_features(model: torch.nn.Module, features: torch.Tensor) -> torch.Tensor:
    return model.head(_move_features_to_head(model, features))


def _paired_logits_from_features(
    model: torch.nn.Module,
    ref_features: torch.Tensor,
    alt_features: torch.Tensor,
) -> torch.Tensor:
    features = torch.cat(
        [
            ref_features,
            alt_features,
            alt_features - ref_features,
        ],
        dim=-1,
    )
    return _logits_from_features(model, features)


def _diff_payload(left: torch.Tensor, right: torch.Tensor) -> dict[str, Any]:
    left_cpu = left.detach().float().cpu()
    right_cpu = right.detach().float().cpu()
    diff = (left_cpu - right_cpu).abs()
    return {
        "shape": list(left_cpu.shape),
        "max_abs_diff": float(diff.max().item()) if diff.numel() else 0.0,
        "mean_abs_diff": float(diff.mean().item()) if diff.numel() else 0.0,
    }


def main() -> None:
    args = build_parser().parse_args()
    started_at = time.perf_counter()
    log_root = Path(args.log_root)
    output_path = log_root / args.model_name / "batched_encode_equivalence.json"

    if str(args.device).startswith("cuda") and torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()

    progress = tqdm(
        total=6,
        desc=f"{args.model_name} equivalence",
        disable=not args.progress,
    )
    try:
        task_meta = load_processed_task_meta(args.processed_download, args.task_name)
        pair_mode = task_meta.get("input_schema") == "pair_sequence"
        dataframe = load_processed_task_dataframe(args.processed_download, args.task_name, args.split_name)
        dataframe = dataframe.iloc[: int(args.batch_size)].reset_index(drop=True)
        dataset = SequenceBenchmarkDataset(dataframe)
        loader = build_sequence_dataloader(
            dataset,
            batch_size=int(args.batch_size),
            num_workers=int(args.num_workers),
            shuffle=False,
            is_distributed=False,
            pin_memory=False,
        )
        batch = next(iter(loader))
        progress.update()

        spec = get_model_spec(args.model_name)
        resolved_finetune_method = validate_finetune_method(spec, args.finetune_method)
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
        model = model.to(args.device).eval()
        progress.update()

        if not hasattr(model.backbone, "encode_serial"):
            raise AttributeError(f"{args.model_name} backbone 缺少 encode_serial。")
        if not hasattr(model.backbone, "encode_batched"):
            raise AttributeError(f"{args.model_name} backbone 缺少 encode_batched。")

        with torch.no_grad():
            if pair_mode:
                serial_ref = model.backbone.encode_serial(
                    batch["sequences_ref"],
                    token_readout=model.token_readout,
                )
                serial_alt = model.backbone.encode_serial(
                    batch["sequences_alt"],
                    token_readout=model.token_readout,
                )
                serial_logits = _paired_logits_from_features(model, serial_ref, serial_alt)
                serial_features = torch.cat([serial_ref, serial_alt], dim=-1)
            else:
                serial_features = model.backbone.encode_serial(
                    batch["sequences"],
                    token_readout=model.token_readout,
                )
                serial_logits = _logits_from_features(model, serial_features)
        progress.update()

        with torch.no_grad():
            if pair_mode:
                batched_ref = model.backbone.encode_batched(
                    batch["sequences_ref"],
                    token_readout=model.token_readout,
                )
                batched_alt = model.backbone.encode_batched(
                    batch["sequences_alt"],
                    token_readout=model.token_readout,
                )
                batched_logits = _paired_logits_from_features(model, batched_ref, batched_alt)
                batched_features = torch.cat([batched_ref, batched_alt], dim=-1)
            else:
                batched_features = model.backbone.encode_batched(
                    batch["sequences"],
                    token_readout=model.token_readout,
                )
                batched_logits = _logits_from_features(model, batched_features)
        progress.update()

        embedding_diff = _diff_payload(serial_features, batched_features)
        logit_diff = _diff_payload(serial_logits, batched_logits)
        passed = (
            embedding_diff["max_abs_diff"] <= float(args.embedding_atol)
            and logit_diff["max_abs_diff"] <= float(args.logit_atol)
        )
        progress.update()

        payload = {
            "status": "ok" if passed else "mismatch",
            "model_name": args.model_name,
            "task_name": args.task_name,
            "split_name": args.split_name,
            "batch_size": int(args.batch_size),
            "chunk_forward_batch_size": int(args.chunk_forward_batch_size),
            "token_readout": model.token_readout,
            "finetune_method": resolved_finetune_method,
            "input_schema": task_meta.get("input_schema", "single_sequence"),
            "embedding_atol": float(args.embedding_atol),
            "logit_atol": float(args.logit_atol),
            "embedding_diff": embedding_diff,
            "logit_diff": logit_diff,
            "elapsed_seconds": round(time.perf_counter() - started_at, 6),
            **_max_memory_payload(args.device),
        }
        write_json(output_path, payload)
        progress.update()
    except Exception as exc:
        payload = {
            "status": "failed",
            "model_name": args.model_name,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "traceback": traceback.format_exc(),
            "elapsed_seconds": round(time.perf_counter() - started_at, 6),
            **_max_memory_payload(args.device),
        }
        write_json(output_path, payload)
        progress.close()
        raise

    progress.close()
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if args.fail_on_mismatch and payload["status"] != "ok":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
