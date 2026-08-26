from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import torch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import (
    SequenceBenchmarkDataset,
    available_compared_models,
    build_linear_probe_model,
    build_sequence_dataloader,
    get_model_spec,
    load_processed_task_dataframe,
    load_processed_task_meta,
    validate_finetune_method,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="对 compared_models 执行单批次前向 smoke test，验证给定数据集是否能跑通。",
    )
    parser.add_argument("--task-name", required=True, help="处理后的 benchmark 任务名。")
    parser.add_argument(
        "--processed-download",
        required=True,
        help="统一 benchmark 数据根目录，例如 /zengxiangxiang/mps/ood_imageDNA/data/low_similarity_sequence_csv_original_parquet/nt。",
    )
    parser.add_argument("--log-root", required=True, help="验证结果输出目录。")
    parser.add_argument("--model-name", required=True, help="模型名，或使用 all。")
    parser.add_argument(
        "--split-name",
        default="train",
        choices=["train", "valid", "test"],
        help="从哪个 split 抽取一个 batch 做前向。",
    )
    parser.add_argument("--batch-size", type=int, default=1, help="smoke batch size。")
    parser.add_argument("--num-workers", type=int, default=0, help="DataLoader worker 数。")
    parser.add_argument(
        "--max-samples",
        type=int,
        default=2,
        help="从指定 split 最多取多少个样本构建 smoke dataset。",
    )
    parser.add_argument(
        "--token-readout",
        default="auto",
        choices=["auto", "mean", "last"],
        help="token 级读出方式。",
    )
    parser.add_argument(
        "--finetune-method",
        default="frozen_linear_probe",
        help="下游微调方式：frozen_linear_probe / full / ia3；具体以当前模型允许方式为准。",
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
        default=1,
        help="单个样本内部 chunk 微批大小。",
    )
    parser.add_argument(
        "--allow-remote-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="本地 checkpoint 不完整时，是否回退到 Hugging Face。",
    )
    parser.add_argument("--checkpoint-override", default=None, help="覆盖默认 checkpoint。")
    parser.add_argument("--device", default="cpu", help="运行设备，默认 cpu。")
    parser.add_argument(
        "--continue-on-error",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="某个模型失败时是否继续剩余模型。",
    )
    parser.add_argument("--cache-dir", default=None, help="Hugging Face 缓存目录。")
    return parser


def resolve_model_names(model_name: str) -> list[str]:
    if model_name != "all":
        return [model_name]
    return available_compared_models(include_unsupported=True)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    args = build_parser().parse_args()
    task_meta = load_processed_task_meta(args.processed_download, args.task_name)
    pair_mode = task_meta.get("input_schema") == "pair_sequence"

    dataframe = load_processed_task_dataframe(args.processed_download, args.task_name, args.split_name)
    if args.max_samples is not None and args.max_samples > 0:
        dataframe = dataframe.iloc[: int(args.max_samples)].reset_index(drop=True)
    dataset = SequenceBenchmarkDataset(dataframe)
    loader = build_sequence_dataloader(
        dataset,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        shuffle=False,
        is_distributed=False,
        pin_memory=False,
    )
    batch = next(iter(loader))
    log_root = Path(args.log_root)
    results: list[dict] = []

    for model_name in resolve_model_names(args.model_name):
        spec = get_model_spec(model_name)
        model_log_dir = log_root / model_name
        if spec.loader_kind == "unsupported":
            payload = {
                "task_name": args.task_name,
                "split_name": args.split_name,
                "model_name": model_name,
                "status": "skipped",
                "reason": spec.unsupported_reason or "当前脚本未接入该模型。",
            }
            write_json(model_log_dir / "metrics.json", payload)
            results.append(payload)
            continue

        started_at = time.perf_counter()
        try:
            resolved_finetune_method = validate_finetune_method(spec, args.finetune_method)
            model = build_linear_probe_model(
                model_name=model_name,
                num_labels=int(task_meta["num_labels"]),
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
            model = model.to(args.device)
            model.eval()
            with torch.no_grad():
                if pair_mode:
                    logits = model(batch["sequences_ref"], batch["sequences_alt"])
                else:
                    logits = model(batch["sequences"])
            payload = {
                "task_name": args.task_name,
                "split_name": args.split_name,
                "model_name": model_name,
                "status": "ok",
                "logits_shape": list(logits.shape),
                "num_labels": int(task_meta["num_labels"]),
                "input_schema": task_meta.get("input_schema", "single_sequence"),
                "finetune_method": resolved_finetune_method,
                "elapsed_seconds": round(time.perf_counter() - started_at, 6),
            }
        except Exception as exc:
            payload = {
                "task_name": args.task_name,
                "split_name": args.split_name,
                "model_name": model_name,
                "status": "failed",
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
                "elapsed_seconds": round(time.perf_counter() - started_at, 6),
            }
            if not args.continue_on_error:
                write_json(model_log_dir / "metrics.json", payload)
                raise

        write_json(model_log_dir / "metrics.json", payload)
        results.append(payload)

    write_json(
        log_root / "all_results.json",
        {
            "task_name": args.task_name,
            "split_name": args.split_name,
            "input_schema": task_meta.get("input_schema", "single_sequence"),
            "results": results,
        },
    )


if __name__ == "__main__":
    main()
