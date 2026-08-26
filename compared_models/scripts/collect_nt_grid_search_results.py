from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from tqdm.auto import tqdm
except ModuleNotFoundError:  # pragma: no cover
    def tqdm(iterable, **_: Any):
        return iterable


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="汇总 compared_models NT 网格搜索产生的 metrics/profile/runtime JSON。",
    )
    parser.add_argument("--log-root", required=True, help="网格搜索输出根目录。")
    parser.add_argument(
        "--output-csv",
        default=None,
        help="输出 CSV 路径；默认写到 <log-root>/grid_search_summary.csv。",
    )
    parser.add_argument(
        "--output-json",
        default=None,
        help="输出 JSON 路径；默认写到 <log-root>/grid_search_summary.json。",
    )
    include_failed_group = parser.add_mutually_exclusive_group()
    include_failed_group.add_argument(
        "--include-failed",
        dest="include_failed",
        action="store_true",
        default=True,
        help="是否把 failed/skipped 的任务也写入汇总表。",
    )
    include_failed_group.add_argument(
        "--no-include-failed",
        dest="include_failed",
        action="store_false",
        help="只汇总 status=ok 的任务。",
    )
    return parser


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"_read_error": f"{type(exc).__name__}: {exc}"}
    return payload if isinstance(payload, dict) else {"_payload": payload}


def get_nested(payload: dict[str, Any], *keys: str) -> Any:
    current: Any = payload
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def flatten_metric_block(
    row: dict[str, Any],
    metrics: dict[str, Any],
    *,
    split_name: str,
    prefix: str,
) -> None:
    block = get_nested(metrics, "best", split_name)
    if not isinstance(block, dict):
        return
    for key, value in block.items():
        if isinstance(value, (int, float, str, bool)) or value is None:
            row[f"{prefix}_{key}"] = value


def discover_metric_paths(log_root: Path) -> list[Path]:
    return sorted(
        path
        for path in log_root.rglob("metrics.json")
        if "launcher_logs" not in path.parts
    )


def build_row(metrics_path: Path, log_root: Path) -> dict[str, Any]:
    run_dir = metrics_path.parent
    metrics = read_json(metrics_path)
    run_config = read_json(run_dir / "run_config.json")
    profile = read_json(run_dir / "model_profile.json")
    runtime = read_json(run_dir / "runtime.json")
    dataset_meta = read_json(run_dir / "dataset_meta.json")

    row: dict[str, Any] = {
        "task_name": metrics.get("task_name") or run_config.get("task_name"),
        "model_name": metrics.get("model_name") or run_config.get("model_name"),
        "model_display_name": metrics.get("model_display_name") or run_config.get("model_display_name"),
        "status": metrics.get("status", "unknown"),
        "finetune_method": metrics.get("finetune_method") or run_config.get("resolved_finetune_method"),
        "token_readout": metrics.get("token_readout") or run_config.get("token_readout"),
        "lr": run_config.get("lr"),
        "batch_size": run_config.get("batch_size"),
        "weight_decay": run_config.get("weight_decay"),
        "optimizer_name": run_config.get("optimizer_name"),
        "epochs": run_config.get("epochs"),
        "epochs_ran": metrics.get("epochs_ran"),
        "best_epoch": metrics.get("best_epoch"),
        "early_stopping_patience": metrics.get("early_stopping_patience"),
        "stopped_early": metrics.get("stopped_early"),
        "seed": run_config.get("seed"),
        "num_labels": metrics.get("num_labels") or run_config.get("num_labels") or dataset_meta.get("num_labels"),
        "input_schema": metrics.get("input_schema") or run_config.get("input_schema") or dataset_meta.get("input_schema"),
        "chunk_length": metrics.get("chunk_length") or get_nested(run_config, "model_spec", "chunk_length"),
        "chunk_unit": metrics.get("chunk_unit") or get_nested(run_config, "model_spec", "chunk_unit"),
        "total_params": profile.get("total_params"),
        "total_params_m": profile.get("total_params_m"),
        "trainable_params": profile.get("trainable_params"),
        "trainable_params_m": profile.get("trainable_params_m"),
        "gmacs": profile.get("gmacs"),
        "gflops": profile.get("gflops"),
        "activation_memory_mib": profile.get("activation_memory_mib"),
        "forward_backward_activation_mib": profile.get("forward_backward_activation_mib"),
        "profile_status": profile.get("status"),
        "profile_errors": json.dumps(profile.get("errors", {}), ensure_ascii=False),
        "wall_time_seconds": runtime.get("wall_time_seconds"),
        "max_cuda_memory_allocated_mib": runtime.get("max_cuda_memory_allocated_mib"),
        "max_cuda_memory_reserved_mib": runtime.get("max_cuda_memory_reserved_mib"),
        "cuda_device_name": runtime.get("cuda_device_name"),
        "cuda_device_total_memory_mib": runtime.get("cuda_device_total_memory_mib"),
        "metrics_path": str(metrics_path),
        "run_dir": str(run_dir),
        "relative_run_dir": str(run_dir.relative_to(log_root)),
    }
    for split_name in ("train", "valid", "test"):
        flatten_metric_block(row, metrics, split_name=split_name, prefix=f"best_{split_name}")
    if "error_type" in metrics:
        row["error_type"] = metrics.get("error_type")
        row["error_message"] = metrics.get("error_message")
    return row


def main() -> None:
    args = build_parser().parse_args()
    log_root = Path(args.log_root)
    if not log_root.exists():
        raise FileNotFoundError(f"日志根目录不存在: {log_root}")

    metric_paths = discover_metric_paths(log_root)
    rows = [
        build_row(path, log_root)
        for path in tqdm(metric_paths, desc="collect grid results", unit="run")
    ]
    if not args.include_failed:
        rows = [row for row in rows if row.get("status") == "ok"]

    output_csv = Path(args.output_csv) if args.output_csv else log_root / "grid_search_summary.csv"
    output_json = Path(args.output_json) if args.output_json else log_root / "grid_search_summary.json"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    dataframe = pd.DataFrame(rows)
    dataframe.to_csv(output_csv, index=False)
    output_json.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(rows)} rows to {output_csv}")
    print(f"wrote {len(rows)} rows to {output_json}")


if __name__ == "__main__":
    main()
