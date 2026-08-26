from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VISUALDNA_ROOT = PROJECT_ROOT / "visualdna"
for path in (PROJECT_ROOT, VISUALDNA_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from compared_models.benchmark import available_genomic_benchmark_tasks, prepare_genomic_benchmark_task
from visualdna.utils.seed_utils import set_seed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="把 genomic_benchmarks 原始文本序列预处理成统一 benchmark 格式。",
    )
    parser.add_argument("--task-name", required=True, help="任务名，例如 human_nontata_promoters；或使用 all。")
    parser.add_argument("--raw-download", required=True, help="原始 genomic_benchmarks 根目录。")
    parser.add_argument("--processed-download", required=True, help="统一 benchmark 输出根目录。")
    parser.add_argument("--valid-ratio", type=float, default=0.1, help="从 train 划出 valid 的比例。")
    parser.add_argument("--split-seed", type=int, default=42, help="train/valid 切分随机种子。")
    parser.add_argument(
        "--stratified-split",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否按标签分层切分 train/valid。",
    )
    parser.add_argument(
        "--max-samples-per-split",
        type=int,
        default=None,
        help="可选：每个 split 只保留前 N 个样本，便于冒烟测试。",
    )
    parser.add_argument(
        "--overwrite",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="若目标目录已存在，是否覆盖写入。",
    )
    parser.add_argument("--seed", type=int, default=42, help="运行随机种子。")
    return parser


def resolve_task_names(task_name: str, raw_root: str | Path) -> list[str]:
    if task_name != "all":
        return [task_name]
    task_names = available_genomic_benchmark_tasks(raw_root)
    if not task_names:
        raise RuntimeError(f"在 {raw_root} 下未找到任何可处理的 genomic_benchmarks 任务。")
    return task_names


def main() -> None:
    args = build_parser().parse_args()
    set_seed(args.seed)
    raw_root = Path(args.raw_download)
    processed_root = Path(args.processed_download)
    summaries = []
    for task_name in resolve_task_names(args.task_name, raw_root):
        meta = prepare_genomic_benchmark_task(
            task_name=task_name,
            raw_root=raw_root,
            output_root=processed_root,
            valid_ratio=args.valid_ratio,
            split_seed=args.split_seed,
            stratified=args.stratified_split,
            overwrite=args.overwrite,
            max_samples_per_split=args.max_samples_per_split,
        )
        summaries.append(meta)
        print(
            json.dumps(
                {
                    "task_name": task_name,
                    "output_dir": str(processed_root / task_name),
                    "num_labels": meta["num_labels"],
                    "splits": meta["splits"],
                },
                indent=2,
                ensure_ascii=False,
            )
        )

    if len(summaries) > 1:
        processed_root.mkdir(parents=True, exist_ok=True)
        (processed_root / "all_tasks_meta.json").write_text(
            json.dumps({"mode": "all", "tasks": summaries}, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
