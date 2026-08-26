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

from compared_models.benchmark import available_eqtl_tasks, prepare_eqtl_task
from visualdna.utils.seed_utils import set_seed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="把 eQTL 原始数据预处理成统一 benchmark 双序列格式。",
    )
    parser.add_argument("--task-name", required=True, help="组织名，例如 Adipose_Subcutaneous；或使用 all。")
    parser.add_argument("--raw-download", required=True, help="原始 eQTL 根目录，例如 /zengxiangxiang/mps/visualdna/data/raw_download/eqtl。")
    parser.add_argument("--processed-download", required=True, help="统一 benchmark 输出根目录。")
    parser.add_argument(
        "--max-samples-per-split",
        type=int,
        default=None,
        help="可选：每个 split 最多保留前 N 个样本，便于冒烟测试。",
    )
    parser.add_argument(
        "--sequence-length-cutoff-override",
        type=int,
        default=None,
        help="可选：覆盖原始 seq_len_cutoff，生成更短的 smoke 版 eQTL 序列。",
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
    task_names = available_eqtl_tasks(raw_root)
    if not task_names:
        raise RuntimeError(f"在 {raw_root} 下未找到任何可处理的 eQTL 任务。")
    return task_names


def main() -> None:
    args = build_parser().parse_args()
    set_seed(args.seed)
    raw_root = Path(args.raw_download)
    processed_root = Path(args.processed_download)
    summaries = []
    for task_name in resolve_task_names(args.task_name, raw_root):
        meta = prepare_eqtl_task(
            task_name=task_name,
            raw_root=raw_root,
            output_root=processed_root,
            overwrite=args.overwrite,
            max_samples_per_split=args.max_samples_per_split,
            sequence_length_cutoff_override=args.sequence_length_cutoff_override,
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
