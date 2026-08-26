from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import tabix

from visualdna.data.initializer.helpers.dnalongbench.utils import (
    FastaStringExtractor,
    Interval,
    parse_config,
    rcDNA,
)


def available_eqtl_tasks(raw_root: str | Path) -> list[str]:
    config_root = Path(raw_root) / "eQTL" / "config"
    if not config_root.exists():
        raise FileNotFoundError(f"eQTL config 目录不存在: {config_root}")
    task_names = []
    for config_path in sorted(config_root.glob("gtex_hg38.*.config")):
        stem = config_path.stem
        if stem == "gtex_hg38":
            continue
        task_names.append(stem.replace("gtex_hg38.", ""))
    return task_names


def _parse_eqtl_split(
    *,
    eqtl_root: Path,
    config: dict[str, Any],
    split_name: str,
    max_records: int | None,
    sequence_length_cutoff_override: int | None,
) -> pd.DataFrame:
    fasta_reader = FastaStringExtractor(str(eqtl_root / config["genome_fa"]))
    blacklist_tabix = tabix.open(str(eqtl_root / config["eQTL_tabix_file"]))
    df = pd.read_csv(eqtl_root / config["eQTL_file"], sep="\t", header=0)
    seq_len_cutoff = int(
        sequence_length_cutoff_override
        if sequence_length_cutoff_override is not None
        else config.get("seq_len_cutoff", 450000)
    )
    tss_flank_upstream = int(config.get("tss_flank_upstream", 3000))
    tss_flank_downstream = int(config.get("tss_flank_downstream", 3000))
    region_flank_upstream = int(config.get("region_flank_upstream", 500))
    region_flank_downstream = int(config.get("region_flank_downstream", 500))

    records: list[dict[str, Any]] = []
    try:
        for _, row in df.iterrows():
            if row["subset"] != split_name:
                continue
            if row["gene_chrom"] != row["region_chrom"]:
                continue

            if row["gene_strand"] == "+":
                tss_start = int(row["gene_start"])
            elif row["gene_strand"] == "-":
                tss_start = int(row["gene_end"]) - 1
            else:
                continue

            tss_end = tss_start + 1
            tss_start -= tss_flank_upstream
            tss_end += tss_flank_downstream

            variant_start = int(row["region_start"])
            variant_end = int(row["region_end"])
            region_start = variant_start - region_flank_upstream
            region_end = variant_end + region_flank_downstream
            distance = max(0, max(tss_start, region_start) - min(tss_end, region_end))
            if distance > seq_len_cutoff:
                continue

            sequence_start = min(tss_start, region_start)
            sequence_end = max(tss_end, region_end)
            region_interval = Interval(row["region_chrom"], sequence_start, sequence_end)
            region_seq = fasta_reader.extract(region_interval)

            variant_rel_start = variant_start - sequence_start
            variant_rel_end = variant_end - sequence_start
            if region_seq[variant_rel_start:variant_rel_end] != row["allele1"]:
                continue

            if distance > 0:
                if tss_start > region_end:
                    query_start, query_end = region_end, tss_start
                else:
                    query_start, query_end = tss_end, region_start
                if query_start < query_end:
                    for overlap in blacklist_tabix.query(row["region_chrom"], query_start, query_end):
                        overlap_start = int(overlap[1]) - sequence_start
                        overlap_end = int(overlap[2]) - sequence_start
                        if 0 < overlap_start < len(region_seq) and 0 < overlap_end < len(region_seq):
                            region_seq = (
                                region_seq[:overlap_start]
                                + "N" * (overlap_end - overlap_start)
                                + region_seq[overlap_end:]
                            )

            region_seq_alt = region_seq[:variant_rel_start] + row["allele2"] + region_seq[variant_rel_end:]
            if int(row["gene_start"]) > region_end:
                region_seq = rcDNA(region_seq)
                region_seq_alt = rcDNA(region_seq_alt)

            if len(region_seq) <= seq_len_cutoff:
                region_seq = region_seq + "N" * (seq_len_cutoff - len(region_seq))
            else:
                region_seq = region_seq[:seq_len_cutoff]
            if len(region_seq_alt) <= seq_len_cutoff:
                region_seq_alt = region_seq_alt + "N" * (seq_len_cutoff - len(region_seq_alt))
            else:
                region_seq_alt = region_seq_alt[:seq_len_cutoff]

            label_original = str(row["target"])
            records.append(
                {
                    "sequence_ref": region_seq,
                    "sequence_alt": region_seq_alt,
                    "raw_length_ref": len(region_seq),
                    "raw_length_alt": len(region_seq_alt),
                    "label_original": label_original,
                    "label": 1 if label_original == "positive" else 0,
                    "name": str(row["gene_id"]),
                }
            )
            if max_records is not None and len(records) >= int(max_records):
                break
    finally:
        fasta_reader.close()

    if not records:
        raise RuntimeError(f"未从 eQTL split={split_name} 解析出任何样本。")
    return pd.DataFrame.from_records(records)


def _canonicalize_split(
    df: pd.DataFrame,
    *,
    task_name: str,
    split_name: str,
    num_labels: int,
) -> pd.DataFrame:
    canonical = df.copy().reset_index(drop=True)
    canonical.insert(0, "sample_id", [f"{task_name}-{split_name}-{idx:08d}" for idx in range(len(canonical))])
    canonical["split"] = split_name
    canonical["task_name"] = task_name
    canonical["task_type"] = "classification"
    canonical["input_schema"] = "pair_sequence"
    canonical["num_labels"] = int(num_labels)
    keep_columns = [
        "sample_id",
        "task_name",
        "task_type",
        "input_schema",
        "split",
        "sequence_ref",
        "sequence_alt",
        "raw_length_ref",
        "raw_length_alt",
        "label",
        "label_original",
        "num_labels",
        "name",
    ]
    return canonical.loc[:, keep_columns]


def _length_summary(df: pd.DataFrame) -> dict[str, float | int]:
    if df.empty:
        return {
            "count": 0,
            "min_ref": 0,
            "max_ref": 0,
            "mean_ref": 0.0,
            "min_alt": 0,
            "max_alt": 0,
            "mean_alt": 0.0,
        }
    ref_lengths = df["raw_length_ref"].astype(int)
    alt_lengths = df["raw_length_alt"].astype(int)
    return {
        "count": int(len(df)),
        "min_ref": int(ref_lengths.min()),
        "max_ref": int(ref_lengths.max()),
        "mean_ref": float(ref_lengths.mean()),
        "min_alt": int(alt_lengths.min()),
        "max_alt": int(alt_lengths.max()),
        "mean_alt": float(alt_lengths.mean()),
    }


def prepare_eqtl_task(
    *,
    task_name: str,
    raw_root: str | Path,
    output_root: str | Path,
    overwrite: bool = False,
    max_samples_per_split: int | None = None,
    sequence_length_cutoff_override: int | None = None,
) -> dict[str, Any]:
    raw_root_path = Path(raw_root)
    output_root_path = Path(output_root)
    eqtl_root = raw_root_path / "eQTL"
    task_output_dir = output_root_path / task_name
    meta_path = task_output_dir / "meta.json"
    if meta_path.exists() and not overwrite:
        return json.loads(meta_path.read_text(encoding="utf-8"))

    config_path = eqtl_root / "config" / f"gtex_hg38.{task_name}.config"
    if not config_path.exists():
        raise FileNotFoundError(f"缺少 eQTL config: {config_path}")
    config = parse_config(str(config_path))

    train_df = _parse_eqtl_split(
        eqtl_root=eqtl_root,
        config=config,
        split_name="train",
        max_records=max_samples_per_split,
        sequence_length_cutoff_override=sequence_length_cutoff_override,
    )
    valid_df = _parse_eqtl_split(
        eqtl_root=eqtl_root,
        config=config,
        split_name="valid",
        max_records=max_samples_per_split,
        sequence_length_cutoff_override=sequence_length_cutoff_override,
    )
    test_df = _parse_eqtl_split(
        eqtl_root=eqtl_root,
        config=config,
        split_name="test",
        max_records=max_samples_per_split,
        sequence_length_cutoff_override=sequence_length_cutoff_override,
    )

    num_labels = 2
    train_canonical = _canonicalize_split(train_df, task_name=task_name, split_name="train", num_labels=num_labels)
    valid_canonical = _canonicalize_split(valid_df, task_name=task_name, split_name="valid", num_labels=num_labels)
    test_canonical = _canonicalize_split(test_df, task_name=task_name, split_name="test", num_labels=num_labels)

    meta = {
        "task_name": task_name,
        "task_family": "eqtl",
        "task_type": "classification",
        "input_schema": "pair_sequence",
        "num_labels": 2,
        "label_mapping": {"negative": 0, "positive": 1},
        "inverse_label_mapping": {"0": "negative", "1": "positive"},
        "raw_root": str(raw_root_path),
        "source_config": str(config_path),
        "max_samples_per_split": None if max_samples_per_split is None else int(max_samples_per_split),
        "sequence_length_cutoff": int(
            sequence_length_cutoff_override
            if sequence_length_cutoff_override is not None
            else config.get("seq_len_cutoff", 450000)
        ),
        "splits": {
            "train": _length_summary(train_canonical),
            "valid": _length_summary(valid_canonical),
            "test": _length_summary(test_canonical),
        },
    }

    task_output_dir.mkdir(parents=True, exist_ok=True)
    train_canonical.to_parquet(task_output_dir / "train.parquet", index=False)
    valid_canonical.to_parquet(task_output_dir / "valid.parquet", index=False)
    test_canonical.to_parquet(task_output_dir / "test.parquet", index=False)
    train_canonical.to_csv(task_output_dir / "train.csv", index=False)
    valid_canonical.to_csv(task_output_dir / "valid.csv", index=False)
    test_canonical.to_csv(task_output_dir / "test.csv", index=False)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return meta
