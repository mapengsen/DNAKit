"""Backward-compatible reverse-complement bridge to shared file processing."""

from __future__ import annotations

from inspect import signature
from pathlib import Path

from dnakit._processing.sequence import run_sequence
from dnakit.core import DNARecord
from dnakit.io import ReadConfig, WriteConfig, WriteResult

from .records import FeatureOverlapPolicy, RecordOperationResult, reverse_complement_record


def reverse_complement(
    value: DNARecord, *, feature_policy: FeatureOverlapPolicy
) -> RecordOperationResult:
    return reverse_complement_record(value, feature_policy=feature_policy)


def reverse_complement_file(
    source: Path,
    target: Path | None,
    *,
    feature_policy: FeatureOverlapPolicy,
    read_config: ReadConfig | None,
    write_config: WriteConfig | None,
    sequence_column: str | None,
    show_progress: bool,
    workers: int = -1,
) -> WriteResult:
    from ._fasta_direction import try_reverse_fasta
    from .records import _feature_policy

    _feature_policy(feature_policy)
    result = try_reverse_fasta(
        source,
        target,
        read_config=read_config,
        write_config=write_config,
        sequence_column=sequence_column,
        workers=workers,
        show_progress=show_progress,
    )
    if result is not None:
        return result
    bound = signature(reverse_complement).bind(source, feature_policy=feature_policy)
    return run_sequence(
        reverse_complement,
        bound,
        "value",
        "record",
        {
            "output_path": target,
            "read_config": read_config,
            "write_config": write_config,
            "sequence_column": sequence_column,
            "show_progress": show_progress,
            "workers": workers,
        },
    )
