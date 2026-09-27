"""File-based dataset partitioning without materializing DNA sequence objects."""

from __future__ import annotations

import hashlib
import random
from collections import Counter
from collections.abc import Callable, Mapping
from contextlib import ExitStack
from pathlib import Path
from typing import Any

from dnakit._dependency_hints import dependency_hint
from dnakit.datasets.config import SplitConfig
from dnakit.datasets.split import _assign_individuals
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.io import RowWriter, iter_rows, write_rows


def export_grouped_split(
    source: str | Path,
    output_dir: str | Path,
    *,
    group_field: str,
    config: SplitConfig | None = None,
    sample_sizes: Mapping[str, int] | None = None,
    sample_filter: Callable[[Mapping[str, Any]], bool] | None = None,
    sample_digest_field: str | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a row file by group and export an assignment CSV and one CSV per split.

    The first occurrence of each group sets its input order. All rows sharing
    that value remain in the same split. Optionally write seeded reservoir
    samples under ``analysis/``; a digest field adds a SHA-256 column named
    ``<field>_sha256`` to those samples. The full row file is streamed twice.
    """
    config = config or SplitConfig()
    if config.method != "random":
        raise ConfigurationError(
            "export_grouped_split requires a random SplitConfig.",
            code="INVALID_TABULAR_SPLIT_METHOD",
        )
    if sample_sizes is not None and any(
        not isinstance(size, int) or isinstance(size, bool) or size < 0
        for size in sample_sizes.values()
    ):
        raise ConfigurationError(
            "sample_sizes must contain non-negative integers.",
            code="INVALID_TABULAR_SAMPLE_SIZE",
        )
    source, output_dir = Path(source), Path(output_dir)
    names = tuple(config.ratios)
    groups = list(dict.fromkeys(row[group_field] for row in iter_rows(source)))
    by_index = _assign_individuals(
        range(len(groups)),
        config.ratios,
        generator=random.Random(config.seed),
        shuffle=config.shuffle,
    )
    assignments = {group: by_index[i] for i, group in enumerate(groups)}
    del by_index
    write_rows(
        ({group_field: group, "split": assignments[group]} for group in groups),
        output_dir / "assignment.csv",
        columns=[group_field, "split"],
        overwrite=overwrite,
    )

    samples: dict[str, list[dict[str, Any]]] = {name: [] for name in names}
    seen = dict.fromkeys(names, 0)
    generator = random.Random(config.seed)
    fields = list(next(iter_rows(source)))
    with ExitStack() as stack:
        writers = {
            name: stack.enter_context(
                RowWriter(output_dir / f"{name}.csv", columns=fields, overwrite=overwrite)
            )
            for name in names
        }
        for row in iter_rows(source):
            name = assignments[row[group_field]]
            writers[name].write(row)
            if (
                sample_sizes is None
                or name not in sample_sizes
                or (sample_filter is not None and not sample_filter(row))
            ):
                continue
            size = sample_sizes.get(name, 0)
            seen[name] += 1
            if len(samples[name]) < size:
                samples[name].append(row)
            else:
                index = generator.randrange(seen[name])
                if index < size:
                    samples[name][index] = row

    if sample_sizes is not None:
        digest_column = f"{sample_digest_field}_sha256" if sample_digest_field else None
        for name in names:
            rows = samples[name]
            if digest_column is not None:
                rows = [
                    dict(
                        row,
                        **{
                            digest_column: hashlib.sha256(
                                str(row[sample_digest_field]).encode()
                            ).hexdigest()
                        },
                    )
                    for row in rows
                ]
            write_rows(
                rows,
                output_dir / "analysis" / f"{name}.csv",
                columns=[*fields, *([digest_column] if digest_column else [])],
                overwrite=overwrite,
            )

    counts = dict(Counter(assignments.values()))
    return {
        "seed": config.seed,
        "group_counts": {name: counts.get(name, 0) for name in names},
        "analysis_counts": {name: len(samples[name]) for name in names},
    }


def export_kfold_table(
    source: str | Path,
    target: str | Path,
    *,
    config: SplitConfig | None = None,
    label_field: str | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """Append a zero-based ``fold`` column and return per-fold label counts.

    The returned list has ``fold``, ``train`` and ``valid`` keys when
    ``label_field`` is provided; otherwise it is empty. Source row order is
    retained. A shuffled 10-fold split with seed 42 is the default.
    """
    config = config or SplitConfig(method="kfold")
    if config.method != "kfold":
        raise ConfigurationError(
            "export_kfold_table requires a kfold SplitConfig.",
            code="INVALID_TABULAR_SPLIT_METHOD",
        )
    rows = list(iter_rows(source))
    assert config.n_splits is not None
    if config.n_splits > len(rows):
        raise ConfigurationError(
            "n_splits cannot exceed the number of input records.",
            code="KFOLD_SIZE_LIMIT",
            context={"n_splits": config.n_splits, "record_count": len(rows)},
        )
    try:
        from sklearn.model_selection import KFold  # type: ignore[import-untyped]
    except ImportError as exc:
        raise BackendUnavailableError(
            "K-fold splitting requires scikit-learn.",
            code="MISSING_KFOLD_DEPENDENCY",
            hint=dependency_hint("bio-analysis"),
        ) from exc

    folds = [0] * len(rows)
    splitter = KFold(
        n_splits=config.n_splits,
        shuffle=config.shuffle,
        random_state=config.seed if config.shuffle else None,
    )
    for fold, (_, valid) in enumerate(splitter.split(rows)):
        for index in valid:
            folds[int(index)] = fold
    write_rows(
        (dict(row, fold=fold) for row, fold in zip(rows, folds, strict=True)),
        target,
        columns=[*rows[0], "fold"],
        overwrite=overwrite,
    )
    if label_field is None:
        return []
    return [
        {
            "fold": fold,
            "train": dict(
                Counter(row[label_field] for row, k in zip(rows, folds, strict=True) if k != fold)
            ),
            "valid": dict(
                Counter(row[label_field] for row, k in zip(rows, folds, strict=True) if k == fold)
            ),
        }
        for fold in range(config.n_splits)
    ]
