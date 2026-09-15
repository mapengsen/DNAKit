"""Spool random/K-fold sequence payloads to disk; retain only IDs and offsets."""

from __future__ import annotations

import inspect
import tempfile
from collections.abc import Iterator
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._parallel import pack, unpack
from dnakit.core import DNARecord, DNASequence
from dnakit.datasets.file_results import FileRecords, FileSplitResult, FileSplitSubset

from .inputs import input_config, iter_inputs, table_config
from .progress import SplitProgress


@dataclass(frozen=True)
class _SpoolRecords:
    path: Path
    locations: tuple[tuple[int, int], ...]

    def __iter__(self) -> Iterator[DNARecord]:
        with self.path.open("rb") as handle:
            for offset, length in self.locations:
                handle.seek(offset)
                yield unpack(handle.read(length))


def random_file_split(
    function: Any,
    bound: Any,
    input_name: str,
    target: Path,
    sources: tuple[Path, ...],
    writer_config: Any,
    read_config: Any,
    column: str | None,
    show: bool,
) -> FileSplitResult:
    from dataclasses import replace

    from .datasets import _export
    from .indexed_records import try_index

    if writer_config.create_parents:
        target.parent.mkdir(parents=True, exist_ok=True)
    with (
        tempfile.TemporaryDirectory(prefix=".dnakit-split-spool-", dir=target.parent) as directory,
        SplitProgress(sources, show) as progress,
    ):
        spool = Path(directory) / "records.bin"
        locations: list[tuple[int, int]] = []
        headers: list[DNARecord] = []
        with progress.reading():
            indexed = (
                try_index(bound.arguments[input_name], sources)
                if read_config is None and column is None
                else None
            )
        if indexed is not None:
            headers, indexed_records = indexed
        else:
            with (
                spool.open("wb") as handle,
                closing(
                    iter_inputs(
                        bound.arguments[input_name], read_config=read_config, sequence_column=column
                    )
                ) as records,
                progress.reading(),
            ):
                for record in records:
                    payload = pack(record)
                    locations.append((handle.tell(), len(payload)))
                    handle.write(payload)
                    headers.append(DNARecord(DNASequence(""), record.id))
                    del payload, record
        progress.read_done()
        bound.arguments[input_name] = headers
        # The existing algorithm supplies exactly the same assignments and order.
        native = inspect.unwrap(function)
        kwargs = {
            key: value
            for key, value in bound.kwargs.items()
            if key in inspect.signature(native).parameters
        }
        plan = native(*bound.args, **kwargs)
        indices: dict[str, list[int]] = {subset.name: [] for subset in plan.subsets}
        for assignment in plan.assignments:
            indices[assignment.split].append(assignment.input_index)
        if not plan.preserve_order and plan.method == "random":
            # Assignment order in the audit is input order; regenerate the RNG ranking.
            import random

            from dnakit.datasets.split import _assign_individuals

            order = _assign_individuals(
                tuple(range(len(headers))),
                plan.ratios,
                generator=random.Random(plan.seed),
                shuffle=plan.shuffle,
            )
            rank = {index: position for position, index in enumerate(order)}
            for values in indices.values():
                values.sort(key=rank.__getitem__)
        staged = SimpleNamespace(
            subsets=tuple(
                SimpleNamespace(
                    name=name,
                    records=(
                        indexed_records.select(values)
                        if indexed is not None
                        else _SpoolRecords(spool, tuple(locations[index] for index in values))
                    ),
                )
                for name, values in indices.items()
            ),
            to_dict=plan.to_dict,
            validate_sources=indexed_records.validate_sources
            if indexed is not None
            else lambda: None,
        )
        _export(
            staged,
            target,
            function.__name__,
            sources,
            writer_config,
            read_config,
            column,
            progress=progress,
        )
    csv_input = any(
        suffix.lower() in {".csv", ".tsv"} for source in sources for suffix in source.suffixes
    )
    subsets = []
    for subset in plan.subsets:
        extension = (".csv" if csv_input else ".jsonl") + (
            ".gz" if writer_config.compression == "gzip" else ""
        )
        path = target / f"{subset.name}{extension}"
        selected = table_config(sources, path, writer_config, read_config, column)
        reader = replace(
            input_config(path, None), compression=selected.compression, encoding=selected.encoding
        )
        if csv_input:
            reader = replace(
                reader,
                sequence_column=selected.sequence_column,
                id_column=selected.id_column,
                description_column=selected.description_column,
                metadata_column=selected.metadata_column,
                features_column=selected.features_column,
                parts_column=selected.parts_column,
                letter_annotations_column=selected.letter_annotations_column,
                delimiter=selected.delimiter,
            )
        subsets.append(FileSplitSubset(subset.name, FileRecords(path, subset.records.ids, reader)))
    return FileSplitResult(tuple(subsets), target, plan)
