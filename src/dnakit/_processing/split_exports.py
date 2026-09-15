"""Balance indexed split exports across record groups, including the large train subset."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

from dnakit._parallel import available_memory, current_workers, ordered_map

from .fragments import write_fragment
from .indexed_records import IndexedRecords

_SHARD_BYTES = 8 * 1024**2
_MIN_PARALLEL_BYTES = 16 * 1024**2


def _write_group(job: Any) -> tuple[Path, Path, int]:
    records, fragment, config, target = job
    write_fragment(records, fragment, config)
    return target, fragment, sum(end - start for _, start, end in records.locations)


def write_indexed_outputs(jobs: list[Any]) -> bool:
    if not jobs or any(
        not isinstance(records, IndexedRecords)
        or path.suffix != ".jsonl"
        or config.compression == "gzip"
        or config.max_output_bytes is not None
        or config.encoding.lower().replace("_", "-") not in {"utf-8", "utf8"}
        for records, path, config, *_ in jobs
    ):
        return False
    largest = max(
        (end - start for records, *_ in jobs for _, start, end in records.locations), default=0
    )
    count = min(
        current_workers(), 8, max(1, available_memory() // max(512 * 1024**2, largest * 12))
    )
    total = sum(end - start for records, *_ in jobs for _, start, end in records.locations)
    if total < _MIN_PARALLEL_BYTES:
        count = 1
    if count == 1:
        for records, target, config, *counters in jobs:
            if counters:
                from .progress import tracked_records

                stream = tracked_records(records, counters[0])
            else:
                stream = records
            write_fragment(stream, target, config)
            with target.open("r+b") as handle:
                os.fsync(handle.fileno())
            if counters:
                total = sum(end - start for _, start, end in records.locations)
                counters[0].write_text(f"{total} done")
        return True
    from .positioned_split import write_positioned

    if write_positioned(jobs, count):
        return True
    groups: list[Any] = []
    fragments: dict[Path, list[Path]] = {}
    totals: dict[Path, int] = {}
    counters_by_target = {path: counters[0] for _, path, _, *counters in jobs if counters}

    def submit(records: IndexedRecords, target: Path, config: Any, locations: Any) -> None:
        fragment = target.with_name(f".{target.stem}.{len(groups)}.jsonl")
        groups.append((IndexedRecords(tuple(locations), records.stamps), fragment, config, target))

    for records, target, config, *_ in jobs:
        fragments[target] = []
        totals[target] = 0
        locations, size = [], 0

        for location in records.locations:
            locations.append(location)
            size += location[2] - location[1]
            if size >= _SHARD_BYTES:
                submit(records, target, config, locations)
                locations, size = [], 0
        if locations:
            submit(records, target, config, locations)
    for target, fragment, size in ordered_map(_write_group, groups, workers=count, batch_size=1):
        fragments[target].append(fragment)
        totals[target] += size
        if target in counters_by_target:
            counters_by_target[target].write_text(f"{totals[target]} writing")
    for target, pieces in fragments.items():
        if target in counters_by_target:
            counters_by_target[target].write_text(f"{totals[target]} checking")
        with target.open("wb") as handle:
            for fragment in pieces:
                with fragment.open("rb") as source:
                    shutil.copyfileobj(source, handle, length=1024**2)
                fragment.unlink()
            handle.flush()
            os.fsync(handle.fileno())
        if target in counters_by_target:
            counters_by_target[target].write_text(f"{totals[target]} done")
    return True
