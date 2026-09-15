"""Write indexed FASTA records directly into disjoint JSONL output regions."""

from __future__ import annotations

import json
import os
from typing import Any

from dnakit._parallel import ordered_map
from dnakit.core import DNARecord, DNASequence
from dnakit.io._formats import _header, record_to_mapping

_GROUP_BYTES = 4 * 1024**2


def _encode(record: DNARecord) -> bytes:
    return (
        json.dumps(
            record_to_mapping(record), ensure_ascii=False, sort_keys=True, separators=(",", ":")
        )
        + "\n"
    ).encode("utf-8")


def _json_size(layout: tuple[int, bytes]) -> int:
    bases, header = layout
    identifier, description = _header(
        header.decode("utf-8").rstrip("\n"), line_number=1, marker="FASTA"
    )
    probe = DNARecord(
        DNASequence("A" if bases else "", alphabet="iupac"), identifier, description=description
    )
    # JSONL stores symbols in both `sequence` and `parts`. Valid DNA is ASCII.
    size = len(_encode(probe))
    return size + 2 * (bases - 1) if bases else size


def _write_positions(job: Any) -> tuple[Any, int, bool]:
    records, target, positions = job
    complete = 0
    with target.open("r+b") as output:
        for record, (offset, size), (_, start, end) in zip(
            records, positions, records.locations, strict=True
        ):
            data = _encode(record)
            if len(data) != size:
                return target, 0, False  # Irregular layout or Unicode expansion: use the codec.
            output.seek(offset)
            output.write(data)
            complete += end - start
    return target, complete, True


def write_positioned(jobs: list[Any], workers: int) -> bool:
    if any(len(records.layouts) != len(records.locations) for records, *_ in jobs):
        return False
    groups = []
    counters = {path: extra[0] for _, path, _, *extra in jobs if extra}
    completed = {path: 0 for _, path, *_ in jobs}
    for records, target, _config, *_ in jobs:
        positions = []
        offset = 0
        for layout in records.layouts:
            size = _json_size(layout)
            positions.append((offset, size))
            offset += size
        with target.open("wb") as stream:
            stream.truncate(offset)
        first, weight = 0, 0
        for index, (_, start, end) in enumerate(records.locations):
            weight += end - start
            if weight >= _GROUP_BYTES or index == len(records.locations) - 1:
                groups.append(
                    (
                        records.select(list(range(first, index + 1))),
                        target,
                        tuple(positions[first : index + 1]),
                    )
                )
                first, weight = index + 1, 0
    compatible = True
    for path, weight, valid in ordered_map(_write_positions, groups, workers=workers, batch_size=1):
        compatible &= valid
        completed[path] += weight
        if path in counters:
            counters[path].write_text(f"{completed[path]} writing")
    if not compatible:
        for _, target, *_ in jobs:
            target.unlink()
        return False
    for _, target, *_ in jobs:
        if target in counters:
            counters[target].write_text(f"{completed[target]} checking")
        with target.open("r+b") as handle:
            os.fsync(handle.fileno())
        if target in counters:
            counters[target].write_text(f"{completed[target]} done")
    return True
