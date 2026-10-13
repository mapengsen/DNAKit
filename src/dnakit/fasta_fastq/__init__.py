"""Dependency-free FASTA/FASTQ candidate implementations for local staging."""

from __future__ import annotations

import random
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence, Tuple, Union

from dnakit.sequence import (
    normalize_dna_sequence,
    reverse_complement_dna,
    split_sequence_windows,
)


@dataclass(frozen=True)
class FastqRecord:
    """One normalized FASTQ record."""

    identifier: str
    sequence: str
    quality: str

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.quality):
            raise ValueError("FASTQ sequence and quality lengths must match")


Record = Tuple[str, str]
FastqInput = Union[FastqRecord, Tuple[str, str, str]]


def _read_source(source: str | Path) -> str:
    if isinstance(source, Path) or (isinstance(source, str) and "\n" not in source and Path(source).exists()):
        return Path(source).read_text(encoding="utf-8")
    if not isinstance(source, str):
        raise TypeError("source must be text or a filesystem path")
    return source


def _fasta_records(source: str | Path | Sequence[Record]) -> list[Record]:
    if isinstance(source, Sequence) and not isinstance(source, (str, bytes, Path)):
        return [(str(identifier), normalize_dna_sequence(sequence)) for identifier, sequence in source]
    text = _read_source(source)
    records: list[Record] = []
    identifier: str | None = None
    sequence_parts: list[str] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, normalize_dna_sequence("".join(sequence_parts))))
            identifier = line[1:].strip().split(None, 1)[0]
            if not identifier:
                raise ValueError(f"FASTA header at line {line_number} is empty")
            sequence_parts = []
        elif identifier is None:
            raise ValueError("FASTA sequence appears before the first header")
        else:
            sequence_parts.append(line)
    if identifier is not None:
        records.append((identifier, normalize_dna_sequence("".join(sequence_parts))))
    return records


def _fastq_records(source: str | Path | Sequence[FastqInput]) -> list[FastqRecord]:
    if isinstance(source, Sequence) and not isinstance(source, (str, bytes, Path)):
        result = []
        for item in source:
            result.append(item if isinstance(item, FastqRecord) else FastqRecord(str(item[0]), normalize_dna_sequence(item[1]), str(item[2])))
        return result
    lines = [line.rstrip("\r") for line in _read_source(source).splitlines() if line.strip()]
    if len(lines) % 4:
        raise ValueError("FASTQ input must contain complete four-line records")
    result = []
    for offset in range(0, len(lines), 4):
        header, sequence, plus, quality = lines[offset : offset + 4]
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"invalid FASTQ record at line {offset + 1}")
        result.append(FastqRecord(header[1:].split(None, 1)[0], normalize_dna_sequence(sequence), quality))
    return result


def _as_fastq(records: Sequence[FastqInput]) -> list[FastqRecord]:
    return _fastq_records(records)


def _write_fasta(records: Sequence[Record], *, width: int = 80) -> str:
    if width < 1:
        raise ValueError("width must be positive")
    lines = []
    for identifier, sequence in records:
        lines.append(f">{identifier}")
        lines.extend(sequence[index : index + width] for index in range(0, len(sequence), width))
    return "\n".join(lines) + ("\n" if lines else "")


def _write_fastq(records: Sequence[FastqRecord]) -> str:
    lines = []
    for record in records:
        lines.extend((f"@{record.identifier}", record.sequence, "+", record.quality))
    return "\n".join(lines) + ("\n" if lines else "")


def build_fasta_index(source: str | Path | Sequence[Record]) -> dict[str, dict[str, int]]:
    """Build an in-memory FASTA index with record order and lengths."""

    records = _fasta_records(source)
    return {identifier: {"record_index": index, "length": len(sequence)} for index, (identifier, sequence) in enumerate(records)}


def build_fastq_index(source: str | Path | Sequence[FastqInput]) -> dict[str, dict[str, int]]:
    """Build an in-memory FASTQ index with record order and lengths."""

    records = _fastq_records(source)
    return {record.identifier: {"record_index": index, "length": len(record.sequence), "quality_length": len(record.quality)} for index, record in enumerate(records)}


def calculate_gc_profile(source: str | Path | Sequence[Record]) -> list[dict[str, object]]:
    """Return per-record GC fractions."""

    result = []
    for identifier, sequence in _fasta_records(source):
        gc = sum(base in "GC" for base in sequence) / len(sequence) if sequence else 0.0
        result.append({"id": identifier, "length": len(sequence), "gc_fraction": gc, "gc_percent": gc * 100})
    return result


def calculate_length_distribution(source: str | Path | Sequence[Record] | Sequence[FastqInput]) -> dict[str, object]:
    """Summarize record lengths for FASTA or FASTQ."""

    if isinstance(source, Sequence) and source and isinstance(source[0], FastqRecord):
        lengths = [len(record.sequence) for record in _fastq_records(source)]
    else:
        try:
            records = _fasta_records(source)  # type: ignore[arg-type]
            lengths = [len(sequence) for _, sequence in records]
        except (ValueError, TypeError):
            lengths = [len(record.sequence) for record in _fastq_records(source)]  # type: ignore[arg-type]
    return _distribution(lengths)


def _distribution(values: Sequence[int]) -> dict[str, object]:
    if not values:
        return {"count": 0, "min": None, "max": None, "mean": None, "counts": {}}
    counts = Counter(values)
    return {"count": len(values), "min": min(values), "max": max(values), "mean": sum(values) / len(values), "counts": dict(sorted(counts.items()))}


def calculate_n_content_profile(source: str | Path | Sequence[Record]) -> list[dict[str, object]]:
    """Return N counts and fractions for each FASTA record."""

    return [{"id": identifier, "n_count": sequence.count("N"), "n_fraction": sequence.count("N") / len(sequence) if sequence else 0.0} for identifier, sequence in _fasta_records(source)]


def calculate_quality_distribution(source: str | Path | Sequence[FastqInput]) -> dict[str, object]:
    """Summarize Phred+33 quality characters."""

    qualities = [ord(character) - 33 for record in _fastq_records(source) for character in record.quality]
    return _distribution(qualities)


def compare_fasta_files(left: str | Path | Sequence[Record], right: str | Path | Sequence[Record]) -> dict[str, object]:
    """Compare FASTA IDs and normalized sequences."""

    left_map = dict(_fasta_records(left))
    right_map = dict(_fasta_records(right))
    common = sorted(set(left_map) & set(right_map))
    differences = [identifier for identifier in common if left_map[identifier] != right_map[identifier]]
    return {"left_only": sorted(set(left_map) - set(right_map)), "right_only": sorted(set(right_map) - set(left_map)), "common": common, "different_sequences": differences, "identical": not differences and set(left_map) == set(right_map)}


def convert_from_fasta(source: str | Path | Sequence[Record], *, quality: str = "I") -> list[FastqRecord]:
    """Convert FASTA records to FASTQ with a constant quality character."""

    if len(quality) != 1 or ord(quality) < 33:
        raise ValueError("quality must be one printable Phred+33 character")
    return [FastqRecord(identifier, sequence, quality * len(sequence)) for identifier, sequence in _fasta_records(source)]


def convert_from_tabular(source: str | Path | Sequence[Mapping[str, object]]) -> list[Record]:
    """Convert rows with id/name and sequence columns to FASTA records."""

    if isinstance(source, Sequence) and not isinstance(source, (str, bytes, Path)):
        return [(str(row.get("id", row.get("name"))), normalize_dna_sequence(str(row["sequence"]))) for row in source]
    rows = []
    for line in _read_source(source).splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 2:
            raise ValueError("tabular sequence input requires id and sequence columns")
        rows.append((parts[0], normalize_dna_sequence(parts[1])))
    return rows


def convert_to_fasta(source: str | Path | Sequence[FastqInput]) -> list[Record]:
    """Drop quality strings from FASTQ records."""

    return [(record.identifier, record.sequence) for record in _fastq_records(source)]


def convert_to_tabular(source: str | Path | Sequence[Record]) -> list[dict[str, object]]:
    """Convert FASTA records to serializable rows."""

    return [{"id": identifier, "sequence": sequence} for identifier, sequence in _fasta_records(source)]


def count_reads(source: str | Path | Sequence[FastqInput]) -> int:
    """Count FASTQ records."""

    return len(_fastq_records(source))


def count_records(source: str | Path | Sequence[Record]) -> int:
    """Count FASTA records."""

    return len(_fasta_records(source))


def deinterleave_pairs(source: str | Path | Sequence[FastqInput]) -> tuple[list[FastqRecord], list[FastqRecord]]:
    """Split alternating R1/R2 records into two lists."""

    records = _fastq_records(source)
    return records[::2], records[1::2]


def estimate_duplication(source: str | Path | Sequence[FastqInput]) -> dict[str, object]:
    """Estimate exact sequence duplication in FASTQ records."""

    records = _fastq_records(source)
    counts = Counter(record.sequence for record in records)
    duplicate_reads = sum(count - 1 for count in counts.values() if count > 1)
    return {"records": len(records), "unique_sequences": len(counts), "duplicate_records": duplicate_reads, "duplicate_fraction": duplicate_reads / len(records) if records else 0.0}


def estimate_error_rate(source: str | Path | Sequence[FastqInput], *, expected_quality: int = 30) -> float:
    """Estimate the fraction of bases below an expected Phred score."""

    records = _fastq_records(source)
    values = [ord(character) - 33 for record in records for character in record.quality]
    return sum(value < expected_quality for value in values) / len(values) if values else 0.0


def _by_id(source: str | Path | Sequence[Record]) -> dict[str, str]:
    return dict(_fasta_records(source))


def _extract(source: str | Path | Sequence[Record], identifier: str, start: int, end: int) -> Record:
    records = _by_id(source)
    if identifier not in records:
        raise KeyError(identifier)
    sequence = records[identifier]
    if start < 0 or end < start or end > len(sequence):
        raise ValueError("region is outside the record")
    return f"{identifier}:{start}-{end}", sequence[start:end]


def extract_cds(source: str | Path | Sequence[Record], *, identifier: str, intervals: Sequence[Sequence[int]]) -> Record:
    sequence = _by_id(source)[identifier]
    pieces = [sequence[start:end] for start, end in intervals]
    return f"{identifier}:cds", "".join(pieces)


def extract_exon(source: str | Path | Sequence[Record], *, identifier: str, start: int, end: int) -> Record:
    return _extract(source, identifier, start, end)


def extract_flanking(source: str | Path | Sequence[Record], *, identifier: str, start: int, end: int, flank: int = 100) -> Record:
    sequence = _by_id(source)[identifier]
    return f"{identifier}:flank", sequence[max(0, start - flank) : min(len(sequence), end + flank)]


def extract_introns(source: str | Path | Sequence[Record], *, identifier: str, exons: Sequence[Sequence[int]]) -> list[Record]:
    sequence = _by_id(source)[identifier]
    ordered = sorted(exons)
    boundaries = [ordered[index][1] for index in range(len(ordered) - 1)]
    return [(f"{identifier}:intron:{index}", sequence[start:end]) for index, (start, end) in enumerate(zip(boundaries, [item[0] for item in ordered[1:]]), start=1) if end >= start]


def extract_promoter(source: str | Path | Sequence[Record], *, identifier: str, tss: int, upstream: int = 1000, downstream: int = 100) -> Record:
    sequence = _by_id(source)[identifier]
    return f"{identifier}:promoter", sequence[max(0, tss - upstream) : min(len(sequence), tss + downstream)]


def extract_read_names(source: str | Path | Sequence[FastqInput]) -> list[str]:
    return [record.identifier for record in _fastq_records(source)]


def extract_region(source: str | Path | Sequence[Record], *, identifier: str, start: int, end: int) -> Record:
    return _extract(source, identifier, start, end)


def extract_regions_from_bed(source: str | Path | Sequence[Record], bed: str | Path | Sequence[Sequence[object]]) -> list[Record]:
    records = _by_id(source)
    if isinstance(bed, Sequence) and not isinstance(bed, (str, bytes, Path)):
        rows = bed
    else:
        rows = [line.split("\t") for line in _read_source(bed).splitlines() if line.strip() and not line.startswith("#")]
    return [
        (str(row[3]) if len(row) > 3 else f"{row[0]}:{row[1]}-{row[2]}", records[str(row[0])][int(row[1]) : int(row[2])])
        for row in rows
        if str(row[0]) in records
    ]


def filter(source: str | Path | Sequence[Record], *, min_length: int | None = None, max_length: int | None = None, ids: set[str] | None = None) -> list[Record]:
    return [(identifier, sequence) for identifier, sequence in _fasta_records(source) if (min_length is None or len(sequence) >= min_length) and (max_length is None or len(sequence) <= max_length) and (ids is None or identifier in ids)]


def filter_by_length(source: str | Path | Sequence[FastqInput], *, min_length: int | None = None, max_length: int | None = None) -> list[FastqRecord]:
    return [record for record in _fastq_records(source) if (min_length is None or len(record.sequence) >= min_length) and (max_length is None or len(record.sequence) <= max_length)]


def filter_by_quality(source: str | Path | Sequence[FastqInput], *, min_quality: float) -> list[FastqRecord]:
    return [record for record in _fastq_records(source) if (sum(ord(character) - 33 for character in record.quality) / len(record.quality) if record.quality else 0.0) >= min_quality]


def get_ids(source: str | Path | Sequence[Record]) -> list[str]:
    return [identifier for identifier, _ in _fasta_records(source)]


def get_lengths(source: str | Path | Sequence[Record]) -> dict[str, int]:
    return {identifier: len(sequence) for identifier, sequence in _fasta_records(source)}


def get_sequence_by_id(source: str | Path | Sequence[Record], *, identifier: str) -> str:
    return _by_id(source)[identifier]


def interleave_pairs(first: str | Path | Sequence[FastqInput], second: str | Path | Sequence[FastqInput]) -> list[FastqRecord]:
    left, right = _fastq_records(first), _fastq_records(second)
    if len(left) != len(right):
        raise ValueError("paired FASTQ inputs must contain the same number of records")
    result = []
    for first_record, second_record in zip(left, right):
        result.extend((first_record, second_record))
    return result


def iterate(source: str | Path | Sequence[Record] | Sequence[FastqInput]) -> list[Record] | list[FastqRecord]:
    """Materialize records in input order for a streaming-compatible facade."""

    try:
        return _fasta_records(source)  # type: ignore[arg-type]
    except (ValueError, TypeError):
        return _fastq_records(source)  # type: ignore[arg-type]


def mask_regions(source: str | Path | Sequence[Record], regions: Mapping[str, Sequence[Sequence[int]]], *, symbol: str = "N") -> list[Record]:
    result = []
    for identifier, sequence in _fasta_records(source):
        chars = list(sequence)
        for start, end in regions.get(identifier, ()):
            chars[start:end] = [symbol] * max(0, end - start)
        result.append((identifier, "".join(chars)))
    return result


def merge(sources: Sequence[str | Path | Sequence[Record]]) -> list[Record]:
    result = []
    for source in sources:
        result.extend(_fasta_records(source))
    return result


def merge_lanes(sources: Sequence[str | Path | Sequence[FastqInput]]) -> list[FastqRecord]:
    result = []
    for source in sources:
        result.extend(_fastq_records(source))
    return result


def query_fasta_index(index: Mapping[str, Mapping[str, int]], *, identifier: str) -> Mapping[str, int]:
    return index[identifier]


def remove_duplicate_ids(source: str | Path | Sequence[Record]) -> list[Record]:
    result = []
    seen = set()
    for record in _fasta_records(source):
        if record[0] not in seen:
            seen.add(record[0])
            result.append(record)
    return result


def remove_duplicate_sequences(source: str | Path | Sequence[Record]) -> list[Record]:
    result = []
    seen = set()
    for record in _fasta_records(source):
        if record[1] not in seen:
            seen.add(record[1])
            result.append(record)
    return result


def remove_poly_a(source: str | Path | Sequence[FastqInput], *, minimum: int = 5) -> list[FastqRecord]:
    return [record for record in _fastq_records(source) if "A" * minimum not in record.sequence]


def remove_poly_g(source: str | Path | Sequence[FastqInput], *, minimum: int = 5) -> list[FastqRecord]:
    return [record for record in _fastq_records(source) if "G" * minimum not in record.sequence]


def rename_headers(source: str | Path | Sequence[Record], *, prefix: str = "record") -> list[Record]:
    return [(f"{prefix}_{index}", sequence) for index, (_, sequence) in enumerate(_fasta_records(source), start=1)]


def repair_fasta_headers(source: str | Path | Sequence[Record]) -> list[Record]:
    return [(identifier.split()[0] or f"record_{index}", sequence) for index, (identifier, sequence) in enumerate(_fasta_records(source), start=1)]


def repair_pairs(source: str | Path | Sequence[FastqInput]) -> list[FastqRecord]:
    records = _fastq_records(source)
    grouped: dict[str, list[FastqRecord]] = {}
    for record in records:
        grouped.setdefault(re.sub(r"/[12]$", "", record.identifier), []).append(record)
    return [record for pair in grouped.values() if len(pair) == 2 for record in pair]


def reverse_complement_records(source: str | Path | Sequence[Record]) -> list[Record]:
    return [(identifier, reverse_complement_dna(sequence)) for identifier, sequence in _fasta_records(source)]


def sample(source: str | Path | Sequence[FastqInput], *, size: int, seed: int | None = None) -> list[FastqRecord]:
    records = _fastq_records(source)
    if size < 0:
        raise ValueError("size must be non-negative")
    return random.Random(seed).sample(records, min(size, len(records)))


def sample_records(source: str | Path | Sequence[Record], *, size: int, seed: int | None = None) -> list[Record]:
    records = _fasta_records(source)
    return random.Random(seed).sample(records, min(size, len(records)))


def sanitize_headers(source: str | Path | Sequence[Record]) -> list[Record]:
    return [(re.sub(r"[^A-Za-z0-9_.:-]+", "_", identifier), sequence) for identifier, sequence in _fasta_records(source)]


def shuffle_reads(source: str | Path | Sequence[FastqInput], *, seed: int | None = None) -> list[FastqRecord]:
    result = _fastq_records(source)
    random.Random(seed).shuffle(result)
    return result


def shuffle_records(source: str | Path | Sequence[Record], *, seed: int | None = None) -> list[Record]:
    result = _fasta_records(source)
    random.Random(seed).shuffle(result)
    return result


def sort(source: str | Path | Sequence[Record], *, by: str = "id", reverse: bool = False) -> list[Record]:
    records = _fasta_records(source)
    if by not in {"id", "length", "sequence"}:
        raise ValueError("by must be id, length or sequence")
    key = {"id": lambda item: item[0], "length": lambda item: len(item[1]), "sequence": lambda item: item[1]}[by]
    return sorted(records, key=key, reverse=reverse)


def split_by_barcode(source: str | Path | Sequence[FastqInput], barcodes: Mapping[str, str]) -> dict[str, list[FastqRecord]]:
    result = {name: [] for name in barcodes}
    result["unmatched"] = []
    for record in _fastq_records(source):
        matched = False
        for name, barcode in barcodes.items():
            if record.sequence.startswith(normalize_dna_sequence(barcode)):
                result[name].append(record)
                matched = True
                break
        if not matched:
            result["unmatched"].append(record)
    return result


def split_by_chromosome(source: str | Path | Sequence[Record]) -> dict[str, list[Record]]:
    result: dict[str, list[Record]] = {}
    for record in _fasta_records(source):
        chromosome = re.split(r"[:|]", record[0], maxsplit=1)[0]
        result.setdefault(chromosome, []).append(record)
    return result


def split_by_length(source: str | Path | Sequence[Record], *, threshold: int) -> dict[str, list[Record]]:
    return {"short": [record for record in _fasta_records(source) if len(record[1]) < threshold], "long": [record for record in _fasta_records(source) if len(record[1]) >= threshold]}


def split_by_record(source: str | Path | Sequence[Record], *, records_per_part: int) -> list[list[Record]]:
    records = _fasta_records(source)
    if records_per_part < 1:
        raise ValueError("records_per_part must be positive")
    return [records[index : index + records_per_part] for index in range(0, len(records), records_per_part)]


def split_by_size(source: str | Path | Sequence[Record], *, bases_per_part: int) -> list[list[Record]]:
    if bases_per_part < 1:
        raise ValueError("bases_per_part must be positive")
    parts: list[list[Record]] = [[]]
    size = 0
    for record in _fasta_records(source):
        if parts[-1] and size + len(record[1]) > bases_per_part:
            parts.append([])
            size = 0
        parts[-1].append(record)
        size += len(record[1])
    return parts if parts != [[]] else []


def split_by_window(source: str | Path | Sequence[Record], *, window: int, step: int | None = None) -> list[Record]:
    result = []
    for identifier, sequence in _fasta_records(source):
        result.extend((f"{identifier}:{item['start']}-{item['end']}", str(item["sequence"])) for item in split_sequence_windows(sequence, window=window, step=step))
    return result


def stats(source: str | Path | Sequence[Record]) -> dict[str, object]:
    records = _fasta_records(source)
    lengths = [len(sequence) for _, sequence in records]
    return {"records": len(records), "bases": sum(lengths), "lengths": _distribution(lengths)}


def subsample(source: str | Path | Sequence[FastqInput], *, fraction: float, seed: int | None = None) -> list[FastqRecord]:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be between 0 and 1")
    rng = random.Random(seed)
    return [record for record in _fastq_records(source) if rng.random() < fraction]


def subsample_records(source: str | Path | Sequence[Record], *, fraction: float, seed: int | None = None) -> list[Record]:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be between 0 and 1")
    rng = random.Random(seed)
    return [record for record in _fasta_records(source) if rng.random() < fraction]


def synchronize_pairs(first: str | Path | Sequence[FastqInput], second: str | Path | Sequence[FastqInput]) -> tuple[list[FastqRecord], list[FastqRecord]]:
    left = {re.sub(r"/[12]$", "", record.identifier): record for record in _fastq_records(first)}
    right = {re.sub(r"/[12]$", "", record.identifier): record for record in _fastq_records(second)}
    shared = sorted(set(left) & set(right))
    return [left[key] for key in shared], [right[key] for key in shared]


def trim_adapters(source: str | Path | Sequence[FastqInput], *, adapters: Sequence[str]) -> list[FastqRecord]:
    normalized = [normalize_dna_sequence(adapter) for adapter in adapters]
    result = []
    for record in _fastq_records(source):
        sequence, quality = record.sequence, record.quality
        for adapter in normalized:
            if sequence.startswith(adapter):
                sequence, quality = sequence[len(adapter) :], quality[len(adapter) :]
            if sequence.endswith(adapter):
                sequence, quality = sequence[: -len(adapter)], quality[: -len(adapter)]
        result.append(FastqRecord(record.identifier, sequence, quality))
    return result


def trim_low_quality_ends(source: str | Path | Sequence[FastqInput], *, minimum_quality: int = 20) -> list[FastqRecord]:
    result = []
    for record in _fastq_records(source):
        left = 0
        right = len(record.quality)
        while left < right and ord(record.quality[left]) - 33 < minimum_quality:
            left += 1
        while right > left and ord(record.quality[right - 1]) - 33 < minimum_quality:
            right -= 1
        result.append(FastqRecord(record.identifier, record.sequence[left:right], record.quality[left:right]))
    return result


def unmask_regions(source: str | Path | Sequence[Record], regions: Mapping[str, Sequence[Sequence[int]]]) -> list[Record]:
    result = []
    for identifier, sequence in _fasta_records(source):
        chars = list(sequence)
        for start, end in regions.get(identifier, ()):
            chars[start:end] = [base.upper() for base in chars[start:end]]
        result.append((identifier, "".join(chars)))
    return result


def unwrap_fasta_lines(source: str | Path | Sequence[Record]) -> str:
    return _write_fasta(_fasta_records(source), width=10**9)


def validate(source: str | Path | Sequence[Record] | Sequence[FastqInput]) -> dict[str, object]:
    try:
        records = _fasta_records(source)  # type: ignore[arg-type]
        ids = [identifier for identifier, _ in records]
        return {"valid": len(ids) == len(set(ids)), "format": "fasta", "records": len(records), "duplicate_ids": [identifier for identifier, count in Counter(ids).items() if count > 1]}
    except (ValueError, TypeError):
        records = _fastq_records(source)  # type: ignore[arg-type]
        ids = [record.identifier for record in records]
        return {"valid": len(ids) == len(set(ids)), "format": "fastq", "records": len(records), "duplicate_ids": [identifier for identifier, count in Counter(ids).items() if count > 1]}


def validate_sequence_ids(source: str | Path | Sequence[Record]) -> dict[str, object]:
    ids = get_ids(source)
    return {"valid": len(ids) == len(set(ids)), "ids": ids, "duplicates": [identifier for identifier, count in Counter(ids).items() if count > 1]}


def wrap_fasta_lines(source: str | Path | Sequence[Record], *, width: int = 80) -> str:
    return _write_fasta(_fasta_records(source), width=width)


from .external import cluster_vsearch, convert_seqmagick


__all__ = [
    "FastqRecord", "build_fasta_index", "build_fastq_index", "calculate_gc_profile", "cluster_vsearch",
    "calculate_length_distribution", "calculate_n_content_profile", "calculate_quality_distribution",
    "compare_fasta_files", "convert_from_fasta", "convert_from_tabular", "convert_to_fasta",
    "convert_to_tabular", "count_reads", "count_records", "deinterleave_pairs", "estimate_duplication",
    "estimate_error_rate", "extract_cds", "extract_exon", "extract_flanking", "extract_introns",
    "extract_promoter", "extract_read_names", "extract_region", "extract_regions_from_bed", "filter",
    "filter_by_length", "filter_by_quality", "get_ids", "get_lengths", "get_sequence_by_id",
    "interleave_pairs", "iterate", "mask_regions", "merge", "merge_lanes", "query_fasta_index",
    "remove_duplicate_ids", "remove_duplicate_sequences", "remove_poly_a", "remove_poly_g",
    "rename_headers", "repair_fasta_headers", "repair_pairs", "reverse_complement_records", "sample",
    "sample_records", "sanitize_headers", "shuffle_reads", "shuffle_records", "sort", "split_by_barcode",
    "split_by_chromosome", "split_by_length", "split_by_record", "split_by_size", "split_by_window",
    "stats", "subsample", "subsample_records", "synchronize_pairs", "trim_adapters",
    "trim_low_quality_ends", "unmask_regions", "unwrap_fasta_lines", "validate", "validate_sequence_ids",
    "wrap_fasta_lines", "convert_seqmagick",
]


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
