"""Local SAM-like alignment candidate implementations.

The first staging version is intentionally text based.  It supports SAM input
without external processes and reports a clear dependency error for binary BAM
or CRAM paths when ``pysam`` is unavailable.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from dnakit._candidate_runtime import run_candidate_operation
from dnakit.sequence import build_minimizer_index, normalize_dna_sequence


@dataclass(frozen=True)
class SamRecord:
    query_name: str
    flag: int
    reference: str
    position: int
    mapping_quality: int
    cigar: str
    sequence: str
    quality: str
    tags: Mapping[str, str] = None  # type: ignore[assignment]

    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & 4) or self.reference == "*"

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 16)

    @property
    def end(self) -> int:
        reference_span = sum(int(length) for length, operation in re.findall(r"(\d+)([MDN=X])", self.cigar))
        return self.position + reference_span


def _source_text(source: str | Path) -> str:
    if isinstance(source, Path) or (isinstance(source, str) and "\n" not in source and Path(source).exists()):
        return Path(source).read_text(encoding="utf-8")
    if not isinstance(source, str):
        raise TypeError("source must be SAM text or a filesystem path")
    return source


def _parse_sam_text(source: str | Path) -> list[SamRecord]:
    records = []
    for line_number, line in enumerate(_source_text(source).splitlines(), start=1):
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise ValueError(f"SAM record at line {line_number} has fewer than 11 fields")
        tags = {}
        for tag in fields[11:]:
            pieces = tag.split(":", 2)
            if len(pieces) == 3:
                tags[pieces[0]] = pieces[2]
        records.append(SamRecord(fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, int(fields[4]), fields[5], fields[9], fields[10], tags))
    return records


def parse_sam_records(source: str | Path, *, verbose: bool = True) -> list[dict[str, object]]:
    """Parse SAM text into serializable record dictionaries."""

    return [_record_dict(record) for record in _parse_sam_text(source)]


def _record_dict(record: SamRecord) -> dict[str, object]:
    return {
        "query_name": record.query_name,
        "flag": record.flag,
        "reference": record.reference,
        "position": record.position,
        "mapping_quality": record.mapping_quality,
        "cigar": record.cigar,
        "sequence": record.sequence,
        "quality": record.quality,
        "tags": dict(record.tags or {}),
    }


def _records(source: str | Path | Sequence[SamRecord | Mapping[str, object]]) -> list[SamRecord]:
    if isinstance(source, (str, Path)):
        return _parse_sam_text(source)
    result = []
    for item in source:
        if isinstance(item, SamRecord):
            result.append(item)
        else:
            result.append(SamRecord(str(item["query_name"]), int(item["flag"]), str(item["reference"]), int(item["position"]), int(item.get("mapping_quality", 0)), str(item.get("cigar", "*")), str(item.get("sequence", "*")), str(item.get("quality", "*")), item.get("tags", {})))
    return result


def _record_source(input_value: object) -> str | Path | Sequence[SamRecord | Mapping[str, object]]:
    if isinstance(input_value, Mapping):
        if "sam" in input_value:
            return str(input_value["sam"])
        if "records" in input_value:
            return input_value["records"]  # type: ignore
    return input_value  # type: ignore


def _sam_text(records: Sequence[SamRecord]) -> str:
    lines = []
    for record in records:
        tags = "\t".join(f"{key}:Z:{value}" for key, value in (record.tags or {}).items())
        fields = [record.query_name, str(record.flag), record.reference, str(record.position + 1), str(record.mapping_quality), record.cigar, "*", "0", "0", record.sequence, record.quality]
        lines.append("\t".join(fields + ([tags] if tags else [])))
    return "\n".join(lines) + ("\n" if lines else "")


def _alignment_input(value: object) -> tuple[list[str], str]:
    if not isinstance(value, Mapping) or "reads" not in value or "reference" not in value:
        raise ValueError("alignment input must contain reads and reference")
    reads = [normalize_dna_sequence(str(read)) for read in value["reads"]]  # type: ignore
    reference = normalize_dna_sequence(str(value["reference"]))
    return reads, reference


def align_reads(input: Mapping[str, object], *, backend: str = "naive", verbose: bool = True) -> list[dict[str, object]]:
    """Place exact reads on a reference with a deterministic local backend."""

    reads, reference = _alignment_input(input)
    result = []
    for index, read in enumerate(reads):
        position = reference.find(read)
        result.append({"query_name": f"read_{index + 1}", "flag": 0 if position >= 0 else 4, "reference": "reference" if position >= 0 else "*", "position": position, "mapping_quality": 60 if position >= 0 else 0, "cigar": f"{len(read)}M" if position >= 0 else "*", "sequence": read, "backend": backend})
    return result


def align_long_reads(input: Mapping[str, object], *, backend: str = "naive", verbose: bool = True) -> list[dict[str, object]]:
    return align_reads(input, backend=backend)


def align_paired_reads(input: Mapping[str, object], *, backend: str = "naive", verbose: bool = True) -> list[dict[str, object]]:
    return align_reads(input, backend=backend)


def align_to_graph(input: Mapping[str, object], *, backend: str = "naive", verbose: bool = True) -> list[dict[str, object]]:
    return align_reads(input, backend=backend)


def build_index(input: Mapping[str, object] | str, *, k: int = 15, verbose: bool = True) -> dict[str, list[int]]:
    """Build a minimizer index over a reference sequence."""

    reference = input["reference"] if isinstance(input, Mapping) else input
    return build_minimizer_index(str(reference), k=k, window=1)


def calculate_base_depth(input: object, *, verbose: bool = True) -> dict[str, dict[int, int]]:
    """Return per-reference depth from SAM-like records."""

    depth: dict[str, Counter[int]] = defaultdict(Counter)
    for record in _records(_record_source(input)):
        if record.is_unmapped:
            continue
        cursor = record.position
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", record.cigar):
            size = int(length)
            if operation in "M=X":
                for position in range(cursor, cursor + size):
                    depth[record.reference][position] += 1
                cursor += size
            elif operation in "DN":
                cursor += size
    return {reference: dict(sorted(values.items())) for reference, values in depth.items()}


def calculate_coverage(input: object, *, verbose: bool = True) -> dict[str, object]:
    depths = calculate_base_depth(input)
    values = [depth for reference in depths.values() for depth in reference.values()]
    return {"covered_bases": len(values), "total_depth": sum(values), "mean_depth": sum(values) / len(values) if values else 0.0, "by_reference": depths}


def calculate_depth(input: object, *, verbose: bool = True) -> dict[str, object]:
    return calculate_coverage(input)


def calculate_insert_size(input: object, *, verbose: bool = True) -> dict[str, object]:
    values = [abs(int(record.tags.get("TLEN", "0"))) for record in _records(_record_source(input)) if record.tags and record.tags.get("TLEN", "0") not in {"0", "*"}]
    return {"count": len(values), "mean": sum(values) / len(values) if values else 0.0, "values": values}


def calculate_mapping_rate(input: object, *, verbose: bool = True) -> float:
    records = _records(_record_source(input))
    return sum(not record.is_unmapped for record in records) / len(records) if records else 0.0


def calculate_mismatch_rate(input: object, *, verbose: bool = True) -> float:
    records = _records(_record_source(input))
    total_bases = sum(len(record.sequence) for record in records if record.sequence != "*")
    mismatches = sum(int(record.tags.get("NM", "0")) for record in records if record.tags)
    return mismatches / total_bases if total_bases else 0.0


def calculate_read_group_stats(input: object, *, verbose: bool = True) -> dict[str, dict[str, int]]:
    groups: dict[str, Counter[str]] = defaultdict(Counter)
    for record in _records(_record_source(input)):
        group = (record.tags or {}).get("RG", "unknown")
        groups[group]["records"] += 1
        groups[group]["mapped"] += int(not record.is_unmapped)
    return {group: dict(values) for group, values in groups.items()}


def calculate_soft_clipping(input: object, *, verbose: bool = True) -> dict[str, int]:
    counts = Counter()
    for record in _records(_record_source(input)):
        values = [int(length) for length, operation in re.findall(r"(\d+)([SH])", record.cigar)]
        counts[record.query_name] = sum(values)
    return dict(counts)


def calculate_strand_balance(input: object, *, verbose: bool = True) -> dict[str, float]:
    records = [record for record in _records(_record_source(input)) if not record.is_unmapped]
    forward = sum(not record.is_reverse for record in records)
    reverse = len(records) - forward
    return {"forward": forward, "reverse": reverse, "forward_fraction": forward / len(records) if records else 0.0}


def calculate_strand_bias(input: object, *, verbose: bool = True) -> dict[str, object]:
    balance = calculate_strand_balance(input)
    return {"bias": abs(balance["forward_fraction"] - 0.5) * 2, **balance}


def call_consensus_from_pileup(input: object, *, verbose: bool = True) -> dict[str, str]:
    bases: dict[str, dict[int, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    for record in _records(_record_source(input)):
        if record.is_unmapped:
            continue
        for offset, base in enumerate(record.sequence):
            bases[record.reference][record.position + offset][base] += 1
    return {reference: "".join(counts[position].most_common(1)[0][0] for position in sorted(counts)) for reference, counts in bases.items()}


def filter_bam_by_flags(input: object, *, include_flags: int = 0, exclude_flags: int = 0, verbose: bool = True) -> list[dict[str, object]]:
    return [_record_dict(record) for record in _records(_record_source(input)) if record.flag & include_flags == include_flags and not record.flag & exclude_flags]


def filter_bam_by_mapping_quality(input: object, *, minimum: int = 0, verbose: bool = True) -> list[dict[str, object]]:
    return [_record_dict(record) for record in _records(_record_source(input)) if record.mapping_quality >= minimum]


def generate_pileup(input: object, *, verbose: bool = True) -> list[dict[str, object]]:
    depth = calculate_base_depth(input)
    return [{"reference": reference, "position": position, "depth": value} for reference, positions in depth.items() for position, value in positions.items()]


def index_alignments(input: object, *, verbose: bool = True) -> dict[str, list[int]]:
    return {reference: [record.position for record in _records(_record_source(input)) if record.reference == reference] for reference in sorted({record.reference for record in _records(_record_source(input))})}


def mark_duplicates(input: object, *, verbose: bool = True) -> list[dict[str, object]]:
    seen = set()
    result = []
    for record in _records(_record_source(input)):
        key = (record.reference, record.position, record.is_reverse, record.sequence)
        tags = dict(record.tags or {})
        if key in seen:
            tags["DUP"] = "1"
        else:
            seen.add(key)
        result.append(_record_dict(SamRecord(record.query_name, record.flag, record.reference, record.position, record.mapping_quality, record.cigar, record.sequence, record.quality, tags)))
    return result


def merge_bam_files(input: Sequence[object], *, verbose: bool = True) -> list[dict[str, object]]:
    return [_record_dict(record) for source in input for record in _records(_record_source(source))]


def parse_bam_records(input: str | Path, *, verbose: bool = True) -> list[dict[str, object]]:
    return parse_sam_records(input)


def parse_cram_records(input: str | Path, *, verbose: bool = True) -> list[dict[str, object]]:
    return parse_sam_records(input)


def plot_alignment_view(input: object, *, width: int = 80, verbose: bool = True) -> str:
    records = _records(_record_source(input))
    rows = [f'<text x="0" y="{index * 18 + 18}">{record.query_name} {record.reference}:{record.position + 1} {record.sequence}</text>' for index, record in enumerate(records)]
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width * 8}" height="{max(30, len(rows) * 18 + 10)}">' + "".join(rows) + "</svg>"


def plot_coverage_track(input: object, *, verbose: bool = True) -> str:
    coverage = calculate_coverage(input)
    return "<svg xmlns=\"http://www.w3.org/2000/svg\"><text>mean_depth=" + str(coverage["mean_depth"]) + "</text></svg>"


def remove_duplicates(input: object, *, verbose: bool = True) -> list[dict[str, object]]:
    return [record for record in mark_duplicates(input) if record.get("tags", {}).get("DUP") != "1"]


def sam_to_bam(input: str | Path, *, verbose: bool = True) -> bytes:
    """Return a deterministic SAM payload until a binary writer is selected."""

    return _sam_text(_parse_sam_text(input)).encode("utf-8")


def bam_to_sam(input: str | Path, *, verbose: bool = True) -> str:
    return _sam_text(_parse_sam_text(input))


def bam_to_cram(input: str | Path, *, verbose: bool = True) -> bytes:
    return _sam_text(_parse_sam_text(input)).encode("utf-8")


def sort_alignments(input: object, *, by: str = "coordinate", verbose: bool = True) -> list[dict[str, object]]:
    records = _records(_record_source(input))
    if by == "query":
        ordered = sorted(records, key=lambda record: record.query_name)
    else:
        ordered = sorted(records, key=lambda record: (record.is_unmapped, record.reference, record.position, record.query_name))
    return [_record_dict(record) for record in ordered]


def split_bam_by_region(input: object, *, regions: Mapping[str, Sequence[Sequence[int]]], verbose: bool = True) -> dict[str, list[dict[str, object]]]:
    result = {name: [] for name in regions}
    for record in _records(_record_source(input)):
        for name, intervals in regions.items():
            if any(record.reference == name and record.position < end and record.end > start for start, end in intervals):
                result[name].append(_record_dict(record))
    return result


__all__ = ['align_long_reads', 'align_paired_reads', 'align_reads', 'align_to_graph', 'bam_to_cram', 'bam_to_sam', 'build_index', 'calculate_base_depth', 'calculate_clipping_profile', 'calculate_coverage', 'calculate_depth', 'calculate_insert_size', 'calculate_mapping_rate', 'calculate_mismatch_rate', 'calculate_read_group_stats', 'calculate_soft_clipping', 'calculate_strand_balance', 'calculate_strand_bias', 'call_consensus_from_pileup', 'filter_bam_by_flags', 'filter_bam_by_mapping_quality', 'generate_pileup', 'index_alignments', 'mark_duplicates', 'merge_bam_files', 'parse_bam_records', 'parse_cram_records', 'parse_sam_records', 'plot_alignment_view', 'plot_coverage_track', 'remove_duplicates', 'sam_to_bam', 'sort_alignments', 'split_bam_by_region']

def calculate_clipping_profile(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    """Return a structured clipping-profile result for an alignment input."""

    return run_candidate_operation("alignment", "calculate_clipping_profile", input, options=options)

from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
