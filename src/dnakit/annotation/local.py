"""Dependency-free interval operations for explicit gene annotations.

These functions operate on materialized, validated feature mappings.  They do
not infer gene structure, query an annotation database, or inspect a reference
genome.  Coordinate conventions are explicit so that a paper can reproduce
the exact interval arithmetic.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any

from dnakit.exceptions import ConfigurationError

_ZERO_BASED = "0-based-half-open"
_ONE_BASED = "1-based-closed"
_COORDINATE_SYSTEMS = {_ZERO_BASED, _ONE_BASED}


def _envelope(operation: str, result: Mapping[str, Any], *, coordinate_system: str | None = None) -> dict[str, Any]:
    provenance: dict[str, Any] = {
        "implementation": "dnakit.native_annotation_intervals",
        "backend": "local_dependency_free",
        "scientific_boundary": "deterministic interval arithmetic; no reference or annotation database",
    }
    if coordinate_system is not None:
        provenance["coordinate_system"] = coordinate_system
    return {
        "status": "completed",
        "operation": operation,
        "result": dict(result),
        "provenance": provenance,
    }


def _coordinate_system(options: Mapping[str, object], *, default: str = _ONE_BASED) -> str:
    value = options.get("coordinate_system", default)
    if value not in _COORDINATE_SYSTEMS:
        raise ConfigurationError(
            "options.coordinate_system must be '0-based-half-open' or '1-based-closed'."
        )
    return str(value)


def _records(value: object, *, key_names: tuple[str, ...] = ("records", "features", "genes", "items")) -> list[dict[str, Any]]:
    if isinstance(value, Mapping):
        for key in key_names:
            if key in value:
                value = value[key]
                break
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("input must be a sequence of annotation mappings or a mapping containing records.")
    records: list[dict[str, Any]] = []
    for index, record in enumerate(value):
        if not isinstance(record, Mapping):
            raise ConfigurationError(f"annotation record {index} must be a mapping.")
        records.append(dict(record))
    return records


def _coordinate(record: Mapping[str, Any], coordinate_system: str) -> tuple[str, int, int]:
    chrom = record.get("chrom", record.get("seqid", record.get("chromosome")))
    start = record.get("start")
    end = record.get("end")
    if not isinstance(chrom, str) or not chrom.strip():
        raise ConfigurationError("annotation records require a non-empty chrom/seqid/chromosome.")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (start, end)):
        raise ConfigurationError("annotation start and end must be integers.")
    assert isinstance(start, int) and isinstance(end, int)
    minimum = 0 if coordinate_system == _ZERO_BASED else 1
    if start < minimum or end < start or (coordinate_system == _ONE_BASED and end < 1):
        raise ConfigurationError(f"invalid {coordinate_system} interval: {chrom}:{start}-{end}.")
    return chrom.strip(), start, end


def _strand(record: Mapping[str, Any]) -> str:
    value = record.get("strand", "+")
    if value not in {"+", "-", ".", "?"}:
        raise ConfigurationError("strand must be '+', '-', '.', or '?'.")
    return str(value)


def _contig_lengths(input: object, options: Mapping[str, object]) -> dict[str, int]:
    source: object = options.get("contig_lengths")
    if source is None and isinstance(input, Mapping):
        source = input.get("contig_lengths", input.get("contigs"))
    if source is None:
        raise ConfigurationError("contig_lengths are required for intergenic-region generation.")
    if isinstance(source, Mapping):
        pairs = source.items()
    elif isinstance(source, Sequence) and not isinstance(source, (str, bytes)):
        pairs = []
        for record in source:
            if not isinstance(record, Mapping):
                raise ConfigurationError("contig length records must be mappings.")
            pairs.append((record.get("chrom", record.get("seqid")), record.get("length")))
    else:
        raise ConfigurationError("contig_lengths must be a mapping or sequence of mappings.")
    result: dict[str, int] = {}
    for chrom, length in pairs:
        if not isinstance(chrom, str) or not chrom.strip() or isinstance(length, bool) or not isinstance(length, int) or length < 1:
            raise ConfigurationError("contig lengths require a name and a positive integer length.")
        result[chrom.strip()] = length
    if not result:
        raise ConfigurationError("contig_lengths cannot be empty.")
    return result


def convert_annotation_coordinates(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Convert explicit intervals between 0-based half-open and 1-based closed coordinates."""

    settings = dict(options or {})
    source = _coordinate_system(settings, default=_ONE_BASED)
    target = settings.get("target_coordinate_system", _ZERO_BASED if source == _ONE_BASED else _ONE_BASED)
    if target not in _COORDINATE_SYSTEMS:
        raise ConfigurationError(
            "options.target_coordinate_system must be '0-based-half-open' or '1-based-closed'."
        )
    target = str(target)
    records = _records(input)
    converted: list[dict[str, Any]] = []
    for record in records:
        chrom, start, end = _coordinate(record, source)
        row = deepcopy(record)
        row["chrom"] = chrom
        if source == _ZERO_BASED and target == _ONE_BASED:
            row["start"], row["end"] = start + 1, end
        elif source == _ONE_BASED and target == _ZERO_BASED:
            row["start"], row["end"] = start - 1, end
        else:
            row["start"], row["end"] = start, end
        row["coordinate_system"] = target
        converted.append(row)
    return _envelope(
        "convert_annotation_coordinates",
        {"records": converted, "record_count": len(converted), "source_coordinate_system": source, "target_coordinate_system": target},
        coordinate_system=target,
    )


def make_intergenic_regions(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Create gaps between explicit gene intervals and contig boundaries."""

    settings = dict(options or {})
    coordinate_system = _coordinate_system(settings)
    merge_adjacent = settings.get("merge_adjacent", True)
    if not isinstance(merge_adjacent, bool):
        raise ConfigurationError("options.merge_adjacent must be boolean.")
    contigs = _contig_lengths(input, settings)
    genes = _records(input, key_names=("genes", "records", "features", "items"))
    grouped: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in contigs}
    for record in genes:
        chrom, start, end = _coordinate(record, coordinate_system)
        if chrom not in contigs:
            raise ConfigurationError(f"gene interval uses unknown contig {chrom!r}.")
        length = contigs[chrom]
        if end > length:
            raise ConfigurationError(f"gene interval {chrom}:{start}-{end} exceeds contig length {length}.")
        grouped[chrom].append((start, end))
    result: list[dict[str, Any]] = []
    for chrom, length in contigs.items():
        intervals = sorted(grouped[chrom])
        merged: list[list[int]] = []
        for start, end in intervals:
            overlaps = bool(merged and start <= merged[-1][1])
            adjacent = bool(
                merged
                and merge_adjacent
                and start == merged[-1][1] + (1 if coordinate_system == _ONE_BASED else 0)
            )
            if overlaps or adjacent:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        cursor = 1 if coordinate_system == _ONE_BASED else 0
        index = 1
        for start, end in merged:
            if cursor < start:
                result.append({"id": f"intergenic_{chrom}_{index}", "chrom": chrom, "start": cursor, "end": start - 1 if coordinate_system == _ONE_BASED else start, "feature_type": "intergenic", "coordinate_system": coordinate_system})
                index += 1
            cursor = max(cursor, end + 1 if coordinate_system == _ONE_BASED else end)
        if cursor <= length - (0 if coordinate_system == _ZERO_BASED else 1):
            result.append({"id": f"intergenic_{chrom}_{index}", "chrom": chrom, "start": cursor, "end": length, "feature_type": "intergenic", "coordinate_system": coordinate_system})
    return _envelope(
        "make_intergenic_regions",
        {
            "regions": result,
            "region_count": len(result),
            "contig_lengths": contigs,
            "gene_count": len(genes),
            "merge_adjacent": merge_adjacent,
        },
        coordinate_system=coordinate_system,
    )


def _length_option(settings: Mapping[str, object], name: str, default: int) -> int:
    value = settings.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"options.{name} must be a non-negative integer.")
    return value


def make_promoter_regions(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Generate strand-aware promoter windows around explicit gene intervals."""

    settings = dict(options or {})
    coordinate_system = _coordinate_system(settings)
    upstream = _length_option(settings, "upstream", 1000)
    downstream = _length_option(settings, "downstream", 100)
    records = _records(input, key_names=("genes", "records", "features", "items"))
    contig_lengths = None
    if "contig_lengths" in settings or isinstance(input, Mapping) and ("contig_lengths" in input or "contigs" in input):
        contig_lengths = _contig_lengths(input, settings)
    promoters: list[dict[str, Any]] = []
    for index, gene in enumerate(records, 1):
        chrom, start, end = _coordinate(gene, coordinate_system)
        strand = _strand(gene)
        if strand == "-":
            raw_start = end - downstream + 1 if coordinate_system == _ONE_BASED else end - downstream
            raw_end = end + upstream
        else:
            raw_start = start - upstream if coordinate_system == _ZERO_BASED else start - upstream
            raw_end = start + downstream if coordinate_system == _ZERO_BASED else start + downstream - 1
        if coordinate_system == _ONE_BASED:
            raw_start = max(1, raw_start)
        else:
            raw_start = max(0, raw_start)
        if contig_lengths is not None:
            if chrom not in contig_lengths:
                raise ConfigurationError(f"gene interval uses unknown contig {chrom!r}.")
            raw_end = min(contig_lengths[chrom], raw_end)
        if raw_end < raw_start:
            raise ConfigurationError(f"promoter window is empty for gene {gene.get('id', index)!r}.")
        promoters.append({
            "id": f"promoter_{gene.get('id', index)}",
            "gene_id": gene.get("id", gene.get("gene_id", index)),
            "chrom": chrom,
            "start": raw_start,
            "end": raw_end,
            "strand": strand,
            "feature_type": "promoter",
            "coordinate_system": coordinate_system,
            "upstream": upstream,
            "downstream": downstream,
        })
    return _envelope(
        "make_promoter_regions",
        {"promoters": promoters, "promoter_count": len(promoters), "upstream": upstream, "downstream": downstream},
        coordinate_system=coordinate_system,
    )


def merge_gene_models(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Merge explicit exon/segment intervals by gene identity."""

    settings = dict(options or {})
    coordinate_system = _coordinate_system(settings)
    merge_adjacent = settings.get("merge_adjacent", True)
    if not isinstance(merge_adjacent, bool):
        raise ConfigurationError("options.merge_adjacent must be boolean.")
    records = _records(input, key_names=("segments", "exons", "records", "features", "items"))
    grouped: dict[tuple[str, str, str], list[tuple[int, int]]] = {}
    names: dict[tuple[str, str, str], str] = {}
    for record in records:
        chrom, start, end = _coordinate(record, coordinate_system)
        gene_id = record.get("gene_id", record.get("gene", record.get("id")))
        if not isinstance(gene_id, str) or not gene_id.strip():
            raise ConfigurationError("gene model records require gene_id, gene, or id.")
        strand = _strand(record)
        key = (chrom, gene_id.strip(), strand)
        grouped.setdefault(key, []).append((start, end))
        names[key] = str(record.get("gene_name", gene_id)).strip()
    models: list[dict[str, Any]] = []
    for (chrom, gene_id, strand), intervals in sorted(grouped.items()):
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            touches = start <= merged[-1][1] if merged else False
            if merged and (touches or (merge_adjacent and start == merged[-1][1] + 1)):
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        models.append({
            "gene_id": gene_id,
            "gene_name": names[(chrom, gene_id, strand)],
            "chrom": chrom,
            "strand": strand,
            "start": min(interval[0] for interval in merged),
            "end": max(interval[1] for interval in merged),
            "segments": [{"start": start, "end": end} for start, end in merged],
            "segment_count": len(merged),
            "coordinate_system": coordinate_system,
        })
    return _envelope(
        "merge_gene_models",
        {"models": models, "model_count": len(models), "input_segment_count": len(records), "merge_adjacent": merge_adjacent},
        coordinate_system=coordinate_system,
    )


__all__ = [
    "convert_annotation_coordinates",
    "make_intergenic_regions",
    "make_promoter_regions",
    "merge_gene_models",
]
