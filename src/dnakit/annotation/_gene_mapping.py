"""CIGAR-aware projection of GFF3 gene intervals onto genomic query DNA."""

from __future__ import annotations

import re
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from dnakit.core import Interval
from dnakit.core._json import FrozenDict
from dnakit.exceptions import CoordinateError, InputFormatError

from .gene_models import GeneAnnotationConfig, GeneRecord, GeneSegment


@dataclass(frozen=True, slots=True)
class ReferenceGene:
    sequence_id: str
    feature_id: str
    gene_id: str
    gene_name: str | None
    locus_tag: str | None
    biotype: str | None
    strand: str
    parts: tuple[Interval, ...]


@dataclass(frozen=True, slots=True)
class MappingBlock:
    query: Interval
    target: Interval


@dataclass(frozen=True, slots=True)
class Mapping:
    id: str
    query: Interval
    target_name: str
    target: Interval
    strand: str
    mapq: int
    identity: float
    primary: bool
    accepted: bool
    blocks: tuple[MappingBlock, ...]

    def to_dict(self) -> FrozenDict:
        return FrozenDict(
            {
                "mapping_id": self.id,
                "query_start": self.query.start + 1,
                "query_end": self.query.end,
                "reference_sequence": self.target_name,
                "reference_start": self.target.start + 1,
                "reference_end": self.target.end,
                "strand": self.strand,
                "mapq": self.mapq,
                "identity": self.identity,
                "primary": self.primary,
                "accepted": self.accepted,
            }
        )


def read_mappings(
    path: Path,
    query_length: int,
    reference_lengths: dict[str, int],
    config: GeneAnnotationConfig,
) -> tuple[Mapping, ...]:
    """Require exact CIGAR consumption; PAF coordinates alone are insufficient."""

    if path.stat().st_size > config.max_alignment_bytes:
        raise InputFormatError("PAF exceeds max_alignment_bytes.")
    results: list[Mapping] = []
    seen: set[tuple[str, ...]] = set()
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            columns = line.rstrip("\r\n").split("\t")
            try:
                if len(columns) < 12 or columns[0] != "query":
                    raise ValueError("Unexpected query or missing columns")
                length, qs, qe = map(int, columns[1:4])
                target_length, ts, te, matches, span, mapq = map(int, columns[6:12])
                strand, target = columns[4:6]
                if (
                    length != query_length
                    or not 0 <= qs < qe <= length
                    or target_length != reference_lengths.get(target)
                    or not 0 <= ts < te <= target_length
                    or strand not in {"+", "-"}
                    or not 0 <= matches <= span
                    or span < 1
                    or not 0 <= mapq <= 255
                ):
                    raise ValueError("Invalid PAF coordinates, reference, or score")
                tags: dict[str, str] = {}
                for raw in columns[12:]:
                    key, value = raw.rsplit(":", 1)
                    if key in tags:
                        raise ValueError("Duplicate PAF tag")
                    tags[key] = value
                cigar = tags["cg:Z"]
                if tags["tp:A"] not in {"P", "S", "I", "i"}:
                    raise ValueError("Invalid minimap2 alignment type")
                primary = tags["tp:A"] in {"P", "I"}
                signature = (*columns[:12], cigar, tags["tp:A"])
                if signature in seen:
                    continue
                seen.add(signature)
                blocks: list[MappingBlock] = []
                qpos, tpos, offset = (qs if strand == "+" else qe), ts, 0
                consumed = 0
                for match in re.finditer(r"([1-9][0-9]*)([MIDN=X])", cigar):
                    if match.start() != offset:
                        raise ValueError("Unsupported or malformed CIGAR")
                    offset = match.end()
                    size, operation = int(match[1]), match[2]
                    if operation != "N":
                        consumed += size
                    if operation in "M=X":
                        qinterval = (
                            Interval(qpos, qpos + size)
                            if strand == "+"
                            else Interval(qpos - size, qpos)
                        )
                        blocks.append(MappingBlock(qinterval, Interval(tpos, tpos + size)))
                    if operation in "MI=X":
                        qpos += size if strand == "+" else -size
                    if operation in "MDN=X":
                        tpos += size
                if (
                    offset != len(cigar)
                    or not blocks
                    or consumed != span
                    or tpos != te
                    or qpos != (qe if strand == "+" else qs)
                ):
                    raise ValueError("CIGAR does not match PAF spans")
            except (ValueError, KeyError, CoordinateError) as exc:
                raise InputFormatError(
                    f"Invalid minimap2 PAF at line {line_number}: {exc}.",
                    code="INVALID_GENE_MAPPING",
                ) from exc
            identity = matches / span
            accepted = (
                primary
                and config.min_mapq <= mapq < 255
                and identity >= config.min_identity
                and sum(len(block.query) for block in blocks) >= config.min_alignment_length
            )
            results.append(
                Mapping(
                    f"mapping_{len(results) + 1}",
                    Interval(qs, qe),
                    target,
                    Interval(ts, te),
                    strand,
                    mapq,
                    identity,
                    primary,
                    accepted,
                    tuple(blocks),
                )
            )
    return tuple(results)


def _project(block: MappingBlock, part: Interval, strand: str) -> Interval | None:
    start, end = max(block.target.start, part.start), min(block.target.end, part.end)
    if end <= start:
        return None
    if strand == "+":
        return Interval(
            block.query.start + start - block.target.start,
            block.query.start + end - block.target.start,
        )
    return Interval(
        block.query.end - (end - block.target.start),
        block.query.end - (start - block.target.start),
    )


def _merge(parts: list[Interval]) -> tuple[Interval, ...]:
    merged: list[Interval] = []
    for part in sorted(parts, key=lambda item: (item.start, item.end)):
        if merged and part.start <= merged[-1].end:
            merged[-1] = Interval(merged[-1].start, max(part.end, merged[-1].end))
        else:
            merged.append(part)
    return tuple(merged)


def annotate_mappings(
    length: int,
    mappings: tuple[Mapping, ...],
    genes: tuple[ReferenceGene, ...],
    annotated_contigs: set[str],
) -> tuple[tuple[GeneRecord, ...], tuple[GeneSegment, ...]]:
    """Sweep mapping/gene boundaries, retaining overlapping genes on a unique mapping."""

    # Events are (kind, id, delta). Sets alone cannot represent touching CIGAR blocks.
    events: dict[int, list[tuple[str, str, int]]] = defaultdict(list)
    events[0] = []
    events[length] = []

    def event(part: Interval, kind: str, identifier: str) -> None:
        events[part.start].append((kind, identifier, 1))
        events[part.end].append((kind, identifier, -1))

    by_contig: dict[str, list[ReferenceGene]] = defaultdict(list)
    for gene in genes:
        by_contig[gene.sequence_id].append(gene)
    candidates: dict[str, GeneRecord] = {}
    mapping_by_id = {mapping.id: mapping for mapping in mappings}
    for mapping in mappings:
        event(mapping.query, "observed", mapping.id)
        if not mapping.accepted:
            continue
        for block in mapping.blocks:
            event(block.query, "mapped", mapping.id)
        block_ends = [block.target.end for block in mapping.blocks]
        for gene in by_contig[mapping.target_name]:
            parts = []
            for ref_part in gene.parts:
                if ref_part.end <= mapping.target.start or ref_part.start >= mapping.target.end:
                    continue
                first = bisect_right(block_ends, ref_part.start)
                for block_index in range(first, len(mapping.blocks)):
                    block = mapping.blocks[block_index]
                    if block.target.start >= ref_part.end:
                        break
                    projected = _project(block, ref_part, mapping.strand)
                    if projected is not None:
                        parts.append(projected)
            if not parts:
                continue
            identifier = f"gene_{len(candidates) + 1}"
            strand = gene.strand
            if mapping.strand == "-" and strand in {"+", "-"}:
                strand = "+" if strand == "-" else "-"
            candidates[identifier] = GeneRecord(
                identifier,
                gene.gene_id,
                gene.gene_name,
                gene.locus_tag,
                gene.feature_id,
                gene.biotype,
                strand,
                gene.strand,
                gene.sequence_id,
                gene.parts,
                _merge(parts),
                mapping.id,
            )
            for part in candidates[identifier].parts:
                event(part, "gene", identifier)

    counts: dict[str, dict[str, int]] = {"mapped": {}, "observed": {}, "gene": {}}
    used_parts: dict[str, list[Interval]] = defaultdict(list)
    segments: list[GeneSegment] = []
    positions = sorted(events)
    for index, position in enumerate(positions[:-1]):
        for kind, identifier, delta in events[position]:
            counts[kind][identifier] = counts[kind].get(identifier, 0) + delta
            if counts[kind][identifier] == 0:
                del counts[kind][identifier]
        location = Interval(position, positions[index + 1])
        active = counts["mapped"]
        mapping_id = next(iter(active)) if len(active) == 1 else None
        if mapping_id is None:
            reason = (
                "ambiguous_mapping"
                if len(active) > 1
                else "unaligned_or_low_confidence"
                if counts["observed"]
                else "unmapped"
            )
            segment = GeneSegment(location, "unknown", reason=reason)
        elif mapping_by_id[mapping_id].target_name not in annotated_contigs:
            segment = GeneSegment(
                location, "unknown", mapping_id=mapping_id, reason="missing_contig_annotation"
            )
        else:
            ids = tuple(
                sorted(key for key in counts["gene"] if candidates[key].mapping_id == mapping_id)
            )
            for key in ids:
                used_parts[key].append(location)
            segment = GeneSegment(
                location,
                "gene" if ids else "intergenic",
                ids,
                tuple(sorted({candidates[key].gene_id for key in ids})),
                tuple(
                    sorted({name for key in ids if (name := candidates[key].gene_name) is not None})
                ),
                mapping_id,
            )
        if segments and replace(segment, location=segments[-1].location) == segments[-1]:
            segments[-1] = replace(
                segment, location=Interval(segments[-1].location.start, location.end)
            )
        else:
            segments.append(segment)
    records = tuple(
        sorted(
            (replace(candidates[key], parts=_merge(parts)) for key, parts in used_parts.items()),
            key=lambda gene: (gene.location.start, gene.location.end, gene.record_id),
        )
    )
    return records, tuple(segments)
