"""Configuration and coordinate-explicit results for reference-based gene annotation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal

from dnakit.core import Interval
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class GeneAnnotationConfig:
    """Limits for genomic DNA mapping; imports never download or run minimap2.

    A cache directory pins the first successful assembly/annotation selection for
    each species and assembly request. Use a new directory to refresh that pin.
    ``asm5`` is intended for long, accurate genomic DNA, not spliced cDNA.
    """

    cache_dir: str = ".dnakit/gene_annotation"
    minimap2_executable: str = "minimap2"
    preset: str = "asm5"
    threads: int = 1
    min_mapq: int = 20
    min_identity: float = 0.90
    min_alignment_length: int = 100
    timeout_seconds: float = 3600.0
    allow_download: bool = True
    api_key: str | None = None
    max_query_bases: int = 100_000_000
    max_reference_bases: int = 20_000_000_000
    max_download_bytes: int = 20_000_000_000
    max_gff_records: int = 10_000_000
    max_alignment_bytes: int = 100_000_000
    max_index_bytes: int = 10_000_000_000

    def __post_init__(self) -> None:
        for name in ("cache_dir", "minimap2_executable"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or "\x00" in value:
                raise ConfigurationError(f"{name} must be non-empty text.")
        if self.preset not in {"asm5", "asm10", "asm20", "map-hifi", "map-ont", "sr"}:
            raise ConfigurationError("Unsupported genomic minimap2 preset.")
        for name, maximum, minimum in (
            ("threads", 256, 1),
            ("min_mapq", 254, 1),
            ("min_alignment_length", 100_000_000, 1),
            ("max_query_bases", 1_000_000_000, 1),
            ("max_reference_bases", 100_000_000_000, 1),
            ("max_download_bytes", 1_000_000_000_000, 1),
            ("max_gff_records", 100_000_000, 1),
            ("max_alignment_bytes", 100_000_000, 1),
            ("max_index_bytes", 10_000_000_000, 1),
        ):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or not minimum <= value <= maximum
            ):
                raise ConfigurationError(f"{name} must be in [{minimum}, {maximum}].")
        for name, float_maximum in (("min_identity", 1.0), ("timeout_seconds", 86_400.0)):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or not 0 < value <= float_maximum
            ):
                raise ConfigurationError(f"{name} must be in (0, {float_maximum}].")
        if not isinstance(self.allow_download, bool):
            raise ConfigurationError("allow_download must be boolean.")
        if self.api_key is not None and (
            not isinstance(self.api_key, str) or not self.api_key.strip()
        ):
            raise ConfigurationError("api_key must be None or non-empty text.")


def _interval(interval: Interval) -> dict[str, int]:
    return {"start": interval.start + 1, "end": interval.end}


@dataclass(frozen=True, slots=True)
class GeneRecord:
    """One projected gene instance; internal intervals are 0-based half-open.

    ``location`` is the envelope; ``parts`` contains only reliably mapped bases.
    Insertions and ambiguous mappings can leave holes inside the envelope.
    """

    record_id: str
    gene_id: str
    gene_name: str | None
    locus_tag: str | None
    feature_id: str
    biotype: str | None
    strand: str
    reference_strand: str
    reference_sequence: str
    reference_parts: tuple[Interval, ...]
    parts: tuple[Interval, ...]
    mapping_id: str

    @property
    def location(self) -> Interval:
        return Interval(self.parts[0].start, self.parts[-1].end)

    def to_dict(self) -> dict[str, object]:
        return {
            "record_id": self.record_id,
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "display_name": self.gene_name or self.locus_tag or self.gene_id,
            "locus_tag": self.locus_tag,
            "feature_id": self.feature_id,
            "biotype": self.biotype,
            "strand": self.strand,
            "reference_strand": self.reference_strand,
            "reference_sequence": self.reference_sequence,
            "reference_parts": [_interval(part) for part in self.reference_parts],
            **_interval(self.location),
            "parts": [_interval(part) for part in self.parts],
            "mapping_id": self.mapping_id,
        }


@dataclass(frozen=True, slots=True)
class GeneSegment:
    """A disjoint input interval; an unknown interval is never called intergenic."""

    location: Interval
    kind: Literal["gene", "intergenic", "unknown"]
    gene_record_ids: tuple[str, ...] = ()
    gene_ids: tuple[str, ...] = ()
    gene_names: tuple[str, ...] = ()
    mapping_id: str | None = None
    reason: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            **_interval(self.location),
            "type": self.kind,
            "gene_record_ids": list(self.gene_record_ids),
            "gene_ids": list(self.gene_ids),
            "gene_names": list(self.gene_names),
            "mapping_id": self.mapping_id,
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class GeneAnnotationResult:
    """Named genes, complete input segmentation, and pinned reference provenance."""

    sequence_id: str
    sequence_length: int
    reference: FrozenDict
    gene_records: tuple[GeneRecord, ...]
    segments: tuple[GeneSegment, ...]
    alignments: tuple[FrozenDict, ...]
    metadata: FrozenDict

    def to_dict(self) -> dict[str, Any]:
        """Export all coordinates as 1-based closed, including reverse mappings."""

        return {
            "sequence_id": self.sequence_id,
            "sequence_length": self.sequence_length,
            "coordinate_system": "1-based-closed",
            "reference": to_json_compatible(self.reference),
            "gene_records": [record.to_dict() for record in self.gene_records],
            "segments": [segment.to_dict() for segment in self.segments],
            "alignments": [to_json_compatible(item) for item in self.alignments],
            "metadata": to_json_compatible(self.metadata),
        }


__all__ = ["GeneAnnotationConfig", "GeneAnnotationResult", "GeneRecord", "GeneSegment"]
