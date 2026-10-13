"""umi_tools functions added by review group G18o (sequencing_qc area).

Read-level helpers take plain values: tag dictionaries, CIGAR strings or read names. Records are
built as in-memory pysam ``AlignedSegment`` objects, so no BAM file is needed for them. The two
metacontig functions open a local SAM or BAM file (the tag-level one needs a BAM index). Nothing
here writes an alignment file or changes global state. umi_tools and pysam are imported only when
a function runs.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_UMI_TOOLS_EXTRA = "bio-umi-tools"
_PYSAM_EXTRA = "bio-reads"
_PACKAGE = "umi-tools"
_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9]$")
_CIGAR = re.compile(r"^(?:[0-9]+[MIDNSHP=X])+$")
_MAX_TAGS = 32
_MAX_NAME = 4_096
_MAX_BARCODES = 1_000_000
_MAX_ERRORS = 3
_MAX_POSITION = 2**31 - 1
_RESOLUTIONS = ("discard", "correct")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _text(value: object, name: str, *, limit: int = _MAX_NAME) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} is longer than {limit} characters.")
    return value


def _tag_name(value: object, name: str) -> str:
    if not isinstance(value, str) or not _TAG.match(value):
        raise ConfigurationError(f"{name} must be a two-letter SAM tag name.")
    return value


def _optional_text(value: object, name: str) -> str | None:
    if value is None:
        return None
    return _text(value, name)


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or float(value) < minimum:
        raise ConfigurationError(f"{name} must be a number >= {minimum}.")
    return float(value)


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _decode(value: object) -> str | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _segment(pysam: Any) -> Any:
    record = pysam.AlignedSegment()
    record.query_name = "read"
    return record


def barcode_from_bam_tag(
    tags: Mapping[str, str],
    *,
    umi_tag: str = "RX",
    umi_tag_split: str | None = None,
    umi_tag_delim: str | None = None,
    cell_barcode: bool = False,
    cell_tag: str | None = None,
    cell_tag_split: str | None = "-",
    cell_tag_delim: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract a UMI, and optionally a cell barcode, from the string tags of one alignment.

    ``tags`` maps two-letter SAM tag names to string values. The UMI is read from ``umi_tag``
    (default ``RX``), cut at ``umi_tag_split`` (keeping the first part) and stripped of the
    characters in ``umi_tag_delim``. The cell barcode is read from ``cell_tag`` in the same way,
    with ``cell_tag_split`` defaulting to ``-`` so that a 10x GEM suffix is dropped. A missing tag
    raises ``ConfigurationError``.
    """

    if not isinstance(tags, Mapping) or not 1 <= len(tags) <= _MAX_TAGS:
        raise ConfigurationError(f"tags must be a mapping of 1 to {_MAX_TAGS} SAM tags.")
    for key, value in tags.items():
        _tag_name(key, "tag name")
        if not isinstance(value, str):
            raise ConfigurationError(f"tag {key} must have a string value.")
    umi_name = _tag_name(umi_tag, "umi_tag")
    flag = _flag(cell_barcode, "cell_barcode")
    cell_name = _tag_name(cell_tag, "cell_tag") if flag else None
    required = [umi_name, cell_name] if cell_name else [umi_name]
    missing = [key for key in required if key not in tags]
    if missing:
        raise ConfigurationError(f"tags do not contain {missing}.")
    pysam = require("pysam", _PYSAM_EXTRA)
    sam = require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)
    record = _segment(pysam)
    for key, value in tags.items():
        record.set_tag(key, value, value_type="Z")
    try:
        umi_bytes, cell_bytes = sam.get_barcode_tag(
            record,
            cell_barcode=flag,
            umi_tag=umi_name,
            cell_tag=cell_name,
            umi_tag_split=_optional_text(umi_tag_split, "umi_tag_split"),
            umi_tag_delim=_optional_text(umi_tag_delim, "umi_tag_delim"),
            cell_tag_split=_optional_text(cell_tag_split, "cell_tag_split"),
            cell_tag_delim=_optional_text(cell_tag_delim, "cell_tag_delim"),
        )
    except ValueError as exc:
        raise ConfigurationError(f"could not extract the barcode from the tags: {exc}") from exc
    return {
        "status": "completed",
        "operation": "barcode_from_bam_tag",
        "umi_tag": umi_name,
        "cell_tag": cell_name,
        "umi": _decode(umi_bytes),
        "cell_barcode": _decode(cell_bytes),
        "provenance": _provenance(method="umi_tools.sam_methods.get_barcode_tag"),
    }


def metacontig_contig_map(
    bam_path: str | Path,
    gene_transcript_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map each gene (metacontig) to the references of a local SAM or BAM that belong to it.

    The gene-transcript file holds gene and transcript names in two tab-separated columns. A
    transcript is kept only when the alignment header lists it as a reference. As in umi_tools,
    lines starting with ``#`` are skipped and the first empty line ends the file.
    """

    bam = input_path(bam_path)
    gene_map = input_path(gene_transcript_path)
    pysam = require("pysam", _PYSAM_EXTRA)
    umi_map = require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)
    with pysam.AlignmentFile(str(bam), "r") as alignment:
        mapped = _gene_map(umi_map, alignment, gene_map)
    metacontigs = {gene: sorted(mapped[gene]) for gene in sorted(mapped)}
    return {
        "status": "completed",
        "operation": "metacontig_contig_map",
        "bam_path": str(bam),
        "gene_transcript_path": str(gene_map),
        "metacontig_count": len(metacontigs),
        "metacontigs": metacontigs,
        "provenance": _provenance(method="umi_tools.sam_methods.getMetaContig2contig"),
    }


def _gene_map(umi_map: Any, alignment: Any, gene_map: Path) -> dict[str, set[str]]:
    try:
        return umi_map.getMetaContig2contig(alignment, str(gene_map))
    except ValueError as exc:
        raise ConfigurationError(
            f"each gene-transcript line must hold two tab-separated fields: {exc}"
        ) from exc


def tag_reads_by_metacontig(
    bam_path: str | Path,
    gene_transcript_path: str | Path,
    *,
    metatag: str = "XG",
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the alignments of an indexed BAM in gene (metacontig) order.

    Alignments of every reference that belongs to a gene are fetched and tagged in memory with
    ``metatag``; no file is written. Counts are per gene, so a reference shared by two genes is
    counted for both. The BAM needs an index (``.bai`` or ``.csi``) for fetching by reference.
    """

    tag = _tag_name(metatag, "metatag")
    bam = input_path(bam_path)
    gene_map = input_path(gene_transcript_path)
    pysam = require("pysam", _PYSAM_EXTRA)
    umi_map = require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        if not alignment.has_index():
            raise ConfigurationError(
                f"{bam.name} needs an index (.bai or .csi) to fetch by reference."
            )
        mapping = _gene_map(umi_map, alignment, gene_map)
        ordered = {gene: sorted(mapping[gene]) for gene in sorted(mapping)}
        counts: Counter[str] = Counter()
        for read in umi_map.metafetcher(alignment, ordered, tag):
            counts[str(read.get_tag(tag))] += 1
    per_gene = {gene: counts.get(gene, 0) for gene in sorted(ordered)}
    return {
        "status": "completed",
        "operation": "tag_reads_by_metacontig",
        "bam_path": str(bam),
        "metatag": tag,
        "metacontig_count": len(per_gene),
        "read_count": sum(per_gene.values()),
        "reads_per_metacontig": per_gene,
        "provenance": _provenance(method="umi_tools.sam_methods.metafetcher"),
    }


def soft_clip_adjusted_position(
    cigar: str,
    *,
    reference_start: int = 0,
    is_reverse: bool = False,
    soft_clip_threshold: float = 4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the UMI deduplication position of an alignment from its CIGAR, after soft clipping.

    For a forward read, ``position`` is the start minus any leading soft clip. For a reverse
    read, ``position`` is the end plus any trailing soft clip, and ``start`` is the unclipped
    start. A read counts as spliced when its CIGAR has an ``N``, or when a soft clip longer than
    ``soft_clip_threshold`` sits at the 5' end. ``splice_offset`` is the offset of the first
    junction from the start of the read, or ``None`` when there is no junction.
    """

    if not isinstance(cigar, str) or not _CIGAR.match(cigar):
        raise ConfigurationError("cigar must be a CIGAR string such as 3S5M2N2M.")
    start = _integer(reference_start, "reference_start", minimum=0, maximum=_MAX_POSITION)
    flag = _flag(is_reverse, "is_reverse")
    threshold = _number(soft_clip_threshold, "soft_clip_threshold", minimum=0.0)
    pysam = require("pysam", _PYSAM_EXTRA)
    umi_map = require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)
    record = _segment(pysam)
    record.cigarstring = cigar
    record.reference_start = start
    record.is_reverse = flag
    try:
        read_start, position, spliced = umi_map.get_read_position(record, threshold)
    except (ValueError, IndexError) as exc:
        raise ConfigurationError(f"CIGAR cannot be used for a position: {exc}") from exc
    offset = None if spliced is False else int(spliced)
    return {
        "status": "completed",
        "operation": "soft_clip_adjusted_position",
        "cigar": cigar,
        "reference_start": start,
        "is_reverse": flag,
        "start": int(read_start),
        "position": int(position),
        "is_spliced": offset is not None,
        "splice_offset": offset,
        "soft_clip_threshold": threshold,
        "provenance": _provenance(method="umi_tools.sam_methods.get_read_position"),
    }


def parse_prefixed_umi_name(
    read_name: str,
    *,
    cell_barcode: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a UMI, and optionally a cell barcode, from a name written by umi_tools extract.

    The name is split at ``:``. A field ``UMI_<sequence>`` gives the UMI, and ``CELL_<sequence>``
    gives the cell barcode when ``cell_barcode`` is True. A name without a UMI field raises
    ``ConfigurationError``.
    """

    name = _text(read_name, "read_name")
    flag = _flag(cell_barcode, "cell_barcode")
    umi_map = require("umi_tools.sam_methods", _UMI_TOOLS_EXTRA)
    try:
        umi_bytes, cell_bytes = umi_map.get_barcode_umis(SimpleNamespace(qname=name), flag)
    except ValueError as exc:
        raise ConfigurationError(f"read name has no UMI_ field: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_prefixed_umi_name",
        "read_name": name,
        "cell_barcode_requested": flag,
        "umi": _decode(umi_bytes),
        "cell_barcode": _decode(cell_bytes),
        "provenance": _provenance(method="umi_tools.sam_methods.get_barcode_umis"),
    }


def resolve_above_knee_barcodes(
    cell_barcode_counts: Mapping[str, int],
    cell_whitelist: Sequence[str],
    true_to_false_map: Mapping[str, Sequence[str]],
    *,
    errors: int = 1,
    resolution_method: str = "discard",
    verbose: bool = True,
) -> dict[str, Any]:
    """Discard, or correct, whitelist barcodes that lie within ``errors`` edits of a commoner one.

    Barcodes are visited from the lowest count upward, and a barcode with a near miss among the
    more common barcodes is discarded. With ``resolution_method="correct"``, a single substitution
    near miss is also recorded in ``true_to_false_map`` (the common barcode gains the erroneous
    one). Counts missing from ``cell_barcode_counts`` count as zero. The inputs are not modified.
    """

    if not isinstance(cell_barcode_counts, Mapping) or not isinstance(cell_whitelist, Sequence):
        raise ConfigurationError(
            "cell_barcode_counts must be a mapping and cell_whitelist a sequence."
        )
    if isinstance(cell_whitelist, str) or len(cell_barcode_counts) > _MAX_BARCODES:
        raise ConfigurationError(f"at most {_MAX_BARCODES} barcodes are supported.")
    counts: Counter[str] = Counter()
    for barcode, value in cell_barcode_counts.items():
        counts[_text(barcode, "barcode count key", limit=256)] = _integer(
            value, "barcode count", minimum=0, maximum=2**53
        )
    texts = [_text(item, "cell_whitelist barcode", limit=256) for item in cell_whitelist]
    if len(set(texts)) != len(texts):
        raise ConfigurationError("cell_whitelist must not repeat barcodes.")
    whitelist = sorted(texts)
    mapping: defaultdict[str, set[str]] = defaultdict(set)
    for true_barcode, false_barcodes in true_to_false_map.items():
        if isinstance(false_barcodes, str) or not isinstance(false_barcodes, Sequence):
            raise ConfigurationError("each true_to_false_map value must be a sequence of barcodes.")
        key = _text(true_barcode, "true_to_false_map key", limit=256)
        mapping[key].update(_text(item, "false barcode", limit=256) for item in false_barcodes)
    level = _integer(errors, "errors", minimum=1, maximum=_MAX_ERRORS)
    method = _choice(resolution_method, _RESOLUTIONS, "resolution_method")
    umi_map = require("umi_tools.whitelist_methods", _UMI_TOOLS_EXTRA)
    kept, corrected = umi_map.errorDetectAboveThreshold(
        counts, whitelist, mapping, errors=level, resolution_method=method
    )
    kept_sorted = sorted(kept)
    return {
        "status": "completed",
        "operation": "resolve_above_knee_barcodes",
        "resolution_method": method,
        "errors": level,
        "kept": kept_sorted,
        "discarded": sorted(set(whitelist) - set(kept_sorted)),
        "true_to_false_map": {key: sorted(corrected[key]) for key in sorted(corrected)},
        "provenance": _provenance(method="umi_tools.whitelist_methods.errorDetectAboveThreshold"),
    }


__all__ = [
    "barcode_from_bam_tag",
    "metacontig_contig_map",
    "parse_prefixed_umi_name",
    "resolve_above_knee_barcodes",
    "soft_clip_adjusted_position",
    "tag_reads_by_metacontig",
]

_register_public_verbose(globals(), tuple(__all__))
