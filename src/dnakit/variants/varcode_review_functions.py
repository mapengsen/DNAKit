"""Offline varcode wrappers for nucleotide predicates, reference names, breakends and layouts.

Every function here runs on local inputs. Breakend pairing and VCF rewriting read a local VCF
with ``varcode.load_vcf``; the genome name is used only to infer contig naming, and no reference
sequence or pyensembl cache is read. Effect prediction, RNA evidence and transcript-level fusion
models need pyensembl data and are not wrapped here.

Coordinates are 1-based and inclusive. Breakend sides use ``keeps='left'`` (the bases up to and
including ``position`` are retained) or ``keeps='right'`` (the bases from ``position`` on).

Packages are imported only when a function runs; a missing ``varcode`` raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-varcode]"`` hint.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_VARCODE_EXTRA = "bio-varcode"
_DEFAULT_GENOME = "GRCh38"
_MAX_TEXT = 10_000
_MAX_SEQUENCE = 1_000_000
_MAX_POSITION = 2_000_000_000
_MAX_RECORDS = 100_000
_MAX_SEGMENTS = 1_000
_MAX_INTERVALS = 10_000
_KEEPS = (None, "left", "right")
_STRANDS = ("+", "-")


def _provenance(reference: str = "none: offline inputs only") -> dict[str, str]:
    return {"package": f"varcode {metadata.version('varcode')}", "reference": reference}


def _text(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > _MAX_TEXT
        or any(character.isspace() or character == "\x00" for character in value)
    ):
        raise ConfigurationError(f"{name} must be a non-empty string without whitespace.")
    return value


def _int(value: object, name: str, *, minimum: int = 0, maximum: int = _MAX_POSITION) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _mapping(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping.")
    return dict(value)


def _listed(value: object, name: str, limit: int) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(value) > limit:
        raise ConfigurationError(f"{name} may hold at most {limit} items.")
    return value


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(str(exc)) from exc


def _nucleotide(value: object, name: str) -> str:
    if not isinstance(value, str) or len(value) != 1 or not value.isalpha():
        raise ConfigurationError(f"{name} must be a single letter nucleotide code.")
    return value.upper()


def _sequence(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or len(value) > _MAX_SEQUENCE:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {_MAX_SEQUENCE}.")
    if not value.isalpha():
        raise ConfigurationError(f"{name} must contain letters only.")
    return value.upper()


def _segments(value: object, name: str, module: Any) -> tuple[Any, ...]:
    """Reference-derived layout segments from mappings with contig, start, end and strand."""

    items = _listed(value, name, _MAX_SEGMENTS)
    segments = []
    for index, item in enumerate(items):
        label = f"{name}[{index}]"
        row = _mapping(item, label)
        missing = sorted({"contig", "start", "end"} - row.keys())
        if missing:
            raise ConfigurationError(f"{label} is missing fields: {', '.join(missing)}.")
        contig = _text(row["contig"], f"{label}.contig")
        start = _int(row["start"], f"{label}.start", minimum=1)
        end = _int(row["end"], f"{label}.end", minimum=1)
        if end < start:
            raise ConfigurationError(f"{label}.end must be at least {label}.start.")
        strand = row.get("strand", "+")
        if strand not in _STRANDS:
            raise ConfigurationError(f"{label}.strand must be '+' or '-'.")
        segments.append(_call(module.LayoutSegment, contig, start, end, strand=strand))
    return tuple(segments)


def _breakend(value: object, name: str, module: Any) -> Any:
    row = _mapping(value, name)
    missing = sorted({"contig", "position"} - row.keys())
    if missing:
        raise ConfigurationError(f"{name} is missing fields: {', '.join(missing)}.")
    keeps = row.get("keeps")
    if keeps not in _KEEPS:
        raise ConfigurationError(f"{name}.keeps must be 'left', 'right' or None.")
    contig = _text(row["contig"], f"{name}.contig")
    position = _int(row["position"], f"{name}.position", minimum=1)
    return module.Breakend(contig, position, keeps)


def _segment_rows(segments: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "contig": str(segment.contig),
            "start": int(segment.start),
            "end": int(segment.end),
            "strand": str(segment.strand),
            "length": int(segment.end) - int(segment.start) + 1,
        }
        for segment in segments
    ]


def _load_collection(source: Path, genome: str, **options: Any) -> Any:
    module = require("varcode", _VARCODE_EXTRA)
    collection = _call(
        module.load_vcf,
        str(source),
        genome=genome,
        include_info=True,
        distinct=False,
        max_variants=_MAX_RECORDS + 1,
        **options,
    )
    if len(collection) > _MAX_RECORDS:
        raise ConfigurationError(f"VCF has more than {_MAX_RECORDS} records.")
    return collection


def _fresh_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"{name} already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {target.parent}")
    return target


def is_purine_base(
    nucleotide: str,
    *,
    allow_extended: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Whether one nucleotide code is a purine (A or G); other codes need ``allow_extended``."""

    base = _nucleotide(nucleotide, "nucleotide")
    if not isinstance(allow_extended, bool):
        raise ConfigurationError("allow_extended must be True or False.")
    module = require("varcode.nucleotides", _VARCODE_EXTRA)
    answer = _call(module.is_purine, base, allow_extended_nucleotides=allow_extended)
    return {
        "status": "completed",
        "operation": "is_purine_base",
        "nucleotide": base,
        "allow_extended": allow_extended,
        "is_purine": bool(answer),
        "provenance": _provenance(),
    }


def has_only_standard_nucleotides(nucleotides: str, *, verbose: bool = True) -> dict[str, Any]:
    """Whether every base of a sequence is A, C, G or T, with the 1-based positions of the rest."""

    text = _sequence(nucleotides, "nucleotides")
    module = require("varcode.nucleotides", _VARCODE_EXTRA)
    standard = set(module.STANDARD_NUCLEOTIDES)
    positions = [index for index, base in enumerate(text, start=1) if base not in standard]
    return {
        "status": "completed",
        "operation": "has_only_standard_nucleotides",
        "length": len(text),
        "standard": bool(module.all_standard_nucleotides(text)),
        "non_standard_positions": positions,
        "provenance": _provenance(),
    }


def is_mitochondrial_contig(contig: str | None, *, verbose: bool = True) -> dict[str, Any]:
    """Whether a contig name (chrM, MT, mtDNA and similar) refers to the mitochondrial genome."""

    if contig is not None:
        _text(contig, "contig")
    module = require("varcode.effects.codon_tables", _VARCODE_EXTRA)
    return {
        "status": "completed",
        "operation": "is_mitochondrial_contig",
        "contig": contig,
        "mitochondrial": bool(module.is_mitochondrial_contig(contig)),
        "provenance": _provenance(),
    }


def is_ucsc_reference_name(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Whether a genome name is one of varcode's UCSC assembly names (for example hg19 or mm10)."""

    text = _text(name, "name")
    module = require("varcode.reference", _VARCODE_EXTRA)
    return {
        "status": "completed",
        "operation": "is_ucsc_reference_name",
        "name": text,
        "ucsc": bool(module.is_ucsc_reference_name(text)),
        "provenance": _provenance(),
    }


def infer_reference_assembly(reference: str, *, verbose: bool = True) -> dict[str, Any]:
    """Infer the canonical assembly name from a reference name or FASTA path, without reading it."""

    text = _text(reference, "reference")
    module = require("varcode.reference", _VARCODE_EXTRA)
    assembly = _call(module.infer_reference_name, text)
    return {
        "status": "completed",
        "operation": "infer_reference_assembly",
        "reference": text,
        "assembly": str(assembly),
        "provenance": _provenance("assembly names are matched against varcode's table only"),
    }


def breakend_junction_span(
    near: Mapping[str, Any],
    far: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """The DEL, DUP or INV span one same-contig novel junction describes, or None."""

    structural = require("varcode.structural_variant", _VARCODE_EXTRA)
    first = _breakend(near, "near", structural)
    second = _breakend(far, "far", structural)
    details = _call(structural.junction_span_details, (first, second))
    span = None
    if details is not None:
        sv_type, start, end, affected_start, affected_end = details
        span = {
            "sv_type": str(sv_type),
            "start": int(start),
            "end": int(end),
            "affected_start": int(affected_start),
            "affected_end": int(affected_end),
        }
    return {
        "status": "completed",
        "operation": "breakend_junction_span",
        "span": span,
        "provenance": _provenance(),
    }


def join_breakend_layouts(
    first_segments: Sequence[Mapping[str, Any]],
    first_breakend: Mapping[str, Any],
    second_segments: Sequence[Mapping[str, Any]],
    second_breakend: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join two retained breakend sides into one oriented layout with an explicit novel junction."""

    module = require("varcode.genomic_layout", _VARCODE_EXTRA)
    structural = require("varcode.structural_variant", _VARCODE_EXTRA)
    left = module.GenomicLayout(_segments(first_segments, "first_segments", module))
    right = module.GenomicLayout(_segments(second_segments, "second_segments", module))
    breakend_a = _breakend(first_breakend, "first_breakend", structural)
    breakend_b = _breakend(second_breakend, "second_breakend", structural)
    if breakend_a.keeps is None or breakend_b.keeps is None:
        raise ConfigurationError("Both breakends must state which side they keep.")
    joined = _call(module.join_breakends, left, breakend_a, right, breakend_b)
    return {
        "status": "completed",
        "operation": "join_breakend_layouts",
        "segments": _segment_rows(joined.segments),
        "length": int(joined.length),
        "provenance": _provenance(),
    }


def project_intervals_onto_layout(
    segments: Sequence[Mapping[str, Any]],
    intervals: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project reference intervals onto a layout; each maximal matching run becomes one piece."""

    module = require("varcode.genomic_layout", _VARCODE_EXTRA)
    layout = module.GenomicLayout(_segments(segments, "segments", module))
    items = _listed(intervals, "intervals", _MAX_INTERVALS)
    projected = []
    for index, item in enumerate(items):
        label = f"intervals[{index}]"
        row = _mapping(item, label)
        missing = sorted({"contig", "start", "end"} - row.keys())
        if missing:
            raise ConfigurationError(f"{label} is missing fields: {', '.join(missing)}.")
        start = _int(row["start"], f"{label}.start", minimum=1)
        end = _int(row["end"], f"{label}.end", minimum=1)
        if end < start:
            raise ConfigurationError(f"{label}.end must be at least {label}.start.")
        projected.append((_text(row["contig"], f"{label}.contig"), start, end))
    pieces = _call(module.project_intervals, layout, projected)
    rows = [
        {"segments": _segment_rows(piece.segments), "length": int(piece.length)} for piece in pieces
    ]
    rows.sort(key=lambda row: (row["segments"][0]["contig"], row["segments"][0]["start"]))
    return {
        "status": "completed",
        "operation": "project_intervals_onto_layout",
        "count": len(rows),
        "pieces": rows,
        "provenance": _provenance(),
    }


def pair_breakend_mates(
    path: str,
    *,
    genome: str = _DEFAULT_GENOME,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge MATEID-paired BND records of a local VCF into one combined record per event."""

    source = input_path(path)
    genome_name = _text(genome, "genome")
    transforms = require("varcode.transforms", _VARCODE_EXTRA)
    collection = _load_collection(source, genome_name, parse_structural_variants=True)
    paired = _call(transforms.pair_breakends, collection)
    rows = []
    for variant in paired:
        sources = tuple(getattr(variant, "source_variants", ()) or ())
        # A combined record does not keep the source contig spelling, so read it from a half.
        origin = sources[0] if sources else variant
        rows.append(
            {
                "contig": str(origin.original_contig),
                "start": int(variant.original_start),
                "ref": str(variant.original_ref),
                "alt": str(variant.alt),
                "source_starts": sorted(int(item.original_start) for item in sources),
            }
        )
    rows.sort(key=lambda row: (row["contig"], row["start"], row["alt"]))
    return {
        "status": "completed",
        "operation": "pair_breakend_mates",
        "path": str(source),
        "input_count": len(collection),
        "output_count": len(rows),
        "paired_count": sum(1 for row in rows if row["source_starts"]),
        "records": rows,
        "provenance": _provenance(f"genome name {genome_name} is used only to infer contigs"),
    }


def rewrite_vcf_records(
    path: str,
    output: str,
    *,
    genome: str = _DEFAULT_GENOME,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the variants of a local VCF back out as a new VCF using varcode rules."""

    source = input_path(path)
    genome_name = _text(genome, "genome")
    target = _fresh_file(output, "output")
    writer = require("varcode.vcf_output", _VARCODE_EXTRA)
    collection = _load_collection(source, genome_name)
    handle = target.open("x", encoding="utf-8")
    try:
        with handle:
            _call(writer.variants_to_vcf, list(collection), collection.metadata, out=handle)
    except Exception:
        target.unlink(missing_ok=True)
        raise
    return {
        "status": "completed",
        "operation": "rewrite_vcf_records",
        "source": str(source),
        "output": str(target),
        "variants": len(collection),
        "provenance": _provenance(f"genome name {genome_name} is used only to infer contigs"),
    }


__all__ = [
    "breakend_junction_span",
    "has_only_standard_nucleotides",
    "infer_reference_assembly",
    "is_mitochondrial_contig",
    "is_purine_base",
    "is_ucsc_reference_name",
    "join_breakend_layouts",
    "pair_breakend_mates",
    "project_intervals_onto_layout",
    "rewrite_vcf_records",
]

_register_public_verbose(globals(), tuple(__all__))
