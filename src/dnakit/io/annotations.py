"""Strict, bounded annotation codecs for GFF3, BED3-BED6, and AGP 2.1."""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal, TextIO, TypeAlias, cast
from urllib.parse import quote, unquote_to_bytes

from dnakit.core._json import to_json_compatible
from dnakit.core.coordinates import Interval
from dnakit.core.enums import GapKind, Strand
from dnakit.core.feature import DNAFeature
from dnakit.core.gap import Gap
from dnakit.exceptions import ConfigurationError, InputFormatError

from ._advanced_common import checked_lines, open_text_source, write_text_path

AnnotationSource: TypeAlias = str | os.PathLike[str] | TextIO
AnnotationTarget: TypeAlias = str | os.PathLike[str] | TextIO
AnnotationFormat = Literal["gff3", "bed"]


def _validate_limit(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code="INVALID_ANNOTATION_LIMIT",
            context={"field": name, "value": value},
        )


def _bounded_append(values: list[object], value: object, *, max_records: int) -> None:
    if len(values) >= max_records:
        raise InputFormatError(
            "Annotation input exceeds max_records.",
            code="ANNOTATION_RECORD_LIMIT_EXCEEDED",
            context={"max_records": max_records},
        )
    values.append(value)


def _bounded_header(values: list[str], value: str, *, max_records: int) -> None:
    if len(values) >= max_records:
        raise InputFormatError(
            "Annotation headers exceed the configured resource limit.",
            code="ANNOTATION_HEADER_LIMIT_EXCEEDED",
            context={"max_header_lines": max_records},
        )
    values.append(value)


def _strand_from_symbol(value: str, *, line_number: int, allow_unknown: bool = True) -> Strand:
    mapping = {"+": Strand.FORWARD, "-": Strand.REVERSE, ".": Strand.UNKNOWN}
    if allow_unknown:
        mapping["?"] = Strand.UNKNOWN
    try:
        return mapping[value]
    except KeyError as exc:
        raise InputFormatError(
            "Invalid annotation strand symbol.",
            code="INVALID_ANNOTATION_STRAND",
            context={"line_number": line_number, "strand": value},
        ) from exc


def _strand_symbol(strand: Strand) -> str:
    if strand is Strand.FORWARD:
        return "+"
    if strand is Strand.REVERSE:
        return "-"
    if strand is Strand.UNKNOWN:
        return "."
    raise InputFormatError(
        "The annotation codec cannot serialize a both-strands feature.",
        code="UNSUPPORTED_ANNOTATION_STRAND",
    )


@dataclass(frozen=True, slots=True)
class AnnotationEntry:
    """A sequence identifier paired with one immutable DNA feature."""

    sequence_id: str
    feature: DNAFeature

    def __post_init__(self) -> None:
        if (
            not isinstance(self.sequence_id, str)
            or not self.sequence_id
            or any(symbol.isspace() for symbol in self.sequence_id)
        ):
            raise ConfigurationError(
                "AnnotationEntry sequence_id must be non-empty text.",
                code="INVALID_ANNOTATION_ENTRY",
            )
        if not isinstance(self.feature, DNAFeature):
            raise ConfigurationError(
                "AnnotationEntry feature must be a DNAFeature.",
                code="INVALID_ANNOTATION_ENTRY",
            )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class AnnotationDocument:
    """Bounded annotation document retaining safe header directives."""

    format: AnnotationFormat
    entries: tuple[AnnotationEntry, ...]
    headers: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.entries, tuple) or not isinstance(self.headers, tuple):
            raise ConfigurationError(
                "Annotation document entries and headers must be tuples.",
                code="INVALID_ANNOTATION_DOCUMENT",
            )
        if self.format not in {"gff3", "bed"}:
            raise ConfigurationError("Unknown annotation format.", code="INVALID_ANNOTATION_FORMAT")
        if any(not isinstance(item, AnnotationEntry) for item in self.entries):
            raise ConfigurationError(
                "entries must contain AnnotationEntry objects.", code="INVALID_ANNOTATION_DOCUMENT"
            )
        if any(not isinstance(item, str) or "\n" in item or "\r" in item for item in self.headers):
            raise ConfigurationError(
                "Annotation headers must be single-line strings.",
                code="INVALID_ANNOTATION_DOCUMENT",
            )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def _decode_gff_value(value: str, *, line_number: int) -> str:
    for match in re.finditer("%", value):
        if not re.fullmatch(r"[0-9A-Fa-f]{2}", value[match.start() + 1 : match.start() + 3]):
            raise InputFormatError(
                "GFF3 attribute contains malformed percent encoding.",
                code="INVALID_GFF3_ATTRIBUTE",
                context={"line_number": line_number},
            )
    try:
        return unquote_to_bytes(value).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "GFF3 attribute percent encoding is not valid UTF-8.",
            code="INVALID_GFF3_ATTRIBUTE",
            context={"line_number": line_number},
        ) from exc


def _parse_gff_attributes(value: str, *, line_number: int) -> dict[str, object]:
    if value == ".":
        return {}
    result: dict[str, object] = {}
    for item in value.split(";"):
        if not item or "=" not in item:
            raise InputFormatError(
                "GFF3 attributes must be semicolon-separated key=value pairs.",
                code="INVALID_GFF3_ATTRIBUTE",
                context={"line_number": line_number},
            )
        raw_key, raw_values = item.split("=", 1)
        key = _decode_gff_value(raw_key, line_number=line_number)
        if not key or key in result:
            raise InputFormatError(
                "GFF3 attribute names must be non-empty and unique.",
                code="INVALID_GFF3_ATTRIBUTE",
                context={"line_number": line_number, "attribute": key},
            )
        values = tuple(
            _decode_gff_value(raw, line_number=line_number) for raw in raw_values.split(",")
        )
        if any(value == "" for value in values):
            raise InputFormatError(
                "GFF3 attribute values must be non-empty.",
                code="INVALID_GFF3_ATTRIBUTE",
                context={"line_number": line_number, "attribute": key},
            )
        result[key] = values[0] if len(values) == 1 else values
    return result


def read_gff3(
    source: AnnotationSource,
    *,
    max_records: int = 1_000_000,
    max_line_length: int = 1_000_000,
    max_header_lines: int = 10_000,
) -> AnnotationDocument:
    """Read GFF3 features using 1-based closed to 0-based half-open conversion."""

    _validate_limit(max_records, "max_records")
    _validate_limit(max_line_length, "max_line_length")
    _validate_limit(max_header_lines, "max_header_lines")
    entries: list[object] = []
    headers: list[str] = []
    seen_version = False
    with open_text_source(source) as handle:
        for line_number, line in checked_lines(handle, max_line_length=max_line_length):
            if not line:
                continue
            if line == "##FASTA":
                raise InputFormatError(
                    "Embedded FASTA is outside the DNAKit GFF3 annotation subset.",
                    code="UNSUPPORTED_GFF3_FASTA",
                    context={"line_number": line_number},
                )
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    version_match = re.fullmatch(r"##gff-version 3(?:\.\d+){0,2}", line)
                    if version_match is None or seen_version or line_number != 1:
                        raise InputFormatError(
                            "GFF3 requires one supported version directive.",
                            code="INVALID_GFF3_VERSION",
                            context={"line_number": line_number},
                        )
                    seen_version = True
                _bounded_header(headers, line, max_records=max_header_lines)
                continue
            if not seen_version:
                raise InputFormatError(
                    "GFF3 data appears before '##gff-version 3'.",
                    code="GFF3_VERSION_REQUIRED",
                    context={"line_number": line_number},
                )
            columns = line.split("\t")
            if len(columns) != 9:
                raise InputFormatError(
                    "GFF3 records must contain exactly nine tab-separated columns.",
                    code="INVALID_GFF3_COLUMNS",
                    context={"line_number": line_number, "column_count": len(columns)},
                )
            (
                seqid,
                source_name,
                feature_type,
                start_raw,
                end_raw,
                score_raw,
                strand_raw,
                phase_raw,
                attributes_raw,
            ) = columns
            if not seqid or seqid == "." or not feature_type or feature_type == ".":
                raise InputFormatError(
                    "GFF3 seqid and type must be present.",
                    code="INVALID_GFF3_FIELD",
                    context={"line_number": line_number},
                )
            if not source_name:
                raise InputFormatError(
                    "GFF3 source must be non-empty or '.'.",
                    code="INVALID_GFF3_FIELD",
                    context={"line_number": line_number},
                )
            try:
                start = int(start_raw)
                end = int(end_raw)
            except ValueError as exc:
                raise InputFormatError(
                    "GFF3 coordinates must be integers.",
                    code="INVALID_GFF3_COORDINATE",
                    context={"line_number": line_number},
                ) from exc
            if start < 1 or end < start:
                raise InputFormatError(
                    "GFF3 coordinates must form a positive 1-based closed interval.",
                    code="INVALID_GFF3_COORDINATE",
                    context={"line_number": line_number},
                )
            try:
                score = None if score_raw == "." else float(score_raw)
            except ValueError as exc:
                raise InputFormatError(
                    "GFF3 score must be numeric or '.'.",
                    code="INVALID_GFF3_SCORE",
                    context={"line_number": line_number},
                ) from exc
            if score is not None and not math.isfinite(score):
                raise InputFormatError(
                    "GFF3 score must be finite.",
                    code="INVALID_GFF3_SCORE",
                    context={"line_number": line_number},
                )
            if phase_raw not in {".", "0", "1", "2"}:
                raise InputFormatError(
                    "GFF3 phase must be 0, 1, 2, or '.'.",
                    code="INVALID_GFF3_PHASE",
                    context={"line_number": line_number},
                )
            if (feature_type == "CDS" and phase_raw == ".") or (
                feature_type != "CDS" and phase_raw != "."
            ):
                raise InputFormatError(
                    "GFF3 phase is required only for CDS features.",
                    code="INVALID_GFF3_PHASE",
                    context={"line_number": line_number, "feature_type": feature_type},
                )
            attributes = _parse_gff_attributes(attributes_raw, line_number=line_number)
            identifier = attributes.pop("ID", None)
            name = attributes.pop("Name", None)
            if identifier is not None and not isinstance(identifier, str):
                raise InputFormatError(
                    "GFF3 ID must have exactly one value.",
                    code="INVALID_GFF3_ATTRIBUTE",
                    context={"line_number": line_number},
                )
            if name is not None and not isinstance(name, str):
                raise InputFormatError(
                    "GFF3 Name must have exactly one value.",
                    code="INVALID_GFF3_ATTRIBUTE",
                    context={"line_number": line_number},
                )
            feature = DNAFeature(
                feature_type,
                Interval(start - 1, end),
                id=identifier if isinstance(identifier, str) else None,
                label=name if isinstance(name, str) else None,
                strand=_strand_from_symbol(strand_raw, line_number=line_number),
                score=score,
                phase=None if phase_raw == "." else int(phase_raw),
                qualifiers=attributes,
                source=None if source_name == "." else source_name,
            )
            _bounded_append(entries, AnnotationEntry(seqid, feature), max_records=max_records)
    if not seen_version:
        raise InputFormatError(
            "GFF3 input is missing '##gff-version 3'.", code="GFF3_VERSION_REQUIRED"
        )
    return AnnotationDocument(
        "gff3", cast(tuple[AnnotationEntry, ...], tuple(entries)), tuple(headers)
    )


_GFF_SAFE = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:^*$@!+_?-|"


def _encode_gff(value: str) -> str:
    return quote(value, safe=_GFF_SAFE)


def _gff_attributes(feature: DNAFeature) -> str:
    values: list[tuple[str, object]] = []
    if feature.id is not None:
        values.append(("ID", feature.id))
    if feature.label is not None:
        values.append(("Name", feature.label))
    if "ID" in feature.qualifiers or "Name" in feature.qualifiers:
        raise InputFormatError(
            "Feature qualifiers conflict with reserved GFF3 ID or Name attributes.",
            code="RESERVED_GFF3_ATTRIBUTE",
        )
    values.extend(feature.qualifiers.items())
    rendered: list[str] = []
    for key, value in values:
        if not isinstance(key, str) or not key:
            raise InputFormatError(
                "GFF3 attribute names must be non-empty strings.", code="INVALID_GFF3_ATTRIBUTE"
            )
        sequence = value if isinstance(value, tuple) else (value,)
        if not sequence or any(not isinstance(item, (str, int, float, bool)) for item in sequence):
            raise InputFormatError(
                "GFF3 attributes must contain non-empty scalar values.",
                code="UNSUPPORTED_GFF3_ATTRIBUTE",
            )
        rendered.append(
            f"{_encode_gff(key)}="
            + ",".join(
                _encode_gff(str(item).lower() if isinstance(item, bool) else str(item))
                for item in sequence
            )
        )
    return ";".join(rendered) if rendered else "."


def write_gff3(
    document: AnnotationDocument | Iterable[AnnotationEntry],
    target: AnnotationTarget,
    *,
    overwrite: bool = False,
    create_parents: bool = False,
) -> int:
    """Write deterministic GFF3 with internal coordinate conversion."""

    entries: Iterable[AnnotationEntry]
    if isinstance(document, AnnotationDocument):
        if document.format != "gff3":
            raise ConfigurationError(
                "AnnotationDocument format must be gff3.", code="INVALID_ANNOTATION_FORMAT"
            )
        entries = document.entries
        extra_headers = tuple(
            header for header in document.headers if not header.startswith("##gff-version")
        )
    else:
        entries = document
        extra_headers = ()

    def writer(handle: TextIO) -> int:
        handle.write("##gff-version 3\n")
        for header in extra_headers:
            if not header.startswith("#") or header == "##FASTA":
                raise InputFormatError(
                    "Unsafe or unsupported GFF3 header.", code="INVALID_GFF3_HEADER"
                )
            handle.write(header + "\n")
        count = 0
        for entry in entries:
            if not isinstance(entry, AnnotationEntry):
                raise TypeError("GFF3 entries must be AnnotationEntry objects.")
            feature = entry.feature
            if not isinstance(feature.location, Interval) or len(feature.location) == 0:
                raise InputFormatError(
                    "GFF3 subset output requires one non-empty interval per feature.",
                    code="UNSUPPORTED_GFF3_LOCATION",
                )
            source_name = feature.source or "."
            score = "." if feature.score is None else format(feature.score, ".17g")
            phase = "." if feature.phase is None else str(feature.phase)
            if (feature.type == "CDS" and feature.phase is None) or (
                feature.type != "CDS" and feature.phase is not None
            ):
                raise InputFormatError(
                    "GFF3 phase is required only for CDS features.",
                    code="INVALID_GFF3_PHASE",
                    context={"feature_type": feature.type},
                )
            fields = (
                entry.sequence_id,
                source_name,
                feature.type,
                str(feature.location.start + 1),
                str(feature.location.end),
                score,
                _strand_symbol(feature.strand),
                phase,
                _gff_attributes(feature),
            )
            if any("\t" in value or "\n" in value or "\r" in value for value in fields):
                raise InputFormatError(
                    "GFF3 fields must not contain tabs or line breaks.", code="INVALID_GFF3_FIELD"
                )
            handle.write("\t".join(fields) + "\n")
            count += 1
        return count

    return write_text_path(target, writer, overwrite=overwrite, create_parents=create_parents)


def read_bed(
    source: AnnotationSource,
    *,
    max_records: int = 1_000_000,
    max_line_length: int = 1_000_000,
    max_header_lines: int = 10_000,
) -> AnnotationDocument:
    """Read the BED3-BED6 subset into immutable features."""

    _validate_limit(max_records, "max_records")
    _validate_limit(max_line_length, "max_line_length")
    _validate_limit(max_header_lines, "max_header_lines")
    entries: list[object] = []
    headers: list[str] = []
    with open_text_source(source) as handle:
        for line_number, line in checked_lines(handle, max_line_length=max_line_length):
            if (
                not line
                or line.startswith("#")
                or line.startswith("track ")
                or line.startswith("browser ")
            ):
                if line:
                    _bounded_header(headers, line, max_records=max_header_lines)
                continue
            columns = line.split("\t")
            if not 3 <= len(columns) <= 6:
                raise InputFormatError(
                    "The DNAKit BED subset accepts exactly 3 through 6 columns.",
                    code="UNSUPPORTED_BED_COLUMNS" if len(columns) > 6 else "INVALID_BED_COLUMNS",
                    context={"line_number": line_number, "column_count": len(columns)},
                )
            chrom, start_raw, end_raw = columns[:3]
            if not chrom:
                raise InputFormatError(
                    "BED chromosome must be non-empty.",
                    code="INVALID_BED_FIELD",
                    context={"line_number": line_number},
                )
            try:
                start, end = int(start_raw), int(end_raw)
            except ValueError as exc:
                raise InputFormatError(
                    "BED coordinates must be integers.",
                    code="INVALID_BED_COORDINATE",
                    context={"line_number": line_number},
                ) from exc
            if start < 0 or end <= start:
                raise InputFormatError(
                    "BED requires a non-empty 0-based half-open interval.",
                    code="INVALID_BED_COORDINATE",
                    context={"line_number": line_number},
                )
            name = columns[3] if len(columns) >= 4 and columns[3] != "." else None
            if name == "":
                raise InputFormatError(
                    "BED name must be non-empty or '.'.",
                    code="INVALID_BED_FIELD",
                    context={"line_number": line_number},
                )
            score: float | None = None
            if len(columns) >= 5 and columns[4] != ".":
                try:
                    raw_score = int(columns[4])
                except ValueError as exc:
                    raise InputFormatError(
                        "BED score must be an integer from 0 through 1000.",
                        code="INVALID_BED_SCORE",
                        context={"line_number": line_number},
                    ) from exc
                if not 0 <= raw_score <= 1000:
                    raise InputFormatError(
                        "BED score must be in [0, 1000].",
                        code="INVALID_BED_SCORE",
                        context={"line_number": line_number},
                    )
                score = float(raw_score)
            strand = (
                Strand.UNKNOWN
                if len(columns) < 6
                else _strand_from_symbol(columns[5], line_number=line_number, allow_unknown=False)
            )
            feature = DNAFeature(
                "region",
                Interval(start, end),
                id=name,
                label=name,
                score=score,
                strand=strand,
                source="BED",
            )
            _bounded_append(entries, AnnotationEntry(chrom, feature), max_records=max_records)
    return AnnotationDocument(
        "bed", cast(tuple[AnnotationEntry, ...], tuple(entries)), tuple(headers)
    )


def write_bed(
    document: AnnotationDocument | Iterable[AnnotationEntry],
    target: AnnotationTarget,
    *,
    overwrite: bool = False,
    create_parents: bool = False,
) -> int:
    """Write the deterministic BED6 subset."""

    entries: Iterable[AnnotationEntry]
    headers: tuple[str, ...]
    if isinstance(document, AnnotationDocument):
        if document.format != "bed":
            raise ConfigurationError(
                "AnnotationDocument format must be bed.", code="INVALID_ANNOTATION_FORMAT"
            )
        entries = document.entries
        headers = document.headers
    else:
        entries, headers = document, ()

    def writer(handle: TextIO) -> int:
        for header in headers:
            if not (
                header.startswith("#")
                or header.startswith("track ")
                or header.startswith("browser ")
            ):
                raise InputFormatError("Unsafe BED header.", code="INVALID_BED_HEADER")
            handle.write(header + "\n")
        count = 0
        for entry in entries:
            if not isinstance(entry, AnnotationEntry):
                raise TypeError("BED entries must be AnnotationEntry objects.")
            feature = entry.feature
            if not isinstance(feature.location, Interval) or len(feature.location) == 0:
                raise InputFormatError(
                    "BED6 output requires one non-empty interval.", code="UNSUPPORTED_BED_LOCATION"
                )
            if (
                feature.type != "region"
                or feature.phase is not None
                or feature.qualifiers
                or feature.source not in {None, "BED"}
                or (
                    feature.id is not None
                    and feature.label is not None
                    and feature.id != feature.label
                )
            ):
                raise InputFormatError(
                    "The BED3-BED6 subset would lose feature semantics.",
                    code="BED_FEATURE_LOSS_NOT_ALLOWED",
                )
            name = feature.id or feature.label
            if feature.score is None:
                score: str | None = None
            elif not feature.score.is_integer() or not 0 <= feature.score <= 1000:
                raise InputFormatError(
                    "BED score must be an integer from 0 through 1000.", code="INVALID_BED_SCORE"
                )
            else:
                score = str(int(feature.score))
            fields = [
                entry.sequence_id,
                str(feature.location.start),
                str(feature.location.end),
            ]
            if name is not None or score is not None or feature.strand is not Strand.UNKNOWN:
                fields.append(name or ".")
            if score is not None or feature.strand is not Strand.UNKNOWN:
                fields.append(score or "0")
            if feature.strand is not Strand.UNKNOWN:
                fields.append(_strand_symbol(feature.strand))
            if any("\t" in item or "\n" in item or "\r" in item for item in fields):
                raise InputFormatError(
                    "BED fields must not contain tabs or line breaks.", code="INVALID_BED_FIELD"
                )
            handle.write("\t".join(fields) + "\n")
            count += 1
        return count

    return write_text_path(target, writer, overwrite=overwrite, create_parents=create_parents)


_COMPONENT_TYPES = frozenset("ADFGOPW")
_GAP_TYPES = frozenset("NU")
_GAP_KIND_BY_AGP = {
    "scaffold": GapKind.SCAFFOLD,
    "contig": GapKind.CONTIG,
    "centromere": GapKind.CENTROMERE,
    "short_arm": GapKind.SHORT_ARM,
    "heterochromatin": GapKind.HETEROCHROMATIN,
    "telomere": GapKind.TELOMERE,
    "repeat": GapKind.REPEAT,
    "contamination": GapKind.CONTAMINATION,
}
_AGP_GAP_TYPES = frozenset(_GAP_KIND_BY_AGP)
_AGP_LINKAGE_EVIDENCE = frozenset(
    {
        "paired-ends",
        "align_genus",
        "align_xgenus",
        "align_trnscpt",
        "within_clone",
        "clone_contig",
        "map",
        "proximity_ligation",
        "pcr",
        "strobe",
        "unspecified",
    }
)


@dataclass(frozen=True, slots=True)
class AGPComponent:
    object_id: str
    object_interval: Interval
    part_number: int
    component_type: str
    component_id: str
    component_interval: Interval
    orientation: str

    def __post_init__(self) -> None:
        if (
            not isinstance(self.object_id, str)
            or not self.object_id
            or any(symbol.isspace() for symbol in self.object_id)
        ):
            raise ConfigurationError("AGP object_id must be non-empty.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.object_interval, Interval) or len(self.object_interval) == 0:
            raise ConfigurationError(
                "AGP object interval must be non-empty.", code="INVALID_AGP_ENTRY"
            )
        if (
            isinstance(self.part_number, bool)
            or not isinstance(self.part_number, int)
            or self.part_number < 1
        ):
            raise ConfigurationError("AGP part_number must be positive.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.component_type, str) or self.component_type not in _COMPONENT_TYPES:
            raise ConfigurationError("Invalid AGP component type.", code="INVALID_AGP_ENTRY")
        if (
            not isinstance(self.component_id, str)
            or not self.component_id
            or any(symbol.isspace() for symbol in self.component_id)
        ):
            raise ConfigurationError(
                "AGP component_id must be non-empty.", code="INVALID_AGP_ENTRY"
            )
        if not isinstance(self.component_interval, Interval) or len(self.component_interval) == 0:
            raise ConfigurationError(
                "AGP component interval must be non-empty.", code="INVALID_AGP_ENTRY"
            )
        if len(self.component_interval) != len(self.object_interval):
            raise ConfigurationError(
                "AGP component and object spans must match.", code="INVALID_AGP_ENTRY"
            )
        if not isinstance(self.orientation, str) or self.orientation not in {
            "+",
            "-",
            "?",
            "0",
            "na",
        }:
            raise ConfigurationError("Invalid AGP orientation.", code="INVALID_AGP_ENTRY")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class AGPGap:
    object_id: str
    object_interval: Interval
    part_number: int
    component_type: str
    gap: Gap
    gap_type: str
    linkage: bool
    linkage_evidence: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.object_id, str)
            or not self.object_id
            or any(symbol.isspace() for symbol in self.object_id)
        ):
            raise ConfigurationError("AGP object_id must be non-empty.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.object_interval, Interval) or len(self.object_interval) == 0:
            raise ConfigurationError(
                "AGP object interval must be non-empty.", code="INVALID_AGP_ENTRY"
            )
        if (
            isinstance(self.part_number, bool)
            or not isinstance(self.part_number, int)
            or self.part_number < 1
        ):
            raise ConfigurationError("AGP part_number must be positive.", code="INVALID_AGP_ENTRY")
        if (
            not isinstance(self.component_type, str)
            or self.component_type not in _GAP_TYPES
            or not isinstance(self.gap, Gap)
        ):
            raise ConfigurationError("Invalid AGP gap entry.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.gap_type, str) or self.gap_type not in _AGP_GAP_TYPES:
            raise ConfigurationError("Invalid AGP gap type.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.linkage, bool):
            raise ConfigurationError("AGP linkage must be a boolean.", code="INVALID_AGP_ENTRY")
        if not isinstance(self.linkage_evidence, tuple) or any(
            not isinstance(item, str) or item not in _AGP_LINKAGE_EVIDENCE
            for item in self.linkage_evidence
        ):
            raise ConfigurationError("Invalid AGP linkage evidence.", code="INVALID_AGP_ENTRY")
        if not self.linkage and self.linkage_evidence:
            raise ConfigurationError(
                "Unlinked AGP gaps cannot carry linkage evidence.", code="INVALID_AGP_ENTRY"
            )
        if self.linkage and not self.linkage_evidence:
            raise ConfigurationError(
                "Linked AGP gaps require linkage evidence.", code="INVALID_AGP_ENTRY"
            )
        if self.component_type == "N" and self.gap.length != len(self.object_interval):
            raise ConfigurationError(
                "Known AGP gap span and length must match.", code="INVALID_AGP_ENTRY"
            )
        if self.component_type == "U" and self.gap.length is not None:
            raise ConfigurationError(
                "Unknown AGP gaps must use Gap(length=None).", code="INVALID_AGP_ENTRY"
            )
        if self.component_type == "U" and len(self.object_interval) != 100:
            raise ConfigurationError(
                "Unknown AGP gaps must occupy exactly 100 bases.",
                code="INVALID_AGP_ENTRY",
            )
        if self.gap.crossable is not self.linkage or self.gap.evidence != self.linkage_evidence:
            raise ConfigurationError(
                "AGP Gap linkage fields must agree with the embedded Gap.",
                code="INVALID_AGP_ENTRY",
            )
        if self.gap.kind is not _GAP_KIND_BY_AGP[self.gap_type]:
            raise ConfigurationError(
                "AGP gap_type must agree with the embedded Gap kind.",
                code="INVALID_AGP_ENTRY",
            )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


AGPEntry = AGPComponent | AGPGap


@dataclass(frozen=True, slots=True)
class AGPDocument:
    entries: tuple[AGPEntry, ...]
    headers: tuple[str, ...] = ("##agp-version 2.1",)

    def __post_init__(self) -> None:
        if not isinstance(self.entries, tuple) or not isinstance(self.headers, tuple):
            raise ConfigurationError(
                "AGP document entries and headers must be tuples.",
                code="INVALID_AGP_DOCUMENT",
            )
        if any(not isinstance(entry, (AGPComponent, AGPGap)) for entry in self.entries):
            raise ConfigurationError("Invalid AGP document entry.", code="INVALID_AGP_DOCUMENT")
        if any(
            not isinstance(header, str)
            or not header.startswith("#")
            or "\n" in header
            or "\r" in header
            for header in self.headers
        ):
            raise ConfigurationError("Invalid AGP document header.", code="INVALID_AGP_DOCUMENT")
        versions = tuple(header for header in self.headers if header.startswith("##agp-version"))
        if len(versions) > 1 or versions not in {(), ("##agp-version 2.1",)}:
            raise ConfigurationError(
                "AGP version pragma must be absent or exactly '##agp-version 2.1'.",
                code="INVALID_AGP_DOCUMENT",
            )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def _positive_integer(value: str, *, field: str, line_number: int) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise InputFormatError(
            "AGP coordinate fields must be integers.",
            code="INVALID_AGP_COORDINATE",
            context={"line_number": line_number, "field": field},
        ) from exc
    if result < 1:
        raise InputFormatError(
            "AGP coordinate fields must be positive.",
            code="INVALID_AGP_COORDINATE",
            context={"line_number": line_number, "field": field},
        )
    return result


def read_agp(
    source: AnnotationSource,
    *,
    max_records: int = 1_000_000,
    max_line_length: int = 1_000_000,
    max_header_lines: int = 10_000,
) -> AGPDocument:
    """Read AGP 2.1 components and explicit gaps with continuity validation."""

    _validate_limit(max_records, "max_records")
    _validate_limit(max_line_length, "max_line_length")
    _validate_limit(max_header_lines, "max_header_lines")
    entries: list[object] = []
    headers: list[str] = []
    state: dict[str, tuple[int, int]] = {}
    active_object: str | None = None
    completed_objects: set[str] = set()
    seen_body = False
    with open_text_source(source) as handle:
        for line_number, line in checked_lines(handle, max_line_length=max_line_length):
            if not line:
                continue
            if line.startswith("#"):
                if seen_body:
                    raise InputFormatError(
                        "AGP comments are allowed only before the body.",
                        code="INVALID_AGP_HEADER_POSITION",
                        context={"line_number": line_number},
                    )
                if line.startswith("##agp-version") and line != "##agp-version 2.1":
                    raise InputFormatError(
                        "Only AGP version 2.1 is supported.",
                        code="INVALID_AGP_VERSION",
                        context={"line_number": line_number},
                    )
                _bounded_header(headers, line, max_records=max_header_lines)
                continue
            seen_body = True
            columns = line.split("\t")
            if len(columns) != 9:
                raise InputFormatError(
                    "AGP rows must contain exactly nine tab-separated columns.",
                    code="INVALID_AGP_COLUMNS",
                    context={"line_number": line_number, "column_count": len(columns)},
                )
            object_id, begin_raw, end_raw, part_raw, component_type = columns[:5]
            if not object_id:
                raise InputFormatError(
                    "AGP object ID must be non-empty.",
                    code="INVALID_AGP_FIELD",
                    context={"line_number": line_number},
                )
            if active_object != object_id:
                if object_id in completed_objects:
                    raise InputFormatError(
                        "AGP object rows must form one contiguous block.",
                        code="INVALID_AGP_CONTINUITY",
                        context={"line_number": line_number, "object_id": object_id},
                    )
                if active_object is not None:
                    completed_objects.add(active_object)
                active_object = object_id
            begin = _positive_integer(begin_raw, field="object_beg", line_number=line_number)
            end = _positive_integer(end_raw, field="object_end", line_number=line_number)
            part_number = _positive_integer(part_raw, field="part_number", line_number=line_number)
            if end < begin:
                raise InputFormatError(
                    "AGP object interval end precedes begin.",
                    code="INVALID_AGP_COORDINATE",
                    context={"line_number": line_number},
                )
            expected = state.get(object_id, (1, 1))
            if (begin, part_number) != expected:
                raise InputFormatError(
                    "AGP rows for each object must be contiguous with sequential part numbers.",
                    code="INVALID_AGP_CONTINUITY",
                    context={
                        "line_number": line_number,
                        "expected_begin": expected[0],
                        "expected_part_number": expected[1],
                    },
                )
            state[object_id] = (end + 1, part_number + 1)
            object_interval = Interval(begin - 1, end)
            if component_type in _COMPONENT_TYPES:
                component_id, component_begin_raw, component_end_raw, orientation = columns[5:]
                if not component_id:
                    raise InputFormatError(
                        "AGP component ID must be non-empty.",
                        code="INVALID_AGP_FIELD",
                        context={"line_number": line_number},
                    )
                component_begin = _positive_integer(
                    component_begin_raw, field="component_beg", line_number=line_number
                )
                component_end = _positive_integer(
                    component_end_raw, field="component_end", line_number=line_number
                )
                if component_end < component_begin or component_end - component_begin + 1 != len(
                    object_interval
                ):
                    raise InputFormatError(
                        "AGP component span must equal its object span.",
                        code="INVALID_AGP_COMPONENT_SPAN",
                        context={"line_number": line_number},
                    )
                if orientation not in {"+", "-", "?", "0", "na"}:
                    raise InputFormatError(
                        "Invalid AGP component orientation.",
                        code="INVALID_AGP_ORIENTATION",
                        context={"line_number": line_number},
                    )
                entry: AGPEntry = AGPComponent(
                    object_id,
                    object_interval,
                    part_number,
                    component_type,
                    component_id,
                    Interval(component_begin - 1, component_end),
                    orientation,
                )
            elif component_type in _GAP_TYPES:
                gap_length_raw, gap_type, linkage_raw, evidence_raw = columns[5:]
                gap_length = _positive_integer(
                    gap_length_raw, field="gap_length", line_number=line_number
                )
                if gap_length != len(object_interval):
                    raise InputFormatError(
                        "AGP gap length must equal its object span.",
                        code="INVALID_AGP_GAP_SPAN",
                        context={"line_number": line_number},
                    )
                if component_type == "U" and gap_length != 100:
                    raise InputFormatError(
                        "AGP U gaps must use length 100.",
                        code="INVALID_AGP_GAP_SPAN",
                        context={"line_number": line_number},
                    )
                if linkage_raw not in {"yes", "no"}:
                    raise InputFormatError(
                        "AGP linkage must be yes or no.",
                        code="INVALID_AGP_LINKAGE",
                        context={"line_number": line_number},
                    )
                evidence = () if evidence_raw == "na" else tuple(evidence_raw.split(";"))
                if any(not item for item in evidence):
                    raise InputFormatError(
                        "AGP linkage evidence is malformed.",
                        code="INVALID_AGP_EVIDENCE",
                        context={"line_number": line_number},
                    )
                if gap_type not in _AGP_GAP_TYPES:
                    raise InputFormatError(
                        "Invalid AGP gap type.",
                        code="INVALID_AGP_GAP_TYPE",
                        context={"line_number": line_number, "gap_type": gap_type},
                    )
                if any(item not in _AGP_LINKAGE_EVIDENCE for item in evidence):
                    raise InputFormatError(
                        "Invalid AGP linkage evidence.",
                        code="INVALID_AGP_EVIDENCE",
                        context={"line_number": line_number},
                    )
                if linkage_raw == "no" and evidence:
                    raise InputFormatError(
                        "Unlinked AGP gaps must use 'na' linkage evidence.",
                        code="INVALID_AGP_EVIDENCE",
                        context={"line_number": line_number},
                    )
                if linkage_raw == "yes" and not evidence:
                    raise InputFormatError(
                        "Linked AGP gaps require linkage evidence.",
                        code="INVALID_AGP_EVIDENCE",
                        context={"line_number": line_number},
                    )
                valid_linkage = (
                    (gap_type in {"scaffold", "contamination"} and linkage_raw == "yes")
                    or (
                        gap_type
                        in {"contig", "centromere", "short_arm", "heterochromatin", "telomere"}
                        and linkage_raw == "no"
                    )
                    or gap_type == "repeat"
                )
                if not valid_linkage or ("unspecified" in evidence and gap_type != "contamination"):
                    raise InputFormatError(
                        "AGP gap type, linkage, and evidence combination is invalid.",
                        code="INVALID_AGP_LINKAGE_COMBINATION",
                        context={"line_number": line_number, "gap_type": gap_type},
                    )
                gap = Gap(
                    None if component_type == "U" else gap_length,
                    kind=_GAP_KIND_BY_AGP.get(gap_type, GapKind.UNKNOWN),
                    crossable=linkage_raw == "yes",
                    evidence=evidence,
                    metadata={"agp_gap_type": gap_type},
                )
                entry = AGPGap(
                    object_id,
                    object_interval,
                    part_number,
                    component_type,
                    gap,
                    gap_type,
                    linkage_raw == "yes",
                    evidence,
                )
            else:
                raise InputFormatError(
                    "Unknown AGP component type.",
                    code="INVALID_AGP_COMPONENT_TYPE",
                    context={"line_number": line_number, "component_type": component_type},
                )
            _bounded_append(entries, entry, max_records=max_records)
    return AGPDocument(
        cast(tuple[AGPEntry, ...], tuple(entries)), tuple(headers) or ("##agp-version 2.1",)
    )


def write_agp(
    document: AGPDocument | Iterable[AGPEntry],
    target: AnnotationTarget,
    *,
    overwrite: bool = False,
    create_parents: bool = False,
) -> int:
    """Write deterministic AGP rows after validating object continuity."""

    entries: Iterable[AGPEntry]
    headers: tuple[str, ...]
    if isinstance(document, AGPDocument):
        entries, headers = document.entries, document.headers
    else:
        entries, headers = document, ("##agp-version 2.1",)

    def writer(handle: TextIO) -> int:
        for header in headers:
            if not header.startswith("#"):
                raise InputFormatError(
                    "AGP headers must begin with '#'.", code="INVALID_AGP_HEADER"
                )
            handle.write(header + "\n")
        state: dict[str, tuple[int, int]] = {}
        count = 0
        for entry in entries:
            if not isinstance(entry, (AGPComponent, AGPGap)):
                raise TypeError("AGP entries must be AGPComponent or AGPGap objects.")
            begin, end = entry.object_interval.start + 1, entry.object_interval.end
            expected = state.get(entry.object_id, (1, 1))
            if (begin, entry.part_number) != expected:
                raise InputFormatError(
                    "AGP entries must be contiguous and sequential.", code="INVALID_AGP_CONTINUITY"
                )
            state[entry.object_id] = (end + 1, entry.part_number + 1)
            prefix = [
                entry.object_id,
                str(begin),
                str(end),
                str(entry.part_number),
                entry.component_type,
            ]
            if isinstance(entry, AGPComponent):
                if entry.component_type not in _COMPONENT_TYPES or len(
                    entry.component_interval
                ) != len(entry.object_interval):
                    raise InputFormatError(
                        "Invalid AGP component entry.", code="INVALID_AGP_COMPONENT_SPAN"
                    )
                suffix = [
                    entry.component_id,
                    str(entry.component_interval.start + 1),
                    str(entry.component_interval.end),
                    entry.orientation,
                ]
            else:
                if entry.component_type not in _GAP_TYPES:
                    raise InputFormatError(
                        "Invalid AGP gap component type.", code="INVALID_AGP_COMPONENT_TYPE"
                    )
                gap_length = len(entry.object_interval)
                if entry.gap.length is not None and entry.gap.length != gap_length:
                    raise InputFormatError(
                        "AGP gap length does not equal object span.", code="INVALID_AGP_GAP_SPAN"
                    )
                evidence = ";".join(entry.linkage_evidence) if entry.linkage_evidence else "na"
                suffix = [
                    str(gap_length),
                    entry.gap_type,
                    "yes" if entry.linkage else "no",
                    evidence,
                ]
            fields = prefix + suffix
            if any(
                not value or "\t" in value or "\n" in value or "\r" in value for value in fields
            ):
                raise InputFormatError(
                    "AGP fields must be non-empty and single-line.", code="INVALID_AGP_FIELD"
                )
            handle.write("\t".join(fields) + "\n")
            count += 1
        return count

    return write_text_path(target, writer, overwrite=overwrite, create_parents=create_parents)


__all__ = [
    "AGPComponent",
    "AGPDocument",
    "AGPEntry",
    "AGPGap",
    "AnnotationDocument",
    "AnnotationEntry",
    "AnnotationSource",
    "AnnotationTarget",
    "read_agp",
    "read_bed",
    "read_gff3",
    "write_agp",
    "write_bed",
    "write_gff3",
]
