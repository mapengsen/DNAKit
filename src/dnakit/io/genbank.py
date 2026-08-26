"""A deliberately bounded, documented subset of the GenBank flat-file format.

The codec supports common LOCUS, DEFINITION, ACCESSION, VERSION, FEATURES and
ORIGIN records. It does not claim complete INSDC compatibility and rejects
fuzzy, remote and other unsupported feature-location syntax explicitly.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Mapping
from typing import TextIO, cast

from dnakit.core._json import JSONValue
from dnakit.core.coordinates import CompoundLocation, Interval, Location, UnresolvedLocation
from dnakit.core.enums import DNAAlphabet, Strand, Strandedness, Topology
from dnakit.core.feature import DNAFeature
from dnakit.core.record import DNARecord
from dnakit.core.sequence import STRICT_SYMBOLS, DNASequence
from dnakit.exceptions import InputFormatError

from .config import ReadConfig, WriteConfig

_DATE_RE = re.compile(r"^\d{2}-[A-Z]{3}-\d{4}$")
_INTERVAL_RE = re.compile(r"^(\d+)(?:\.\.(\d+))?$")
_ORIGIN_TOKEN_RE = re.compile(r"[A-Za-z]+")
_RESERVED = {
    "dnakit_id",
    "dnakit_label",
    "dnakit_score",
    "dnakit_phase",
    "dnakit_source",
}


def _error(message: str, code: str, line_number: int | None = None) -> InputFormatError:
    context: dict[str, object] = {}
    if line_number is not None:
        context["line_number"] = line_number
    return InputFormatError(message, code=code, context=context)


def _parse_location(text: str, *, line_number: int) -> tuple[Location, Strand]:
    value = text.strip()
    strand = Strand.FORWARD
    if value.startswith("complement(") and value.endswith(")"):
        value = value[11:-1].strip()
        strand = Strand.REVERSE
    if value.startswith("join(") and value.endswith(")"):
        tokens = tuple(item.strip() for item in value[5:-1].split(","))
    else:
        tokens = (value,)
    parts: list[Interval] = []
    for token in tokens:
        match = _INTERVAL_RE.fullmatch(token)
        if match is None:
            raise _error(
                "The GenBank subset does not support this feature location.",
                "UNSUPPORTED_GENBANK_LOCATION",
                line_number,
            )
        start = int(match.group(1))
        end = int(match.group(2) or match.group(1))
        if start < 1 or end < start:
            raise _error(
                "GenBank feature coordinates are invalid.", "INVALID_GENBANK_LOCATION", line_number
            )
        parts.append(Interval(start - 1, end))
    if not parts:
        raise _error("GenBank join() must not be empty.", "INVALID_GENBANK_LOCATION", line_number)
    return (parts[0] if len(parts) == 1 else CompoundLocation(parts)), strand


def _qualifier_scalar(value: str) -> JSONValue:
    if value == "true":
        return True
    if value == "false":
        return False
    if value == "null":
        return None
    if re.fullmatch(r"-?(?:0|[1-9]\d*)", value):
        return int(value)
    if re.fullmatch(r"-?(?:(?:0|[1-9]\d*)\.\d+|(?:0|[1-9]\d*)(?:\.\d+)?[eE][+-]?\d+)", value):
        return float(value)
    return value


def _parse_qualifier(text: str, *, line_number: int) -> tuple[str, JSONValue]:
    body = text.strip()
    if not body.startswith("/"):
        raise _error(
            "Malformed GenBank feature qualifier.", "INVALID_GENBANK_QUALIFIER", line_number
        )
    body = body[1:]
    value: JSONValue
    if "=" not in body:
        key, value = body, True
    else:
        key, raw_value = body.split("=", 1)
        if raw_value.startswith('"'):
            if not raw_value.endswith('"') or len(raw_value) < 2:
                raise _error(
                    "Multiline GenBank qualifiers are outside the supported subset.",
                    "UNSUPPORTED_GENBANK_QUALIFIER",
                    line_number,
                )
            value = raw_value[1:-1].replace('""', '"')
        else:
            value = _qualifier_scalar(raw_value)
    if not key or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", key):
        raise _error("Invalid GenBank qualifier name.", "INVALID_GENBANK_QUALIFIER", line_number)
    return key, value


def _append_qualifier(values: dict[str, object], key: str, value: JSONValue) -> None:
    previous = values.get(key)
    if previous is None:
        values[key] = value
    elif isinstance(previous, list):
        previous.append(value)
    else:
        values[key] = [previous, value]


def _feature_from_parts(
    key: str,
    location_text: str,
    qualifier_items: list[tuple[str, JSONValue]],
    *,
    line_number: int,
) -> DNAFeature:
    location, strand = _parse_location(location_text, line_number=line_number)
    qualifiers: dict[str, object] = {}
    reserved: dict[str, JSONValue] = {}
    for name, value in qualifier_items:
        if name in _RESERVED:
            if name in reserved:
                raise _error(
                    "A DNAKit GenBank qualifier is duplicated.",
                    "INVALID_GENBANK_QUALIFIER",
                    line_number,
                )
            reserved[name] = value
        else:
            _append_qualifier(qualifiers, name, value)
    try:
        identifier = reserved.get("dnakit_id")
        label = reserved.get("dnakit_label")
        score = reserved.get("dnakit_score")
        phase = reserved.get("dnakit_phase")
        source = reserved.get("dnakit_source")
        for name, value in (
            ("dnakit_id", identifier),
            ("dnakit_label", label),
            ("dnakit_source", source),
        ):
            if value is not None and not isinstance(value, str):
                raise TypeError(f"{name} must be text")
        if score is not None and (isinstance(score, bool) or not isinstance(score, (int, float))):
            raise TypeError("dnakit_score must be numeric")
        if phase is not None and (isinstance(phase, bool) or not isinstance(phase, int)):
            raise TypeError("dnakit_phase must be an integer")
        return DNAFeature(
            key,
            location,
            id=identifier if isinstance(identifier, str) else None,
            strand=strand,
            label=label if isinstance(label, str) else None,
            score=float(score)
            if isinstance(score, (int, float)) and not isinstance(score, bool)
            else None,
            phase=phase if isinstance(phase, int) and not isinstance(phase, bool) else None,
            qualifiers=qualifiers,
            source=source if isinstance(source, str) else None,
        )
    except (TypeError, ValueError) as exc:
        raise _error(
            "Invalid typed DNAKit GenBank qualifier.", "INVALID_GENBANK_QUALIFIER", line_number
        ) from exc


def iter_genbank(handle: TextIO, config: ReadConfig) -> Iterator[DNARecord]:
    """Parse the auditable GenBank subset lazily."""

    lines: list[tuple[int, str]] = []
    for line_number, raw in enumerate(handle, start=1):
        line = raw.rstrip("\r\n")
        if len(line) > config.max_field_size:
            raise _error(
                "A GenBank line exceeds max_field_size.", "GENBANK_LINE_TOO_LONG", line_number
            )
        if line == "//":
            yield _parse_record(lines, config)
            lines = []
        else:
            if len(lines) >= config.max_record_lines:
                raise _error(
                    "A GenBank record exceeds max_record_lines.",
                    "GENBANK_RECORD_LINE_LIMIT_EXCEEDED",
                    line_number,
                )
            lines.append((line_number, line))
    if lines:
        raise _error(
            "GenBank record is missing the // terminator.", "TRUNCATED_GENBANK_RECORD", lines[-1][0]
        )


def _parse_record(lines: list[tuple[int, str]], config: ReadConfig) -> DNARecord:
    if not lines or not lines[0][1].startswith("LOCUS"):
        raise _error(
            "GenBank record must begin with LOCUS.",
            "GENBANK_MISSING_LOCUS",
            lines[0][0] if lines else None,
        )
    locus_fields = lines[0][1].split()
    if len(locus_fields) < 4 or locus_fields[0] != "LOCUS" or locus_fields[2].isdigit() is False:
        raise _error("Malformed GenBank LOCUS line.", "INVALID_GENBANK_LOCUS", lines[0][0])
    locus = locus_fields[1]
    expected_length = int(locus_fields[2])
    if expected_length < 0 or locus_fields[3].lower() != "bp":
        raise _error("Malformed GenBank LOCUS length.", "INVALID_GENBANK_LOCUS", lines[0][0])
    if expected_length > config.max_sequence_symbols:
        raise _error(
            "GenBank record exceeds max_sequence_symbols.",
            "SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
            lines[0][0],
        )
    lower_fields = [item.lower() for item in locus_fields]
    topology = Topology.CIRCULAR if "circular" in lower_fields else Topology.LINEAR
    strandedness = (
        Strandedness.DOUBLE
        if any("ds-dna" in item for item in lower_fields)
        else config.strandedness
    )
    division = (
        locus_fields[-2]
        if len(locus_fields) >= 2 and _DATE_RE.fullmatch(locus_fields[-1])
        else None
    )
    date = locus_fields[-1] if _DATE_RE.fullmatch(locus_fields[-1]) else None

    definition_parts: list[str] = []
    accession: str | None = None
    version: str | None = None
    feature_rows: list[tuple[str, str, int, list[tuple[str, JSONValue]]]] = []
    sequence_chunks: list[str] = []
    section = "header"
    seen_origin = False
    sequence_symbol_count = 0
    active_header_tag: str | None = None
    for line_number, line in lines[1:]:
        if line.startswith("FEATURES"):
            section = "features"
            continue
        if line.startswith("ORIGIN"):
            section = "origin"
            seen_origin = True
            continue
        if section == "header":
            tag = line[:12].strip()
            value = line[12:].strip()
            if tag:
                active_header_tag = tag
            if tag == "DEFINITION" or (not tag and active_header_tag == "DEFINITION"):
                definition_parts.append(value)
            elif tag == "ACCESSION":
                fields = value.split()
                accession = fields[0] if fields else None
            elif tag == "VERSION":
                fields = value.split()
                version = fields[0] if fields else None
        elif section == "features":
            key = line[5:21].strip() if len(line) >= 21 else ""
            value = line[21:].strip() if len(line) >= 21 else ""
            if key:
                feature_rows.append((key, value, line_number, []))
            elif value.startswith("/"):
                if not feature_rows:
                    raise _error(
                        "GenBank qualifier appears before a feature.",
                        "INVALID_GENBANK_QUALIFIER",
                        line_number,
                    )
                name, qualifier = _parse_qualifier(value, line_number=line_number)
                feature_rows[-1][3].append((name, qualifier))
            elif value:
                raise _error(
                    "Wrapped GenBank feature locations are outside the supported subset.",
                    "UNSUPPORTED_GENBANK_LOCATION",
                    line_number,
                )
            elif line.strip():
                raise _error(
                    "Malformed GenBank feature line.",
                    "INVALID_GENBANK_FEATURE",
                    line_number,
                )
        else:
            stripped = line.strip()
            if stripped and re.fullmatch(r"\d+(?:\s+[A-Za-z]+)*", stripped) is None:
                raise _error(
                    "Malformed GenBank ORIGIN sequence line.",
                    "INVALID_GENBANK_ORIGIN",
                    line_number,
                )
            tokens = _ORIGIN_TOKEN_RE.findall(line)
            sequence_symbol_count += sum(len(token) for token in tokens)
            if sequence_symbol_count > config.max_sequence_symbols:
                raise _error(
                    "GenBank record exceeds max_sequence_symbols.",
                    "SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
                    line_number,
                )
            sequence_chunks.extend(tokens)

    sequence_text = "".join(sequence_chunks).upper()
    if not seen_origin:
        raise _error("GenBank record is missing ORIGIN.", "GENBANK_MISSING_ORIGIN", lines[0][0])
    if len(sequence_text) != expected_length:
        raise _error(
            "GenBank LOCUS length does not equal the ORIGIN sequence length.",
            "GENBANK_LENGTH_MISMATCH",
            lines[0][0],
        )
    alphabet = DNAAlphabet.STRICT if set(sequence_text) <= STRICT_SYMBOLS else config.alphabet
    try:
        sequence = DNASequence(
            sequence_text,
            alphabet=alphabet,
            topology=topology,
            strandedness=strandedness,
        )
        features = tuple(
            _feature_from_parts(key, location, qualifiers, line_number=line_number)
            for key, location, line_number, qualifiers in feature_rows
        )
        metadata = {
            "genbank": {
                "locus": locus,
                "accession": accession,
                "version": version,
                "division": division,
                "date": date,
                "codec": "dnakit.genbank.subset.v1",
            }
        }
        description = " ".join(definition_parts).strip()
        return DNARecord(
            sequence,
            version or accession or locus,
            description="" if description == "." else description,
            features=features,
            metadata=metadata,
        )
    except Exception as exc:
        if isinstance(exc, InputFormatError):
            raise
        raise _error(
            "GenBank record violates DNAKit object constraints.",
            "INVALID_GENBANK_RECORD",
            lines[0][0],
        ) from exc


def _format_location(location: Location, strand: Strand) -> str:
    if isinstance(location, UnresolvedLocation):
        raise InputFormatError(
            "Unresolved locations cannot be represented by the GenBank subset.",
            code="UNSUPPORTED_GENBANK_LOCATION",
        )
    parts = (location,) if isinstance(location, Interval) else location.parts
    tokens = [
        f"{part.start + 1}..{part.end}" if len(part) != 1 else str(part.start + 1) for part in parts
    ]
    value = tokens[0] if len(tokens) == 1 else f"join({','.join(tokens)})"
    if strand is Strand.REVERSE:
        return f"complement({value})"
    if strand not in {Strand.FORWARD, Strand.UNKNOWN}:
        raise InputFormatError(
            "The GenBank subset represents only forward, reverse, or unknown feature strands.",
            code="UNSUPPORTED_GENBANK_STRAND",
        )
    return value


def _serialize_qualifier(name: str, value: object) -> str:
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
        raise InputFormatError("Invalid GenBank qualifier name.", code="INVALID_GENBANK_QUALIFIER")
    if value is True:
        return f"/{name}"
    if value is False or value is None or isinstance(value, (int, float)):
        return (
            f"/{name}={str(value).lower() if isinstance(value, bool) or value is None else value}"
        )
    if not isinstance(value, str) or "\n" in value or "\r" in value:
        raise InputFormatError(
            "GenBank qualifier values must be JSON scalars or arrays of scalars.",
            code="UNSUPPORTED_GENBANK_QUALIFIER",
        )
    return f'/{name}="{value.replace(chr(34), chr(34) * 2)}"'


def _feature_qualifiers(feature: DNAFeature) -> Iterator[tuple[str, object]]:
    for key in feature.qualifiers:
        if key in _RESERVED:
            raise InputFormatError(
                "Feature qualifiers conflict with reserved DNAKit GenBank names.",
                code="RESERVED_GENBANK_QUALIFIER",
                context={"qualifier": key},
            )
    if feature.id is not None:
        yield "dnakit_id", feature.id
    if feature.label is not None:
        yield "dnakit_label", feature.label
    if feature.score is not None:
        yield "dnakit_score", feature.score
    if feature.phase is not None:
        yield "dnakit_phase", feature.phase
    if feature.source is not None:
        yield "dnakit_source", feature.source
    for key, value in feature.qualifiers.items():
        if isinstance(value, tuple):
            for item in value:
                yield key, item
        else:
            yield key, value


def write_genbank(handle: TextIO, records: Iterable[DNARecord], config: WriteConfig) -> int:
    """Serialize records in the deterministic DNAKit GenBank subset."""

    count = 0
    for record in records:
        if record.sequence.is_gapped:
            raise InputFormatError(
                "The GenBank subset cannot serialize explicit Gap objects.",
                code="GAP_LOSS_NOT_ALLOWED",
                context={"record_id": record.id},
            )
        if any(symbol.isspace() for symbol in record.id) or len(record.id) > 32:
            raise InputFormatError(
                "GenBank record IDs must be 1-32 characters without whitespace.",
                code="INVALID_GENBANK_ID",
                context={"record_id": record.id},
            )
        metadata_raw: object = record.metadata.get("genbank", {})
        metadata = metadata_raw if isinstance(metadata_raw, Mapping) else {}
        locus = metadata.get("locus", record.id)
        accession = metadata.get("accession", record.id)
        version = metadata.get("version", record.id)
        division = metadata.get("division", "UNK") or "UNK"
        date = metadata.get("date", "01-JAN-1980") or "01-JAN-1980"
        for label, value in (
            ("locus", locus),
            ("accession", accession),
            ("version", version),
            ("division", division),
            ("date", date),
        ):
            if not isinstance(value, str) or not value or any(ch in value for ch in "\r\n"):
                raise InputFormatError(
                    "GenBank metadata fields must be non-empty single-line strings.",
                    code="INVALID_GENBANK_METADATA",
                    context={"field": label},
                )
        if (
            len(cast(str, locus)) > 32
            or any(symbol.isspace() for symbol in cast(str, locus))
            or any(symbol.isspace() for symbol in cast(str, accession))
            or any(symbol.isspace() for symbol in cast(str, version))
            or re.fullmatch(r"[A-Z]{3}", cast(str, division)) is None
        ):
            raise InputFormatError(
                "GenBank locus/accession/version/division metadata is not representable.",
                code="INVALID_GENBANK_METADATA",
            )
        if not _DATE_RE.fullmatch(cast(str, date)):
            raise InputFormatError(
                "GenBank date must use DD-MMM-YYYY.", code="INVALID_GENBANK_METADATA"
            )
        sequence = record.sequence.to_string()
        topology = record.sequence.topology.value
        stranded = "ds-DNA" if record.sequence.strandedness is Strandedness.DOUBLE else "DNA"
        locus_prefix = f"LOCUS       {locus:<32} {len(sequence):>11} bp    "
        handle.write(f"{locus_prefix}{stranded:<6} {topology:<8} {division:>3} {date}\n")
        definition = record.description or "."
        if "\n" in definition or "\r" in definition:
            raise InputFormatError(
                "GenBank descriptions must be single-line text.", code="INVALID_GENBANK_METADATA"
            )
        handle.write(f"DEFINITION  {definition}\n")
        handle.write(f"ACCESSION   {accession}\n")
        handle.write(f"VERSION     {version}\n")
        handle.write("FEATURES             Location/Qualifiers\n")
        for feature in record.features:
            location = _format_location(feature.location, feature.strand)
            if len(feature.type) > 15 or any(symbol.isspace() for symbol in feature.type):
                raise InputFormatError(
                    "GenBank feature keys must fit the 15-character key field.",
                    code="INVALID_GENBANK_FEATURE",
                )
            handle.write(f"     {feature.type:<15} {location}\n")
            for name, value in _feature_qualifiers(feature):
                rendered = _serialize_qualifier(name, value)
                if len(rendered) > 58:
                    raise InputFormatError(
                        "Multiline GenBank qualifiers are outside the supported subset.",
                        code="UNSUPPORTED_GENBANK_QUALIFIER",
                    )
                handle.write(f"                     {rendered}\n")
        handle.write("ORIGIN\n")
        lower = sequence.lower()
        for offset in range(0, len(lower), 60):
            block = lower[offset : offset + 60]
            grouped = " ".join(block[index : index + 10] for index in range(0, len(block), 10))
            handle.write(f"{offset + 1:>9} {grouped}\n")
        handle.write("//\n")
        count += 1
    return count


__all__ = ["iter_genbank", "write_genbank"]
