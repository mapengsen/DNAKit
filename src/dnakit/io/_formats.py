"""Streaming FASTA, FASTQ, and structured table codecs."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Iterator, Mapping
from typing import Any, TextIO, cast

from dnakit.core._json import to_json_compatible
from dnakit.core.coordinates import CompoundLocation, Interval, Location, UnresolvedLocation
from dnakit.core.enums import DNAAlphabet, Strandedness, Topology
from dnakit.core.feature import DNAFeature
from dnakit.core.gap import Gap
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError

from .config import ReadConfig, WriteConfig

_PHRED_KEY = "phred_quality"
_RECORD_SCHEMA_VERSION = "dnakit.record.v1"


def _json_object_without_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> object:
    raise ValueError(f"non-standard JSON constant: {value}")


def _validate_json_structure(
    value: object,
    config: ReadConfig,
    *,
    context: Mapping[str, object] | None = None,
) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(value, 1)]
    while stack:
        current, depth = stack.pop()
        nodes += 1
        if depth > config.max_json_depth or nodes > config.max_json_nodes:
            raise InputFormatError(
                "JSON structure exceeds its configured depth or node limit.",
                code="JSON_STRUCTURE_LIMIT_EXCEEDED",
                context={
                    **dict(context or {}),
                    "max_json_depth": config.max_json_depth,
                    "max_json_nodes": config.max_json_nodes,
                },
            )
        if isinstance(current, Mapping):
            nodes += len(current)
            if nodes > config.max_json_nodes:
                raise InputFormatError(
                    "JSON structure exceeds its configured depth or node limit.",
                    code="JSON_STRUCTURE_LIMIT_EXCEEDED",
                    context={
                        **dict(context or {}),
                        "max_json_depth": config.max_json_depth,
                        "max_json_nodes": config.max_json_nodes,
                    },
                )
            stack.extend((item, depth + 1) for item in current.values())
        elif isinstance(current, list):
            stack.extend((item, depth + 1) for item in current)


def _strict_json_loads(
    text: str,
    config: ReadConfig,
    *,
    message: str,
    code: str,
    context: Mapping[str, object] | None = None,
) -> object:
    try:
        value = json.loads(
            text,
            object_pairs_hook=_json_object_without_duplicates,
            parse_constant=_reject_json_constant,
        )
    except RecursionError as exc:
        raise InputFormatError(
            "JSON structure exceeds its configured depth or node limit.",
            code="JSON_STRUCTURE_LIMIT_EXCEEDED",
            context={
                **dict(context or {}),
                "max_json_depth": config.max_json_depth,
                "max_json_nodes": config.max_json_nodes,
            },
        ) from exc
    except (json.JSONDecodeError, ValueError) as exc:
        details = dict(context or {})
        if isinstance(exc, json.JSONDecodeError):
            details.setdefault("line_number", exc.lineno)
            details["column_number"] = exc.colno
        raise InputFormatError(message, code=code, context=details) from exc
    _validate_json_structure(value, config, context=context)
    return value


def _input_error(
    message: str,
    *,
    code: str,
    line_number: int | None = None,
    record_id: str | None = None,
) -> InputFormatError:
    context: dict[str, object] = {}
    if line_number is not None:
        context["line_number"] = line_number
    if record_id is not None:
        context["record_id"] = record_id
    return InputFormatError(message, code=code, context=context)


def _sequence(
    text: str,
    config: ReadConfig,
    *,
    record_id: str,
    line_number: int,
    alphabet: DNAAlphabet | str | None = None,
    topology: Topology | str | None = None,
    strandedness: Strandedness | str | None = None,
) -> DNASequence:
    if len(text) > config.max_sequence_symbols:
        raise _input_error(
            "Record sequence exceeds max_sequence_symbols.",
            code="SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
            line_number=line_number,
            record_id=record_id,
        )
    value = text.upper() if config.uppercase else text
    try:
        return DNASequence(
            value,
            alphabet=config.alphabet if alphabet is None else alphabet,
            topology=config.topology if topology is None else topology,
            strandedness=config.strandedness if strandedness is None else strandedness,
        )
    except DNAKitError as exc:
        raise _input_error(
            "Record contains symbols incompatible with the configured DNA alphabet.",
            code="INVALID_SEQUENCE_CONTENT",
            line_number=line_number,
            record_id=record_id,
        ) from exc


def _header(value: str, *, line_number: int, marker: str) -> tuple[str, str]:
    body = value[1:].strip()
    if not body:
        raise _input_error(
            f"{marker} header must contain a record ID.",
            code=f"EMPTY_{marker}_HEADER",
            line_number=line_number,
        )
    fields = body.split(maxsplit=1)
    return fields[0], fields[1] if len(fields) == 2 else ""


def _serialized_header(record: DNARecord, *, format: str) -> str:
    if any(symbol.isspace() for symbol in record.id):
        raise InputFormatError(
            f"{format} record IDs must not contain whitespace.",
            code="INVALID_SEQUENCE_HEADER",
            context={"record_id": record.id, "format": format},
        )
    if "\n" in record.description or "\r" in record.description:
        raise InputFormatError(
            f"{format} descriptions must not contain line breaks.",
            code="INVALID_SEQUENCE_HEADER",
            context={"record_id": record.id, "format": format},
        )
    return record.id if not record.description else f"{record.id} {record.description}"


def iter_fasta(handle: TextIO, config: ReadConfig) -> Iterator[DNARecord]:
    """Parse multiline FASTA records without materializing the file."""

    record_id: str | None = None
    description = ""
    chunks: list[str] = []
    header_line = 0
    record_line_count = 0
    symbol_count = 0
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.rstrip("\r\n")
        if len(line) > config.max_field_size:
            raise _input_error(
                "FASTA line exceeds max_field_size.",
                code="SEQUENCE_LINE_TOO_LONG",
                line_number=line_number,
            )
        if line.startswith(">"):
            if record_id is not None:
                yield DNARecord(
                    _sequence(
                        "".join(chunks),
                        config,
                        record_id=record_id,
                        line_number=header_line,
                    ),
                    record_id,
                    description=description,
                )
            record_id, description = _header(line, line_number=line_number, marker="FASTA")
            chunks = []
            symbol_count = 0
            header_line = line_number
            record_line_count = 1
        else:
            record_line_count += 1
            if record_line_count > config.max_record_lines:
                raise _input_error(
                    "FASTA record exceeds max_record_lines.",
                    code="SEQUENCE_RECORD_LINE_LIMIT_EXCEEDED",
                    line_number=line_number,
                    record_id=record_id,
                )
            if not line:
                continue
            if record_id is None:
                raise _input_error(
                    "FASTA sequence data appeared before the first header.",
                    code="FASTA_MISSING_HEADER",
                    line_number=line_number,
                )
            symbol_count += len(line)
            if symbol_count > config.max_sequence_symbols:
                raise _input_error(
                    "FASTA record exceeds max_sequence_symbols.",
                    code="SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
                    line_number=line_number,
                    record_id=record_id,
                )
            chunks.append(line)
    if record_id is not None:
        yield DNARecord(
            _sequence("".join(chunks), config, record_id=record_id, line_number=header_line),
            record_id,
            description=description,
        )


def iter_fastq(handle: TextIO, config: ReadConfig) -> Iterator[DNARecord]:
    """Parse strict four-line FASTQ records and validate Phred values."""

    line_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        line_number += 1
        if len(header.rstrip("\r\n")) > config.max_field_size:
            raise _input_error(
                "FASTQ header exceeds max_field_size.",
                code="SEQUENCE_LINE_TOO_LONG",
                line_number=line_number,
            )
        if not header.startswith("@"):
            raise _input_error(
                "FASTQ record must begin with '@'.",
                code="FASTQ_MISSING_HEADER",
                line_number=line_number,
            )
        record_id, description = _header(
            header.rstrip("\r\n"), line_number=line_number, marker="FASTQ"
        )
        sequence_line = handle.readline()
        plus_line = handle.readline()
        quality_line = handle.readline()
        if not sequence_line or not plus_line or not quality_line:
            raise _input_error(
                "FASTQ record ended before all four lines were present.",
                code="TRUNCATED_FASTQ_RECORD",
                line_number=line_number,
                record_id=record_id,
            )
        sequence_line_number = line_number + 1
        line_number += 3
        for offset, value in enumerate((sequence_line, plus_line, quality_line), start=1):
            if len(value.rstrip("\r\n")) > config.max_field_size:
                raise _input_error(
                    "FASTQ line exceeds max_field_size.",
                    code="SEQUENCE_LINE_TOO_LONG",
                    line_number=sequence_line_number + offset - 1,
                    record_id=record_id,
                )
        sequence_text = sequence_line.rstrip("\r\n")
        plus = plus_line.rstrip("\r\n")
        quality_text = quality_line.rstrip("\r\n")
        if len(sequence_text) > config.max_sequence_symbols:
            raise _input_error(
                "FASTQ record exceeds max_sequence_symbols.",
                code="SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
                line_number=sequence_line_number,
                record_id=record_id,
            )
        if not plus.startswith("+"):
            raise _input_error(
                "FASTQ separator line must begin with '+'.",
                code="FASTQ_MISSING_SEPARATOR",
                line_number=line_number - 1,
                record_id=record_id,
            )
        repeated_header = plus[1:].strip()
        if repeated_header and repeated_header.split(maxsplit=1)[0] != record_id:
            raise _input_error(
                "FASTQ '+' header ID does not match the '@' header ID.",
                code="FASTQ_HEADER_MISMATCH",
                line_number=line_number - 1,
                record_id=record_id,
            )
        sequence = _sequence(
            sequence_text,
            config,
            record_id=record_id,
            line_number=sequence_line_number,
        )
        if len(quality_text) != sequence.symbol_length:
            raise _input_error(
                "FASTQ quality length must equal sequence symbol_length.",
                code="FASTQ_QUALITY_LENGTH_MISMATCH",
                line_number=line_number,
                record_id=record_id,
            )
        quality_codepoints = tuple(ord(symbol) for symbol in quality_text)
        quality = tuple(value - config.phred_offset for value in quality_codepoints)
        if any(
            codepoint > 126 or value < 0 or value > 93
            for codepoint, value in zip(quality_codepoints, quality, strict=True)
        ):
            raise _input_error(
                "FASTQ quality characters must be printable ASCII and encode Phred values "
                "from 0 through 93.",
                code="FASTQ_QUALITY_OUT_OF_RANGE",
                line_number=line_number,
                record_id=record_id,
            )
        yield DNARecord(
            sequence,
            record_id,
            description=description,
            letter_annotations={_PHRED_KEY: quality},
        )


def _require_object(value: object, *, record_index: int) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise InputFormatError(
            "Each structured record must be a JSON object.",
            code="INVALID_STRUCTURED_RECORD",
            context={"record_index": record_index},
        )
    return cast(Mapping[str, object], value)


def _record_from_mapping(
    row: Mapping[str, object], config: ReadConfig, *, record_index: int
) -> DNARecord:
    sequence_value = row.get(config.sequence_column)
    parts_value = row.get(config.parts_column)
    if parts_value is None and not isinstance(sequence_value, str):
        raise InputFormatError(
            "Structured input requires a string sequence column or a parts array.",
            code="MISSING_SEQUENCE_COLUMN",
            context={"record_index": record_index, "column": config.sequence_column},
        )
    if sequence_value is not None and not isinstance(sequence_value, str):
        raise InputFormatError(
            "Structured input sequence must be text or null.",
            code="INVALID_SEQUENCE_COLUMN",
            context={"record_index": record_index, "column": config.sequence_column},
        )
    if isinstance(sequence_value, str) and len(sequence_value) > config.max_sequence_symbols:
        raise InputFormatError(
            "Structured record exceeds max_sequence_symbols.",
            code="SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
            context={"record_index": record_index, "column": config.sequence_column},
        )
    schema_version = row.get("schema_version")
    if schema_version not in {None, "", _RECORD_SCHEMA_VERSION}:
        raise InputFormatError(
            "Structured input uses an unsupported DNAKit record schema.",
            code="UNSUPPORTED_RECORD_SCHEMA",
            context={"record_index": record_index, "schema_version": schema_version},
        )
    id_value = row.get(config.id_column)
    if id_value is None or id_value == "":
        record_id = f"{config.anonymous_id_prefix}{record_index + 1}"
    elif isinstance(id_value, str):
        record_id = id_value
    else:
        raise InputFormatError(
            "Structured input record ID must be text.",
            code="INVALID_RECORD_ID_COLUMN",
            context={"record_index": record_index, "column": config.id_column},
        )
    description_value = row.get(config.description_column, "")
    if description_value is None:
        description = ""
    elif isinstance(description_value, str):
        description = description_value
    else:
        raise InputFormatError(
            "Structured input description must be text.",
            code="INVALID_DESCRIPTION_COLUMN",
            context={"record_index": record_index, "column": config.description_column},
        )
    metadata_value = row.get(config.metadata_column, {})
    annotations_value = row.get(config.letter_annotations_column, {})
    features_value = row.get(config.features_column, ())
    alphabet_value = row.get("alphabet")
    topology_value = row.get("topology")
    strandedness_value = row.get("strandedness")
    alphabet_value = None if alphabet_value == "" else alphabet_value
    topology_value = None if topology_value == "" else topology_value
    strandedness_value = None if strandedness_value == "" else strandedness_value
    for name, value in (
        ("alphabet", alphabet_value),
        ("topology", topology_value),
        ("strandedness", strandedness_value),
    ):
        if value is not None and not isinstance(value, str):
            raise InputFormatError(
                f"Structured input {name} must be text or null.",
                code="INVALID_SEQUENCE_TYPE_COLUMN",
                context={"record_index": record_index, "column": name},
            )
    if metadata_value is None:
        metadata_value = {}
    if annotations_value is None:
        annotations_value = {}
    if features_value is None:
        features_value = ()
    if isinstance(parts_value, str):
        parts_value = (
            _strict_json_loads(
                parts_value,
                config,
                message=f"Column {config.parts_column!r} must contain a JSON array.",
                code="INVALID_JSON_COLUMN",
                context={"record_index": record_index, "column": config.parts_column},
            )
            if parts_value
            else None
        )
    for column, value in (
        (config.metadata_column, metadata_value),
        (config.letter_annotations_column, annotations_value),
        (config.features_column, features_value),
    ):
        if isinstance(value, str):
            value = (
                _strict_json_loads(
                    value,
                    config,
                    message=f"Column {column!r} must contain valid JSON.",
                    code="INVALID_JSON_COLUMN",
                    context={"record_index": record_index, "column": column},
                )
                if value
                else ([] if column == config.features_column else {})
            )
            if column == config.metadata_column:
                metadata_value = value
            elif column == config.letter_annotations_column:
                annotations_value = value
            else:
                features_value = value
    if not isinstance(metadata_value, Mapping) or not isinstance(annotations_value, Mapping):
        raise InputFormatError(
            "Metadata and letter_annotations must be objects.",
            code="INVALID_STRUCTURED_METADATA",
            context={"record_index": record_index},
        )
    if parts_value is None:
        assert isinstance(sequence_value, str)
        sequence = _sequence(
            sequence_value,
            config,
            record_id=record_id,
            line_number=record_index + 1,
            alphabet=cast(str | None, alphabet_value),
            topology=cast(str | None, topology_value),
            strandedness=cast(str | None, strandedness_value),
        )
    else:
        sequence = _sequence_from_parts(
            parts_value,
            sequence_text=sequence_value,
            config=config,
            record_id=record_id,
            record_index=record_index,
            alphabet=cast(str | None, alphabet_value),
            topology=cast(str | None, topology_value),
            strandedness=cast(str | None, strandedness_value),
        )
    features = _features_from_value(features_value, record_index=record_index)
    try:
        return DNARecord(
            sequence,
            record_id,
            description=description,
            features=features,
            metadata=metadata_value,
            letter_annotations=cast(Mapping[str, Iterable[int | float]], annotations_value),
        )
    except DNAKitError as exc:
        raise InputFormatError(
            "Structured record violates the DNARecord schema.",
            code="INVALID_STRUCTURED_RECORD",
            context={"record_index": record_index, "record_id": record_id},
        ) from exc


def _sequence_from_parts(
    value: object,
    *,
    sequence_text: str | None,
    config: ReadConfig,
    record_id: str,
    record_index: int,
    alphabet: str | None,
    topology: str | None,
    strandedness: str | None,
) -> DNASequence:
    if not isinstance(value, list):
        raise InputFormatError(
            "Structured input parts must be an array.",
            code="INVALID_SEQUENCE_PARTS_COLUMN",
            context={"record_index": record_index},
        )
    parts: list[str | Gap] = []
    symbol_count = 0
    try:
        for item in value:
            if not isinstance(item, Mapping):
                raise TypeError("each sequence part must be an object")
            kind = item.get("kind")
            if kind == "symbols":
                symbols = item.get("symbols")
                if not isinstance(symbols, str):
                    raise TypeError("a symbol part requires string symbols")
                symbol_count += len(symbols)
                if symbol_count > config.max_sequence_symbols:
                    raise InputFormatError(
                        "Structured record exceeds max_sequence_symbols.",
                        code="SEQUENCE_SYMBOL_LIMIT_EXCEEDED",
                        context={"record_index": record_index, "record_id": record_id},
                    )
                parts.append(symbols.upper() if config.uppercase else symbols)
            elif kind == "gap":
                raw_evidence = item.get("evidence", [])
                raw_metadata = item.get("metadata", {})
                if not isinstance(raw_evidence, list) or not isinstance(raw_metadata, Mapping):
                    raise TypeError("gap evidence and metadata have invalid types")
                parts.append(
                    Gap(
                        cast(int | None, item.get("length")),
                        kind=cast(str, item.get("gap_kind", "unknown")),
                        crossable=cast(bool | None, item.get("crossable")),
                        evidence=cast(list[str], raw_evidence),
                        metadata=cast(Mapping[str, object], raw_metadata),
                    )
                )
            else:
                raise TypeError("sequence part kind must be symbols or gap")
        resolved = DNASequence(
            parts,
            alphabet=config.alphabet if alphabet is None else alphabet,
            topology=config.topology if topology is None else topology,
            strandedness=config.strandedness if strandedness is None else strandedness,
        )
    except InputFormatError:
        raise
    except (DNAKitError, TypeError, ValueError) as exc:
        raise InputFormatError(
            "Structured input contains invalid sequence parts.",
            code="INVALID_SEQUENCE_PARTS_COLUMN",
            context={"record_index": record_index, "record_id": record_id},
        ) from exc
    if sequence_text is not None:
        expected = sequence_text.upper() if config.uppercase else sequence_text
        if expected != resolved.symbols:
            raise InputFormatError(
                "Structured sequence text does not match the symbols in parts.",
                code="SEQUENCE_PARTS_MISMATCH",
                context={"record_index": record_index, "record_id": record_id},
            )
    return resolved


def _interval_from_value(value: object) -> Interval:
    if not isinstance(value, Mapping):
        raise TypeError("interval must be an object")
    return Interval(cast(int, value.get("start")), cast(int, value.get("end")))


def _location_from_value(value: object) -> Location:
    if not isinstance(value, Mapping):
        raise TypeError("feature location must be an object")
    kind = value.get("kind")
    if kind == "interval":
        return _interval_from_value(value)
    raw_parts = value.get("parts")
    if kind == "compound":
        if not isinstance(raw_parts, list):
            raise TypeError("compound location parts must be an array")
        return CompoundLocation(_interval_from_value(part) for part in raw_parts)
    if kind == "unresolved":
        reason = value.get("reason")
        raw_anchors = value.get("anchors", [])
        if not isinstance(reason, str) or not isinstance(raw_anchors, list):
            raise TypeError("unresolved location requires a reason and anchor array")
        return UnresolvedLocation(
            reason,
            (_interval_from_value(anchor) for anchor in raw_anchors),
        )
    raise TypeError("feature location kind must be interval, compound, or unresolved")


def _features_from_value(value: object, *, record_index: int) -> tuple[DNAFeature, ...]:
    if not isinstance(value, (list, tuple)):
        raise InputFormatError(
            "Structured input features must be an array.",
            code="INVALID_FEATURES_COLUMN",
            context={"record_index": record_index},
        )
    features: list[DNAFeature] = []
    try:
        for raw_feature in value:
            if not isinstance(raw_feature, Mapping):
                raise TypeError("each feature must be an object")
            qualifiers = raw_feature.get("qualifiers", {})
            if not isinstance(qualifiers, Mapping):
                raise TypeError("feature qualifiers must be an object")
            features.append(
                DNAFeature(
                    cast(str, raw_feature.get("type")),
                    _location_from_value(raw_feature.get("location")),
                    id=cast(str | None, raw_feature.get("id")),
                    strand=cast(str, raw_feature.get("strand", "unknown")),
                    label=cast(str | None, raw_feature.get("label")),
                    score=cast(float | None, raw_feature.get("score")),
                    phase=cast(int | None, raw_feature.get("phase")),
                    qualifiers=cast(Mapping[str, object], qualifiers),
                    source=cast(str | None, raw_feature.get("source")),
                )
            )
    except (DNAKitError, TypeError, ValueError) as exc:
        raise InputFormatError(
            "Structured input contains an invalid feature.",
            code="INVALID_FEATURES_COLUMN",
            context={"record_index": record_index},
        ) from exc
    return tuple(features)


def iter_delimited(handle: TextIO, config: ReadConfig, *, delimiter: str) -> Iterator[DNARecord]:
    previous_limit = csv.field_size_limit()
    try:
        csv.field_size_limit(config.max_field_size)
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames
        if fieldnames is None or (
            config.sequence_column not in fieldnames and config.parts_column not in fieldnames
        ):
            raise InputFormatError(
                "Delimited input must have a header containing sequence or parts.",
                code="MISSING_SEQUENCE_COLUMN",
                context={"column": config.sequence_column},
            )
        if any(not isinstance(name, str) or not name for name in fieldnames) or len(
            set(fieldnames)
        ) != len(fieldnames):
            raise InputFormatError(
                "Delimited input header names must be non-empty and unique.",
                code="INVALID_TABLE_HEADER",
                context={"fieldnames": fieldnames},
            )
        for index, row in enumerate(reader):
            if None in row:
                raise InputFormatError(
                    "Delimited input row contains more fields than its header.",
                    code="EXTRA_TABLE_FIELDS",
                    context={"record_index": index, "line_number": reader.line_num},
                )
            yield _record_from_mapping(row, config, record_index=index)
    except csv.Error as exc:
        raise InputFormatError(
            "Delimited input exceeds the configured field-size limit or is malformed.",
            code="CSV_FIELD_ERROR",
            context={"max_field_size": config.max_field_size, "reason": str(exc)},
        ) from exc
    finally:
        csv.field_size_limit(previous_limit)


def iter_json_lines(handle: TextIO, config: ReadConfig) -> Iterator[DNARecord]:
    record_index = 0
    for line_number, line in enumerate(handle, start=1):
        if len(line.rstrip("\r\n")) > config.max_field_size:
            raise _input_error(
                "JSON Lines record exceeds max_field_size.",
                code="STRUCTURED_LINE_TOO_LONG",
                line_number=line_number,
            )
        if not line.strip():
            continue
        value = _strict_json_loads(
            line,
            config,
            message="JSON Lines input contains invalid JSON.",
            code="INVALID_JSON",
            context={"line_number": line_number},
        )
        yield _record_from_mapping(
            _require_object(value, record_index=record_index),
            config,
            record_index=record_index,
        )
        record_index += 1


def iter_json_array(handle: TextIO, config: ReadConfig) -> Iterator[DNARecord]:
    """Parse the standard JSON array representation.

    The parser materializes at most ``max_input_bytes`` and separately bounds
    JSON depth and node count. Use JSONL, CSV, TSV, FASTA, or FASTQ for streaming.
    """

    value = _strict_json_loads(
        handle.read(),
        config,
        message="JSON input contains invalid JSON.",
        code="INVALID_JSON",
    )
    if not isinstance(value, list):
        raise InputFormatError(
            "JSON record input must be an array of objects.",
            code="INVALID_JSON_ROOT",
        )
    for index, item in enumerate(value):
        yield _record_from_mapping(
            _require_object(item, record_index=index), config, record_index=index
        )


def record_to_mapping(record: DNARecord) -> dict[str, Any]:
    """Return the portable MVP record schema used by JSON and delimited I/O."""

    return {
        "schema_version": _RECORD_SCHEMA_VERSION,
        "id": record.id,
        "sequence": record.sequence.symbols,
        "parts": [_sequence_part_to_mapping(part) for part in record.sequence.parts],
        "alphabet": record.sequence.alphabet.value,
        "topology": record.sequence.topology.value,
        "strandedness": record.sequence.strandedness.value,
        "description": record.description,
        "features": [_feature_to_mapping(feature) for feature in record.features],
        "metadata": to_json_compatible(record.metadata),
        "letter_annotations": to_json_compatible(record.letter_annotations),
    }


def _sequence_part_to_mapping(part: str | Gap) -> dict[str, Any]:
    if isinstance(part, str):
        return {"kind": "symbols", "symbols": part}
    return {
        "kind": "gap",
        "length": part.length,
        "gap_kind": part.kind.value,
        "crossable": part.crossable,
        "evidence": list(part.evidence),
        "metadata": to_json_compatible(part.metadata),
    }


def _location_to_mapping(location: Location) -> dict[str, Any]:
    if isinstance(location, Interval):
        return {"kind": "interval", "start": location.start, "end": location.end}
    if isinstance(location, CompoundLocation):
        return {
            "kind": "compound",
            "parts": [{"start": part.start, "end": part.end} for part in location.parts],
        }
    return {
        "kind": "unresolved",
        "reason": location.reason,
        "anchors": [{"start": anchor.start, "end": anchor.end} for anchor in location.anchors],
    }


def _feature_to_mapping(feature: DNAFeature) -> dict[str, Any]:
    return {
        "type": feature.type,
        "location": _location_to_mapping(feature.location),
        "id": feature.id,
        "strand": feature.strand.value,
        "label": feature.label,
        "score": feature.score,
        "phase": feature.phase,
        "qualifiers": to_json_compatible(feature.qualifiers),
        "source": feature.source,
    }


def _reject_implicit_feature_loss(record: DNARecord, *, format: str, config: WriteConfig) -> None:
    if record.features and config.feature_policy == "error":
        raise InputFormatError(
            f"{format} cannot represent DNARecord features without information loss.",
            code="FEATURE_LOSS_NOT_ALLOWED",
            context={"record_id": record.id, "feature_count": len(record.features)},
            hint="Use JSON/JSONL/CSV/TSV or set feature_policy='drop' explicitly.",
        )


def _reject_implicit_gap_loss(record: DNARecord, *, format: str) -> None:
    if record.sequence.is_gapped:
        raise InputFormatError(
            f"{format} cannot represent explicit Gap objects without information loss.",
            code="GAP_LOSS_NOT_ALLOWED",
            context={"record_id": record.id},
            hint="Use JSON, JSONL, CSV, or TSV for lossless Gap persistence.",
        )


def write_fasta(handle: TextIO, records: Iterable[DNARecord], config: WriteConfig) -> int:
    count = 0
    for record in records:
        _reject_implicit_gap_loss(record, format="FASTA")
        _reject_implicit_feature_loss(record, format="FASTA", config=config)
        header = _serialized_header(record, format="FASTA")
        sequence = record.sequence.to_string()
        handle.write(f">{header}\n")
        for offset in range(0, len(sequence), config.line_width):
            handle.write(sequence[offset : offset + config.line_width] + "\n")
        count += 1
    return count


def write_fastq(handle: TextIO, records: Iterable[DNARecord], config: WriteConfig) -> int:
    count = 0
    for record in records:
        _reject_implicit_gap_loss(record, format="FASTQ")
        _reject_implicit_feature_loss(record, format="FASTQ", config=config)
        sequence = record.sequence.to_string()
        values = record.letter_annotations.get(_PHRED_KEY)
        if values is None:
            raise InputFormatError(
                "FASTQ output requires letter_annotations['phred_quality'].",
                code="FASTQ_QUALITY_MISSING",
                context={"record_id": record.id},
            )
        if len(values) != len(sequence):
            raise InputFormatError(
                "FASTQ quality length must equal sequence symbol_length.",
                code="FASTQ_QUALITY_LENGTH_MISMATCH",
                context={"record_id": record.id},
            )
        quality_characters: list[str] = []
        for value in values:
            if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= 93:
                raise InputFormatError(
                    "FASTQ output Phred values must be integers from 0 through 93.",
                    code="FASTQ_QUALITY_OUT_OF_RANGE",
                    context={"record_id": record.id, "value": value},
                )
            codepoint = value + config.phred_offset
            if codepoint > 126:
                raise InputFormatError(
                    "FASTQ quality value cannot be encoded with the configured Phred offset.",
                    code="FASTQ_QUALITY_OUT_OF_RANGE",
                    context={
                        "record_id": record.id,
                        "value": value,
                        "phred_offset": config.phred_offset,
                    },
                )
            quality_characters.append(chr(codepoint))
        header = _serialized_header(record, format="FASTQ")
        handle.write(f"@{header}\n{sequence}\n+\n{''.join(quality_characters)}\n")
        count += 1
    return count


def write_delimited(
    handle: TextIO,
    records: Iterable[DNARecord],
    config: WriteConfig,
    *,
    delimiter: str,
) -> int:
    del config
    fieldnames = [
        "schema_version",
        "id",
        "sequence",
        "parts",
        "alphabet",
        "topology",
        "strandedness",
        "description",
        "features",
        "metadata",
        "letter_annotations",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    count = 0
    for record in records:
        row = record_to_mapping(record)
        row["parts"] = json.dumps(row["parts"], ensure_ascii=False, sort_keys=True)
        row["features"] = json.dumps(row["features"], ensure_ascii=False, sort_keys=True)
        row["metadata"] = json.dumps(row["metadata"], ensure_ascii=False, sort_keys=True)
        row["letter_annotations"] = json.dumps(
            row["letter_annotations"], ensure_ascii=False, sort_keys=True
        )
        writer.writerow(row)
        count += 1
    return count


def write_json_lines(handle: TextIO, records: Iterable[DNARecord], config: WriteConfig) -> int:
    if config.json_indent is not None:
        raise ConfigurationError(
            "JSONL requires one compact JSON object per physical line.",
            code="JSONL_INDENT_NOT_ALLOWED",
            hint="Use format='json' when pretty-printed multi-line output is required.",
        )
    count = 0
    for record in records:
        payload = record_to_mapping(record)
        handle.write(
            json.dumps(
                payload,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
        )
        count += 1
    return count


def write_json_array(handle: TextIO, records: Iterable[DNARecord], config: WriteConfig) -> int:
    handle.write("[\n" if config.json_indent is not None else "[")
    count = 0
    for record in records:
        if count:
            handle.write(",\n" if config.json_indent is not None else ",")
        payload = record_to_mapping(record)
        encoded = json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=None if config.json_indent is not None else (",", ":"),
            indent=config.json_indent,
        )
        if config.json_indent is not None:
            encoded = "\n".join(f"  {line}" for line in encoded.splitlines())
        handle.write(encoded)
        count += 1
    if config.json_indent is not None:
        handle.write("\n")
    handle.write("]\n")
    return count


__all__ = [
    "iter_delimited",
    "iter_fasta",
    "iter_fastq",
    "iter_json_array",
    "iter_json_lines",
    "write_delimited",
    "write_fasta",
    "write_fastq",
    "write_json_array",
    "write_json_lines",
]
