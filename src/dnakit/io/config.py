"""Configuration value objects for DNA record I/O."""

from __future__ import annotations

import codecs
import sys
from dataclasses import dataclass
from typing import Literal

from dnakit.core.enums import DNAAlphabet, Strandedness, Topology
from dnakit.exceptions import ConfigurationError

CompressionMode = Literal["auto", "none", "gzip"]
AnonymousIDPolicy = Literal["generate", "error"]
FeatureWritePolicy = Literal["error", "drop"]


def _non_empty(value: str, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"{field_name} must be a non-empty string.",
            code="INVALID_IO_CONFIG",
            context={"field": field_name},
        )


def _validate_encoding(value: str) -> None:
    _non_empty(value, "encoding")
    try:
        codecs.lookup(value)
    except LookupError as exc:
        raise ConfigurationError(
            "encoding must name a registered Python text codec.",
            code="INVALID_IO_ENCODING",
            context={"encoding": value},
        ) from exc


@dataclass(frozen=True, slots=True)
class ReadConfig:
    """Control decoding without silently relaxing sequence validation.

    ``close_source=None`` applies ownership-aware behavior: paths are owned and
    closed by :class:`~dnakit.io.RecordSource`, while caller-provided streams
    are borrowed and left open.
    """

    encoding: str = "utf-8"
    alphabet: DNAAlphabet = DNAAlphabet.IUPAC
    topology: Topology = Topology.LINEAR
    strandedness: Strandedness = Strandedness.SINGLE
    uppercase: bool = True
    phred_offset: int = 33
    compression: CompressionMode = "auto"
    close_source: bool | None = None
    delimiter: str | None = None
    id_column: str = "id"
    sequence_column: str = "sequence"
    description_column: str = "description"
    metadata_column: str = "metadata"
    letter_annotations_column: str = "letter_annotations"
    features_column: str = "features"
    parts_column: str = "parts"
    max_field_size: int = 10_000_000
    max_records: int = 1_000_000
    max_record_lines: int = 1_000_000
    max_sequence_symbols: int = 100_000_000
    max_input_bytes: int = 1_000_000_000
    max_json_depth: int = 100
    max_json_nodes: int = 1_000_000
    anonymous_id_prefix: str = "sequence_"

    def __post_init__(self) -> None:
        try:
            alphabet = (
                self.alphabet
                if isinstance(self.alphabet, DNAAlphabet)
                else DNAAlphabet(self.alphabet)
            )
            topology = (
                self.topology if isinstance(self.topology, Topology) else Topology(self.topology)
            )
            strandedness = (
                self.strandedness
                if isinstance(self.strandedness, Strandedness)
                else Strandedness(self.strandedness)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Unknown DNA alphabet, topology, or strandedness for input.",
                code="INVALID_IO_SEQUENCE_TYPE",
                context={
                    "alphabet": self.alphabet,
                    "topology": self.topology,
                    "strandedness": self.strandedness,
                },
            ) from exc
        if self.compression not in {"auto", "none", "gzip"}:
            raise ConfigurationError(
                "compression must be 'auto', 'none', or 'gzip'.",
                code="INVALID_IO_COMPRESSION",
            )
        if (
            isinstance(self.phred_offset, bool)
            or not isinstance(self.phred_offset, int)
            or not 33 <= self.phred_offset <= 64
        ):
            raise ConfigurationError(
                "phred_offset must be an integer between 33 and 64.",
                code="INVALID_PHRED_OFFSET",
            )
        if self.delimiter is not None and (
            not isinstance(self.delimiter, str)
            or len(self.delimiter) != 1
            or self.delimiter in {"\x00", "\r", "\n", '"'}
        ):
            raise ConfigurationError(
                "delimiter must be None or one character distinct from "
                "NUL, line breaks, and the quote character.",
                code="INVALID_TABLE_DELIMITER",
            )
        if not isinstance(self.uppercase, bool):
            raise ConfigurationError(
                "uppercase must be a boolean.",
                code="INVALID_IO_CONFIG",
                context={"field": "uppercase"},
            )
        if self.close_source is not None and not isinstance(self.close_source, bool):
            raise ConfigurationError(
                "close_source must be a boolean or None.",
                code="INVALID_IO_CONFIG",
                context={"field": "close_source"},
            )
        for field_name in (
            "encoding",
            "id_column",
            "sequence_column",
            "description_column",
            "metadata_column",
            "letter_annotations_column",
            "features_column",
            "parts_column",
            "anonymous_id_prefix",
        ):
            _non_empty(getattr(self, field_name), field_name)
        columns = (
            self.id_column,
            self.sequence_column,
            self.description_column,
            self.metadata_column,
            self.letter_annotations_column,
            self.features_column,
            self.parts_column,
        )
        if len(set(columns)) != len(columns):
            raise ConfigurationError(
                "Configured table column names must be unique.",
                code="DUPLICATE_TABLE_COLUMN",
            )
        _validate_encoding(self.encoding)
        if any(symbol.isspace() for symbol in self.anonymous_id_prefix):
            raise ConfigurationError(
                "anonymous_id_prefix must not contain whitespace.",
                code="INVALID_ANONYMOUS_ID_PREFIX",
            )
        if (
            isinstance(self.max_field_size, bool)
            or not isinstance(self.max_field_size, int)
            or not 1 <= self.max_field_size <= sys.maxsize
        ):
            raise ConfigurationError(
                "max_field_size must be an integer between 1 and sys.maxsize.",
                code="INVALID_CSV_FIELD_SIZE_LIMIT",
            )
        for field_name in (
            "max_records",
            "max_record_lines",
            "max_sequence_symbols",
            "max_input_bytes",
            "max_json_depth",
            "max_json_nodes",
        ):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ConfigurationError(
                    f"{field_name} must be a positive integer.",
                    code="INVALID_IO_RESOURCE_LIMIT",
                    context={"field": field_name, "value": value},
                )
        object.__setattr__(self, "alphabet", alphabet)
        object.__setattr__(self, "topology", topology)
        object.__setattr__(self, "strandedness", strandedness)


@dataclass(frozen=True, slots=True)
class WriteConfig:
    """Control deterministic record serialization and destination ownership."""

    encoding: str = "utf-8"
    overwrite: bool = False
    create_parents: bool = False
    compression: CompressionMode = "auto"
    compression_level: int = 6
    close_target: bool | None = None
    line_width: int = 80
    phred_offset: int = 33
    delimiter: str | None = None
    anonymous_id_policy: AnonymousIDPolicy = "generate"
    anonymous_id_prefix: str = "sequence_"
    feature_policy: FeatureWritePolicy = "error"
    json_indent: int | None = None
    max_output_bytes: int | None = None

    def __post_init__(self) -> None:
        if self.compression not in {"auto", "none", "gzip"}:
            raise ConfigurationError(
                "compression must be 'auto', 'none', or 'gzip'.",
                code="INVALID_IO_COMPRESSION",
            )
        if self.anonymous_id_policy not in {"generate", "error"}:
            raise ConfigurationError(
                "anonymous_id_policy must be 'generate' or 'error'.",
                code="INVALID_ANONYMOUS_ID_POLICY",
            )
        if self.feature_policy not in {"error", "drop"}:
            raise ConfigurationError(
                "feature_policy must be 'error' or 'drop'.",
                code="INVALID_FEATURE_WRITE_POLICY",
            )
        if (
            isinstance(self.compression_level, bool)
            or not isinstance(self.compression_level, int)
            or not 0 <= self.compression_level <= 9
        ):
            raise ConfigurationError(
                "compression_level must be an integer between 0 and 9.",
                code="INVALID_COMPRESSION_LEVEL",
            )
        if (
            isinstance(self.line_width, bool)
            or not isinstance(self.line_width, int)
            or self.line_width < 1
        ):
            raise ConfigurationError(
                "line_width must be a positive integer.",
                code="INVALID_FASTA_LINE_WIDTH",
            )
        if (
            isinstance(self.phred_offset, bool)
            or not isinstance(self.phred_offset, int)
            or not 33 <= self.phred_offset <= 64
        ):
            raise ConfigurationError(
                "phred_offset must be an integer between 33 and 64.",
                code="INVALID_PHRED_OFFSET",
            )
        if self.delimiter is not None and (
            not isinstance(self.delimiter, str)
            or len(self.delimiter) != 1
            or self.delimiter in {"\x00", "\r", "\n", '"'}
        ):
            raise ConfigurationError(
                "delimiter must be None or one character distinct from "
                "NUL, line breaks, and the quote character.",
                code="INVALID_TABLE_DELIMITER",
            )
        for field_name in ("overwrite", "create_parents"):
            value = getattr(self, field_name)
            if not isinstance(value, bool):
                raise ConfigurationError(
                    f"{field_name} must be a boolean.",
                    code="INVALID_IO_CONFIG",
                    context={"field": field_name},
                )
        if self.close_target is not None and not isinstance(self.close_target, bool):
            raise ConfigurationError(
                "close_target must be a boolean or None.",
                code="INVALID_IO_CONFIG",
                context={"field": "close_target"},
            )
        if self.json_indent is not None and (
            isinstance(self.json_indent, bool)
            or not isinstance(self.json_indent, int)
            or self.json_indent < 0
        ):
            raise ConfigurationError(
                "json_indent must be a non-negative integer or None.",
                code="INVALID_JSON_INDENT",
            )
        if self.max_output_bytes is not None and (
            isinstance(self.max_output_bytes, bool)
            or not isinstance(self.max_output_bytes, int)
            or self.max_output_bytes < 1
        ):
            raise ConfigurationError(
                "max_output_bytes must be a positive integer or None.",
                code="INVALID_IO_RESOURCE_LIMIT",
                context={"field": "max_output_bytes", "value": self.max_output_bytes},
            )
        _validate_encoding(self.encoding)
        _non_empty(self.anonymous_id_prefix, "anonymous_id_prefix")
        if any(symbol.isspace() for symbol in self.anonymous_id_prefix):
            raise ConfigurationError(
                "anonymous_id_prefix must not contain whitespace.",
                code="INVALID_ANONYMOUS_ID_PREFIX",
            )


__all__ = [
    "AnonymousIDPolicy",
    "CompressionMode",
    "FeatureWritePolicy",
    "ReadConfig",
    "WriteConfig",
]
