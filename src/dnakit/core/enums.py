"""Stable string enums used by DNAKit core value objects."""

from __future__ import annotations

from enum import Enum


class StringEnum(str, Enum):
    """Python 3.10-compatible string enum."""

    def __str__(self) -> str:
        return str(self.value)


class DNAAlphabet(StringEnum):
    STRICT = "strict"
    IUPAC = "iupac"


class Topology(StringEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


class Strandedness(StringEnum):
    SINGLE = "single"
    DOUBLE = "double"


class Strand(StringEnum):
    FORWARD = "forward"
    REVERSE = "reverse"
    BOTH = "both"
    UNKNOWN = "unknown"


class GapKind(StringEnum):
    UNKNOWN = "unknown"
    SCAFFOLD = "scaffold"
    CONTIG = "contig"
    CENTROMERE = "centromere"
    SHORT_ARM = "short_arm"
    HETEROCHROMATIN = "heterochromatin"
    TELOMERE = "telomere"
    REPEAT = "repeat"
    CONTAMINATION = "contamination"


class CoordinateSystem(StringEnum):
    ZERO_BASED_HALF_OPEN = "0-based-half-open"
    ZERO_BASED_CLOSED = "0-based-closed"
    ONE_BASED_CLOSED = "1-based-closed"
    ONE_BASED_HALF_OPEN = "1-based-half-open"


class IssueSeverity(StringEnum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class ImplementationLabel(StringEnum):
    NATIVE = "native"
    ADAPTER = "adapter"
    REIMPLEMENTATION = "reimplementation"
    NOVEL = "novel"


class ExecutionMode(StringEnum):
    INTERNAL = "internal"
    EXTERNAL = "external"
    HYBRID = "hybrid"


class OriginClass(StringEnum):
    DNAKIT = "dnakit"
    STANDARD = "standard"
    PUBLISHED_ALGORITHM = "published_algorithm"
    INTEGRATION = "integration"
    NOVEL = "novel"


__all__ = [
    "CoordinateSystem",
    "DNAAlphabet",
    "ExecutionMode",
    "GapKind",
    "ImplementationLabel",
    "IssueSeverity",
    "OriginClass",
    "Strand",
    "Strandedness",
    "StringEnum",
    "Topology",
]
