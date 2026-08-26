"""DNAKit exception hierarchy.

Exceptions are reserved for violated object invariants, invalid configuration,
and failed infrastructure.  Recoverable sequence quality findings belong in an
``Issue``-based report instead.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType


class DNAKitError(Exception):
    """Base class carrying a stable code and actionable diagnostic context."""

    default_code = "DNAKIT_ERROR"

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        context: Mapping[str, object] | None = None,
        hint: str | None = None,
    ) -> None:
        if not isinstance(message, str) or not message.strip():
            raise ValueError("An exception message must be a non-empty string.")
        resolved_code = self.default_code if code is None else code
        if not isinstance(resolved_code, str) or not resolved_code.strip():
            raise ValueError("An exception code must be a non-empty string.")
        if hint is not None and (not isinstance(hint, str) or not hint.strip()):
            raise ValueError("An exception hint must be None or a non-empty string.")

        self.message = message
        self.code = resolved_code
        self.context = MappingProxyType(dict(context or {}))
        self.hint = hint
        super().__init__(message)

    @property
    def error_code(self) -> str:
        """Alias retained for callers that prefer the explicit name."""

        return self.code

    def __str__(self) -> str:
        text = f"[{self.code}] {self.message}"
        if self.hint is not None:
            text += f" Hint: {self.hint}"
        return text


class SequenceError(DNAKitError):
    """Base class for invalid sequence objects or unsupported sequence states."""

    default_code = "SEQUENCE_ERROR"


class InvalidAlphabetError(SequenceError):
    """Raised when symbols do not satisfy the declared DNA alphabet."""

    default_code = "INVALID_ALPHABET"


class UnknownLengthError(SequenceError):
    """Raised when an exact length is requested across an unknown-length gap."""

    default_code = "UNKNOWN_LENGTH"


class UnsupportedGapOperationError(SequenceError):
    """Raised when an operation cannot safely process sequence gaps."""

    default_code = "UNSUPPORTED_GAP_OPERATION"


class CoordinateError(DNAKitError):
    """Raised for invalid coordinates or impossible coordinate conversions."""

    default_code = "COORDINATE_ERROR"


class FeatureError(DNAKitError):
    """Raised when a feature violates its structural invariants."""

    default_code = "FEATURE_ERROR"


class DuplicateIDError(DNAKitError):
    """Raised when an operation requires unique record identifiers."""

    default_code = "DUPLICATE_ID"


class InputFormatError(DNAKitError):
    """Raised when serialized input does not match its declared format."""

    default_code = "INPUT_FORMAT_ERROR"


class RecordSourceClosedError(DNAKitError):
    """Raised when consuming a manually closed record source."""

    default_code = "RECORD_SOURCE_CLOSED"


class ConfigurationError(DNAKitError):
    """Raised for invalid or internally inconsistent configuration."""

    default_code = "CONFIGURATION_ERROR"


class CacheError(DNAKitError):
    """Raised for cache integrity or persistence failures."""

    default_code = "CACHE_ERROR"


class BackendError(DNAKitError):
    """Base class for optional backend failures."""

    default_code = "BACKEND_ERROR"


class BackendUnavailableError(BackendError):
    """Raised when a requested optional backend cannot be located."""

    default_code = "BACKEND_UNAVAILABLE"


class BackendVersionError(BackendError):
    """Raised when a backend version is incompatible or unparseable."""

    default_code = "BACKEND_VERSION_ERROR"


class BackendTimeoutError(BackendError):
    """Raised when a backend exceeds its configured timeout."""

    default_code = "BACKEND_TIMEOUT"


class BackendExecutionError(BackendError):
    """Raised when a backend process or library call fails."""

    default_code = "BACKEND_EXECUTION_ERROR"


class DownloadError(BackendError):
    """Raised when an external reference-data download fails."""

    default_code = "DOWNLOAD_ERROR"


class QueryError(BackendError):
    """Raised when a public biological-database query fails."""

    default_code = "QUERY_ERROR"


class BatchExecutionError(DNAKitError):
    """Raised when a batch cannot continue under its error policy."""

    default_code = "BATCH_EXECUTION_ERROR"


__all__ = [
    "BackendError",
    "BackendExecutionError",
    "BackendTimeoutError",
    "BackendUnavailableError",
    "BackendVersionError",
    "BatchExecutionError",
    "CacheError",
    "ConfigurationError",
    "CoordinateError",
    "DNAKitError",
    "DownloadError",
    "DuplicateIDError",
    "FeatureError",
    "InputFormatError",
    "InvalidAlphabetError",
    "QueryError",
    "RecordSourceClosedError",
    "SequenceError",
    "UnknownLengthError",
    "UnsupportedGapOperationError",
]
