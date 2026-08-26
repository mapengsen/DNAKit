"""Single-consumption record source lifecycle."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import AbstractContextManager
from enum import Enum
from types import TracebackType

from dnakit.core.collection import DNASet
from dnakit.core.record import DNARecord
from dnakit.exceptions import RecordSourceClosedError


class _SourceState(Enum):
    ACTIVE = "active"
    EXHAUSTED = "exhausted"
    CLOSED = "closed"


class RecordSource(Iterator[DNARecord], AbstractContextManager["RecordSource"]):
    """A lazy, single-use iterator that deterministically releases resources.

    Natural exhaustion remains ordinary iterator exhaustion.  Manual closure is
    distinct and any later attempt to consume records raises
    :class:`~dnakit.exceptions.RecordSourceClosedError`.
    """

    def __init__(
        self,
        iterator: Iterator[DNARecord],
        *,
        close_callback: Callable[[], None] | None = None,
        source_name: str | None = None,
        format: str | None = None,
    ) -> None:
        self._iterator = iterator
        self._close_callback = close_callback
        self._source_name = source_name
        self._format = format
        self._state = _SourceState.ACTIVE
        self._resource_released = False

    @property
    def source_name(self) -> str | None:
        return self._source_name

    @property
    def format(self) -> str | None:
        return self._format

    @property
    def closed(self) -> bool:
        """Whether underlying resources have been released."""

        return self._state is not _SourceState.ACTIVE

    @property
    def exhausted(self) -> bool:
        return self._state is _SourceState.EXHAUSTED

    def __iter__(self) -> RecordSource:
        return self

    def __next__(self) -> DNARecord:
        if self._state is _SourceState.CLOSED:
            raise RecordSourceClosedError(
                "Cannot consume a manually closed RecordSource.",
                context={"source": self._source_name, "format": self._format},
                hint="Call read() again to create a new single-use source.",
            )
        if self._state is _SourceState.EXHAUSTED:
            raise StopIteration
        try:
            return next(self._iterator)
        except StopIteration:
            self._state = _SourceState.EXHAUSTED
            self._release_resource()
            raise
        except BaseException:
            self._state = _SourceState.CLOSED
            self._release_resource()
            raise

    def collect(self) -> DNASet:
        """Materialize only the unconsumed records, then release resources."""

        if self._state is _SourceState.CLOSED:
            raise RecordSourceClosedError(
                "Cannot collect from a manually closed RecordSource.",
                context={"source": self._source_name, "format": self._format},
                hint="Call read() again to create a new source.",
            )
        try:
            records = tuple(self)
        finally:
            if self._state is _SourceState.ACTIVE:
                self.close()
        return DNASet.from_records(records, source=self._source_name)

    def close(self) -> None:
        """Close an active source; repeated calls are safe."""

        if self._state is _SourceState.ACTIVE:
            self._state = _SourceState.CLOSED
            self._release_resource()

    def _release_resource(self) -> None:
        if self._resource_released:
            return
        self._resource_released = True
        iterator_close = getattr(self._iterator, "close", None)
        try:
            if callable(iterator_close):
                iterator_close()
        finally:
            if self._close_callback is not None:
                self._close_callback()

    def __enter__(self) -> RecordSource:
        if self._state is _SourceState.CLOSED:
            raise RecordSourceClosedError(
                "Cannot enter a context with a manually closed RecordSource.",
                context={"source": self._source_name, "format": self._format},
            )
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        del exc_type, exc_value, traceback
        self.close()


__all__ = ["RecordSource"]
