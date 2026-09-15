"""Private context-local observation of file reads, independent of parsing rules."""

from collections.abc import Callable
from contextvars import ContextVar

read_observer: ContextVar[Callable[[str, int], None] | None] = ContextVar(
    "dnakit_read_observer", default=None
)
