"""Bounded retries for transient read-only network operations.

Transfers restart in staging files; this is not HTTP Range resume support.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar

from dnakit.exceptions import DownloadError, QueryError

P = ParamSpec("P")
T = TypeVar("T")
_TRANSIENT_STATUS = {408, 429, 500, 502, 503, 504}
_NETWORK_CODES = {
    "DOWNLOAD_NETWORK_ERROR",
    "NCBI_NETWORK_ERROR",
    "QUERY_NETWORK_ERROR",
    "DOWNLOAD_LENGTH_MISMATCH",
}


def retry_transient(function: Callable[P, T]) -> Callable[P, T]:
    """Try at most three times; never retry validation, checksum or HTTP 404 errors."""

    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> T:
        if kwargs.get("method", "GET") != "GET":
            return function(*args, **kwargs)
        for attempt in range(3):
            try:
                return function(*args, **kwargs)
            except (DownloadError, QueryError) as exc:
                retryable = exc.code in _NETWORK_CODES or (
                    exc.code in {"DOWNLOAD_HTTP_ERROR", "NCBI_HTTP_ERROR", "QUERY_HTTP_ERROR"}
                    and exc.context.get("status") in _TRANSIENT_STATUS
                )
                if not retryable or isinstance(exc.__cause__, PermissionError):
                    raise
                if attempt == 2:
                    raise type(exc)(
                        exc.message,
                        code=exc.code,
                        context={**exc.context, "attempts": attempt + 1},
                        hint=exc.hint
                        or (
                            "Failed after 3 attempts. Retry when the network or server recovers; "
                            "partial transfers restart from the beginning."
                        ),
                    ) from exc
                time.sleep(0.25 * 2**attempt)
        raise AssertionError("unreachable")

    return wrapped
