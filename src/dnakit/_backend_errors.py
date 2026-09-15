"""Consistent model failure diagnostics without importing a tensor framework."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit.exceptions import BackendExecutionError


def backend_failure(
    message: str, *, code: str, error: Exception, context: Mapping[str, object]
) -> BackendExecutionError:
    cause: BaseException | None = error
    for _ in range(10):
        if cause is None:
            break
        if (
            isinstance(cause, MemoryError)
            or "out of memory" in str(cause).lower()
            or type(cause).__name__ == "OutOfMemoryError"
        ):
            return BackendExecutionError(
                "The model ran out of memory; no prediction result was returned.",
                code="MODEL_OUT_OF_MEMORY",
                context={**context, "error_type": type(cause).__name__},
                hint=(
                    "Reduce batch_size or input length, free device memory, "
                    "or select device='cpu' when supported."
                ),
            )
        cause = cause.__cause__
    if isinstance(error, BackendExecutionError):
        return BackendExecutionError(
            error.message,
            code=error.code,
            context={**context, **error.context},
            hint=error.hint,
        )
    return BackendExecutionError(
        message, code=code, context={**context, "error_type": type(error).__name__}
    )
