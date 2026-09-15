"""Consistent, context-local feedback for public Python operations."""

from __future__ import annotations

import inspect
import sys
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import is_dataclass, replace
from functools import wraps
from typing import Any

_enabled: ContextVar[bool] = ContextVar("dnakit_verbose", default=True)
_active: ContextVar[Any] = ContextVar("dnakit_operation", default=None)


def progress_enabled(requested: bool = True) -> bool:
    """Whether DNAKit may render progress in the current operation."""
    return requested and _enabled.get()


def status_enabled(requested: bool = True) -> bool:
    """Legacy adapters defer start/end messages to the outer public operation."""
    return progress_enabled(requested) and _active.get() is None


@contextmanager
def _scope(operation: Any) -> Iterator[None]:
    flag = _enabled.set(operation.enabled and _enabled.get())
    parent = _active.set(_active.get() or operation)
    try:
        yield
    finally:
        _active.reset(parent)
        _enabled.reset(flag)


class _Operation:
    def __init__(self, function: Callable[..., Any], arguments: Mapping[str, Any], enabled: bool):
        self.enabled = enabled
        self.visible = enabled and _active.get() is None
        self.done = False
        self.name = function.__name__.replace("_", " ")
        if function.__module__.startswith(("dnakit.download", "dnakit.references")) and (
            "download" in function.__name__
            or (
                function.__module__.startswith("dnakit.download")
                and function.__name__ not in {"build_index", "metadata"}
            )
        ):
            self.name = "download"
        self.target = next(
            (
                arguments[key]
                for key in ("output_path", "output_dir", "target")
                if arguments.get(key) is not None
            ),
            None,
        )

    def start(self) -> None:
        if self.visible:
            suffix = f": {self.target}" if self.target is not None else ": processing"
            print(f"[DNAKit] Starting {self.name}{suffix}", file=sys.stderr, flush=True)

    def finish(self, result: Any = None, *, state: str = "complete") -> None:
        if self.done:
            return
        self.done = True
        if not self.visible:
            return
        details: list[str] = []
        if state == "complete":
            counts = getattr(result, "counts", None)
            if isinstance(counts, Mapping) and all(isinstance(v, int) for v in counts.values()):
                details.append(", ".join(f"{name}={count}" for name, count in counts.items()))
            elif isinstance(getattr(result, "record_count", None), int):
                details.append(f"{result.record_count} records")
            destination = getattr(result, "fasta_path", None) or self.target
            if destination is None:
                artifact = getattr(result, "target_artifact", None)
                destination = getattr(artifact, "relative_path", None)
            if destination is not None:
                details.append(str(destination))
        if state != "complete" and self.target is not None:
            details.append(str(self.target))
        suffix = ": " + "; ".join(details) if details else ""
        print(f"[DNAKit] {self.name.capitalize()} {state}{suffix}", file=sys.stderr, flush=True)


class _OperationIterator(Iterator[Any]):
    """Resume feedback context per next(), never across the caller's yield boundary."""

    def __init__(self, iterator: Iterator[Any], operation: _Operation) -> None:
        self.iterator = iterator
        self.operation = operation

    def __iter__(self) -> _OperationIterator:
        return self

    def __next__(self) -> Any:
        with _scope(self.operation):
            try:
                return next(self.iterator)
            except StopIteration:
                self.operation.finish()
                raise
            except BaseException:
                self.operation.finish(state="failed or interrupted")
                raise

    def close(self) -> None:
        with _scope(self.operation):
            try:
                close = getattr(self.iterator, "close", None)
                if close is not None:
                    close()
            finally:
                self.operation.finish(state="closed")


def _defer_completion(result: Any, operation: _Operation) -> Any:
    # Preserve RecordSource's public type, collect(), ownership, and context manager.
    if type(result).__module__ == "dnakit.io.source" and type(result).__name__ == "RecordSource":
        result._iterator = _OperationIterator(result._iterator, operation)
        return result
    if isinstance(result, Iterator):
        return _OperationIterator(result, operation)
    operation.finish(result)
    return result


def with_verbose(function: Callable[..., Any]) -> Any:
    """Expose verbose=True without altering computation or existing return contracts."""
    if getattr(function, "__dnakit_verbose_wrapper__", None) is function:
        return function
    signature = inspect.signature(function)
    native_verbose = "verbose" in signature.parameters

    def prepare(args: Any, kwargs: Any, caller: Mapping[str, Any]) -> tuple[Any, _Operation]:
        verbose = kwargs.get("verbose", True) if native_verbose else kwargs.pop("verbose", True)
        if not isinstance(verbose, bool):
            from dnakit.exceptions import ConfigurationError

            raise ConfigurationError("verbose must be boolean.")
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        legacy = bound.arguments.get("show_progress", True)
        config = bound.arguments.get("config")
        callbacks = (bound.arguments.get("progress"), bound.arguments.get("query_progress"))
        importing = getattr(caller.get("__spec__"), "_initializing", False)
        enabled = verbose and _enabled.get() and not importing
        # Explicit legacy quiet switches and caller-owned callbacks remain compatible.
        enabled = (
            enabled and legacy is not False and getattr(config, "show_progress", True) is not False
        )
        visible = enabled and not any(callback is not None for callback in callbacks)
        operation = _Operation(function, bound.arguments, enabled)
        operation.visible = operation.visible and visible
        operation.visible = operation.visible and not caller.get("__name__", "").startswith(
            "dnakit."
        )
        if not enabled:
            if "show_progress" in signature.parameters:
                bound.arguments["show_progress"] = False
            if (
                is_dataclass(config)
                and not isinstance(config, type)
                and hasattr(config, "show_progress")
            ):
                bound.arguments["config"] = replace(config, show_progress=False)
        operation.start()
        return bound, operation

    if inspect.iscoroutinefunction(function):

        @wraps(function)
        async def wrapped(*args: Any, **kwargs: Any) -> Any:
            bound, operation = prepare(args, kwargs, sys._getframe(1).f_globals)
            with _scope(operation):
                try:
                    result = await function(*bound.args, **bound.kwargs)
                except BaseException:
                    operation.finish(state="failed or interrupted")
                    raise
                return _defer_completion(result, operation)

    else:

        @wraps(function)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            bound, operation = prepare(args, kwargs, sys._getframe(1).f_globals)
            with _scope(operation):
                try:
                    result = function(*bound.args, **bound.kwargs)
                except BaseException:
                    operation.finish(state="failed or interrupted")
                    raise
                return _defer_completion(result, operation)

    parameters = list(signature.parameters.values())
    if native_verbose:
        parameters = [
            parameter.replace(annotation=bool) if parameter.name == "verbose" else parameter
            for parameter in parameters
        ]
    else:
        insertion = next(
            (
                index
                for index, parameter in enumerate(parameters)
                if parameter.kind == parameter.VAR_KEYWORD
            ),
            len(parameters),
        )
        parameters.insert(
            insertion,
            inspect.Parameter(
                "verbose", inspect.Parameter.KEYWORD_ONLY, default=True, annotation=bool
            ),
        )
    wrapped.__dict__["__signature__"] = signature.replace(parameters=parameters)
    wrapped.__annotations__ = {**getattr(function, "__annotations__", {}), "verbose": bool}
    wrapped.__dict__["__dnakit_verbose_wrapper__"] = wrapped
    for name in ("cache_clear", "cache_info", "cache_parameters"):
        if hasattr(function, name):
            wrapped.__dict__[name] = getattr(function, name)
    return wrapped


def register_public(namespace: dict[str, Any], names: tuple[str, ...] | None = None) -> None:
    """Register exported functions; classes, constants, and private helpers are untouched."""
    for name in namespace.get("__all__", ()) if names is None else names:
        value = namespace.get(name)
        if (
            callable(value)
            and inspect.isfunction(inspect.unwrap(value))
            and value.__module__.startswith("dnakit")
        ):
            namespace[name] = with_verbose(value)
