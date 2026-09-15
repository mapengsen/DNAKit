"""Bounded, ordered CPU processing without changing scientific results."""

from __future__ import annotations

import copyreg
import importlib
import inspect
import io
import multiprocessing
import os
import pickle
from collections import deque
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from types import MappingProxyType
from typing import Any

from dnakit.exceptions import ConfigurationError

_workers: ContextVar[int] = ContextVar("dnakit_workers", default=-1)


def available_cpus() -> int:
    process_count = getattr(os, "process_cpu_count", None)
    if process_count is not None:
        return process_count() or 1
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return os.cpu_count() or 1


def resolve_workers(workers: int) -> int:
    if isinstance(workers, bool) or not isinstance(workers, int) or workers == 0 or workers < -1:
        raise ConfigurationError("workers must be -1 (all available CPUs) or a positive integer.")
    return available_cpus() if workers == -1 else min(workers, available_cpus())


@contextmanager
def worker_scope(workers: int) -> Iterator[None]:
    resolve_workers(workers)
    token = _workers.set(_workers.get() if workers == -1 else workers)
    try:
        yield
    finally:
        _workers.reset(token)


def current_workers() -> int:
    return resolve_workers(_workers.get())


def available_memory() -> int:
    """Account for Linux container limits when available."""
    available = 512 * 1024**2
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
                break
        limit = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if limit != "max":
            used = int(Path("/sys/fs/cgroup/memory.current").read_text())
            available = min(available, max(0, int(limit) - used))
    except (OSError, ValueError):
        pass
    return max(1, available)


def memory_budget() -> int:
    return max(1, min(256 * 1024**2, available_memory() // 8))


def _proxy(values: dict[Any, Any]) -> Any:
    return MappingProxyType(values)


def _reduce_proxy(value: Any) -> Any:
    return _proxy, (dict(value),)


def _reduce_frozen(value: Any) -> Any:
    return type(value), (dict(value),)


def pack(value: Any) -> bytes:
    # Local reducers: do not alter pickle behaviour globally in the user's application.
    from dnakit.core._json import FrozenDict

    stream = io.BytesIO()
    encoder = pickle.Pickler(stream, protocol=pickle.HIGHEST_PROTOCOL)
    encoder.dispatch_table = {
        **copyreg.dispatch_table,
        MappingProxyType: _reduce_proxy,
        FrozenDict: _reduce_frozen,
    }
    encoder.dump(value)
    return stream.getvalue()


def unpack(value: bytes) -> Any:
    # Only private, locally generated IPC/spool data reaches this function.
    return pickle.loads(value)


def _function(module: str, name: str) -> Any:
    return inspect.unwrap(getattr(importlib.import_module(module), name))


def _execute(module: str, name: str, payloads: list[bytes]) -> bytes:
    from dnakit._verbosity import _enabled

    token = _enabled.set(False)
    try:
        with worker_scope(1):
            function = _function(module, name)
            try:
                return pack((True, [function(unpack(payload)) for payload in payloads]))
            except Exception as exc:
                return pack((False, exc))
    finally:
        _enabled.reset(token)


def ordered_map(
    function: Any, items: Iterable[Any], *, workers: int | None = None, batch_size: int = 64
) -> Iterator[Any]:
    """Map importable unary work, limiting both queued bytes and worker count.

    Tiny workloads and oversized records run locally. Results are yielded in input
    order; worker scheduling never controls random seeds or merge order.
    """
    count = current_workers() if workers is None else resolve_workers(workers)
    if count == 1:
        yield from map(function, items)
        return
    budget = memory_budget()
    # Allow at least ~128 MiB of free memory per child in addition to queued payloads.
    count = min(count, max(1, available_memory() // (128 * 1024**2)))
    if count == 1:
        yield from map(function, items)
        return
    iterator = iter(items)
    pending: deque[Any] = deque()
    queued = 0
    pool = None
    batch: list[bytes] = []
    size = 0

    def consume() -> Iterator[Any]:
        nonlocal queued
        future, amount = pending.popleft()
        queued -= amount
        success, value = unpack(future.result())
        if not success:
            raise value
        yield from value

    def submit() -> Iterator[Any]:
        nonlocal pool, queued, batch, size
        while pending and (len(pending) >= count or queued + size > budget):
            yield from consume()
        if pool is None:
            pool = ProcessPoolExecutor(
                max_workers=count, mp_context=multiprocessing.get_context("spawn")
            )
        pending.append((pool.submit(_execute, function.__module__, function.__name__, batch), size))
        queued += size
        batch, size = [], 0

    try:
        for item in iterator:
            payload = pack(item)
            if len(payload) > budget // count:
                if batch:
                    yield from submit()
                while pending:
                    yield from consume()
                yield function(item)
                continue
            batch.append(payload)
            size += len(payload)
            if len(batch) >= batch_size or size >= min(1024**2, budget // count):
                yield from submit()
        if batch:
            if pool is None:
                for payload in batch:
                    yield function(unpack(payload))
            else:
                yield from submit()
        while pending:
            yield from consume()
    finally:
        for future, _ in pending:
            future.cancel()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)


def invoke(payload: Any) -> Any:
    """Invoke a native public operation by importable name; used by record adapters."""
    module, name, args, kwargs, materialize = payload
    function = _function(module, name)
    parameters = inspect.signature(function).parameters
    kwargs = {
        key: value
        for key, value in kwargs.items()
        if key not in {"verbose", "workers", "show_progress"} or key in parameters
    }
    result = function(*args, **kwargs)
    return list(result) if materialize else result


def _direction_block(payload: tuple[str, bool, bool]) -> str:
    from dnakit.core.sequence import _COMPLEMENT

    text, reverse, complement = payload
    text = text[::-1] if reverse else text
    return text.translate(_COMPLEMENT) if complement else text


def transform_text(text: str, *, reverse: bool, complement: bool) -> str:
    """Use C-level string kernels; copying these strings through IPC costs more.

    Parallel file transforms distribute file offsets instead, outside this kernel.
    """
    return _direction_block((text, reverse, complement))
