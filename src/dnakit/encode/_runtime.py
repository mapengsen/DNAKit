"""Isolated synchronous execution of DNAKit's ENCODE adapters."""

from __future__ import annotations

import asyncio
import importlib.util
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
from pathlib import Path


def invoke(operation: str, arguments: dict[str, object]) -> str:
    """Give each call its own async client state and close all opened resources."""
    path = Path(__file__).parent / "_backend/server/main.py"
    spec = importlib.util.spec_from_file_location("dnakit.encode._invocation", path)
    if spec is None or spec.loader is None:
        raise ImportError("DNAKit ENCODE backend is missing")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    async def execute() -> str:
        async with module.lifespan():
            result = await getattr(module, "encode_" + operation)(**arguments)
            if not isinstance(result, str):
                raise TypeError("ENCODE operation must return JSON text.")
            return result

    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(execute())
    with ThreadPoolExecutor(max_workers=1, thread_name_prefix="dnakit-encode") as pool:
        return pool.submit(copy_context().run, lambda: asyncio.run(execute())).result()
