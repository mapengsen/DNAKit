"""Refresh one repaired function without reloading unrelated module state.

Existing classes, caches and configuration objects keep their identity. New
registry snapshots and module exports use the reviewed function definition.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import linecache
import sys
from pathlib import Path
from typing import Any

from dnakit._verbosity import with_verbose
from dnakit.exceptions import ConfigurationError


def refresh_repaired_function(tool_name: str, source_path: str, function_name: str) -> None:
    """Activate a source-level repair in already imported DNAKit modules."""
    if tool_name.startswith("dnakit.generated."):
        from dnakit.generated import _refresh

        _refresh()
        return
    public_module_name, export_name = tool_name.rsplit(".", 1)
    public_module = importlib.import_module(public_module_name)
    old_export = getattr(public_module, export_name)
    old_function = inspect.unwrap(old_export)
    if not inspect.isfunction(old_function):
        raise ConfigurationError(
            "The repaired export is not a Python function.", code="TOOL_REPAIR_RUNTIME_INVALID"
        )
    if old_function.__name__ != function_name:
        raise ConfigurationError(
            "The repaired definition must retain the original function name.",
            code="TOOL_REPAIR_RUNTIME_INVALID",
        )
    path = Path(source_path).resolve()
    if Path(inspect.getsourcefile(old_function) or "").resolve() != path:
        raise ConfigurationError(
            "The repaired function source location changed.", code="TOOL_REPAIR_RUNTIME_INVALID"
        )
    tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
    definitions = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == function_name
    ]
    if len(definitions) != 1:
        raise ConfigurationError(
            "Cannot locate one repaired definition.", code="TOOL_REPAIR_RUNTIME_INVALID"
        )
    futures = [
        node
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "__future__"
    ]
    selected = ast.Module(body=[*futures, definitions[0]], type_ignores=[])
    # Use the live module namespace: copying it freezes reassigned configuration
    # and sends ``global`` writes into an invisible duplicate module state.
    context = old_function.__globals__
    missing = object()
    previous = context.get(function_name, missing)
    try:
        exec(compile(selected, str(path), "exec", dont_inherit=True), context)
        updated = context[function_name]
        if getattr(old_export, "__dnakit_verbose_wrapper__", None) is old_export:
            updated = with_verbose(updated)
    except BaseException:
        if previous is missing:
            context.pop(function_name, None)
        else:
            context[function_name] = previous
        raise
    # Export modules commonly re-export the same verbosity-wrapped callable.
    replacements: list[tuple[dict[str, Any], str]] = []
    for module_name, module in tuple(sys.modules.items()):
        if module is None or not (module_name == "dnakit" or module_name.startswith("dnakit.")):
            continue
        for name, value in tuple(vars(module).items()):
            if value is old_export or value is old_function:
                replacements.append((vars(module), name))
    for namespace, name in replacements:
        namespace[name] = updated
    setattr(public_module, export_name, updated)
    context[function_name] = updated
    linecache.checkcache(str(path))
