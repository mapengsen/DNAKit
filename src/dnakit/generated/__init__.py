"""Validated generated tools, exported from the local activation manifest.

Only activated functions appear in ``__all__``. Internal refresh replaces the
exported set after generation or rollback; it never installs dependencies.
"""

from __future__ import annotations

import builtins as _builtins
import hashlib as _hashlib
import json as _json
import keyword as _keyword
import re as _re
import sys as _sys
import threading as _threading
import types as _types
from pathlib import Path as _Path
from typing import Any as _Any

from dnakit._verbosity import with_verbose as _with_verbose
from dnakit.exceptions import ConfigurationError as _ConfigurationError

_DIRECTORY = _Path(__file__).resolve().parent
_LOCK = _threading.RLock()
_MODULES: dict[str, _types.ModuleType] = {}
__all__: list[str] = []


def _entries(directory: _Path | None = None) -> dict[str, dict[str, _Any]]:
    """Read metadata without executing generated modules."""
    manifest = (directory or _DIRECTORY) / "manifest.json"
    if not manifest.exists():
        return {}
    try:
        with manifest.open("rb") as stream:
            raw = stream.read(4_000_001)
        if len(raw) > 4_000_000:
            raise ValueError("Manifest exceeds 4 MB.")
        data = _json.loads(raw)
        if not isinstance(data, dict) or data.get("version") != 1:
            raise ValueError("Unsupported manifest version.")
        entries = data.get("tools")
        if not isinstance(entries, dict):
            raise ValueError("Manifest tools must be an object.")
        for name, entry in entries.items():
            if not isinstance(name, str) or not _re.fullmatch(r"[a-z][a-z0-9_]{0,79}", name):
                raise ValueError("Invalid generated function name.")
            if _keyword.iskeyword(name) or hasattr(_builtins, name):
                raise ValueError("Generated function names cannot shadow Python builtins.")
            if not isinstance(entry, dict):
                raise ValueError("Invalid generated entry.")
            digest = entry.get("sha256")
            if not isinstance(digest, str) or not _re.fullmatch(r"[0-9a-f]{64}", digest):
                raise ValueError("Invalid module digest.")
            if entry.get("module") != f"{name}_{digest[:16]}.py":
                raise ValueError("Invalid generated module path.")
            contract = entry.get("contract")
            if not isinstance(contract, dict) or not isinstance(contract.get("effects"), list):
                raise ValueError("Generated tools must include an effect contract.")
            effects = contract["effects"]
            if not effects or any(
                effect not in {"local", "network", "filesystem", "model", "external_process"}
                for effect in effects
            ):
                raise ValueError("Invalid effect contract.")
            order = ("local", "network", "filesystem", "model", "external_process")
            if entry.get("effect") != max(effects, key=order.index):
                raise ValueError(
                    "The manifest effect must match the most conservative declared effect."
                )
            if "local" in effects and len(set(effects)) != 1:
                raise ValueError("A local effect cannot be mixed with side effects.")
        return entries
    except (OSError, UnicodeError, ValueError, TypeError, RecursionError) as exc:
        raise _ConfigurationError(
            "The generated-tool manifest is invalid.",
            code="GENERATED_MANIFEST_INVALID",
            context={"path": str(manifest), "detail": str(exc)},
        ) from None


def _refresh() -> None:
    """Reload validated code and atomically replace public exports after activation/rollback."""
    with _LOCK:
        entries = _entries()
        exports: dict[str, _Any] = {}
        modules: dict[str, _types.ModuleType] = {}
        for name, entry in entries.items():
            path = _DIRECTORY / entry["module"]
            if path.is_symlink() or path.resolve().parent != _DIRECTORY.resolve():
                raise _ConfigurationError(
                    "Generated module paths must stay inside their package.",
                    code="GENERATED_MODULE_INVALID",
                )
            try:
                with path.open("rb") as stream:
                    source = stream.read(200_001)
            except OSError:
                raise _ConfigurationError(
                    "An activated generated module is missing or unreadable.",
                    code="GENERATED_MODULE_INVALID",
                    context={"path": str(path)},
                ) from None
            digest = _hashlib.sha256(source).hexdigest()
            if len(source) > 200_000 or digest != entry["sha256"]:
                raise _ConfigurationError(
                    "An activated generated module failed its content hash check.",
                    code="GENERATED_MODULE_INVALID",
                    context={"path": str(path)},
                )
            key = f"{path.resolve()}:{digest}"
            module = _MODULES.get(key)
            if module is None:
                module_name = f"dnakit.generated.{path.stem}"
                module = _types.ModuleType(module_name)
                module.__file__ = str(path)
                module.__package__ = "dnakit.generated"
                previous = _sys.modules.get(module_name)
                _sys.modules[module_name] = module
                try:
                    exec(compile(source, str(path), "exec"), module.__dict__)
                except BaseException:
                    if previous is None:
                        _sys.modules.pop(module_name, None)
                    else:
                        _sys.modules[module_name] = previous
                    raise
            function = getattr(module, name, None)
            if not isinstance(function, _types.FunctionType):
                raise _ConfigurationError(
                    "An activated generated module does not expose its declared function.",
                    code="GENERATED_MODULE_INVALID",
                    context={"name": name},
                )
            function = _with_verbose(function)
            setattr(module, name, function)
            exports[name] = function
            modules[key] = module
        for name in __all__:
            globals().pop(name, None)
        globals().update(exports)
        __all__[:] = sorted(exports)
        for key, old_module in _MODULES.items():
            if key not in modules and _sys.modules.get(old_module.__name__) is old_module:
                _sys.modules.pop(old_module.__name__, None)
        _MODULES.clear()
        _MODULES.update(modules)


_refresh()
