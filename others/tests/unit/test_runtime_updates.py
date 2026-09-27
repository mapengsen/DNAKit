"""Refresh repaired functions while preserving live module state and public aliases."""

from __future__ import annotations

import importlib.util
import inspect
import sys
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import Mock

import pytest

import dnakit
import dnakit.generated as generated
from dnakit._verbosity import with_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.tools.runtime_updates import refresh_repaired_function

_INITIAL = """from __future__ import annotations

class State:
    pass

INSTANCE = State()
CACHE = {}
FACTOR = 2
CALLS = 0

def measure(value: int) -> int:
    return value - 1

def unrelated(value: int) -> int:
    return value + 10
"""

_UPDATED = """from __future__ import annotations

raise RuntimeError("Do not reexecute module initialization")

def measure(value: int, *, offset: int = 0) -> int:
    global CALLS
    CALLS += 1
    CACHE["last"] = value
    return value * FACTOR + offset
"""


def _module(monkeypatch: pytest.MonkeyPatch, path: Path, *, encoding: str = "utf-8") -> ModuleType:
    name = "dnakit._runtime_updates_test"
    path.write_text(_INITIAL, encoding=encoding)
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("wrapped", [False, True])
def test_refresh_updates_aliases_and_optional_arguments_without_reloading_module_state(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    wrapped: bool,
) -> None:
    path = tmp_path / "implementation.py"
    module = _module(monkeypatch, path)
    original = module.measure
    exported = with_verbose(original) if wrapped else original
    vars(module)["measure"] = exported
    aliases = ModuleType("dnakit._runtime_aliases_test")
    vars(aliases)["public_alias"] = exported
    vars(aliases)["raw_alias"] = original
    monkeypatch.setitem(sys.modules, aliases.__name__, aliases)
    monkeypatch.setattr(dnakit, "_runtime_root_alias", exported, raising=False)
    saved = (module.State, module.INSTANCE, module.CACHE, module.unrelated)
    path.write_text(_UPDATED, encoding="utf-8")

    refresh_repaired_function(f"{module.__name__}.measure", str(path), "measure")

    updated = module.measure
    assert saved == (module.State, module.INSTANCE, module.CACHE, module.unrelated)
    assert (
        aliases.public_alias is aliases.raw_alias is vars(dnakit)["_runtime_root_alias"] is updated
    )
    assert "offset" in inspect.signature(updated).parameters
    kwargs: dict[str, Any] = {"offset": 3}
    if wrapped:
        assert updated.__dnakit_verbose_wrapper__ is updated
        kwargs["verbose"] = False
    assert updated(4, **kwargs) == 11
    assert module.CALLS == 1 and module.CACHE["last"] == 4
    vars(module)["FACTOR"] = 5
    assert updated(4, **kwargs) == 23
    assert module.CALLS == 2
    assert inspect.unwrap(updated).__globals__ is vars(module)


@pytest.mark.parametrize("wrapped", [False, True])
def test_bom_source_refresh_and_rollback_update_the_live_function(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, wrapped: bool
) -> None:
    path = tmp_path / "implementation.py"
    module = _module(monkeypatch, path, encoding="utf-8-sig")
    original_bytes = path.read_bytes()
    assert original_bytes.startswith(b"\xef\xbb\xbf")
    original = module.measure
    if wrapped:
        vars(module)["measure"] = with_verbose(original)
    saved = (module.State, module.INSTANCE, module.CACHE)
    path.write_text(_UPDATED, encoding="utf-8-sig")

    refresh_repaired_function(f"{module.__name__}.measure", str(path), "measure")

    updated = module.measure
    assert updated(4, offset=3) == 11
    assert module.CALLS == 1 and module.CACHE["last"] == 4
    assert saved == (module.State, module.INSTANCE, module.CACHE)

    path.write_bytes(original_bytes)
    refresh_repaired_function(f"{module.__name__}.measure", str(path), "measure")

    restored = module.measure
    assert restored is not updated and restored(4) == 3
    assert "offset" not in inspect.signature(restored).parameters
    assert saved == (module.State, module.INSTANCE, module.CACHE)
    assert module.CALLS == 1 and module.CACHE["last"] == 4
    assert path.read_bytes() == original_bytes
    if wrapped:
        assert updated.__dnakit_verbose_wrapper__ is updated
        assert restored.__dnakit_verbose_wrapper__ is restored


def test_failed_definition_refresh_preserves_old_callable_and_state(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    path = tmp_path / "implementation.py"
    module = _module(monkeypatch, path)
    original, state = module.measure, module.INSTANCE
    path.write_text(
        "def measure(value: int, *, offset=missing_default) -> int:\n    return value\n",
        encoding="utf-8",
    )

    with pytest.raises(NameError):
        refresh_repaired_function(f"{module.__name__}.measure", str(path), "measure")

    assert module.measure is original and module.INSTANCE is state
    assert module.measure(4) == 3


def test_refresh_rejects_changed_source_or_function_identity(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    path = tmp_path / "implementation.py"
    module = _module(monkeypatch, path)
    original = module.measure
    other = tmp_path / "other.py"
    other.write_text(_UPDATED, encoding="utf-8")

    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_RUNTIME_INVALID"):
        refresh_repaired_function(f"{module.__name__}.measure", str(other), "measure")
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_RUNTIME_INVALID"):
        refresh_repaired_function(f"{module.__name__}.measure", str(path), "unrelated")
    assert module.measure is original


def test_generated_repair_uses_manifest_refresh_instead_of_native_module_execution(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    refresh = Mock()
    monkeypatch.setattr(generated, "_refresh", refresh)

    refresh_repaired_function(
        "dnakit.generated.generated_probe", "missing-source.py", "generated_probe"
    )

    refresh.assert_called_once_with()
