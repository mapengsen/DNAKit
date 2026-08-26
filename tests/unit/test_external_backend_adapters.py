"""Tests for passive, non-redistributed external CLI adapter metadata."""

from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path

import pytest

from dnakit.backends import EXTERNAL_CLI_ADAPTERS, backend_registry, execute_bounded_command
from dnakit.backends.external import ExternalCLIAdapter
from dnakit.exceptions import BackendExecutionError, BackendTimeoutError, ConfigurationError


def test_external_cli_adapters_are_registered_without_execution(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    calls: list[str] = []

    def locate(executable: str) -> None:
        calls.append(executable)
        return None

    monkeypatch.setattr(shutil, "which", locate)
    expected = {"blastn", "mmseqs2", "sourmash", "dashing"}

    assert expected <= set(backend_registry)
    for adapter in EXTERNAL_CLI_ADAPTERS:
        info = backend_registry.probe(adapter.backend_id)
        assert info.available is False
        assert info.version is None
        assert info.metadata["probe_mode"] == "path-only"
        assert info.metadata["redistributed"] is False
    assert calls == [adapter.executable for adapter in EXTERNAL_CLI_ADAPTERS]


def test_dashing_adapter_declares_external_gpl_boundary() -> None:
    adapter = next(item for item in EXTERNAL_CLI_ADAPTERS if item.backend_id == "dashing")
    info = adapter.probe()

    assert info.license_expression == "GPL-3.0-only"
    assert info.package_location is None
    assert info.metadata["execution"] == "explicit opt-in only"


def _fake_adapter(path: Path) -> ExternalCLIAdapter:
    return ExternalCLIAdapter("fake", str(path), frozenset({"test"}), "MIT")


def test_explicit_version_query_parses_fake_executable(tmp_path: Path) -> None:
    executable = tmp_path / "fake-version"
    executable.write_text("#!/bin/sh\necho 'fake 1.2.3'\n", encoding="utf-8")
    executable.chmod(0o700)

    info = _fake_adapter(executable).version()

    assert info.available is True
    assert info.version == "1.2.3"
    assert info.metadata["probe_mode"] == "explicit-version-command"


def test_explicit_version_query_reports_nonzero_exit(tmp_path: Path) -> None:
    executable = tmp_path / "fake-failure"
    executable.write_text("#!/bin/sh\necho 'fake 2.0.0'\nexit 3\n", encoding="utf-8")
    executable.chmod(0o700)

    info = _fake_adapter(executable).version()

    assert info.available is False
    assert info.version == "2.0.0"
    assert info.metadata["return_code"] == 3


def test_explicit_version_query_bounds_time_and_output(tmp_path: Path) -> None:
    sleepy = tmp_path / "sleepy"
    sleepy.write_text("#!/bin/sh\nsleep 1\n", encoding="utf-8")
    sleepy.chmod(0o700)
    noisy = tmp_path / "noisy"
    noisy.write_text("#!/bin/sh\nprintf '%0200d' 1\n", encoding="utf-8")
    noisy.chmod(0o700)

    with pytest.raises(BackendTimeoutError):
        _fake_adapter(sleepy).version(timeout_seconds=0.01)
    with pytest.raises(BackendExecutionError) as output_error:
        _fake_adapter(noisy).version(max_output_bytes=10)
    assert output_error.value.code == "BACKEND_VERSION_OUTPUT_LIMIT"
    observed = output_error.value.context["observed_bytes"]
    assert isinstance(observed, int) and observed <= 11


def test_explicit_version_timeout_terminates_descendants(tmp_path: Path) -> None:
    executable = tmp_path / "descendant"
    executable.write_text("#!/bin/sh\nsleep 3 &\nwait\n", encoding="utf-8")
    executable.chmod(0o700)

    started = time.monotonic()
    with pytest.raises(BackendTimeoutError):
        _fake_adapter(executable).version(timeout_seconds=0.01)
    assert time.monotonic() - started < 1.0


def test_explicit_version_query_wraps_spawn_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    executable = tmp_path / "exists"
    executable.write_text("ignored", encoding="utf-8")
    executable.chmod(0o700)

    def fail(*args: object, **kwargs: object) -> subprocess.Popen[bytes]:
        del args, kwargs
        raise OSError("cannot execute")

    monkeypatch.setattr(subprocess, "Popen", fail)
    with pytest.raises(BackendExecutionError) as error:
        _fake_adapter(executable).version()
    assert error.value.code == "BACKEND_VERSION_EXECUTION_FAILED"


def test_bounded_executor_rejects_unbounded_arguments_before_execution(tmp_path: Path) -> None:
    executable = tmp_path / "must-not-run"
    executable.write_text("#!/bin/sh\nexit 99\n", encoding="utf-8")
    executable.chmod(0o700)

    with pytest.raises(ConfigurationError) as error:
        execute_bounded_command(
            executable,
            (f"argument-{index}" for index in range(300)),
            backend_id="fake",
            cwd=tmp_path,
        )

    assert error.value.code == "EXTERNAL_ARGUMENT_LIMIT_EXCEEDED"


def test_bounded_executor_confines_monitored_outputs_to_cwd(tmp_path: Path) -> None:
    executable = tmp_path / "must-not-run"
    executable.write_text("#!/bin/sh\nexit 99\n", encoding="utf-8")
    executable.chmod(0o700)

    with pytest.raises(ConfigurationError) as error:
        execute_bounded_command(
            executable,
            ("fixed-subcommand",),
            backend_id="fake",
            cwd=tmp_path,
            monitored_output_paths=(tmp_path.parent / "outside.tsv",),
        )

    assert error.value.code == "UNSAFE_EXTERNAL_OUTPUT_PATH"
