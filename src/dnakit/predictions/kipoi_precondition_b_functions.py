"""Kipoi model operations that run inside Singularity or in the current Python environment (kipoi 0.8.6).

Three operations are offered:

- ``predict_kipoi_model_in_singularity``: the ``kipoi predict --singularity`` path. It runs ``kipoi predict`` for a
  local model directory inside an explicit local Singularity image (a ``.sif`` file), with bind mounts derived from the
  model directory, the dataloader's local file arguments and the output file.
- ``test_kipoi_model_in_singularity``: the ``kipoi test`` path of ``kipoi test-source --singularity`` for a local model
  directory inside an explicit local Singularity image.
- ``test_kipoi_model_in_current_env``: ``kipoi test`` (kipoi's own ``cli_test``) for a local model directory in the
  running Python environment. The check fails when kipoi logs a WARNING, as ``kipoi.cli.source_test.test_model`` does
  for its subprocess output. No conda environment is created.

The singularity runtime is an external program: the user installs Singularity or Apptainer and DNAKit calls ``singularity
exec``. A missing program raises ``BackendUnavailableError`` with an install hint. Images are never downloaded here;
kipoi's remote container lookup is not used. ``dry_run=True`` returns the command line without running it.

Prediction output is written as TSV, as in ``predict_kipoi_model``: kipoi 0.8.6 cannot write HDF5 with NumPy 2.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import shutil
from collections.abc import Mapping, Sequence
from io import StringIO
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_KIPOI_HINT = "Install kipoi with: python -m pip install kipoi"
_SINGULARITY_HINT = (
    "Install Singularity or Apptainer (for example: conda install -c conda-forge singularity) "
    "and make sure the singularity executable is on PATH."
)
__all__ = [
    "predict_kipoi_model_in_singularity",
    "test_kipoi_model_in_current_env",
    "test_kipoi_model_in_singularity",
]

_TSV = ".tsv"
_SOURCE = "dir"
_PREDICT_TIMEOUT = 3600.0
_TEST_TIMEOUT = 1800.0
_LOG_TAIL_CHARS = 4000


def _kipoi() -> Any:
    try:
        import kipoi
        import kipoi.cli.main
    except ImportError as exc:
        raise BackendUnavailableError(
            f"kipoi could not be imported: {exc}", hint=_KIPOI_HINT
        ) from exc
    return kipoi


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value.strip()


def _count(value: object, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _existing_dir(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def _existing_file(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def _new_tsv(value: object, name: str) -> Path:
    path = Path(_name(value, name)).expanduser().resolve()
    if path.suffix.lower() != _TSV:
        raise ConfigurationError(f"{name} must end with {_TSV}: only TSV output is offered.")
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def _mapping(value: object, name: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of argument names to values.")
    return {str(key): item for key, item in value.items()}


def _bind_directories(paths: Sequence[Path], arguments: Mapping[str, Any]) -> list[str]:
    """Directories that must be visible inside the container: the given paths and the folders of local file arguments."""
    directories: set[str] = set()
    for path in paths:
        directories.add(str(path if path.is_dir() else path.parent))
    for value in arguments.values():
        if isinstance(value, str) and os.path.exists(os.path.expanduser(value)):
            candidate = Path(os.path.expanduser(value)).resolve()
            directories.add(str(candidate if candidate.is_dir() else candidate.parent))
    return sorted(directories)


def _singularity() -> str:
    found = shutil.which("singularity")
    if found is None:
        raise BackendUnavailableError(
            "The singularity executable was not found on PATH.", hint=_SINGULARITY_HINT
        )
    return found


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


def _singularity_version(executable: str) -> str | None:
    result = execute_bounded_command(
        executable, ["--version"], backend_id="singularity", cwd=os.getcwd(), verbose=False
    )
    if result.return_code != 0:
        return None
    return result.output.strip() or None


def _run_in_container(
    *,
    operation: str,
    container: Path,
    inside: list[str],
    bind_dirs: list[str],
    outputs: list[Path],
    timeout: float,
    dry_run: bool,
) -> dict[str, Any]:
    """Build ``singularity exec -B <dirs> <image> <kipoi command>``; run it unless ``dry_run`` is set."""
    binds = ["-B", ",".join(bind_dirs)] if bind_dirs else []
    if dry_run:
        return {
            "status": "dry_run",
            "operation": operation,
            "container_image": str(container),
            "command": ["singularity", "exec", *binds, str(container), *inside],
            "bind_directories": bind_dirs,
            "outputs": [],
            "return_code": None,
            "log_tail": "",
            "provenance": {"backend": "singularity", "version": None},
        }
    executable = _singularity()
    version = _singularity_version(executable)
    arguments = ["exec", *binds, str(container), *inside]
    cwd = os.getcwd()
    result = execute_bounded_command(
        executable,
        arguments,
        backend_id="singularity",
        cwd=cwd,
        timeout_seconds=timeout,
        verbose=False,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"{operation}: singularity exec exited with code {result.return_code}: {_tail(result.output)}"
        )
    written = [str(path) for path in outputs if path.is_file()]
    return {
        "status": "completed",
        "operation": operation,
        "container_image": str(container),
        "command": [executable, *arguments],
        "bind_directories": bind_dirs,
        "outputs": written,
        "return_code": result.return_code,
        "log_tail": _tail(result.output),
        "provenance": {"backend": "singularity", "version": version},
    }


def predict_kipoi_model_in_singularity(
    container_image: str,
    model_dir: str,
    output_file: str,
    *,
    dataloader_args: Mapping[str, Any] | None = None,
    batch_size: int = 32,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run ``kipoi predict`` for a local model directory inside an explicit Singularity image and write TSV predictions.

    This is the ``kipoi predict --singularity`` operation with the container given explicitly: ``container_image`` is a
    local ``.sif`` file that contains kipoi. The command is ``singularity exec -B <dirs> <image> kipoi predict <model_dir>
    --source dir -o <output_file>``; ``dataloader_args`` is passed to ``--dataloader_args`` as JSON, and local file
    arguments are bind-mounted. ``output_file`` must not exist yet and must end with ``.tsv``. A missing singularity
    executable raises ``BackendUnavailableError``; with ``dry_run=True`` the command is returned and nothing runs.
    """
    image = _existing_file(container_image, "container_image")
    directory = _existing_dir(model_dir, "model_dir")
    out = _new_tsv(output_file, "output_file")
    size = _count(batch_size, "batch_size", 1)
    arguments = _mapping(dataloader_args, "dataloader_args")
    inside = [
        "kipoi",
        "predict",
        str(directory),
        "--source",
        _SOURCE,
        "--batch_size",
        str(size),
        "-o",
        str(out),
    ]
    if arguments:
        inside += ["--dataloader_args", json.dumps(arguments)]
    bind_dirs = _bind_directories([directory, out.parent], arguments)
    return _run_in_container(
        operation="predict_kipoi_model_in_singularity",
        container=image,
        inside=inside,
        bind_dirs=bind_dirs,
        outputs=[out],
        timeout=_PREDICT_TIMEOUT,
        dry_run=_flag(dry_run, "dry_run"),
    )


def test_kipoi_model_in_singularity(
    container_image: str,
    model_dir: str,
    *,
    batch_size: int = 32,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run ``kipoi test`` on a local model directory inside an explicit Singularity image.

    The command is ``singularity exec -B <model dir> <image> kipoi test --batch_size <n> --source dir <model_dir>``, the
    same ``kipoi test`` that ``kipoi test-source --singularity`` runs for each model. kipoi's own test compares the
    predictions with ``test.expect`` in ``model.yaml`` when that file is local. A non-zero exit raises
    ``BackendExecutionError`` with the log tail; a missing singularity executable raises ``BackendUnavailableError``.
    """
    image = _existing_file(container_image, "container_image")
    directory = _existing_dir(model_dir, "model_dir")
    size = _count(batch_size, "batch_size", 1)
    inside = ["kipoi", "test", "--batch_size", str(size), "--source", _SOURCE, str(directory)]
    return _run_in_container(
        operation="test_kipoi_model_in_singularity",
        container=image,
        inside=inside,
        bind_dirs=_bind_directories([directory], {}),
        outputs=[],
        timeout=_TEST_TIMEOUT,
        dry_run=_flag(dry_run, "dry_run"),
    )


class _WarningCollector(logging.Handler):
    """Collects the messages of WARNING and higher records emitted while a kipoi test runs."""

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def test_kipoi_model_in_current_env(
    model_dir: str, *, batch_size: int = 32, verbose: bool = True
) -> dict[str, Any]:
    """Run ``kipoi test`` on a local model directory in the current Python environment and fail on logged warnings.

    kipoi's own ``cli_test`` runs with ``--source dir --batch_size <n>``: the example inputs from ``dataloader.yaml`` are
    predicted, and the expected predictions in ``model.yaml`` (``test.expect``, when it is a local file) are compared
    with the model's precision. As ``kipoi.cli.source_test.test_model`` does for its subprocess output, any WARNING record
    logged during the test makes the call raise ``BackendExecutionError``. Nothing is installed and no conda environment is
    created; the model's dependencies must already be importable.
    """
    kipoi = _kipoi()
    directory = _existing_dir(model_dir, "model_dir")
    size = _count(batch_size, "batch_size", 1)
    collector = _WarningCollector()
    root = logging.getLogger()
    root.addHandler(collector)
    captured = StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            kipoi.cli.main.cli_test(
                "test", ["--batch_size", str(size), "--source", _SOURCE, str(directory)]
            )
    except SystemExit as exc:
        if exc.code not in (None, 0):
            raise BackendExecutionError(
                f"kipoi test exited with code {exc.code}: {_tail(captured.getvalue())}"
            ) from exc
    except Exception as exc:
        raise BackendExecutionError(f"kipoi test failed for {directory}: {exc}") from exc
    finally:
        root.removeHandler(collector)
    if collector.messages:
        raise BackendExecutionError(
            f"kipoi test logged {len(collector.messages)} warning(s) for {directory}: "
            + " | ".join(collector.messages[:5])
        )
    return {
        "status": "completed",
        "operation": "test_kipoi_model_in_current_env",
        "model_dir": str(directory),
        "batch_size": size,
        "warnings": 0,
        "log_tail": _tail(captured.getvalue()),
        "provenance": {"backend": "kipoi", "version": getattr(kipoi, "__version__", None)},
    }


_register_public_verbose(globals())
