"""Private helpers shared by the Picard command-line facades of integration pack P08.

A facade runs one fixed Picard tool with an explicit ``KEY=value`` argument list (no shell) through
``dnakit.backends.external.execute_bounded_command``. Picard is found from an explicit
``picard_path`` (a ``picard`` launcher script or a ``picard.jar``), then from ``picard`` on PATH,
then from the ``PICARD_JAR`` environment variable. A jar is started with ``java -jar``; ``java``
must be on PATH. A missing program raises ``BackendUnavailableError`` carrying an install hint.
Results are JSON-serialisable dictionaries: status, backend, tool, version, command, return code,
existing outputs and a bounded log tail.
"""

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DNAKitError

TESTED_VERSION = "3.5.0"
LOG_TAIL_CHARS = 4000
DEFAULT_TIMEOUT_SECONDS = 3600.0
MAX_OUTPUT_BYTES = 5_000_000
MAX_ARTIFACT_BYTES = 2_000_000_000
INSTALL_HINT = (
    "Install Picard 3.x (for example `conda install -c bioconda picard`, or download "
    "picard.jar from "
    "https://github.com/broadinstitute/picard/releases and install a Java 17 runtime), "
    "then put `picard` on PATH, "
    "pass picard_path (a launcher or a picard.jar) or set PICARD_JAR."
)
STRINGENCY = ("STRICT", "LENIENT", "SILENT")
_VERSION = re.compile(r"Version:\s*(\d+(?:\.\d+)+)")
_PLAIN = re.compile(r"[^\x00-\x1f\x7f]+")


# ---------------------------------------------------------------------------
# Locating and running Picard
# ---------------------------------------------------------------------------


def launcher(picard_path: str | None) -> tuple[str, tuple[str, ...]]:
    """The executable and its leading arguments (``-jar picard.jar`` for a jar)."""
    jar: Path | None = None
    if picard_path is not None:
        candidate = Path(picard_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(
                f"picard_path does not exist: {candidate}", hint=INSTALL_HINT
            )
        if candidate.suffix == ".jar":
            jar = candidate.resolve()
        else:
            return str(candidate.resolve()), ()
    else:
        found = shutil.which("picard")
        if found is not None:
            return str(Path(found).resolve()), ()
        env_jar = os.environ.get("PICARD_JAR")
        if env_jar and Path(env_jar).expanduser().is_file():
            jar = Path(env_jar).expanduser().resolve()
    if jar is None:
        raise BackendUnavailableError(
            "Picard was not found (no `picard` on PATH and no PICARD_JAR).", hint=INSTALL_HINT
        )
    java = shutil.which("java")
    if java is None:
        raise BackendUnavailableError(
            "Java was not found on PATH (needed to run picard.jar).", hint=INSTALL_HINT
        )
    return str(Path(java).resolve()), ("-jar", str(jar))


@lru_cache(maxsize=16)
def _version_of(executable: str, prefix: tuple[str, ...]) -> str | None:
    try:
        result = execute_bounded_command(
            executable,
            [*prefix, "CheckTerminatorBlock", "--version"],
            backend_id="picard",
            cwd=Path(executable).parent,
            timeout_seconds=120.0,
            max_output_bytes=100_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    match = _VERSION.search(result.output or "")
    return match.group(1) if match else None


def run(
    tool: str,
    arguments: Sequence[str],
    *,
    outputs: Mapping[str, Path],
    picard_path: str | None,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    ok_codes: Sequence[int] = (0,),
    cwd: Path | None = None,
    required: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Run one Picard tool. ``status`` is completed when the exit code is accepted and the required
    outputs exist.

    ``outputs`` maps a result key to a file or directory path; every key in ``required`` (default:
    all keys) must exist after the run. ``ok_codes`` lists the exit codes that are results rather
    than failures.
    """
    executable, prefix = launcher(picard_path)
    if cwd is None:
        parents = [str(path.parent) for path in outputs.values()]
        cwd = Path(os.path.commonpath(parents)) if parents else Path.cwd()
    monitored = [
        path
        for path in outputs.values()
        if path != cwd and cwd in path.parents and not path.is_dir()
    ]
    result = execute_bounded_command(
        executable,
        [*prefix, tool, *arguments],
        backend_id="picard",
        cwd=cwd,
        timeout_seconds=timeout_seconds,
        max_output_bytes=MAX_OUTPUT_BYTES,
        monitored_output_paths=monitored[:16],
        max_monitored_output_bytes=MAX_ARTIFACT_BYTES,
        verbose=False,
    )
    must = list(outputs) if required is None else list(required)
    present = {key: str(path) for key, path in outputs.items() if path.exists()}
    complete = result.return_code in ok_codes and all(key in present for key in must)
    return {
        "status": "completed" if complete else "failed",
        "backend": "picard",
        "tool": tool,
        "version": _version_of(executable, prefix),
        "tested_version": TESTED_VERSION,
        "command": ["picard", tool, *arguments],
        "return_code": result.return_code,
        "outputs": present,
        "log_tail": (result.output or "")[-LOG_TAIL_CHARS:],
    }


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------


def existing_file(value: str | Path, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def existing_files(values: Sequence[str | Path], name: str) -> list[Path]:
    if isinstance(values, (str, Path)) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of existing files.")
    return [existing_file(value, name) for value in values]


def existing_directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def optional_file(value: str | Path | None, name: str) -> Path | None:
    return None if value is None else existing_file(value, name)


def new_file(value: str | Path, name: str) -> Path:
    """A path that does not exist yet and whose folder exists; facades never overwrite outputs."""
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def optional_new_file(value: str | Path | None, name: str) -> Path | None:
    return None if value is None else new_file(value, name)


def new_directory(value: str | Path, name: str) -> Path:
    """A directory path that does not exist yet (Picard creates it); its parent must exist."""
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or _PLAIN.fullmatch(value) is None:
        raise ConfigurationError(f"{name} must be a non-empty string without control characters.")
    return value


def integer(
    value: object, name: str, *, minimum: int | None = None, maximum: int | None = None
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def number(
    value: object, name: str, *, minimum: float | None = None, maximum: float | None = None
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return float(value)


def boolean(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def choice(value: object, options: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {', '.join(options)}; got {value!r}.")
    return value


def choices(values: Sequence[str] | None, options: Sequence[str], name: str) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        raise ConfigurationError(f"{name} must be a list of strings.")
    return [choice(value, options, name) for value in values]


def texts(values: Sequence[str] | None, name: str) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        raise ConfigurationError(f"{name} must be a list of strings.")
    return [text(value, name) for value in values]


def stringency(value: str | None) -> str | None:
    return None if value is None else choice(value, STRINGENCY, "validation_stringency")


# ---------------------------------------------------------------------------
# Argument building
# ---------------------------------------------------------------------------


def _render(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def kv(key: str, value: object | None) -> list[str]:
    """``KEY=value`` for a given value, nothing for ``None`` (so Picard's own default applies)."""
    return [] if value is None else [f"{key}={_render(value)}"]


def kv_many(key: str, values: Sequence[object] | None) -> list[str]:
    """``KEY=a KEY=b`` for a list-valued Picard argument."""
    return [] if not values else [f"{key}={_render(value)}" for value in values]


# ---------------------------------------------------------------------------
# Metrics files
# ---------------------------------------------------------------------------


def _scalar(value: str) -> object:
    if value == "":
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def read_metrics(path: Path) -> dict[str, Any]:
    """Parse a Picard metrics file: the metrics class, its rows and the optional histogram rows."""
    sections: dict[str, tuple[str, list[list[str]]]] = {}
    current: str | None = None
    label = ""
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("## METRICS CLASS") or line.startswith("## HISTOGRAM"):
            if current is not None:
                sections[current] = (label, rows)
            current = "metrics" if "METRICS CLASS" in line else "histogram"
            label = line.split("\t", 1)[1] if "\t" in line else ""
            rows = []
        elif current is not None and not line.startswith("#"):
            if line.strip():
                rows.append(line.split("\t"))
    if current is not None:
        sections[current] = (label, rows)

    def table(name: str) -> list[dict[str, object]]:
        if name not in sections or not sections[name][1]:
            return []
        header, *body = sections[name][1]
        return [{h: _scalar(v) for h, v in zip(header, row, strict=False)} for row in body]

    return {
        "metrics_class": sections.get("metrics", ("", []))[0] or None,
        "metrics": table("metrics"),
        "histogram": table("histogram"),
    }


def count_lines(path: Path, prefix: str | None = None, *, skip_prefix: str | None = None) -> int:
    """Number of lines in a text file, optionally only those starting with ``prefix``.

    Lines starting with ``skip_prefix`` are not counted.
    """
    total = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if prefix is not None and not line.startswith(prefix):
                continue
            if skip_prefix is not None and line.startswith(skip_prefix):
                continue
            total += 1
    return total


# ---------------------------------------------------------------------------
# Optional-argument shortcuts: validate when given, pass nothing when ``None``
# ---------------------------------------------------------------------------


def a_int(
    key: str,
    value: int | None,
    name: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> list[str]:
    return [] if value is None else kv(key, integer(value, name, minimum=minimum, maximum=maximum))


def a_num(
    key: str,
    value: float | None,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> list[str]:
    return [] if value is None else kv(key, number(value, name, minimum=minimum, maximum=maximum))


def a_bool(key: str, value: bool | None, name: str) -> list[str]:
    return [] if value is None else kv(key, boolean(value, name))


def a_choice(key: str, value: str | None, options: Sequence[str], name: str) -> list[str]:
    return [] if value is None else kv(key, choice(value, options, name))


def a_choices(
    key: str, values: Sequence[str] | None, options: Sequence[str], name: str
) -> list[str]:
    return kv_many(key, choices(values, options, name))


def a_text(key: str, value: str | None, name: str) -> list[str]:
    return [] if value is None else kv(key, text(value, name))


def a_texts(key: str, values: Sequence[str] | None, name: str) -> list[str]:
    return kv_many(key, texts(values, name))


def a_file(key: str, value: str | Path | None, name: str) -> list[str]:
    return [] if value is None else kv(key, existing_file(value, name))


def a_files(key: str, values: Sequence[str | Path] | None, name: str) -> list[str]:
    return [] if not values else [f"{key}={path}" for path in existing_files(values, name)]
