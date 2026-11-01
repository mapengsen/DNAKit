"""Private helpers shared by the GATK command-line facades of pack P00.

One GATK tool is run with an explicit argument list (no shell) through
``dnakit.backends.external.execute_bounded_command``. The launcher is resolved from an explicit
path, the ``DNAKIT_GATK`` environment variable or PATH. Beyond the typed parameters of a
function, a caller may pass ``options``: a mapping of long option names to values that is
checked against the options that ``gatk <Tool> --help`` lists for that tool (table in
``_gatk_p00_options``); options that a parameter of the function already sets are refused. A
missing launcher raises ``BackendUnavailableError`` with an install hint, a failed run raises
``BackendExecutionError`` with the tail of the log.
"""

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from ._gatk_p00_options import TOOL_OPTIONS

PROGRAM = "gatk"
TESTED_VERSION = "4.6.2.0"
CITATION = "https://doi.org/10.1101/gr.107524.110"
INSTALL_HINT = (
    "Install GATK 4.6.2.0 or later (for example `conda install -c bioconda gatk4`, or the release "
    "zip from https://github.com/broadinstitute/gatk/releases) together with Java 17, then put the "
    "`gatk` launcher on PATH, set DNAKIT_GATK or pass gatk_path."
)
TIMEOUT_SECONDS = 86_400.0
MAX_OUTPUT_BYTES = 5_000_000
MAX_ARTIFACT_BYTES = 2_000_000_000
LOG_TAIL_CHARS = 4000
_VERSION = re.compile(r"v?(\d+\.\d+\.\d+(?:\.\d+)?)")
_MEMORY = re.compile(r"^[1-9]\d{0,5}[mMgG]$")
_OPTION_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
_VERSION_CACHE: dict[str, str | None] = {}


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def text(value: object, name: str) -> str:
    if isinstance(value, os.PathLike):
        value = os.fspath(value)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if value.startswith("-") or "\x00" in value:
        raise ConfigurationError(f"{name} must not start with '-' or contain NUL characters.")
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


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def choice(value: object, options: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {list(options)}; got {value!r}.")
    return value


def existing_file(value: object, name: str) -> Path:
    path = Path(text(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def existing_files(values: object, name: str) -> tuple[Path, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of file paths.")
    return tuple(existing_file(value, name) for value in values)


def existing_path(value: object, name: str) -> Path:
    """An existing file or directory (some GATK inputs accept either)."""
    path = Path(text(value, name)).expanduser().resolve()
    if not path.exists():
        raise ConfigurationError(f"{name} must exist: {path}")
    return path


def existing_directory(value: object, name: str) -> Path:
    path = Path(text(value, name)).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def new_file(value: object, name: str) -> Path:
    """A path that does not exist yet and whose folder exists; nothing is overwritten."""
    path = Path(text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def new_directory(value: object, name: str) -> Path:
    """A directory that does not exist yet and whose parent exists; the facade creates it."""
    path = Path(text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be reused: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The parent folder of {name} does not exist: {path.parent}")
    return path


def new_prefix(value: object, name: str) -> Path:
    """A file-name prefix (not a file): its folder must exist; GATK appends its own suffixes."""
    path = Path(text(value, name)).expanduser().resolve()
    if path.is_dir():
        raise ConfigurationError(f"{name} is a folder; give a file-name prefix inside it: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def texts(values: object, name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    return tuple(text(value, name) for value in values)


def integers(values: object, name: str, *, minimum: int | None = None) -> tuple[int, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of integers.")
    return tuple(integer(value, name, minimum=minimum) for value in values)


def choices(values: object, options: Sequence[str], name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    return tuple(choice(value, options, name) for value in values)


def interval(value: object, name: str) -> str:
    """A locus such as ``chr1:100-200`` or an interval file (resolved to an absolute path)."""
    raw = text(value, name)
    candidate = Path(raw).expanduser()
    return str(candidate.resolve()) if candidate.is_file() else raw


def intervals(values: object, name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of intervals or interval files.")
    return tuple(interval(value, name) for value in values)


def working_directory(expected: Sequence[Path], outputs: Mapping[str, Path]) -> Path:
    """The folder that holds the outputs (monitored below the cwd); the current folder if none."""
    paths = [*expected, *(path for path in outputs.values() if path not in expected)]
    if not paths:
        return Path.cwd()
    return common_directory(paths)


def common_directory(paths: Sequence[Path]) -> Path:
    """The folder that contains every expected output (outputs are monitored below the cwd)."""
    folders = [str(path if path.is_dir() else path.parent) for path in paths]
    return Path(os.path.commonpath(folders))


# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------


def _scalar(value: object, name: str) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, str):
        if not value or "\x00" in value:
            raise ConfigurationError(f"option {name!r} has an empty or NUL-containing value.")
        if value.startswith("-") and not re.fullmatch(r"-\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", value):
            raise ConfigurationError(
                f"option {name!r}: values must not look like command-line flags."
            )
        return value
    raise ConfigurationError(f"option {name!r} must be a string, number, bool or a list of them.")


def _option_value(tool: str, name: str, kind: str, value: object) -> str:
    """One option value; paths become absolute, ``out`` paths must be new files."""
    if kind in {"path", "out", "dir"} and isinstance(value, (str, Path)):
        raw = str(value)
        if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", raw):
            return _scalar(raw, name)
        path = Path(text(raw, name)).expanduser().resolve()
        if kind == "out":
            if path.exists():
                raise ConfigurationError(
                    f"option {name!r} must be a new file, not an existing one: {path}"
                )
            if not path.parent.is_dir():
                raise ConfigurationError(
                    f"The folder for option {name!r} does not exist: {path.parent}"
                )
        return str(path)
    return _scalar(value, name)


def option_arguments(
    tool: str, options: Mapping[str, Any] | None, *, reserved: Sequence[str] = ()
) -> list[str]:
    """Turn ``{"MAX_RECORDS_IN_RAM": 5000}`` into ``["--MAX_RECORDS_IN_RAM", "5000"]``.

    Only the options listed for the tool are accepted. A list value repeats the option;
    ``True`` / ``False`` become ``true`` / ``false``. Options that the calling function already
    sets (``reserved``) are refused; output-file options must name a new file.
    """
    if options is None:
        return []
    if not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping of GATK long option names to values.")
    allowed = TOOL_OPTIONS[tool]
    arguments: list[str] = []
    for raw_name, value in options.items():
        name = raw_name[2:] if isinstance(raw_name, str) and raw_name.startswith("--") else raw_name
        if not isinstance(name, str) or not _OPTION_NAME.fullmatch(name):
            raise ConfigurationError(f"invalid option name {raw_name!r}.")
        if name in reserved:
            raise ConfigurationError(
                f"option {name!r} is set by a parameter of this function; "
                "use that parameter instead."
            )
        if name not in allowed:
            raise ConfigurationError(f"gatk {tool} option {name!r} is not allowed by this facade.")
        values = value if isinstance(value, (list, tuple)) else [value]
        if not values:
            raise ConfigurationError(f"option {name!r} has an empty list.")
        for item in values:
            arguments += [f"--{name}", _option_value(tool, name, allowed[name], item)]
    return arguments


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def resolve_gatk(gatk_path: str | None) -> str:
    """Explicit path, then ``DNAKIT_GATK``, then PATH."""
    candidate_text = gatk_path if gatk_path else os.environ.get("DNAKIT_GATK") or None
    if candidate_text is not None:
        if os.sep not in candidate_text:
            found = shutil.which(candidate_text)
            if found is None:
                raise BackendUnavailableError(
                    f"External executable {candidate_text!r} was not found on PATH.",
                    hint=INSTALL_HINT,
                    code="GATK_UNAVAILABLE",
                )
            return str(Path(found).absolute())
        candidate = Path(candidate_text).expanduser()
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise BackendUnavailableError(
                f"GATK launcher is missing or not executable: {candidate_text}",
                hint=INSTALL_HINT,
                code="GATK_UNAVAILABLE",
            )
        return str(candidate.absolute())
    found = shutil.which(PROGRAM)
    if found is None:
        raise BackendUnavailableError(
            "External executable 'gatk' was not found on PATH.",
            hint=INSTALL_HINT,
            code="GATK_UNAVAILABLE",
        )
    return str(Path(found).absolute())


def _version(executable: str, cwd: Path, verbose: bool) -> str | None:
    if executable in _VERSION_CACHE:
        return _VERSION_CACHE[executable]
    result = execute_bounded_command(
        executable,
        ("--version",),
        backend_id=PROGRAM,
        cwd=cwd,
        timeout_seconds=300.0,
        max_output_bytes=MAX_OUTPUT_BYTES,
        verbose=verbose,
    )
    found: str | None = None
    for line in result.output.splitlines():
        if "GATK" in line or "Genome Analysis Toolkit" in line:
            match = _VERSION.search(line)
            if match:
                found = match.group(1)
                break
    _VERSION_CACHE[executable] = found
    return found


def tail(value: str) -> str:
    return value[-LOG_TAIL_CHARS:]


def run_tool(
    operation: str,
    tool: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    gatk_path: str | None,
    java_max_memory: str | None,
    verbose: bool,
    timeout_seconds: float = TIMEOUT_SECONDS,
    success_codes: Sequence[int] = (0,),
) -> dict[str, Any]:
    """Run ``gatk [--java-options -Xmx..] <tool> <arguments>``; ``expected`` files must exist after.

    ``expected`` paths must be inside ``cwd`` (the executor monitors them there).
    ``success_codes`` lists the exit codes that are not errors (a few GATK checkers report a
    finding with a non-zero code).
    """
    if tool not in TOOL_OPTIONS:
        raise ConfigurationError(f"unknown GATK tool {tool!r}.")
    launcher: list[str] = []
    if java_max_memory is not None:
        if not isinstance(java_max_memory, str) or not _MEMORY.fullmatch(java_max_memory):
            raise ConfigurationError("java_max_memory must look like '4g' or '2048m'.")
        launcher = ["--java-options", f"-Xmx{java_max_memory}"]
    executable = resolve_gatk(gatk_path)
    version = _version(executable, cwd, verbose)
    command = (*launcher, tool, *arguments)
    result = execute_bounded_command(
        executable,
        command,
        backend_id=PROGRAM,
        cwd=cwd,
        timeout_seconds=timeout_seconds,
        max_output_bytes=MAX_OUTPUT_BYTES,
        monitored_output_paths=tuple(expected),
        max_monitored_output_bytes=MAX_ARTIFACT_BYTES,
        verbose=verbose,
    )
    if result.return_code not in success_codes:
        raise BackendExecutionError(
            f"gatk {tool} exited with code {result.return_code}.",
            code="GATK_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": tail(result.output)},
        )
    missing = [str(path) for path in expected if not path.exists()]
    if missing:
        raise BackendExecutionError(
            f"gatk {tool} did not write the expected files: {missing}",
            code="GATK_OUTPUT_MISSING",
            context={"operation": operation, "log_tail": tail(result.output)},
        )
    return {
        "status": "completed",
        "backend": PROGRAM,
        "operation": operation,
        "tool": tool,
        "version": version,
        "command": [PROGRAM, *command],
        "return_code": result.return_code,
        "log_tail": tail(result.output),
        "provenance": {
            "tool": "GATK",
            "version": version,
            "tested_version": TESTED_VERSION,
            "executable": executable,
            "citation": CITATION,
        },
        "_log": result.output,
    }


def finish(
    result: dict[str, Any], *, outputs: Mapping[str, Path | None], **extra: object
) -> dict[str, Any]:
    """Drop the private log, add the output paths and the extra fields."""
    result.pop("_log", None)
    result["outputs"] = {name: str(path) for name, path in outputs.items() if path is not None}
    result.update(extra)
    return result


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------


def read_text(path: Path, limit: int = 2_000_000) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return handle.read(limit)


def _convert(value: str) -> object:
    for kind in (int, float):
        try:
            return kind(value)
        except ValueError:
            continue
    return value


def picard_metrics(path: Path) -> dict[str, Any]:
    """Parse a Picard metrics file: the ``## METRICS CLASS`` rows and the ``## HISTOGRAM`` rows."""
    sections: dict[str, list[list[str]]] = {"metrics": [], "histogram": []}
    current: str | None = None
    for line in read_text(path).splitlines():
        if line.startswith("## METRICS CLASS"):
            current = "metrics"
        elif line.startswith("## HISTOGRAM"):
            current = "histogram"
        elif not line.strip():
            if current is not None and sections[current]:
                current = None
        elif line.startswith("#"):
            continue
        elif current is not None:
            sections[current].append(line.split("\t"))
    parsed: dict[str, Any] = {}
    for key, rows in sections.items():
        parsed[key] = (
            [
                {name: _convert(value) for name, value in zip(rows[0], row, strict=False)}
                for row in rows[1:]
            ]
            if rows
            else []
        )
    return parsed


def tsv_rows(path: Path, *, comment: str = "#", limit: int = 1000) -> list[dict[str, object]]:
    """Rows of a tab-separated table; ``comment`` and ``@`` lines are skipped, then a header."""
    lines = [
        line
        for line in read_text(path).splitlines()
        if line.strip() and not line.startswith(comment) and not line.startswith("@")
    ]
    if not lines:
        return []
    names = lines[0].split("\t")
    return [
        {name: _convert(value) for name, value in zip(names, line.split("\t"), strict=False)}
        for line in lines[1 : limit + 1]
    ]


def first_integer(path: Path) -> int | None:
    """The first integer found in a small text file (for example a CountReads result)."""
    match = re.search(r"-?\d+", read_text(path, 10_000))
    return int(match.group(0)) if match else None


def delimited_rows(path: Path, *, limit: int = 1000) -> list[dict[str, object]]:
    """Rows of a comma- or tab-separated table (the delimiter is taken from the header line)."""
    lines = [line for line in read_text(path).splitlines() if line.strip()]
    if not lines:
        return []
    delimiter = "\t" if "\t" in lines[0] else ","
    names = lines[0].split(delimiter)
    return [
        {name: _convert(value) for name, value in zip(names, line.split(delimiter), strict=False)}
        for line in lines[1 : limit + 1]
    ]


def count_data_lines(path: Path) -> int | None:
    """Data lines (not ``#`` or ``@``) of an uncompressed text file; ``None`` if compressed."""
    if path.suffix in {".gz", ".bgz", ".bam", ".cram"}:
        return None
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and line[0] not in "#@")


def names(joined: str) -> tuple[str, ...]:
    """The option names of a space-separated string (keeps generated source lines short)."""
    return tuple(joined.split())
