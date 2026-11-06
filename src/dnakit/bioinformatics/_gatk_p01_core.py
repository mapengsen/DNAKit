"""Private helpers shared by the GATK command-line facades of pack P01.

One GATK tool is run with an explicit argument list (no shell) through
``dnakit.backends.external.execute_bounded_command``. The launcher is resolved from an explicit path, the
``DNAKIT_GATK`` environment variable or PATH. The arguments a caller may add are limited to the long options
that ``gatk <Tool> --help`` lists for that tool (table in ``_gatk_p01_options``); the options a function sets
itself cannot be repeated. A missing launcher raises ``BackendUnavailableError`` with an install hint, a failed
run raises ``BackendExecutionError`` with the tail of the log.
"""

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping, Sequence
from contextlib import suppress
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from ._gatk_p01_options import TOOL_OPTIONS

PROGRAM = "gatk"
TESTED_VERSION = "4.7.0.0"
CITATION = "https://doi.org/10.1101/gr.107524.110"
INSTALL_HINT = (
    "Install GATK 4.7.0.0 or later (for example `conda install -c bioconda gatk4`, or the release zip from "
    "https://github.com/broadinstitute/gatk/releases) together with Java 17, then put the `gatk` launcher on "
    "PATH, set DNAKIT_GATK or pass gatk_path."
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


def flag_value(value: object, name: str) -> bool:
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


def new_file(value: object, name: str) -> Path:
    """A path that does not exist yet and whose folder exists; nothing is overwritten."""
    path = Path(text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def new_directory(value: object, name: str) -> Path:
    """A directory that does not exist yet; it is created together with missing parents' check."""
    path = Path(text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be reused: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The parent folder of {name} does not exist: {path.parent}")
    return path


def existing_directory(value: object, name: str) -> Path:
    path = Path(text(value, name)).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


# ---------------------------------------------------------------------------
# Options
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


def option_arguments(
    tool: str, options: Mapping[str, Any] | None, *, reserved: Sequence[str] = ()
) -> list[str]:
    """Turn ``{"max-records-in-ram": 5000}`` into ``["--max-records-in-ram", "5000"]`` for the listed options only.

    A list value repeats the option; ``True`` / ``False`` become ``true`` / ``false``. Options that the calling
    function already sets (``reserved``) are refused.
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
                f"option {name!r} is set by a parameter of this function; use that parameter instead."
            )
        if name not in allowed:
            raise ConfigurationError(f"gatk {tool} option {name!r} is not allowed by this facade.")
        values = value if isinstance(value, (list, tuple)) else [value]
        if not values:
            raise ConfigurationError(f"option {name!r} has an empty list.")
        for item in values:
            arguments += [f"--{name}", _scalar(item, name)]
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
    max_output_bytes: int = MAX_OUTPUT_BYTES,
) -> dict[str, Any]:
    """Run ``gatk [--java-options -Xmx..] <tool> <arguments>``; every ``expected`` file must exist afterwards.

    ``expected`` paths must be inside ``cwd`` (the executor monitors them there).
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
        max_output_bytes=max_output_bytes,
        monitored_output_paths=tuple(expected),
        max_monitored_output_bytes=MAX_ARTIFACT_BYTES,
        verbose=verbose,
    )
    if result.return_code != 0:
        problems = [
            line.strip()
            for line in result.output.splitlines()
            if re.search(r"Exception|USER ERROR|Error:|^Argument '", line.strip())
            and not line.lstrip().startswith(("at ", "-"))
        ]
        raise BackendExecutionError(
            f"gatk {tool} exited with code {result.return_code}."
            + (f" {problems[0][:300]}" if problems else ""),
            code="GATK_EXECUTION_FAILED",
            context={
                "operation": operation,
                "log_tail": tail(result.output),
                "errors": problems[:5],
            },
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
    """Drop the private log, add the output paths (``None`` outputs are omitted) and extra fields."""
    result.pop("_log", None)
    result["outputs"] = {name: str(path) for name, path in outputs.items() if path is not None}
    result.update(extra)
    return result


# ---------------------------------------------------------------------------
# Declarative operation specifications
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Param:
    """One named parameter of a public function and the GATK option it becomes.

    ``kind`` is ``in`` (existing file), ``in_dir`` (existing directory), ``ins`` (list of existing files, option
    repeated), ``src`` (existing file or ``gendb://<existing directory>``), ``out`` (new file), ``outd`` (new
    directory created here), ``outd_tool`` (new directory created by the tool), ``str``, ``strs`` (repeated
    option), ``int``, ``float``, ``bool`` or ``choice``. ``flag`` may differ per variant of a merged function.
    ``v`` lists the variants the parameter applies to; ``req`` is True, False or the variants that require it.
    """

    name: str
    kind: str
    flag: str | Mapping[str, str] | None
    v: tuple[str, ...] | None = None
    req: bool | tuple[str, ...] = False
    choices: tuple[str, ...] | None = None
    minimum: float | None = None


@dataclass(frozen=True)
class Spec:
    operation: str
    tools: Mapping[str, str]
    params: tuple[Param, ...]
    selector: str | None = None
    any_of: tuple[tuple[tuple[str, ...] | None, tuple[str, ...]], ...] = ()
    # (output key, input parameter, suffixes, all_created): files the tool writes next to its input
    derived: tuple[tuple[str, str, tuple[str, ...], bool], ...] = ()
    tables: Mapping[str, str] = field(default_factory=dict)
    hook: str | None = None
    big_output: bool = False


def _flag_of(param: Param, variant: str) -> str | None:
    flag = param.flag
    if flag is None or isinstance(flag, str):
        return flag
    return flag.get(variant)


def _used(param: Param, value: object, variant: str) -> bool:
    """True when the value is to be passed; raises when a value is given for a variant it does not apply to."""
    if param.v is not None and variant not in param.v:
        if value is not None:
            raise ConfigurationError(f"{param.name} is not used when the mode is {variant!r}.")
        return False
    needed = param.req is True or (isinstance(param.req, tuple) and variant in param.req)
    if value is None:
        if needed:
            raise ConfigurationError(f"{param.name} is required when the mode is {variant!r}.")
        return False
    return True


def _source(value: object, name: str) -> tuple[str, Path]:
    if isinstance(value, str) and value.startswith("gendb://"):
        directory = existing_directory(value[len("gendb://") :], name)
        return f"gendb://{directory}", directory
    path = existing_file(value, name)
    return str(path), path


def table_rows(path: Path, kind: str, *, limit: int = 500) -> list[dict[str, str]]:
    """Rows of a tab-separated results file: ``tsv`` (``#`` lines skipped) or Picard ``metrics``."""
    lines = read_text(path).splitlines()
    if kind == "metrics":
        for index, line in enumerate(lines):
            if line.startswith("## METRICS CLASS") and index + 2 <= len(lines):
                body: list[str] = []
                for row in lines[index + 2 :]:
                    if not row.strip():
                        break
                    body.append(row)
                names = lines[index + 1].split("\t")
                return [dict(zip(names, row.split("\t"), strict=False)) for row in body[:limit]]
        return []
    data = [line for line in lines if line.strip() and not line.startswith("#")]
    if not data:
        return []
    names = data[0].split("\t")
    return [dict(zip(names, row.split("\t"), strict=False)) for row in data[1:][:limit]]


_BLOCK = re.compile(
    r"Block #(\d+) at file offset (\d+)\s+- compressed size: (\d+)\s+- uncompressed size: (\d+)"
)
_BASE_LINE = re.compile(r"^([A-Za-z]) : (\d+)\s*$", re.MULTILINE)


def _hook_results(hook: str | None, log: str) -> dict[str, object]:
    if hook == "bgzf_blocks":
        blocks = [
            {
                "block": int(m.group(1)),
                "offset": int(m.group(2)),
                "compressed_size": int(m.group(3)),
                "uncompressed_size": int(m.group(4)),
            }
            for m in _BLOCK.finditer(log)
        ]
        return {"block_count": len(blocks), "blocks": blocks[:10000]}
    if hook == "base_counts":
        counts = {m.group(1): int(m.group(2)) for m in _BASE_LINE.finditer(log)}
        return {"base_counts": counts, "total_bases": sum(counts.values())}
    if hook == "tab_table":
        lines = log.splitlines()
        for index, line in enumerate(lines):
            if line.startswith("MD5\t"):
                rows = []
                for row in lines[index:]:
                    if "\t" not in row:
                        break
                    rows.append(row.split("\t"))
                return {"comparison_table": rows}
        return {"comparison_table": []}
    return {}


def execute(
    spec: Spec,
    values: Mapping[str, Any],
    *,
    mode: str | None,
    options: Mapping[str, Any] | None,
    gatk_path: str | None,
    java_max_memory: str | None,
    verbose: bool,
) -> dict[str, Any]:
    """Validate the arguments of one public function, run its GATK tool and collect the results."""
    variant = "default" if spec.selector is None else choice(mode, tuple(spec.tools), spec.selector)
    tool = spec.tools[variant]
    arguments: list[str] = []
    reserved: list[str] = []
    outputs: dict[str, Path] = {}
    expected: list[Path] = []
    first_input: Path | None = None
    created: list[Path] = []
    for param in spec.params:
        value = values.get(param.name)
        if not _used(param, value, variant):
            continue
        flag = _flag_of(param, variant)
        kind = param.kind
        pieces: list[str]
        if kind == "in":
            path = existing_file(value, param.name)
            first_input = first_input or path
            pieces = [str(path)]
        elif kind == "in_dir":
            path = existing_directory(value, param.name)
            first_input = first_input or path
            pieces = [str(path)]
        elif kind == "ins":
            files = existing_files(value, param.name)
            first_input = first_input or files[0]
            pieces = [str(item) for item in files]
        elif kind == "src":
            text_value, path = _source(value, param.name)
            first_input = first_input or path
            pieces = [text_value]
        elif kind == "out":
            path = new_file(value, param.name)
            outputs[param.name] = path
            expected.append(path)
            pieces = [str(path)]
        elif kind in ("outd", "outd_tool"):
            path = new_directory(value, param.name)
            outputs[param.name] = path
            expected.append(path)
            if kind == "outd":
                created.append(path)
            pieces = [str(path)]
        elif kind == "str":
            pieces = [text(value, param.name)]
        elif kind == "strs":
            if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
                raise ConfigurationError(f"{param.name} must be a non-empty list of strings.")
            pieces = [text(item, param.name) for item in value]
        elif kind == "int":
            pieces = [
                str(
                    integer(
                        value,
                        param.name,
                        minimum=None if param.minimum is None else int(param.minimum),
                    )
                )
            ]
        elif kind == "float":
            pieces = [str(number(value, param.name, minimum=param.minimum))]
        elif kind == "bool":
            pieces = ["true" if flag_value(value, param.name) else "false"]
        elif kind == "choice":
            pieces = [choice(value, param.choices or (), param.name)]
        else:  # pragma: no cover - specification error
            raise ConfigurationError(f"unknown parameter kind {kind!r}.")
        if flag is not None:
            reserved.append(flag)
            for piece in pieces:
                arguments += [f"--{flag}", piece]
    for scope, names in spec.any_of:
        if (scope is None or variant in scope) and all(values.get(name) is None for name in names):
            raise ConfigurationError(f"at least one of {', '.join(names)} is required.")
    derived_paths: dict[str, list[Path]] = {}
    for key, source_name, suffixes, all_created in spec.derived:
        source_value = values.get(source_name)
        if source_value is None or f"{key}_path" in outputs:
            continue
        base = Path(str(source_value).removeprefix("gendb://")).expanduser().resolve()
        candidates = [base.with_name(base.name + suffix) for suffix in suffixes]
        for candidate in candidates:
            if candidate.exists():
                raise ConfigurationError(f"{candidate} already exists and will not be overwritten.")
        derived_paths[key] = candidates
        if all_created:
            expected.extend(candidates)
    arguments += option_arguments(tool, options, reserved=reserved)
    if expected:
        parents = [str(path.parent) for path in expected]
        cwd = Path(os.path.commonpath(parents))
    elif first_input is not None:
        cwd = first_input.parent if first_input.is_file() else first_input
    else:
        raise ConfigurationError("this operation needs at least one input or output path.")
    for directory in created:
        directory.mkdir()
    try:
        result = run_tool(
            spec.operation,
            tool,
            arguments,
            cwd=cwd,
            expected=expected,
            gatk_path=gatk_path,
            java_max_memory=java_max_memory,
            verbose=verbose,
            max_output_bytes=100_000_000 if spec.big_output else MAX_OUTPUT_BYTES,
        )
    except BaseException:
        for directory in created:
            with suppress(OSError):
                directory.rmdir()
        raise
    extra: dict[str, Any] = {}
    if spec.selector is not None:
        extra["mode"] = variant
    produced: dict[str, Path | None] = dict(outputs)
    for key, candidates in derived_paths.items():
        found = [path for path in candidates if path.exists()]
        if len(found) == 1:
            produced[key] = found[0]
        elif found:
            produced[key] = found[0]
            extra[f"{key}_files"] = [str(path) for path in found]
    listings = {
        name: sorted(item.name for item in path.iterdir())
        for name, path in outputs.items()
        if path.is_dir()
    }
    if listings:
        extra["files"] = listings
    if spec.tables:
        extra["tables"] = {
            name: table_rows(outputs[name], kind)
            for name, kind in spec.tables.items()
            if name in outputs and outputs[name].is_file()
        }
    extra.update(_hook_results(spec.hook, result["_log"]))
    return finish(result, outputs=produced, **extra)


def read_text(path: Path, limit: int = 1_000_000) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return handle.read(limit)


def tsv_rows(path: Path, *, comment: str = "#") -> list[dict[str, str]]:
    """Rows of a GATK/Picard metrics-style table: ``#`` lines skipped, first remaining line is the header."""
    lines = [
        line
        for line in read_text(path).splitlines()
        if line.strip() and not line.startswith(comment)
    ]
    if not lines:
        return []
    names = lines[0].split("\t")
    return [dict(zip(names, line.split("\t"), strict=False)) for line in lines[1:]]


def count_vcf_records(path: Path) -> int | None:
    """Number of data lines of an uncompressed VCF (``None`` for compressed files)."""
    if path.suffix == ".gz":
        return None
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))
