"""Private runner shared by the MMseqs2 command-line facades of integration pack P07.

Each public facade (``mmseqs2_p07_workflow_cli_functions``, ``mmseqs2_p07_taxonomy_cli_functions``
and ``mmseqs2_p07_align_cli_functions``) runs one fixed MMseqs2 module with an explicit argument
list (no shell) through ``dnakit.backends.external.execute_bounded_command``. MMseqs2 databases are
addressed by their path prefix (``prefix``, ``prefix.index``, ``prefix.dbtype``, ...). Every option
is optional and is passed only when the caller gives it, so MMseqs2's own defaults apply
otherwise. Outputs are never overwritten.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics import _gap_cli
from dnakit.exceptions import BackendExecutionError, ConfigurationError, DNAKitError

PROGRAM = "mmseqs"
BACKEND = "mmseqs2"
TESTED_VERSION = "76de5d8ae114750adc59288516d947b6a72a4b05"
CITATION = "https://doi.org/10.1038/nbt.3988"
INSTALL_HINT = (
    "Install MMseqs2 (for example `conda install -c bioconda mmseqs2`, or a release binary from "
    "https://github.com/soedinglab/MMseqs2/releases), then put mmseqs on PATH or pass mmseqs_path."
)
TIMEOUT_SECONDS = 86_400.0
LOG_TAIL_CHARS = 4000
MAX_OUTPUT_BYTES = 5_000_000
MAX_ARTIFACT_BYTES = 2_000_000_000
PREVIEW_ROWS = 1000

_FRAMES = re.compile(r"^[1-3](,[1-3])*$")
_COLUMNS = re.compile(r"^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$")
_TAXON_EXPRESSION = re.compile(r"^[0-9&|!,()]+$")


# ---------------------------------------------------------------------------
# Input and output paths
# ---------------------------------------------------------------------------


def input_file(value: str | Path, name: str) -> Path:
    return _gap_cli.existing_file(value, name)


def input_directory(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def fasta_files(value: str | Path | Sequence[str | Path], name: str) -> list[Path]:
    """One or several existing files (FASTA/FASTQ, plain or compressed; GFF3 for gff2db)."""
    items = [value] if isinstance(value, (str, Path)) else list(value)
    if not items:
        raise ConfigurationError(f"{name} must name at least one file.")
    return [input_file(item, name) for item in items]


def _db_exists(prefix: Path) -> bool:
    return any(Path(f"{prefix}{suffix}").exists() for suffix in (".dbtype", ".index"))


def input_database(value: str | Path, name: str) -> Path:
    """An existing MMseqs2 database, given as its path prefix (``<prefix>.dbtype`` must exist)."""
    prefix = Path(value).expanduser().resolve()
    if not Path(f"{prefix}.dbtype").is_file():
        raise ConfigurationError(
            f"{name} must be the path prefix of an existing MMseqs2 database: {prefix}"
        )
    return prefix


def output_database(value: str | Path, name: str) -> Path:
    """The path prefix of a database that does not exist yet; its folder must exist."""
    prefix = Path(value).expanduser().resolve()
    if prefix.exists() or _db_exists(prefix):
        raise ConfigurationError(f"{name} must not exist yet: {prefix}")
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {prefix.parent}")
    return prefix


def output_file(value: str | Path, name: str) -> Path:
    return _gap_cli.new_file(value, name)


def choose_kind(value: object) -> str:
    """``fasta`` (easy-* workflows, FASTA in and text out) or ``database`` (database in and out)."""
    return _gap_cli.choice(value, ("fasta", "database"), "input_kind")


def fasta_or_database(value: str | Path, name: str) -> Path:
    """A FASTA file or the prefix of an existing database (easy-* workflows accept both)."""
    path = Path(value).expanduser().resolve()
    if path.is_file():
        return path
    if Path(f"{path}.dbtype").is_file():
        return path
    raise ConfigurationError(f"{name} must be an existing FASTA file or database prefix: {path}")


def output_family(
    value: str | Path, suffixes: Mapping[str, str], name: str
) -> tuple[Path, dict[str, Path]]:
    """A new file prefix and the ``<prefix><suffix>`` files an easy-* workflow writes."""
    prefix = Path(value).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {prefix.parent}")
    family = {key: Path(f"{prefix}{suffix}") for key, suffix in suffixes.items()}
    clashes = [str(path) for path in family.values() if path.exists()]
    if clashes:
        raise ConfigurationError(f"{name} would overwrite existing files: {clashes}")
    return prefix, family


# ---------------------------------------------------------------------------
# Option builders: nothing is passed unless the caller gives a value
# ---------------------------------------------------------------------------


def opt_int(
    flag: str, value: object, name: str, minimum: int | None = None, maximum: int | None = None
) -> list[str]:
    if value is None:
        return []
    return [flag, str(_gap_cli.integer(value, name, minimum=minimum, maximum=maximum))]


def opt_float(
    flag: str, value: object, name: str, minimum: float | None = None, maximum: float | None = None
) -> list[str]:
    if value is None:
        return []
    return [flag, repr(_gap_cli.number(value, name, minimum=minimum, maximum=maximum))]


def opt_bool(flag: str, value: object, name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return [flag, "1" if value else "0"]


def opt_choice(flag: str, value: object, options: Sequence[int], name: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, bool) or not isinstance(value, int) or value not in options:
        raise ConfigurationError(
            f"{name} must be one of {', '.join(str(option) for option in options)}."
        )
    return [flag, str(value)]


def opt_text(flag: str, value: object, name: str, pattern: str | None = None) -> list[str]:
    if value is None:
        return []
    cleaned = _gap_cli.text(value, name)
    if pattern is not None and re.fullmatch(pattern, cleaned) is None:
        raise ConfigurationError(f"{name} has an unsupported form: {cleaned!r}.")
    return [flag, cleaned]


def opt_file(flag: str, value: str | Path | None, name: str) -> list[str]:
    return [] if value is None else [flag, str(input_file(value, name))]


def opt_directory(flag: str, value: str | Path | None, name: str) -> list[str]:
    return [] if value is None else [flag, str(input_directory(value, name))]


def frames(flag: str, value: str | None, name: str) -> list[str]:
    return opt_text(flag, value, name, _FRAMES.pattern)


def columns(flag: str, value: str | None, name: str) -> list[str]:
    return opt_text(flag, value, name, _COLUMNS.pattern)


def taxon_list(flag: str, value: str | None, name: str) -> list[str]:
    return opt_text(flag, value, name, _TAXON_EXPRESSION.pattern)


def reject_with(condition: bool, message: str) -> None:
    if condition:
        raise ConfigurationError(message)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


@lru_cache(maxsize=16)
def version_of(executable: str) -> str | None:
    """MMseqs2 prints its release or commit identifier for ``mmseqs version``."""
    try:
        result = execute_bounded_command(
            executable,
            ["version"],
            backend_id=BACKEND,
            cwd=Path(executable).parent,
            timeout_seconds=60.0,
            max_output_bytes=100_000,
            verbose=False,
        )
    except DNAKitError:
        return None
    for line in (result.output or "").splitlines():
        line = line.strip()
        if line:
            return line
    return None


def _common_folder(paths: Sequence[Path]) -> Path:
    folders = [path.parent for path in paths]
    first = folders[0]
    for folder in folders[1:]:
        while folder != first and first not in folder.parents:
            first = first.parent
    return first


def _entry_count(prefix: Path) -> int | None:
    index = Path(f"{prefix}.index")
    if not index.is_file():
        return None
    with index.open("rb") as handle:
        return sum(1 for _ in handle)


def _file_preview(path: Path) -> tuple[int, list[list[str]]]:
    rows: list[list[str]] = []
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            count += 1
            if len(rows) < PREVIEW_ROWS:
                rows.append(line.rstrip("\n").split("\t"))
    return count, rows


def run(
    operation: str,
    subcommand: str,
    positionals: Sequence[str | Path],
    flags: Sequence[str],
    *,
    databases: Mapping[str, Path] | None = None,
    files: Mapping[str, Path] | None = None,
    extra_outputs: Mapping[str, Path] | None = None,
    needs_tmp: bool = False,
    threads: int | None = 1,
    mmseqs_path: str | None = None,
    anchor: Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run ``mmseqs <subcommand> <positionals> <flags>`` and collect the new databases and files.

    ``databases`` and ``files`` name the outputs this call must create (database prefixes and plain
    files); ``extra_outputs`` are additional plain files that must exist afterwards. A
    temporary folder
    next to the outputs is passed as the final argument when ``needs_tmp`` and removed afterwards.
    """
    databases = dict(databases or {})
    files = dict(files or {})
    extra_outputs = dict(extra_outputs or {})
    expected_paths = [*databases.values(), *files.values(), *extra_outputs.values()]
    if not expected_paths and anchor is None:
        raise ConfigurationError("An MMseqs2 facade needs at least one output.")
    cwd = _common_folder(expected_paths) if expected_paths else anchor
    assert cwd is not None
    if not cwd.is_dir():
        raise ConfigurationError(f"The output folder does not exist: {cwd}")
    thread_flags = (
        []
        if threads is None
        else ["--threads", str(_gap_cli.integer(threads, "threads", minimum=1, maximum=1024))]
    )
    executable = _gap_cli.program(PROGRAM, explicit=mmseqs_path, hint=INSTALL_HINT)
    tmp_folder = Path(tempfile.mkdtemp(prefix=".mmseqs_tmp_", dir=cwd)) if needs_tmp else None
    try:
        arguments = [subcommand, *(str(item) for item in positionals)]
        if tmp_folder is not None:
            arguments.append(str(tmp_folder))
        arguments += [*flags, *thread_flags]
        result = execute_bounded_command(
            executable,
            arguments,
            backend_id=BACKEND,
            cwd=cwd,
            timeout_seconds=TIMEOUT_SECONDS,
            max_output_bytes=MAX_OUTPUT_BYTES,
            monitored_output_paths=expected_paths,
            max_monitored_output_bytes=MAX_ARTIFACT_BYTES,
            verbose=verbose,
        )
    finally:
        if tmp_folder is not None:
            shutil.rmtree(tmp_folder, ignore_errors=True)
    log_tail = (result.output or "")[-LOG_TAIL_CHARS:]
    if result.return_code != 0:
        raise BackendExecutionError(
            f"MMseqs2 {subcommand} exited with code {result.return_code}.",
            code="MMSEQS2_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": log_tail},
        )
    missing = [
        str(path) for path in [*files.values(), *extra_outputs.values()] if not path.is_file()
    ]
    missing += [str(prefix) for prefix in databases.values() if not _db_exists(prefix)]
    if missing:
        raise BackendExecutionError(
            f"MMseqs2 {subcommand} did not write the expected outputs: {missing}",
            code="MMSEQS2_OUTPUT_MISSING",
            context={"operation": operation, "log_tail": log_tail},
        )
    version = version_of(executable)
    outputs: dict[str, str] = {}
    entry_counts: dict[str, int | None] = {}
    for key, prefix in databases.items():
        outputs[key] = str(prefix)
        entry_counts[key] = _entry_count(prefix)
    line_counts: dict[str, int] = {}
    for key, path in {**files, **extra_outputs}.items():
        outputs[key] = str(path)
        line_counts[key] = (
            _file_preview(path)[0] if path.stat().st_size <= MAX_ARTIFACT_BYTES else -1
        )
    summary: dict[str, Any] = {
        "status": "completed",
        "backend": BACKEND,
        "operation": operation,
        "version": version,
        "command": [PROGRAM, *arguments],
        "outputs": outputs,
        "log_tail": log_tail,
        "provenance": {
            "tool": "MMseqs2",
            "version": version,
            "tested_version": TESTED_VERSION,
            "executable": executable,
            "citation": CITATION,
        },
    }
    if entry_counts:
        summary["entry_counts"] = entry_counts
    if line_counts:
        summary["line_counts"] = line_counts
    return summary


def table(path: Path, header: Sequence[str] | None = None) -> list[Any]:
    """Up to ``PREVIEW_ROWS`` rows of a tab-separated file; dicts when column names are given."""
    _, rows = _file_preview(path)
    if header is None:
        return rows
    return [dict(zip(header, row, strict=False)) for row in rows]
