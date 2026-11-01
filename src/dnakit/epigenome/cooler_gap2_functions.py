"""cooler command-line operations as typed DNAKit facades (balance, cload, coarsen, csort, digest,
dump, fileops, info, load, makebins, merge, show, zoomify) plus headless counterparts of the
interactive display.

The user installs cooler (``python -m pip install cooler`` or ``conda install -c bioconda cooler``).
Each function runs one ``cooler`` subcommand with an explicit, typed keyword list that maps to fixed
flags; no other option can be passed and no shell is used. Programs are looked up as an explicit
``cooler_path``, then on ``PATH``, then next to the running Python interpreter (where a virtual
environment keeps its console scripts). Outputs are written to explicit new paths; modes that modify
a file in place are named in the docstring. Results are JSON-serialisable dictionaries with status,
backend, version, command, outputs and a bounded log tail.

The four ``cooler_csort_*`` step functions expose the pieces of ``cooler csort`` one by one (read,
flip, sort, index). The flip step runs in-process, the others call the system programs (gzip or
pigz, sort, tabix or pairix) with argument lists. ``cooler_show_view`` and ``cooler_pairs_header``
use the cooler and matplotlib Python packages in-process.

Verified against cooler 0.10.4. Known upstream defects are noted in the docstrings concerned.
"""

from __future__ import annotations

import functools
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import require
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)

_COOLER_HINT = (
    "Install cooler 0.10 or later "
    "(python -m pip install cooler, or conda install -c bioconda cooler)."
)
_SORT_HINT = (
    "sort is part of GNU coreutils (apt install coreutils); it is present on Linux and macOS."
)
_HTSLIB_HINT = "Install htslib, which provides bgzip and tabix (conda install -c bioconda htslib)."
_PAIRIX_HINT = (
    "Install pairix (conda install -c bioconda pairix) "
    "and pypairix (python -m pip install pypairix)."
)
_GZIP_HINT = "Install pigz (conda install -c conda-forge pigz) or gzip (apt install gzip)."
_PYTHON_HINT = 'python -m pip install "dnakit[bio-cooler,bio-plots,bio-analysis]"'
_MAX_TAIL = 4000
_MAX_CAPTURE = 20_000_000
_DATA_CELL_LIMIT = 4_000_000
_BIG_TAIL = 500_000
_DEFAULT_TIMEOUT = 3600.0
_MAX_MATRIX_CELLS = 40_000
_SCALES = ("linear", "log2", "log10")
_IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".svg", ".pdf", ".ps", ".eps")

_FIELD_WITH_COLUMN = re.compile(r"^[A-Za-z_]\w*(=[1-9]\d*)?(:(dtype|agg)=\w+(,(dtype|agg)=\w+)*)?$")
_FIELD_NAME_ONLY = re.compile(r"^[A-Za-z_]\w*(:(dtype|agg)=\w+(,(dtype|agg)=\w+)*)?$")
_STORAGE_OPTIONS = re.compile(r"^\w+=[\w.]+(,\w+=[\w.]+)*$")
_ASSEMBLY = re.compile(r"^[A-Za-z0-9_][\w.\-]*$")
_ENZYME = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
_FLOAT_FORMAT = re.compile(r"^[-+ #0]*\d*(\.\d+)?[eEfFgG]$")
_RESOLUTION_TOKEN = re.compile(r"^(\d+|\d*[nNbB]|4[dD][nN])$")
_COLUMN_NAME = re.compile(r"^[A-Za-z_][\w.\-]*$")
_SORT_BUFFER = re.compile(r"^\d+[%bkKMGTPEZYRQ]?$")
_COLORMAP = re.compile(r"^[A-Za-z0-9_]+$")
_VERSION = re.compile(r"version\s+(\S+)")


# ------------------------------------------------------------------ program lookup and running


def _find_program(name: str, explicit: str | None, hint: str) -> str:
    """Explicit path, then PATH, then the folder of the running Python interpreter."""

    if explicit is not None:
        candidate = Path(explicit).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name}_path does not exist: {candidate}", hint=hint)
        return str(candidate)
    found = shutil.which(name)
    if found is not None:
        return found
    sibling = Path(sys.executable).parent / name
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)
    raise BackendUnavailableError(
        f"{name} was not found on PATH or next to the Python interpreter ({sibling.parent}).",
        hint=hint,
    )


def _path_with(program: str) -> str:
    """PATH with the program's own folder first, so sibling tools (bgzip, tabix) are found."""

    return f"{Path(program).absolute().parent}{os.pathsep}{os.environ.get('PATH', '')}"


def _tool_on_path(name: str, program: str, hint: str) -> str:
    found = shutil.which(name, path=_path_with(program))
    if found is None:
        raise BackendUnavailableError(
            f"cooler csort needs {name}, which was not found on PATH.", hint=hint
        )
    return found


@functools.lru_cache(maxsize=16)
def _cooler_version(program: str) -> str | None:
    try:
        result = execute_bounded_command(
            program,
            ["--version"],
            backend_id="cooler:version",
            cwd=Path(program).absolute().parent,
            timeout_seconds=120,
            verbose=False,
        )
    except (
        BackendExecutionError,
        BackendTimeoutError,
        BackendUnavailableError,
        ConfigurationError,
    ):
        return None
    match = _VERSION.search(result.output or "")
    return match.group(1) if match else None


def _timeout(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value <= 86_400:
        raise ConfigurationError("timeout_seconds must be a number in (0, 86400].")
    return float(value)


def _run_cooler(
    subcommand: str,
    arguments: list[str],
    *,
    cwd: Path,
    outputs: Sequence[Path],
    cooler_path: str | None,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
    tail_chars: int = _MAX_TAIL,
) -> dict[str, Any]:
    timeout = _timeout(timeout_seconds)
    program = _find_program("cooler", cooler_path, _COOLER_HINT)
    env = {"PATH": _path_with(program)}
    env.update(environment or {})
    result = execute_bounded_command(
        program,
        [subcommand, *arguments],
        backend_id=f"cooler:{subcommand}",
        cwd=cwd,
        timeout_seconds=timeout,
        max_output_bytes=_MAX_CAPTURE,
        environment_overrides=env,
        verbose=False,
    )
    ok = result.return_code == 0 and all(path.exists() for path in outputs)
    return {
        "status": "completed" if ok else "failed",
        "backend": "cooler",
        "backend_version": _cooler_version(program),
        "subcommand": subcommand,
        "command": ["cooler", subcommand, *arguments],
        "return_code": result.return_code,
        "outputs": [str(path) for path in outputs if path.exists()],
        "log_tail": (result.output or "")[-tail_chars:],
    }


def _run_to_file(
    program: str, arguments: list[str], output: Path, timeout_seconds: float
) -> dict[str, Any]:
    """Run ``program arguments`` with its standard output written to ``output`` (argument list, no
    shell).
    """

    timeout = _timeout(timeout_seconds)
    try:
        with output.open("wb") as handle:
            completed = subprocess.run(
                [program, *arguments],
                stdin=subprocess.DEVNULL,
                stdout=handle,
                stderr=subprocess.PIPE,
                shell=False,
                timeout=timeout,
                check=False,
            )
    except subprocess.TimeoutExpired as exc:
        output.unlink(missing_ok=True)
        raise BackendTimeoutError(
            "External backend command exceeded its timeout.", context={"program": program}
        ) from exc
    except OSError as exc:
        output.unlink(missing_ok=True)
        raise BackendExecutionError(
            "External backend command could not be started.", context={"program": program}
        ) from exc
    if completed.returncode != 0:
        output.unlink(missing_ok=True)
    return {
        "return_code": completed.returncode,
        "log_tail": completed.stderr.decode("utf-8", errors="replace")[-_MAX_TAIL:],
    }


# --------------------------------------------------------------------------- argument checks


def _integer(value: object, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError(f"{name} must be an integer >= {minimum}.")
    return value


def _number(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if positive and value <= 0:
        raise ConfigurationError(f"{name} must be positive.")
    return float(value)


def _text(value: object, name: str, *, option_value: bool = False) -> str:
    """A non-empty string; positional values must not start with '-' (read as options)."""

    if not isinstance(value, str) or not value.strip() or "\x00" in value or len(value) > 4096:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if value.startswith("-") and not option_value:
        raise ConfigurationError(f"{name} must not start with '-'.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}.")
    return value


def _existing(path: object, name: str) -> Path:
    if not isinstance(path, (str, os.PathLike)):
        raise ConfigurationError(f"{name} must be a path.")
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _existing_folder(path: object, name: str) -> Path:
    if not isinstance(path, (str, os.PathLike)):
        raise ConfigurationError(f"{name} must be a path.")
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise ConfigurationError(f"{name} must be an existing folder: {resolved}")
    return resolved


def _new(path: object, name: str) -> Path:
    if not isinstance(path, (str, os.PathLike)):
        raise ConfigurationError(f"{name} must be a path.")
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _split_uri(uri: object, name: str) -> tuple[str, str | None]:
    if not isinstance(uri, (str, os.PathLike)):
        raise ConfigurationError(f"{name} must be a path or a cooler URI.")
    parts = os.fspath(uri).split("::")
    if len(parts) > 2 or not parts[0]:
        raise ConfigurationError(f"{name} is not a valid cooler URI: {uri!r}")
    return parts[0], (parts[1] if len(parts) == 2 else None)


def _cool_input(uri: object, name: str) -> tuple[Path, str]:
    """Existing cooler file with an optional ``::group``; returns the file and the absolute URI."""

    file_text, group = _split_uri(uri, name)
    path = _existing(file_text, name)
    return path, (f"{path}::{group}" if group else str(path))


def _cool_output(uri: object, name: str, append: bool) -> tuple[Path, str]:
    """New cooler file (or an existing one when appending) with an optional ``::group``."""

    file_text, group = _split_uri(uri, name)
    path = Path(file_text).expanduser().resolve()
    if path.exists() and not append:
        raise ConfigurationError(
            f"{name} must not exist yet: {path} "
            "(pass append=True to add a cooler to an existing file)."
        )
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path, (f"{path}::{group}" if group else str(path))


def _fields(values: Sequence[str] | None, name: str, *, column_number: bool) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        raise ConfigurationError(f"{name} must be a list of strings.")
    pattern = _FIELD_WITH_COLUMN if column_number else _FIELD_NAME_ONLY
    shape = (
        "<name>=<column number>[:dtype=<dtype>][,agg=<agg>]"
        if column_number
        else "<name>[:dtype=<dtype>][,agg=<agg>]"
    )
    result: list[str] = []
    for value in values:
        if not isinstance(value, str) or not pattern.match(value):
            raise ConfigurationError(f"Each entry of {name} must look like {shape}; got {value!r}.")
        result.append(value)
    return result


def _bins_spec(value: object, name: str) -> str:
    """``<chromsizes file>:<bin size in bp>`` or the path of a BED file with the bins."""

    text = _text(value, name)
    if ":" in os.path.splitdrive(text)[1]:
        sizes_text, _, size_text = text.rpartition(":")
        sizes = _existing(sizes_text, f"the chromsizes file of {name}")
        try:
            size = int(size_text)
        except ValueError as exc:
            raise ConfigurationError(
                f"{name} must end with an integer bin size in bp, got {size_text!r}."
            ) from exc
        if size < 1:
            raise ConfigurationError(f"The bin size in {name} must be positive.")
        return f"{sizes}:{size}"
    return str(_existing(text, name))


def _region(value: object, name: str) -> str:
    text = _text(value, name)
    if any(character.isspace() for character in text):
        raise ConfigurationError(f"{name} must not contain whitespace: {text!r}")
    return text


def _single_character(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 8
        or any(character in value for character in "\"'\\\r\n\x00")
    ):
        raise ConfigurationError(
            f"{name} must be a short string without quotes, backslashes or newlines."
        )
    return value


def _optional_metadata(metadata_path: str | None, assembly: str | None) -> list[str]:
    arguments: list[str] = []
    if metadata_path is not None:
        arguments += ["--metadata", str(_existing(metadata_path, "metadata_path"))]
    if assembly is not None:
        if not isinstance(assembly, str) or not _ASSEMBLY.match(assembly):
            raise ConfigurationError("assembly must be a name such as hg38 or mm10.")
        arguments += ["--assembly", assembly]
    return arguments


def _storage_options(value: str | None) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, str) or not _STORAGE_OPTIONS.match(value):
        raise ConfigurationError(
            "storage_options must look like 'compression=gzip,compression_opts=6'."
        )
    return ["--storage-options", value]


def _comma_names(values: Sequence[str] | None, name: str) -> str | None:
    if values is None:
        return None
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of column names.")
    for value in values:
        if not isinstance(value, str) or not _COLUMN_NAME.match(value):
            raise ConfigurationError(f"{name} must contain plain column names; got {value!r}.")
    return ",".join(values)


def _text_preview(
    path: Path, lines: int = 5, count_limit: int = 256_000_000
) -> tuple[int | None, list[str]]:
    """Line count (None above ``count_limit`` bytes) and first lines of a text or .gz file."""

    opener: Any = gzip.open if path.suffix == ".gz" else open
    count_all = path.stat().st_size <= count_limit
    preview: list[str] = []
    count = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            count += 1
            if len(preview) < lines:
                preview.append(line.rstrip("\n"))
            elif not count_all:
                break
    return (count if count_all else None), preview


def _plain(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


# --------------------------------------------------------------------------- balance


def cooler_balance(
    cool_uri: str,
    *,
    cis_only: bool = False,
    trans_only: bool = False,
    ignore_diags: int | None = None,
    ignore_dist: int | None = None,
    mad_max: int | None = None,
    min_nnz: int | None = None,
    min_count: int | None = None,
    blacklist: str | None = None,
    nproc: int | None = None,
    chunksize: int | None = None,
    tol: float | None = None,
    max_iters: int | None = None,
    name: str | None = None,
    force: bool = False,
    convergence_policy: Literal["store_final", "store_nan", "discard", "error"] | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Balance a cooler by iterative correction and store the weights (cooler balance).

    The cooler file is modified in place: a ``weight`` column (or ``name``) is added to ``bins``.
    Options left as None keep cooler's own defaults (ignore_diags 2, mad_max 5, min_nnz 10,
    min_count 0, tol 1e-5, max_iters 200, convergence_policy store_final, 8 processes). A
    ``blacklist`` BED file needs at least two lines: cooler 0.10.4 reads a one-line file as a header
    and fails with "need at least one array to concatenate". ``cooler balance --stdout`` is not
    exposed: read the weights afterwards with ``cooler_dump(table="bins", columns=["weight"])``. The
    result adds ``weights_stored``, the outcome of ``cooler balance --check`` after the run.
    """

    path, uri = _cool_input(cool_uri, "cool_uri")
    if cis_only and trans_only:
        raise ConfigurationError("Provide at most one of cis_only and trans_only.")
    arguments: list[str] = []
    if cis_only:
        arguments.append("--cis-only")
    if trans_only:
        arguments.append("--trans-only")
    for flag, value, label, minimum in (
        ("--ignore-diags", ignore_diags, "ignore_diags", 0),
        ("--ignore-dist", ignore_dist, "ignore_dist", 0),
        ("--mad-max", mad_max, "mad_max", 0),
        ("--min-nnz", min_nnz, "min_nnz", 0),
        ("--min-count", min_count, "min_count", 0),
        ("--nproc", nproc, "nproc", 1),
        ("--chunksize", chunksize, "chunksize", 1),
        ("--max-iters", max_iters, "max_iters", 1),
    ):
        if value is not None:
            arguments += [flag, str(_integer(value, label, minimum))]
    if tol is not None:
        arguments += ["--tol", repr(_number(tol, "tol", positive=True))]
    if blacklist is not None:
        arguments += ["--blacklist", str(_existing(blacklist, "blacklist"))]
    column = "weight" if name is None else _text(name, "name")
    if name is not None:
        arguments += ["--name", column]
    if force:
        arguments.append("--force")
    if convergence_policy is not None:
        arguments += [
            "--convergence-policy",
            _choice(
                convergence_policy,
                "convergence_policy",
                ("store_final", "store_nan", "discard", "error"),
            ),
        ]
    arguments.append(uri)
    result = _run_cooler(
        "balance",
        arguments,
        cwd=path.parent,
        outputs=[path],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    stored: bool | None = None
    if result["status"] == "completed":
        check = cooler_balance_check(
            uri, name=column, cooler_path=cooler_path, timeout_seconds=300.0, verbose=False
        )
        stored = check["balanced"]
    result.update({"cool_uri": uri, "weights_column": column, "weights_stored": stored})
    return result


def cooler_balance_check(
    cool_uri: str,
    *,
    name: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report whether a cooler already has a balancing weight column (cooler balance --check)."""

    path, uri = _cool_input(cool_uri, "cool_uri")
    column = "weight" if name is None else _text(name, "name")
    arguments = ["--check"]
    if name is not None:
        arguments += ["--name", column]
    arguments.append(uri)
    result = _run_cooler(
        "balance",
        arguments,
        cwd=path.parent,
        outputs=[],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    tail = result["log_tail"]
    balanced: bool | None = None
    if "is balanced" in tail:
        balanced = True
    elif "column found" in tail:
        balanced = False
    # cooler exits 1 both for 'not balanced' and for unreadable files; only the message tells which.
    result["status"] = "completed" if balanced is not None else "failed"
    result.update({"cool_uri": uri, "weights_column": column, "balanced": balanced})
    return result


# --------------------------------------------------------------------------- cload


def _cload(
    mode: str,
    bins: str,
    pairs_path: Path,
    output_path: str,
    extra: list[str],
    *,
    metadata_path: str | None,
    assembly: str | None,
    append: bool,
    cooler_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    spec = _bins_spec(bins, "bins")
    target, uri = _cool_output(output_path, "output_path", append)
    arguments = [
        mode,
        *_optional_metadata(metadata_path, assembly),
        *extra,
        spec,
        str(pairs_path),
        uri,
    ]
    result = _run_cooler(
        "cload",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result["subcommand"] = f"cload {mode}"
    result["output_uri"] = uri
    return result


def cooler_cload_pairs(
    bins: str,
    pairs_path: str,
    output_path: str,
    *,
    chrom1: int,
    pos1: int,
    chrom2: int,
    pos2: int,
    zero_based: bool = False,
    comment_char: str | None = None,
    no_symmetric_upper: bool = False,
    input_copy_status: Literal["unique", "duplex"] | None = None,
    fields: Sequence[str] | None = None,
    chunksize: int | None = None,
    mergebuf: int | None = None,
    max_merge: int | None = None,
    temp_dir: str | None = None,
    no_delete_temp: bool = False,
    storage_options: str | None = None,
    metadata_path: str | None = None,
    assembly: str | None = None,
    append: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin an unsorted text file of read pairs (optionally gzipped) into a cooler (cload pairs).

    ``bins`` is ``<chromsizes file>:<bin size in bp>`` or the path of a BED file with the bins.
    Field numbers are one-based. ``fields`` entries look like ``score=7:dtype=float,agg=mean``. With
    ``append=True`` the cooler is added to an existing file (``output_path`` may be
    ``file.mcool::resolutions/10000``). Options left as None keep cooler's defaults (comment
    character ``#``, input_copy_status unique, chunksize 15,000,000). Reading from standard input is
    not supported.
    """

    source = _existing(pairs_path, "pairs_path")
    extra: list[str] = []
    for flag, number, label in (
        ("--chrom1", chrom1, "chrom1"),
        ("--pos1", pos1, "pos1"),
        ("--chrom2", chrom2, "chrom2"),
        ("--pos2", pos2, "pos2"),
    ):
        extra += [flag, str(_integer(number, label))]
    if zero_based:
        extra.append("--zero-based")
    if comment_char is not None:
        extra += ["--comment-char", _single_character(comment_char, "comment_char")]
    if no_symmetric_upper:
        extra.append("--no-symmetric-upper")
    if input_copy_status is not None:
        extra += [
            "--input-copy-status",
            _choice(input_copy_status, "input_copy_status", ("unique", "duplex")),
        ]
    for field in _fields(fields, "fields", column_number=True):
        extra += ["--field", field]
    for flag, value, label in (
        ("--chunksize", chunksize, "chunksize"),
        ("--mergebuf", mergebuf, "mergebuf"),
        ("--max-merge", max_merge, "max_merge"),
    ):
        if value is not None:
            extra += [flag, str(_integer(value, label))]
    if temp_dir is not None:
        extra += ["--temp-dir", str(_existing_folder(temp_dir, "temp_dir"))]
    if no_delete_temp:
        extra.append("--no-delete-temp")
    extra += _storage_options(storage_options)
    if append:
        extra.append("--append")
    return _cload(
        "pairs",
        bins,
        source,
        output_path,
        extra,
        metadata_path=metadata_path,
        assembly=assembly,
        append=append,
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )


def cooler_cload_tabix(
    bins: str,
    pairs_path: str,
    output_path: str,
    *,
    chrom2: int | None = None,
    pos2: int | None = None,
    zero_based: bool = False,
    nproc: int | None = None,
    max_split: int | None = None,
    metadata_path: str | None = None,
    assembly: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin a bgzip-compressed, tabix-indexed pairs file into a cooler (cooler cload tabix).

    ``pairs_path`` needs its ``.tbi`` (or ``.csi``) index next to it, built on chrom1 and pos1 (see
    ``cooler_csort``). The chromosome order of ``bins`` must agree with the order in the file. Field
    numbers are one-based; cooler's defaults are chrom2 4 and pos2 5, nproc 8 and max_split 2 when
    left as None.
    """

    source = _existing(pairs_path, "pairs_path")
    if not (Path(str(source) + ".tbi").is_file() or Path(str(source) + ".csi").is_file()):
        raise ConfigurationError(f"No tabix index (.tbi or .csi) found next to {source}.")
    extra: list[str] = []
    for flag, value, label in (
        ("--chrom2", chrom2, "chrom2"),
        ("--pos2", pos2, "pos2"),
        ("--nproc", nproc, "nproc"),
        ("--max-split", max_split, "max_split"),
    ):
        if value is not None:
            extra += [flag, str(_integer(value, label))]
    if zero_based:
        extra.append("--zero-based")
    return _cload(
        "tabix",
        bins,
        source,
        output_path,
        extra,
        metadata_path=metadata_path,
        assembly=assembly,
        append=False,
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )


def cooler_cload_pairix(
    bins: str,
    pairs_path: str,
    output_path: str,
    *,
    zero_based: bool = False,
    nproc: int | None = None,
    max_split: int | None = None,
    block_char: str | None = None,
    metadata_path: str | None = None,
    assembly: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin a bgzip-compressed, pairix-indexed pairs file into a cooler (cooler cload pairix).

    ``pairs_path`` needs its ``.px2`` index next to it (see ``cooler_csort`` with index pairix).
    cooler needs the ``pypairix`` Python package for this mode; without it the run is reported as
    failed with the import error in ``log_tail``. Defaults when left as None: nproc 8, max_split 2,
    block_char ``|``.
    """

    source = _existing(pairs_path, "pairs_path")
    if not Path(str(source) + ".px2").is_file():
        raise ConfigurationError(f"No pairix index (.px2) found next to {source}.")
    extra: list[str] = []
    for flag, value, label in (
        ("--nproc", nproc, "nproc"),
        ("--max-split", max_split, "max_split"),
    ):
        if value is not None:
            extra += [flag, str(_integer(value, label))]
    if zero_based:
        extra.append("--zero-based")
    if block_char is not None:
        extra += ["--block-char", _single_character(block_char, "block_char")]
    return _cload(
        "pairix",
        bins,
        source,
        output_path,
        extra,
        metadata_path=metadata_path,
        assembly=assembly,
        append=False,
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )


def cooler_cload_hiclib(
    bins: str,
    pairs_path: str,
    output_path: str,
    *,
    chunksize: int | None = None,
    metadata_path: str | None = None,
    assembly: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin a hiclib HDF5 contact list (chrms1, cuts1, chrms2, cuts2) into a cooler (cload hiclib).

    The pairs must be sorted on the first chromosome and upper triangular with respect to the
    chromosome order of ``bins``; positions are zero-based. chunksize defaults to 100,000,000 when
    left as None.
    """

    source = _existing(pairs_path, "pairs_path")
    extra: list[str] = []
    if chunksize is not None:
        extra += ["--chunksize", str(_integer(chunksize, "chunksize"))]
    return _cload(
        "hiclib",
        bins,
        source,
        output_path,
        extra,
        metadata_path=metadata_path,
        assembly=assembly,
        append=False,
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )


# ------------------------------------------------------------------ coarsen, merge, zoomify, load


def cooler_coarsen(
    cool_uri: str,
    output_path: str,
    *,
    factor: int | None = None,
    nproc: int | None = None,
    chunksize: int | None = None,
    fields: Sequence[str] | None = None,
    append: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coarsen a cooler to a lower resolution by pooling k-by-k pixel tiles (cooler coarsen).

    ``fields`` names the value columns to aggregate, e.g. ``["count",
    "score:dtype=float,agg=mean"]``; the default is the ``count`` column summed. factor defaults to
    2, nproc to 1 and chunksize to 10,000,000 when left as None. ``append=True`` adds the result to
    an existing file (``output_path`` may carry ``::group``).
    """

    _, source = _cool_input(cool_uri, "cool_uri")
    target, uri = _cool_output(output_path, "output_path", append)
    arguments: list[str] = []
    for flag, value, label in (
        ("--factor", factor, "factor"),
        ("--nproc", nproc, "nproc"),
        ("--chunksize", chunksize, "chunksize"),
    ):
        if value is not None:
            arguments += [flag, str(_integer(value, label, 2 if label == "factor" else 1))]
    for field in _fields(fields, "fields", column_number=False):
        arguments += ["--field", field]
    if append:
        arguments.append("--append")
    arguments += ["--out", uri, source]
    result = _run_cooler(
        "coarsen",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result.update({"input_uri": source, "output_uri": uri})
    return result


def cooler_merge(
    output_path: str,
    input_paths: Sequence[str],
    *,
    chunksize: int | None = None,
    fields: Sequence[str] | None = None,
    append: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge coolers that share identical bins into one by summing their pixels (cooler merge).

    ``input_paths`` are cooler files or URIs with identical chromosomes and bins. ``fields`` are as
    in ``cooler_coarsen``. chunksize (merge buffer in pixel rows) defaults to 20,000,000 when left
    as None.
    """

    if isinstance(input_paths, str) or len(input_paths) < 2:
        raise ConfigurationError("input_paths must list at least two coolers to merge.")
    sources = [_cool_input(item, "input_paths")[1] for item in input_paths]
    target, uri = _cool_output(output_path, "output_path", append)
    arguments: list[str] = []
    if chunksize is not None:
        arguments += ["--chunksize", str(_integer(chunksize, "chunksize"))]
    for field in _fields(fields, "fields", column_number=False):
        arguments += ["--field", field]
    if append:
        arguments.append("--append")
    arguments += [uri, *sources]
    result = _run_cooler(
        "merge",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result.update({"output_uri": uri, "input_uris": sources})
    return result


_BALANCE_FLAGS: dict[str, tuple[str, str]] = {
    "cis_only": ("--cis-only", "flag"),
    "trans_only": ("--trans-only", "flag"),
    "ignore_diags": ("--ignore-diags", "int"),
    "ignore_dist": ("--ignore-dist", "int"),
    "mad_max": ("--mad-max", "int"),
    "min_nnz": ("--min-nnz", "int"),
    "min_count": ("--min-count", "int"),
    "nproc": ("--nproc", "int"),
    "max_iters": ("--max-iters", "int"),
    "tol": ("--tol", "float"),
}


def _balance_arguments(options: Mapping[str, Any] | None) -> str | None:
    if not options:
        return None
    parts: list[str] = []
    for key, value in options.items():
        if key not in _BALANCE_FLAGS:
            raise ConfigurationError(
                f"Unsupported balance option {key!r}; use one of {sorted(_BALANCE_FLAGS)}."
            )
        flag, kind = _BALANCE_FLAGS[key]
        if kind == "flag":
            if value is True:
                parts.append(flag)
            elif value is not False:
                raise ConfigurationError(f"balance option {key} must be true or false.")
        elif kind == "int":
            parts += [flag, str(_integer(value, key, 0))]
        else:
            parts += [flag, repr(_number(value, key, positive=True))]
    return " ".join(parts) or None


def cooler_zoomify(
    cool_uri: str,
    output_path: str,
    *,
    resolutions: str | None = None,
    nproc: int | None = None,
    chunksize: int | None = None,
    balance: bool = False,
    balance_options: Mapping[str, Any] | None = None,
    base_uris: Sequence[str] | None = None,
    fields: Sequence[str] | None = None,
    legacy: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a multi-resolution cooler (.mcool) by recursive coarsening (cooler zoomify).

    ``resolutions`` is a comma-separated list such as ``1000,5000,25000`` or a progression: ``B``
    (binary, steps of 2), ``N`` (nice, 1-2-5 steps), ``<start>N`` or ``4DN``; the default is ``B``.
    The progression ``<start>B`` (e.g. ``1000B``) is documented by cooler but fails in 0.10.4
    (zoomify tests ``res.endswith("n")`` twice, so ``1000b`` reaches ``int("1000b")``; reproduced
    with ``-r 100B``): use ``B`` or an explicit list. ``legacy=True`` writes the integer-labelled
    zoom levels of the legacy layout. ``balance=True`` balances every zoom level;
    ``balance_options`` passes cooler balance settings (keys cis_only, trans_only, ignore_diags,
    ignore_dist, mad_max, min_nnz, min_count, nproc, max_iters, tol). ``base_uris`` are extra base
    coolers to aggregate from.
    """

    _, source = _cool_input(cool_uri, "cool_uri")
    target = _new(output_path, "output_path")
    arguments: list[str] = []
    if resolutions is not None:
        text = _text(resolutions, "resolutions").replace(" ", "")
        if not all(_RESOLUTION_TOKEN.match(token) for token in text.split(",")):
            raise ConfigurationError(
                "resolutions must be integers or progressions (B, N, <start>N, <start>B, 4DN) "
                "separated by commas."
            )
        arguments += ["--resolutions", text]
    for flag, value, label in (
        ("--nproc", nproc, "nproc"),
        ("--chunksize", chunksize, "chunksize"),
    ):
        if value is not None:
            arguments += [flag, str(_integer(value, label))]
    if balance:
        arguments.append("--balance")
        balance_text = _balance_arguments(balance_options)
        if balance_text is not None:
            arguments.append(f"--balance-args={balance_text}")
    elif balance_options:
        raise ConfigurationError("balance_options needs balance=True.")
    if base_uris is not None:
        if isinstance(base_uris, str):
            raise ConfigurationError("base_uris must be a list of cooler paths or URIs.")
        for item in base_uris:
            arguments += ["--base-uri", _cool_input(item, "base_uris")[1]]
    for field in _fields(fields, "fields", column_number=False):
        arguments += ["--field", field]
    if legacy:
        arguments.append("--legacy")
    arguments += ["--out", str(target), source]
    result = _run_cooler(
        "zoomify",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result.update(
        {
            "input_uri": source,
            "output_path": str(target),
            "legacy": bool(legacy),
            "balanced": bool(balance),
        }
    )
    return result


def cooler_load(
    bins: str,
    pixels_path: str,
    output_path: str,
    *,
    input_format: Literal["coo", "bg2"],
    metadata_path: str | None = None,
    assembly: str | None = None,
    fields: Sequence[str] | None = None,
    count_as_float: bool = False,
    one_based: bool = False,
    comment_char: str | None = None,
    no_symmetric_upper: bool = False,
    input_copy_status: Literal["unique", "duplex"] | None = None,
    chunksize: int | None = None,
    mergebuf: int | None = None,
    temp_dir: str | None = None,
    max_merge: int | None = None,
    no_delete_temp: bool = False,
    storage_options: str | None = None,
    append: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a cooler from a pre-binned matrix in COO or BG2 text format (cooler load).

    ``input_format`` is ``coo`` (bin1_id, bin2_id, count) or ``bg2`` (chrom1, start1, end1, chrom2,
    start2, end2, count). ``bins`` is ``<chromsizes file>:<bin size in bp>`` or a BED file.
    ``fields`` entries look like ``score=4:dtype=float``; once any field is given the count column
    is loaded only if listed too (``count=7`` or ``count:dtype=float``). ``count_as_float`` stores
    counts as floats. The file may be gzipped; standard input is not supported. Options left as None
    keep cooler's defaults.
    """

    spec = _bins_spec(bins, "bins")
    source = _existing(pixels_path, "pixels_path")
    target, uri = _cool_output(output_path, "output_path", append)
    arguments = ["--format", _choice(input_format, "input_format", ("coo", "bg2"))]
    arguments += _optional_metadata(metadata_path, assembly)
    for field in _fields(fields, "fields", column_number=True):
        arguments += ["--field", field]
    if count_as_float:
        arguments.append("--count-as-float")
    if one_based:
        arguments.append("--one-based")
    if comment_char is not None:
        arguments += ["--comment-char", _single_character(comment_char, "comment_char")]
    if no_symmetric_upper:
        arguments.append("--no-symmetric-upper")
    if input_copy_status is not None:
        arguments += [
            "--input-copy-status",
            _choice(input_copy_status, "input_copy_status", ("unique", "duplex")),
        ]
    for flag, value, label in (
        ("--chunksize", chunksize, "chunksize"),
        ("--mergebuf", mergebuf, "mergebuf"),
        ("--max-merge", max_merge, "max_merge"),
    ):
        if value is not None:
            arguments += [flag, str(_integer(value, label))]
    if temp_dir is not None:
        arguments += ["--temp-dir", str(_existing_folder(temp_dir, "temp_dir"))]
    if no_delete_temp:
        arguments.append("--no-delete-temp")
    arguments += _storage_options(storage_options)
    if append:
        arguments.append("--append")
    arguments += [spec, str(source), uri]
    result = _run_cooler(
        "load",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result["output_uri"] = uri
    return result


# ------------------------------------------------------- text outputs: bins, digest, dump, info


def _bin_table_options(header: bool, rel_ids: int | None) -> list[str]:
    arguments: list[str] = []
    if header:
        arguments.append("--header")
    if rel_ids is not None:
        if isinstance(rel_ids, bool) or rel_ids not in (0, 1):
            raise ConfigurationError("rel_ids must be 0, 1 or None.")
        arguments += ["--rel-ids", str(rel_ids)]
    return arguments


def cooler_makebins(
    chromsizes_path: str,
    bin_size: int,
    output_path: str,
    *,
    header: bool = False,
    rel_ids: Literal[0, 1] | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write fixed-width genomic bins as a BED file from a chromsizes file (cooler makebins).

    ``rel_ids`` adds a column of per-chromosome bin numbers, zero-based (0) or one-based (1).
    """

    sizes = _existing(chromsizes_path, "chromsizes_path")
    size = _integer(bin_size, "bin_size")
    target = _new(output_path, "output_path")
    arguments = [*_bin_table_options(header, rel_ids), "--out", str(target), str(sizes), str(size)]
    result = _run_cooler(
        "makebins",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    if result["status"] == "completed":
        count, preview = _text_preview(target)
        result.update({"line_count": count, "preview": preview})
    return result


def cooler_digest(
    chromsizes_path: str,
    fasta_path: str,
    enzyme: str,
    output_path: str,
    *,
    header: bool = False,
    rel_ids: Literal[0, 1] | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write restriction-fragment bins of a genome as a BED file (cooler digest).

    ``fasta_path`` is a genome FASTA file or a folder of ``*.fa`` / ``*.fasta`` files
    (uncompressed); pyfaidx may create a ``.fai`` index next to the FASTA. ``enzyme`` is a Biopython
    restriction enzyme name such as HindIII or DpnII. Chromosomes follow the order of the chromsizes
    file.
    """

    sizes = _existing(chromsizes_path, "chromsizes_path")
    fasta = (
        Path(fasta_path).expanduser().resolve()
        if isinstance(fasta_path, (str, os.PathLike))
        else None
    )
    if fasta is None or not (fasta.is_file() or fasta.is_dir()):
        raise ConfigurationError(
            f"fasta_path must be an existing FASTA file or folder: {fasta_path}"
        )
    if not isinstance(enzyme, str) or not _ENZYME.match(enzyme):
        raise ConfigurationError("enzyme must be a restriction enzyme name such as HindIII.")
    target = _new(output_path, "output_path")
    arguments = [
        *_bin_table_options(header, rel_ids),
        "--out",
        str(target),
        str(sizes),
        str(fasta),
        enzyme,
    ]
    result = _run_cooler(
        "digest",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    if result["status"] == "completed":
        count, preview = _text_preview(target)
        result.update({"line_count": count, "preview": preview})
    return result


def cooler_dump(
    cool_uri: str,
    output_path: str,
    *,
    table: Literal["chroms", "bins", "pixels"] = "pixels",
    columns: Sequence[str] | None = None,
    header: bool = False,
    na_rep: str | None = None,
    float_format: str | None = None,
    region: str | None = None,
    region2: str | None = None,
    fill_lower: bool = False,
    balanced: bool = False,
    join: bool = False,
    annotate: Sequence[str] | None = None,
    one_based_ids: bool = False,
    one_based_starts: bool = False,
    chunksize: int | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dump a cooler table (chroms, bins or pixels) to tab-separated text (cooler dump).

    For pixels, ``region`` / ``region2`` select a UCSC-style query box (``chr1:10,000-20,000``);
    ``balanced`` adds a ``balanced`` column (weight1 * weight2 * count), ``join`` replaces bin ids
    with chromosome coordinates, ``annotate`` joins extra bin-table columns (suffixed 1 and 2),
    ``fill_lower`` generates the implicit lower triangle of symmetric-upper coolers.
    ``one_based_ids`` and ``one_based_starts`` shift ids or starts, but cooler applies them only
    together with ``join``, ``balanced`` or ``annotate``. ``columns`` applies to chroms and bins.
    float_format defaults to ``g``; ``na_rep`` must be non-empty (cooler's default is empty). The
    result reports ``line_count`` and a short ``preview``.
    """

    _, source = _cool_input(cool_uri, "cool_uri")
    target = _new(output_path, "output_path")
    arguments = ["--table", _choice(table, "table", ("chroms", "bins", "pixels"))]
    column_text = _comma_names(columns, "columns")
    if column_text is not None:
        arguments += ["--columns", column_text]
    if header:
        arguments.append("--header")
    if na_rep is not None:
        arguments += ["--na-rep", _text(na_rep, "na_rep", option_value=True)]
    if float_format is not None:
        if not isinstance(float_format, str) or not _FLOAT_FORMAT.match(float_format):
            raise ConfigurationError("float_format must look like 'g', '.12g' or '03.2f'.")
        arguments += ["--float-format", float_format]
    if region is not None:
        arguments += ["--range", _region(region, "region")]
    if region2 is not None:
        if region is None:
            raise ConfigurationError("region2 needs region.")
        arguments += ["--range2", _region(region2, "region2")]
    if fill_lower:
        arguments.append("--fill-lower")
    if balanced:
        arguments.append("--balanced")
    if join:
        arguments.append("--join")
    annotate_text = _comma_names(annotate, "annotate")
    if annotate_text is not None:
        arguments += ["--annotate", annotate_text]
    if one_based_ids:
        arguments.append("--one-based-ids")
    if one_based_starts:
        arguments.append("--one-based-starts")
    if chunksize is not None:
        arguments += ["--chunksize", str(_integer(chunksize, "chunksize"))]
    arguments += ["--out", str(target), source]
    result = _run_cooler(
        "dump",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result.update({"input_uri": source, "table": table, "line_count": None, "preview": []})
    if result["status"] == "completed":
        count, preview = _text_preview(target)
        result.update({"line_count": count, "preview": preview})
    return result


def cooler_info(
    cool_uri: str,
    *,
    field: str | None = None,
    metadata: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the info block of a cooler: bin size, bin and chromosome counts, pixel sum.

    ``field`` returns one info value; ``metadata=True`` returns the user metadata instead. The
    parsed JSON is in ``info``, ``value`` or ``metadata`` of the result.
    """

    _, uri = _cool_input(cool_uri, "cool_uri")
    arguments: list[str] = []
    if metadata and field is not None:
        raise ConfigurationError("Give either field or metadata, not both.")
    if metadata:
        arguments.append("--metadata")
    if field is not None:
        arguments += ["--field", _text(field, "field")]
    with tempfile.TemporaryDirectory(prefix="dnakit-cooler-info-") as scratch:
        out = Path(scratch) / "info.txt"
        result = _run_cooler(
            "info",
            [*arguments, "--out", str(out), uri],
            cwd=Path(scratch),
            outputs=[],
            cooler_path=cooler_path,
            timeout_seconds=timeout_seconds,
        )
        text = out.read_text(encoding="utf-8") if out.is_file() else ""
    result.update({"cool_uri": uri, "info": None, "metadata": None, "value": None})
    if result["status"] == "completed":
        if field is not None:
            stripped = text.strip()
            try:
                result["value"] = json.loads(stripped)
            except ValueError:
                result["value"] = stripped
        else:
            try:
                parsed = json.loads(text)
            except ValueError as exc:
                raise BackendExecutionError("cooler info did not return valid JSON.") from exc
            result["metadata" if metadata else "info"] = parsed
    return result


# --------------------------------------------------------------------------- file operations


def cooler_ls(
    cool_path: str,
    *,
    long: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """List all coolers inside an HDF5 file, e.g. the resolutions of an .mcool (cooler ls).

    ``long=True`` also reports the bin size of each cooler (None for variable-size bins).
    """

    path = _existing(cool_path, "cool_path")
    arguments = ["--long", str(path)] if long else [str(path)]
    result = _run_cooler(
        "ls",
        arguments,
        cwd=path.parent,
        outputs=[],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
        tail_chars=_BIG_TAIL,
    )
    coolers: list[str] = []
    bin_sizes: dict[str, int | None] = {}
    for line in result["log_tail"].splitlines():
        uri, _, size_text = line.partition("\t")
        if not uri.startswith(f"{path}::"):
            continue
        coolers.append(uri)
        if long:
            size_text = size_text.strip().replace(",", "")
            bin_sizes[uri] = int(size_text) if size_text.isdigit() else None
    result.update({"cool_path": str(path), "coolers": coolers})
    if long:
        result["bin_sizes"] = bin_sizes
    return result


def _transfer(
    subcommand: str,
    source_uri: str,
    destination_uri: str,
    extra: list[str],
    cooler_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    source_path, source = _cool_input(source_uri, "source_uri")
    file_text, group = _split_uri(destination_uri, "destination_uri")
    destination = Path(file_text).expanduser().resolve()
    if not destination.parent.is_dir():
        raise ConfigurationError(
            f"The folder for destination_uri does not exist: {destination.parent}"
        )
    target = f"{destination}::{group}" if group else str(destination)
    result = _run_cooler(
        subcommand,
        [*extra, source, target],
        cwd=destination.parent,
        outputs=[],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    if destination.exists():
        result["outputs"] = [str(destination)]
    result.update(
        {"source_uri": source, "destination_uri": target, "source_file": str(source_path)}
    )
    return result


def cooler_cp(
    source_uri: str,
    destination_uri: str,
    *,
    overwrite: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy a cooler from one file to another or to another group of the same file (cooler cp).

    URIs are ``file`` or ``file::/group``. ``overwrite=True`` passes ``-w``, which truncates and
    replaces the whole destination file if it exists. When source and destination are in the same
    file, cooler 0.10.4 cannot honour it and fails with "unable to truncate a file which is already
    open" (reproduced with ``cooler cp -w f.cool::/a f.cool::/b``).
    """

    return _transfer(
        "cp",
        source_uri,
        destination_uri,
        ["--overwrite"] if overwrite else [],
        cooler_path,
        timeout_seconds,
    )


def cooler_mv(
    source_uri: str,
    destination_uri: str,
    *,
    overwrite: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rename a cooler group within the same HDF5 file (cooler mv).

    ``overwrite=True`` passes ``-w``, which cooler 0.10.4 cannot honour inside one file (see
    ``cooler_cp``); it is accepted for completeness.
    """

    return _transfer(
        "mv",
        source_uri,
        destination_uri,
        ["--overwrite"] if overwrite else [],
        cooler_path,
        timeout_seconds,
    )


def cooler_ln(
    source_uri: str,
    destination_uri: str,
    *,
    overwrite: bool = False,
    soft: bool = False,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Link a cooler instead of copying it (cooler ln).

    A hard link within one file by default; ``soft=True`` makes a soft link in the same file, or an
    external link when the destination is another file. ``overwrite=True`` passes ``-w``, which
    truncates the destination file (and fails inside one file, see ``cooler_cp``). Do not link the
    root group ``file::/`` into one of its own sub-groups: that makes a cycle and cooler then
    fails to list the file (RecursionError, reproduced with ``cooler ln f.cool::/ f.cool::/x``).
    """

    extra: list[str] = []
    if overwrite:
        extra.append("--overwrite")
    if soft:
        extra.append("--soft")
    return _transfer("ln", source_uri, destination_uri, extra, cooler_path, timeout_seconds)


def _tree_command(
    subcommand: str, uri: str, level: int | None, cooler_path: str | None, timeout_seconds: float
) -> dict[str, Any]:
    path, full = _cool_input(uri, "uri")
    arguments: list[str] = []
    if level is not None:
        arguments += ["--level", str(_integer(level, "level", 0))]
    arguments.append(full)
    result = _run_cooler(
        subcommand,
        arguments,
        cwd=path.parent,
        outputs=[],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
        tail_chars=_BIG_TAIL,
    )
    result["uri"] = full
    return result


def cooler_tree(
    uri: str,
    *,
    level: int | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the HDF5 data hierarchy of a cooler file or group as text (cooler tree).

    ``level`` limits the depth. The text is the value of ``cooler.fileops.pprint_data_tree``; it is
    in ``tree``.
    """

    result = _tree_command("tree", uri, level, cooler_path, timeout_seconds)
    result["tree"] = result["log_tail"].rstrip("\n") if result["status"] == "completed" else None
    return result


def cooler_attrs(
    uri: str,
    *,
    level: int | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the HDF5 attributes of a cooler file or group as YAML and a dict (cooler attrs).

    ``level`` limits the depth. The text is the value of ``cooler.fileops.pprint_attr_tree``; it is
    in ``attributes_yaml`` and, when PyYAML is installed, parsed into ``attributes``.
    """

    result = _tree_command("attrs", uri, level, cooler_path, timeout_seconds)
    result["attributes_yaml"] = None
    result["attributes"] = None
    if result["status"] == "completed":
        text = result["log_tail"]
        result["attributes_yaml"] = text
        try:
            import yaml  # type: ignore[import-untyped]

            result["attributes"] = yaml.safe_load(text)
        except Exception:
            result["attributes"] = None
    return result


# --------------------------------------------------------------------------- csort


def _field_numbers(chrom1: int, pos1: int, chrom2: int, pos2: int) -> dict[str, int]:
    numbers = {
        "C1": _integer(chrom1, "chrom1"),
        "P1": _integer(pos1, "pos1"),
        "C2": _integer(chrom2, "chrom2"),
        "P2": _integer(pos2, "pos2"),
    }
    if len(set(numbers.values())) != 4:
        raise ConfigurationError(
            "chrom1, pos1, chrom2 and pos2 must be four different field numbers."
        )
    return numbers


def _sort_keys(index: str, fields: Mapping[str, int]) -> list[str]:
    """The sort keys of ``cooler csort``: positional for tabix, block for pairix."""

    c1, p1, c2, p2 = fields["C1"], fields["P1"], fields["C2"], fields["P2"]
    if index == "tabix":
        return [f"-k{c1},{c1}", f"-k{p1},{p1}n", f"-k{c2},{c2}", f"-k{p2},{p2}n"]
    return [f"-k{c1},{c1}", f"-k{c2},{c2}", f"-k{p1},{p1}n", f"-k{p2},{p2}n"]


def _index_arguments(
    index: str, fields: Mapping[str, int], zero_based: bool, output: str
) -> list[str]:
    c1, p1, c2, p2 = fields["C1"], fields["P1"], fields["C2"], fields["P2"]
    if index == "tabix":
        arguments = ["-f", f"-s{c1}", f"-b{p1}", f"-e{p1}"]
    else:
        arguments = ["-f", f"-s{c1}", f"-d{c2}", f"-b{p1}", f"-e{p1}", f"-u{p2}", f"-v{p2}"]
    if zero_based:
        arguments.append("-0")
    return [*arguments, output]


def _read_command(infile: str, pigz: str | None, gzip_program: str | None) -> list[str]:
    """The decompression command of ``cooler csort``: pigz, else gzip, for ``.gz`` files."""

    if pigz is not None:
        return [pigz, "-dc", infile]
    if gzip_program is not None:
        return [gzip_program, "-dc", infile]
    raise BackendUnavailableError(
        "No gzip decompressor (pigz or gzip) was found on PATH.", hint=_GZIP_HINT
    )


def cooler_csort(
    pairs_path: str,
    chromosomes_path: str,
    output_path: str,
    *,
    chrom1: int,
    pos1: int,
    chrom2: int,
    pos2: int,
    index: Literal["tabix", "pairix"] = "pairix",
    nproc: int | None = None,
    zero_based: bool = False,
    sep: str | None = None,
    comment_char: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = _DEFAULT_TIMEOUT,
    verbose: bool = True,
) -> dict[str, Any]:
    """Orient, sort, bgzip and index a contact list for cload tabix or pairix (cooler csort).

    Mates are ordered so every contact is upper triangular for the chromosome order of
    ``chromosomes_path`` (pairs on other chromosomes are dropped), rows are sorted (positional for
    tabix, block for pairix), compressed with bgzip and indexed. Field numbers are one-based. Needs
    sort, bgzip and tabix or pairix on PATH (or next to cooler); missing programs raise
    BackendUnavailableError. sep defaults to a tab and comment_char to ``#``. The result lists the
    sorted file and its index (``.tbi`` or ``.px2``). Reading standard input and the free
    ``--sort-options`` text are not supported (see ``cooler_csort_sort_pairs`` for sort tuning).
    """

    source = _existing(pairs_path, "pairs_path")
    chromosomes = _existing(chromosomes_path, "chromosomes_path")
    fields = _field_numbers(chrom1, pos1, chrom2, pos2)
    choice = _choice(index, "index", ("tabix", "pairix"))
    target = _new(output_path, "output_path")
    if target == source:
        raise ConfigurationError("output_path must differ from pairs_path.")
    arguments = [
        "--chrom1",
        str(fields["C1"]),
        "--chrom2",
        str(fields["C2"]),
        "--pos1",
        str(fields["P1"]),
        "--pos2",
        str(fields["P2"]),
        "--index",
        choice,
    ]
    if nproc is not None:
        arguments += ["--nproc", str(_integer(nproc, "nproc"))]
    if zero_based:
        arguments.append("--zero-based")
    if sep is not None:
        arguments += ["--sep", _single_character(sep, "sep")]
    if comment_char is not None:
        arguments += ["--comment-char", _single_character(comment_char, "comment_char")]
    arguments += ["--out", str(target), str(source), str(chromosomes)]
    program = _find_program("cooler", cooler_path, _COOLER_HINT)
    _tool_on_path("sort", program, _SORT_HINT)
    _tool_on_path("bgzip", program, _HTSLIB_HINT)
    _tool_on_path(choice, program, _HTSLIB_HINT if choice == "tabix" else _PAIRIX_HINT)
    suffix = ".tbi" if choice == "tabix" else ".px2"
    index_file = Path(str(target) + suffix)
    result = _run_cooler(
        "csort",
        arguments,
        cwd=target.parent,
        outputs=[target, index_file],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
    )
    result.update({"pairs_path": str(source), "index": choice, "index_file": str(index_file)})
    return result


def cooler_csort_read_pairs(
    pairs_path: str,
    output_path: str,
    *,
    pigz_path: str | None = None,
    gzip_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a contact list as plain text, decompressing .gz with pigz or gzip (csort read step).

    As in cooler csort, a file ending in ``.gz`` is read with ``pigz -dc`` when pigz is available
    and with ``gzip -dc`` otherwise; any other file is copied unchanged (cooler pipes it through
    ``cat``). The command is an argument list; its output goes straight to ``output_path``.
    ``command`` is None for a copied file.
    """

    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    if not source.name.endswith(".gz"):
        shutil.copyfile(source, target)
        return {
            "status": "completed",
            "backend": "copy",
            "command": None,
            "return_code": 0,
            "outputs": [str(target)],
            "bytes_written": target.stat().st_size,
            "log_tail": "",
        }
    pigz = pigz_path if pigz_path is not None else shutil.which("pigz")
    gzip_program = gzip_path if gzip_path is not None else shutil.which("gzip")
    for label, value in (("pigz_path", pigz_path), ("gzip_path", gzip_path)):
        if value is not None and not Path(value).is_file():
            raise BackendUnavailableError(f"{label} does not exist: {value}", hint=_GZIP_HINT)
    command = _read_command(str(source), pigz, gzip_program)
    outcome = _run_to_file(command[0], command[1:], target, timeout_seconds)
    ok = outcome["return_code"] == 0 and target.is_file()
    return {
        "status": "completed" if ok else "failed",
        "backend": Path(command[0]).name,
        "command": command,
        "return_code": outcome["return_code"],
        "outputs": [str(target)] if ok else [],
        "bytes_written": target.stat().st_size if ok else 0,
        "log_tail": outcome["log_tail"],
    }


def cooler_csort_flip_pairs(
    pairs_path: str,
    chromosomes_path: str,
    output_path: str,
    *,
    chrom1: int,
    pos1: int,
    chrom2: int,
    pos2: int,
    sep: str = "\t",
    comment_char: str = "#",
    verbose: bool = True,
) -> dict[str, Any]:
    """Orient each contact upper triangular for a chromosome order (csort flip step).

    ``chromosomes_path`` lists the chromosomes in the wanted order (first tab-separated column; a
    chromsizes file works); contacts on other chromosomes are dropped and lines starting with
    ``comment_char`` are skipped. The mates swap when chrom1 ranks after chrom2, or on the same
    chromosome when pos1 > pos2. Field numbers are one-based; a ``.gz`` input is read through gzip.
    The algorithm is the one cooler runs as a ``python -c`` pipeline stage; here it runs in-process
    on the files and returns counts. Output columns are joined with tabs.
    """

    source = _existing(pairs_path, "pairs_path")
    chromosomes = _existing(chromosomes_path, "chromosomes_path")
    fields = _field_numbers(chrom1, pos1, chrom2, pos2)
    target = _new(output_path, "output_path")
    separator = _single_character(sep, "sep").encode()
    comment = _single_character(comment_char, "comment_char").encode()
    c1, c2, p1, p2 = fields["C1"] - 1, fields["C2"] - 1, fields["P1"] - 1, fields["P2"] - 1
    rank: dict[bytes, int] = {}
    with chromosomes.open("rb") as handle:
        for number, line in enumerate(handle, 1):
            name = line.split(b"\t")[0].strip()
            if name:
                rank[name] = number
    if not rank:
        raise ConfigurationError(f"No chromosome names found in {chromosomes}.")
    opener: Any = gzip.open if source.name.endswith(".gz") else open
    counts = {"lines_read": 0, "comment_lines": 0, "dropped": 0, "flipped": 0, "written": 0}
    needed = max(c1, c2, p1, p2) + 1
    try:
        with opener(source, "rb") as instream, target.open("wb") as outstream:
            for number, line in enumerate(instream, 1):
                counts["lines_read"] += 1
                if line.startswith(comment):
                    counts["comment_lines"] += 1
                    continue
                parts = line.strip().split(separator)
                if len(parts) < needed:
                    raise ConfigurationError(
                        f"Line {number} of {source} has {len(parts)} columns; "
                        f"field number {needed} is needed."
                    )
                name1, name2 = parts[c1], parts[c2]
                if name1 in rank and name2 in rank:
                    id1, id2 = rank[name1], rank[name2]
                    try:
                        position1, position2 = int(parts[p1]), int(parts[p2])
                    except ValueError as exc:
                        raise ConfigurationError(
                            f"Line {number} of {source} has a non-integer position."
                        ) from exc
                    if id1 > id2 or (id1 == id2 and position1 > position2):
                        parts[c1], parts[c2] = parts[c2], parts[c1]
                        parts[p1], parts[p2] = parts[p2], parts[p1]
                        counts["flipped"] += 1
                    outstream.write(b"\t".join(parts) + b"\n")
                    counts["written"] += 1
                else:
                    counts["dropped"] += 1
    except (OSError, EOFError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"Could not read the pairs file {source}: {exc}") from exc
    except ConfigurationError:
        target.unlink(missing_ok=True)
        raise
    return {
        "status": "completed",
        "backend": "python",
        "outputs": [str(target)],
        "chromosomes": len(rank),
        **counts,
    }


def cooler_csort_sort_pairs(
    pairs_path: str,
    output_path: str,
    *,
    chrom1: int,
    pos1: int,
    chrom2: int,
    pos2: int,
    index: Literal["tabix", "pairix"] = "pairix",
    nproc: int | None = None,
    buffer_size: str | None = None,
    temp_dir: str | None = None,
    sort_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort a contact list with the system sort, cooler csort's keys and C locale (csort sort step).

    ``index="tabix"`` sorts by chrom1, pos1, chrom2, pos2 (positional); ``index="pairix"`` by
    chrom1, chrom2, pos1, pos2 (block); positions sort numerically, names byte-wise (LC_ALL=C, set
    for the child process only). The input should already be flipped (``cooler_csort_flip_pairs``).
    ``nproc`` adds ``--parallel``, ``buffer_size`` (for example ``50%`` or ``2G``) adds
    ``--buffer-size`` and ``temp_dir`` adds ``-T``; cooler uses ``--parallel=8 --buffer-size=50%``
    when the sort supports it. The result is written with ``sort -o``.
    """

    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    fields = _field_numbers(chrom1, pos1, chrom2, pos2)
    choice = _choice(index, "index", ("tabix", "pairix"))
    arguments = _sort_keys(choice, fields)
    if nproc is not None:
        arguments.append(f"--parallel={_integer(nproc, 'nproc')}")
    if buffer_size is not None:
        if not isinstance(buffer_size, str) or not _SORT_BUFFER.match(buffer_size):
            raise ConfigurationError("buffer_size must look like '50%' or '2G'.")
        arguments.append(f"--buffer-size={buffer_size}")
    if temp_dir is not None:
        arguments += ["-T", str(_existing_folder(temp_dir, "temp_dir"))]
    program = _find_program("sort", sort_path, _SORT_HINT)
    timeout = _timeout(timeout_seconds)
    arguments += ["-o", str(target), str(source)]
    outcome = execute_bounded_command(
        program,
        arguments,
        backend_id="sort",
        cwd=target.parent,
        timeout_seconds=timeout,
        environment_overrides={"LC_ALL": "C"},
        verbose=False,
    )
    ok = outcome.return_code == 0 and target.is_file()
    return {
        "status": "completed" if ok else "failed",
        "backend": "sort",
        "command": ["sort", *arguments],
        "return_code": outcome.return_code,
        "outputs": [str(target)] if ok else [],
        "log_tail": (outcome.output or "")[-_MAX_TAIL:],
    }


def cooler_csort_index_pairs(
    pairs_path: str,
    *,
    chrom1: int,
    pos1: int,
    chrom2: int,
    pos2: int,
    index: Literal["tabix", "pairix"] = "pairix",
    zero_based: bool = False,
    index_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index a sorted, bgzip-compressed contact list with tabix or pairix (csort index step).

    Writes ``<pairs_path>.tbi`` (tabix, on chrom1 and pos1) or ``<pairs_path>.px2`` (pairix, on
    chrom1, chrom2, pos1 and pos2). The file must be bgzip-compressed and sorted with the matching
    order (``cooler_csort_sort_pairs``). ``index_path`` is the tabix or pairix executable when it is
    not on PATH or next to the interpreter.
    """

    source = _existing(pairs_path, "pairs_path")
    with source.open("rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            raise ConfigurationError(
                f"{source} is not bgzip-compressed; compress it with bgzip before indexing."
            )
    fields = _field_numbers(chrom1, pos1, chrom2, pos2)
    choice = _choice(index, "index", ("tabix", "pairix"))
    program = _find_program(choice, index_path, _HTSLIB_HINT if choice == "tabix" else _PAIRIX_HINT)
    arguments = _index_arguments(choice, fields, bool(zero_based), str(source))
    timeout = _timeout(timeout_seconds)
    outcome = execute_bounded_command(
        program,
        arguments,
        backend_id=choice,
        cwd=source.parent,
        timeout_seconds=timeout,
        verbose=False,
    )
    index_file = Path(str(source) + (".tbi" if choice == "tabix" else ".px2"))
    ok = outcome.return_code == 0 and index_file.is_file()
    return {
        "status": "completed" if ok else "failed",
        "backend": choice,
        "command": [choice, *arguments],
        "return_code": outcome.return_code,
        "outputs": [str(index_file)] if ok else [],
        "index_file": str(index_file),
        "log_tail": (outcome.output or "")[-_MAX_TAIL:],
    }


# --------------------------------------------------------------------------- pairs header


def cooler_pairs_header(
    pairs_path: str,
    *,
    comment_char: str = "#",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the leading comment lines of a contact list, plain or gzipped (cload get_header).

    Lines are stripped of trailing spaces and newlines and keep their comment character. Reading
    stops at the first line that does not start with ``comment_char``. The result also reports
    whether data lines follow.
    """

    source = _existing(pairs_path, "pairs_path")
    comment = _single_character(comment_char, "comment_char")
    opener: Any = gzip.open if source.name.endswith(".gz") else open
    header: list[str] = []
    has_data = False
    try:
        with opener(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(comment):
                    header.append(line.strip())
                else:
                    has_data = True
                    break
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        raise ConfigurationError(f"Could not read the pairs file {source}: {exc}") from exc
    return {
        "status": "completed",
        "pairs_path": str(source),
        "header": header,
        "header_line_count": len(header),
        "has_data_lines": has_data,
    }


# --------------------------------------------------------------------------- plots


def _matrix_payload(matrix: Any, np: Any) -> dict[str, Any]:
    """The values behind an image: the matrix when small, else a summary (NaN becomes None)."""

    values = np.asarray(matrix, dtype=float)
    finite = values[np.isfinite(values)]
    summary = {
        "shape": [int(values.shape[0]), int(values.shape[1])]
        if values.ndim == 2
        else list(values.shape),
        "finite_cells": int(finite.size),
        "nan_cells": int(values.size - finite.size),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
        "sum": float(finite.sum()) if finite.size else None,
        "mean": float(finite.mean()) if finite.size else None,
    }
    rows = None
    if values.ndim == 2 and values.size <= _MAX_MATRIX_CELLS:
        rows = [[_plain(float(cell)) for cell in row] for row in values.tolist()]
    return {"summary": summary, "matrix": rows, "matrix_omitted": rows is None}


def _scaled(values: Any, scale: str, np: Any) -> Any:
    out = np.asarray(values, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        if scale == "log2":
            return np.log2(out)
        if scale == "log10":
            return np.log10(out)
    return out


_SCALE_LABELS = {
    "linear": "relative contact frequency",
    "log2": "log 2 ( relative contact frequency )",
    "log10": "log 10 ( relative contact frequency )",
}


def _python_modules() -> tuple[Any, Any, Any, Any, Any]:
    cooler = require("cooler", "bio-cooler")
    np = require("numpy", "bio-analysis")
    figure_module = require("matplotlib.figure", "bio-plots")
    canvas_module = require("matplotlib.backends.backend_agg", "bio-plots")
    util = require("cooler.util", "bio-cooler")
    return cooler, np, util, figure_module, canvas_module


def cooler_show(
    cool_uri: str,
    region: str,
    output_path: str,
    *,
    region2: str | None = None,
    balanced: bool = False,
    dpi: int | None = None,
    scale: Literal["linear", "log2", "log10"] | None = None,
    force: bool = False,
    zmin: float | None = None,
    zmax: float | None = None,
    cmap: str | None = None,
    field: str | None = None,
    cooler_path: str | None = None,
    timeout_seconds: float = 1800.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a contact-matrix region of a cooler to an image file (cooler show --out, headless).

    ``region`` is UCSC-style (``chr1:10,000,000-11,000,000``); ``region2`` gives the column range
    (trans or asymmetric views). The image format follows the extension of ``output_path`` (png,
    jpg, svg, pdf, ps, eps). The colour scale is applied for the image only (cooler defaults: log10,
    cmap YlOrRd, field count). The result returns the image path and the values behind it (``data``:
    summary, and the raw matrix when it has at most 40,000 cells), read with the cooler package when
    it is installed.
    """

    _, source = _cool_input(cool_uri, "cool_uri")
    row_region = _region(region, "region")
    target = _new(output_path, "output_path")
    if target.suffix.lower() not in _IMAGE_SUFFIXES:
        raise ConfigurationError(f"output_path must end with one of {list(_IMAGE_SUFFIXES)}.")
    arguments = ["--out", str(target)]
    if region2 is not None:
        arguments += ["--range2", _region(region2, "region2")]
    if balanced:
        arguments.append("--balanced")
    if dpi is not None:
        arguments += ["--dpi", str(_integer(dpi, "dpi", 10))]
    if scale is not None:
        arguments += ["--scale", _choice(scale, "scale", _SCALES)]
    if force:
        arguments.append("--force")
    if zmin is not None:
        arguments.append(f"--zmin={_number(zmin, 'zmin')!r}")
    if zmax is not None:
        arguments.append(f"--zmax={_number(zmax, 'zmax')!r}")
    if cmap is not None:
        if not isinstance(cmap, str) or not _COLORMAP.match(cmap):
            raise ConfigurationError("cmap must be a matplotlib colormap name such as YlOrRd.")
        arguments += ["--cmap", cmap]
    if field is not None:
        arguments += ["--field", _text(field, "field")]
    arguments += [source, row_region]
    result = _run_cooler(
        "show",
        arguments,
        cwd=target.parent,
        outputs=[target],
        cooler_path=cooler_path,
        timeout_seconds=timeout_seconds,
        environment={"MPLBACKEND": "Agg"},
    )
    result.update(
        {
            "image_path": str(target) if target.is_file() else None,
            "image_bytes": target.stat().st_size if target.is_file() else 0,
            "cool_uri": source,
            "region": row_region,
            "region2": region2 or row_region,
            "balanced": bool(balanced),
            "scale": scale or "log10",
            "field": field or "count",
            "data": None,
        }
    )
    if result["status"] == "completed":
        try:
            cooler, np, _, _, _ = _python_modules()
            clr = cooler.Cooler(source)
            second = region2 or row_region
            rows = clr.extent(row_region)[1] - clr.extent(row_region)[0]
            columns = clr.extent(second)[1] - clr.extent(second)[0]
            if rows * columns > _DATA_CELL_LIMIT:
                result["data"] = {
                    "summary": None,
                    "matrix": None,
                    "matrix_omitted": True,
                    "reason": f"the region has {rows * columns} cells (limit {_DATA_CELL_LIMIT})",
                }
            else:
                matrix = clr.matrix(balance=bool(balanced), field=field or "count").fetch(
                    row_region, second
                )
                result["data"] = {
                    "binsize": clr.info.get("bin-size"),
                    **_matrix_payload(matrix, np),
                }
        except (BackendUnavailableError, OSError, KeyError, ValueError):
            result["data"] = None
    return result


def cooler_show_view(
    cool_uri: str,
    region: str,
    output_path: str,
    *,
    region2: str | None = None,
    view: Sequence[float] | None = None,
    balanced: bool = False,
    field: str = "count",
    scale: Literal["linear", "log2", "log10"] = "log10",
    cmap: str = "YlOrRd",
    zmin: float | None = None,
    zmax: float | None = None,
    dpi: int | None = None,
    max_matrix_size: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the heat map the interactive ``cooler show`` window shows after a zoom or pan.

    The interactive window redraws the matrix for the visible window: it rounds the window outwards
    to whole bins, trims it to the chromosome (the last, partial bin is excluded), and shows a
    placeholder text instead of the matrix when the window holds ``max_matrix_size`` cells or more.
    This function does the same without a display. ``view`` is the window as ``[column start, column
    end, row start, row end]`` in bp on the chromosomes of ``region2`` (default ``region``) and
    ``region``; without ``view`` the whole of ``region`` is drawn. Needs fixed-size bins. Returns
    the image path and the values behind it (``data``: summary, and the matrix when it has at most
    40,000 cells), the rounded window and whether the placeholder was drawn.
    """

    _, source = _cool_input(cool_uri, "cool_uri")
    row_text = _region(region, "region")
    col_text = row_text if region2 is None else _region(region2, "region2")
    target = _new(output_path, "output_path")
    if target.suffix.lower() not in _IMAGE_SUFFIXES:
        raise ConfigurationError(f"output_path must end with one of {list(_IMAGE_SUFFIXES)}.")
    choice = _choice(scale, "scale", _SCALES)
    if not isinstance(field, str) or not field.strip():
        raise ConfigurationError("field must be a non-empty column name.")
    if not isinstance(cmap, str) or not _COLORMAP.match(cmap):
        raise ConfigurationError("cmap must be a matplotlib colormap name such as YlOrRd.")
    limit = _integer(max_matrix_size, "max_matrix_size")
    pixel_dpi = None if dpi is None else _integer(dpi, "dpi", 10)
    window: list[float] | None = None
    if view is not None:
        if (
            isinstance(view, str)  # type: ignore[unreachable]
            or len(view) != 4
            or any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in view)
        ):
            raise ConfigurationError(
                "view must be four numbers: column start, column end, row start, row end."
            )
        if view[0] >= view[1] or view[2] >= view[3]:
            raise ConfigurationError("view must have start < end for the columns and for the rows.")
        window = [view[0], view[1], view[2], view[3]]
    cooler, np, util, figure_module, canvas_module = _python_modules()
    try:
        clr = cooler.Cooler(source)
        try:
            binsize = int(clr.info["bin-size"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ConfigurationError(
                "cooler_show_view needs a cooler with fixed-size bins."
            ) from exc
        sizes = clr.chromsizes
        row_chrom, row_lo, row_hi = util.parse_region(row_text, sizes)
        col_chrom, col_lo, col_hi = util.parse_region(col_text, sizes)
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Could not read the cooler or its regions: {exc}") from exc
    if window is None:
        x_lo, x_hi, y_lo, y_hi = col_lo, col_hi, row_lo, row_hi
        row_region: Any = row_text
        col_region: Any = col_text
    else:
        x_lo = max(0, int(math.floor(window[0] / binsize) * binsize))
        x_hi = min(
            int(math.ceil(window[1] / binsize) * binsize),
            int(math.floor(sizes[col_chrom] / binsize) * binsize),
        )
        y_lo = max(0, int(math.floor(window[2] / binsize) * binsize))
        y_hi = min(
            int(math.ceil(window[3] / binsize) * binsize),
            int(math.floor(sizes[row_chrom] / binsize) * binsize),
        )
        if x_hi <= x_lo or y_hi <= y_lo:
            raise ConfigurationError(
                "view is empty after rounding to whole bins and trimming to the chromosome."
            )
        row_region = (row_chrom, y_lo, y_hi)
        col_region = (col_chrom, x_lo, x_hi)
    try:
        rows = clr.extent(row_region)[1] - clr.extent(row_region)[0]
        columns = clr.extent(col_region)[1] - clr.extent(col_region)[0]
        placeholder = rows * columns >= limit
        matrix = None
        if not placeholder:
            if balanced and "weight" not in clr.bins()[:1].columns:
                raise ConfigurationError(
                    "This cooler has no balancing weights; run cooler_balance first."
                )
            matrix = np.asarray(
                clr.matrix(balance=bool(balanced), field=field).fetch(row_region, col_region)
            )
    except (KeyError, ValueError) as exc:
        raise ConfigurationError(f"Region is not valid for this cooler: {exc}") from exc
    figure = figure_module.Figure(figsize=(11, 10))
    canvas_module.FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    if placeholder:
        axes.plot(
            [x_lo, x_hi, x_hi, x_lo, x_lo, x_hi],
            [y_lo, y_hi, y_lo, y_lo, y_hi, y_hi],
            c="k",
            lw=0.5,
        )
        axes.set_xlim(x_lo, x_hi)
        axes.set_ylim(y_hi, y_lo)
        axes.text(
            0.5,
            0.5,
            "The requested region is too large\nto display at this resolution.",
            horizontalalignment="center",
            verticalalignment="center",
            transform=axes.transAxes,
        )
        data = {"summary": None, "matrix": None, "matrix_omitted": True}
    else:
        image = axes.imshow(
            _scaled(matrix, choice, np),
            interpolation="none",
            extent=[x_lo, x_hi, y_hi, y_lo],
            vmin=zmin,
            vmax=zmax,
            cmap=cmap,
        )
        figure.colorbar(image, ax=axes).set_label(_SCALE_LABELS[choice])
        data = _matrix_payload(matrix, np)
    axes.set_ylabel(f"{row_chrom} coordinate")
    axes.set_xlabel(f"{col_chrom} coordinate")
    try:
        figure.savefig(str(target), dpi=pixel_dpi)
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Could not save the image: {exc}") from exc
    return {
        "status": "completed",
        "backend": "cooler+matplotlib",
        "image_path": str(target),
        "image_bytes": target.stat().st_size,
        "cool_uri": source,
        "row_region": [row_chrom, int(y_lo), int(y_hi)],
        "column_region": [col_chrom, int(x_lo), int(x_hi)],
        "view_applied": window is not None,
        "placeholder": bool(placeholder),
        "binsize": binsize,
        "balanced": bool(balanced),
        "scale": choice,
        "field": field,
        "data": data,
    }


__all__ = [
    "cooler_attrs",
    "cooler_balance",
    "cooler_balance_check",
    "cooler_cload_hiclib",
    "cooler_cload_pairix",
    "cooler_cload_pairs",
    "cooler_cload_tabix",
    "cooler_coarsen",
    "cooler_cp",
    "cooler_csort",
    "cooler_csort_flip_pairs",
    "cooler_csort_index_pairs",
    "cooler_csort_read_pairs",
    "cooler_csort_sort_pairs",
    "cooler_digest",
    "cooler_dump",
    "cooler_info",
    "cooler_ln",
    "cooler_load",
    "cooler_ls",
    "cooler_makebins",
    "cooler_merge",
    "cooler_mv",
    "cooler_pairs_header",
    "cooler_show",
    "cooler_show_view",
    "cooler_tree",
    "cooler_zoomify",
]
_register_public_verbose(globals())
