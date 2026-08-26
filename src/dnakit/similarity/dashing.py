"""Explicit, bounded adapter for Dashing file-level Jaccard computation."""

from __future__ import annotations

import hashlib
import os
import re
import stat
import tempfile
from collections.abc import Iterable
from contextlib import suppress
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Literal, TypeAlias

from dnakit.backends import execute_bounded_command
from dnakit.core import (
    DNA,
    BackendInfo,
    Citation,
    DNARecord,
    DNASequence,
    DNASet,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    Topology,
)
from dnakit.core._json import freeze_mapping
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.similarity.results import (
    DashingJaccardMatrixResult,
    DashingNeighborHit,
    DashingNeighborRow,
    DashingTopKResult,
)

DashingMode: TypeAlias = Literal["sketch", "exact"]
DashingInput: TypeAlias = DNA | DNASequence | DNARecord | str | Path

DEFAULT_MAX_DASHING_ITEMS = 1_000
DEFAULT_MAX_DASHING_INPUT_BYTES = 1_000_000_000
DEFAULT_MAX_DASHING_OUTPUT_BYTES = 100_000_000
DEFAULT_MAX_DASHING_CAPTURE_BYTES = 1_000_000
DEFAULT_MAX_DASHING_SKETCH_MEMORY_BYTES = 512_000_000

_DASHING_VERSION = re.compile(r"Dashing version:\s*(?P<version>[-+._A-Za-z0-9]{1,128})")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_COMMAND_PROFILE = "dashing-dist-jaccard-v1"
_FIXED_SUBCOMMANDS = frozenset({"dist"})
_FIXED_FLAGS = frozenset(
    {
        "--avoid-sorting",
        "--use-full-khash-sets",
        "-C",
        "-F",
        "-O",
        "-S",
        "-k",
        "-o",
        "-p",
    }
)
_JACCARD_TOLERANCE = 1e-6


@dataclass(frozen=True, slots=True)
class _FileSnapshot:
    device: int
    inode: int
    size: int
    modified_ns: int


@dataclass(frozen=True, slots=True)
class _PreparedInput:
    label: str
    path: Path
    sha256: str
    snapshot: _FileSnapshot


def _snapshot(path: Path) -> _FileSnapshot:
    try:
        details = path.stat()
    except OSError as exc:
        raise ConfigurationError(
            "Dashing input file could not be inspected.",
            code="DASHING_INPUT_UNAVAILABLE",
            context={"path": str(path)},
        ) from exc
    if not stat.S_ISREG(details.st_mode):
        raise ConfigurationError(
            "Dashing path inputs must be regular files.",
            code="INVALID_DASHING_INPUT_PATH",
            context={"path": str(path)},
        )
    return _FileSnapshot(
        device=details.st_dev,
        inode=details.st_ino,
        size=details.st_size,
        modified_ns=details.st_mtime_ns,
    )


def _hash_file(path: Path, *, max_bytes: int) -> str:
    digest = hashlib.sha256()
    observed = 0
    try:
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(min(1_048_576, max_bytes + 1 - observed))
                if not chunk:
                    break
                observed += len(chunk)
                if observed > max_bytes:
                    raise ConfigurationError(
                        "Dashing input exceeds max_input_bytes while being hashed.",
                        code="DASHING_INPUT_LIMIT_EXCEEDED",
                        context={"path": str(path), "max_input_bytes": max_bytes},
                    )
                digest.update(chunk)
    except OSError as exc:
        raise ConfigurationError(
            "Dashing input file could not be read.",
            code="DASHING_INPUT_UNAVAILABLE",
            context={"path": str(path)},
        ) from exc
    return digest.hexdigest()


def _stable_snapshot_and_hash(path: Path, *, max_bytes: int) -> tuple[_FileSnapshot, str]:
    before = _snapshot(path)
    if before.size > max_bytes:
        raise ConfigurationError(
            "Dashing input exceeds max_input_bytes.",
            code="DASHING_INPUT_LIMIT_EXCEEDED",
            context={"path": str(path), "max_input_bytes": max_bytes, "size": before.size},
        )
    digest = _hash_file(path, max_bytes=max_bytes)
    after = _snapshot(path)
    if before != after:
        raise ConfigurationError(
            "Dashing input changed while it was being prepared.",
            code="DASHING_INPUT_CHANGED",
            context={"path": str(path)},
        )
    return after, digest


def _validate_positive_int(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
    return value


def _validate_bool(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be boolean.")
    return value


def _materialize_inputs(
    inputs: DashingInput | DNASet | Iterable[DashingInput], *, max_items: int
) -> tuple[DashingInput, ...]:
    if isinstance(inputs, DNASet):
        resolved: tuple[DashingInput, ...] = inputs.records[: max_items + 1]
    elif isinstance(inputs, DNA):
        resolved = inputs.records[: max_items + 1]
    elif isinstance(inputs, (DNASequence, DNARecord, str, Path)):
        resolved = (inputs,)
    else:
        try:
            resolved = tuple(islice(iter(inputs), max_items + 1))
        except TypeError as exc:
            raise ConfigurationError(
                "Dashing inputs must be sequence objects, a DNASet, or an iterable of paths."
            ) from exc
    if len(resolved) > max_items:
        raise ConfigurationError(
            "Dashing input count exceeds max_items.",
            code="DASHING_ITEM_LIMIT_EXCEEDED",
            context={"max_items": max_items, "item_count_is_lower_bound": True},
        )
    if len(resolved) < 2:
        raise ConfigurationError("Dashing Jaccard requires at least two inputs.")
    if any(not isinstance(item, (DNA, DNASequence, DNARecord, str, Path)) for item in resolved):
        raise ConfigurationError(
            "Every Dashing input must be DNASequence, DNARecord, or a FASTA/FASTQ path."
        )
    return resolved


def _sequence_and_label(item: DNA | DNASequence | DNARecord, index: int) -> tuple[DNASequence, str]:
    if isinstance(item, DNA):
        return item.sequence, item.id
    if isinstance(item, DNARecord):
        return item.sequence, item.id
    return item, f"sequence_{index + 1}"


def _write_sequence_fasta(path: Path, sequence: DNASequence, index: int) -> None:
    symbols = sequence.to_string()
    try:
        with path.open("x", encoding="ascii", newline="\n") as handle:
            handle.write(f">dnakit_input_{index + 1}\n")
            for offset in range(0, len(symbols), 80):
                handle.write(symbols[offset : offset + 80])
                handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    except (OSError, UnicodeEncodeError) as exc:
        raise ConfigurationError(
            "DNAKit could not stage a Dashing FASTA input.",
            code="DASHING_INPUT_STAGING_FAILED",
            context={"input_index": index},
        ) from exc


def _prepare_inputs(
    values: tuple[DashingInput, ...],
    workspace: Path,
    *,
    k: int,
    max_input_bytes: int,
) -> tuple[_PreparedInput, ...]:
    prepared: list[_PreparedInput] = []
    total_bytes = 0
    for index, item in enumerate(values):
        if isinstance(item, (DNA, DNASequence, DNARecord)):
            sequence, label = _sequence_and_label(item, index)
            if sequence.topology is Topology.CIRCULAR:
                raise ConfigurationError(
                    "Dashing adapter requires linear sequence objects.",
                    code="DASHING_CIRCULAR_INPUT_UNSUPPORTED",
                    context={"input_index": index},
                    hint=(
                        "Linearize the sequence explicitly so omitted origin-spanning "
                        "k-mers are auditable."
                    ),
                )
            if sequence.is_gapped:
                raise ConfigurationError(
                    "Dashing adapter does not silently discard sequence gaps.",
                    code="DASHING_GAPPED_INPUT_UNSUPPORTED",
                    context={"input_index": index},
                )
            if sequence.symbol_length < k:
                raise ConfigurationError(
                    "Dashing k cannot exceed an in-memory sequence length.",
                    code="DASHING_K_EXCEEDS_SEQUENCE",
                    context={
                        "input_index": index,
                        "k": k,
                        "sequence_length": sequence.symbol_length,
                    },
                )
            predicted_bytes = sequence.symbol_length + (sequence.symbol_length + 79) // 80 + 32
            if total_bytes + predicted_bytes > max_input_bytes:
                raise ConfigurationError(
                    "Combined Dashing input size exceeds max_input_bytes.",
                    code="DASHING_INPUT_LIMIT_EXCEEDED",
                    context={
                        "max_input_bytes": max_input_bytes,
                        "observed_bytes_is_lower_bound": total_bytes + predicted_bytes,
                    },
                )
            path = workspace / f"input_{index + 1:06d}.fasta"
            _write_sequence_fasta(path, sequence, index)
        else:
            raw_path = Path(item).expanduser()
            if any(character in str(raw_path) for character in ("\x00", "\r", "\n", "\t")):
                raise ConfigurationError(
                    "Dashing input paths cannot contain NUL, tabs, or line breaks.",
                    code="UNSAFE_DASHING_INPUT_PATH",
                    context={"input_index": index},
                )
            try:
                path = raw_path.resolve(strict=True)
            except (OSError, RuntimeError) as exc:
                raise ConfigurationError(
                    "Dashing path input is unavailable.",
                    code="DASHING_INPUT_UNAVAILABLE",
                    context={"input_index": index, "path": str(raw_path)},
                ) from exc
            try:
                str(path).encode("utf-8", errors="strict")
            except UnicodeEncodeError as exc:
                raise ConfigurationError(
                    "Dashing input paths must be valid UTF-8.",
                    code="UNSAFE_DASHING_INPUT_PATH",
                    context={"input_index": index},
                ) from exc
            label = path.name
        snapshot, digest = _stable_snapshot_and_hash(path, max_bytes=max_input_bytes)
        total_bytes += snapshot.size
        if total_bytes > max_input_bytes:
            raise ConfigurationError(
                "Combined Dashing input size exceeds max_input_bytes.",
                code="DASHING_INPUT_LIMIT_EXCEEDED",
                context={
                    "max_input_bytes": max_input_bytes,
                    "observed_bytes": total_bytes,
                    "input_count_is_lower_bound": True,
                },
            )
        prepared.append(_PreparedInput(label, path, digest, snapshot))
    return tuple(prepared)


def _verify_inputs_unchanged(prepared: tuple[_PreparedInput, ...], *, max_bytes: int) -> None:
    for item in prepared:
        current, digest = _stable_snapshot_and_hash(item.path, max_bytes=max_bytes)
        if current != item.snapshot or digest != item.sha256:
            raise BackendExecutionError(
                "A Dashing input changed during external execution.",
                code="DASHING_INPUT_CHANGED",
                context={"path": str(item.path)},
            )


def _write_path_list(path: Path, prepared: tuple[_PreparedInput, ...]) -> None:
    payload = "".join(f"{item.path}\n" for item in prepared)
    try:
        encoded_size = len(payload.encode("utf-8", errors="strict"))
    except UnicodeEncodeError as exc:  # pragma: no cover - rejected while preparing paths.
        raise ConfigurationError("Dashing path list must be valid UTF-8.") from exc
    if encoded_size > 1_000_000:
        raise ConfigurationError(
            "Encoded Dashing path list exceeds its 1 MB safety limit.",
            code="DASHING_PATH_LIST_LIMIT_EXCEEDED",
        )
    try:
        path.write_text(payload, encoding="utf-8", newline="\n")
    except (OSError, UnicodeError) as exc:
        raise ConfigurationError(
            "DNAKit could not stage the Dashing path list.",
            code="DASHING_INPUT_STAGING_FAILED",
        ) from exc


def _command_arguments(
    *,
    k: int,
    threads: int,
    canonical: bool,
    mode: DashingMode,
    sketch_size_log2: int | None,
    matrix_path: Path,
    sizes_path: Path,
    paths_path: Path,
) -> tuple[str, ...]:
    arguments = [
        "dist",
        "--avoid-sorting",
        "-k",
        str(k),
        "-p",
        str(threads),
        "-O",
        str(matrix_path),
        "-o",
        str(sizes_path),
        "-F",
        str(paths_path),
    ]
    if sketch_size_log2 is not None:
        arguments.extend(("-S", str(sketch_size_log2)))
    if not canonical:
        arguments.append("-C")
    if mode == "exact":
        arguments.append("--use-full-khash-sets")
    _assert_command_whitelist(tuple(arguments))
    return tuple(arguments)


def _assert_command_whitelist(arguments: tuple[str, ...]) -> None:
    if not arguments or arguments[0] not in _FIXED_SUBCOMMANDS:
        raise AssertionError("Dashing command escaped the fixed subcommand whitelist.")
    if any(
        token.startswith("-") and token not in _FIXED_FLAGS
        for token in arguments[1:]
        if not token.lstrip("-").isdigit()
    ):
        raise AssertionError("Dashing command escaped the fixed flag whitelist.")


def _parse_matrix(
    path: Path,
    prepared: tuple[_PreparedInput, ...],
    *,
    max_output_bytes: int,
) -> tuple[tuple[float, ...], ...]:
    if path.is_symlink():
        raise BackendExecutionError(
            "Dashing matrix output cannot be a symbolic link.",
            code="UNSAFE_DASHING_OUTPUT_PATH",
        )
    try:
        details = path.stat()
    except OSError as exc:
        raise BackendExecutionError(
            "Dashing did not create its Jaccard matrix output.",
            code="DASHING_OUTPUT_MISSING",
        ) from exc
    if not stat.S_ISREG(details.st_mode):
        raise BackendExecutionError(
            "Dashing matrix output must be a regular file.",
            code="UNSAFE_DASHING_OUTPUT_PATH",
        )
    size = details.st_size
    if not 1 <= size <= max_output_bytes:
        raise BackendExecutionError(
            "Dashing matrix output violates max_output_bytes.",
            code="DASHING_OUTPUT_LIMIT_EXCEEDED",
            context={"observed_bytes": size, "max_output_bytes": max_output_bytes},
        )
    expected_paths = tuple(str(item.path) for item in prepared)
    line_limit = max(1_000_000, max(len(value.encode("utf-8")) for value in expected_paths) + 256)
    try:
        with path.open("r", encoding="utf-8", errors="strict", newline=None) as handle:
            lines = tuple(islice(handle, len(prepared) + 2))
    except (OSError, UnicodeError) as exc:
        raise BackendExecutionError(
            "Dashing matrix output is not readable UTF-8 TSV.",
            code="DASHING_OUTPUT_PROTOCOL_ERROR",
        ) from exc
    if len(lines) != len(prepared) + 1 or any(
        len(line.encode("utf-8")) > line_limit for line in lines
    ):
        raise BackendExecutionError(
            "Dashing matrix output has an invalid row count or line length.",
            code="DASHING_OUTPUT_PROTOCOL_ERROR",
        )
    header = lines[0].rstrip("\r\n").split("\t")
    if tuple(header) != ("##Names", *expected_paths):
        raise BackendExecutionError(
            "Dashing matrix header does not match the staged input order.",
            code="DASHING_OUTPUT_PROTOCOL_ERROR",
        )

    item_count = len(prepared)
    matrix = [
        [1.0 if row == column else 0.0 for column in range(item_count)] for row in range(item_count)
    ]
    for row_index, line in enumerate(lines[1:]):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != item_count + 1 or fields[0] != expected_paths[row_index]:
            raise BackendExecutionError(
                "Dashing matrix row does not match the staged input order.",
                code="DASHING_OUTPUT_PROTOCOL_ERROR",
                context={"row_index": row_index},
            )
        for column_index, token in enumerate(fields[1:]):
            if column_index <= row_index:
                if token != "-":
                    raise BackendExecutionError(
                        "Dashing upper-triangle placeholder is invalid.",
                        code="DASHING_OUTPUT_PROTOCOL_ERROR",
                        context={"row_index": row_index, "column_index": column_index},
                    )
                continue
            try:
                value = float(token)
            except ValueError as exc:
                raise BackendExecutionError(
                    "Dashing matrix contains a non-numeric Jaccard value.",
                    code="DASHING_OUTPUT_PROTOCOL_ERROR",
                    context={"row_index": row_index, "column_index": column_index},
                ) from exc
            if not -_JACCARD_TOLERANCE <= value <= 1 + _JACCARD_TOLERANCE:
                raise BackendExecutionError(
                    "Dashing matrix contains a Jaccard value outside [0, 1].",
                    code="DASHING_OUTPUT_PROTOCOL_ERROR",
                    context={
                        "row_index": row_index,
                        "column_index": column_index,
                        "value": token,
                    },
                )
            bounded = min(1.0, max(0.0, value))
            matrix[row_index][column_index] = bounded
            matrix[column_index][row_index] = bounded
    return tuple(tuple(row) for row in matrix)


def _resolve_output_path(output_path: str | Path | None, *, overwrite: bool) -> Path | None:
    if output_path is None:
        return None
    raw = Path(output_path).expanduser()
    if not raw.name or any(character in raw.name for character in ("\x00", "\r", "\n")):
        raise ConfigurationError("Dashing output_path has an unsafe filename.")
    try:
        parent = raw.parent.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise ConfigurationError("Dashing output_path parent must already exist.") from exc
    if not parent.is_dir():
        raise ConfigurationError("Dashing output_path parent must be a directory.")
    target = parent / raw.name
    if target.is_symlink():
        raise ConfigurationError(
            "Dashing output_path cannot be a symbolic link.", code="UNSAFE_DASHING_OUTPUT_PATH"
        )
    if target.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing Dashing output: {target}")
    if target.exists() and not target.is_file():
        raise ConfigurationError("Dashing output_path must be a regular file path.")
    return target


def _persist_output(source: Path, target: Path, *, overwrite: bool, max_bytes: int) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    written = 0
    try:
        with source.open("rb") as input_handle, os.fdopen(descriptor, "wb") as output_handle:
            while True:
                chunk = input_handle.read(1_048_576)
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    raise BackendExecutionError(
                        "Dashing output exceeded max_output_bytes during persistence.",
                        code="DASHING_OUTPUT_LIMIT_EXCEEDED",
                    )
                output_handle.write(chunk)
            output_handle.flush()
            os.fsync(output_handle.fileno())
        if overwrite:
            os.replace(temporary, target)
        else:
            os.link(temporary, target)
            temporary.unlink()
    except BaseException:
        with suppress(OSError):
            os.close(descriptor)
        temporary.unlink(missing_ok=True)
        raise


def _backend_version(output: str) -> str:
    match = _DASHING_VERSION.search(output)
    if match is None:
        raise BackendExecutionError(
            "Dashing execution did not report a parseable backend version.",
            code="DASHING_VERSION_UNPARSEABLE",
        )
    return match.group("version")


def _provenance(
    executable: Path,
    version: str,
    *,
    elapsed_seconds: float,
    execution_output_sha256: str,
    hybrid: bool,
) -> Provenance:
    backend = BackendInfo(
        "dashing",
        version=version,
        executable_path=str(executable),
        license_expression="GPL-3.0-only",
        capabilities={"sketch", "similarity"},
        available=True,
        metadata={
            "probe_mode": "explicit-scientific-execution",
            "command_profile": _COMMAND_PROFILE,
            "elapsed_seconds": elapsed_seconds,
            "execution_output_sha256": execution_output_sha256,
            "redistributed": False,
            "shell": False,
        },
    )
    return Provenance(
        dependency_versions={"dashing": version},
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.HYBRID if hybrid else ExecutionMode.EXTERNAL,
            origin_class=OriginClass.INTEGRATION,
            license_expression="GPL-3.0-only",
            citations=(
                Citation(
                    "baker-langmead-2019-dashing",
                    title="Dashing: fast and accurate genomic distances with HyperLogLog",
                    doi="10.1186/s13059-019-1875-0",
                ),
            ),
        ),
        backend=backend,
    )


@dataclass(frozen=True, slots=True, init=False)
class DashingAdapter:
    """Opt-in Dashing adapter that never discovers, installs, or bundles Dashing."""

    executable_path: Path

    def __init__(self, executable_path: str | Path) -> None:
        try:
            resolved = Path(executable_path).expanduser().resolve(strict=True)
            details = resolved.stat()
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            raise BackendUnavailableError(
                "A Dashing executable path must be supplied explicitly.",
                context={"executable": str(executable_path)},
            ) from exc
        if not stat.S_ISREG(details.st_mode) or not os.access(resolved, os.X_OK):
            raise BackendUnavailableError(
                "The supplied Dashing path must be an executable regular file.",
                context={"executable": str(resolved)},
            )
        object.__setattr__(self, "executable_path", resolved)

    def matrix(
        self,
        inputs: DashingInput | DNASet | Iterable[DashingInput],
        *,
        k: int = 31,
        mode: DashingMode = "sketch",
        sketch_size_log2: int | None = None,
        canonical: bool = True,
        threads: int = 1,
        temp_dir: str | Path | None = None,
        output_path: str | Path | None = None,
        overwrite: bool = False,
        timeout_seconds: float = 300.0,
        max_items: int = DEFAULT_MAX_DASHING_ITEMS,
        max_input_bytes: int = DEFAULT_MAX_DASHING_INPUT_BYTES,
        max_output_bytes: int = DEFAULT_MAX_DASHING_OUTPUT_BYTES,
        max_capture_bytes: int = DEFAULT_MAX_DASHING_CAPTURE_BYTES,
        max_sketch_memory_bytes: int = DEFAULT_MAX_DASHING_SKETCH_MEMORY_BYTES,
    ) -> DashingJaccardMatrixResult:
        """Compute a validated file-level exact or HLL-sketch Jaccard matrix."""

        resolved_k = _validate_positive_int(k, "k", 32)
        resolved_threads = _validate_positive_int(threads, "threads", 256)
        resolved_max_items = _validate_positive_int(max_items, "max_items", 1_000)
        resolved_max_input = _validate_positive_int(
            max_input_bytes, "max_input_bytes", 10_000_000_000
        )
        resolved_max_output = _validate_positive_int(
            max_output_bytes, "max_output_bytes", 1_000_000_000
        )
        resolved_max_capture = _validate_positive_int(
            max_capture_bytes, "max_capture_bytes", 100_000_000
        )
        resolved_max_sketch_memory = _validate_positive_int(
            max_sketch_memory_bytes, "max_sketch_memory_bytes", 10_000_000_000
        )
        resolved_canonical = _validate_bool(canonical, "canonical")
        resolved_overwrite = _validate_bool(overwrite, "overwrite")
        if mode not in {"sketch", "exact"}:
            raise ConfigurationError("mode must be 'sketch' or 'exact'.")
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or not 0 < timeout_seconds <= 86_400
        ):
            raise ConfigurationError("timeout_seconds must be in (0, 86400].")
        if mode == "sketch":
            resolved_sketch_size = (
                10
                if sketch_size_log2 is None
                else _validate_positive_int(sketch_size_log2, "sketch_size_log2", 30)
            )
        else:
            if sketch_size_log2 is not None:
                raise ConfigurationError("sketch_size_log2 must be None in exact mode.")
            resolved_sketch_size = None

        values = _materialize_inputs(inputs, max_items=resolved_max_items)
        if (
            resolved_sketch_size is not None
            and len(values) * (1 << resolved_sketch_size) > resolved_max_sketch_memory
        ):
            raise ConfigurationError(
                "Requested Dashing sketches exceed max_sketch_memory_bytes.",
                code="DASHING_SKETCH_MEMORY_LIMIT_EXCEEDED",
                context={
                    "item_count": len(values),
                    "sketch_size_log2": resolved_sketch_size,
                    "max_sketch_memory_bytes": resolved_max_sketch_memory,
                },
            )
        output_target = _resolve_output_path(output_path, overwrite=resolved_overwrite)
        if temp_dir is None:
            temp_base: str | None = None
        else:
            try:
                resolved_temp_base = Path(temp_dir).expanduser().resolve(strict=True)
            except (OSError, RuntimeError) as exc:
                raise ConfigurationError("temp_dir must be an existing directory.") from exc
            if not resolved_temp_base.is_dir():
                raise ConfigurationError("temp_dir must be an existing directory.")
            temp_base = str(resolved_temp_base)

        with tempfile.TemporaryDirectory(prefix="dnakit-dashing-", dir=temp_base) as name:
            workspace = Path(name).resolve()
            workspace.chmod(0o700)
            prepared = _prepare_inputs(
                values,
                workspace,
                k=resolved_k,
                max_input_bytes=resolved_max_input,
            )
            if output_target is not None and any(
                output_target.resolve(strict=False) == item.path for item in prepared
            ):
                raise ConfigurationError(
                    "Dashing output_path cannot replace an input file.",
                    code="DASHING_OUTPUT_INPUT_CONFLICT",
                )
            paths_path = workspace / "inputs.paths"
            matrix_path = workspace / "jaccard.tsv"
            sizes_path = workspace / "sizes.tsv"
            _write_path_list(paths_path, prepared)
            arguments = _command_arguments(
                k=resolved_k,
                threads=resolved_threads,
                canonical=resolved_canonical,
                mode=mode,
                sketch_size_log2=resolved_sketch_size,
                matrix_path=matrix_path,
                sizes_path=sizes_path,
                paths_path=paths_path,
            )
            try:
                execution = execute_bounded_command(
                    self.executable_path,
                    arguments,
                    backend_id="dashing",
                    cwd=workspace,
                    timeout_seconds=float(timeout_seconds),
                    max_output_bytes=resolved_max_capture,
                    monitored_output_paths=(matrix_path, sizes_path),
                    max_monitored_output_bytes=resolved_max_output,
                )
            finally:
                _verify_inputs_unchanged(prepared, max_bytes=resolved_max_input)
            if execution.return_code != 0:
                raise BackendExecutionError(
                    "Dashing Jaccard execution returned a non-zero exit status.",
                    code="DASHING_EXECUTION_FAILED",
                    context={
                        "return_code": execution.return_code,
                        "output_excerpt": execution.output[:1_000],
                    },
                )
            version = _backend_version(execution.output)
            matrix = _parse_matrix(
                matrix_path,
                prepared,
                max_output_bytes=resolved_max_output,
            )
            raw_digest = _hash_file(matrix_path, max_bytes=resolved_max_output)
            if not _SHA256.fullmatch(raw_digest):  # pragma: no cover - hashlib guarantee.
                raise AssertionError("hashlib returned an invalid SHA-256 digest")
            if output_target is not None:
                _persist_output(
                    matrix_path,
                    output_target,
                    overwrite=resolved_overwrite,
                    max_bytes=resolved_max_output,
                )
                _, persisted_digest = _stable_snapshot_and_hash(
                    output_target, max_bytes=resolved_max_output
                )
                if persisted_digest != raw_digest:
                    raise BackendExecutionError(
                        "Persisted Dashing output failed its SHA-256 integrity check.",
                        code="DASHING_OUTPUT_INTEGRITY_ERROR",
                    )
            execution_digest = hashlib.sha256(execution.output.encode("utf-8")).hexdigest()
            provenance = _provenance(
                self.executable_path,
                version,
                elapsed_seconds=execution.elapsed_seconds,
                execution_output_sha256=execution_digest,
                hybrid=False,
            )
            parameters = freeze_mapping(
                {
                    "backend": "dashing",
                    "canonical": resolved_canonical,
                    "command_profile": _COMMAND_PROFILE,
                    "deterministic": True,
                    "jaccard_tolerance": _JACCARD_TOLERANCE,
                    "k": resolved_k,
                    "max_capture_bytes": resolved_max_capture,
                    "max_input_bytes": resolved_max_input,
                    "max_items": resolved_max_items,
                    "max_output_bytes": resolved_max_output,
                    "max_sketch_memory_bytes": resolved_max_sketch_memory,
                    "mode": mode,
                    "seed": None,
                    "sketch_size_log2": resolved_sketch_size,
                    "threads": resolved_threads,
                    "timeout_seconds": float(timeout_seconds),
                }
            )
            return DashingJaccardMatrixResult(
                name="dashing_jaccard_matrix",
                method="hll-sketch" if mode == "sketch" else "exact-kmer-set",
                algorithm_version=f"dashing-{version}",
                labels=tuple(item.label for item in prepared),
                values=matrix,
                input_sha256=tuple(item.sha256 for item in prepared),
                parameters=parameters,
                provenance=provenance,
                raw_output_path=None if output_target is None else str(output_target),
                raw_output_sha256=raw_digest,
            )

    def top_k(
        self,
        inputs: DashingInput | DNASet | Iterable[DashingInput],
        *,
        top_k: int,
        k: int = 31,
        mode: DashingMode = "sketch",
        sketch_size_log2: int | None = None,
        canonical: bool = True,
        threads: int = 1,
        temp_dir: str | Path | None = None,
        output_path: str | Path | None = None,
        overwrite: bool = False,
        timeout_seconds: float = 300.0,
        max_items: int = DEFAULT_MAX_DASHING_ITEMS,
        max_input_bytes: int = DEFAULT_MAX_DASHING_INPUT_BYTES,
        max_output_bytes: int = DEFAULT_MAX_DASHING_OUTPUT_BYTES,
        max_capture_bytes: int = DEFAULT_MAX_DASHING_CAPTURE_BYTES,
        max_sketch_memory_bytes: int = DEFAULT_MAX_DASHING_SKETCH_MEMORY_BYTES,
    ) -> DashingTopKResult:
        """Return deterministic per-item Top-k neighbors from one Dashing matrix run."""

        resolved_top_k = _validate_positive_int(top_k, "top_k", 999)
        resolved_max_items = _validate_positive_int(max_items, "max_items", 1_000)
        input_values = _materialize_inputs(inputs, max_items=resolved_max_items)
        if resolved_top_k >= len(input_values):
            raise ConfigurationError("top_k must be smaller than the Dashing input count.")
        matrix = self.matrix(
            input_values,
            k=k,
            mode=mode,
            sketch_size_log2=sketch_size_log2,
            canonical=canonical,
            threads=threads,
            temp_dir=temp_dir,
            output_path=output_path,
            overwrite=overwrite,
            timeout_seconds=timeout_seconds,
            max_items=resolved_max_items,
            max_input_bytes=max_input_bytes,
            max_output_bytes=max_output_bytes,
            max_capture_bytes=max_capture_bytes,
            max_sketch_memory_bytes=max_sketch_memory_bytes,
        )
        rows: list[DashingNeighborRow] = []
        for query_index, row_values in enumerate(matrix.values):
            ranked = sorted(
                (
                    DashingNeighborHit(matrix.labels[index], index, value)
                    for index, value in enumerate(row_values)
                    if index != query_index
                ),
                key=lambda hit: (-hit.jaccard, hit.index),
            )[:resolved_top_k]
            rows.append(
                DashingNeighborRow(
                    query_label=matrix.labels[query_index],
                    query_index=query_index,
                    hits=tuple(ranked),
                )
            )
        backend = matrix.provenance.backend
        if backend is None or backend.version is None:  # pragma: no cover - matrix invariant.
            raise AssertionError("Dashing matrix provenance is missing its backend version")
        elapsed_value = backend.metadata["elapsed_seconds"]
        output_digest_value = backend.metadata["execution_output_sha256"]
        if (
            isinstance(elapsed_value, bool)
            or not isinstance(elapsed_value, (int, float))
            or not isinstance(output_digest_value, str)
        ):  # pragma: no cover - created by _provenance.
            raise AssertionError("Dashing matrix provenance metadata is invalid")
        parameters = freeze_mapping(
            {
                **dict(matrix.parameters),
                "top_k": resolved_top_k,
                "top_k_selection": "dnakit-deterministic-score-index-order",
            }
        )
        provenance = _provenance(
            self.executable_path,
            backend.version,
            elapsed_seconds=float(elapsed_value),
            execution_output_sha256=output_digest_value,
            hybrid=True,
        )
        return DashingTopKResult(
            name="dashing_jaccard_top_k",
            method=matrix.method,
            algorithm_version=matrix.algorithm_version,
            labels=matrix.labels,
            rows=tuple(rows),
            top_k=resolved_top_k,
            input_sha256=matrix.input_sha256,
            parameters=parameters,
            provenance=provenance,
            raw_output_path=matrix.raw_output_path,
            raw_output_sha256=matrix.raw_output_sha256,
        )


__all__ = [
    "DEFAULT_MAX_DASHING_CAPTURE_BYTES",
    "DEFAULT_MAX_DASHING_INPUT_BYTES",
    "DEFAULT_MAX_DASHING_ITEMS",
    "DEFAULT_MAX_DASHING_OUTPUT_BYTES",
    "DEFAULT_MAX_DASHING_SKETCH_MEMORY_BYTES",
    "DashingAdapter",
    "DashingInput",
    "DashingMode",
]
