"""Explicit bounded adapter for caller-installed ``primer3_core``."""

from __future__ import annotations

import math
import os
import re
import stat
import tempfile
from collections.abc import Mapping
from pathlib import Path

from dnakit.backends import execute_bounded_command
from dnakit.core import (
    BackendInfo,
    Citation,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    Issue,
    IssueSeverity,
    OriginClass,
    Provenance,
)
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from ._shared import freeze_parameters, reverse_complement_text, validate_positive_int
from .results import PrimerDesignCandidate, PrimerDesignInterfaceResult

_BOULDER_KEY = re.compile(r"^[A-Z][A-Z0-9_]{0,255}$")
_POSITION_KEY = re.compile(r"^PRIMER_(?:LEFT|RIGHT)_\d+$")
_INTEGER_KEY = re.compile(r"(?:_NUM_RETURNED|_PRODUCT_SIZE)$")
_FLOAT_KEY = re.compile(r"(?:_TM|_GC_PERCENT|_PENALTY)$")


def _explicit_executable(value: str | Path | None) -> tuple[str | None, str | None]:
    if value is None:
        return None, None
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(
            "primer3_core_path must be an explicit filesystem path or None.",
            code="INVALID_PRIMER3_EXECUTABLE",
        )
    requested = str(value)
    if (
        not requested.strip()
        or any(character in requested for character in ("\x00", "\r", "\n"))
        or (not Path(requested).is_absolute() and "/" not in requested and "\\" not in requested)
    ):
        raise ConfigurationError(
            "primer3_core_path must be explicit, not a PATH command name.",
            code="INVALID_PRIMER3_EXECUTABLE",
        )
    try:
        resolved = Path(requested).expanduser().resolve(strict=True)
        details = resolved.stat()
    except (OSError, RuntimeError, ValueError):
        return requested, None
    if not stat.S_ISREG(details.st_mode) or not os.access(resolved, os.X_OK):
        return requested, None
    return requested, str(resolved)


def _explicit_parameter_directory(value: str | Path | None) -> tuple[str | None, str | None]:
    if value is None:
        return None, None
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(
            "thermodynamic_parameters_path must be an explicit directory path or None.",
            code="INVALID_PRIMER3_PARAMETER_PATH",
        )
    requested = str(value)
    try:
        resolved = Path(requested).expanduser().resolve(strict=True)
    except (OSError, RuntimeError, ValueError):
        return requested, None
    return requested, str(resolved) if resolved.is_dir() else None


def _provenance(info: BackendInfo) -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.EXTERNAL,
            origin_class=OriginClass.INTEGRATION,
            license_expression=info.license_expression,
            citations=(
                Citation(
                    "primer3",
                    title="Primer3--new capabilities and interfaces",
                    doi="10.1093/nar/gks596",
                ),
            ),
        ),
        backend=info,
    )


def _bounded_text(value: object, name: str, *, max_length: int) -> str:
    if not isinstance(value, str) or len(value) > max_length:
        raise BackendExecutionError(
            f"Primer3 output {name} must be bounded text.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"field": name},
        )
    return value


def _output_int(value: object, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise BackendExecutionError(
            f"Primer3 output {name} must be an integer >= {minimum}.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"field": name},
        )
    return value


def _output_float(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise BackendExecutionError(
            f"Primer3 output {name} must be finite numeric data.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"field": name},
        )
    return float(value)


def _position_length(value: object, name: str) -> tuple[int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise BackendExecutionError(
            f"Primer3 output {name} must contain position and length.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"field": name},
        )
    return (
        _output_int(value[0], f"{name}[0]"),
        _output_int(value[1], f"{name}[1]", minimum=1),
    )


def _condition_number(conditions: Mapping[str, object], name: str) -> float:
    value = conditions[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"Primer design condition {name} must be numeric.")
    return float(value)


def _boulder_value(key: str, value: object) -> str:
    if key == "SEQUENCE_TARGET":
        if not isinstance(value, list) or len(value) != 2:
            raise AssertionError("SEQUENCE_TARGET must contain start and length.")
        text = f"{value[0]},{value[1]}"
    elif key == "SEQUENCE_EXCLUDED_REGION":
        if not isinstance(value, list):
            raise AssertionError("SEQUENCE_EXCLUDED_REGION must be a list.")
        text = " ".join(f"{item[0]},{item[1]}" for item in value)
    elif key == "PRIMER_PRODUCT_SIZE_RANGE":
        if not isinstance(value, list) or len(value) != 1:
            raise AssertionError("PRIMER_PRODUCT_SIZE_RANGE must contain one interval.")
        text = f"{value[0][0]}-{value[0][1]}"
    elif isinstance(value, bool):
        text = "1" if value else "0"
    elif isinstance(value, (int, float, str)):
        text = str(value)
    else:
        raise ConfigurationError(
            "Primer3 request contains an unsupported Boulder-IO value.",
            code="INVALID_PRIMER3_DESIGN_INPUT",
            context={"key": key, "value_type": type(value).__name__},
        )
    if any(character in text for character in ("\x00", "\r", "\n")):
        raise ConfigurationError(
            "Primer3 Boulder-IO values cannot contain line breaks or NUL.",
            code="INVALID_PRIMER3_DESIGN_INPUT",
            context={"key": key},
        )
    return text


def _boulder_record(sequence_args: Mapping[str, object], global_args: Mapping[str, object]) -> str:
    lines = [
        f"{key}={_boulder_value(key, value)}"
        for arguments in (sequence_args, global_args)
        for key, value in arguments.items()
    ]
    return "\n".join((*lines, "=", ""))


def _parse_boulder_value(key: str, value: str) -> object:
    try:
        if _POSITION_KEY.fullmatch(key):
            parts = value.split(",")
            if len(parts) != 2:
                raise ValueError
            return [int(parts[0]), int(parts[1])]
        if _INTEGER_KEY.search(key):
            return int(value)
        if _FLOAT_KEY.search(key):
            parsed = float(value)
            if not math.isfinite(parsed):
                raise ValueError
            return parsed
    except ValueError as exc:
        raise BackendExecutionError(
            "primer3_core returned malformed numeric Boulder-IO output.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"key": key},
        ) from exc
    return value


def _parse_boulder_output(
    text: str, *, max_result_keys: int, max_text_length: int
) -> dict[str, object]:
    result: dict[str, object] = {}
    terminated = False
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line == "=":
            terminated = True
            if any(trailing.strip() for trailing in lines[index + 1 :]):
                raise BackendExecutionError(
                    "primer3_core returned data after its Boulder-IO record terminator.",
                    code="INVALID_PRIMER3_DESIGN_OUTPUT",
                )
            break
        if "=" not in line:
            raise BackendExecutionError(
                "primer3_core returned malformed Boulder-IO output.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
            )
        key, value = line.split("=", 1)
        if not _BOULDER_KEY.fullmatch(key) or key in result:
            raise BackendExecutionError(
                "primer3_core returned an invalid or duplicate Boulder-IO key.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
                context={"key": key[:256]},
            )
        if len(value) > max_text_length:
            raise BackendExecutionError(
                "primer3_core returned an oversized Boulder-IO value.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
                context={"key": key},
            )
        result[key] = _parse_boulder_value(key, value)
        if len(result) > max_result_keys:
            raise BackendExecutionError(
                "primer3_core output exceeded max_result_keys.",
                code="PRIMER3_DESIGN_OUTPUT_LIMIT_EXCEEDED",
            )
    if not terminated:
        raise BackendExecutionError(
            "primer3_core output did not terminate its Boulder-IO record.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
        )
    return result


def _read_bounded(path: Path, max_bytes: int) -> str:
    try:
        if not path.is_file():
            raise OSError
        with path.open("rb") as handle:
            payload = handle.read(max_bytes + 1)
        if len(payload) > max_bytes:
            raise OSError
        return payload.decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise BackendExecutionError(
            "primer3_core output artifact is missing, oversized, or invalid UTF-8.",
            code="INVALID_PRIMER3_DESIGN_OUTPUT",
            context={"artifact": path.name, "max_bytes": max_bytes},
        ) from exc


class Primer3CLIDesignAdapter:
    """Run primer design through an explicit user-installed ``primer3_core`` path."""

    def __init__(
        self,
        primer3_core_path: str | Path | None = None,
        *,
        thermodynamic_parameters_path: str | Path | None = None,
        timeout_seconds: float = 60.0,
        max_output_bytes: int = 10_000_000,
    ) -> None:
        self._requested_executable, self._path = _explicit_executable(primer3_core_path)
        self._requested_parameters, self._parameters_path = _explicit_parameter_directory(
            thermodynamic_parameters_path
        )
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or not 0 < timeout_seconds <= 300
        ):
            raise ConfigurationError("timeout_seconds must be in (0, 300].")
        if (
            isinstance(max_output_bytes, bool)
            or not isinstance(max_output_bytes, int)
            or not 1 <= max_output_bytes <= 100_000_000
        ):
            raise ConfigurationError("max_output_bytes must be in [1, 100000000].")
        self._timeout_seconds = float(timeout_seconds)
        self._max_output_bytes = max_output_bytes

    @property
    def info(self) -> BackendInfo:
        parameters_valid = self._requested_parameters is None or self._parameters_path is not None
        available = self._path is not None and parameters_valid
        return BackendInfo(
            "primer3-cli",
            executable_path=self._path,
            license_expression="GPL-2.0-or-later",
            capabilities=("primer_design",) if available else (),
            available=available,
            metadata={
                "adapter_status": "explicit-cli-execution",
                "requested_primer3_core_path": self._requested_executable,
                "resolved_primer3_core_path": self._path,
                "requested_thermodynamic_parameters_path": self._requested_parameters,
                "resolved_thermodynamic_parameters_path": self._parameters_path,
                "probe_mode": "explicit-path-stat-only",
                "automatic_path_search": False,
                "automatic_install": False,
                "automatic_download": False,
                "redistributed": False,
            },
        )

    def ensure_available(self) -> None:
        if self._requested_parameters is not None and self._parameters_path is None:
            raise BackendUnavailableError(
                "The explicit Primer3 thermodynamic parameter directory is unavailable.",
                code="PRIMER3_PARAMETER_PATH_UNAVAILABLE",
                context={"requested_path": self._requested_parameters},
            )
        if self._path is None:
            raise BackendUnavailableError(
                "The explicit primer3_core executable is unavailable.",
                code="PRIMER3_CLI_UNAVAILABLE",
                context={"requested_executable": self._requested_executable},
                hint="Install Primer3 separately and pass primer3_core_path explicitly.",
            )

    @staticmethod
    def _arguments(
        result: PrimerDesignInterfaceResult,
    ) -> tuple[dict[str, object], dict[str, object]]:
        request = result.request
        conditions = request.thermodynamic_conditions
        primer_min, primer_max = request.primer_length_range
        tm_min, tm_max = request.tm_range_celsius
        gc_min, gc_max = request.gc_range
        product_min, product_max = request.product_length_range
        sequence_args: dict[str, object] = {
            "SEQUENCE_ID": "dnakit-template",
            "SEQUENCE_TEMPLATE": request.template.symbols,
            "SEQUENCE_TARGET": [request.target_start, request.target_end - request.target_start],
            "SEQUENCE_EXCLUDED_REGION": [
                [start, end - start] for start, end in request.excluded_regions
            ],
        }
        global_args: dict[str, object] = {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_NUM_RETURN": request.candidate_count,
            "PRIMER_OPT_SIZE": (primer_min + primer_max) // 2,
            "PRIMER_MIN_SIZE": primer_min,
            "PRIMER_MAX_SIZE": primer_max,
            "PRIMER_OPT_TM": (tm_min + tm_max) / 2.0,
            "PRIMER_MIN_TM": tm_min,
            "PRIMER_MAX_TM": tm_max,
            "PRIMER_MIN_GC": gc_min * 100.0,
            "PRIMER_MAX_GC": gc_max * 100.0,
            "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
            "PRIMER_SALT_MONOVALENT": (
                _condition_number(conditions, "sodium_molar")
                + _condition_number(conditions, "potassium_molar")
            )
            * 1_000.0,
            "PRIMER_SALT_DIVALENT": _condition_number(conditions, "magnesium_molar") * 1_000.0,
            "PRIMER_DNTP_CONC": _condition_number(conditions, "dntp_molar") * 1_000.0,
            "PRIMER_DNA_CONC": _condition_number(conditions, "strand_concentration_molar")
            * 1_000_000_000.0,
            "PRIMER_DMSO_CONC": _condition_number(conditions, "dmso_percent"),
            "PRIMER_DMSO_FACTOR": _condition_number(conditions, "dmso_factor_celsius_per_percent"),
            "PRIMER_FORMAMIDE_CONC": _condition_number(conditions, "formamide_molar"),
        }
        return sequence_args, global_args

    @staticmethod
    def _candidate(
        raw: Mapping[str, object],
        *,
        rank: int,
        template: str,
        max_text_length: int,
    ) -> PrimerDesignCandidate:
        template_length = len(template)
        prefix = f"PRIMER_PAIR_{rank}"
        left_position, left_length = _position_length(
            raw.get(f"PRIMER_LEFT_{rank}"), f"left_{rank}"
        )
        right_anchor, right_length = _position_length(
            raw.get(f"PRIMER_RIGHT_{rank}"), f"right_{rank}"
        )
        left_start = left_position
        left_end = left_start + left_length
        right_start = right_anchor - right_length + 1
        right_end = right_anchor + 1
        if not 0 <= left_start < left_end <= right_start < right_end <= template_length:
            raise BackendExecutionError(
                "Primer3 returned candidate coordinates outside the template.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
                context={"candidate_rank": rank},
            )
        left_sequence = _bounded_text(
            raw.get(f"PRIMER_LEFT_{rank}_SEQUENCE"),
            f"left_sequence_{rank}",
            max_length=max_text_length,
        ).upper()
        right_sequence = _bounded_text(
            raw.get(f"PRIMER_RIGHT_{rank}_SEQUENCE"),
            f"right_sequence_{rank}",
            max_length=max_text_length,
        ).upper()
        expected_left = template[left_start:left_end]
        expected_right = reverse_complement_text(template[right_start:right_end])
        if (
            len(left_sequence) != left_length
            or len(right_sequence) != right_length
            or left_sequence != expected_left
            or right_sequence != expected_right
        ):
            raise BackendExecutionError(
                "Primer3 returned candidate sequences inconsistent with template coordinates.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
                context={"candidate_rank": rank},
            )
        penalty_value = raw.get(f"{prefix}_PENALTY")
        try:
            product_size = _output_int(raw.get(f"{prefix}_PRODUCT_SIZE"), "product_size", minimum=1)
            if product_size != right_end - left_start:
                raise ConfigurationError(
                    "Primer3 product size does not match candidate coordinates."
                )
            return PrimerDesignCandidate(
                rank=rank,
                left_primer_sequence=left_sequence,
                right_primer_sequence=right_sequence,
                left_start=left_start,
                left_end=left_end,
                right_start=right_start,
                right_end=right_end,
                product_size=product_size,
                left_tm_celsius=_output_float(raw.get(f"PRIMER_LEFT_{rank}_TM"), f"left_tm_{rank}"),
                right_tm_celsius=_output_float(
                    raw.get(f"PRIMER_RIGHT_{rank}_TM"), f"right_tm_{rank}"
                ),
                left_gc_fraction=_output_float(
                    raw.get(f"PRIMER_LEFT_{rank}_GC_PERCENT"), f"left_gc_{rank}"
                )
                / 100.0,
                right_gc_fraction=_output_float(
                    raw.get(f"PRIMER_RIGHT_{rank}_GC_PERCENT"), f"right_gc_{rank}"
                )
                / 100.0,
                pair_penalty=(
                    None
                    if penalty_value is None
                    else _output_float(penalty_value, f"pair_penalty_{rank}")
                ),
            )
        except ConfigurationError as exc:
            raise BackendExecutionError(
                "Primer3 returned an internally inconsistent candidate.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
                context={"candidate_rank": rank},
            ) from exc

    def _execute(
        self,
        record: str,
        *,
        max_result_keys: int,
        max_text_length: int,
    ) -> dict[str, object]:
        self.ensure_available()
        assert self._path is not None
        with tempfile.TemporaryDirectory(prefix="dnakit-primer3-design-") as workspace_text:
            workspace = Path(workspace_text)
            input_path = workspace / "input.boulder"
            output_path = workspace / "output.boulder"
            error_path = workspace / "error.txt"
            try:
                with input_path.open("x", encoding="utf-8", newline="\n") as handle:
                    handle.write(record)
                    handle.flush()
                    os.fsync(handle.fileno())
            except OSError as exc:
                raise BackendExecutionError(
                    "DNAKit could not stage the primer3_core input.",
                    code="PRIMER3_DESIGN_INPUT_STAGING_FAILED",
                ) from exc
            arguments = (
                "--strict_tags",
                "--io_version=4",
                "--output=output.boulder",
                "--error=error.txt",
                "input.boulder",
            )
            result = execute_bounded_command(
                self._path,
                arguments,
                backend_id="primer3-cli",
                cwd=workspace,
                timeout_seconds=self._timeout_seconds,
                max_output_bytes=min(self._max_output_bytes, 1_000_000),
                monitored_output_paths=(output_path, error_path),
                max_monitored_output_bytes=self._max_output_bytes,
            )
            error_text = ""
            if error_path.exists():
                error_text = _read_bounded(error_path, self._max_output_bytes)
            if result.return_code != 0:
                raise BackendExecutionError(
                    "primer3_core primer design failed.",
                    code="PRIMER3_CLI_EXECUTION_FAILED",
                    context={
                        "return_code": result.return_code,
                        "process_output_excerpt": result.output[:1_000],
                        "error_excerpt": error_text[:1_000],
                    },
                )
            output_text = _read_bounded(output_path, self._max_output_bytes)
        return _parse_boulder_output(
            output_text,
            max_result_keys=max_result_keys,
            max_text_length=max_text_length,
        )

    def design(
        self,
        prepared: PrimerDesignInterfaceResult,
        *,
        max_returned_candidates: int = 10_000,
        max_result_keys: int = 100_000,
        max_text_length: int = 100_000,
        max_template_length: int = 1_000_000,
    ) -> PrimerDesignInterfaceResult:
        """Execute one prepared request; no discovery or execution happens earlier."""

        if not isinstance(prepared, PrimerDesignInterfaceResult):
            raise ConfigurationError("prepared must be a PrimerDesignInterfaceResult.")
        if prepared.execution_performed:
            raise ConfigurationError("A completed primer-design result cannot be executed again.")
        validate_positive_int(max_returned_candidates, "max_returned_candidates", maximum=10_000)
        validate_positive_int(max_result_keys, "max_result_keys", maximum=1_000_000)
        validate_positive_int(max_text_length, "max_text_length", maximum=1_000_000)
        validate_positive_int(
            max_template_length,
            "max_template_length",
            maximum=10_000_000,
        )
        template_length = prepared.request.template.symbol_length
        if template_length > max_template_length:
            raise ConfigurationError(
                "Primer design template exceeds max_template_length.",
                code="PRIMER3_DESIGN_TEMPLATE_LIMIT_EXCEEDED",
                context={
                    "template_length": template_length,
                    "max_template_length": max_template_length,
                },
            )
        output_limit = min(prepared.request.candidate_count, max_returned_candidates)
        sequence_args, global_args = self._arguments(prepared)
        global_args["PRIMER_NUM_RETURN"] = output_limit
        if self._parameters_path is not None:
            global_args["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = self._parameters_path + os.sep
        raw = self._execute(
            _boulder_record(sequence_args, global_args),
            max_result_keys=max_result_keys,
            max_text_length=max_text_length,
        )
        if len(raw) > max_result_keys or any(
            not isinstance(key, str) or not key or len(key) > 256 for key in raw
        ):
            raise BackendExecutionError(
                "primer3_core returned an invalid or oversized result mapping.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
            )
        if any(isinstance(value, str) and len(value) > max_text_length for value in raw.values()):
            raise BackendExecutionError(
                "primer3_core returned oversized text output.",
                code="INVALID_PRIMER3_DESIGN_OUTPUT",
            )
        backend_error = _bounded_text(
            raw.get("PRIMER_ERROR", ""), "PRIMER_ERROR", max_length=max_text_length
        )
        if backend_error:
            raise BackendExecutionError(
                "primer3_core rejected the primer-design request.",
                code="PRIMER3_CLI_DESIGN_REJECTED",
                context={"backend_error": backend_error},
            )
        returned = _output_int(raw.get("PRIMER_PAIR_NUM_RETURNED", 0), "pair_count")
        if returned > output_limit:
            raise BackendExecutionError(
                "primer3_core returned more candidates than the configured limit.",
                code="PRIMER3_DESIGN_OUTPUT_LIMIT_EXCEEDED",
                context={"returned": returned, "output_limit": output_limit},
            )
        candidates = tuple(
            self._candidate(
                raw,
                rank=rank,
                template=prepared.request.template.symbols,
                max_text_length=max_text_length,
            )
            for rank in range(returned)
        )
        warning = _bounded_text(
            raw.get("PRIMER_WARNING", ""), "PRIMER_WARNING", max_length=max_text_length
        )
        issues = (
            (Issue("PRIMER3_DESIGN_WARNING", IssueSeverity.WARNING, warning),) if warning else ()
        )
        info = self.info
        return PrimerDesignInterfaceResult(
            request=prepared.request,
            backend_name=info.name,
            status=(
                "execution-complete-with-candidates"
                if candidates
                else "execution-complete-no-candidates"
            ),
            execution_performed=True,
            candidates=candidates,
            reason="The user explicitly invoked a separately installed primer3_core executable.",
            method="primer3-cli-core-design-adapter",
            algorithm_version="dnakit-primer3-cli-design-adapter-v1",
            parameters=freeze_parameters(
                {
                    "explicit_execution": True,
                    "automatic_path_search": False,
                    "automatic_install": False,
                    "automatic_probe": False,
                    "automatic_execution": False,
                    "requested_candidates": prepared.request.candidate_count,
                    "returned_candidates": len(candidates),
                    "max_returned_candidates": max_returned_candidates,
                    "max_result_keys": max_result_keys,
                    "max_text_length": max_text_length,
                    "max_template_length": max_template_length,
                    "primer3_core_path": self._path,
                    "thermodynamic_parameters_path": self._parameters_path,
                }
            ),
            provenance=_provenance(info),
            issues=issues,
        )


__all__ = ["Primer3CLIDesignAdapter"]
