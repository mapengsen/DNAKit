"""Explicit bounded adapters for caller-installed Primer3 CLI programs."""

from __future__ import annotations

import math
import os
import re
import stat
import tempfile
from pathlib import Path
from typing import Protocol, runtime_checkable

from dnakit.backends import execute_bounded_command
from dnakit.core import (
    BackendInfo,
    Citation,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    Issue,
    IssueSeverity,
    OriginClass,
    Provenance,
)
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from ._shared import canonical_linear_symbols
from .config import ThermodynamicConditions
from .results import ConditionalCapability, Primer3ThermodynamicResult

_PRIMER3_CAPABILITIES = frozenset({"hairpin", "self_dimer", "heterodimer", "tm"})
_PRIMER3_STRUCTURE_CAPABILITIES = frozenset({"hairpin", "self_dimer", "heterodimer"})
_THERMODYNAMIC_LINE = re.compile(
    r"^Calculated thermodynamical parameters for dimer:\s*"
    r"(?:\d+\s+)?dS\s*=\s*(?P<ds>[-+0-9.eE]+)\s+"
    r"dH\s*=\s*(?P<dh>[-+0-9.eE]+)\s+"
    r"dG\s*=\s*(?P<dg>[-+0-9.eE]+)\s+"
    r"t\s*=\s*(?P<tm>[-+0-9.eE]+)\s*$"
)
_NO_STRUCTURE_LINE = re.compile(r"^0\s+dS\s*=.*(?:inf|nan)", re.IGNORECASE)


def _primer3_number(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{field} must be finite numeric output")
    return float(value)


def validate_primer3_result(
    result: object,
    *,
    capability: str,
    sequences_5to3: tuple[str, ...],
    conditions: ThermodynamicConditions,
    max_loop: int | None,
    output_structure: bool,
    error_code: str = "INVALID_PRIMER3_RESULT",
) -> Primer3ThermodynamicResult:
    """Validate that an adapter result is bound to the exact resolved request."""

    if not isinstance(result, Primer3ThermodynamicResult):
        raise BackendExecutionError(
            "Primer3 adapter returned an unsupported result object.",
            code=error_code,
            context={"result_type": type(result).__name__, "capability": capability},
        )
    mismatches: list[str] = []
    if result.capability != capability:
        mismatches.append("capability")
    if result.sequences_5to3 != sequences_5to3:
        mismatches.append("sequences_5to3")
    if result.conditions != conditions:
        mismatches.append("conditions")
    if result.max_loop != max_loop:
        mismatches.append("max_loop")
    if not output_structure and result.ascii_structure is not None:
        mismatches.append("ascii_structure")
    if result.backend != result.provenance.backend:
        mismatches.append("provenance.backend")
    if result.provenance.implementation.label is not ImplementationLabel.ADAPTER:
        mismatches.append("provenance.implementation")
    if not result.backend.available or capability not in result.backend.capabilities:
        mismatches.append("backend.capabilities")
    if mismatches:
        raise BackendExecutionError(
            "Primer3 adapter result does not match the resolved request contract.",
            code=error_code,
            context={"capability": capability, "mismatched_fields": tuple(mismatches)},
        )
    return result


def _explicit_executable(value: str | Path | None, name: str) -> tuple[str | None, str | None]:
    if value is None:
        return None, None
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(
            f"{name} must be an explicit filesystem path or None.",
            code="INVALID_PRIMER3_EXECUTABLE",
        )
    requested = str(value)
    if (
        not requested.strip()
        or any(character in requested for character in ("\x00", "\r", "\n"))
        or (not Path(requested).is_absolute() and "/" not in requested and "\\" not in requested)
    ):
        raise ConfigurationError(
            f"{name} must be an explicit filesystem path, not a PATH command name.",
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


def _explicit_parameter_directory(
    value: str | Path | None,
) -> tuple[str | None, str | None]:
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


def _conditions(value: ThermodynamicConditions | None) -> ThermodynamicConditions:
    resolved = ThermodynamicConditions() if value is None else value
    if not isinstance(resolved, ThermodynamicConditions):
        raise ConfigurationError("conditions must be ThermodynamicConditions or None.")
    return resolved


def _max_loop(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 30:
        raise ConfigurationError(
            "max_loop must be an integer in [1, 30].",
            code="INVALID_PRIMER3_MAX_LOOP",
        )
    return value


def _number_argument(value: float) -> str:
    return format(value, ".15g")


def _adapter_provenance(info: BackendInfo) -> Provenance:
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


@runtime_checkable
class ThermodynamicsBackend(Protocol):
    """Minimal backend-neutral interface exposed to higher-level workflows."""

    @property
    def info(self) -> BackendInfo: ...

    def supports(self, capability: str) -> bool: ...

    def ensure_available(self, capability: str) -> None: ...


class Primer3CLIAdapter:
    """Call only explicit user-installed ``oligotm`` and ``ntthal`` paths.

    Construction and :attr:`info` are passive. No PATH search, import,
    download, installation, or subprocess execution occurs until a calculation
    method is explicitly called.
    """

    def __init__(
        self,
        *,
        oligotm_path: str | Path | None = None,
        ntthal_path: str | Path | None = None,
        thermodynamic_parameters_path: str | Path | None = None,
        timeout_seconds: float = 30.0,
        max_output_bytes: int = 1_000_000,
    ) -> None:
        self._requested_oligotm, self._oligotm_path = _explicit_executable(
            oligotm_path, "oligotm_path"
        )
        self._requested_ntthal, self._ntthal_path = _explicit_executable(ntthal_path, "ntthal_path")
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
            or not 1 <= max_output_bytes <= 10_000_000
        ):
            raise ConfigurationError("max_output_bytes must be in [1, 10000000].")
        self._timeout_seconds = float(timeout_seconds)
        self._max_output_bytes = max_output_bytes

    @property
    def info(self) -> BackendInfo:
        paths_valid = self._requested_parameters is None or self._parameters_path is not None
        capabilities: set[str] = set()
        if self._oligotm_path is not None:
            capabilities.add("tm")
        if paths_valid and self._ntthal_path is not None:
            capabilities.update(_PRIMER3_STRUCTURE_CAPABILITIES)
        executable_path = (
            self._ntthal_path
            if paths_valid and self._ntthal_path is not None
            else self._oligotm_path
        )
        return BackendInfo(
            "primer3-cli",
            version=None,
            executable_path=executable_path,
            license_expression="GPL-2.0-or-later",
            capabilities=tuple(sorted(capabilities)),
            available=bool(capabilities),
            metadata={
                "adapter_status": "explicit-cli-execution",
                "requested_oligotm_path": self._requested_oligotm,
                "resolved_oligotm_path": self._oligotm_path,
                "requested_ntthal_path": self._requested_ntthal,
                "resolved_ntthal_path": self._ntthal_path,
                "requested_thermodynamic_parameters_path": self._requested_parameters,
                "resolved_thermodynamic_parameters_path": self._parameters_path,
                "probe_mode": "explicit-path-stat-only",
                "version_probe_executed": False,
                "automatic_path_search": False,
                "automatic_install": False,
                "automatic_download": False,
                "redistributed": False,
            },
        )

    def supports(self, capability: str) -> bool:
        return capability in _PRIMER3_CAPABILITIES

    def ensure_available(self, capability: str) -> None:
        if not self.supports(capability):
            raise BackendExecutionError(
                "Primer3 does not declare the requested thermodynamic capability.",
                code="BACKEND_CAPABILITY_UNSUPPORTED",
                context={"backend": "primer3-cli", "capability": capability},
            )
        if (
            capability in _PRIMER3_STRUCTURE_CAPABILITIES
            and self._requested_parameters is not None
            and self._parameters_path is None
        ):
            raise BackendUnavailableError(
                "The explicit Primer3 thermodynamic parameter directory is unavailable.",
                code="PRIMER3_PARAMETER_PATH_UNAVAILABLE",
                context={"requested_path": self._requested_parameters},
            )
        path = self._oligotm_path if capability == "tm" else self._ntthal_path
        requested = self._requested_oligotm if capability == "tm" else self._requested_ntthal
        if path is None:
            raise BackendUnavailableError(
                "The required explicit Primer3 CLI executable is unavailable.",
                code="PRIMER3_CLI_UNAVAILABLE",
                context={"capability": capability, "requested_executable": requested},
                hint="Install Primer3 separately and pass the executable's explicit path.",
            )

    def _run(self, executable: str, arguments: tuple[str, ...], capability: str) -> str:
        with tempfile.TemporaryDirectory(prefix="dnakit-primer3-") as workspace:
            result = execute_bounded_command(
                executable,
                arguments,
                backend_id="primer3-cli",
                cwd=workspace,
                timeout_seconds=self._timeout_seconds,
                max_output_bytes=self._max_output_bytes,
            )
        if result.return_code != 0:
            raise BackendExecutionError(
                "Primer3 CLI calculation failed.",
                code="PRIMER3_CLI_EXECUTION_FAILED",
                context={
                    "capability": capability,
                    "return_code": result.return_code,
                    "output_excerpt": result.output[:1_000],
                },
            )
        return result.output

    @staticmethod
    def _oligotm_arguments(symbols: str, conditions: ThermodynamicConditions) -> tuple[str, ...]:
        # oligotm currently caps argc at 14, so values equal to its documented
        # defaults are omitted and only effective overrides are passed.
        condition_options: list[tuple[str, float]] = []
        candidates = (
            ("-mv", conditions.monovalent_molar * 1_000.0, 50.0),
            ("-dv", conditions.magnesium_molar * 1_000.0, 1.5),
            ("-n", conditions.dntp_molar * 1_000.0, 0.6),
            ("-d", conditions.strand_concentration_molar * 1_000_000_000.0, 50.0),
            ("-dm", conditions.dmso_percent, 0.0),
            ("-df", conditions.dmso_factor_celsius_per_percent, 0.6),
            ("-fo", conditions.formamide_molar, 0.0),
        )
        for flag, value, default in candidates:
            if value != default:
                condition_options.append((flag, value))
        if len(condition_options) > 4:
            raise ConfigurationError(
                "The installed oligotm CLI cannot represent all requested condition overrides.",
                code="PRIMER3_OLIGOTM_ARGUMENT_LIMIT",
                context={"override_count": len(condition_options), "maximum": 4},
            )
        # Pin both scientific method selectors instead of relying on the
        # installed executable's defaults. oligotm caps argc at 14, leaving
        # room for four non-default condition flags plus the sequence.
        options = [("-tp", 1.0), ("-sc", 1.0), *condition_options]
        arguments = tuple(
            item for flag, value in options for item in (flag, _number_argument(value))
        )
        return (*arguments, symbols)

    def tm(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
    ) -> Primer3ThermodynamicResult:
        self.ensure_available("tm")
        symbols = canonical_linear_symbols(
            sequence, operation="Primer3 oligotm", min_length=2, max_length=36
        )
        resolved = _conditions(conditions)
        assert self._oligotm_path is not None
        output = self._run(
            self._oligotm_path,
            self._oligotm_arguments(symbols, resolved),
            "tm",
        )
        try:
            tm_celsius = float(output.strip())
            if not math.isfinite(tm_celsius):
                raise ValueError
        except ValueError as exc:
            raise BackendExecutionError(
                "oligotm returned invalid Tm output.",
                code="INVALID_PRIMER3_CLI_OUTPUT",
                context={"capability": "tm", "output_excerpt": output[:1_000]},
            ) from exc
        info = self.info
        return Primer3ThermodynamicResult(
            capability="tm",
            sequences_5to3=(symbols,),
            structure_found=None,
            tm_celsius=tm_celsius,
            delta_g_kcal_per_mol=None,
            delta_h_kcal_per_mol=None,
            delta_s_cal_per_k_mol=None,
            ascii_structure=None,
            conditions=resolved,
            max_loop=None,
            method="primer3-cli-oligotm-santalucia",
            algorithm_version="dnakit-primer3-cli-adapter-v1",
            backend=info,
            provenance=_adapter_provenance(info),
            issues=(),
        )

    def _structure(
        self,
        capability: str,
        sequence_a: DNASequence,
        sequence_b: DNASequence | None,
        *,
        conditions: ThermodynamicConditions | None,
        max_loop: int,
        output_structure: bool,
    ) -> Primer3ThermodynamicResult:
        self.ensure_available(capability)
        if not isinstance(output_structure, bool):
            raise ConfigurationError("output_structure must be boolean.")
        symbols_a = canonical_linear_symbols(
            sequence_a, operation=f"Primer3 {capability}", min_length=1, max_length=60
        )
        symbols_b = (
            None
            if sequence_b is None
            else canonical_linear_symbols(
                sequence_b, operation=f"Primer3 {capability}", min_length=1, max_length=60
            )
        )
        resolved = _conditions(conditions)
        if resolved.dmso_percent != 0.0 or resolved.formamide_molar != 0.0:
            raise ConfigurationError(
                "ntthal does not expose DMSO or formamide condition flags.",
                code="PRIMER3_NTTHAL_COSOLVENT_UNSUPPORTED",
            )
        loop = _max_loop(max_loop)
        mode = "HAIRPIN" if capability == "hairpin" else "ANY"
        arguments = [
            "-a",
            mode,
            "-mv",
            _number_argument(resolved.monovalent_molar * 1_000.0),
            "-dv",
            _number_argument(resolved.magnesium_molar * 1_000.0),
            "-n",
            _number_argument(resolved.dntp_molar * 1_000.0),
            "-d",
            _number_argument(resolved.strand_concentration_molar * 1_000_000_000.0),
            "-t",
            _number_argument(resolved.temperature_celsius),
            "-maxloop",
            str(loop),
        ]
        if self._parameters_path is not None:
            arguments.extend(("-path", self._parameters_path))
        arguments.extend(("-s1", symbols_a))
        if capability == "self_dimer":
            arguments.extend(("-s2", symbols_a))
        elif capability == "heterodimer":
            if symbols_b is None:
                raise AssertionError("Heterodimer calculations require a second sequence.")
            arguments.extend(("-s2", symbols_b))
        assert self._ntthal_path is not None
        output = self._run(self._ntthal_path, tuple(arguments), capability)
        lines = output.splitlines()
        match = _THERMODYNAMIC_LINE.fullmatch(lines[0].strip()) if lines else None
        issues: tuple[Issue, ...] = ()
        if match is None and lines and _NO_STRUCTURE_LINE.match(lines[0].strip()):
            structure_found = False
            tm_celsius = delta_g = delta_h = delta_s = 0.0
            ascii_structure = None
            issues = (
                Issue(
                    "PRIMER3_NO_STRUCTURE",
                    IssueSeverity.INFO,
                    "ntthal reported no predicted secondary structure.",
                ),
            )
        elif match is not None:
            try:
                tm_celsius = float(match.group("tm"))
                delta_g = float(match.group("dg")) / 1_000.0
                delta_h = float(match.group("dh")) / 1_000.0
                delta_s = float(match.group("ds"))
                if not all(math.isfinite(item) for item in (tm_celsius, delta_g, delta_h, delta_s)):
                    raise ValueError
            except ValueError as exc:
                raise BackendExecutionError(
                    "ntthal returned non-finite thermodynamic output.",
                    code="INVALID_PRIMER3_CLI_OUTPUT",
                    context={"capability": capability},
                ) from exc
            structure_found = True
            ascii_structure = "\n".join(lines[1:]) if output_structure and len(lines) > 1 else None
        else:
            raise BackendExecutionError(
                "ntthal returned unrecognized thermodynamic output.",
                code="INVALID_PRIMER3_CLI_OUTPUT",
                context={"capability": capability, "output_excerpt": output[:1_000]},
            )
        info = self.info
        return Primer3ThermodynamicResult(
            capability=capability,
            sequences_5to3=(symbols_a,) if symbols_b is None else (symbols_a, symbols_b),
            structure_found=structure_found,
            tm_celsius=tm_celsius,
            delta_g_kcal_per_mol=delta_g,
            delta_h_kcal_per_mol=delta_h,
            delta_s_cal_per_k_mol=delta_s,
            ascii_structure=ascii_structure,
            conditions=resolved,
            max_loop=loop,
            method=f"primer3-cli-ntthal-{capability.replace('_', '-')}",
            algorithm_version="dnakit-primer3-cli-adapter-v1",
            backend=info,
            provenance=_adapter_provenance(info),
            issues=issues,
        )

    def hairpin(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._structure(
            "hairpin",
            sequence,
            None,
            conditions=conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )

    def self_dimer(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._structure(
            "self_dimer",
            sequence,
            None,
            conditions=conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )

    def heterodimer(
        self,
        sequence_a: DNASequence,
        sequence_b: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult:
        return self._structure(
            "heterodimer",
            sequence_a,
            sequence_b,
            conditions=conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )


def probe_primer3(
    *,
    oligotm_path: str | Path | None = None,
    ntthal_path: str | Path | None = None,
    thermodynamic_parameters_path: str | Path | None = None,
) -> BackendInfo:
    """Inspect explicit Primer3 paths without searching, importing, or executing."""

    return Primer3CLIAdapter(
        oligotm_path=oligotm_path,
        ntthal_path=ntthal_path,
        thermodynamic_parameters_path=thermodynamic_parameters_path,
    ).info


def conditional_capabilities(
    *,
    primer3_ntthal_path: str | Path | None = None,
    primer3_thermodynamic_parameters_path: str | Path | None = None,
) -> tuple[ConditionalCapability, ...]:
    """Describe conditional structure capabilities without automatic discovery."""

    primer3 = probe_primer3(
        ntthal_path=primer3_ntthal_path,
        thermodynamic_parameters_path=primer3_thermodynamic_parameters_path,
    )
    from dnakit.secondary_structure import probe_nupack

    nupack = probe_nupack()
    execution_supported = all(
        capability in primer3.capabilities for capability in _PRIMER3_STRUCTURE_CAPABILITIES
    )
    status = "available-user-supplied-cli" if execution_supported else "conditional-unavailable"
    return (
        ConditionalCapability(
            "THERMO-008",
            "hairpin",
            status,
            execution_supported,
            False,
            False,
            ("Primer3 CLI", "NUPACK-manual-only"),
            primer3,
            "Requires an independently installed Primer3 ntthal executable and explicit path.",
        ),
        ConditionalCapability(
            "THERMO-009",
            "self_dimer",
            status,
            execution_supported,
            False,
            False,
            ("Primer3 CLI", "NUPACK-manual-only"),
            primer3,
            "Requires an independently installed Primer3 ntthal executable and explicit path.",
        ),
        ConditionalCapability(
            "THERMO-010",
            "heterodimer",
            status,
            execution_supported,
            False,
            False,
            ("Primer3 CLI", "NUPACK-manual-only"),
            primer3,
            "Requires an independently installed Primer3 ntthal executable and explicit path.",
        ),
        ConditionalCapability(
            "THERMO-011",
            "secondary_structure",
            ("available-user-installed-adapter" if nupack.available else "conditional-unavailable"),
            nupack.available,
            False,
            False,
            ("NUPACK",),
            nupack,
            (
                "A bounded explicit adapter is available, but NUPACK must be separately "
                "licensed and installed; DNAKit never installs or downloads it automatically."
            ),
        ),
    )


__all__ = [
    "Primer3CLIAdapter",
    "ThermodynamicsBackend",
    "conditional_capabilities",
    "probe_primer3",
    "validate_primer3_result",
]
