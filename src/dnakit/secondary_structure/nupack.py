"""Explicit adapter for a separately licensed and installed NUPACK 4 module."""

from __future__ import annotations

import importlib
import importlib.util
import math
from collections.abc import Iterable, Mapping
from importlib.metadata import PackageNotFoundError, distribution, version
from itertools import islice
from typing import Any, Literal, cast

from dnakit.core import (
    BackendInfo,
    Citation,
    DNASequence,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
)
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.thermodynamics import ThermodynamicConditions

from .dotbracket import _strands, analyze_dot_bracket, pair_probability_metrics
from .results import (
    ComplexConcentration,
    NupackComplexResult,
    NupackTubeResult,
    PredictedSecondaryStructure,
)

_NUPACK_CAPABILITIES = frozenset(
    {
        "pfunc",
        "mfe",
        "pairs",
        "subopt",
        "sample",
        "ensemble_size",
        "structure_probability",
        "ensemble_defect",
        "tube_concentrations",
    }
)
_MAX_SAMPLES = 100_000
_MAX_TUBE_STRANDS = 20
_MAX_COMPLEX_SIZE = 4


def probe_nupack() -> BackendInfo:
    """Passively locate NUPACK without importing, installing, or downloading it."""

    available = importlib.util.find_spec("nupack") is not None
    package_version: str | None = None
    package_location: str | None = None
    if available:
        try:
            package_version = version("nupack")
            package_location = str(distribution("nupack").locate_file("nupack"))
        except PackageNotFoundError:
            available = False
    return BackendInfo(
        "nupack",
        version=package_version,
        package_location=package_location,
        license_expression="LicenseRef-NUPACK",
        capabilities=_NUPACK_CAPABILITIES,
        available=available,
        metadata={
            "adapter_status": "execution-enabled-only-when-user-licensed-and-installed",
            "automatic_install": False,
            "automatic_download": False,
            "import_executed": False,
            "separate_license_or_subscription_required": True,
        },
    )


def _provenance(info: BackendInfo) -> Provenance:
    dependencies = {"nupack": info.version} if info.version is not None else {}
    return Provenance(
        dependency_versions=dependencies,
        implementation=ImplementationInfo(
            label=ImplementationLabel.ADAPTER,
            execution_mode=ExecutionMode.HYBRID,
            origin_class=OriginClass.INTEGRATION,
            license_expression=info.license_expression,
            citations=(
                Citation(
                    "nupack2011",
                    title="NUPACK: analysis and design of nucleic acid systems",
                    doi="10.1002/jcc.21596",
                ),
            ),
        ),
        backend=info,
    )


def _finite(value: object, field: str) -> float:
    if isinstance(value, (bool, str, bytes)):
        raise BackendExecutionError(
            f"NUPACK output {field} must be finite numeric data.",
            code="INVALID_NUPACK_OUTPUT",
        )
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError, OverflowError) as exc:
        raise BackendExecutionError(
            f"NUPACK output {field} must be finite numeric data.",
            code="INVALID_NUPACK_OUTPUT",
        ) from exc
    if not math.isfinite(result):
        raise BackendExecutionError(
            f"NUPACK output {field} must be finite numeric data.",
            code="INVALID_NUPACK_OUTPUT",
        )
    return result


def _input_finite(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"{field} must be finite numeric data.",
            code="INVALID_NUPACK_CONFIGURATION",
            context={"field": field, "value": value},
        )
    return float(value)


def _structure_text(value: object) -> str:
    text = str(value)
    if not text or len(text) > 100_000:
        raise BackendExecutionError(
            "NUPACK returned an invalid structure string.", code="INVALID_NUPACK_OUTPUT"
        )
    return text


class NupackAdapter:
    """NUPACK 4 utilities and tube-analysis adapter with bounded outputs."""

    def __init__(self) -> None:
        self._info = probe_nupack()

    @property
    def info(self) -> BackendInfo:
        return self._info

    def ensure_available(self) -> None:
        if not self.info.available:
            raise BackendUnavailableError(
                "NUPACK is not installed in this environment.",
                code="NUPACK_UNAVAILABLE",
                hint=(
                    "Obtain an appropriate NUPACK license/subscription, install the downloaded "
                    "package locally, then call this adapter explicitly."
                ),
            )

    def _module(self) -> Any:
        self.ensure_available()
        try:
            return importlib.import_module("nupack")
        except Exception as exc:
            raise BackendExecutionError(
                "The installed NUPACK module could not be imported.",
                code="NUPACK_IMPORT_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc

    @staticmethod
    def _conditions(
        value: ThermodynamicConditions | None,
    ) -> ThermodynamicConditions:
        resolved = ThermodynamicConditions(sodium_molar=1.0) if value is None else value
        if not isinstance(resolved, ThermodynamicConditions):
            raise ConfigurationError(
                "conditions must be ThermodynamicConditions or None.",
                code="INVALID_NUPACK_CONDITIONS",
            )
        if not 0.05 <= resolved.monovalent_molar <= 1.1:
            raise ConfigurationError(
                "NUPACK DNA monovalent concentration must be in [0.05, 1.1] M.",
                code="NUPACK_SALT_OUT_OF_DOMAIN",
            )
        if not 0.0 <= resolved.magnesium_molar <= 0.2:
            raise ConfigurationError(
                "NUPACK DNA magnesium concentration must be in [0, 0.2] M.",
                code="NUPACK_SALT_OUT_OF_DOMAIN",
            )
        if (
            resolved.dntp_molar != 0.0
            or resolved.dmso_percent != 0.0
            or resolved.formamide_molar != 0.0
        ):
            raise ConfigurationError(
                "The NUPACK adapter does not map dNTP, DMSO, or formamide fields.",
                code="NUPACK_CONDITION_UNSUPPORTED",
            )
        return resolved

    @staticmethod
    def _model(module: Any, conditions: ThermodynamicConditions, ensemble: str) -> object:
        if ensemble not in {"stacking", "dangle-stacking", "coaxial-stacking", "nostacking"}:
            raise ConfigurationError(
                "Unknown NUPACK ensemble.",
                code="INVALID_NUPACK_ENSEMBLE",
                context={"ensemble": ensemble},
            )
        try:
            return module.Model(
                material="dna",
                ensemble=ensemble,
                celsius=conditions.temperature_celsius,
                sodium=conditions.monovalent_molar,
                magnesium=conditions.magnesium_molar,
            )
        except Exception as exc:
            raise BackendExecutionError(
                "NUPACK model construction failed.",
                code="NUPACK_MODEL_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc

    @staticmethod
    def _predicted_structures(
        raw_values: object,
        strands: tuple[DNASequence, ...],
        *,
        maximum: int,
    ) -> tuple[PredictedSecondaryStructure, ...]:
        if isinstance(raw_values, (str, bytes)) or not isinstance(raw_values, Iterable):
            raise BackendExecutionError(
                "NUPACK returned a non-iterable structure collection.",
                code="INVALID_NUPACK_OUTPUT",
            )
        items = tuple(islice(iter(raw_values), maximum + 1))
        if len(items) > maximum:
            raise BackendExecutionError(
                "NUPACK structure output exceeded the configured bound.",
                code="NUPACK_OUTPUT_LIMIT_EXCEEDED",
            )
        results: list[PredictedSecondaryStructure] = []
        for item in items:
            typed_item = cast(Any, item)
            try:
                structure_text = _structure_text(typed_item.structure)
                energy = _finite(typed_item.energy, "energy")
                raw_stack_energy = typed_item.stack_energy
                stack_energy = (
                    None if raw_stack_energy is None else _finite(raw_stack_energy, "stack_energy")
                )
            except AttributeError as exc:
                raise BackendExecutionError(
                    "NUPACK structure output omitted a required field.",
                    code="INVALID_NUPACK_OUTPUT",
                ) from exc
            results.append(
                PredictedSecondaryStructure(
                    summary=analyze_dot_bracket(strands, structure_text),
                    free_energy_kcal_per_mol=energy,
                    stack_free_energy_kcal_per_mol=stack_energy,
                )
            )
        return tuple(results)

    def analyze_complex(
        self,
        strands: Iterable[DNASequence],
        *,
        conditions: ThermodynamicConditions | None = None,
        ensemble: Literal[
            "stacking", "dangle-stacking", "coaxial-stacking", "nostacking"
        ] = "stacking",
        target_structure: str | None = None,
        suboptimal_energy_gap_kcal_per_mol: float = 1.0,
        num_samples: int = 100,
        accessibility_window_size: int = 4,
    ) -> NupackComplexResult:
        """Run bounded NUPACK single-complex ensemble utilities."""

        strand_symbols = _strands(strands)
        strand_objects = tuple(DNASequence(item) for item in strand_symbols)
        resolved = self._conditions(conditions)
        if (
            isinstance(num_samples, bool)
            or not isinstance(num_samples, int)
            or not 0 <= num_samples <= _MAX_SAMPLES
        ):
            raise ConfigurationError(
                f"num_samples must be an integer in [0, {_MAX_SAMPLES}].",
                code="INVALID_NUPACK_SAMPLE_COUNT",
            )
        gap = _input_finite(suboptimal_energy_gap_kcal_per_mol, "suboptimal_energy_gap")
        if not 0.0 <= gap <= 20.0:
            raise ConfigurationError(
                "suboptimal_energy_gap_kcal_per_mol must be in [0, 20].",
                code="INVALID_NUPACK_ENERGY_GAP",
            )
        target = (
            None
            if target_structure is None
            else analyze_dot_bracket(strand_objects, target_structure)
        )
        module = self._module()
        model = self._model(module, resolved, ensemble)
        sequence_list = list(strand_symbols)
        try:
            pfunc_raw = module.pfunc(strands=sequence_list, model=model)
            ensemble_free_energy = _finite(pfunc_raw[1], "ensemble_free_energy")
            partition_function = pfunc_raw[0]
            if hasattr(partition_function, "log"):
                partition_log = _finite(partition_function.log(), "partition_function_log")
            else:
                partition_log = math.log(_finite(partition_function, "partition_function"))
            mfe_raw = module.mfe(strands=sequence_list, model=model)
            pairs_raw = module.pairs(strands=sequence_list, model=model)
            array = pairs_raw.to_array()
            matrix = array.tolist() if hasattr(array, "tolist") else array
            suboptimal_raw = module.subopt(
                strands=sequence_list,
                energy_gap=gap,
                model=model,
            )
            samples_raw = (
                []
                if num_samples == 0
                else module.sample(
                    strands=sequence_list,
                    num_sample=num_samples,
                    model=model,
                )
            )
            ensemble_size = int(module.ensemble_size(strands=sequence_list, model=model))
            target_probability = (
                None
                if target_structure is None
                else _finite(
                    module.structure_probability(
                        strands=sequence_list,
                        structure=target_structure,
                        model=model,
                    ),
                    "target_structure_probability",
                )
            )
            target_defect = (
                None
                if target_structure is None
                else _finite(
                    module.defect(
                        strands=sequence_list,
                        structure=target_structure,
                        model=model,
                    ),
                    "target_ensemble_defect",
                )
            )
        except (BackendExecutionError, ConfigurationError):
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "NUPACK complex analysis failed.",
                code="NUPACK_EXECUTION_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc
        if ensemble_size < 1:
            raise BackendExecutionError(
                "NUPACK ensemble_size must be positive.", code="INVALID_NUPACK_OUTPUT"
            )
        mfe_structures = self._predicted_structures(mfe_raw, strand_objects, maximum=100)
        if not mfe_structures:
            raise BackendExecutionError(
                "NUPACK returned no MFE structure.", code="INVALID_NUPACK_OUTPUT"
            )
        suboptimal = self._predicted_structures(
            suboptimal_raw,
            strand_objects,
            maximum=100_000,
        )
        if isinstance(samples_raw, (str, bytes)) or not isinstance(samples_raw, Iterable):
            raise BackendExecutionError(
                "NUPACK samples output is invalid.", code="INVALID_NUPACK_OUTPUT"
            )
        sample_items = tuple(islice(iter(samples_raw), num_samples + 1))
        if len(sample_items) != num_samples:
            raise BackendExecutionError(
                "NUPACK sample count does not match the request.",
                code="INVALID_NUPACK_OUTPUT",
            )
        samples = tuple(_structure_text(item) for item in sample_items)
        probability_result = pair_probability_metrics(
            strand_objects,
            cast(Iterable[Iterable[float]], matrix),
            accessibility_window_size=accessibility_window_size,
        )
        if target_probability is not None and not 0.0 <= target_probability <= 1.0:
            raise BackendExecutionError(
                "NUPACK target probability is outside [0, 1].",
                code="INVALID_NUPACK_OUTPUT",
            )
        if target_defect is not None and not 0.0 <= target_defect <= 1.0:
            raise BackendExecutionError(
                "NUPACK normalized ensemble defect is outside [0, 1].",
                code="INVALID_NUPACK_OUTPUT",
            )
        return NupackComplexResult(
            strands_5to3=strand_symbols,
            material="dna",
            temperature_celsius=resolved.temperature_celsius,
            monovalent_molar=resolved.monovalent_molar,
            magnesium_molar=resolved.magnesium_molar,
            ensemble=ensemble,
            partition_function_log=partition_log,
            ensemble_free_energy_kcal_per_mol=ensemble_free_energy,
            mfe_structures=mfe_structures,
            pair_probabilities=probability_result,
            suboptimal_structures=suboptimal,
            boltzmann_samples=samples,
            ensemble_size=ensemble_size,
            target_structure=target,
            target_structure_probability=target_probability,
            target_ensemble_defect=target_defect,
            method="nupack-4-utilities-adapter-v1",
            backend=self.info,
            provenance=_provenance(self.info),
        )

    def analyze_tube(
        self,
        strands: Mapping[str, DNASequence],
        concentrations_molar: Mapping[str, float],
        *,
        target_strand_names: Iterable[str],
        conditions: ThermodynamicConditions | None = None,
        max_complex_size: int = 2,
    ) -> NupackTubeResult:
        """Calculate target and off-target equilibrium complex concentrations."""

        if not isinstance(strands, Mapping) or not strands:
            raise ConfigurationError(
                "strands must be a non-empty name-to-DNASequence mapping.",
                code="INVALID_NUPACK_TUBE_STRANDS",
            )
        if len(strands) > _MAX_TUBE_STRANDS or any(
            not isinstance(name, str) or not name.strip() or not isinstance(sequence, DNASequence)
            for name, sequence in strands.items()
        ):
            raise ConfigurationError(
                f"NUPACK tube strand names and sequences are invalid or exceed "
                f"the {_MAX_TUBE_STRANDS}-strand limit.",
                code="INVALID_NUPACK_TUBE_STRANDS",
            )
        if set(concentrations_molar) != set(strands):
            raise ConfigurationError(
                "concentrations_molar keys must exactly match strands.",
                code="NUPACK_TUBE_CONCENTRATION_KEYS",
            )
        names = tuple(strands)
        sequences = _strands(tuple(strands[name] for name in names))
        concentrations = tuple(_input_finite(concentrations_molar[name], name) for name in names)
        if any(value <= 0.0 or value > 1.0 for value in concentrations):
            raise ConfigurationError(
                "Every strand concentration must be in (0, 1] M.",
                code="INVALID_NUPACK_TUBE_CONCENTRATION",
            )
        if isinstance(target_strand_names, (str, bytes)) or not isinstance(
            target_strand_names, Iterable
        ):
            raise ConfigurationError(
                "target_strand_names must be an iterable.",
                code="INVALID_NUPACK_TARGET_COMPLEX",
            )
        target_names = tuple(islice(iter(target_strand_names), _MAX_COMPLEX_SIZE + 1))
        if not target_names or any(name not in strands for name in target_names):
            raise ConfigurationError(
                "Target complex contains an unknown strand name.",
                code="INVALID_NUPACK_TARGET_COMPLEX",
            )
        if (
            isinstance(max_complex_size, bool)
            or not isinstance(max_complex_size, int)
            or not 1 <= max_complex_size <= _MAX_COMPLEX_SIZE
            or len(target_names) > max_complex_size
        ):
            raise ConfigurationError(
                f"max_complex_size must be in [1, {_MAX_COMPLEX_SIZE}] and include the target.",
                code="INVALID_NUPACK_COMPLEX_SIZE",
            )
        resolved = self._conditions(conditions)
        module = self._module()
        model = self._model(module, resolved, "stacking")
        try:
            strand_objects = {
                name: module.Strand(sequence, name=name)
                for name, sequence in zip(names, sequences, strict=True)
            }
            target_complex = module.Complex(
                [strand_objects[name] for name in target_names], name="dnakit-target"
            )
            tube = module.Tube(
                strands={
                    strand_objects[name]: concentration
                    for name, concentration in zip(names, concentrations, strict=True)
                },
                complexes=module.SetSpec(max_size=max_complex_size, include=[target_complex]),
                name="dnakit-tube",
            )
            analysis = module.tube_analysis(tubes=[tube], model=model, compute=["pfunc", "pairs"])
            tube_result = analysis[tube]
            raw_concentrations = tube_result.complex_concentrations
            fraction_unpaired = _finite(
                tube_result.fraction_bases_unpaired, "fraction_bases_unpaired"
            )
        except (BackendExecutionError, ConfigurationError):
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "NUPACK tube analysis failed.",
                code="NUPACK_EXECUTION_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc
        complex_results: list[ComplexConcentration] = []
        target_concentration = 0.0
        for raw_complex, raw_concentration in raw_concentrations.items():
            concentration = _finite(raw_concentration, "complex_concentration")
            if concentration < 0.0 or concentration > math.fsum(concentrations) * (1.0 + 1e-6):
                raise BackendExecutionError(
                    "NUPACK returned an invalid complex concentration.",
                    code="INVALID_NUPACK_OUTPUT",
                )
            raw_name = getattr(raw_complex, "name", None)
            name = str(raw_name) if raw_name else str(raw_complex)
            raw_complex_strands = getattr(raw_complex, "strands", ())
            complex_strand_names = tuple(
                str(getattr(item, "name", item)) for item in raw_complex_strands
            )
            is_target = raw_complex is target_complex or name == "dnakit-target"
            if is_target:
                target_concentration += concentration
            complex_results.append(
                ComplexConcentration(
                    name=name,
                    strand_names=complex_strand_names,
                    concentration_molar=concentration,
                    is_target=is_target,
                )
            )
        total_complex_concentration = math.fsum(
            item.concentration_molar for item in complex_results
        )
        target_fraction = (
            0.0
            if total_complex_concentration == 0.0
            else target_concentration / total_complex_concentration
        )
        if not 0.0 <= target_fraction <= 1.0:
            raise BackendExecutionError(
                "NUPACK target complex fraction is outside [0, 1].",
                code="INVALID_NUPACK_OUTPUT",
            )
        if not 0.0 <= fraction_unpaired <= 1.0:
            raise BackendExecutionError(
                "NUPACK fraction_bases_unpaired is outside [0, 1].",
                code="INVALID_NUPACK_OUTPUT",
            )
        return NupackTubeResult(
            strand_names=names,
            sequences_5to3=sequences,
            input_concentrations_molar=concentrations,
            complex_concentrations=tuple(complex_results),
            target_complex_name="dnakit-target",
            target_strand_names=target_names,
            target_complex_concentration_molar=target_concentration,
            complex_fraction_denominator_molar=total_complex_concentration,
            target_complex_fraction=target_fraction,
            non_target_complex_fraction=1.0 - target_fraction,
            fraction_bases_unpaired=fraction_unpaired,
            max_complex_size=max_complex_size,
            method="nupack-4-tube-analysis-adapter-v1",
            backend=self.info,
            provenance=_provenance(self.info),
        )


__all__ = ["NupackAdapter", "probe_nupack"]
