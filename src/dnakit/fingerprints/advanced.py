"""Interpretable fixed-schema advanced DNA fingerprints and preprocessing."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import islice
from typing import Literal, Protocol, cast

from dnakit.core import DNASequence, Provenance, Topology
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.descriptors import codon_statistics, exact_repeat_fraction, gc_at_content
from dnakit.exceptions import (
    BackendExecutionError,
    ConfigurationError,
    UnsupportedGapOperationError,
)
from dnakit.patterns import (
    RestrictionEnzyme,
    scan_codon_sites,
    scan_motif,
    scan_restriction_sites,
)
from dnakit.thermodynamics import (
    Primer3ThermodynamicResult,
    ThermodynamicConditions,
    nearest_neighbor,
    validate_primer3_result,
)

from ._shared import SequenceInput, sequence_and_id
from .kmer import kmer as kmer_fingerprint

AdvancedVector = tuple[float, ...]
PreprocessMode = Literal["none", "standard", "minmax", "l1", "l2"]
ThermodynamicMissingStrategy = Literal["zero", "sentinel", "error"]
DEFAULT_MAX_PREPROCESS_ROWS = 1_000_000
DEFAULT_MAX_MULTISCALE_LEVELS = 16


class ThermodynamicStructureAdapter(Protocol):
    """Explicit structure adapter used by the thermodynamic fingerprint."""

    def hairpin(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult: ...

    def self_dimer(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult: ...

    def heterodimer(
        self,
        sequence_a: DNASequence,
        sequence_b: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
        max_loop: int = 30,
        output_structure: bool = False,
    ) -> Primer3ThermodynamicResult: ...


@dataclass(frozen=True)
class AdvancedFingerprintResult:
    """One fixed-order, JSON-compatible interpretable feature vector."""

    name: str
    method: str
    schema_version: str
    sequence_id: str | None
    feature_names: tuple[str, ...]
    values: AdvancedVector
    parameters: FrozenDict
    algorithm_version: str
    provenance: Provenance

    def __post_init__(self) -> None:
        for name in ("name", "method", "schema_version", "algorithm_version"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ConfigurationError(f"Advanced fingerprint {name} must be non-empty.")
        if len(self.feature_names) != len(self.values) or len(set(self.feature_names)) != len(
            self.feature_names
        ):
            raise ConfigurationError("Advanced fingerprint schema and values do not align.")
        if any(not isinstance(name, str) or not name for name in self.feature_names):
            raise ConfigurationError("Advanced fingerprint names must be non-empty strings.")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            for value in self.values
        ):
            raise ConfigurationError("Advanced fingerprint values must be finite numbers.")
        if self.sequence_id is not None and (
            not isinstance(self.sequence_id, str) or not self.sequence_id.strip()
        ):
            raise ConfigurationError("Advanced fingerprint sequence_id must be non-empty or None.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("Advanced fingerprint parameters must be FrozenDict.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("Advanced fingerprint provenance must be Provenance.")

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], to_json_compatible(self))


@dataclass(frozen=True)
class FeaturePreprocessor:
    """Training-set-fitted deterministic column transform."""

    feature_names: tuple[str, ...]
    mode: PreprocessMode
    centers: tuple[float, ...]
    scales: tuple[float, ...]
    imputation_values: tuple[float, ...]
    kept_indices: tuple[int, ...]
    variance_threshold: float
    missing_strategy: Literal["error", "mean", "zero"]
    max_rows: int = DEFAULT_MAX_PREPROCESS_ROWS
    schema_version: str = "dnakit.feature-preprocessor.v2"

    def __post_init__(self) -> None:
        width = len(self.feature_names)
        if not self.feature_names or any(
            not isinstance(name, str) or not name.strip() for name in self.feature_names
        ):
            raise ConfigurationError("Preprocessor feature names must be non-empty strings.")
        if len(set(self.feature_names)) != width:
            raise ConfigurationError("Preprocessor feature names must be unique.")
        if self.mode not in {"none", "standard", "minmax", "l1", "l2"}:
            raise ConfigurationError("Unknown preprocessing mode.")
        if self.missing_strategy not in {"error", "mean", "zero"}:
            raise ConfigurationError("Unknown missing strategy.")
        if not all(
            isinstance(values, tuple) and len(values) == width
            for values in (self.centers, self.scales, self.imputation_values)
        ):
            raise ConfigurationError("Preprocessor numeric vectors must match its feature schema.")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            for values in (self.centers, self.scales, self.imputation_values)
            for value in values
        ):
            raise ConfigurationError("Preprocessor numeric vectors must contain finite numbers.")
        if any(value <= 0 for value in self.scales):
            raise ConfigurationError("Preprocessor scales must be positive.")
        if (
            any(
                isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < width
                for index in self.kept_indices
            )
            or tuple(sorted(set(self.kept_indices))) != self.kept_indices
        ):
            raise ConfigurationError("Preprocessor kept_indices must be sorted unique indices.")
        if (
            isinstance(self.variance_threshold, bool)
            or not isinstance(self.variance_threshold, (int, float))
            or not math.isfinite(self.variance_threshold)
            or self.variance_threshold < 0
        ):
            raise ConfigurationError("variance_threshold must be finite and non-negative.")
        if (
            isinstance(self.max_rows, bool)
            or not isinstance(self.max_rows, int)
            or not 0 < self.max_rows <= DEFAULT_MAX_PREPROCESS_ROWS
        ):
            raise ConfigurationError(f"max_rows must be in [1, {DEFAULT_MAX_PREPROCESS_ROWS}].")
        if self.schema_version != "dnakit.feature-preprocessor.v2":
            raise ConfigurationError("Unsupported feature-preprocessor schema version.")

    def transform(self, rows: Iterable[Iterable[int | float | None]]) -> tuple[AdvancedVector, ...]:
        return _transform_rows(rows, self)

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], to_json_compatible(self))


def _result(
    value: SequenceInput,
    *,
    name: str,
    method: str,
    feature_names: Sequence[str],
    values: Sequence[int | float],
    parameters: Mapping[str, object],
    schema_version: str | None = None,
    provenance: Provenance | None = None,
) -> AdvancedFingerprintResult:
    _, sequence_id = sequence_and_id(value)
    return AdvancedFingerprintResult(
        name,
        method,
        schema_version or f"dnakit.{name}.v1",
        sequence_id,
        tuple(feature_names),
        tuple(float(item) for item in values),
        FrozenDict(parameters),
        "1.0",
        provenance or Provenance(),
    )


def motif_fingerprint(
    value: SequenceInput,
    motifs: Mapping[str, str],
    *,
    binary: bool = False,
) -> AdvancedFingerprintResult:
    """Encode ordered exact/IUPAC motif counts or presence."""

    if not isinstance(motifs, Mapping) or not motifs:
        raise ConfigurationError("motifs must be a non-empty ordered mapping.")
    if any(not isinstance(name, str) or not name.strip() for name in motifs):
        raise ConfigurationError("motif names must be non-empty strings.")
    if not isinstance(binary, bool):
        raise ConfigurationError("binary must be a boolean.")
    names = tuple(sorted(motifs))
    results = [
        scan_motif(value, motifs[name], mode="iupac", name=name, merge_strands=True)
        for name in names
    ]
    counts = tuple(len(result.hits) for result in results)
    return _result(
        value,
        name="motif_fingerprint",
        method="iupac-motif-presence" if binary else "iupac-motif-count",
        feature_names=tuple(f"motif:{name}" for name in names),
        values=tuple(int(count > 0) for count in counts) if binary else counts,
        parameters={"motifs": {name: motifs[name] for name in names}, "binary": binary},
    )


def restriction_fingerprint(
    value: SequenceInput,
    enzymes: Iterable[str | RestrictionEnzyme],
    *,
    binary: bool = False,
    max_enzymes: int = 10_000,
) -> AdvancedFingerprintResult:
    """Encode recognition-site counts for a versioned/user-defined panel."""

    if isinstance(enzymes, (str, bytes)):
        raise ConfigurationError("enzymes must be an iterable of enzyme definitions.")
    if not isinstance(binary, bool):
        raise ConfigurationError("binary must be a boolean.")
    scan = scan_restriction_sites(value, enzymes, max_enzymes=max_enzymes)
    requested = scan.parameters.get("enzymes")
    if not isinstance(requested, tuple):
        raise AssertionError("Restriction scan must audit its resolved enzyme panel.")
    resolved: list[str] = []
    for item in requested:
        if not isinstance(item, Mapping):
            raise AssertionError("Restriction scan enzyme audit entries must be mappings.")
        name = item.get("name")
        if not isinstance(name, str):
            raise AssertionError("Restriction scan enzyme audit entries must contain names.")
        resolved.append(name)
    resolved_names = tuple(resolved)
    counts = {name: 0 for name in resolved_names}
    for hit in scan.hits:
        counts[hit.enzyme] += 1
    return _result(
        value,
        name="restriction_fingerprint",
        method="restriction-site-presence" if binary else "restriction-site-count",
        feature_names=tuple(f"restriction:{name}" for name in resolved_names),
        values=tuple(int(counts[name] > 0) if binary else counts[name] for name in resolved_names),
        parameters={"enzymes": resolved_names, "binary": binary, "max_enzymes": max_enzymes},
    )


def gc_spatial_fingerprint(
    value: SequenceInput,
    *,
    bins: int = 10,
) -> AdvancedFingerprintResult:
    """Pool canonical GC fraction over a fixed number of normalized-position bins."""

    if isinstance(bins, bool) or not isinstance(bins, int) or bins <= 0 or bins > 10_000:
        raise ConfigurationError("bins must be an integer in [1, 10000].")
    sequence, _ = sequence_and_id(value)
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Spatial fingerprints cannot silently omit explicit Gap objects.",
            code="SPATIAL_FINGERPRINT_GAP_NOT_ALLOWED",
        )
    text = sequence.symbols
    fractions: list[float] = []
    for index in range(bins):
        start = len(text) * index // bins
        end = len(text) * (index + 1) // bins
        window = text[start:end]
        canonical = tuple(symbol for symbol in window if symbol in "ACGT")
        fractions.append(
            sum(symbol in "GC" for symbol in canonical) / len(canonical) if canonical else 0.0
        )
    return _result(
        value,
        name="gc_spatial_fingerprint",
        method="equal-symbol-bin-gc-fraction",
        feature_names=tuple(f"gc_bin:{index}" for index in range(bins)),
        values=fractions,
        parameters={
            "bins": bins,
            "empty_bin_value": 0.0,
            "ambiguity_policy": "ignore",
            "topology": sequence.topology.value,
            "circular_origin": "stored-origin" if sequence.topology is Topology.CIRCULAR else None,
        },
    )


def repeat_fingerprint(value: SequenceInput) -> AdvancedFingerprintResult:
    """Encode homopolymer/tandem-repeat summary values."""

    repeat = exact_repeat_fraction(value)
    from dnakit.descriptors import homopolymer_runs

    homopolymer = homopolymer_runs(value, ambiguity_policy="ignore")
    return _result(
        value,
        name="repeat_fingerprint",
        method="interpretable-exact-repeat-summary",
        feature_names=(
            "repeat_fraction",
            "repeat_run_count",
            "longest_homopolymer",
            "homopolymer_run_count",
        ),
        values=(
            repeat.repeat_fraction,
            len(repeat.runs),
            homopolymer.longest_length,
            len(homopolymer.runs),
        ),
        parameters={"repeat_method": repeat.method, "homopolymer_method": homopolymer.method},
    )


def coding_fingerprint(
    value: SequenceInput,
    *,
    frame: int = 0,
    genetic_code: int = 1,
) -> AdvancedFingerprintResult:
    """Encode codon/start/stop/ORF summaries without activity prediction."""

    codons = codon_statistics(value, frame=frame, genetic_code=genetic_code)
    sites = scan_codon_sites(value, genetic_code=genetic_code)
    return _result(
        value,
        name="coding_fingerprint",
        method="codon-and-site-summary",
        feature_names=("codon_count", "start_count", "stop_count", "six_frame_site_count"),
        values=(codons.codon_count, codons.start_count, codons.stop_count, len(sites.hits)),
        parameters={"frame": frame, "genetic_code": genetic_code},
    )


def thermodynamic_fingerprint(
    value: SequenceInput,
    *,
    conditions: ThermodynamicConditions | None = None,
    structure_adapter: ThermodynamicStructureAdapter | None = None,
    paired_value: SequenceInput | None = None,
    missing_strategy: ThermodynamicMissingStrategy = "zero",
    missing_value: float = -999.0,
    max_loop: int = 30,
) -> AdvancedFingerprintResult:
    """Encode a fixed v2 duplex/hairpin/dimer schema with explicit missing data."""

    sequence, _ = sequence_and_id(value)
    paired_sequence = None if paired_value is None else sequence_and_id(paired_value)[0]
    if missing_strategy not in {"zero", "sentinel", "error"}:
        raise ConfigurationError("Unknown thermodynamic fingerprint missing strategy.")
    if (
        isinstance(missing_value, bool)
        or not isinstance(missing_value, (int, float))
        or not math.isfinite(missing_value)
    ):
        raise ConfigurationError("missing_value must be finite numeric data.")
    if isinstance(max_loop, bool) or not isinstance(max_loop, int) or not 1 <= max_loop <= 30:
        raise ConfigurationError("max_loop must be an integer in [1, 30].")
    if structure_adapter is not None and any(
        not callable(getattr(structure_adapter, name, None))
        for name in ("hairpin", "self_dimer", "heterodimer")
    ):
        raise ConfigurationError(
            "structure_adapter must implement hairpin, self_dimer, and heterodimer.",
            code="INVALID_THERMODYNAMIC_STRUCTURE_ADAPTER",
        )
    nearest = nearest_neighbor(sequence, conditions=conditions)
    fill = 0.0 if missing_strategy == "zero" else float(missing_value)
    if missing_strategy == "error" and structure_adapter is None:
        raise ConfigurationError(
            "The hairpin fingerprint feature is unavailable under missing_strategy=error.",
            code="THERMODYNAMIC_FINGERPRINT_MISSING_FEATURE",
            context={"capability": "hairpin"},
        )
    if missing_strategy == "error" and paired_sequence is None:
        raise ConfigurationError(
            "The heterodimer fingerprint feature is unavailable under missing_strategy=error.",
            code="THERMODYNAMIC_FINGERPRINT_MISSING_FEATURE",
            context={"capability": "heterodimer"},
        )

    def missing_features(capability: str) -> tuple[float, float, float, float]:
        if missing_strategy == "error":
            raise ConfigurationError(
                f"The {capability} fingerprint feature is unavailable under "
                "missing_strategy=error.",
                code="THERMODYNAMIC_FINGERPRINT_MISSING_FEATURE",
                context={"capability": capability},
            )
        return (0.0, fill, fill, fill)

    validated_structure_results: list[Primer3ThermodynamicResult] = []

    def structure_features(
        result: object,
        *,
        capability: str,
        expected_sequences: tuple[str, ...],
    ) -> tuple[float, float, float, float]:
        validated = validate_primer3_result(
            result,
            capability=capability,
            sequences_5to3=expected_sequences,
            conditions=nearest.conditions,
            max_loop=max_loop,
            output_structure=False,
            error_code="MISMATCHED_THERMODYNAMIC_FINGERPRINT_RESULT",
        )
        if validated.structure_found is None or validated.delta_g_kcal_per_mol is None:
            raise BackendExecutionError(
                "Structure adapter omitted required fingerprint values.",
                code="INVALID_THERMODYNAMIC_FINGERPRINT_RESULT",
                context={"capability": capability},
            )
        validated_structure_results.append(validated)
        return (
            1.0,
            float(validated.structure_found),
            validated.tm_celsius,
            validated.delta_g_kcal_per_mol,
        )

    adapter_provenance: Provenance | None = None
    backend_name: str | None = None
    backend_version: str | None = None
    if structure_adapter is None:
        hairpin_values = missing_features("hairpin")
        self_dimer_values = missing_features("self_dimer")
        heterodimer_values = missing_features("heterodimer")
    else:
        hairpin_result = structure_adapter.hairpin(
            sequence,
            conditions=conditions,
            max_loop=max_loop,
            output_structure=False,
        )
        self_dimer_result = structure_adapter.self_dimer(
            sequence,
            conditions=conditions,
            max_loop=max_loop,
            output_structure=False,
        )
        hairpin_values = structure_features(
            hairpin_result,
            capability="hairpin",
            expected_sequences=(sequence.symbols,),
        )
        self_dimer_values = structure_features(
            self_dimer_result,
            capability="self_dimer",
            expected_sequences=(sequence.symbols,),
        )
        if paired_sequence is None:
            heterodimer_values = missing_features("heterodimer")
        else:
            heterodimer_values = structure_features(
                structure_adapter.heterodimer(
                    sequence,
                    paired_sequence,
                    conditions=conditions,
                    max_loop=max_loop,
                    output_structure=False,
                ),
                capability="heterodimer",
                expected_sequences=(sequence.symbols, paired_sequence.symbols),
            )
        adapter_provenance = hairpin_result.provenance
        backend_name = hairpin_result.backend.name
        backend_version = hairpin_result.backend.version
        if any(
            result.backend != hairpin_result.backend
            or result.provenance != hairpin_result.provenance
            for result in validated_structure_results[1:]
        ):
            raise BackendExecutionError(
                "Structure adapter returned inconsistent backend or provenance metadata.",
                code="INCONSISTENT_THERMODYNAMIC_FINGERPRINT_RESULTS",
            )
    return _result(
        value,
        name="thermodynamic_fingerprint",
        method="santalucia1998-duplex-with-explicit-primer3-structure-features",
        feature_names=(
            "tm_celsius",
            "delta_h",
            "delta_s",
            "delta_g",
            "hairpin_available",
            "hairpin_found",
            "hairpin_tm_celsius",
            "hairpin_delta_g",
            "self_dimer_available",
            "self_dimer_found",
            "self_dimer_tm_celsius",
            "self_dimer_delta_g",
            "heterodimer_available",
            "heterodimer_found",
            "heterodimer_tm_celsius",
            "heterodimer_delta_g",
        ),
        values=(
            nearest.tm_celsius,
            nearest.delta_h_kcal_per_mol,
            nearest.delta_s_cal_per_k_mol,
            nearest.delta_g_kcal_per_mol,
            *hairpin_values,
            *self_dimer_values,
            *heterodimer_values,
        ),
        parameters={
            "conditions": to_json_compatible(nearest.conditions),
            "parameter_set": nearest.parameter_set,
            "hairpin_dimer_included": structure_adapter is not None,
            "heterodimer_included": structure_adapter is not None and paired_sequence is not None,
            "structure_adapter_supplied": structure_adapter is not None,
            "automatic_backend_probe": False,
            "automatic_backend_execution": False,
            "missing_strategy": missing_strategy,
            "missing_value": fill,
            "paired_sequence_supplied": paired_sequence is not None,
            "max_loop": max_loop,
            "backend_name": backend_name,
            "backend_version": backend_version,
        },
        schema_version="dnakit.thermodynamic_fingerprint.v2",
        provenance=adapter_provenance,
    )


def multiscale_fingerprint(
    value: SequenceInput,
    *,
    k_values: Iterable[int] = (1, 2, 3),
) -> AdvancedFingerprintResult:
    """Concatenate frequency k-mer vectors at explicit scales plus global GC."""

    if isinstance(k_values, (str, bytes)):
        raise ConfigurationError("k_values must be an iterable of positive integers.")
    scales_list: list[int] = []
    try:
        iterator = iter(k_values)
    except TypeError as exc:
        raise ConfigurationError("k_values must be an iterable of positive integers.") from exc
    for index, scale in enumerate(iterator):
        if index >= DEFAULT_MAX_MULTISCALE_LEVELS:
            raise ConfigurationError(
                "k_values exceeds the hard multiscale-level limit.",
                code="MULTISCALE_LEVEL_LIMIT",
            )
        scales_list.append(scale)
    scales = tuple(scales_list)
    if not scales or any(isinstance(k, bool) or not isinstance(k, int) or k <= 0 for k in scales):
        raise ConfigurationError("k_values must contain positive integers.")
    if len(scales) != len(set(scales)):
        raise ConfigurationError("k_values must be unique.")
    names: list[str] = []
    values: list[float] = []
    for k in scales:
        result = kmer_fingerprint(value, k=k, mode="frequency")
        names.extend(f"k{k}:{name}" for name in result.feature_names)
        values.extend(float(item) for item in result.dense_values())
    gc = gc_at_content(value, ambiguity_policy="ignore")
    names.append("global_gc")
    values.append(0.0 if gc.gc_fraction is None else gc.gc_fraction)
    return _result(
        value,
        name="multiscale_fingerprint",
        method="concatenated-kmer-frequency-and-global-gc",
        feature_names=names,
        values=values,
        parameters={"k_values": scales},
    )


def hybrid_fingerprint(
    components: Mapping[str, AdvancedFingerprintResult],
    *,
    weights: Mapping[str, int | float] | None = None,
) -> AdvancedFingerprintResult:
    """Concatenate named compatible component results with optional scalar weights."""

    if not isinstance(components, Mapping) or not components:
        raise ConfigurationError("components must be a non-empty mapping.")
    if any(not isinstance(name, str) or not name.strip() for name in components):
        raise ConfigurationError("Component names must be non-empty strings.")
    if any(not isinstance(value, AdvancedFingerprintResult) for value in components.values()):
        raise ConfigurationError("Component values must be AdvancedFingerprintResult objects.")
    if weights is not None and not isinstance(weights, Mapping):
        raise ConfigurationError("weights must be a mapping or None.")
    ordered = tuple(sorted(components))
    resolved_weights = {name: 1.0 for name in ordered}
    for name, raw in (weights or {}).items():
        if name not in resolved_weights:
            raise ConfigurationError("Weight names must match component names.")
        if (
            isinstance(raw, bool)
            or not isinstance(raw, (int, float))
            or not math.isfinite(raw)
            or raw < 0
        ):
            raise ConfigurationError("Component weights must be finite non-negative numbers.")
        resolved_weights[name] = float(raw)
    sequence_ids = {components[name].sequence_id for name in ordered}
    if len(sequence_ids) > 1:
        raise ConfigurationError("Hybrid components must refer to the same sequence ID.")
    names = tuple(
        f"{component}:{feature}"
        for component in ordered
        for feature in components[component].feature_names
    )
    values = tuple(
        value * resolved_weights[component]
        for component in ordered
        for value in components[component].values
    )
    return AdvancedFingerprintResult(
        "hybrid_fingerprint",
        "weighted-component-concatenation",
        "dnakit.hybrid_fingerprint.v1",
        next(iter(sequence_ids)),
        names,
        values,
        FrozenDict(
            {
                "components": {name: components[name].schema_version for name in ordered},
                "weights": resolved_weights,
            }
        ),
        "1.0",
        Provenance(),
    )


def _materialize_rows(
    rows: Iterable[Iterable[int | float | None]],
    *,
    max_rows: int,
    width: int,
) -> tuple[tuple[float | None, ...], ...]:
    output: list[tuple[float | None, ...]] = []
    for index, row in enumerate(rows):
        if index >= max_rows:
            raise ConfigurationError("Feature rows exceed max_rows.")
        try:
            values = tuple(islice(iter(row), width + 1))
        except TypeError as exc:
            raise ConfigurationError("Each feature row must be an iterable.") from exc
        if len(values) != width:
            raise ConfigurationError("Feature row width does not match the fitted schema.")
        for value in values:
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ConfigurationError("Feature values must be finite numbers or None.")
        output.append(tuple(None if item is None else float(item) for item in values))
    if not output:
        raise ConfigurationError("Feature matrix must contain at least one row.")
    return tuple(output)


def fit_preprocessor(
    rows: Iterable[Iterable[int | float | None]],
    *,
    feature_names: Sequence[str],
    mode: PreprocessMode = "standard",
    missing_strategy: Literal["error", "mean", "zero"] = "error",
    variance_threshold: float = 0.0,
    max_rows: int = DEFAULT_MAX_PREPROCESS_ROWS,
) -> FeaturePreprocessor:
    """Fit imputation/scaling/low-variance parameters on training rows only."""

    if isinstance(feature_names, (str, bytes)):
        raise ConfigurationError("feature_names must be a sequence of names, not text.")
    names = tuple(feature_names)
    if (
        not names
        or any(not isinstance(name, str) or not name.strip() for name in names)
        or len(set(names)) != len(names)
    ):
        raise ConfigurationError("feature_names must be unique non-empty strings.")
    if mode not in {"none", "standard", "minmax", "l1", "l2"}:
        raise ConfigurationError("Unknown preprocessing mode.")
    if missing_strategy not in {"error", "mean", "zero"}:
        raise ConfigurationError("Unknown missing strategy.")
    if (
        isinstance(variance_threshold, bool)
        or not isinstance(variance_threshold, (int, float))
        or not math.isfinite(variance_threshold)
        or variance_threshold < 0
    ):
        raise ConfigurationError("variance_threshold must be finite and non-negative.")
    if (
        isinstance(max_rows, bool)
        or not isinstance(max_rows, int)
        or not 0 < max_rows <= DEFAULT_MAX_PREPROCESS_ROWS
    ):
        raise ConfigurationError(f"max_rows must be in [1, {DEFAULT_MAX_PREPROCESS_ROWS}].")
    matrix = _materialize_rows(rows, max_rows=max_rows, width=len(names))
    columns = tuple(tuple(row[index] for row in matrix) for index in range(len(names)))
    centers: list[float] = []
    scales: list[float] = []
    imputations: list[float] = []
    variances: list[float] = []
    for column in columns:
        observed = tuple(value for value in column if value is not None)
        if missing_strategy == "error" and len(observed) != len(column):
            raise ConfigurationError("Missing feature encountered under error policy.")
        center = math.fsum(observed) / len(observed) if observed else 0.0
        imputation = center if missing_strategy == "mean" else 0.0
        imputed = tuple(imputation if value is None else value for value in column)
        mean = math.fsum(imputed) / len(imputed)
        variance = math.fsum((value - mean) ** 2 for value in imputed) / len(imputed)
        if mode == "standard":
            centers.append(mean)
            scales.append(math.sqrt(variance) or 1.0)
        elif mode == "minmax":
            minimum, maximum = min(imputed), max(imputed)
            centers.append(minimum)
            scales.append(maximum - minimum or 1.0)
        else:
            centers.append(0.0)
            scales.append(1.0)
        imputations.append(imputation)
        variances.append(variance)
    kept = tuple(
        index for index, variance in enumerate(variances) if variance >= variance_threshold
    )
    return FeaturePreprocessor(
        names,
        mode,
        tuple(centers),
        tuple(scales),
        tuple(imputations),
        kept,
        float(variance_threshold),
        missing_strategy,
        max_rows,
    )


def _transform_rows(
    rows: Iterable[Iterable[int | float | None]], preprocessor: FeaturePreprocessor
) -> tuple[AdvancedVector, ...]:
    if not isinstance(preprocessor, FeaturePreprocessor):
        raise ConfigurationError("preprocessor must be FeaturePreprocessor.")
    matrix = _materialize_rows(
        rows,
        max_rows=preprocessor.max_rows,
        width=len(preprocessor.feature_names),
    )
    transformed: list[AdvancedVector] = []
    for row in matrix:
        values: list[float] = []
        for index, raw in enumerate(row):
            if raw is None:
                if preprocessor.missing_strategy == "error":
                    raise ConfigurationError("Missing feature encountered under error policy.")
                raw = preprocessor.imputation_values[index]
            values.append((raw - preprocessor.centers[index]) / preprocessor.scales[index])
        selected = [values[index] for index in preprocessor.kept_indices]
        if preprocessor.mode in {"l1", "l2"}:
            norm = (
                math.fsum(abs(value) for value in selected)
                if preprocessor.mode == "l1"
                else math.sqrt(math.fsum(value * value for value in selected))
            )
            if norm:
                selected = [value / norm for value in selected]
        transformed.append(tuple(selected))
    return tuple(transformed)


__all__ = [
    "DEFAULT_MAX_MULTISCALE_LEVELS",
    "DEFAULT_MAX_PREPROCESS_ROWS",
    "AdvancedFingerprintResult",
    "FeaturePreprocessor",
    "ThermodynamicMissingStrategy",
    "ThermodynamicStructureAdapter",
    "coding_fingerprint",
    "fit_preprocessor",
    "gc_spatial_fingerprint",
    "hybrid_fingerprint",
    "motif_fingerprint",
    "multiscale_fingerprint",
    "repeat_fingerprint",
    "restriction_fingerprint",
    "thermodynamic_fingerprint",
]
