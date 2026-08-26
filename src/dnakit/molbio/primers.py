"""Primer matching, PCR simulation, properties, and design request contracts."""

from __future__ import annotations

from collections.abc import Iterable
from itertools import islice
from typing import Protocol, cast

from dnakit.core import DNAAlphabet, DNASequence, Issue, IssueSeverity, Strand, Topology
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.thermodynamics import (
    MeltingTemperatureResult,
    NearestNeighborConfig,
    Primer3ThermodynamicResult,
    ThermodynamicConditions,
    melting_temperature,
    validate_primer3_result,
)

from ._shared import (
    adapter_provenance,
    circular_slice,
    finite_fraction,
    freeze_parameters,
    iupac_compatible,
    native_provenance,
    require_sequence,
    reverse_complement_text,
    validate_positive_int,
)
from .results import (
    Amplicon,
    ConditionalAnalysis,
    PCRResult,
    PrimerBindingHit,
    PrimerDesignInterfaceResult,
    PrimerDesignRequest,
    PrimerMatchResult,
    PrimerPropertiesResult,
)


class PrimerStructureAdapter(Protocol):
    """Explicit structure adapter accepted by :func:`primer_properties`."""

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


class _PrimerTmAdapter(Protocol):
    def tm(
        self,
        sequence: DNASequence,
        *,
        conditions: ThermodynamicConditions | None = None,
    ) -> Primer3ThermodynamicResult: ...


def _coerce_strands(value: Strand | str) -> tuple[Strand, ...]:
    try:
        resolved = value if isinstance(value, Strand) else Strand(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("strand must be forward, reverse, or both.") from exc
    if resolved is Strand.BOTH:
        return (Strand.FORWARD, Strand.REVERSE)
    if resolved not in (Strand.FORWARD, Strand.REVERSE):
        raise ConfigurationError("strand must be forward, reverse, or both.")
    return (resolved,)


def _comparison_positions(
    primer: str,
    template_window: str,
    strand: Strand,
    *,
    allow_iupac: bool,
) -> tuple[int, ...]:
    oriented = primer if strand is Strand.FORWARD else reverse_complement_text(primer)
    mismatches: list[int] = []
    for index, (query, target) in enumerate(zip(oriented, template_window, strict=True)):
        matches = iupac_compatible(query, target) if allow_iupac else query == target
        if not matches:
            primer_index = index if strand is Strand.FORWARD else len(primer) - index - 1
            mismatches.append(primer_index)
    return tuple(sorted(mismatches))


def match_primer(
    primer: DNASequence,
    template: DNASequence,
    *,
    primer_name: str = "primer",
    strand: Strand | str = Strand.BOTH,
    max_mismatches: int = 0,
    strict_three_prime_bases: int = 0,
    allow_iupac: bool = False,
    max_hits: int = 100_000,
    max_template_length: int = 10_000_000,
    max_comparison_cells: int = 100_000_000,
) -> PrimerMatchResult:
    """Find ungapped primer matches and report mismatches in primer 5'-to-3' order."""

    primer_symbols = require_sequence(
        primer,
        operation="primer matching",
        max_length=1_000,
        canonical=not allow_iupac,
        allow_circular=False,
    )
    template_symbols = require_sequence(
        template,
        operation="primer matching",
        max_length=max_template_length,
        canonical=not allow_iupac,
    )
    if not isinstance(primer_name, str) or not primer_name.strip():
        raise ConfigurationError("primer_name must be non-empty.")
    if not isinstance(allow_iupac, bool):
        raise ConfigurationError("allow_iupac must be boolean.")
    if isinstance(max_mismatches, bool) or not isinstance(max_mismatches, int):
        raise ConfigurationError("max_mismatches must be a non-negative integer.")
    if not 0 <= max_mismatches <= len(primer_symbols):
        raise ConfigurationError("max_mismatches lies outside the primer length.")
    if isinstance(strict_three_prime_bases, bool) or not isinstance(strict_three_prime_bases, int):
        raise ConfigurationError("strict_three_prime_bases must be a non-negative integer.")
    if not 0 <= strict_three_prime_bases <= len(primer_symbols):
        raise ConfigurationError("strict_three_prime_bases lies outside the primer length.")
    validate_positive_int(max_hits, "max_hits", maximum=1_000_000)
    validate_positive_int(max_comparison_cells, "max_comparison_cells", maximum=1_000_000_000)
    strands = _coerce_strands(strand)
    circular = template.topology is Topology.CIRCULAR
    if len(primer_symbols) > len(template_symbols):
        window_count = 0
    else:
        window_count = (
            len(template_symbols) if circular else len(template_symbols) - len(primer_symbols) + 1
        )
    estimated_cells = window_count * len(primer_symbols) * len(strands)
    if estimated_cells > max_comparison_cells:
        raise ConfigurationError(
            "Primer scan exceeds max_comparison_cells.",
            code="PRIMER_MATCH_LIMIT_EXCEEDED",
            context={
                "estimated_cells": estimated_cells,
                "max_comparison_cells": max_comparison_cells,
            },
        )
    hits: list[PrimerBindingHit] = []
    truncated = False
    for resolved_strand in strands:
        for start in range(window_count):
            window = (
                circular_slice(template_symbols, start, len(primer_symbols))
                if circular
                else template_symbols[start : start + len(primer_symbols)]
            )
            mismatches = _comparison_positions(
                primer_symbols,
                window,
                resolved_strand,
                allow_iupac=allow_iupac,
            )
            three_prime = sum(
                position >= len(primer_symbols) - strict_three_prime_bases
                for position in mismatches
            )
            if len(mismatches) > max_mismatches or three_prime:
                continue
            if len(hits) >= max_hits:
                truncated = True
                break
            raw_end = start + len(primer_symbols)
            hits.append(
                PrimerBindingHit(
                    primer_name=primer_name,
                    strand=resolved_strand,
                    start=start,
                    end=(
                        raw_end - len(template_symbols)
                        if circular and raw_end > len(template_symbols)
                        else raw_end
                    ),
                    matched_template=window,
                    mismatch_positions_5to3=mismatches,
                    three_prime_mismatch_count=three_prime,
                    wraps_origin=circular and raw_end > len(template_symbols),
                )
            )
        if truncated:
            break
    return PrimerMatchResult(
        primer_sequence=primer_symbols,
        template_length=len(template_symbols),
        template_topology=template.topology.value,
        hits=tuple(hits),
        truncated=truncated,
        method="bounded_ungapped_hamming_primer_match",
        algorithm_version="dnakit-primer-match-v1",
        parameters=freeze_parameters(
            {
                "strand": str(strand),
                "max_mismatches": max_mismatches,
                "strict_three_prime_bases": strict_three_prime_bases,
                "allow_indels": False,
                "allow_iupac": allow_iupac,
                "max_hits": max_hits,
                "max_template_length": max_template_length,
                "estimated_comparison_cells": estimated_cells,
                "max_comparison_cells": max_comparison_cells,
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


def simulate_pcr(
    template: DNASequence,
    forward_primer: DNASequence,
    reverse_primer: DNASequence,
    *,
    max_mismatches: int = 0,
    strict_three_prime_bases: int = 1,
    min_product_length: int = 1,
    max_product_length: int = 100_000,
    allow_iupac: bool = False,
    max_binding_hits: int = 10_000,
    max_pair_checks: int = 1_000_000,
    max_products: int = 10_000,
    max_total_product_bases: int = 100_000_000,
    max_template_length: int = 10_000_000,
) -> PCRResult:
    """Predict inward-facing linear or single-traversal circular amplicons."""

    template_symbols = require_sequence(
        template,
        operation="PCR simulation",
        max_length=max_template_length,
        canonical=not allow_iupac,
    )
    forward_symbols = require_sequence(
        forward_primer,
        operation="PCR simulation",
        max_length=1_000,
        canonical=not allow_iupac,
        allow_circular=False,
    )
    reverse_symbols = require_sequence(
        reverse_primer,
        operation="PCR simulation",
        max_length=1_000,
        canonical=not allow_iupac,
        allow_circular=False,
    )
    validate_positive_int(min_product_length, "min_product_length")
    validate_positive_int(max_product_length, "max_product_length", maximum=100_000_000)
    if min_product_length > max_product_length:
        raise ConfigurationError("min_product_length cannot exceed max_product_length.")
    validate_positive_int(max_pair_checks, "max_pair_checks", maximum=100_000_000)
    validate_positive_int(max_products, "max_products", maximum=1_000_000)
    validate_positive_int(
        max_total_product_bases,
        "max_total_product_bases",
        maximum=10_000_000_000,
    )
    forward_matches = match_primer(
        forward_primer,
        template,
        primer_name="forward",
        strand=Strand.FORWARD,
        max_mismatches=max_mismatches,
        strict_three_prime_bases=strict_three_prime_bases,
        allow_iupac=allow_iupac,
        max_hits=max_binding_hits,
        max_template_length=max_template_length,
    )
    reverse_matches = match_primer(
        reverse_primer,
        template,
        primer_name="reverse",
        strand=Strand.REVERSE,
        max_mismatches=max_mismatches,
        strict_three_prime_bases=strict_three_prime_bases,
        allow_iupac=allow_iupac,
        max_hits=max_binding_hits,
        max_template_length=max_template_length,
    )
    if forward_matches.truncated or reverse_matches.truncated:
        raise ConfigurationError(
            "Primer binding hits were truncated; increase max_binding_hits before PCR pairing.",
            code="PCR_BINDING_HITS_TRUNCATED",
        )
    pair_checks = len(forward_matches.hits) * len(reverse_matches.hits)
    if pair_checks > max_pair_checks:
        raise ConfigurationError(
            "PCR primer pairing exceeds max_pair_checks.",
            code="PCR_PAIR_LIMIT_EXCEEDED",
            context={"pair_checks": pair_checks, "max_pair_checks": max_pair_checks},
        )
    circular = template.topology is Topology.CIRCULAR
    reverse_product_end = reverse_complement_text(reverse_symbols)
    products: list[Amplicon] = []
    total_product_bases = 0
    truncated = False
    for forward_hit in forward_matches.hits:
        for reverse_hit in reverse_matches.hits:
            if circular:
                product_length = (
                    reverse_hit.start + len(reverse_symbols) - forward_hit.start
                ) % len(template_symbols)
                if product_length == 0:
                    product_length = len(template_symbols)
                if product_length < len(forward_symbols) + len(reverse_symbols):
                    continue
                interior_length = product_length - len(forward_symbols) - len(reverse_symbols)
                interior = circular_slice(
                    template_symbols,
                    forward_hit.start + len(forward_symbols),
                    interior_length,
                )
                template_end = (forward_hit.start + product_length) % len(template_symbols)
                wraps = forward_hit.start + product_length > len(template_symbols)
            else:
                if forward_hit.end > reverse_hit.start:
                    continue
                product_length = reverse_hit.end - forward_hit.start
                interior = template_symbols[forward_hit.end : reverse_hit.start]
                template_end = reverse_hit.end
                wraps = False
            if not min_product_length <= product_length <= max_product_length:
                continue
            if len(products) >= max_products:
                truncated = True
                break
            if total_product_bases + product_length > max_total_product_bases:
                truncated = True
                break
            product_symbols = forward_symbols + interior + reverse_product_end
            alphabet = (
                DNAAlphabet.IUPAC if set(product_symbols) - set("ACGT") else DNAAlphabet.STRICT
            )
            products.append(
                Amplicon(
                    sequence=DNASequence(
                        product_symbols,
                        alphabet=alphabet,
                        topology=Topology.LINEAR,
                        strandedness=template.strandedness,
                    ),
                    template_start=forward_hit.start,
                    template_end=template_end,
                    wraps_origin=wraps,
                    forward_binding=forward_hit,
                    reverse_binding=reverse_hit,
                )
            )
            total_product_bases += product_length
        if truncated:
            break
    return PCRResult(
        template_length=len(template_symbols),
        template_topology=template.topology.value,
        amplicons=tuple(products),
        truncated=truncated,
        method="inward_facing_ungapped_primer_pairing",
        algorithm_version="dnakit-pcr-simulation-v1",
        parameters=freeze_parameters(
            {
                "max_mismatches": max_mismatches,
                "strict_three_prime_bases": strict_three_prime_bases,
                "min_product_length": min_product_length,
                "max_product_length": max_product_length,
                "allow_iupac": allow_iupac,
                "allow_indels": False,
                "circular_traversals": 1 if circular else 0,
                "pair_checks": pair_checks,
                "max_pair_checks": max_pair_checks,
                "max_products": max_products,
                "total_product_bases": total_product_bases,
                "max_total_product_bases": max_total_product_bases,
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


def _conditional_structure(capability: str) -> ConditionalAnalysis:
    return ConditionalAnalysis(
        capability=capability,
        status="conditional-unavailable",
        available=False,
        backend_required="Primer3-compatible user-installed adapter",
        automatic_probe=False,
        automatic_execution=False,
        reason=(
            "The native module does not claim a validated secondary-structure model; "
            "prepare a backend request and execute it explicitly outside this function."
        ),
    )


def _completed_structure(
    result: object,
    *,
    capability: str,
    expected_sequences: tuple[str, ...],
    conditions: ThermodynamicConditions,
    max_loop: int,
    output_structure: bool,
) -> tuple[ConditionalAnalysis, Primer3ThermodynamicResult]:
    validated = validate_primer3_result(
        result,
        capability=capability,
        sequences_5to3=expected_sequences,
        conditions=conditions,
        max_loop=max_loop,
        output_structure=output_structure,
        error_code="MISMATCHED_PRIMER_STRUCTURE_RESULT",
    )
    return (
        ConditionalAnalysis(
            capability=capability,
            status="execution-complete",
            available=True,
            backend_required=validated.backend.name,
            automatic_probe=False,
            automatic_execution=False,
            execution_performed=True,
            reason="The user explicitly supplied and executed a structure adapter.",
            structure_found=validated.structure_found,
            tm_celsius=validated.tm_celsius,
            delta_g_kcal_per_mol=validated.delta_g_kcal_per_mol,
            delta_h_kcal_per_mol=validated.delta_h_kcal_per_mol,
            delta_s_cal_per_k_mol=validated.delta_s_cal_per_k_mol,
            ascii_structure=validated.ascii_structure,
            backend_name=validated.backend.name,
            backend_version=validated.backend.version,
        ),
        validated,
    )


def _validate_consistent_adapter_results(
    results: Iterable[Primer3ThermodynamicResult],
) -> None:
    resolved = tuple(results)
    if not resolved:
        return
    first = resolved[0]
    if any(
        result.backend != first.backend or result.provenance != first.provenance
        for result in resolved[1:]
    ):
        raise BackendExecutionError(
            "Primer adapter returned inconsistent backend or provenance metadata.",
            code="INCONSISTENT_PRIMER_ADAPTER_RESULTS",
        )


def primer_properties(
    primer: DNASequence,
    *,
    paired_primer: DNASequence | None = None,
    tm_method: str = "nearest_neighbor",
    conditions: ThermodynamicConditions | None = None,
    nn_config: NearestNeighborConfig | None = None,
    structure_adapter: PrimerStructureAdapter | None = None,
    max_loop: int = 30,
    output_structure: bool = False,
) -> PrimerPropertiesResult:
    """Calculate GC/Tm and explicitly execute an optional structure adapter."""

    symbols = require_sequence(
        primer,
        operation="primer property calculation",
        max_length=60,
        canonical=True,
        allow_circular=False,
    )
    paired_symbols: str | None = None
    paired_tm: float | None = None
    paired_gc: float | None = None
    if paired_primer is not None:
        paired_symbols = require_sequence(
            paired_primer,
            operation="primer pair property calculation",
            max_length=60,
            canonical=True,
            allow_circular=False,
        )
    if tm_method not in ("wallace", "nearest_neighbor", "primer3-cli"):
        raise ConfigurationError("tm_method must be wallace, nearest_neighbor, or primer3-cli.")
    if isinstance(max_loop, bool) or not isinstance(max_loop, int) or not 1 <= max_loop <= 30:
        raise ConfigurationError("max_loop must be an integer in [1, 30].")
    if not isinstance(output_structure, bool):
        raise ConfigurationError("output_structure must be boolean.")
    if structure_adapter is not None and any(
        not callable(getattr(structure_adapter, name, None))
        for name in ("hairpin", "self_dimer", "heterodimer")
    ):
        raise ConfigurationError(
            "structure_adapter must implement hairpin, self_dimer, and heterodimer.",
            code="INVALID_PRIMER_STRUCTURE_ADAPTER",
        )
    resolved_conditions = ThermodynamicConditions() if conditions is None else conditions
    if not isinstance(resolved_conditions, ThermodynamicConditions):
        raise ConfigurationError(
            "conditions must be ThermodynamicConditions or None.",
            code="INVALID_PRIMER_CONDITIONS",
        )
    tm: Primer3ThermodynamicResult | MeltingTemperatureResult
    if tm_method == "primer3-cli":
        if structure_adapter is None or not callable(getattr(structure_adapter, "tm", None)):
            raise ConfigurationError(
                "tm_method='primer3-cli' requires an explicit adapter with a callable tm method.",
                code="PRIMER3_TM_ADAPTER_REQUIRED",
            )
        if nn_config is not None:
            raise ConfigurationError(
                "nn_config is not applicable to tm_method='primer3-cli'.",
                code="UNUSED_PRIMER_NATIVE_TM_CONFIG",
            )
        tm_adapter = cast(_PrimerTmAdapter, structure_adapter)
        tm = validate_primer3_result(
            tm_adapter.tm(primer, conditions=resolved_conditions),
            capability="tm",
            sequences_5to3=(symbols,),
            conditions=resolved_conditions,
            max_loop=None,
            output_structure=False,
            error_code="MISMATCHED_PRIMER_TM_RESULT",
        )
    else:
        tm = melting_temperature(
            primer,
            method=tm_method,  # type: ignore[arg-type]
            conditions=resolved_conditions,
            config=nn_config,
        )
    if paired_primer is not None and paired_symbols is not None:
        paired_result: Primer3ThermodynamicResult | MeltingTemperatureResult
        if tm_method == "primer3-cli":
            assert structure_adapter is not None
            tm_adapter = cast(_PrimerTmAdapter, structure_adapter)
            paired_result = validate_primer3_result(
                tm_adapter.tm(paired_primer, conditions=resolved_conditions),
                capability="tm",
                sequences_5to3=(paired_symbols,),
                conditions=resolved_conditions,
                max_loop=None,
                output_structure=False,
                error_code="MISMATCHED_PRIMER_TM_RESULT",
            )
        else:
            paired_result = melting_temperature(
                paired_primer,
                method=tm_method,  # type: ignore[arg-type]
                conditions=resolved_conditions,
                config=nn_config,
            )
        paired_tm = paired_result.tm_celsius
        paired_gc = sum(base in "GC" for base in paired_symbols) / len(paired_symbols)
    issues: tuple[Issue, ...]
    if structure_adapter is None:
        hairpin = _conditional_structure("hairpin")
        self_dimer = _conditional_structure("self_dimer")
        paired_hairpin = _conditional_structure("hairpin") if paired_primer is not None else None
        paired_self_dimer = (
            _conditional_structure("self_dimer") if paired_primer is not None else None
        )
        heterodimer = _conditional_structure("heterodimer") if paired_primer is not None else None
        issues = (
            Issue(
                "PRIMER_STRUCTURE_BACKEND_REQUIRED",
                IssueSeverity.INFO,
                "Hairpin and dimer values require an explicitly supplied validated adapter.",
            ),
        )
        provenance = native_provenance()
    else:
        primary_hairpin_raw = structure_adapter.hairpin(
            primer,
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        primary_self_raw = structure_adapter.self_dimer(
            primer,
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        hairpin, primary_hairpin_result = _completed_structure(
            primary_hairpin_raw,
            capability="hairpin",
            expected_sequences=(symbols,),
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        self_dimer, primary_self_result = _completed_structure(
            primary_self_raw,
            capability="self_dimer",
            expected_sequences=(symbols,),
            conditions=resolved_conditions,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        paired_hairpin = None
        paired_self_dimer = None
        heterodimer = None
        completed_results = [primary_hairpin_result, primary_self_result]
        if paired_primer is not None and paired_symbols is not None:
            paired_hairpin, paired_hairpin_result = _completed_structure(
                structure_adapter.hairpin(
                    paired_primer,
                    conditions=resolved_conditions,
                    max_loop=max_loop,
                    output_structure=output_structure,
                ),
                capability="hairpin",
                expected_sequences=(paired_symbols,),
                conditions=resolved_conditions,
                max_loop=max_loop,
                output_structure=output_structure,
            )
            paired_self_dimer, paired_self_result = _completed_structure(
                structure_adapter.self_dimer(
                    paired_primer,
                    conditions=resolved_conditions,
                    max_loop=max_loop,
                    output_structure=output_structure,
                ),
                capability="self_dimer",
                expected_sequences=(paired_symbols,),
                conditions=resolved_conditions,
                max_loop=max_loop,
                output_structure=output_structure,
            )
            heterodimer, heterodimer_result = _completed_structure(
                structure_adapter.heterodimer(
                    primer,
                    paired_primer,
                    conditions=resolved_conditions,
                    max_loop=max_loop,
                    output_structure=output_structure,
                ),
                capability="heterodimer",
                expected_sequences=(symbols, paired_symbols),
                conditions=resolved_conditions,
                max_loop=max_loop,
                output_structure=output_structure,
            )
            completed_results.extend(
                (paired_hairpin_result, paired_self_result, heterodimer_result)
            )
        if tm_method == "primer3-cli":
            if not isinstance(tm, Primer3ThermodynamicResult):
                raise AssertionError("Primer3 Tm branch must produce a Primer3 result.")
            completed_results.append(tm)
            if paired_primer is not None and paired_symbols is not None:
                if not isinstance(paired_result, Primer3ThermodynamicResult):
                    raise AssertionError("Primer3 paired Tm must produce a Primer3 result.")
                completed_results.append(paired_result)
        _validate_consistent_adapter_results(completed_results)
        issues = ()
        provenance = primary_hairpin_raw.provenance
    return PrimerPropertiesResult(
        primer_sequence=symbols,
        gc_fraction=sum(base in "GC" for base in symbols) / len(symbols),
        tm_celsius=tm.tm_celsius,
        tm_method=tm.method,
        hairpin=hairpin,
        self_dimer=self_dimer,
        heterodimer=heterodimer,
        paired_primer_sequence=paired_symbols,
        paired_gc_fraction=paired_gc,
        paired_tm_celsius=paired_tm,
        paired_hairpin=paired_hairpin,
        paired_self_dimer=paired_self_dimer,
        method="native_gc_tm_with_explicit_structure_adapter",
        algorithm_version="dnakit-primer-properties-v2",
        parameters=freeze_parameters(
            {
                "tm_method": tm_method,
                "conditions": tm.conditions.to_dict()
                if hasattr(tm.conditions, "to_dict")
                else {
                    "temperature_celsius": tm.conditions.temperature_celsius,
                    "sodium_molar": tm.conditions.sodium_molar,
                    "potassium_molar": tm.conditions.potassium_molar,
                    "magnesium_molar": tm.conditions.magnesium_molar,
                    "dntp_molar": tm.conditions.dntp_molar,
                    "strand_concentration_molar": tm.conditions.strand_concentration_molar,
                    "dmso_percent": tm.conditions.dmso_percent,
                    "dmso_factor_celsius_per_percent": (
                        tm.conditions.dmso_factor_celsius_per_percent
                    ),
                    "formamide_molar": tm.conditions.formamide_molar,
                    "salt_model": tm.conditions.salt_model,
                },
                "hairpin_backend_executed": structure_adapter is not None,
                "dimer_backend_executed": structure_adapter is not None,
                "structure_adapter_supplied": structure_adapter is not None,
                "automatic_backend_probe": False,
                "automatic_backend_execution": False,
                "max_loop": max_loop,
                "output_structure": output_structure,
            }
        ),
        provenance=provenance,
        issues=issues,
    )


def prepare_primer_design(
    template: DNASequence,
    *,
    target_start: int,
    target_end: int,
    primer_length_range: tuple[int, int] = (18, 25),
    tm_range_celsius: tuple[float, float] = (57.0, 63.0),
    gc_range: tuple[float, float] = (0.4, 0.6),
    product_length_range: tuple[int, int] = (100, 1_000),
    excluded_regions: object = (),
    candidate_count: int = 5,
    backend_name: str = "primer3",
    conditions: ThermodynamicConditions | None = None,
    max_excluded_regions: int = 10_000,
    max_template_length: int = 10_000_000,
) -> PrimerDesignInterfaceResult:
    """Validate and serialize a backend-neutral primer-design request only."""

    symbols = require_sequence(
        template,
        operation="primer design request",
        max_length=max_template_length,
        canonical=True,
        allow_circular=False,
    )
    if (
        isinstance(target_start, bool)
        or not isinstance(target_start, int)
        or isinstance(target_end, bool)
        or not isinstance(target_end, int)
        or not 0 <= target_start < target_end <= len(symbols)
    ):
        raise ConfigurationError("Target must be a valid 0-based half-open template interval.")

    def pair(value: object, name: str) -> tuple[object, object]:
        if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
            raise ConfigurationError(f"{name} must contain exactly two values.")
        items = tuple(islice(iter(value), 3))
        if len(items) != 2:
            raise ConfigurationError(f"{name} must contain exactly two values.")
        return items[0], items[1]

    raw_primer_min, raw_primer_max = pair(primer_length_range, "primer_length_range")
    raw_product_min, raw_product_max = pair(product_length_range, "product_length_range")
    primer_min = validate_positive_int(raw_primer_min, "minimum primer length", maximum=1_000)
    primer_max = validate_positive_int(raw_primer_max, "maximum primer length", maximum=1_000)
    product_min = validate_positive_int(
        raw_product_min, "minimum product length", maximum=100_000_000
    )
    product_max = validate_positive_int(
        raw_product_max, "maximum product length", maximum=100_000_000
    )
    if primer_min > primer_max or product_min > product_max:
        raise ConfigurationError("Length range lower bounds cannot exceed upper bounds.")
    tm_min, tm_max = pair(tm_range_celsius, "tm_range_celsius")
    if (
        isinstance(tm_min, bool)
        or not isinstance(tm_min, (int, float))
        or isinstance(tm_max, bool)
        or not isinstance(tm_max, (int, float))
    ):
        raise ConfigurationError("Tm range bounds must be numeric.")
    resolved_tm_min, resolved_tm_max = float(tm_min), float(tm_max)
    if not -100.0 <= resolved_tm_min <= resolved_tm_max <= 200.0:
        raise ConfigurationError("Tm range is invalid or outside the request boundary.")
    raw_gc_min, raw_gc_max = pair(gc_range, "gc_range")
    gc_min, gc_max = (
        finite_fraction(raw_gc_min, "minimum GC"),
        finite_fraction(raw_gc_max, "maximum GC"),
    )
    if gc_min > gc_max:
        raise ConfigurationError("Minimum GC cannot exceed maximum GC.")
    validate_positive_int(candidate_count, "candidate_count", maximum=10_000)
    validate_positive_int(max_excluded_regions, "max_excluded_regions", maximum=1_000_000)
    if not isinstance(backend_name, str) or not backend_name.strip():
        raise ConfigurationError("backend_name must be non-empty.")
    if isinstance(excluded_regions, (str, bytes)) or not isinstance(excluded_regions, Iterable):
        raise ConfigurationError("excluded_regions must be an iterable of intervals.")
    region_values = tuple(islice(iter(excluded_regions), max_excluded_regions + 1))
    if len(region_values) > max_excluded_regions:
        raise ConfigurationError("excluded_regions exceeds max_excluded_regions.")
    regions: list[tuple[int, int]] = []
    for region in region_values:
        if (
            not isinstance(region, tuple)
            or len(region) != 2
            or any(isinstance(item, bool) or not isinstance(item, int) for item in region)
            or not 0 <= region[0] < region[1] <= len(symbols)
        ):
            raise ConfigurationError("Each excluded region must be a valid template interval.")
        regions.append(region)
    resolved_conditions = ThermodynamicConditions() if conditions is None else conditions
    if not isinstance(resolved_conditions, ThermodynamicConditions):
        raise ConfigurationError(
            "conditions must be ThermodynamicConditions or None.",
            code="INVALID_PRIMER_CONDITIONS",
        )
    condition_values = freeze_parameters(
        {
            "temperature_celsius": resolved_conditions.temperature_celsius,
            "sodium_molar": resolved_conditions.sodium_molar,
            "potassium_molar": resolved_conditions.potassium_molar,
            "magnesium_molar": resolved_conditions.magnesium_molar,
            "dntp_molar": resolved_conditions.dntp_molar,
            "strand_concentration_molar": resolved_conditions.strand_concentration_molar,
            "dmso_percent": resolved_conditions.dmso_percent,
            "dmso_factor_celsius_per_percent": (
                resolved_conditions.dmso_factor_celsius_per_percent
            ),
            "formamide_molar": resolved_conditions.formamide_molar,
            "salt_model": resolved_conditions.salt_model,
        }
    )
    request = PrimerDesignRequest(
        template=template,
        target_start=target_start,
        target_end=target_end,
        primer_length_range=(primer_min, primer_max),
        tm_range_celsius=(resolved_tm_min, resolved_tm_max),
        gc_range=(gc_min, gc_max),
        product_length_range=(product_min, product_max),
        excluded_regions=tuple(regions),
        candidate_count=candidate_count,
        thermodynamic_conditions=condition_values,
    )
    return PrimerDesignInterfaceResult(
        request=request,
        backend_name=backend_name,
        status="request-ready-execution-not-performed",
        execution_performed=False,
        candidates=(),
        reason=(
            "DNAKit validates a neutral request but does not install, probe, or execute "
            "an external primer-design backend in this interface."
        ),
        method="backend_neutral_primer_design_request",
        algorithm_version="dnakit-primer-design-interface-v1",
        parameters=freeze_parameters(
            {
                "backend_name": backend_name,
                "automatic_install": False,
                "automatic_probe": False,
                "automatic_execution": False,
                "max_excluded_regions": max_excluded_regions,
                "max_template_length": max_template_length,
            }
        ),
        provenance=adapter_provenance(reference_name="Primer3-compatible interface"),
        issues=(),
    )


__all__ = [
    "match_primer",
    "prepare_primer_design",
    "primer_properties",
    "simulate_pcr",
]
