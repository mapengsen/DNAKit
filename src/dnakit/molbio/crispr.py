"""CRISPR candidate scanning with bounded sequence-only off-target enumeration."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit.core import DNASequence, Issue, IssueSeverity, Strand
from dnakit.exceptions import ConfigurationError
from dnakit.patterns import PAMRule, scan_pam_candidates
from dnakit.patterns.results import GuideCandidate

from ._shared import freeze_parameters, native_provenance, require_sequence, validate_positive_int
from .primers import match_primer
from .results import CrisprScanResult, OffTargetHit


def scan_crispr_candidates(
    target: DNASequence,
    rule: str | PAMRule,
    *,
    guide_length: int | None = None,
    min_gc: float = 0.0,
    max_gc: float = 1.0,
    references: Mapping[str, DNASequence] | None = None,
    max_off_target_mismatches: int = 0,
    max_candidates: int = 10_000,
    max_off_targets: int = 100_000,
    max_reference_sequences: int = 10_000,
    max_reference_length: int = 10_000_000,
    max_total_comparison_cells: int = 500_000_000,
) -> CrisprScanResult:
    """Reuse PAM rules and optionally enumerate ungapped sequence matches.

    This function intentionally does not predict editing efficiency, bulges,
    chromatin effects, or biological off-target risk.
    """

    require_sequence(
        target,
        operation="CRISPR candidate scanning",
        max_length=max_reference_length,
    )
    validate_positive_int(max_candidates, "max_candidates", maximum=1_000_000)
    validate_positive_int(max_off_targets, "max_off_targets", maximum=10_000_000)
    validate_positive_int(max_reference_sequences, "max_reference_sequences", maximum=1_000_000)
    validate_positive_int(
        max_total_comparison_cells,
        "max_total_comparison_cells",
        maximum=10_000_000_000,
    )
    if (
        isinstance(max_off_target_mismatches, bool)
        or not isinstance(max_off_target_mismatches, int)
        or max_off_target_mismatches < 0
    ):
        raise ConfigurationError("max_off_target_mismatches must be non-negative.")
    candidate_scan = scan_pam_candidates(
        target,
        rule,
        guide_length=guide_length,
        min_gc=min_gc,
        max_gc=max_gc,
        max_matches=max_candidates,
        max_scan_length=max_reference_length,
    )
    resolved_references = references or {}
    if len(resolved_references) > max_reference_sequences:
        raise ConfigurationError("Reference collection exceeds max_reference_sequences.")
    reference_items = tuple(resolved_references.items())
    total_reference_length = 0
    for reference_id, sequence in reference_items:
        if not isinstance(reference_id, str) or not reference_id.strip():
            raise ConfigurationError("Reference identifiers must be non-empty strings.")
        total_reference_length += len(
            require_sequence(
                sequence,
                operation="CRISPR off-target enumeration",
                max_length=max_reference_length,
                canonical=True,
            )
        )
    estimated_cells = sum(
        total_reference_length * len(candidate.guide_sequence) * 2
        for candidate in candidate_scan.hits
    )
    if estimated_cells > max_total_comparison_cells:
        raise ConfigurationError(
            "CRISPR off-target enumeration exceeds max_total_comparison_cells.",
            code="CRISPR_OFF_TARGET_LIMIT_EXCEEDED",
            context={
                "estimated_cells": estimated_cells,
                "max_total_comparison_cells": max_total_comparison_cells,
            },
        )
    off_targets: list[OffTargetHit] = []
    off_target_truncated = False
    for candidate in candidate_scan.hits:
        if not isinstance(candidate, GuideCandidate):
            raise AssertionError("Unexpected PAM candidate result type.")
        guide = DNASequence(candidate.guide_sequence)
        for reference_id, reference in reference_items:
            matches = match_primer(
                guide,
                reference,
                primer_name="guide",
                strand=Strand.BOTH,
                max_mismatches=max_off_target_mismatches,
                max_hits=max_off_targets,
                max_template_length=max_reference_length,
                max_comparison_cells=max_total_comparison_cells,
            )
            for hit in matches.hits:
                if len(off_targets) >= max_off_targets:
                    off_target_truncated = True
                    break
                off_targets.append(
                    OffTargetHit(
                        guide_sequence=candidate.guide_sequence,
                        reference_id=reference_id,
                        strand=hit.strand,
                        start=hit.start,
                        end=hit.end,
                        mismatch_positions_5to3=hit.mismatch_positions_5to3,
                        wraps_origin=hit.wraps_origin,
                    )
                )
            if matches.truncated:
                off_target_truncated = True
            if off_target_truncated:
                break
        if off_target_truncated:
            break
    issues: list[Issue] = []
    if not reference_items:
        issues.append(
            Issue(
                "CRISPR_REFERENCE_NOT_PROVIDED",
                IssueSeverity.INFO,
                "Candidates were returned without sequence off-target enumeration.",
            )
        )
    issues.append(
        Issue(
            "CRISPR_BIOLOGICAL_RISK_NOT_PREDICTED",
            IssueSeverity.INFO,
            "Sequence matches are not editing-efficiency or biological off-target predictions.",
        )
    )
    return CrisprScanResult(
        candidates=tuple(candidate_scan.hits),
        off_targets=tuple(off_targets),
        candidate_truncated=candidate_scan.truncated,
        off_target_truncated=off_target_truncated,
        efficiency_prediction_performed=False,
        method="pam_scan_plus_bounded_ungapped_sequence_matching",
        algorithm_version="dnakit-crispr-candidates-v1",
        parameters=freeze_parameters(
            {
                "guide_length": guide_length,
                "min_gc": min_gc,
                "max_gc": max_gc,
                "reference_count": len(reference_items),
                "max_off_target_mismatches": max_off_target_mismatches,
                "allow_indels_or_bulges": False,
                "efficiency_prediction": False,
                "estimated_comparison_cells": estimated_cells,
                "max_total_comparison_cells": max_total_comparison_cells,
                "max_candidates": max_candidates,
                "max_off_targets": max_off_targets,
            }
        ),
        provenance=native_provenance(
            reimplementation=True,
            reference_name="DNAKit common PAM rule definitions",
            reference_version="common-nucleases-v1",
        ),
        issues=tuple(issues),
    )


__all__ = ["scan_crispr_candidates"]
