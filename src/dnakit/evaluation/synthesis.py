"""Rule-based synthesis-risk screening with no experimental success claim."""

from __future__ import annotations

import hashlib
import math

from dnakit.core import IssueSeverity, Topology
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.descriptors import gc_at_content, homopolymer_runs, window_descriptors
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.patterns import find_inverted_repeats, find_tandem_repeats

from ._shared import EvaluationInput, aggregate_numeric, as_float, issue, materialize_input, report
from .config import SynthesisRiskConfig
from .results import EvaluationEntry, EvaluationReport

_DISCLAIMER = (
    "DNAKit transparent sequence-rule screen only; this score is not a vendor acceptance "
    "criterion, thermodynamic structure prediction, or experimental synthesis-success probability."
)


def _window_count(
    sequence_length: int,
    *,
    window_size: int,
    step: int,
    include_partial: bool,
) -> int:
    """Return the exact number of rows produced by ``window_descriptors``."""

    last_start = sequence_length if include_partial else sequence_length - window_size + 1
    stop = max(0, last_start)
    return (stop + step - 1) // step if stop else 0


def _risk_entry(item: object, index: int, config: SynthesisRiskConfig) -> EvaluationEntry:
    from ._shared import InputItem

    if not isinstance(item, InputItem):
        raise TypeError("Internal synthesis input must be InputItem.")
    sequence = item.sequence
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Synthesis-risk rules require gap-free input.",
            code="SYNTHESIS_RISK_GAP_UNSUPPORTED",
            context={"subject_id": item.subject_id},
        )
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Synthesis-risk rules currently require linear topology.",
            code="SYNTHESIS_RISK_LINEAR_REQUIRED",
            context={"subject_id": item.subject_id},
        )
    canonical_symbols = sequence.symbols
    canonical_sequence_sha256 = hashlib.sha256(canonical_symbols.encode("ascii")).hexdigest()
    window_count = _window_count(
        sequence.symbol_length,
        window_size=config.window_size,
        step=config.window_step,
        include_partial=True,
    )
    if window_count > config.max_windows_per_sequence:
        raise ConfigurationError(
            "Synthesis-risk calculation exceeds max_windows_per_sequence.",
            code="SYNTHESIS_RISK_WINDOW_LIMIT",
            context={
                "window_count": window_count,
                "max_windows_per_sequence": config.max_windows_per_sequence,
                "sequence_length": sequence.symbol_length,
                "window_size": config.window_size,
                "window_step": config.window_step,
                "include_partial": True,
            },
        )
    global_gc = gc_at_content(sequence, ambiguity_policy="ignore").gc_fraction
    windows = window_descriptors(
        sequence,
        ("gc",),
        window_size=config.window_size,
        step=config.window_step,
        include_partial=True,
        ambiguity_policy="ignore",
    )
    risky_windows = tuple(
        {
            "start": window.symbol_start,
            "end": window.symbol_end,
            "gc_fraction": window.values["gc_fraction"],
        }
        for window in windows.windows
        if isinstance(window.values["gc_fraction"], (int, float))
        and not config.local_gc_min <= float(window.values["gc_fraction"]) <= config.local_gc_max
    )
    homopolymer = homopolymer_runs(
        sequence,
        min_run_length=config.homopolymer_threshold,
        ambiguity_policy="ignore",
    )
    tandem = find_tandem_repeats(
        sequence,
        min_unit_length=config.tandem_min_unit,
        max_unit_length=config.tandem_max_unit,
        min_repeats=config.tandem_min_repeats,
        max_comparisons=config.max_pattern_comparisons_per_sequence,
        max_comparison_cells=config.max_pattern_comparisons_per_sequence,
        max_matches=config.max_matches_per_sequence,
    )
    inverted = find_inverted_repeats(
        sequence,
        min_arm_length=config.inverted_min_arm,
        max_arm_length=config.inverted_max_arm,
        max_loop_length=config.inverted_max_loop,
        max_comparisons=config.max_pattern_comparisons_per_sequence,
        max_comparison_cells=config.max_pattern_comparisons_per_sequence,
        max_matches=config.max_matches_per_sequence,
    )
    components = {
        "global_gc": float(
            global_gc is None or not config.global_gc_min <= global_gc <= config.global_gc_max
        ),
        "local_gc": min(1.0, len(risky_windows) / max(1, len(windows.windows))),
        "homopolymer": min(
            1.0,
            max(0, homopolymer.longest_length - config.homopolymer_threshold + 1)
            / config.homopolymer_threshold,
        ),
        "tandem_repeat": min(1.0, len(tandem.hits) / 5.0),
        "inverted_repeat_structure_proxy": min(1.0, len(inverted.hits) / 5.0),
    }
    score = math.fsum(components.values()) / len(components)
    findings = []
    if any(value > 0 for value in components.values()):
        findings.append(
            issue(
                "EVAL_SYNTHESIS_RULE_FINDING",
                IssueSeverity.WARNING,
                "One or more transparent sequence-risk rules were triggered.",
                triggered=tuple(name for name, value in components.items() if value > 0),
            )
        )
    if sequence.ambiguity_count:
        findings.append(
            issue(
                "EVAL_SYNTHESIS_AMBIGUITY_IGNORED",
                IssueSeverity.WARNING,
                "Ambiguous IUPAC symbols were excluded from composition-based calculations.",
                ambiguity_count=sequence.ambiguity_count,
            )
        )
    return EvaluationEntry(
        item.subject_id,
        index,
        FrozenDict(
            {
                "score": score,
                "risk_score": score,
                "risk_level": "high" if score >= 0.5 else ("medium" if score >= 0.2 else "low"),
                "canonical_sequence_sha256": canonical_sequence_sha256,
                "canonical_sequence_length": sequence.symbol_length,
                "components": components,
                "global_gc_fraction": global_gc,
                "risky_gc_windows": risky_windows,
                "window_count": window_count,
                "longest_homopolymer": homopolymer.longest_length,
                "homopolymer_runs": to_json_compatible(homopolymer.runs),
                "tandem_repeat_hits": to_json_compatible(tandem.hits),
                "inverted_repeat_hits": to_json_compatible(inverted.hits),
                "tandem_scan_truncated": tandem.truncated,
                "inverted_repeat_scan_truncated": inverted.truncated,
                "tandem_comparisons": tandem.parameters["comparisons"],
                "inverted_repeat_comparisons": inverted.parameters["comparisons"],
                "disclaimer": _DISCLAIMER,
            }
        ),
        tuple(findings),
    )


def evaluate_synthesis_risk(
    value: EvaluationInput,
    *,
    config: SynthesisRiskConfig | None = None,
) -> EvaluationReport:
    """Screen transparent GC, run, repeat, and inverted-repeat proxy rules."""

    resolved = SynthesisRiskConfig() if config is None else config
    if not isinstance(resolved, SynthesisRiskConfig):
        raise TypeError("config must be SynthesisRiskConfig or None.")
    items = materialize_input(value, limits=resolved.limits)
    entries = tuple(_risk_entry(item, index, resolved) for index, item in enumerate(items))
    return report(
        name="synthesis_risk",
        method="equal-weight-transparent-sequence-rule-screen",
        version="eval-synthesis-risk-v2",
        parameters={
            "global_gc_range": (resolved.global_gc_min, resolved.global_gc_max),
            "local_gc_range": (resolved.local_gc_min, resolved.local_gc_max),
            "window_size": resolved.window_size,
            "window_step": resolved.window_step,
            "partial_windows": True,
            "homopolymer_threshold": resolved.homopolymer_threshold,
            "tandem": {
                "min_unit": resolved.tandem_min_unit,
                "max_unit": resolved.tandem_max_unit,
                "min_repeats": resolved.tandem_min_repeats,
            },
            "inverted_repeat": {
                "min_arm": resolved.inverted_min_arm,
                "max_arm": resolved.inverted_max_arm,
                "max_loop": resolved.inverted_max_loop,
            },
            "aggregation": "arithmetic mean of five component risk values",
            "structure_method": "inverted-repeat count proxy; no folding backend",
            "sequence_binding": {
                "digest": "SHA-256",
                "payload": "DNASequence.symbols encoded as ASCII",
                "length": "DNASequence.symbol_length",
            },
            "vendor_rule_set": None,
            "disclaimer": _DISCLAIMER,
            "resource_limits": {
                "max_windows_per_sequence": resolved.max_windows_per_sequence,
                "max_pattern_comparisons_per_sequence": (
                    resolved.max_pattern_comparisons_per_sequence
                ),
                "max_matches_per_sequence": resolved.max_matches_per_sequence,
                "evaluation_limits": resolved.limits,
            },
        },
        metrics={
            "score": (
                math.fsum(as_float(entry.metrics, "score") for entry in entries) / len(entries)
                if entries
                else None
            ),
            "record_count": len(entries),
            "risk_score_summary": aggregate_numeric(entries, "risk_score"),
            "high_risk_count": sum(entry.metrics["risk_level"] == "high" for entry in entries),
            "medium_risk_count": sum(entry.metrics["risk_level"] == "medium" for entry in entries),
            "low_risk_count": sum(entry.metrics["risk_level"] == "low" for entry in entries),
            "disclaimer": _DISCLAIMER,
        },
        entries=entries,
    )


__all__ = ["evaluate_synthesis_risk"]
