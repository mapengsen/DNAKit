"""Known-cis point-edit annotation and explicit transcript-model prediction (pack P32).

These functions use the caller's local GTF/FASTA reference. They do not infer phase:
``variants`` must describe alleles already known to be in cis. Results preserve
unresolved outcomes from varcode instead of manufacturing a successful annotation.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab.varcode_gap_functions import (
    _VARCODE_EXTRA,
    _call,
    _candidate_record,
    _document,
    _effect_record,
    _genome,
    _int,
    _mutant_record,
    _transcript,
    _variant,
    _variant_record,
)


def _inputs(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    transcript_id: str,
    cache_dir: str | None,
) -> tuple[Any, Any, list[Any], dict[str, Any]]:
    if isinstance(variants, (str, bytes)) or not isinstance(variants, Sequence) or not variants:
        raise ConfigurationError("variants must be a non-empty list of known-cis variant records.")
    if len(variants) > 1000:
        raise ConfigurationError("variants may contain at most 1000 records.")
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    return (
        genome,
        _transcript(genome, transcript_id),
        [_variant(genome, item, f"variants[{index}]") for index, item in enumerate(variants)],
        info,
    )


def varcode_annotate_cis_haplotype(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    transcript_id: str,
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate the joint effect of known-cis point edits with FastEffectAnnotator.annotate_haplotype.

    At least two variants are required; each has contig, 1-based start, ref and alt. Structural variants and
    unsupported transcript models raise ConfigurationError, matching varcode's
    explicit NotImplemented outcome. No germline or RNA evidence is inferred.
    """
    if len(variants) < 2:
        raise ConfigurationError("Haplotype annotation requires at least two known-cis variants.")
    _, transcript, items, info = _inputs(gtf, genome_fasta, variants, transcript_id, cache_dir)
    annotator = require("varcode.annotators.fast", _VARCODE_EXTRA).FastEffectAnnotator()
    effect = _call(annotator.annotate_haplotype, items, transcript)
    if effect is NotImplemented:
        raise ConfigurationError(
            "varcode fast annotator does not support this haplotype/transcript."
        )
    mutant = getattr(effect, "mutant_transcript", None)
    return _document(
        "annotate_cis_haplotype",
        variants=[_variant_record(item) for item in items],
        transcript_id=transcript_id,
        effect=_effect_record(effect),
        mutant_transcript=None if mutant is None else _mutant_record(mutant),
        cache_directory=info["cache_directory"],
    )


def varcode_predict_transcript_model(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    transcript_id: str,
    *,
    germline_variants: Sequence[Mapping[str, Any]] = (),
    max_hypotheses: int = 64,
    use_genomic_sequence: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict transcript-model effects with splice alternatives via predict_transcript_model_effect.

    Known-cis variants are evaluated jointly; germline_variants retain upstream
    phase uncertainty. Genomic sequence enables intron-retention/cryptic-splice
    realisation. HypothesisLimit and Unresolved remain explicit effect classes.
    """
    genome, transcript, items, info = _inputs(gtf, genome_fasta, variants, transcript_id, cache_dir)
    if isinstance(germline_variants, (str, bytes)) or not isinstance(germline_variants, Sequence):
        raise ConfigurationError("germline_variants must be a list of variant records.")
    if not isinstance(use_genomic_sequence, bool):
        raise ConfigurationError("use_genomic_sequence must be a boolean.")
    varcode = require("varcode", _VARCODE_EXTRA)
    germline = [_variant(genome, item, "germline_variant") for item in germline_variants]
    cap = _int(max_hypotheses, "max_hypotheses", minimum=1, maximum=10000)
    provider = (
        (lambda contig, start, end: varcode.reference_range(genome, contig, start, end))
        if use_genomic_sequence
        else None
    )
    effect = _call(
        varcode.predict_transcript_model_effect,
        items,
        transcript,
        germline_variants=germline,
        sequence_provider=provider,
        max_hypotheses=cap,
    )
    candidates = list(getattr(effect, "candidates", ()) or ())
    return _document(
        "predict_transcript_model",
        variants=[_variant_record(item) for item in items],
        transcript_id=transcript_id,
        effect=_effect_record(effect),
        candidates=[_candidate_record(candidate) for candidate in candidates],
        candidate_count=len(candidates),
        genomic_sequence_used=use_genomic_sequence,
        cache_directory=info["cache_directory"],
    )


__all__ = ["varcode_annotate_cis_haplotype", "varcode_predict_transcript_model"]
_register_public_verbose(globals())
