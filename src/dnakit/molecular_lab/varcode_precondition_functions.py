"""varcode operations that need RNA-seq evidence or a splice outcome set, on a local GTF and genome FASTA.

The RNA functions take the observed junctions or fusion models as explicit inputs: a sequence in transcript
orientation with its model ID, an optional partner transcript and a read count, or ExACTO structure and integration
tables written by an RNA tool. The annotation is the same local pyensembl index that ``varcode_gap_functions`` builds
from a GTF and a genome FASTA. No RNA alignment is read and nothing is downloaded. Read counts are kept as evidence
and never become likelihoods, as in varcode.

``varcode_enumerate_splice_outcomes`` predicts a variant's effect on one transcript and wraps a splice-signal
disruption in varcode's splice outcome set (exon skipping, intron retention, cryptic sites), using the genome
sequence of the local reference when ``use_genomic_sequence`` is True.

Results are JSON-serialisable dictionaries, built like those of the other varcode functions. Missing varcode raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-varcode]"`` hint.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab.varcode_gap_functions import (
    _SV_TYPES,
    _VARCODE_EXTRA,
    _bases,
    _call,
    _candidate_record,
    _document,
    _effect_record,
    _genome,
    _identifier,
    _int,
    _mapping,
    _mutant_record,
    _path,
    _text,
    _transcript,
    _variant,
    _variant_record,
    _vc,
)

_ENDS = ("five_prime", "three_prime")
_MAX_FUSIONS = 200
_MAX_EFFECTS = 1_000
_MAX_CANDIDATES = 1_000
_MAX_SEQUENCE = 100_000
_MAX_CALLS = 500


# --- validation ------------------------------------------------------------------------------


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _structural(genome: Any, record: object, name: str) -> tuple[Any, dict[str, Any]]:
    """Build a varcode StructuralVariant on the local reference from contig, start, end and sv_type."""
    item = _mapping(record, name)
    missing = sorted({"contig", "start", "end", "sv_type"} - set(item))
    if missing:
        raise ConfigurationError(f"{name} is missing fields: {', '.join(missing)}.")
    contig = _identifier(item["contig"], f"{name}.contig")
    start = _int(item["start"], f"{name}.start", minimum=1)
    end = _int(item["end"], f"{name}.end", minimum=start)
    sv_type = _text(item["sv_type"], f"{name}.sv_type").upper()
    if sv_type not in _SV_TYPES:
        raise ConfigurationError(f"{name}.sv_type must be one of: {', '.join(_SV_TYPES)}.")
    structural = _call(
        _vc("structural_variant").StructuralVariant, contig, start, sv_type, end=end, genome=genome
    )
    return structural, {"contig": contig, "start": start, "end": end, "sv_type": sv_type}


def _fusion_candidate(genome: Any, structural: Any, record: object, name: str) -> Any:
    """Build one observed fusion outcome (varcode.rna_evidence.make_fusion_outcome) from a mapping."""
    item = _mapping(record, name)
    missing = sorted({"transcript_id", "transcript_model_id", "sequence"} - set(item))
    if missing:
        raise ConfigurationError(f"{name} is missing fields: {', '.join(missing)}.")
    sequence = _bases(item["sequence"], f"{name}.sequence", allow_empty=False)
    if len(sequence) > _MAX_SEQUENCE:
        raise ConfigurationError(f"{name}.sequence may hold at most {_MAX_SEQUENCE} bases.")
    anchor = _transcript(genome, item["transcript_id"])
    partner_id = item.get("partner_transcript_id")
    partner = None if partner_id is None else _transcript(genome, partner_id)
    cds_start = item.get("cds_start")
    read_count = item.get("read_count")
    return _call(
        _vc("rna_evidence").make_fusion_outcome,
        structural,
        anchor,
        sequence=sequence,
        transcript_model_id=_identifier(item["transcript_model_id"], f"{name}.transcript_model_id"),
        partner_transcript=partner,
        cds_start=None if cds_start is None else _int(cds_start, f"{name}.cds_start", minimum=0),
        source=_identifier(item.get("source", "rna"), f"{name}.source"),
        read_count=None
        if read_count is None
        else _int(read_count, f"{name}.read_count", minimum=0),
    )


def _fusion_list(values: object, name: str) -> list[Mapping[str, Any]]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of fusion mappings.")
    if not values or len(values) > _MAX_FUSIONS:
        raise ConfigurationError(f"{name} must hold between 1 and {_MAX_FUSIONS} entries.")
    return [_mapping(item, f"{name}[{index}]") for index, item in enumerate(values)]


def _exacto_variants(genome: Any, values: object) -> tuple[dict[str, Any], dict[str, Any]]:
    """Map Exacto DNA call IDs to StructuralVariant objects; returns the objects and their records."""
    if not isinstance(values, Mapping) or not values or len(values) > _MAX_FUSIONS:
        raise ConfigurationError(f"structural_variants must map 1 to {_MAX_FUSIONS} DNA call IDs.")
    variants: dict[str, Any] = {}
    records: dict[str, Any] = {}
    for call_id, record in values.items():
        key = _text(call_id, "structural_variants key")
        variants[key], records[key] = _structural(genome, record, f"structural_variants[{key}]")
    return variants, records


def _transcript_ids(genome: Any, values: object, name: str) -> list[Any] | None:
    if values is None:
        return None
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of transcript IDs or None.")
    if len(values) > _MAX_CALLS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_CALLS} transcript IDs.")
    return [_transcript(genome, item) for item in values]


# --- public operations -----------------------------------------------------------------------


def varcode_make_fusion_outcome(
    gtf: str,
    genome_fasta: str,
    structural_variant: Mapping[str, Any],
    transcript_id: str,
    sequence: str,
    transcript_model_id: str,
    *,
    partner_transcript_id: str | None = None,
    cds_start: int | None = None,
    source: str = "rna",
    read_count: int | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import one observed RNA junction or fusion model on a structural variant (varcode.make_fusion_outcome).

    ``sequence`` is the observed 5'-to-3' sequence of the model (not genomic-forward). ``transcript_id`` is its
    annotated 5' anchor. Give ``partner_transcript_id`` only when the observed path joins that transcript in sense
    orientation; otherwise the outcome is a translocation to an intergenic region. ``cds_start`` (zero-based) asks
    for translation of a complete start-to-stop ORF in the sequence; the protein is then a prediction from RNA.
    ``structural_variant`` has ``contig``, ``start``, ``end`` and ``sv_type`` (DEL, DUP, INV, INS or CNV).
    """
    require("varcode", _VARCODE_EXTRA)
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    structural, sv_record = _structural(genome, structural_variant, "structural_variant")
    candidate = _fusion_candidate(
        genome,
        structural,
        {
            "transcript_id": transcript_id,
            "transcript_model_id": transcript_model_id,
            "sequence": sequence,
            "partner_transcript_id": partner_transcript_id,
            "cds_start": cds_start,
            "source": source,
            "read_count": read_count,
        },
        "fusion",
    )
    return _document(
        "make_fusion_outcome",
        structural_variant=sv_record,
        transcript_id=str(transcript_id),
        partner_transcript_id=partner_transcript_id,
        transcript_model_id=str(transcript_model_id),
        candidate=_candidate_record(candidate),
        mutant_transcript=_mutant_record(candidate.effect.mutant_transcript),
        cache_directory=info["cache_directory"],
    )


def varcode_apply_rna_evidence_to_effects(
    gtf: str,
    genome_fasta: str,
    structural_variant: Mapping[str, Any],
    fusion_outcomes: list[Mapping[str, Any]],
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Attach observed fusion outcomes to the predicted effects of a structural variant (varcode.apply_rna_evidence_to_effects).

    Each entry of ``fusion_outcomes`` has ``transcript_id``, ``transcript_model_id`` and ``sequence``, and may have
    ``partner_transcript_id``, ``cds_start``, ``read_count`` and ``source`` (see ``varcode_make_fusion_outcome``).
    The entries are built into a varcode ``RNAEvidence`` resolver and applied to the effects predicted on the local
    reference. Each effect in the result lists the observed candidates it gained; the DNA-predicted effect is kept.
    """
    varcode = require("varcode", _VARCODE_EXTRA)
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    structural, sv_record = _structural(genome, structural_variant, "structural_variant")
    items = _fusion_list(fusion_outcomes, "fusion_outcomes")
    candidates = [
        _fusion_candidate(genome, structural, item, f"fusion_outcomes[{index}]")
        for index, item in enumerate(items)
    ]
    resolver = _call(_vc("rna_evidence").RNAEvidence, candidates)
    effects = _call(_vc("effects").predict_variant_effects, structural)
    resolved = _call(varcode.apply_rna_evidence_to_effects, effects, resolver)
    rows: list[dict[str, Any]] = []
    for effect in list(resolved)[:_MAX_EFFECTS]:
        row = _effect_record(effect)
        row["candidates"] = [
            _candidate_record(candidate)
            for candidate in list(getattr(effect, "candidates", ()) or ())[:_MAX_CANDIDATES]
        ]
        rows.append(row)
    return _document(
        "apply_rna_evidence_to_effects",
        structural_variant=sv_record,
        fusion_count=len(candidates),
        effects=rows,
        count=len(rows),
        observed_candidate_count=sum(len(row["candidates"]) for row in rows),
        cache_directory=info["cache_directory"],
    )


def varcode_load_exacto_fusions(
    gtf: str,
    genome_fasta: str,
    structures_path: str,
    integrated_path: str,
    structural_variants: Mapping[str, Mapping[str, Any]],
    *,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import ExACTO structure and integration tables linked to structural variants (varcode.load_exacto_fusions).

    ``structures_path`` and ``integrated_path`` are ExACTO's transcript-structure and integrated-variant TSV files
    (gzip is accepted). ``structural_variants`` maps each ExACTO DNA call ID used in the integration table to a
    structural variant (``contig``, ``start``, ``end``, ``sv_type``); only these calls are imported. Sequences are
    already in transcript orientation. Linear two-locus models become fusion outcomes; circular, multi-gene and
    incomplete paths raise an error instead of being guessed. Native primary-structure peptides and explicit ORF
    starts are not exposed by this wrapper.
    """
    require("varcode", _VARCODE_EXTRA)
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    structures = _path(structures_path, "structures_path")
    integrated = _path(integrated_path, "integrated_path")
    variants, records = _exacto_variants(genome, structural_variants)
    resolver = _call(
        _vc("exacto").load_exacto_fusions,
        str(structures),
        str(integrated),
        variants_by_id=variants,
    )
    rows = [
        _candidate_record(candidate) for candidate in list(resolver.candidates)[:_MAX_CANDIDATES]
    ]
    return _document(
        "load_exacto_fusions",
        structural_variants=records,
        candidates=rows,
        count=len(rows),
        structures_path=str(structures),
        integrated_path=str(integrated),
        cache_directory=info["cache_directory"],
    )


def varcode_reference_completion_hypotheses(
    gtf: str,
    genome_fasta: str,
    structures_path: str,
    integrated_path: str,
    structural_variants: Mapping[str, Mapping[str, Any]],
    *,
    ends: list[str] | None = None,
    five_prime_transcript_ids: list[str] | None = None,
    three_prime_transcript_ids: list[str] | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Complete the mapped ends of ExACTO fusions as explicit reference hypotheses (varcode.reference_completion_hypotheses).

    The observed fusions are loaded as in ``varcode_load_exacto_fusions``. For each one, the unmapped 5' and/or 3'
    end is extended with reference bases from an annotated transcript, one hypothesis per choice. ``ends`` selects
    ``five_prime``, ``three_prime`` or both (the default). The transcript lists give the alternatives to consider
    at each end; None uses the annotated transcript of that end. The hypotheses are never reported as RNA
    observations: they are returned with their assumed spans, and the read counts are not used as support.
    """
    require("varcode", _VARCODE_EXTRA)
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    structures = _path(structures_path, "structures_path")
    integrated = _path(integrated_path, "integrated_path")
    chosen = list(_ENDS if ends is None else ends)
    if not chosen or set(chosen) - set(_ENDS):
        raise ConfigurationError("ends must contain five_prime and/or three_prime.")
    five = _transcript_ids(genome, five_prime_transcript_ids, "five_prime_transcript_ids")
    three = _transcript_ids(genome, three_prime_transcript_ids, "three_prime_transcript_ids")
    variants, records = _exacto_variants(genome, structural_variants)
    resolver = _call(
        _vc("exacto").load_exacto_fusions,
        str(structures),
        str(integrated),
        variants_by_id=variants,
    )
    observed = list(resolver.candidates)[:_MAX_CANDIDATES]
    hypotheses: list[dict[str, Any]] = []
    for candidate in observed:
        for hypothesis in _call(
            _vc("reference_completion").reference_completion_hypotheses,
            candidate,
            five_prime_transcripts=five,
            three_prime_transcripts=three,
            ends=tuple(chosen),
        ):
            hypotheses.append(_candidate_record(hypothesis))
    return _document(
        "reference_completion_hypotheses",
        structural_variants=records,
        ends=chosen,
        observed_count=len(observed),
        hypothesis_count=len(hypotheses),
        hypotheses=hypotheses[:_MAX_CANDIDATES],
        cache_directory=info["cache_directory"],
    )


def varcode_enumerate_splice_outcomes(
    gtf: str,
    genome_fasta: str,
    variant: Mapping[str, Any],
    transcript_id: str,
    *,
    use_genomic_sequence: bool = True,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate the splice mechanisms behind a splice-signal disruption on one transcript (varcode.enumerate_splice_outcomes).

    The variant's effect on the transcript is predicted first. A splice donor, acceptor or splice-site effect is then
    wrapped in varcode's outcome set, with the candidate mechanisms (for example exon skipping, intron retention and
    cryptic sites). With ``use_genomic_sequence`` True the sequence of the local reference is supplied, so cryptic
    sites are scanned and intron retention is materialised. Other effects are returned unchanged and have no
    candidates. Use a transcript with several exons; a single-exon transcript has no splice junction.
    """
    varcode = require("varcode", _VARCODE_EXTRA)
    genome, info = _genome(gtf, genome_fasta, cache_dir)
    item = _variant(genome, variant)
    transcript = _transcript(genome, transcript_id)
    effect = _call(_vc("effects").predict_variant_effect_on_transcript, item, transcript)
    provider: Callable[[str, int, int], str] | None = None
    if _flag(use_genomic_sequence, "use_genomic_sequence"):

        def provider(contig: str, start: int, end: int) -> str:
            return varcode.reference_range(genome, contig, start, end)

    outcome = _call(
        _vc("splice_outcomes").enumerate_splice_outcomes, effect, genomic_sequence=provider
    )
    outcome_set = _vc("splice_outcomes").SpliceOutcomeSet
    candidates = list(getattr(outcome, "candidates", ()) or ())[:_MAX_CANDIDATES]
    return _document(
        "enumerate_splice_outcomes",
        variant=_variant_record(item),
        transcript_id=str(transcript_id),
        input_effect=_effect_record(effect),
        outcome_class=type(outcome).__name__,
        is_outcome_set=isinstance(outcome, outcome_set),
        genomic_sequence_used=provider is not None,
        candidates=[_candidate_record(candidate) for candidate in candidates],
        count=len(candidates),
        cache_directory=info["cache_directory"],
    )


_register_public_verbose(globals())
