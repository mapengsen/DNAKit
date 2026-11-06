"""Explicit reference-ORF and observed-RNA annotations via varcode >=10.16.3.

Reference ORF models and taxonomy/assembly context are supplied by the caller.
No catalogue or reference data are downloaded. Annotation preserves upstream
unresolved outcomes and never treats reference peptide evidence as mutant RNA
or translation evidence.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab.varcode_gap_functions import (
    _call,
    _document,
    _genome,
    _transcript,
    _variant,
)

_EXTRA = "bio-varcode-orf"


def _api() -> tuple[Any, Any]:
    varcode = require("varcode", _EXTRA)
    if not all(
        hasattr(varcode, name)
        for name in ("annotate_orf_variants", "annotate_orf_rna_product", "predict_orf_effects")
    ):
        raise BackendUnavailableError(
            "This varcode version has no reference-ORF annotation API.",
            hint='Install "varcode>=10.16.3" and "orfdata>=0.34.1" (dnakit[bio-varcode-orf]).',
        )
    return varcode, require("orfdata", _EXTRA)


def _reference_genome(
    gtf: str, genome_fasta: str, context_spec: Mapping[str, Any], cache_dir: str | None
) -> tuple[Any, Any, Any, Any]:
    varcode, orfdata = _api()
    if not isinstance(context_spec, Mapping):
        raise ConfigurationError(
            "reference_context must supply taxon_id, species, genome_build and annotation_name."
        )
    context = _call(orfdata.ReferenceContext, **dict(context_spec))
    if context.genome_build is None or context.annotation_name is None:
        raise ConfigurationError(
            "Binding local genomic transcripts requires genome_build and annotation_name."
        )
    _, info = _genome(gtf, genome_fasta, cache_dir)
    root = Path(info["cache_directory"])
    pyensembl = require("pyensembl", _EXTRA)
    # Build an independent view with caller-specified real reference identity;
    # never mutate the shared synthetic Genome cached by other DNAKit wrappers.
    genome = pyensembl.Genome(
        reference_name=context.genome_build,
        annotation_name=context.annotation_name,
        annotation_version=context.annotation_version,
        gtf_path_or_url=str(root / "annotation.gtf"),
        transcript_fasta_paths_or_urls=[str(root / "cdna.fa")],
        protein_fasta_paths_or_urls=[str(root / "protein.fa")],
        genome_fasta_path_or_url=str(Path(genome_fasta).expanduser().resolve()),
        cache_directory_path=str(root / "orf_pyensembl"),
    )
    _call(genome.index)
    return varcode, orfdata, genome, context


def _model(orfdata: Any, genome: Any, context: Any, spec: Mapping[str, Any]) -> Any:
    if not isinstance(spec, Mapping):
        raise ConfigurationError("Each model must be an ORF record mapping.")
    allowed = {
        "orf_id",
        "gene_id",
        "gene_name",
        "transcript_id",
        "contig",
        "strand",
        "coding_blocks",
        "protein_sequence",
        "biotype",
        "genome_build",
        "taxon_id",
        "species",
        "translation_table",
        "source",
        "start_codon",
    }
    if set(spec) - allowed:
        raise ConfigurationError(
            "Unsupported ORF model fields: " + ", ".join(sorted(set(spec) - allowed))
        )
    fields = dict(spec)
    try:
        fields["coding_blocks"] = tuple(tuple(block) for block in fields["coding_blocks"])
    except (KeyError, TypeError) as exc:
        raise ConfigurationError(
            "coding_blocks must be a list of 1-based inclusive [start, end] pairs."
        ) from exc
    fields.setdefault("genome_build", context.genome_build)
    fields.setdefault("taxon_id", context.taxon_id)
    fields.setdefault("species", context.species)
    fields.setdefault("biotype", "ORF")
    record = _call(orfdata.ORFRecord, **fields, ribosome_studies=(), peptides=())
    return record


def _variants(genome: Any, variants: Sequence[Mapping[str, Any]]) -> list[Any]:
    if (
        isinstance(variants, (str, bytes))
        or not isinstance(variants, Sequence)
        or not variants
        or len(variants) > 1000
    ):
        raise ConfigurationError("variants must contain 1 to 1000 known-cis variant records.")
    return [_variant(genome, record) for record in variants]


def varcode_annotate_reference_orf(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    model: Mapping[str, Any],
    *,
    reference_context: Mapping[str, Any],
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate known-cis variants against one explicit reference ORF with annotate_orf_variants.

    model supplies orf_id, gene_id, gene_name, transcript_id, contig, strand,
    coding_blocks (1-based inclusive genomic pairs), and protein_sequence
    (without terminal stop). reference_context supplies taxon_id, species,
    genome_build and annotation_name, optionally annotation_version. Species,
    sequence and coding-path mismatches are rejected by orfdata binding.
    """
    varcode, orfdata, genome, context = _reference_genome(
        gtf, genome_fasta, reference_context, cache_dir
    )
    record = _model(orfdata, genome, context, model)
    reference = _call(record.bind, _transcript(genome, record.transcript_id), context=context)
    annotation = _call(varcode.annotate_orf_variants, _variants(genome, variants), reference)
    return _document("annotate_reference_orf", annotation=annotation.to_dict())


def varcode_annotate_orf_rna_product(
    gtf: str,
    genome_fasta: str,
    model: Mapping[str, Any],
    cdna_sequence: str,
    base_origins: Sequence[Sequence[Any] | None],
    *,
    reference_context: Mapping[str, Any],
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate an observed RNA product against an ORF via annotate_orf_rna_product.

    base_origins must contain one [contig, 1-based position, strand] or null per
    cDNA base, supplied by RNA reconstruction. The mapping is never inferred
    from peptide similarity. Missing/duplicated reference start codons remain
    unresolved. No RNA reads are assembled by this function.
    """
    varcode, orfdata, genome, context = _reference_genome(
        gtf, genome_fasta, reference_context, cache_dir
    )
    record = _model(orfdata, genome, context, model)
    transcript = _transcript(genome, record.transcript_id)
    reference = _call(record.bind, transcript, context=context)
    if (
        not isinstance(cdna_sequence, str)
        or not cdna_sequence
        or set(cdna_sequence) - set("ACGTNacgtn")
    ):
        raise ConfigurationError("cdna_sequence must be a non-empty DNA sequence.")
    if len(base_origins) != len(cdna_sequence):
        raise ConfigurationError("base_origins needs one entry per cDNA nucleotide.")
    origins = []
    for origin in base_origins:
        if origin is None:
            origins.append(None)
        elif (
            isinstance(origin, (list, tuple))
            and len(origin) == 3
            and isinstance(origin[0], str)
            and isinstance(origin[1], int)
            and not isinstance(origin[1], bool)
            and origin[1] >= 1
            and origin[2] in ("+", "-")
        ):
            origins.append(tuple(origin))
        else:
            raise ConfigurationError(
                "Each base origin must be null or [contig, positive_position, +/-]."
            )
    products = require("varcode.effect_hypotheses", _EXTRA)
    product = products.RealizedTranscriptProduct(
        transcript, cdna_sequence, None, evidence={"base_origins": origins}
    )
    annotation = _call(varcode.annotate_orf_rna_product, reference, product)
    return _document("annotate_orf_rna_product", annotation=annotation.to_dict())


def varcode_predict_orf_catalogue(
    gtf: str,
    genome_fasta: str,
    variants: Sequence[Mapping[str, Any]],
    models: Sequence[Mapping[str, Any]],
    *,
    reference_context: Mapping[str, Any],
    catalogue_metadata: Mapping[str, Any] | None = None,
    cache_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict effects over explicit ORF catalogue models via predict_orf_effects.

    models uses the same record format as varcode_annotate_reference_orf.
    An empty result means no overlapping reference model, not absence of a
    biological effect. Invalid model bindings remain unresolved annotations.
    """
    varcode, orfdata, genome, context = _reference_genome(
        gtf, genome_fasta, reference_context, cache_dir
    )
    if isinstance(models, (str, bytes)) or not isinstance(models, Sequence) or len(models) > 10000:
        raise ConfigurationError("models must be a list of at most 10000 ORF models.")
    records = [_model(orfdata, genome, context, model) for model in models]
    metadata = dict(catalogue_metadata or {})
    metadata.setdefault("genome_build", context.genome_build)
    catalog = _call(orfdata.ORFCatalog, records, metadata)
    annotations = _call(
        varcode.predict_orf_effects,
        _variants(genome, variants),
        catalog=catalog,
        reference_context=context,
    )
    return _document(
        "predict_orf_catalogue",
        annotations=[item.to_dict() for item in annotations],
        count=len(annotations),
    )


__all__ = [
    "varcode_annotate_orf_rna_product",
    "varcode_annotate_reference_orf",
    "varcode_predict_orf_catalogue",
]
_register_public_verbose(globals())
