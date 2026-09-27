"""Optional Prodigal gene prediction through the official Pyrodigal API.

This module adapts Pyrodigal; it does not reimplement Prodigal. Cite Larralde
(2022), https://doi.org/10.21105/joss.04296, and Hyatt et al. (2010),
https://doi.org/10.1186/1471-2105-11-119. Pyrodigal is distributed separately
under GPL-3.0. Its API and coordinate contract are documented at
https://pyrodigal.readthedocs.io/en/stable/api/genes.html.
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from typing import Any

from dnakit._progress import Progress
from dnakit.bioinformatics._common import require, sequence_text
from dnakit.core import DNA, DNARecord, DNASequence, Topology
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_SequenceInput = str | DNASequence | DNARecord | DNA
_CITATIONS = (
    "https://doi.org/10.21105/joss.04296",
    "https://doi.org/10.1186/1471-2105-11-119",
)


@dataclass(frozen=True)
class PredictedGene:
    """A predicted CDS, with genomic coordinates and strand-oriented sequence.

    ``start`` and ``end`` are 0-based, half-open reference coordinates on
    either strand. ``original_begin``/``original_end`` retain Pyrodigal's
    1-based inclusive coordinates. ``partial_begin``/``partial_end`` retain
    the upstream flags without reinterpreting them as 5-prime/3-prime labels.
    Scores and confidence percentages are model outputs, not experimental
    confirmation of a gene or annotation of its molecular function.
    """

    gene_id: str
    start: int
    end: int
    strand: int
    original_begin: int
    original_end: int
    cds: str
    protein: str
    translation_table: int
    partial_begin: bool
    partial_end: bool
    start_type: str
    score: float
    coding_score: float
    start_score: float
    confidence_percent: float


@dataclass(frozen=True)
class GenePredictionResult:
    """Serializable predictions with method, sequence identity and parameters."""

    sequence_id: str
    sequence_length: int
    input_sha256: str
    ambiguous_bases_to_n: int
    genes: tuple[PredictedGene, ...]
    backend_version: str
    prodigal_version: str
    meta: bool
    closed: bool
    mask: bool
    include_stop: bool
    training_sha256: str | None
    training_translation_table: int | None
    method: str = "Prodigal via Pyrodigal"
    coordinate_system: str = "0-based half-open"
    original_coordinate_system: str = "1-based inclusive"
    citations: tuple[str, ...] = _CITATIONS

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible mapping, without any backend objects."""
        result = asdict(self)
        result["genes"] = [asdict(gene) for gene in self.genes]
        result["citations"] = list(self.citations)
        return result


def _canonical_with_unknowns(sequence: _SequenceInput) -> tuple[str, str, int]:
    value = sequence.record if isinstance(sequence, DNA) else sequence
    value = value.sequence if isinstance(value, DNARecord) else value
    if isinstance(value, DNASequence) and value.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Pyrodigal predicts on linear contigs; explicitly linearize circular DNA "
            "and account for genes crossing the chosen origin."
        )
    text = sequence_text(sequence)
    unknowns = sum(base not in "ACGT" for base in text)
    # Preserve length and all reference positions; Prodigal represents IUPAC
    # ambiguity as N and does not infer the missing base identity.
    normalized = "".join(base if base in "ACGT" else "N" for base in text)
    return normalized, hashlib.sha256(text.encode("ascii")).hexdigest(), unknowns


def predict_prokaryotic_genes(
    sequence: _SequenceInput,
    *,
    sequence_id: str = "sequence",
    meta: bool = True,
    training_sequence: _SequenceInput | None = None,
    translation_table: int | None = None,
    closed: bool = False,
    mask: bool = True,
    include_stop: bool = True,
    verbose: bool = True,
) -> GenePredictionResult:
    """Predict prokaryotic CDSs using an optional, locally installed Pyrodigal.

    Install with ``pip install 'dnakit[genomics-genes]'``. Metagenomic mode
    uses Prodigal's supplied models and selects the translation table from
    the winning model; each gene reports that table. To train a genome model,
    explicitly pass ``meta=False`` and ``training_sequence=...``. Training
    uses table 11 unless ``translation_table`` is provided. Pyrodigal checks
    training length and rejects unsuitable short inputs; no silent fallback
    to metagenomic mode is performed.

    Input must be non-empty, ungapped IUPAC DNA. Ambiguous IUPAC letters are
    represented as N without changing positions. Returned CDS strings are
    oriented 5-prime to 3-prime on the predicted gene strand. Protein strings
    use the gene's prediction table and include a terminal ``*`` when present
    unless ``include_stop=False``. ``mask=True`` uses Pyrodigal's default
    50-base unknown-region masking threshold. ``closed`` controls contig
    edge predictions, not circular topology.

    Progress shows prediction/training stages; the backend does not expose
    nucleotide-level progress. No software, weights or sequence data are
    downloaded by this call. This predicts CDS structure, not gene function.
    """
    for name, value in (
        ("meta", meta),
        ("closed", closed),
        ("mask", mask),
        ("include_stop", include_stop),
        ("verbose", verbose),
    ):
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be a boolean.")
    if not isinstance(sequence_id, str) or not sequence_id.strip():
        raise ConfigurationError("sequence_id must be a non-empty string.")
    if meta and (training_sequence is not None or translation_table is not None):
        raise ConfigurationError(
            "Metagenomic mode selects its own model/table; training_sequence and "
            "translation_table require meta=False."
        )
    if not meta and training_sequence is None:
        raise ConfigurationError("meta=False requires an explicit training_sequence.")
    if translation_table is not None and (
        isinstance(translation_table, bool)
        or not isinstance(translation_table, int)
        or translation_table < 1
    ):
        raise ConfigurationError("translation_table must be a positive integer.")

    text, input_digest, unknowns = _canonical_with_unknowns(sequence)
    training_text: str | None = None
    training_digest: str | None = None
    if training_sequence is not None:
        training_text, training_digest, _ = _canonical_with_unknowns(training_sequence)
    training_table = None if meta else (11 if translation_table is None else translation_table)
    pyrodigal = require("pyrodigal", "genomics-genes")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Predicting prokaryotic genes", total=2 if meta else 3)
        try:
            finder = pyrodigal.GeneFinder(meta=meta, closed=closed, mask=mask)
            if training_text is not None:
                finder.train(training_text, translation_table=training_table)
                progress.advance(task)
            predictions = finder.find_genes(text)
            progress.advance(task)
            genes = tuple(
                PredictedGene(
                    gene_id=f"{sequence_id}_{index}",
                    start=int(gene.begin) - 1,
                    end=int(gene.end),
                    strand=int(gene.strand),
                    original_begin=int(gene.begin),
                    original_end=int(gene.end),
                    cds=gene.sequence(),
                    protein=gene.translate(include_stop=include_stop),
                    translation_table=int(gene.translation_table),
                    partial_begin=bool(gene.partial_begin),
                    partial_end=bool(gene.partial_end),
                    start_type=str(gene.start_type),
                    score=float(gene.score),
                    coding_score=float(gene.cscore),
                    start_score=float(gene.sscore),
                    confidence_percent=float(gene.confidence()),
                )
                for index, gene in enumerate(predictions, 1)
            )
            progress.advance(task)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                f"Pyrodigal rejected the gene prediction input: {exc}"
            ) from exc
        except RuntimeError as exc:
            raise BackendExecutionError(f"Pyrodigal gene prediction failed: {exc}") from exc

    return GenePredictionResult(
        sequence_id=sequence_id,
        sequence_length=len(text),
        input_sha256=input_digest,
        ambiguous_bases_to_n=unknowns,
        genes=genes,
        backend_version=str(pyrodigal.__version__),
        prodigal_version=str(pyrodigal.PRODIGAL_VERSION),
        meta=meta,
        closed=closed,
        mask=mask,
        include_stop=include_stop,
        training_sha256=training_digest,
        training_translation_table=training_table,
    )


__all__ = ["GenePredictionResult", "PredictedGene", "predict_prokaryotic_genes"]
