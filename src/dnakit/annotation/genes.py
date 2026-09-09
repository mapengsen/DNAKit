"""Annotate raw genomic DNA with gene names from a matching reference assembly."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

from dnakit.backends import execute_bounded_command
from dnakit.core._json import FrozenDict
from dnakit.download.files import _write_manifest
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.references import supported_genome_aliases

from ._gene_mapping import annotate_mappings, read_mappings
from ._gene_reference import (
    _ACCESSION,
    _json_object,
    fasta_lengths,
    prepare_reference,
    read_genes,
    sha256_file,
)
from .gene_models import GeneAnnotationConfig, GeneAnnotationResult


def _executable(config: GeneAnnotationConfig) -> Path:
    raw = config.minimap2_executable
    found = shutil.which(raw)
    path = Path(found or raw).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise BackendUnavailableError(
            "minimap2 is required for gene annotation.",
            code="MINIMAP2_UNAVAILABLE",
            hint="Install minimap2 or set minimap2_executable to its executable path.",
        )
    return path


def _index(
    fasta: Path,
    fasta_sha256: str,
    reference_bases: int,
    executable: Path,
    version: str,
    config: GeneAnnotationConfig,
) -> Path:
    # Presets affect the index's k/w values. They must be part of its identity.
    identity = {"fasta_sha256": fasta_sha256, "preset": config.preset, "minimap2_version": version}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    parent = Path(config.cache_dir).expanduser().resolve() / "indexes"
    parent.mkdir(parents=True, exist_ok=True)
    directory = parent / key
    if not directory.exists():
        lock = parent / f"{key}.lock"
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError as exc:
            raise ConfigurationError(
                "The reference index is being built; retry after it finishes."
            ) from exc
        os.close(descriptor)
        try:
            with tempfile.TemporaryDirectory(prefix=".gene-index-", dir=parent) as raw:
                work = Path(raw)
                data = work / "index"
                data.mkdir()
                output = data / "reference.mmi"
                # A single-part index is necessary for trustworthy MAPQ. NCBI
                # plant genomes may exceed minimap2's default indexing batch size.
                command = execute_bounded_command(
                    executable,
                    (
                        "-x",
                        config.preset,
                        "-I",
                        str(reference_bases + 1),
                        "-d",
                        str(output),
                        "-t",
                        str(config.threads),
                        str(fasta),
                    ),
                    backend_id="minimap2",
                    cwd=data,
                    timeout_seconds=config.timeout_seconds,
                    max_output_bytes=10_000_000,
                    monitored_output_paths=(output,),
                    max_monitored_output_bytes=config.max_index_bytes,
                )
                if command.return_code or not output.is_file() or output.stat().st_size == 0:
                    raise BackendExecutionError(
                        "minimap2 failed to build its reference index.",
                        code="GENE_INDEX_FAILED",
                        context={
                            "return_code": command.return_code,
                            "output": command.output[-2000:],
                        },
                    )
                _write_manifest(
                    {**identity, "index_sha256": sha256_file(output)},
                    data / "index.json",
                    overwrite=False,
                )
                os.rename(data, directory)
        finally:
            lock.unlink()
    output = directory / "reference.mmi"
    manifest = _json_object(directory / "index.json")
    if (
        any(manifest.get(key) != value for key, value in identity.items())
        or output.stat().st_size > config.max_index_bytes
        or sha256_file(output) != manifest.get("index_sha256")
    ):
        raise ConfigurationError(
            "Reference index cache failed validation; use a new cache directory."
        )
    return output


def annotate_genes(
    sequence: str,
    *,
    species: str | int,
    assembly: str | None = None,
    sequence_id: str = "query",
    reference_fasta: str | os.PathLike[str] | None = None,
    annotation_gff: str | os.PathLike[str] | None = None,
    config: GeneAnnotationConfig | None = None,
    progress: Callable[[str], None] | None = None,
) -> GeneAnnotationResult:
    """Map genomic DNA and project named genes from a matching GFF3 annotation.

    ``assembly=None`` selects the species' current RefSeq reference assembly on
    first use, then pins that exact reference/annotation pair in ``cache_dir``.
    Explicit assemblies can be versioned accessions, supported aliases, or an
    unambiguous species-scoped assembly name. Multiple/no candidates are errors.

    Pass both local files for offline use (plain or gzip FASTA/GFF3). Their
    species is caller-declared; assembly declarations and sequence IDs/bounds
    are checked when available. Genomic DNA only: this is annotation transfer,
    not de novo gene prediction or a transcript/CDS annotation pipeline.

    Results retain overlapping genes and cover every input base with a gene,
    intergenic, or unknown segment. Unmapped/low-confidence bases, insertions
    and ambiguous primary mappings stay unknown. Python Interval objects use
    0-based half-open coordinates; ``result.to_dict()`` uses 1-based closed.
    """

    resolved = config if config is not None else GeneAnnotationConfig()
    if not isinstance(resolved, GeneAnnotationConfig):
        raise ConfigurationError("config must be GeneAnnotationConfig.")
    if not isinstance(sequence, str):
        raise ConfigurationError("sequence must be raw DNA text.")
    dna = "".join(sequence.split()).upper()
    if (
        not 1 <= len(dna) <= resolved.max_query_bases
        or re.fullmatch(r"[ACGTRYSWKMBDHVN]+", dna) is None
    ):
        raise ConfigurationError(
            "sequence must be non-empty ungapped IUPAC DNA within max_query_bases."
        )
    if (
        isinstance(species, bool)
        or not isinstance(species, (str, int))
        or (isinstance(species, int) and species < 1)
        or not str(species).strip()
        or len(str(species)) > 512
    ):
        raise ConfigurationError("species must be a scientific name, alias, or positive TaxID.")
    organism = str(species).strip()
    if assembly is not None and (not isinstance(assembly, str) or not assembly.strip()):
        raise ConfigurationError("assembly must be None or a non-empty assembly name/accession.")
    requested = None if assembly is None else assembly.strip()
    if (
        not isinstance(sequence_id, str)
        or not sequence_id.strip()
        or any(c in sequence_id for c in "\r\n\x00")
    ):
        raise ConfigurationError("sequence_id must be non-empty single-line text.")
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.")
    if (reference_fasta is None) != (annotation_gff is None):
        raise ConfigurationError("Provide both reference_fasta and annotation_gff together.")

    def notify(message: str) -> None:
        if progress is not None:
            progress(message)

    executable = _executable(resolved)
    notify("Preparing reference FASTA and gene annotation")
    reference: dict[str, Any]
    if reference_fasta is not None and annotation_gff is not None:
        fasta = Path(reference_fasta).expanduser().resolve(strict=True)
        gff = Path(annotation_gff).expanduser().resolve(strict=True)
        if not fasta.is_file() or not gff.is_file():
            raise ConfigurationError("Local reference paths must be files.")
        reference = {
            "assembly": requested,
            "source": "local files",
            "selection": "local",
            "requested_assembly": requested,
            "species_validated": False,
            "species_validation": "caller-declared; no remote taxonomy lookup",
            "fasta_sha256": sha256_file(fasta),
            "gff3_sha256": sha256_file(gff),
        }
    else:
        fasta, gff, reference = prepare_reference(organism, requested, resolved)
    lengths = fasta_lengths(fasta, resolved)
    genes, represented, headers = read_genes(gff, lengths, resolved)
    declared = _ACCESSION.search(headers.get("genome-build-accession", ""))
    if reference["selection"] == "local":
        actual = declared[0].upper() if declared is not None else None
        expected = (
            None
            if requested is None
            else supported_genome_aliases().get(requested.casefold(), requested)
        )
        if expected is not None:
            if _ACCESSION.fullmatch(expected) and actual is not None and expected.upper() != actual:
                raise ConfigurationError("Local GFF3 declares a different assembly accession.")
            if not _ACCESSION.fullmatch(expected) and headers.get("genome-build") not in {
                None,
                expected,
            }:
                raise ConfigurationError("Local GFF3 declares a different assembly name.")
        reference["assembly"] = actual or expected or f"local:{reference['fasta_sha256'][:16]}"
        reference["annotation_headers"] = headers
    reference.update({"input_species": organism, "fasta_path": str(fasta), "gff3_path": str(gff)})
    notify("Preparing minimap2 index")
    with tempfile.TemporaryDirectory(prefix="dnakit-gene-mapping-") as raw:
        work = Path(raw)
        version_result = execute_bounded_command(
            executable,
            ("--version",),
            backend_id="minimap2",
            cwd=work,
            timeout_seconds=10,
        )
        version = version_result.output.strip()
        if version_result.return_code or not version:
            raise BackendExecutionError("Could not determine minimap2 version.")
        index = _index(
            fasta,
            str(reference["fasta_sha256"]),
            sum(lengths.values()),
            executable,
            version,
            resolved,
        )
        query = work / "query.fasta"
        query.write_text(">query\n" + dna + "\n", encoding="ascii")
        paf = work / "alignment.paf"
        notify("Mapping DNA to the selected reference assembly")
        arguments = (
            "-c",
            "-x",
            resolved.preset,
            "--secondary=yes",
            "-N",
            "50",
            "-t",
            str(resolved.threads),
            "-o",
            str(paf),
            str(index),
            str(query),
        )
        command = execute_bounded_command(
            executable,
            arguments,
            backend_id="minimap2",
            cwd=work,
            timeout_seconds=resolved.timeout_seconds,
            max_output_bytes=10_000_000,
            monitored_output_paths=(paf,),
            max_monitored_output_bytes=resolved.max_alignment_bytes,
        )
        if command.return_code or not paf.is_file():
            raise BackendExecutionError(
                "minimap2 did not produce a successful genomic alignment.",
                code="GENE_MAPPING_FAILED",
                context={"return_code": command.return_code, "output": command.output[-2000:]},
            )
        mappings = read_mappings(paf, len(dna), lengths, resolved)
        paf_sha256 = sha256_file(paf)
    notify("Projecting gene intervals and resolving overlaps")
    records, segments = annotate_mappings(len(dna), mappings, genes, represented)
    return GeneAnnotationResult(
        sequence_id,
        len(dna),
        FrozenDict(reference),
        records,
        segments,
        tuple(mapping.to_dict() for mapping in mappings),
        FrozenDict(
            {
                "method": "reference_gene_annotation_transfer",
                "algorithm_version": "dnakit-gene-annotation-v1",
                "query_sha256": hashlib.sha256(dna.encode("ascii")).hexdigest(),
                "minimap2_version": version,
                "preset": resolved.preset,
                "min_mapq": resolved.min_mapq,
                "min_identity": resolved.min_identity,
                "min_alignment_length": resolved.min_alignment_length,
                "paf_sha256": paf_sha256,
                "annotated_fraction": sum(len(s.location) for s in segments if s.kind != "unknown")
                / len(dna),
                "intergenic_definition": (
                    "reliably mapped bases outside gene features in the selected GFF3"
                ),
                "gene_name_policy": (
                    "null when absent; display_name falls back to locus_tag/gene_id"
                ),
            }
        ),
    )


__all__ = ["annotate_genes"]
