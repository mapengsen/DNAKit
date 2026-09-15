"""NCBI Datasets data-package downloads."""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.references import download_genome as genome
from dnakit.search import SearchConfig
from dnakit.search._http import build_url, request_text, require_text
from dnakit.search._shared import query_values, quoted_path_values
from dnakit.search.ncbi import (
    DEFAULT_NCBI_DATASETS_URL,
    DEFAULT_NCBI_EUTILS_URL,
    _eutils_parameters,
)

from ._sequence import format_converter, save_sequence, sequence_accession, sequence_target
from .files import dataset, resolved_config
from .models import DatasetDownloadResult, DownloadConfig, DownloadProgress, RemoteFile

_GENOME_INCLUDE = frozenset(
    {
        "CDS_FASTA",
        "GENOME_FASTA",
        "GENOME_GBFF",
        "GENOME_GFF",
        "GENOME_GTF",
        "NONE",
        "PROT_FASTA",
        "RNA_FASTA",
        "SEQUENCE_REPORT",
    }
)
_GENE_INCLUDE = frozenset(
    {"FASTA_3P_UTR", "FASTA_5P_UTR", "FASTA_CDS", "FASTA_GENE", "FASTA_PROTEIN", "FASTA_RNA"}
)
_VIRUS_INCLUDE = frozenset({"GENOME", "CDS", "PROTEIN", "NONE"})
_VIRUS_AUX = frozenset({"ANNOTATION", "BIOSAMPLE_REPORT"})
_VIRUS_ACCESSION = re.compile(r"^[A-Z]{1,6}_?\d{3,}(?:\.\d+)?$", flags=re.IGNORECASE)


def _validate_include(values: Sequence[str], allowed: frozenset[str]) -> tuple[str, ...]:
    resolved = tuple(value.upper() for value in values)
    if not resolved or any(value not in allowed for value in resolved):
        raise ConfigurationError(
            "include contains an unsupported NCBI package file type.",
            code="INVALID_DOWNLOAD_CONFIG",
            context={"allowed": sorted(allowed)},
        )
    return resolved


def _slug(values: Sequence[str]) -> str:
    joined = "_".join(values[:3])
    if len(values) > 3:
        joined += f"_plus_{len(values) - 3}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", joined)


def genome_package(
    accessions: str | Sequence[str],
    output_dir: str | Path,
    *,
    include: Sequence[str] = ("GENOME_FASTA",),
    chromosomes: Sequence[str] = (),
    hydrated: Literal["FULLY_HYDRATED", "DATA_REPORT_ONLY"] = "FULLY_HYDRATED",
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download an NCBI genome ZIP with selected sequences, annotations, and reports."""

    resolved = resolved_config(config)
    values = query_values(accessions, name="accessions", maximum=100)
    if any(
        re.fullmatch(r"GC[AF]_\d+\.\d+", value, flags=re.IGNORECASE) is None for value in values
    ):
        raise ConfigurationError(
            "genome_package requires versioned GCA_/GCF_ accessions.",
            code="INVALID_ASSEMBLY_ACCESSION",
        )
    includes = _validate_include(include, _GENOME_INCLUDE)
    if hydrated not in {"FULLY_HYDRATED", "DATA_REPORT_ONLY"}:
        raise ConfigurationError("Unsupported NCBI hydration mode.")
    chromosome_values = tuple(str(value).strip() for value in chromosomes)
    if any(not value or len(value) > 128 for value in chromosome_values):
        raise ConfigurationError("chromosomes contains an invalid value.")
    url = build_url(
        api_base_url,
        f"/genome/accession/{quoted_path_values(values)}/download",
        (
            ("include_annotation_type", includes),
            ("chromosomes", chromosome_values),
            ("hydrated", hydrated),
        ),
    )
    filename = f"ncbi_genome_{_slug(values)}.zip"
    headers = {} if resolved.api_key is None else {"api-key": resolved.api_key}
    return dataset(
        (RemoteFile(url, filename=filename),),
        output_dir,
        kind="ncbi_genome_package",
        source="NCBI Datasets Genome",
        metadata={
            "accessions": values,
            "include": includes,
            "chromosomes": chromosome_values,
            "hydrated": hydrated,
        },
        config=resolved,
        progress=progress,
        headers=headers,
    )


def annotation(
    accession: str,
    output_dir: str | Path,
    *,
    formats: Sequence[str] = ("GENOME_GFF", "GENOME_GTF", "GENOME_GBFF", "SEQUENCE_REPORT"),
    chromosomes: Sequence[str] = (),
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download an NCBI genome-annotation package (GFF3/GTF/GBFF/reports)."""

    return genome_package(
        accession,
        output_dir,
        include=formats,
        chromosomes=chromosomes,
        config=config,
        progress=progress,
        api_base_url=api_base_url,
    )


def gene(
    gene_ids: int | str | Sequence[int | str],
    output_dir: str | Path,
    *,
    include: Sequence[str] = (
        "FASTA_GENE",
        "FASTA_RNA",
        "FASTA_CDS",
        "FASTA_PROTEIN",
        "FASTA_5P_UTR",
        "FASTA_3P_UTR",
    ),
    accession_filter: Sequence[str] = (),
    include_product_report: bool = True,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download gene, transcript, CDS, UTR, and protein FASTA from NCBI."""

    resolved = resolved_config(config)
    values = query_values(gene_ids, name="gene_ids", maximum=100)
    if any(not value.isdigit() for value in values):
        raise ConfigurationError("gene requires numeric NCBI GeneIDs.", code="INVALID_QUERY")
    includes = _validate_include(include, _GENE_INCLUDE)
    accessions = tuple(str(value).strip() for value in accession_filter)
    if any(not value or len(value) > 256 for value in accessions):
        raise ConfigurationError("accession_filter contains an invalid value.")
    url = build_url(
        api_base_url,
        f"/gene/id/{quoted_path_values(values)}/download",
        (
            ("include_annotation_type", includes),
            ("accession_filter", accessions),
            ("aux_report", "PRODUCT_REPORT" if include_product_report else None),
            (
                "tabular_reports",
                ("DATASET_REPORT", "PRODUCT_REPORT")
                if include_product_report
                else ("DATASET_REPORT",),
            ),
        ),
    )
    headers = {} if resolved.api_key is None else {"api-key": resolved.api_key}
    return dataset(
        (RemoteFile(url, filename=f"ncbi_gene_{_slug(values)}.zip"),),
        output_dir,
        kind="ncbi_gene_package",
        source="NCBI Datasets Gene",
        metadata={"gene_ids": values, "include": includes, "accession_filter": accessions},
        config=resolved,
        progress=progress,
        headers=headers,
    )


def taxonomy(
    tax_ids: int | str | Sequence[int | str],
    output_dir: str | Path,
    *,
    include_names: bool = True,
    include_summary: bool = True,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download taxonomy JSONL, lineage/names, and optional summary TSV."""

    resolved = resolved_config(config)
    values = query_values(tax_ids, name="tax_ids", maximum=100)
    if any(not value.isdigit() for value in values):
        raise ConfigurationError("taxonomy requires numeric NCBI TaxIDs.", code="INVALID_QUERY")
    reports: list[str] = []
    if include_names:
        reports.append("NAMES_REPORT")
    if include_summary:
        reports.append("TAXONOMY_SUMMARY")
    url = build_url(
        api_base_url,
        f"/taxonomy/taxon/{quoted_path_values(values)}/download",
        (("aux_reports", tuple(reports)),),
    )
    headers = {} if resolved.api_key is None else {"api-key": resolved.api_key}
    return dataset(
        (RemoteFile(url, filename=f"ncbi_taxonomy_{_slug(values)}.zip"),),
        output_dir,
        kind="ncbi_taxonomy_package",
        source="NCBI Datasets Taxonomy",
        metadata={"tax_ids": values, "aux_reports": tuple(reports)},
        config=resolved,
        progress=progress,
        headers=headers,
    )


def virus_package(
    query: str | Sequence[str],
    output_dir: str | Path,
    *,
    by: Literal["auto", "accession", "taxon"] = "auto",
    include: Sequence[str] = ("GENOME", "CDS", "PROTEIN"),
    aux_reports: Sequence[str] = ("ANNOTATION", "BIOSAMPLE_REPORT"),
    refseq_only: bool = False,
    annotated_only: bool = False,
    complete_only: bool = False,
    host: str | None = None,
    geo_location: str | None = None,
    released_since: str | None = None,
    updated_since: str | None = None,
    config: DownloadConfig | None = None,
    progress: Callable[[DownloadProgress], None] | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download an NCBI Virus genome, CDS, protein, annotation, and metadata ZIP."""

    resolved = resolved_config(config)
    values = query_values(query, name="query", maximum=100)
    if by == "auto":
        flags = tuple(_VIRUS_ACCESSION.fullmatch(value) is not None for value in values)
        if any(flags) and not all(flags):
            raise ConfigurationError("Do not mix virus accessions and taxon names.")
        by = "accession" if all(flags) else "taxon"
    if by not in {"accession", "taxon"}:
        raise ConfigurationError("virus_package by must be auto, accession, or taxon.")
    if by == "taxon" and len(values) != 1:
        raise ConfigurationError("Virus taxon downloads accept one taxon at a time.")
    includes = _validate_include(include, _VIRUS_INCLUDE)
    reports = tuple(value.upper() for value in aux_reports)
    if any(value not in _VIRUS_AUX for value in reports):
        raise ConfigurationError(
            "aux_reports contains an unsupported virus report.",
            context={"allowed": sorted(_VIRUS_AUX)},
        )
    for name, value in (
        ("refseq_only", refseq_only),
        ("annotated_only", annotated_only),
        ("complete_only", complete_only),
    ):
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be boolean.")
    path_values = values if by == "accession" else (values[0],)
    url = build_url(
        api_base_url,
        f"/virus/{by}/{quoted_path_values(path_values)}/genome/download",
        (
            ("include_sequence", includes),
            ("aux_report", reports),
            ("refseq_only", refseq_only),
            ("annotated_only", annotated_only),
            ("complete_only", complete_only),
            ("host", None if host is None else require_text(host, "host", max_length=256)),
            (
                "geo_location",
                None
                if geo_location is None
                else require_text(geo_location, "geo_location", max_length=256),
            ),
            (
                "released_since",
                None
                if released_since is None
                else require_text(released_since, "released_since", max_length=64),
            ),
            (
                "updated_since",
                None
                if updated_since is None
                else require_text(updated_since, "updated_since", max_length=64),
            ),
        ),
    )
    headers = {} if resolved.api_key is None else {"api-key": resolved.api_key}
    return dataset(
        (RemoteFile(url, filename=f"ncbi_virus_{_slug(values)}.zip"),),
        output_dir,
        kind="ncbi_virus_package",
        source="NCBI Datasets Virus",
        metadata={
            "query": values,
            "by": by,
            "include": includes,
            "aux_reports": reports,
            "refseq_only": refseq_only,
            "annotated_only": annotated_only,
            "complete_only": complete_only,
            "host": host,
            "geo_location": geo_location,
            "released_since": released_since,
            "updated_since": updated_since,
        },
        config=resolved,
        progress=progress,
        headers=headers,
    )


def sequence(
    accession: str,
    output_path: str | Path,
    *,
    database: Literal["nuccore", "protein"] = "nuccore",
    format: Literal["fasta", "genbank", "embl"] = "fasta",
    seq_start: int | None = None,
    seq_stop: int | None = None,
    strand: Literal[1, 2] = 1,
    config: DownloadConfig | None = None,
    search_config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
    verbose: bool = True,
) -> DatasetDownloadResult:
    """Download a nucleotide/RNA/protein accession via NCBI EFetch, with a manifest.

    FASTA and GenBank/GenPept are native responses. Nucleotide EMBL output is
    converted from GenBank using optional ``dnakit[bio-formats]`` (Biopython).
    Optional seq_start/seq_stop use NCBI's 1-based inclusive coordinates;
    strand=2 retrieves the reverse complement. Protein requests use full records.
    SearchConfig controls contact/API credentials and bounded response size
    (20 MB by default, at most 100 MB); DownloadConfig controls output limits
    and overwrite. DownloadConfig.api_key is used when search_config is absent.
    """
    resolved = resolved_config(config)
    identifier = sequence_accession(accession)
    if database not in {"nuccore", "protein"} or format not in {"fasta", "genbank", "embl"}:
        raise ConfigurationError("Unsupported sequence database or format.")
    if isinstance(strand, bool) or strand not in (1, 2):
        raise ConfigurationError("strand must be 1 or 2.")
    if (seq_start is None) != (seq_stop is None):
        raise ConfigurationError("Provide both seq_start and seq_stop.")
    if seq_start is not None and (
        isinstance(seq_start, bool)
        or isinstance(seq_stop, bool)
        or not isinstance(seq_start, int)
        or not isinstance(seq_stop, int)
        or not 1 <= seq_start <= seq_stop
    ):
        raise ConfigurationError("Sequence coordinates must satisfy 1 <= seq_start <= seq_stop.")
    if database == "protein" and (format == "embl" or seq_start is not None or strand != 1):
        raise ConfigurationError(
            "Protein supports full FASTA or GenPept (format='genbank') records."
        )
    source_format = "fasta" if format == "fasta" else "genbank"
    converter = format_converter(source_format, format)
    target, manifest = sequence_target(output_path, resolved)
    query_config = search_config or SearchConfig(
        timeout=min(resolved.timeout, 300),
        max_response_bytes=min(resolved.max_file_bytes, 20_000_000),
        api_key=resolved.api_key,
    )
    rettype = "fasta" if format == "fasta" else ("gp" if database == "protein" else "gbwithparts")
    url = build_url(
        api_base_url,
        "/efetch.fcgi",
        [
            ("db", database),
            ("id", identifier),
            ("rettype", rettype),
            ("retmode", "text"),
            ("complexity", 1),
            ("seq_start", seq_start),
            ("seq_stop", seq_stop),
            ("strand", strand if database == "nuccore" else None),
            *_eutils_parameters(query_config),
        ],
    )
    text = request_text(url, query_config, provider="NCBI Entrez")
    return save_sequence(
        text,
        accession=identifier,
        database=database,
        provider="NCBI Entrez",
        url=url,
        source_format=source_format,
        format=format,
        converter=converter,
        target=target,
        manifest=manifest,
        config=resolved,
        metadata={
            "seq_start": seq_start,
            "seq_stop": seq_stop,
            "strand": strand,
            "coordinate_system": "1-based inclusive",
        },
    )


__all__ = [
    "annotation",
    "gene",
    "genome",
    "genome_package",
    "sequence",
    "taxonomy",
    "virus_package",
]

_register_public_verbose(globals())
