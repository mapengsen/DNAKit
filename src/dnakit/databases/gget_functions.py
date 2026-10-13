"""gget wrappers: offline local COSMIC lookup plus gene- and sequence-level database queries.

Offline: ``gget_cosmic_local`` searches a COSMIC mutation TSV already on disk.
Network: the remaining functions call public gget services (NCBI, UCSC, OpenTargets, ARCHS4,
Bgee, Enrichr, 8cubeDB). They are explicit adapters; tests replace the upstream call.

Functions already wrapped in ``gene_backends`` (``gget.search``, ``gget.info``, ``gget.seq``,
``gget.mutate``) are not repeated here. Packages import only inside each function.
"""

from __future__ import annotations

import importlib.metadata
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_GGET = "bio-gget"
_COSMIC_PROJECTS = (
    "cancer",
    "cancer_example",
    "census",
    "resistance",
    "cell_line",
    "genome_screen",
    "targeted_screen",
    "other",
)
_BLAST_PROGRAMS = ("default", "blastn", "blastp", "blastx", "tblastn", "tblastx")
_BLAT_SEQTYPES = ("default", "DNA", "protein", "translated RNA", "translated DNA")
_OPENTARGETS_RESOURCES = (
    "depmap",
    "diseases",
    "drugs",
    "expression",
    "interactions",
    "pharmacogenetics",
    "tractability",
)
_ARCHS4_WHICH = ("correlation", "tissue")
_ARCHS4_SPECIES = ("human", "mouse")
_BGEE_TYPES = ("orthologs", "expression")
_ENRICHR_SPECIES = ("human", "mouse")
_ENSEMBL_GENE = re.compile(r"ENS[A-Z]{0,6}G\d{11}")
_GENE_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,63}")
_TEXT_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9 ._:+/-]{0,127}")
_SEQUENCE = re.compile(r"[A-Za-z*-]+")
_MAX_GENES = 500
_MAX_SEQUENCE = 20_000


def _version() -> str:
    try:
        return importlib.metadata.version("gget")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - source checkouts only
        return "unknown"


def _provenance(**notes: object) -> dict[str, Any]:
    return {"package": "gget", "version": _version(), **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "complete", "operation": operation, "provenance": provenance, **fields}


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > limit:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {limit} chars.")
    return value.strip()


def _token(value: object, name: str, pattern: re.Pattern[str]) -> str:
    text = _text(value, name, limit=128)
    if pattern.fullmatch(text) is None:
        raise ConfigurationError(f"{name} contains unsupported characters: {text!r}.")
    return text


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _tokens(values: object, name: str, pattern: re.Pattern[str], *, limit: int) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} accepts at most {limit} entries.")
    return [_token(item, name, pattern) for item in values]


def _sequence(value: object, name: str) -> str:
    text = _text(value, name, limit=_MAX_SEQUENCE)
    if _SEQUENCE.fullmatch(text) is None:
        raise ConfigurationError(f"{name} must contain only sequence letters, '*' or '-'.")
    return text.upper()


def _plain(value: Any) -> Any:
    """Make upstream results JSON-serializable and free of NaN/NumPy scalars."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "to_dict"):
        return _plain(value.to_dict(orient="records"))
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    return str(value)


def _records(value: Any) -> list[Any]:
    plain = _plain(value)
    if plain is None:
        return []
    if isinstance(plain, dict):
        return [plain]
    return list(plain)


def _upstream(name: str, call: Any) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"gget {name} failed: {exc}") from exc


# --- offline ----------------------------------------------------------------------------


def gget_cosmic_local(
    searchterm: str,
    cosmic_tsv_path: str | Path,
    *,
    cosmic_project: str | None = None,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a local COSMIC mutation TSV through gget.cosmic (offline; no login or download)."""

    gget = require("gget", _GGET)
    term = _text(searchterm, "searchterm", limit=128)
    path = input_path(cosmic_tsv_path)
    count = _int(limit, "limit", 1, 10_000)
    flag = _flag(verbose, "verbose")
    if cosmic_project is None:
        # Mirror gget's own inference so the recorded project matches what gget applies.
        project = "cancer" if "CancerMutationCensus_AllData" in str(path) else "other"
    else:
        project = _choice(cosmic_project, "cosmic_project", _COSMIC_PROJECTS)

    def call() -> Any:
        try:
            return gget.cosmic(
                term,
                cosmic_tsv_path=str(path),
                cosmic_project=project,
                limit=count,
                json=True,
                verbose=flag,
            )
        except ValueError as exc:
            if "No results were found" in str(exc):
                return []
            raise

    rows = _records(_upstream("cosmic", call))
    return _document(
        "gget_cosmic_local",
        _provenance(source="local COSMIC mutation TSV", network=False),
        searchterm=term,
        cosmic_project=project,
        source_file=path.name,
        limit=count,
        count=len(rows),
        records=rows,
    )


# --- network: sequence alignment ----------------------------------------------------------


def gget_blast(
    sequence: str,
    *,
    program: str = "default",
    database: str = "default",
    limit: int = 50,
    expect: float = 10.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """BLAST one nucleotide or protein sequence through gget.blast (NCBI network)."""

    gget = require("gget", _GGET)
    query = _sequence(sequence, "sequence")
    chosen_program = _choice(program, "program", _BLAST_PROGRAMS)
    chosen_db = _token(database, "database", re.compile(r"[A-Za-z0-9_.-]{1,64}|default"))
    hits_limit = _int(limit, "limit", 1, 500)
    threshold = _number(expect, "expect", 1e-30, 1000.0)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "blast",
        lambda: gget.blast(
            query,
            program=chosen_program,
            database=chosen_db,
            limit=hits_limit,
            expect=threshold,
            json=True,
            verbose=flag,
        ),
    )
    hits = _records(rows)
    return _document(
        "gget_blast",
        _provenance(source="NCBI BLAST via gget"),
        query_length=len(query),
        program=chosen_program,
        database=chosen_db,
        expect=threshold,
        count=len(hits),
        hits=hits,
    )


def gget_blat(
    sequence: str,
    *,
    seqtype: str = "default",
    assembly: str = "human",
    verbose: bool = True,
) -> dict[str, Any]:
    """BLAT one DNA or protein sequence against a UCSC assembly through gget.blat (network)."""

    gget = require("gget", _GGET)
    query = _sequence(sequence, "sequence")
    kind = _choice(seqtype, "seqtype", _BLAT_SEQTYPES)
    build = _token(assembly, "assembly", re.compile(r"[A-Za-z0-9_]{1,32}"))
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "blat",
        lambda: gget.blat(query, seqtype=kind, assembly=build, json=True, verbose=flag),
    )
    hits = _records(rows)
    return _document(
        "gget_blat",
        _provenance(source="UCSC BLAT via gget"),
        query_length=len(query),
        seqtype=kind,
        assembly=build,
        count=len(hits),
        hits=hits,
    )


# --- network: gene-level annotation -------------------------------------------------------


def gget_opentargets(
    ensembl_id: str,
    *,
    resource: str = "diseases",
    limit: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Query OpenTargets for diseases, drugs or other records of one Ensembl gene (network)."""

    gget = require("gget", _GGET)
    gene = _token(ensembl_id, "ensembl_id", _ENSEMBL_GENE)
    chosen = _choice(resource, "resource", _OPENTARGETS_RESOURCES)
    cap = None if limit is None else _int(limit, "limit", 1, 10_000)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "opentargets",
        lambda: gget.opentargets(gene, resource=chosen, limit=cap, json=True, verbose=flag),
    )
    records = _records(rows)
    return _document(
        "gget_opentargets",
        _provenance(source="OpenTargets Platform GraphQL via gget"),
        ensembl_id=gene,
        resource=chosen,
        count=len(records),
        records=records,
    )


def gget_archs4(
    gene: str,
    *,
    ensembl: bool = False,
    which: str = "correlation",
    gene_count: int = 100,
    species: str = "human",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return ARCHS4 co-expressed genes or tissue expression for one gene (network)."""

    gget = require("gget", _GGET)
    symbol = _token(gene, "gene", _GENE_TOKEN)
    use_ensembl = _flag(ensembl, "ensembl")
    if use_ensembl and _ENSEMBL_GENE.fullmatch(symbol) is None:
        raise ConfigurationError("ensembl=True requires an Ensembl gene ID.")
    mode = _choice(which, "which", _ARCHS4_WHICH)
    count = _int(gene_count, "gene_count", 1, 1000)
    taxon = _choice(species, "species", _ARCHS4_SPECIES)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "archs4",
        lambda: gget.archs4(
            symbol,
            ensembl=use_ensembl,
            which=mode,
            gene_count=count,
            species=taxon,
            json=True,
            verbose=flag,
        ),
    )
    records = _records(rows)
    return _document(
        "gget_archs4",
        _provenance(source="ARCHS4 via gget"),
        gene=symbol,
        which=mode,
        species=taxon,
        count=len(records),
        records=records,
    )


def gget_bgee(
    gene_id: str | Sequence[str],
    *,
    record_type: str = "orthologs",
    verbose: bool = True,
) -> dict[str, Any]:
    """Get Bgee orthologs (one gene) or expression calls (several genes) through gget (network)."""

    gget = require("gget", _GGET)
    kind = _choice(record_type, "record_type", _BGEE_TYPES)
    if isinstance(gene_id, str):
        ids = [_token(gene_id, "gene_id", _ENSEMBL_GENE)]
    else:
        ids = _tokens(gene_id, "gene_id", _ENSEMBL_GENE, limit=100)
    if kind == "orthologs" and len(ids) != 1:
        raise ConfigurationError("record_type='orthologs' accepts exactly one gene_id.")
    flag = _flag(verbose, "verbose")
    query: Any = ids[0] if kind == "orthologs" else ids
    rows = _upstream(
        "bgee",
        lambda: gget.bgee(query, type=kind, json=True, verbose=flag),
    )
    records = _records(rows)
    return _document(
        "gget_bgee",
        _provenance(source="Bgee via gget"),
        gene_ids=ids,
        record_type=kind,
        count=len(records),
        records=records,
    )


def gget_enrichr(
    genes: Sequence[str],
    *,
    database: str,
    species: str = "human",
    background_list: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run Enrichr gene-set enrichment for a gene list against one library (network)."""

    gget = require("gget", _GGET)
    symbols = _tokens(genes, "genes", _GENE_TOKEN, limit=_MAX_GENES)
    library = _token(database, "database", _TEXT_TOKEN)
    taxon = _choice(species, "species", _ENRICHR_SPECIES)
    background = None
    if background_list is not None:
        background = _tokens(background_list, "background_list", _GENE_TOKEN, limit=50_000)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "enrichr",
        lambda: gget.enrichr(
            symbols,
            database=library,
            species=taxon,
            background_list=background,
            json=True,
            verbose=flag,
        ),
    )
    records = _records(rows)
    return _document(
        "gget_enrichr",
        _provenance(source="Enrichr via gget"),
        gene_count=len(symbols),
        database=library,
        species=taxon,
        count=len(records),
        records=records,
    )


# --- network: 8cubeDB expression ----------------------------------------------------------


def gget_specificity(gene_list: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Retrieve gene-level tissue or cell specificity statistics from 8cubeDB (network)."""

    gget = require("gget", _GGET)
    genes = _tokens(gene_list, "gene_list", _GENE_TOKEN, limit=_MAX_GENES)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "specificity",
        lambda: gget.specificity(genes, json=True, verbose=flag),
    )
    records = _records(rows)
    return _document(
        "gget_specificity",
        _provenance(source="8cubeDB via gget"),
        gene_count=len(genes),
        count=len(records),
        records=records,
    )


def gget_gene_expression(
    gene_list: Sequence[str],
    *,
    analysis_level: str,
    analysis_type: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve normalized 8cubeDB expression values for genes in one partition (network)."""

    gget = require("gget", _GGET)
    genes = _tokens(gene_list, "gene_list", _GENE_TOKEN, limit=_MAX_GENES)
    level = _token(analysis_level, "analysis_level", _TEXT_TOKEN)
    partition = _token(analysis_type, "analysis_type", _TEXT_TOKEN)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "gene_expression",
        lambda: gget.gene_expression(
            genes,
            analysis_level=level,
            analysis_type=partition,
            json=True,
            verbose=flag,
        ),
    )
    records = _records(rows)
    return _document(
        "gget_gene_expression",
        _provenance(source="8cubeDB normalized expression via gget"),
        gene_count=len(genes),
        analysis_level=level,
        analysis_type=partition,
        count=len(records),
        records=records,
    )


def gget_psi_block(
    gene_list: Sequence[str],
    *,
    analysis_level: str,
    analysis_type: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve 8cubeDB psi-block specificity scores for genes in one partition (network)."""

    gget = require("gget", _GGET)
    genes = _tokens(gene_list, "gene_list", _GENE_TOKEN, limit=_MAX_GENES)
    level = _token(analysis_level, "analysis_level", _TEXT_TOKEN)
    partition = _token(analysis_type, "analysis_type", _TEXT_TOKEN)
    flag = _flag(verbose, "verbose")
    rows = _upstream(
        "psi_block",
        lambda: gget.psi_block(
            genes,
            analysis_level=level,
            analysis_type=partition,
            json=True,
            verbose=flag,
        ),
    )
    records = _records(rows)
    return _document(
        "gget_psi_block",
        _provenance(source="8cubeDB psi-block scores via gget"),
        gene_count=len(genes),
        analysis_level=level,
        analysis_type=partition,
        count=len(records),
        records=records,
    )


__all__ = [
    "gget_archs4",
    "gget_bgee",
    "gget_blast",
    "gget_blat",
    "gget_cosmic_local",
    "gget_enrichr",
    "gget_gene_expression",
    "gget_opentargets",
    "gget_psi_block",
    "gget_specificity",
]

_register_public_verbose(globals(), tuple(__all__))
