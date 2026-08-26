"""NCBI Datasets and Entrez E-utilities query adapters."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any, Literal
from urllib.parse import quote

from dnakit.exceptions import ConfigurationError, QueryError

from ._http import (
    build_url,
    limited_records,
    mapping_records,
    redact_url,
    request_json,
    require_text,
)
from ._shared import (
    adapter_provenance,
    query_values,
    quoted_path_values,
    resolved_config,
)
from ._shared import (
    page_size as validate_page_size,
)
from .models import QueryResult, SearchConfig

DEFAULT_NCBI_DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$", flags=re.IGNORECASE)
_REFSEQ_ACCESSION = re.compile(r"^[A-Z]{1,6}_\d+(?:\.\d+)?$", flags=re.IGNORECASE)
_VIRUS_ACCESSION = re.compile(r"^[A-Z]{1,6}_?\d{3,}(?:\.\d+)?$", flags=re.IGNORECASE)
_ENTREZ_DATABASES = frozenset(
    {
        "assembly",
        "bioproject",
        "biosample",
        "clinvar",
        "gds",
        "gene",
        "geoprofiles",
        "nuccore",
        "protein",
        "pubmed",
        "snp",
        "sra",
    }
)


def _datasets_headers(config: SearchConfig) -> dict[str, str]:
    return {} if config.api_key is None else {"api-key": config.api_key}


def _page_metadata(payload: Mapping[str, Any]) -> tuple[int | None, str | None]:
    raw_total = payload.get("total_count")
    total = raw_total if isinstance(raw_total, int) and not isinstance(raw_total, bool) else None
    token = payload.get("next_page_token")
    next_token = token if isinstance(token, str) and token.strip() else None
    return total, next_token


def taxonomy(
    query: str | int | Sequence[str | int],
    *,
    report: Literal["summary", "names", "related"] = "summary",
    children: bool = False,
    include_lineage: bool = True,
    ranks: Sequence[str] = (),
    page_size: int = 20,
    page_token: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
) -> QueryResult:
    """Query NCBI Taxonomy by TaxID, scientific name, common name, or alias.

    Use ``report="names"`` for the complete names/aliases report.  Use
    ``report="related"`` with one numeric TaxID for child nodes and lineage.
    """

    resolved = resolved_config(config)
    size = validate_page_size(page_size, resolved)
    values = query_values(query, maximum=100)
    if report not in {"summary", "names", "related"}:
        raise ConfigurationError("Unsupported taxonomy report.", code="INVALID_QUERY")
    params: list[tuple[str, object]] = [("page_size", size)]
    if page_token is not None:
        params.append(("page_token", require_text(page_token, "page_token")))
    if ranks:
        params.append(
            ("ranks", tuple(require_text(rank, "rank", max_length=128) for rank in ranks))
        )
    if report == "related":
        if len(values) != 1 or not values[0].isdigit():
            raise ConfigurationError(
                "taxonomy report='related' requires one numeric TaxID.", code="INVALID_QUERY"
            )
        path = f"/taxonomy/taxon/{values[0]}/related_ids"
        params.append(("include_lineage", include_lineage))
    else:
        suffix = "name_report" if report == "names" else "dataset_report"
        path = f"/taxonomy/taxon/{quoted_path_values(values)}/{suffix}"
        params.append(("children", children))
    url = build_url(api_base_url, path, params)
    payload = request_json(
        url,
        resolved,
        provider="NCBI Datasets",
        headers=_datasets_headers(resolved),
    )
    if not isinstance(payload, Mapping):
        raise QueryError(
            "NCBI Taxonomy returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
        )
    records: tuple[dict[str, object], ...]
    if report == "related":
        raw_ids = payload.get("tax_ids")
        records = (
            tuple(
                {"tax_id": item}
                for item in raw_ids
                if isinstance(item, int) and not isinstance(item, bool)
            )
            if isinstance(raw_ids, list)
            else ()
        )
    else:
        records = mapping_records(payload, key="reports")
    records = limited_records(records, resolved)
    total, next_token = _page_metadata(payload)
    total = max(total, len(records)) if total is not None else len(records)
    return QueryResult(
        "taxonomy",
        "NCBI Datasets",
        redact_url(url),
        records,
        adapter_provenance(
            "NCBI Taxonomy",
            citation_url="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/",
            filters={
                "query": values,
                "report": report,
                "children": children,
                "ranks": tuple(ranks),
            },
        ),
        total_count=total,
        next_page_token=next_token,
        metadata={"coordinate_system": None, "report": report},
    )


def assembly(
    query: str | Sequence[str],
    *,
    by: Literal["auto", "accession", "taxon"] = "auto",
    reference_only: bool = False,
    assembly_source: Literal["all", "refseq", "genbank"] = "all",
    current_only: bool = True,
    page_size: int = 20,
    page_token: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
) -> QueryResult:
    """Query genome assemblies and their statistics from NCBI Datasets."""

    resolved = resolved_config(config)
    size = validate_page_size(page_size, resolved)
    values = query_values(query, maximum=100)
    if by == "auto":
        accession_flags = tuple(
            _ASSEMBLY_ACCESSION.fullmatch(value) is not None for value in values
        )
        if any(accession_flags) and not all(accession_flags):
            raise ConfigurationError(
                "Do not mix assembly accessions and taxon names in one query.", code="INVALID_QUERY"
            )
        by = "accession" if all(accession_flags) else "taxon"
    if by not in {"accession", "taxon"}:
        raise ConfigurationError(
            "assembly by must be auto, accession, or taxon.", code="INVALID_QUERY"
        )
    path = f"/genome/{by}/{quoted_path_values(values)}/dataset_report"
    params: list[tuple[str, object]] = [
        ("filters.reference_only", reference_only),
        ("filters.assembly_source", assembly_source),
        ("filters.assembly_version", "current" if current_only else "all_assemblies"),
        ("page_size", size),
    ]
    if page_token is not None:
        params.append(("page_token", require_text(page_token, "page_token")))
    url = build_url(api_base_url, path, params)
    payload = request_json(
        url, resolved, provider="NCBI Datasets", headers=_datasets_headers(resolved)
    )
    if not isinstance(payload, Mapping):
        raise QueryError(
            "NCBI Assembly returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
        )
    records = limited_records(mapping_records(payload, key="reports"), resolved)
    total, next_token = _page_metadata(payload)
    total = max(total, len(records)) if total is not None else len(records)
    return QueryResult(
        "assembly",
        "NCBI Datasets",
        redact_url(url),
        records,
        adapter_provenance(
            "NCBI Genome",
            citation_url="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/get-genome-metadata/",
            filters={
                "query": values,
                "by": by,
                "reference_only": reference_only,
                "assembly_source": assembly_source,
                "current_only": current_only,
            },
        ),
        total_count=total,
        next_page_token=next_token,
    )


def gene(
    query: str | int | Sequence[str | int],
    *,
    taxon: str | int | None = None,
    by: Literal["auto", "id", "accession", "symbol", "search"] = "auto",
    report: Literal["gene", "product"] = "gene",
    gene_types: Sequence[str] = (),
    page_size: int = 20,
    page_token: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
) -> QueryResult:
    """Query genes by GeneID, RefSeq accession, symbol, alias, or name."""

    resolved = resolved_config(config)
    size = validate_page_size(page_size, resolved)
    values = query_values(query, maximum=100)
    taxon_text = None if taxon is None else require_text(str(taxon), "taxon", max_length=1_024)
    if report not in {"gene", "product"}:
        raise ConfigurationError("gene report must be gene or product.", code="INVALID_QUERY")
    if by == "auto":
        if all(value.isdigit() for value in values):
            by = "id"
        elif all(_REFSEQ_ACCESSION.fullmatch(value) is not None for value in values):
            by = "accession"
        elif taxon_text is not None:
            by = "symbol"
        else:
            raise ConfigurationError(
                "Gene symbols, aliases, and names require taxon; otherwise set by explicitly.",
                code="INVALID_QUERY",
            )
    suffix = "product_report" if report == "product" else "dataset_report"
    params: list[tuple[str, object]] = [("page_size", size)]
    if page_token is not None:
        params.append(("page_token", require_text(page_token, "page_token")))
    if gene_types:
        params.append(
            ("types", tuple(require_text(item, "gene_type", max_length=128) for item in gene_types))
        )
    if by == "id":
        if not all(value.isdigit() for value in values):
            raise ConfigurationError("GeneID queries must be numeric.", code="INVALID_QUERY")
        path = f"/gene/id/{quoted_path_values(values)}/{suffix}"
    elif by == "accession":
        path = f"/gene/accession/{quoted_path_values(values)}/{suffix}"
    elif by == "symbol":
        if taxon_text is None:
            raise ConfigurationError("Gene symbol queries require taxon.", code="INVALID_QUERY")
        path = (
            f"/gene/symbol/{quoted_path_values(values)}/taxon/{quote(taxon_text, safe='')}/{suffix}"
        )
    elif by == "search":
        if taxon_text is None or len(values) != 1:
            raise ConfigurationError(
                "Gene alias/name search requires one query and taxon.", code="INVALID_QUERY"
            )
        path = f"/gene/taxon/{quote(taxon_text, safe='')}/{suffix}"
        params.append(("query", values[0]))
    else:
        raise ConfigurationError("Unsupported gene query mode.", code="INVALID_QUERY")
    url = build_url(api_base_url, path, params)
    payload = request_json(
        url, resolved, provider="NCBI Datasets", headers=_datasets_headers(resolved)
    )
    if not isinstance(payload, Mapping):
        raise QueryError("NCBI Gene returned an unexpected response.", code="QUERY_RESPONSE_ERROR")
    records = limited_records(mapping_records(payload, key="reports"), resolved)
    total, next_token = _page_metadata(payload)
    total = max(total, len(records)) if total is not None else len(records)
    return QueryResult(
        "gene",
        "NCBI Datasets",
        redact_url(url),
        records,
        adapter_provenance(
            "NCBI Gene",
            citation_url="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/get-gene-metadata/",
            filters={"query": values, "taxon": taxon_text, "by": by, "report": report},
        ),
        total_count=total,
        next_page_token=next_token,
        metadata={"report": report},
    )


def ncbi_orthologs(
    gene_id: int | str,
    *,
    returned_content: Literal["complete", "ids_only"] = "complete",
    page_size: int = 20,
    page_token: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
) -> QueryResult:
    """Query the NCBI ortholog set for one numeric GeneID."""

    resolved = resolved_config(config)
    value = require_text(str(gene_id), "gene_id", max_length=32)
    if not value.isdigit():
        raise ConfigurationError("gene_id must be numeric.", code="INVALID_QUERY")
    size = validate_page_size(page_size, resolved)
    if returned_content not in {"complete", "ids_only"}:
        raise ConfigurationError("Unsupported ortholog returned_content.", code="INVALID_QUERY")
    params: list[tuple[str, object]] = [
        ("returned_content", returned_content.upper()),
        ("page_size", size),
    ]
    if page_token is not None:
        params.append(("page_token", require_text(page_token, "page_token")))
    url = build_url(api_base_url, f"/gene/id/{value}/orthologs", params)
    payload = request_json(
        url, resolved, provider="NCBI Datasets", headers=_datasets_headers(resolved)
    )
    if not isinstance(payload, Mapping):
        raise QueryError(
            "NCBI Orthologs returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
        )
    records = limited_records(mapping_records(payload, key="reports"), resolved)
    total, next_token = _page_metadata(payload)
    total = max(total, len(records)) if total is not None else len(records)
    return QueryResult(
        "homology",
        "NCBI Datasets",
        redact_url(url),
        records,
        adapter_provenance(
            "NCBI Gene Orthologs",
            citation_url="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-set/",
            filters={"gene_id": value},
        ),
        total_count=total,
        next_page_token=next_token,
    )


def virus(
    query: str | Sequence[str],
    *,
    by: Literal["auto", "accession", "taxon"] = "auto",
    report: Literal["dataset", "annotation"] = "dataset",
    refseq_only: bool = False,
    annotated_only: bool = False,
    complete_only: bool = False,
    host: str | None = None,
    geo_location: str | None = None,
    released_since: str | None = None,
    updated_since: str | None = None,
    page_size: int = 20,
    page_token: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_DATASETS_URL,
) -> QueryResult:
    """Query NCBI Virus genome and annotation metadata by taxon or accession."""

    resolved = resolved_config(config)
    values = query_values(query, maximum=100)
    size = validate_page_size(page_size, resolved)
    if by == "auto":
        accession_flags = tuple(_VIRUS_ACCESSION.fullmatch(value) is not None for value in values)
        if any(accession_flags) and not all(accession_flags):
            raise ConfigurationError(
                "Do not mix virus accessions and taxon names.", code="INVALID_QUERY"
            )
        by = "accession" if all(accession_flags) else "taxon"
    if by not in {"accession", "taxon"}:
        raise ConfigurationError("virus by must be auto, accession, or taxon.")
    if by == "taxon" and len(values) != 1:
        raise ConfigurationError("Virus taxon queries accept one taxon at a time.")
    if report not in {"dataset", "annotation"}:
        raise ConfigurationError("virus report must be dataset or annotation.")
    for name, value in (
        ("refseq_only", refseq_only),
        ("annotated_only", annotated_only),
        ("complete_only", complete_only),
    ):
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be boolean.")
    suffix = "dataset_report" if report == "dataset" else "annotation_report"
    path_values = values if by == "accession" else (values[0],)
    path = f"/virus/{by}/{quoted_path_values(path_values)}/{suffix}"
    params: list[tuple[str, object]] = [
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
        ("page_size", size),
    ]
    if page_token is not None:
        params.append(("page_token", require_text(page_token, "page_token")))
    url = build_url(api_base_url, path, params)
    payload = request_json(
        url, resolved, provider="NCBI Datasets", headers=_datasets_headers(resolved)
    )
    if not isinstance(payload, Mapping):
        raise QueryError("NCBI Virus returned an unexpected response.", code="QUERY_RESPONSE_ERROR")
    records = limited_records(mapping_records(payload, key="reports"), resolved)
    total, next_token = _page_metadata(payload)
    total = max(total, len(records)) if total is not None else len(records)
    return QueryResult(
        "virus",
        "NCBI Datasets",
        redact_url(url),
        records,
        adapter_provenance(
            "NCBI Virus",
            citation_url="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/virus-metadata/",
            filters={
                "query": values,
                "by": by,
                "report": report,
                "refseq_only": refseq_only,
                "annotated_only": annotated_only,
                "complete_only": complete_only,
                "host": host,
                "geo_location": geo_location,
                "released_since": released_since,
                "updated_since": updated_since,
            },
        ),
        total_count=total,
        next_page_token=next_token,
        metadata={"report": report},
    )


def _eutils_parameters(config: SearchConfig) -> list[tuple[str, object]]:
    parameters: list[tuple[str, object]] = [("tool", config.tool)]
    if config.email is not None:
        parameters.append(("email", config.email))
    if config.api_key is not None:
        parameters.append(("api_key", config.api_key))
    return parameters


def entrez(
    database: str,
    term: str,
    *,
    retmax: int = 20,
    retstart: int = 0,
    sort: str | None = None,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
    query_type: str | None = None,
) -> QueryResult:
    """Run a bounded ESearch followed by ESummary against an Entrez database."""

    resolved = resolved_config(config)
    db = require_text(database, "database", max_length=64).casefold()
    if db not in _ENTREZ_DATABASES:
        raise ConfigurationError(
            "Unsupported Entrez database.",
            code="INVALID_QUERY",
            context={"database": db, "allowed": sorted(_ENTREZ_DATABASES)},
        )
    query = require_text(term, "term")
    size = validate_page_size(retmax, resolved, maximum=10_000)
    if isinstance(retstart, bool) or not isinstance(retstart, int) or retstart < 0:
        raise ConfigurationError("retstart must be a non-negative integer.", code="INVALID_QUERY")
    search_params: list[tuple[str, object]] = [
        ("db", db),
        ("term", query),
        ("retmode", "json"),
        ("retmax", size),
        ("retstart", retstart),
        *_eutils_parameters(resolved),
    ]
    if sort is not None:
        search_params.append(("sort", require_text(sort, "sort", max_length=128)))
    search_url = build_url(api_base_url, "/esearch.fcgi", search_params)
    search_payload = request_json(search_url, resolved, provider="NCBI Entrez")
    if not isinstance(search_payload, Mapping) or not isinstance(
        search_payload.get("esearchresult"), Mapping
    ):
        raise QueryError(
            "NCBI ESearch returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
        )
    search_result = search_payload["esearchresult"]
    raw_ids = search_result.get("idlist")
    ids = tuple(str(item) for item in raw_ids) if isinstance(raw_ids, list) else ()
    if len(ids) > resolved.max_records:
        raise QueryError("NCBI ESearch exceeded max_records.", code="QUERY_RECORD_LIMIT")
    count_raw = search_result.get("count")
    try:
        total = int(count_raw)
    except (TypeError, ValueError):
        total = len(ids)
    records: tuple[dict[str, object], ...] = ()
    summary_url: str | None = None
    if ids:
        summary_url = build_url(
            api_base_url,
            "/esummary.fcgi",
            [
                ("db", db),
                ("id", ",".join(ids)),
                ("retmode", "json"),
                ("version", "2.0"),
                *_eutils_parameters(resolved),
            ],
        )
        summary_payload = request_json(summary_url, resolved, provider="NCBI Entrez")
        result = summary_payload.get("result") if isinstance(summary_payload, Mapping) else None
        if not isinstance(result, Mapping):
            raise QueryError(
                "NCBI ESummary returned an unexpected response.", code="QUERY_RESPONSE_ERROR"
            )
        uid_order = result.get("uids")
        ordered = tuple(str(item) for item in uid_order) if isinstance(uid_order, list) else ids
        records = tuple(
            dict(item) for uid in ordered if isinstance((item := result.get(uid)), Mapping)
        )
        records = limited_records(records, resolved)
    kind = query_type or db
    return QueryResult(
        kind,
        "NCBI Entrez",
        redact_url(search_url),
        records,
        adapter_provenance(
            f"NCBI Entrez {db}",
            citation_url="https://www.ncbi.nlm.nih.gov/books/NBK25501/",
            filters={"database": db, "term": query, "retstart": retstart, "retmax": size},
        ),
        total_count=max(total, len(records)),
        metadata={
            "database": db,
            "uids": ids,
            "summary_url": None if summary_url is None else redact_url(summary_url),
        },
    )


def accession(
    identifier: str,
    *,
    database: Literal["nuccore", "protein", "assembly", "sra"] = "nuccore",
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Check whether an accession exists and retrieve its current summary/version."""

    value = require_text(identifier, "identifier", max_length=256)
    return entrez(
        database,
        f'"{value}"[Accession]',
        retmax=20,
        config=config,
        api_base_url=api_base_url,
        query_type="accession",
    )


def project(
    term: str,
    *,
    retmax: int = 20,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Search BioProject metadata with the Entrez query language."""

    return entrez(
        "bioproject",
        term,
        retmax=retmax,
        config=config,
        api_base_url=api_base_url,
        query_type="project",
    )


def sample(
    term: str,
    *,
    retmax: int = 20,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Search BioSample metadata with the Entrez query language."""

    return entrez(
        "biosample",
        term,
        retmax=retmax,
        config=config,
        api_base_url=api_base_url,
        query_type="sample",
    )


def clinical_variant(
    term: str,
    *,
    retmax: int = 20,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Search ClinVar by rsID, HGVS, gene, condition, or Entrez expression."""

    return entrez(
        "clinvar",
        term,
        retmax=retmax,
        config=config,
        api_base_url=api_base_url,
        query_type="clinical_variant",
    )


def expression(
    term: str,
    *,
    retmax: int = 20,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Search GEO DataSets/Series metadata through NCBI E-utilities."""

    return entrez(
        "gds",
        term,
        retmax=retmax,
        config=config,
        api_base_url=api_base_url,
        query_type="expression",
    )


def literature(
    term: str,
    *,
    retmax: int = 20,
    sort: str = "relevance",
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Search linked biological literature in PubMed."""

    resolved = resolved_config(config)
    return entrez(
        "pubmed",
        term,
        retmax=retmax,
        sort=sort,
        config=resolved,
        api_base_url=api_base_url,
        query_type="literature",
    )


def database_version(
    database: str,
    *,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_NCBI_EUTILS_URL,
) -> QueryResult:
    """Return Entrez database record counts, indexed fields, and last-update metadata."""

    resolved = resolved_config(config)
    db = require_text(database, "database", max_length=64).casefold()
    if db not in _ENTREZ_DATABASES:
        raise ConfigurationError("Unsupported Entrez database.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        "/einfo.fcgi",
        [("db", db), ("retmode", "json"), ("version", "2.0"), *_eutils_parameters(resolved)],
    )
    payload = request_json(url, resolved, provider="NCBI Entrez")
    if not isinstance(payload, Mapping):
        raise QueryError("NCBI EInfo returned an unexpected response.", code="QUERY_RESPONSE_ERROR")
    root = payload.get("einforesult")
    raw_info = root.get("dbinfo") if isinstance(root, Mapping) else None
    if isinstance(raw_info, list):
        records = mapping_records(raw_info)
    elif isinstance(raw_info, Mapping):
        records = (dict(raw_info),)
    else:
        records = ()
    records = limited_records(records, resolved)
    return QueryResult(
        "database_version",
        "NCBI Entrez",
        redact_url(url),
        records,
        adapter_provenance(
            f"NCBI Entrez {db}",
            citation_url="https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EInfo",
            filters={"database": db},
        ),
        total_count=len(records),
    )


__all__ = [
    "DEFAULT_NCBI_DATASETS_URL",
    "DEFAULT_NCBI_EUTILS_URL",
    "accession",
    "assembly",
    "clinical_variant",
    "database_version",
    "entrez",
    "expression",
    "gene",
    "literature",
    "ncbi_orthologs",
    "project",
    "sample",
    "taxonomy",
    "virus",
]
