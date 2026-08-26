"""European Nucleotide Archive Portal API query adapters."""

from __future__ import annotations

import re
from collections.abc import Sequence

from dnakit.exceptions import ConfigurationError, QueryError

from ._http import (
    build_url,
    limited_records,
    mapping_records,
    redact_url,
    request_json,
    require_text,
)
from ._shared import adapter_provenance, resolved_config
from .models import QueryResult, SearchConfig

DEFAULT_ENA_PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api"
_RESULT_TYPES = frozenset(
    {
        "analysis",
        "assembly",
        "read_experiment",
        "read_run",
        "sample",
        "sequence",
        "study",
        "taxon",
        "tsa_set",
        "wgs_set",
    }
)
_FIELD = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
DEFAULT_READ_FIELDS = (
    "accession",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "scientific_name",
    "tax_id",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "nominal_length",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
    "first_public",
    "last_updated",
)


def ena_search(
    result: str,
    query: str,
    *,
    fields: Sequence[str],
    limit: int = 100,
    offset: int = 0,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENA_PORTAL_URL,
    query_type: str | None = None,
) -> QueryResult:
    """Run a bounded ENA Portal query and return selected metadata fields."""

    resolved = resolved_config(config)
    result_type = require_text(result, "result", max_length=64)
    if result_type not in _RESULT_TYPES:
        raise ConfigurationError(
            "Unsupported ENA result type.",
            code="INVALID_QUERY",
            context={"allowed": sorted(_RESULT_TYPES)},
        )
    expression = require_text(query, "query")
    field_values = tuple(require_text(field, "field", max_length=128) for field in fields)
    if (
        not field_values
        or len(field_values) > 128
        or any(_FIELD.fullmatch(field) is None for field in field_values)
    ):
        raise ConfigurationError(
            "ENA fields must contain 1-128 valid field names.", code="INVALID_QUERY"
        )
    if (
        isinstance(limit, bool)
        or not isinstance(limit, int)
        or not 1 <= limit <= resolved.max_records
    ):
        raise ConfigurationError(
            "limit exceeds SearchConfig max_records.", code="QUERY_RECORD_LIMIT"
        )
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise ConfigurationError("offset must be a non-negative integer.", code="INVALID_QUERY")
    url = build_url(
        api_base_url,
        "/search",
        (
            ("result", result_type),
            ("query", expression),
            ("fields", ",".join(field_values)),
            ("format", "json"),
            ("limit", limit),
            ("offset", offset),
        ),
    )
    payload = request_json(url, resolved, provider="ENA")
    if not isinstance(payload, list):
        raise QueryError("ENA returned an unexpected response.", code="QUERY_RESPONSE_ERROR")
    records = limited_records(mapping_records(payload), resolved)
    return QueryResult(
        query_type or result_type,
        "ENA",
        redact_url(url),
        records,
        adapter_provenance(
            "European Nucleotide Archive",
            citation_url="https://www.ebi.ac.uk/ena/portal/api/",
            filters={
                "result": result_type,
                "query": expression,
                "fields": field_values,
                "limit": limit,
                "offset": offset,
            },
        ),
        total_count=len(records),
        metadata={
            "result_type": result_type,
            "offset": offset,
            "page_is_full": len(records) == limit,
        },
    )


def reads(
    query: str,
    *,
    fields: Sequence[str] = DEFAULT_READ_FIELDS,
    result: str = "read_run",
    limit: int = 100,
    offset: int = 0,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENA_PORTAL_URL,
) -> QueryResult:
    """Query sequencing Study/Experiment/Run metadata and public file locations."""

    return ena_search(
        result,
        query,
        fields=fields,
        limit=limit,
        offset=offset,
        config=config,
        api_base_url=api_base_url,
        query_type="reads",
    )


def project(
    query: str,
    *,
    fields: Sequence[str] = (
        "study_accession",
        "secondary_study_accession",
        "study_title",
        "study_type",
        "study_abstract",
        "center_name",
        "first_public",
        "last_updated",
    ),
    limit: int = 100,
    offset: int = 0,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENA_PORTAL_URL,
) -> QueryResult:
    """Search public ENA Study/BioProject records."""

    return ena_search(
        "study",
        query,
        fields=fields,
        limit=limit,
        offset=offset,
        config=config,
        api_base_url=api_base_url,
        query_type="project",
    )


def sample(
    query: str,
    *,
    fields: Sequence[str] = (
        "sample_accession",
        "secondary_sample_accession",
        "scientific_name",
        "tax_id",
        "description",
        "collection_date",
        "country",
        "location",
        "host",
        "host_tax_id",
        "first_public",
        "last_updated",
    ),
    limit: int = 100,
    offset: int = 0,
    config: SearchConfig | None = None,
    api_base_url: str = DEFAULT_ENA_PORTAL_URL,
) -> QueryResult:
    """Search public ENA Sample/BioSample metadata."""

    return ena_search(
        "sample",
        query,
        fields=fields,
        limit=limit,
        offset=offset,
        config=config,
        api_base_url=api_base_url,
        query_type="sample",
    )


__all__ = [
    "DEFAULT_ENA_PORTAL_URL",
    "DEFAULT_READ_FIELDS",
    "ena_search",
    "project",
    "reads",
    "sample",
]
