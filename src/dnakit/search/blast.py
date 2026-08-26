"""Asynchronous NCBI BLAST URL API adapter for sequence identification."""

from __future__ import annotations

import hashlib
import math
import re
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Literal, cast
from urllib.parse import urlencode

from dnakit.core import DNASequence, Provenance
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError, QueryError

from ._http import (
    build_url,
    redact_url,
    request_json,
    request_text,
    require_https_base,
    require_text,
)
from ._shared import adapter_provenance, resolved_config
from .models import QueryResult, SearchConfig

DEFAULT_NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
_RID = re.compile(r"^[A-Za-z0-9-]{8,64}$")
_DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", flags=re.IGNORECASE)
_PROGRAMS = frozenset({"blastn", "blastp", "blastx", "tblastn", "tblastx"})
_STATUS = re.compile(r"Status=(?P<status>[A-Z]+)")
_HITS = re.compile(r"ThereAreHits=(?P<hits>yes|no)", flags=re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class BlastJob:
    """An NCBI BLAST request identifier that can be polled explicitly."""

    rid: str
    estimated_seconds: int
    database: str
    program: str
    query_sha256: str
    request_url: str
    provenance: Provenance

    def __post_init__(self) -> None:
        if not isinstance(self.rid, str) or _RID.fullmatch(self.rid) is None:
            raise ConfigurationError("BlastJob rid is invalid.", code="INVALID_BLAST_JOB")
        if (
            isinstance(self.estimated_seconds, bool)
            or not isinstance(self.estimated_seconds, int)
            or self.estimated_seconds < 0
        ):
            raise ConfigurationError(
                "BlastJob estimated_seconds is invalid.", code="INVALID_BLAST_JOB"
            )
        if not re.fullmatch(r"[0-9a-f]{64}", self.query_sha256):
            raise ConfigurationError("BlastJob query_sha256 is invalid.", code="INVALID_BLAST_JOB")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("BlastJob provenance is invalid.", code="INVALID_BLAST_JOB")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class BlastStatus:
    """Current state of one submitted BLAST request."""

    rid: str
    status: Literal["waiting", "ready", "failed", "unknown"]
    has_hits: bool | None


@dataclass(frozen=True, slots=True)
class BlastProgress:
    """Polling progress emitted by ``wait_for_blast``."""

    rid: str
    status: str
    elapsed_seconds: float


def _endpoint_origin(url: str) -> tuple[str, str]:
    value = require_text(url, "api_url", max_length=2_048)
    if "?" in value or "#" in value:
        raise ConfigurationError(
            "BLAST api_url cannot contain query or fragment.", code="INVALID_QUERY"
        )
    slash = value.rfind("/")
    if slash <= len("https://"):
        raise ConfigurationError(
            "BLAST api_url must include an endpoint path.", code="INVALID_QUERY"
        )
    base = require_https_base(value[:slash])
    path = value[slash:]
    return base, path


def _sequence_text(sequence: str | DNASequence, program: str) -> str:
    if isinstance(sequence, DNASequence):
        if sequence.is_gapped:
            raise ConfigurationError(
                "BLAST query cannot contain structured gaps.", code="INVALID_QUERY"
            )
        value = sequence.symbols
    elif isinstance(sequence, str):
        value = "".join(sequence.split())
    else:
        raise TypeError("sequence must be str or DNASequence.")
    if not value or len(value) > 1_000_000:
        raise ConfigurationError(
            "BLAST query length must be in [1, 1000000].", code="QUERY_REQUEST_SIZE_LIMIT"
        )
    if program in {"blastn", "blastx", "tblastx"} and _DNA.fullmatch(value) is None:
        raise ConfigurationError(
            "Nucleotide BLAST queries must contain IUPAC DNA symbols only.", code="INVALID_QUERY"
        )
    return value.upper()


def _contact_parameters(config: SearchConfig) -> list[tuple[str, str]]:
    if config.email is None:
        raise ConfigurationError(
            "NCBI BLAST requires SearchConfig(email=...) for responsible API use.",
            code="BLAST_CONTACT_REQUIRED",
        )
    return [("TOOL", config.tool), ("EMAIL", config.email)]


def submit_blast(
    sequence: str | DNASequence,
    *,
    database: str = "core_nt",
    program: str = "blastn",
    hitlist_size: int = 50,
    expect: float = 10.0,
    megablast: bool = True,
    config: SearchConfig | None = None,
    api_url: str = DEFAULT_NCBI_BLAST_URL,
) -> BlastJob:
    """Submit one bounded BLAST query; no implicit polling or sleeping occurs."""

    resolved = resolved_config(config)
    program_value = require_text(program, "program", max_length=32).casefold()
    if program_value not in _PROGRAMS:
        raise ConfigurationError("Unsupported BLAST program.", code="INVALID_QUERY")
    database_value = require_text(database, "database", max_length=256)
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", database_value):
        raise ConfigurationError("BLAST database name is invalid.", code="INVALID_QUERY")
    if (
        isinstance(hitlist_size, bool)
        or not isinstance(hitlist_size, int)
        or not 1 <= hitlist_size <= min(5_000, resolved.max_records)
    ):
        raise ConfigurationError(
            "hitlist_size exceeds the configured record limit.", code="QUERY_RECORD_LIMIT"
        )
    if (
        isinstance(expect, bool)
        or not isinstance(expect, (int, float))
        or not math.isfinite(expect)
        or expect <= 0
    ):
        raise ConfigurationError("expect must be a positive finite number.", code="INVALID_QUERY")
    query = _sequence_text(sequence, program_value)
    base, path = _endpoint_origin(api_url)
    endpoint = build_url(base, path)
    parameters: list[tuple[str, object]] = [
        ("CMD", "Put"),
        ("PROGRAM", program_value),
        ("DATABASE", database_value),
        ("QUERY", query),
        ("HITLIST_SIZE", hitlist_size),
        ("EXPECT", expect),
        ("FORMAT_TYPE", "JSON2"),
        *_contact_parameters(resolved),
    ]
    if program_value == "blastn" and megablast:
        parameters.append(("MEGABLAST", "on"))
    body = urlencode(parameters).encode("utf-8")
    response = request_text(
        endpoint,
        resolved,
        provider="NCBI BLAST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
        data=body,
    )
    rid_match = re.search(r"^\s*RID\s*=\s*(\S+)\s*$", response, flags=re.MULTILINE)
    rtoe_match = re.search(r"^\s*RTOE\s*=\s*(\d+)\s*$", response, flags=re.MULTILINE)
    if rid_match is None or rtoe_match is None:
        raise QueryError("NCBI BLAST did not return RID/RTOE.", code="BLAST_SUBMISSION_ERROR")
    rid = rid_match.group(1)
    return BlastJob(
        rid=rid,
        estimated_seconds=int(rtoe_match.group(1)),
        database=database_value,
        program=program_value,
        query_sha256=hashlib.sha256(query.encode("ascii")).hexdigest(),
        request_url=endpoint,
        provenance=adapter_provenance(
            "NCBI BLAST",
            citation_url="https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html",
            filters={
                "database": database_value,
                "program": program_value,
                "hitlist_size": hitlist_size,
                "expect": float(expect),
                "megablast": megablast,
                "query_sha256": hashlib.sha256(query.encode("ascii")).hexdigest(),
            },
        ),
    )


def blast_status(
    job: BlastJob | str,
    *,
    config: SearchConfig | None = None,
    api_url: str = DEFAULT_NCBI_BLAST_URL,
) -> BlastStatus:
    """Check one RID without automatic retrying."""

    resolved = resolved_config(config)
    rid = job.rid if isinstance(job, BlastJob) else require_text(job, "rid", max_length=64)
    if _RID.fullmatch(rid) is None:
        raise ConfigurationError("rid is invalid.", code="INVALID_BLAST_JOB")
    base, path = _endpoint_origin(api_url)
    url = build_url(
        base,
        path,
        (
            ("CMD", "Get"),
            ("FORMAT_OBJECT", "SearchInfo"),
            ("RID", rid),
            *_contact_parameters(resolved),
        ),
    )
    response = request_text(url, resolved, provider="NCBI BLAST")
    status_match = _STATUS.search(response)
    raw_status = status_match.group("status").casefold() if status_match is not None else "unknown"
    status = raw_status if raw_status in {"waiting", "ready", "failed"} else "unknown"
    hits_match = _HITS.search(response)
    has_hits = None if hits_match is None else hits_match.group("hits").casefold() == "yes"
    return BlastStatus(
        rid, cast(Literal["waiting", "ready", "failed", "unknown"], status), has_hits
    )


def _blast_hits(
    payload: object, *, query_length_hint: int | None = None
) -> tuple[dict[str, object], ...]:
    if not isinstance(payload, Mapping):
        return ()
    outputs = payload.get("BlastOutput2")
    if not isinstance(outputs, list):
        return ()
    records: list[dict[str, object]] = []
    for output in outputs:
        if not isinstance(output, Mapping):
            continue
        report = output.get("report")
        results = report.get("results") if isinstance(report, Mapping) else None
        search = results.get("search") if isinstance(results, Mapping) else None
        if not isinstance(search, Mapping):
            continue
        query_length = search.get("query_len")
        if not isinstance(query_length, int) or query_length < 1:
            query_length = query_length_hint
        hits = search.get("hits")
        if not isinstance(hits, list):
            continue
        for hit in hits:
            if not isinstance(hit, Mapping):
                continue
            descriptions = hit.get("description")
            description = (
                descriptions[0]
                if isinstance(descriptions, list)
                and descriptions
                and isinstance(descriptions[0], Mapping)
                else {}
            )
            hsps = hit.get("hsps")
            valid_hsps = (
                tuple(item for item in hsps if isinstance(item, Mapping))
                if isinstance(hsps, list)
                else ()
            )
            best = max(
                valid_hsps,
                key=lambda item: (
                    float(item.get("bit_score", 0.0))
                    if isinstance(item.get("bit_score"), (int, float))
                    else 0.0
                ),
                default={},
            )
            align_len = best.get("align_len") if isinstance(best, Mapping) else None
            identities = best.get("identity") if isinstance(best, Mapping) else None
            identity = (
                float(identities) / align_len
                if isinstance(identities, int) and isinstance(align_len, int) and align_len > 0
                else None
            )
            coverage = (
                float(align_len) / query_length
                if isinstance(align_len, int) and isinstance(query_length, int) and query_length > 0
                else None
            )
            record: dict[str, object] = {
                "accession": description.get("accession")
                if isinstance(description, Mapping)
                else None,
                "title": description.get("title") if isinstance(description, Mapping) else None,
                "tax_id": description.get("taxid") if isinstance(description, Mapping) else None,
                "scientific_name": description.get("sciname")
                if isinstance(description, Mapping)
                else None,
                "identity": identity,
                "query_coverage": coverage,
                "evalue": best.get("evalue") if isinstance(best, Mapping) else None,
                "bit_score": best.get("bit_score") if isinstance(best, Mapping) else None,
                "alignment_length": align_len,
                "query_length": query_length,
                "hsp_count": len(valid_hsps),
            }
            records.append(record)
    return tuple(records)


def blast_results(
    job: BlastJob | str,
    *,
    config: SearchConfig | None = None,
    api_url: str = DEFAULT_NCBI_BLAST_URL,
) -> QueryResult:
    """Retrieve ready BLAST JSON2 results and normalize the best HSP per hit."""

    resolved = resolved_config(config)
    rid = job.rid if isinstance(job, BlastJob) else require_text(job, "rid", max_length=64)
    if _RID.fullmatch(rid) is None:
        raise ConfigurationError("rid is invalid.", code="INVALID_BLAST_JOB")
    base, path = _endpoint_origin(api_url)
    url = build_url(
        base,
        path,
        (("CMD", "Get"), ("RID", rid), ("FORMAT_TYPE", "JSON2"), *_contact_parameters(resolved)),
    )
    payload = request_json(url, resolved, provider="NCBI BLAST")
    hits = _blast_hits(payload)
    if len(hits) > resolved.max_records:
        raise QueryError("BLAST results exceeded max_records.", code="QUERY_RECORD_LIMIT")
    provenance = (
        job.provenance
        if isinstance(job, BlastJob)
        else adapter_provenance(
            "NCBI BLAST",
            citation_url="https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html",
            filters={"rid": rid},
        )
    )
    identity_values = tuple(
        float(value) for hit in hits if isinstance((value := hit.get("identity")), float)
    )
    max_identity = max(
        identity_values,
        default=None,
    )
    coverage_values = tuple(
        float(value) for hit in hits if isinstance((value := hit.get("query_coverage")), float)
    )
    max_coverage = max(
        coverage_values,
        default=None,
    )
    return QueryResult(
        "sequence_similarity",
        "NCBI BLAST",
        redact_url(url),
        hits,
        provenance,
        total_count=len(hits),
        metadata={
            "rid": rid,
            "max_identity": max_identity,
            "max_query_coverage": max_coverage,
            "novelty_score": None if max_identity is None else 1.0 - max_identity,
            "novelty_definition": "1 - maximum best-HSP identity; inspect coverage separately",
        },
    )


def wait_for_blast(
    job: BlastJob,
    *,
    poll_interval: float = 60.0,
    timeout: float = 1_800.0,
    progress: Callable[[BlastProgress], None] | None = None,
    config: SearchConfig | None = None,
    api_url: str = DEFAULT_NCBI_BLAST_URL,
) -> QueryResult:
    """Poll at NCBI-compliant intervals until a submitted BLAST job is ready."""

    if not isinstance(job, BlastJob):
        raise TypeError("job must be BlastJob.")
    if (
        isinstance(poll_interval, bool)
        or not isinstance(poll_interval, (int, float))
        or not math.isfinite(poll_interval)
        or poll_interval < 60
    ):
        raise ConfigurationError(
            "poll_interval must be at least 60 seconds for NCBI BLAST.", code="INVALID_QUERY"
        )
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or not 60 <= timeout <= 86_400
    ):
        raise ConfigurationError("timeout must be in [60, 86400].", code="INVALID_QUERY")
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.", code="INVALID_QUERY")
    resolved = resolved_config(config)
    started = time.monotonic()
    while True:
        elapsed = time.monotonic() - started
        if elapsed + poll_interval > timeout:
            raise QueryError(
                "NCBI BLAST job did not finish before timeout.",
                code="BLAST_TIMEOUT",
                context={"rid": job.rid, "timeout": float(timeout)},
            )
        time.sleep(float(poll_interval))
        status = blast_status(job, config=resolved, api_url=api_url)
        elapsed = time.monotonic() - started
        if progress is not None:
            progress(BlastProgress(job.rid, status.status, elapsed))
        if status.status == "ready":
            if status.has_hits is False:
                return QueryResult(
                    "sequence_similarity",
                    "NCBI BLAST",
                    job.request_url,
                    (),
                    job.provenance,
                    total_count=0,
                    metadata={"rid": job.rid, "max_identity": None, "novelty_score": 1.0},
                )
            return blast_results(job, config=resolved, api_url=api_url)
        if status.status in {"failed", "unknown"}:
            raise QueryError(
                "NCBI BLAST job failed or expired.",
                code="BLAST_JOB_FAILED",
                context={"rid": job.rid, "status": status.status},
            )


def identify(
    sequence: str | DNASequence,
    *,
    wait: bool = False,
    config: SearchConfig | None = None,
    **kwargs: object,
) -> BlastJob | QueryResult:
    """Submit an NCBI nucleotide search for likely sequence source identification.

    The default is asynchronous and returns ``BlastJob``.  Set ``wait=True``
    to poll no faster than once per minute and return normalized hits.
    """

    submit_keys = {"database", "program", "hitlist_size", "expect", "megablast", "api_url"}
    unknown = set(kwargs) - submit_keys - {"poll_interval", "timeout", "progress"}
    if unknown:
        raise TypeError(f"Unexpected identify arguments: {sorted(unknown)}")
    submit_arguments = {key: value for key, value in kwargs.items() if key in submit_keys}
    job = submit_blast(sequence, config=config, **submit_arguments)  # type: ignore[arg-type]
    if not wait:
        return job
    wait_arguments = {
        key: value
        for key, value in kwargs.items()
        if key in {"poll_interval", "timeout", "progress", "api_url"}
    }
    return wait_for_blast(job, config=config, **wait_arguments)  # type: ignore[arg-type]


def novelty(
    sequence: str | DNASequence,
    *,
    wait: bool = False,
    config: SearchConfig | None = None,
    **kwargs: object,
) -> BlastJob | QueryResult:
    """Submit or run the same BLAST search used for remote novelty inspection."""

    return identify(sequence, wait=wait, config=config, **kwargs)


__all__ = [
    "DEFAULT_NCBI_BLAST_URL",
    "BlastJob",
    "BlastProgress",
    "BlastStatus",
    "blast_results",
    "blast_status",
    "identify",
    "novelty",
    "submit_blast",
    "wait_for_blast",
]
