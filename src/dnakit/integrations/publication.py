"""Publication-oriented software evidence and citation records.

The registry distinguishes the right to cite a tool in a paper from the right
to redistribute its source or bundle its runtime.  It is intentionally
offline: every returned record is a frozen evidence snapshot that should be
paired with the exact package, executable, model and database versions used in
an experiment.
"""

from __future__ import annotations

import copy
import json
from datetime import date
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError


@lru_cache(maxsize=1)
def _records() -> list[dict[str, Any]]:
    source = Path(__file__).with_name("publication_records.json")
    payload = json.loads(source.read_text(encoding="utf-8"))
    records = payload.get("records")
    if not isinstance(records, list):
        raise ConfigurationError("Publication registry is missing its records list.")
    return records


def _record_key(value: str) -> str:
    return " ".join(value.casefold().replace("_", "-").split())


def _find(name: str) -> dict[str, Any]:
    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("name must be a non-empty string.")
    wanted = _record_key(name)
    for record in _records():
        names = [record["name"], record["id"], *record.get("aliases", [])]
        if any(_record_key(candidate) == wanted for candidate in names):
            return copy.deepcopy(record)
    raise ConfigurationError(f"No publication record is registered for {name!r}.")


def _journal_citation_status(record: dict[str, Any]) -> str:
    """Classify the source type; a DOI alone does not prove peer review."""

    explicit = record.get("publication_type")
    if explicit in {"journal_article", "preprint", "software_archive", "project_only", "source_check_required"}:
        return explicit
    doi = str(record.get("doi") or "").casefold()
    if "preprint" in doi or "peerj.preprints" in doi:
        return "preprint"
    if doi.startswith("10.5281/") or "zenodo" in doi:
        return "software_archive"
    if record.get("citation_status") == "verified_doi":
        return "journal_article"
    if record.get("citation_status") == "preprint_or_project":
        return "preprint"
    if record.get("citation_status") == "project_citation":
        return "project_only"
    return "source_check_required"


def list_publication_records(
    query: str = "",
    *,
    reuse_route: str | None = None,
    citation_status: str | None = None,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """List frozen software evidence records matching a name or route."""

    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    if reuse_route is not None and not isinstance(reuse_route, str):
        raise ConfigurationError("reuse_route must be a string or None.")
    if citation_status is not None and not isinstance(citation_status, str):
        raise ConfigurationError("citation_status must be a string or None.")
    terms = query.casefold().split()
    result = []
    for record in _records():
        haystack = json.dumps(record, ensure_ascii=False).casefold()
        if terms and not all(term in haystack for term in terms):
            continue
        if reuse_route is not None and record.get("reuse_route") != reuse_route:
            continue
        if citation_status is not None and record.get("citation_status") != citation_status:
            continue
        result.append(copy.deepcopy(record))
    return result


def describe_publication_record(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return citation, license, source and integration-boundary evidence."""

    return _find(name)


def publication_readiness(
    name: str,
    *,
    version: str | None = None,
    backend_version: str | None = None,
    database_version: str | None = None,
    accessed: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assess whether a software use is ready to record in a manuscript.

    This is an evidence checklist, not a legal opinion or a scientific
    validation. A DOI fills the citation field but does not by itself prove
    peer review; project-only records require a version and access date.
    """

    record = _find(name)
    journal_citation_status = _journal_citation_status(record)
    if version is not None and (not isinstance(version, str) or not version.strip()):
        raise ConfigurationError("version must be a non-empty string or None.")
    for field, value in (("backend_version", backend_version), ("database_version", database_version)):
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigurationError(f"{field} must be a non-empty string or None.")
    if accessed is not None and (not isinstance(accessed, str) or not accessed.strip()):
        raise ConfigurationError("accessed must be a non-empty ISO date string or None.")
    if accessed is not None:
        try:
            date.fromisoformat(accessed)
        except ValueError as exc:
            raise ConfigurationError("accessed must use YYYY-MM-DD ISO date format.") from exc
    if accessed is None and record["citation_status"] in {"project_citation", "paper_citation_requires_source_check"}:
        accessed = date.today().isoformat()
    checks = {
        "citation_evidence": record.get("citation_status") in {"verified_doi", "project_citation", "preprint_or_project"},
        "license_evidence": record.get("license_status") == "verified",
        "version_recorded": bool(version),
        "backend_version_recorded": bool(backend_version) or record.get("kind") not in {"pip_or_cli", "software_or_service", "cli"},
        "database_version_recorded": bool(database_version) or record.get("reuse_route") not in {"optional_pip_adapter", "external_cli_adapter", "optional_pip_or_external_cli"},
        "access_date_recorded": bool(accessed),
    }
    if record.get("license_status") == "restricted":
        status = "external_service_terms_review"
    elif not checks["citation_evidence"] or not checks["version_recorded"]:
        status = "incomplete_citation_evidence"
    elif record.get("license_status") != "verified":
        status = "license_review_required"
    elif record.get("reuse_route") == "reference_only_until_license_confirmed":
        status = "reference_only_pending_license"
    else:
        status = "ready_for_manuscript_record"
    return {
        "name": record["name"],
        "status": status,
        "journal_citation_status": journal_citation_status,
        "ready_for_journal": journal_citation_status == "journal_article" and checks["citation_evidence"],
        "article_citation": record["article_citation"],
        "reuse_route": record["reuse_route"],
        "checks": checks,
        "version": version,
        "backend_version": backend_version,
        "database_version": database_version,
        "accessed": accessed,
        "citation": record["citation"],
        "doi": record.get("doi"),
        "license": record.get("license"),
        "sources": record.get("source_urls", []),
        "notes": record.get("notes"),
    }


def build_software_citation(
    name: str,
    version: str,
    *,
    accessed: str | None = None,
    backend_version: str | None = None,
    database_version: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a manuscript-ready software citation record with provenance."""

    if accessed is None:
        accessed = date.today().isoformat()
    readiness = publication_readiness(
        name,
        version=version,
        backend_version=backend_version,
        database_version=database_version,
        accessed=accessed,
        verbose=False,
    )
    record = _find(name)
    citation = f"{record['name']} version {version}; {record['citation']}"
    if record.get("doi"):
        citation += f" DOI: {record['doi']}."
    else:
        citation += f" Available at {record['source_urls'][0]} (accessed {accessed})."
    return {
        "status": readiness["status"],
        "journal_citation_status": readiness["journal_citation_status"],
        "ready_for_journal": readiness["ready_for_journal"],
        "citation": citation,
        "software": record["name"],
        "version": version,
        "backend_version": backend_version,
        "database_version": database_version,
        "accessed": accessed,
        "doi": record.get("doi"),
        "license": record.get("license"),
        "reuse_route": record.get("reuse_route"),
        "sources": record.get("source_urls", []),
        "checks": readiness["checks"],
    }


def publication_registry_summary(*, verbose: bool = True) -> dict[str, Any]:
    """Return counts by citation and reuse evidence state."""

    records = _records()
    by_citation: dict[str, int] = {}
    by_reuse: dict[str, int] = {}
    for record in records:
        by_citation[record["citation_status"]] = by_citation.get(record["citation_status"], 0) + 1
        by_reuse[record["reuse_route"]] = by_reuse.get(record["reuse_route"], 0) + 1
    return {
        "records": len(records),
        "citation_status": by_citation,
        "reuse_route": by_reuse,
        "disclaimer": "Citation eligibility, source licensing and scientific validation are separate checks.",
    }


__all__ = [
    "build_software_citation",
    "describe_publication_record",
    "list_publication_records",
    "publication_readiness",
    "publication_registry_summary",
]

_register_public_verbose(globals(), tuple(__all__))
