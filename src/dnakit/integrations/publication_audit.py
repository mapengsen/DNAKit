"""Reproducible publication-evidence audit for the 320 surveyed tools.

The audit intentionally keeps three questions independent:

* ``journal_citation_status`` says whether a peer-reviewed journal article was
  identified. A DOI by itself is not enough: preprints and software archives
  are labelled separately.
* ``manuscript_category`` says what a Methods section can do next.
* ``license_status`` and ``scientific_validation`` are independent fields;
  citation evidence never grants redistribution rights or validates results.

The source inventory and the curated DOI evidence are frozen JSON snapshots so
the report can be regenerated offline and reviewed in a manuscript audit.
"""

from __future__ import annotations

import copy
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError

_CATEGORIES = (
    "manuscript_ready",
    "needs_version_or_access_date",
    "project_only",
    "missing_reliable_citation",
    "license_restricted",
)
_PUBLICATION_TYPES = {
    "journal_article",
    "preprint",
    "software_archive",
    "project_only",
    "source_check_required",
    "missing",
}


def _key(value: str) -> str:
    """Normalize a display name without throwing away non-Latin characters."""

    return " ".join(value.casefold().replace("_", "-").split())


def _data(name: str) -> dict[str, Any]:
    path = Path(__file__).with_name(name)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ConfigurationError(f"{name} must contain a JSON object.")
    return payload


@lru_cache(maxsize=1)
def _inventory() -> list[dict[str, Any]]:
    records = _data("third_party_tools.json").get("records")
    if not isinstance(records, list):
        raise ConfigurationError("third_party_tools.json is missing records.")
    return records


@lru_cache(maxsize=1)
def _evidence_by_name() -> dict[str, dict[str, Any]]:
    """Index verified local records and manually reviewed DOI records."""

    evidence: dict[str, dict[str, Any]] = {}
    for filename in ("publication_records.json", "publication_audit_records.json"):
        payload = _data(filename)
        for record in payload.get("records", []):
            if not isinstance(record, dict):
                continue
            aliases = [record.get("name"), record.get("id"), *record.get("aliases", [])]
            for alias in aliases:
                if isinstance(alias, str) and alias.strip():
                    # The audit overrides the small publication registry when
                    # it has an explicit publication type for the same name.
                    evidence[_key(alias)] = copy.deepcopy(record)
    return evidence


def _lookup_evidence(name: str) -> dict[str, Any] | None:
    """Resolve exact names and survey labels with parenthetical suffixes."""

    wanted = _key(name)
    evidence = _evidence_by_name()
    exact = evidence.get(wanted)
    if exact is not None:
        return exact
    # The survey uses labels such as ``GATK4（现行源码）`` and
    # ``MMseqs2（当前主仓库）``. Match a reviewed alias only when it is the
    # complete leading token, avoiding accidental BWA/BWA-MEM2 collisions.
    for alias, record in evidence.items():
        if wanted.startswith(alias + "（") or wanted.startswith(alias + " ("):
            return record
    return None


def _publication_type(record: dict[str, Any]) -> tuple[str, dict[str, Any] | None]:
    matched = _lookup_evidence(record["name"])
    if matched:
        explicit = matched.get("publication_type")
        if explicit in _PUBLICATION_TYPES:
            return explicit, matched
        status = matched.get("citation_status")
        if status == "verified_doi":
            # Existing DOI records are journal articles unless their DOI is a
            # preprint or software archive. Those two exceptions are explicit
            # because a DOI alone cannot establish peer review.
            doi = str(matched.get("doi") or "").casefold()
            if "preprint" in doi or "peerj.preprints" in doi:
                return "preprint", matched
            if doi.startswith("10.5281/") or "zenodo" in doi:
                return "software_archive", matched
            return "journal_article", matched
        if status == "preprint_or_project":
            return "preprint", matched
        if status == "project_citation":
            return "project_only", matched
        if status == "paper_citation_requires_source_check":
            return "source_check_required", matched
    # A survey citation reference without a DOI is not promoted to a journal
    # article: it needs a source check before manuscript use.
    if record.get("citation_references"):
        return "source_check_required", None
    if record.get("urls") or "prose_link" in record.get("evidence", []):
        return "project_only", None
    return "missing", None


def _license_status(record: dict[str, Any], matched: dict[str, Any] | None) -> tuple[str, str | None]:
    values = [str(value) for value in record.get("licenses", []) if value]
    if matched and matched.get("license"):
        values.append(str(matched["license"]))
    text = " ".join(values).casefold()
    decision = str(record.get("primary_decision", "")).casefold()
    restrictive_terms = (
        "agpl",
        "gpl",
        "cc by-nc",
        "商业",
        "commercial",
        "academic-purpose",
        "academic purposes",
        "restricted",
        "未找到",
        "未确认",
        "不清晰",
        "缺少",
    )
    if decision in {"external_only", "blocked_or_unverified"} or any(term in text for term in restrictive_terms):
        return "restricted_or_unverified", "Redistribution or linking terms require a separate license review."
    if values:
        return "declared", "; ".join(values)
    return "unknown", "The survey does not provide a license evidence record."


def _category(publication_type: str, license_status: str) -> str:
    if license_status == "restricted_or_unverified":
        return "license_restricted"
    if publication_type == "journal_article":
        return "manuscript_ready"
    if publication_type in {"preprint", "software_archive"}:
        return "needs_version_or_access_date"
    if publication_type == "project_only":
        return "project_only"
    return "missing_reliable_citation"


def _requirements(publication_type: str) -> list[str]:
    if publication_type == "journal_article":
        return ["exact_software_version"]
    if publication_type in {"preprint", "software_archive"}:
        return ["exact_software_version", "access_date", "label_non_journal_source"]
    if publication_type == "project_only":
        return ["exact_software_version", "access_date", "official_project_url"]
    if publication_type == "source_check_required":
        return ["verify_primary_paper_or_doi", "exact_software_version", "access_date"]
    return ["identify_primary_paper_or_doi", "exact_software_version"]


def _audit_one(record: dict[str, Any]) -> dict[str, Any]:
    publication_type, matched = _publication_type(record)
    license_status, license_note = _license_status(record, matched)
    doi = (matched or {}).get("doi")
    citation = (matched or {}).get("citation")
    source_urls = list((matched or {}).get("source_urls", []))
    if not source_urls:
        source_urls = list(record.get("urls", []))
    category = _category(publication_type, license_status)
    return {
        "name": record["name"],
        "key": record.get("key"),
        "source_lines": record.get("source_lines", []),
        "primary_decision": record.get("primary_decision"),
        "integration_status": record.get("integration_status"),
        "local_match": bool(record.get("local_match")),
        "journal_citation_status": publication_type,
        "citation_evidence": {
            "doi": doi,
            "citation": citation,
            "source_urls": source_urls,
            "verified_registry": matched is not None,
        },
        "manuscript_category": category,
        "required_metadata": _requirements(publication_type),
        "license_status": license_status,
        "license_note": license_note,
        "reuse_route": (matched or {}).get("reuse_route"),
        "scientific_validation": "not_assessed_by_publication_audit",
        "notes": (matched or {}).get("notes"),
    }


def audit_publication_evidence(
    query: str = "",
    *,
    manuscript_category: str | None = None,
    journal_citation_status: str | None = None,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """Audit every surveyed project, or filter the frozen audit by name/state."""

    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    if manuscript_category is not None and manuscript_category not in _CATEGORIES:
        raise ConfigurationError(f"Unknown manuscript_category: {manuscript_category!r}.")
    if journal_citation_status is not None and journal_citation_status not in _PUBLICATION_TYPES:
        raise ConfigurationError(f"Unknown journal_citation_status: {journal_citation_status!r}.")
    terms = query.casefold().split()
    result = []
    for source in _inventory():
        name = str(source.get("name", ""))
        if terms and not all(term in name.casefold() for term in terms):
            continue
        item = _audit_one(source)
        if manuscript_category is not None and item["manuscript_category"] != manuscript_category:
            continue
        if journal_citation_status is not None and item["journal_citation_status"] != journal_citation_status:
            continue
        result.append(item)
    return result


def describe_publication_audit(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the single audit row for an exact surveyed project name."""

    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("name must be a non-empty string.")
    rows = audit_publication_evidence(name, verbose=False)
    exact = [row for row in rows if _key(row["name"]) == _key(name)]
    if len(exact) != 1:
        raise ConfigurationError(f"No unique surveyed project named {name!r}.")
    return exact[0]


def publication_audit_summary(*, verbose: bool = True) -> dict[str, Any]:
    """Return five manuscript categories plus independent evidence counts."""

    rows = audit_publication_evidence(verbose=False)
    categories = {name: 0 for name in _CATEGORIES}
    citation_status: dict[str, int] = {name: 0 for name in _PUBLICATION_TYPES}
    licenses: dict[str, int] = {}
    for row in rows:
        categories[row["manuscript_category"]] += 1
        citation_status[row["journal_citation_status"]] += 1
        state = row["license_status"]
        licenses[state] = licenses.get(state, 0) + 1
    return {
        "records": len(rows),
        "survey_records_expected": 320,
        "survey_complete": len(rows) == 320,
        "manuscript_categories": categories,
        "journal_citation_status": citation_status,
        "license_status": licenses,
        "ready_for_journal": sum(row["journal_citation_status"] == "journal_article" and row["manuscript_category"] == "manuscript_ready" for row in rows),
        "scientific_validation": "not_assessed_by_publication_audit",
        "disclaimer": "A journal citation does not grant redistribution rights and does not validate a scientific result; report exact versions, databases, parameters and independent validation separately.",
    }


__all__ = [
    "audit_publication_evidence",
    "describe_publication_audit",
    "publication_audit_summary",
]

_register_public_verbose(globals(), tuple(__all__))
