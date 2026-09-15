"""Offline discovery of curated skill metadata and bounded workflow resources."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import re
import shutil
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path, PurePosixPath
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

_STOP_WORDS = frozenset({"a", "an", "and", "for", "from", "how", "in", "of", "the", "to", "with"})


def _locations() -> tuple[Path, Path]:
    packaged = Path(__file__).resolve().parent / "_data"
    if (packaged / "catalog.json").is_file():
        return packaged, packaged / "catalog.json"
    checkout = Path(__file__).resolve().parents[3]
    return checkout / "skills", checkout / "registry" / "catalog.json"


@lru_cache(maxsize=1)
def _catalog() -> dict[str, Any]:
    _, path = _locations()
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(
            "The bundled DNAKit skill catalog is missing or unreadable.",
            code="SKILL_CATALOG_UNAVAILABLE",
            hint="Reinstall a DNAKit distribution containing the skill catalog.",
        ) from exc
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise ConfigurationError("Unsupported skill catalog schema.", code="SKILL_CATALOG_INVALID")
    return value


def _limit(value: int, name: str, minimum: int, maximum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"{name} must be an integer in [{minimum}, {maximum}].", code="SKILL_ARGUMENT_INVALID"
        )


def _records(category: str | None, source: str | None) -> list[dict[str, Any]]:
    catalog = _catalog()
    if category is not None and category not in catalog["categories"]:
        raise ConfigurationError(
            f"Unknown skill category: {category}", code="SKILL_CATEGORY_INVALID"
        )
    if source is not None and source not in {entry["id"] for entry in catalog["sources"]}:
        raise ConfigurationError(f"Unknown skill source: {source}", code="SKILL_SOURCE_INVALID")
    return [
        entry
        for entry in catalog["skills"]
        if (category is None or entry["category"] == category)
        and (source is None or entry["source_id"] == source)
    ]


def _summary(entry: dict[str, Any]) -> dict[str, Any]:
    fields = (
        "id",
        "description",
        "summary_zh",
        "category",
        "capability",
        "source_id",
        "license",
        "validation",
    )
    return copy.deepcopy({key: entry[key] for key in fields if key in entry})


def list_skill_categories() -> dict[str, Any]:
    """Return category labels and counts without reading workflow bodies."""
    counts = Counter(entry["category"] for entry in _catalog()["skills"])
    return {
        key: {**copy.deepcopy(value), "count": counts[key]}
        for key, value in _catalog()["categories"].items()
    }


def skill_catalog_stats() -> dict[str, Any]:
    """Summarize the bundled release, sources, licenses, and validation states."""
    catalog = _catalog()
    entries = catalog["skills"]
    return {
        "schema_version": catalog["schema_version"],
        "snapshot_date": catalog["snapshot_date"],
        "count": len(entries),
        "sources": dict(Counter(entry["source_id"] for entry in entries)),
        "licenses": dict(Counter(entry["license"] for entry in entries)),
        "execution_validation": dict(
            Counter(entry["validation"]["execution"] for entry in entries)
        ),
        "selection": copy.deepcopy(catalog.get("selection", {})),
    }


def list_skills(
    *, category: str | None = None, source: str | None = None, offset: int = 0, limit: int = 20
) -> dict[str, Any]:
    """Return a bounded page of skill summaries."""
    _limit(offset, "offset", 0, 1_000_000)
    _limit(limit, "limit", 1, 100)
    entries = _records(category, source)
    return {
        "total": len(entries),
        "offset": offset,
        "limit": limit,
        "skills": [_summary(entry) for entry in entries[offset : offset + limit]],
    }


def search_skills(
    query: str, *, category: str | None = None, source: str | None = None, limit: int = 10
) -> dict[str, Any]:
    """Rank skill summaries locally, with curated Chinese-to-English query aliases.

    No embedding API, network access, workflow execution, or installation is performed.
    """
    _limit(limit, "limit", 1, 100)
    if not isinstance(query, str) or not query.strip() or len(query) > 2000:
        raise ConfigurationError(
            "query must contain 1-2000 characters.", code="SKILL_QUERY_INVALID"
        )
    normalized = query.casefold().strip()
    for phrase, replacement in sorted(
        _catalog().get("query_aliases", {}).items(), key=lambda pair: -len(pair[0])
    ):
        normalized = normalized.replace(phrase.casefold(), f" {replacement} ")
    terms = set(re.findall(r"[a-z0-9]+|[\u4e00-\u9fff]+", normalized)) - _STOP_WORDS
    phrase = " ".join(re.findall(r"[a-z0-9]+|[\u4e00-\u9fff]+", normalized))
    if not terms:
        return {"query": query, "total": 0, "skills": []}
    scored: list[tuple[float, dict[str, Any]]] = []
    for entry in _records(category, source):
        identity = " ".join(
            [entry["id"], entry["capability"], *entry.get("aliases", [])]
        ).casefold()
        labels = " ".join([entry.get("summary_zh", ""), *entry.get("keywords", [])]).casefold()
        description = entry["description"].casefold()
        fields = ((identity, 8.0), (labels, 5.0), (description, 2.0))
        score = 0.0
        matched = 0
        for term in terms:
            weights = [weight for value, weight in fields if term in value]
            if weights:
                matched += 1
                score += max(weights)
        if not matched or matched / len(terms) < 0.5:
            continue
        score *= matched / len(terms)
        if matched == len(terms):
            score += 10.0
        if len(terms) > 1 and phrase in identity.replace("-", " "):
            score += 20.0
        if query.casefold() == entry["id"] or query.casefold() in entry.get("aliases", []):
            score += 100.0
        scored.append((score, entry))
    scored.sort(key=lambda pair: (-pair[0], pair[1]["id"]))
    return {
        "query": query,
        "total": len(scored),
        "skills": [
            {**_summary(entry), "score": round(score, 3)} for score, entry in scored[:limit]
        ],
    }


def get_skill(skill_id: str) -> dict[str, Any]:
    """Return provenance, dependencies, resource manifest, and validation metadata."""
    if not isinstance(skill_id, str) or not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", skill_id):
        raise ConfigurationError("Invalid skill identifier.", code="SKILL_ID_INVALID")
    entries = _catalog()["skills"]
    for entry in entries:
        if entry["id"] == skill_id:
            return copy.deepcopy(entry)
    aliases = [entry for entry in entries if skill_id in entry.get("aliases", [])]
    if len(aliases) == 1:
        return copy.deepcopy(aliases[0])
    raise ConfigurationError(
        f"Unknown or ambiguous skill: {skill_id}",
        code="SKILL_NOT_FOUND",
        hint="Use dnakit skills search to find an installed catalog identifier.",
    )


def _resource_bytes(entry: dict[str, Any], resource: str) -> bytes:
    relative = PurePosixPath(resource)
    if (
        not resource
        or relative.is_absolute()
        or ".." in relative.parts
        or "\\" in resource
        or str(relative) != resource
    ):
        raise ConfigurationError(
            "Invalid relative skill resource path.", code="SKILL_RESOURCE_INVALID"
        )
    record = next((item for item in entry["files"] if item["path"] == resource), None)
    if record is None:
        raise ConfigurationError(
            f"Resource is not bundled: {resource}", code="SKILL_RESOURCE_NOT_FOUND"
        )
    root, _ = _locations()
    directory: Path = root / entry["id"]
    path: Path = directory / resource
    if path.is_symlink() or not path.resolve().is_relative_to(directory.resolve()):
        raise ConfigurationError(
            "Skill resources cannot escape their directory.", code="SKILL_RESOURCE_INVALID"
        )
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise InputFormatError(
            f"Bundled resource checksum mismatch: {entry['id']}/{resource}",
            code="SKILL_RESOURCE_CHANGED",
            hint="Rebuild the skill catalog after intentional source changes, or reinstall DNAKit.",
        )
    return data


def read_skill_file(
    skill_id: str, *, resource: str = "SKILL.md", offset: int = 0, limit: int = 24000
) -> dict[str, Any]:
    """Read a checksum-verified UTF-8 resource with character-based pagination."""
    _limit(offset, "offset", 0, 100_000_000)
    _limit(limit, "limit", 1, 100_000)
    entry = get_skill(skill_id)
    try:
        content = _resource_bytes(entry, resource).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "This resource is binary; install the skill to access the original file.",
            code="SKILL_RESOURCE_BINARY",
        ) from exc
    end = min(offset + limit, len(content))
    return {
        "skill": _summary(entry),
        "resource": resource,
        "content": content[offset:end],
        "offset": offset,
        "total_characters": len(content),
        "next_offset": end if end < len(content) else None,
        "upstream_url": entry["upstream_url"],
        "runtime": copy.deepcopy(entry.get("runtime", {})),
        "required_skills": copy.deepcopy(entry.get("required_skills", [])),
        "attribution": entry["attribution"],
        "upstream_revision": entry["upstream_revision"],
        "files": [item["path"] for item in entry["files"]] if resource == "SKILL.md" else [],
    }


def check_skill(skill_id: str) -> dict[str, Any]:
    """Inspect declared local dependencies without importing or running third-party code.

    API access, reference data, model licenses, and scientific output correctness require
    separate verification. A successful local probe is not an execution validation.
    """
    entry = get_skill(skill_id)
    runtime = entry.get("runtime", {})
    checks: list[dict[str, Any]] = []
    if minimum := runtime.get("python_minimum"):
        parts = tuple(int(part) for part in minimum.split("."))
        checks.append(
            {
                "kind": "python_version",
                "name": ">=" + minimum,
                "available": sys.version_info[: len(parts)] >= parts,
                "current": ".".join(map(str, sys.version_info[:3])),
            }
        )
    for module in runtime.get("python_modules", []):
        try:
            present = importlib.util.find_spec(module.split(".")[0]) is not None
        except (ValueError, ImportError):
            present = False
        checks.append({"kind": "python_module", "name": module, "available": present})
    for command in runtime.get("commands", []):
        checks.append(
            {"kind": "command", "name": command, "available": shutil.which(command) is not None}
        )
    for variable in runtime.get("environment", []):
        checks.append(
            {"kind": "environment", "name": variable, "available": bool(os.environ.get(variable))}
        )
    missing = any(not check["available"] for check in checks)
    return {
        "id": entry["id"],
        "status": "missing_dependencies" if missing else "local_checks_passed",
        "checks": checks,
        "manual_checks": copy.deepcopy(runtime.get("manual_checks", [])),
        "execution_validation": entry["validation"]["execution"],
        "note": (
            "Local availability checks do not verify API access, reference "
            "data, or workflow execution."
        ),
    }


def cite_skill(skill_id: str, *, format: str = "text") -> str | dict[str, Any]:
    """Return a version-pinned software citation and the upstream citation guidance."""
    entry = get_skill(skill_id)
    citation = {
        "skill": entry["id"],
        "upstream_name": entry["upstream_name"],
        "source": entry["source_id"],
        "author": entry["attribution"],
        "version": entry["upstream_revision"],
        "url": entry["upstream_url"],
        "license": entry["license"],
        "citation_guidance": entry.get("citation_guidance", ""),
        "note": "Also cite the tools, models, databases, and papers actually used in the analysis.",
    }
    if format == "json":
        return citation
    if format == "text":
        return (
            f"{citation['author']}. {citation['upstream_name']}. "
            f"Revision {citation['version']}. {citation['url']}\n"
            f"License: {citation['license']}. Distributed in DNAKit as {entry['id']}.\n"
            f"{citation['citation_guidance']}\n{citation['note']}"
        )
    if format == "bibtex":

        def escape(value: object) -> str:
            return (
                str(value)
                .replace("\\", "\\textbackslash{}")
                .replace("{", "\\{")
                .replace("}", "\\}")
                .replace("%", "\\%")
                .replace("\n", " ")
            )

        key = str(entry["id"]).replace("-", "_")
        fields = {
            "author": "{" + escape(citation["author"]) + "}",
            "title": escape(citation["upstream_name"]),
            "version": escape(citation["version"]),
            "url": escape(citation["url"]),
            "note": "Distributed in DNAKit; license " + escape(citation["license"]),
        }
        return (
            "@software{"
            + key
            + ",\n"
            + ",\n".join("  " + name + " = {" + value + "}" for name, value in fields.items())
            + "\n}"
        )
    raise ConfigurationError(
        "format must be text, json, or bibtex.", code="SKILL_CITATION_FORMAT_INVALID"
    )
