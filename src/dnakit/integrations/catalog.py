"""Offline discovery of reviewed analyses and their original scientific sources."""

from __future__ import annotations

import copy
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, cast

from dnakit.exceptions import ConfigurationError


@lru_cache(maxsize=2)
def _data(name: str) -> dict[str, Any]:
    return cast(
        dict[str, Any], json.loads(Path(__file__).with_name(name).read_text(encoding="utf-8"))
    )


def describe_analysis(analysis: str, *, verbose: bool = True) -> dict[str, Any]:
    """Describe inputs, parameters, dependencies, outputs and sources without execution."""
    for entry in _data("analyses.json")["analyses"]:
        if entry["id"] == analysis:
            return copy.deepcopy(entry)
    raise ConfigurationError(f"Unknown reviewed analysis: {analysis!r}.")


def list_analyses(
    query: str = "", *, domain: str | None = None, verbose: bool = True
) -> list[dict[str, Any]]:
    """Search available executable adapters; availability does not imply installation."""
    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    terms = query.casefold().split()
    result = []
    for entry in _data("analyses.json")["analyses"]:
        text = json.dumps(entry, ensure_ascii=False).casefold()
        if (domain is None or domain in entry["domains"]) and all(t in text for t in terms):
            result.append(copy.deepcopy(entry))
    return result


def search_projects(
    query: str = "", *, domain: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Search the dated evidence survey, including projects not integrated into DNAKit.

    Project entries distinguish skills, MCP servers, software, models and papers.
    A search hit is evidence of discovery, not a claim of local execution.
    """
    if not isinstance(query, str):
        raise ConfigurationError("query must be a string.")
    data = _data("projects.json")
    terms = query.casefold().split()
    projects = [
        copy.deepcopy(entry)
        for entry in data["projects"]
        if (domain is None or domain in entry["domains"])
        and all(t in json.dumps(entry, ensure_ascii=False).casefold() for t in terms)
    ]
    return {"snapshot_date": data["snapshot_date"], "total": len(projects), "projects": projects}
