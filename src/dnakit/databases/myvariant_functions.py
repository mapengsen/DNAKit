"""myvariant wrapper: batch variant lookup by HGVS ID or rsID (MyVariant.info network).

``MyVariantInfo`` methods ``getvariants``, ``query``, ``format_hgvs`` and ``get_hgvs_from_vcf``
are already wrapped by ``gene_backends`` (``fetch_variant_annotations``,
``query_variant_annotations``, ``format_hgvs_id`` and ``hgvs_ids_from_vcf``).
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_MYVARIANT = "bio-myvariant"
_MAX_TERMS = 1000
_TERM = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:>()+*=-]{0,299}")
_SCOPES = re.compile(r"[A-Za-z0-9_.,]{1,128}")
_FIELDS = re.compile(r"[A-Za-z0-9_.,*]{1,512}")
_DEFAULT_FIELDS = "dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,dbnsfp.genename"


def _version() -> str:
    try:
        return importlib.metadata.version("myvariant")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - source checkouts only
        return "unknown"


def _pattern(value: object, name: str, pattern: re.Pattern[str], *, limit: int) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > limit:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {limit} chars.")
    text = value.strip()
    if pattern.fullmatch(text) is None:
        raise ConfigurationError(f"{name} contains unsupported characters: {text!r}.")
    return text


def _terms(values: object) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("query_terms must be a non-empty list of HGVS IDs or rsIDs.")
    if len(values) > _MAX_TERMS:
        raise ConfigurationError(f"query_terms accepts at most {_MAX_TERMS} entries per call.")
    return [_pattern(item, "query_terms", _TERM, limit=300) for item in values]


def _flag(value: object) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError("verbose must be a boolean.")
    return value


def _plain(item: Any) -> Any:
    if item is None or isinstance(item, (str, bool, int, float)):
        return item
    if isinstance(item, dict):
        return {str(key): _plain(sub) for key, sub in sorted(item.items())}
    if isinstance(item, (list, tuple)):
        return [_plain(sub) for sub in item]
    return str(item)


def myvariant_querymany(
    query_terms: Sequence[str],
    *,
    scopes: str = "dbsnp.rsid",
    fields: str = _DEFAULT_FIELDS,
    assembly: str = "hg38",
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up many variants by HGVS ID or rsID in one MyVariant.info batch (network).

    Results follow the order of ``query_terms``; terms MyVariant does not resolve are listed
    in ``not_found``.
    """

    myvariant = require("myvariant", _MYVARIANT)
    terms = _terms(query_terms)
    scope = _pattern(scopes, "scopes", _SCOPES, limit=128)
    field_text = _pattern(fields, "fields", _FIELDS, limit=512)
    build = _pattern(assembly, "assembly", re.compile(r"hg19|hg38"), limit=8)
    flag = _flag(verbose)
    client = myvariant.MyVariantInfo()
    try:
        rows = client.querymany(
            terms, scopes=scope, fields=field_text, assembly=build, verbose=flag
        )
    except Exception as exc:
        raise BackendExecutionError(f"myvariant querymany failed: {exc}") from exc

    position = {term: index for index, term in enumerate(terms)}
    hits: list[dict[str, Any]] = []
    missing: list[str] = []
    for row in rows or []:
        query = str(row.get("query", ""))
        if row.get("notfound"):
            missing.append(query)
        else:
            hits.append(_plain(row))
    hits.sort(key=lambda row: position.get(str(row.get("query", "")), len(terms)))
    matched = {str(row.get("query", "")) for row in hits}
    return {
        "status": "complete",
        "operation": "myvariant_querymany",
        "provenance": {"package": "myvariant", "version": _version(), "assembly": build},
        "scopes": scope,
        "requested": len(terms),
        "matched_terms": len(matched),
        "not_found": sorted(
            set(missing) - matched, key=lambda item: position.get(item, len(terms))
        ),
        "hit_count": len(hits),
        "records": hits,
    }


__all__ = ["myvariant_querymany"]

_register_public_verbose(globals(), tuple(__all__))
