"""mygene wrapper: batch gene-record lookup by Entrez or Ensembl ID (MyGene.info network).

``mygene.MyGeneInfo`` methods ``query``, ``querymany`` and ``getgene`` are already wrapped by
``gene_backends`` (``query_mygene_genes``, ``map_gene_identifiers``, ``fetch_mygene_gene``).
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_MYGENE = "bio-mygene"
_MAX_IDS = 1000
_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,63}")
_TAXA = re.compile(r"[A-Za-z0-9_,]{1,128}")
_FIELDS = re.compile(r"[A-Za-z0-9_.,*]{1,512}")
_DEFAULT_FIELDS = "symbol,name,entrezgene,ensembl.gene,taxid"


def _version() -> str:
    try:
        return importlib.metadata.version("mygene")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - source checkouts only
        return "unknown"


def _provenance(**notes: object) -> dict[str, Any]:
    return {"package": "mygene", "version": _version(), **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "complete", "operation": operation, "provenance": provenance, **fields}


def _pattern(value: object, name: str, pattern: re.Pattern[str], *, limit: int) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > limit:
        raise ConfigurationError(f"{name} must be a non-empty string of at most {limit} chars.")
    text = value.strip()
    if pattern.fullmatch(text) is None:
        raise ConfigurationError(f"{name} contains unsupported characters: {text!r}.")
    return text


def _ids(values: object) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("gene_ids must be a non-empty list of Entrez or Ensembl IDs.")
    if len(values) > _MAX_IDS:
        raise ConfigurationError(f"gene_ids accepts at most {_MAX_IDS} entries per call.")
    return [_pattern(item, "gene_ids", _ID, limit=64) for item in values]


def _fields(value: object) -> str:
    if isinstance(value, str):
        return _pattern(value, "fields", _FIELDS, limit=512)
    if isinstance(value, Sequence) and value and all(isinstance(item, str) for item in value):
        return _pattern(",".join(value), "fields", _FIELDS, limit=512)
    raise ConfigurationError("fields must be a comma-separated string or a list of strings.")


def mygene_getgenes(
    gene_ids: Sequence[str],
    *,
    fields: str | Sequence[str] = _DEFAULT_FIELDS,
    species: str = "human",
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch MyGene.info gene records for Entrez or Ensembl IDs in one batch (network).

    Records are returned in the order of ``gene_ids``; IDs MyGene does not know are reported
    in ``not_found`` instead of being dropped silently.
    """

    mygene = require("mygene", _MYGENE)
    ids = _ids(gene_ids)
    field_text = _fields(fields)
    taxon = _pattern(species, "species", _TAXA, limit=128)
    client = mygene.MyGeneInfo()
    try:
        rows = client.getgenes(ids, fields=field_text, species=taxon, verbose=verbose)
    except Exception as exc:
        raise BackendExecutionError(f"mygene getgenes failed: {exc}") from exc

    position = {identifier: index for index, identifier in enumerate(ids)}
    found: list[dict[str, Any]] = []
    missing: list[str] = []
    for row in rows or []:
        query = str(row.get("query", row.get("_id", "")))
        if row.get("notfound"):
            missing.append(query)
        else:
            found.append({key: _value(item) for key, item in sorted(row.items())})
    found.sort(key=lambda row: position.get(str(row.get("query", row.get("_id"))), len(ids)))
    return _document(
        "mygene_getgenes",
        _provenance(source="MyGene.info gene batch endpoint", species=taxon),
        requested=len(ids),
        found=len(found),
        not_found=sorted(set(missing), key=lambda item: position.get(item, len(ids))),
        records=found,
    )


def _value(item: Any) -> Any:
    if item is None or isinstance(item, (str, bool, int, float)):
        return item
    if isinstance(item, dict):
        return {str(key): _value(sub) for key, sub in sorted(item.items())}
    if isinstance(item, (list, tuple)):
        return [_value(sub) for sub in item]
    return str(item)


__all__ = ["mygene_getgenes"]

_register_public_verbose(globals(), tuple(__all__))
