"""biothings_client wrapper: orthologs of human or other genes through MyGene.info homologs.

The HGVS helpers (``get_hgvs_from_vcf`` and ``format_hgvs``) are already wrapped by
``gene_backends`` as ``hgvs_ids_from_vcf`` and ``format_hgvs_id``. The remaining public names
are client factories, cache and join utilities, or generic helpers that are not DNA-level.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_BIOTHINGS = "bio-biothings"
_MAX_GENES = 500
_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,63}")
_TAXA = re.compile(r"all|[0-9]{1,12}(,[0-9]{1,12}){0,31}")
_FIELDS = re.compile(r"all|[A-Za-z0-9_.,*]{1,512}")


def _version() -> str:
    try:
        return importlib.metadata.version("biothings_client")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover - source checkouts only
        return "unknown"


def _pattern(value: object, name: str, pattern: re.Pattern[str]) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value.strip()) is None:
        raise ConfigurationError(f"{name} has an unsupported value: {value!r}.")
    return value.strip()


def _plain(item: Any) -> Any:
    if item is None or isinstance(item, (str, bool, int, float)):
        return item
    if isinstance(item, dict):
        return {str(key): _plain(sub) for key, sub in sorted(item.items())}
    if isinstance(item, (list, tuple)):
        return [_plain(sub) for sub in item]
    return str(item)


def _grouped(row: Any) -> dict[str, Any]:
    """Sort homolog groups by taxon ID and genes by ID so output does not depend on dict order."""
    plain = _plain(row)
    groups = sorted(plain.get("homologs", []), key=lambda group: int(group["taxid"]))
    for group in groups:
        group["genes"] = sorted(group["genes"], key=lambda gene: str(gene.get("_id", "")))
    return {"gene": str(plain["gene"]), "homologs": groups}


def biothings_client_get_homologs(
    gene_list: Sequence[str],
    *,
    fields: str = "all",
    species: str = "all",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return homolog gene records for each gene in ``gene_list`` (MyGene.info network).

    ``species`` is ``"all"`` or comma-separated NCBI taxonomy IDs such as ``"10090,7955"``.
    Output keeps the input gene order; each entry lists homologs grouped by taxon ID.
    """

    gene_utils = require("biothings_client.utils.gene", _BIOTHINGS)
    if isinstance(gene_list, (str, bytes)) or not isinstance(gene_list, Sequence) or not gene_list:
        raise ConfigurationError("gene_list must be a non-empty list of Entrez or Ensembl IDs.")
    if len(gene_list) > _MAX_GENES:
        raise ConfigurationError(f"gene_list accepts at most {_MAX_GENES} genes per call.")
    genes = [_pattern(item, "gene_list", _ID) for item in gene_list]
    field_text = _pattern(fields, "fields", _FIELDS)
    taxa = _pattern(species, "species", _TAXA)
    try:
        rows = gene_utils.get_homologs(genes, fields=field_text, species=taxa)
    except Exception as exc:
        raise BackendExecutionError(f"biothings_client get_homologs failed: {exc}") from exc
    results = [_grouped(row) for row in rows]
    return {
        "status": "complete",
        "operation": "biothings_client_get_homologs",
        "provenance": {
            "package": "biothings_client",
            "version": _version(),
            "source": "MyGene.info homologene annotations",
            "species": taxa,
        },
        "requested": len(genes),
        "results": results,
    }


__all__ = ["biothings_client_get_homologs"]

_register_public_verbose(globals(), tuple(__all__))
