"""NCBI E-utilities DNA records, database searches, summaries and typed links.

Original upstream defaults are retained. One request is made per invocation;
callers must obey NCBI's 3 requests/second (10 with API key) across processes.
"""

from __future__ import annotations

import threading
import time
from collections.abc import Mapping, Sequence
from typing import Any
from xml.etree import ElementTree

from dnakit.databases._network_p44 import _request
from dnakit.exceptions import ConfigurationError, QueryError
from dnakit.search._http import build_url, require_text

_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_FETCH = {
    "nuccore": {"fasta", "gbc", "ft", "native", "fasta_cds_na"},
    "gene": {"gene_table", "xml", "asn.1"},
    "snp": {"xml", "docsum"},
    "clinvar": {"vcv"},
    "sra": {"full", "xml"},
}
_SEARCH = {"assembly", "annotinfo", "snp", "dbvar", "grasp", "geoprofiles", "gtr", "seqannot"}
_SUMMARY = {"snp", "sra", "omim"}
_LINKS = {
    "nuccore_sra_wgs": ("nuccore", "sra"),
    "genome_assembly": ("genome", "assembly"),
    "assembly_assembly_diploid": ("assembly", "assembly"),
    "assembly_genome": ("assembly", "genome"),
    "assembly_sra": ("assembly", "sra"),
    "gene_clinvar": ("gene", "clinvar"),
    "gene_clinvar_specific": ("gene", "clinvar"),
    "gene_gene_neighbors": ("gene", "gene"),
    "gene_omim": ("gene", "omim"),
    "snp_clinvar": ("snp", "clinvar"),
    "snp_gene": ("snp", "gene"),
    "clinvar_dbvar": ("clinvar", "dbvar"),
    "clinvar_gene": ("clinvar", "gene"),
    "clinvar_gene_specific": ("clinvar", "gene"),
    "clinvar_gtr": ("clinvar", "gtr"),
    "clinvar_omim": ("clinvar", "omim"),
    "clinvar_snp": ("clinvar", "snp"),
    "dbvar_clinvar": ("dbvar", "clinvar"),
    "sra_assembly": ("sra", "assembly"),
    "sra_gap_all": ("sra", "gap"),
    "sra_gds": ("sra", "gds"),
    "sra_nuccore_wgs": ("sra", "nuccore"),
    "gap_sra_all": ("gap", "sra"),
    "gds_gds_similar": ("gds", "gds"),
    "gds_sra": ("gds", "sra"),
    "omim_clinvar": ("omim", "clinvar"),
    "omim_gene": ("omim", "gene"),
    "omim_omim": ("omim", "omim"),
}
_LOCK = threading.Lock()
_LAST = 0.0


def _call(
    endpoint: str,
    params: Mapping[str, Any],
    email: str | None,
    api_key: str | None,
    timeout: float,
    max_bytes: int,
    json_response: bool,
) -> dict[str, Any]:
    global _LAST
    query = dict(params)
    query["tool"] = "DNAKit"
    if email is not None:
        query["email"] = require_text(email, "email")
    if api_key is not None:
        query["api_key"] = require_text(api_key, "api_key")
    # This rate limit covers all four adapters in this process.
    with _LOCK:
        delay = (0.11 if api_key else 0.35) - (time.monotonic() - _LAST)
        if delay > 0:
            time.sleep(delay)
        _LAST = time.monotonic()
    result = _request(
        build_url(_BASE, "/" + endpoint + ".fcgi", query.items()),
        timeout=timeout,
        max_bytes=max_bytes,
        json_response=json_response,
    )
    if not json_response and str(result["data"]).lstrip().startswith("<?xml"):
        try:
            root = ElementTree.fromstring(result["data"])
        except ElementTree.ParseError as exc:
            raise QueryError("NCBI returned malformed XML.") from exc
        errors = root.findall(".//ERROR")
        if errors:
            raise QueryError("NCBI E-utilities returned an error: " + str(errors[0].text))
    result["provenance"]["service"] = "NCBI E-utilities"
    return result


def _ids(ids: Sequence[str]) -> str:
    if isinstance(ids, (str, bytes)) or not isinstance(ids, Sequence) or not 1 <= len(ids) <= 200:
        raise ConfigurationError("ids must contain 1-200 identifiers.")
    values = [require_text(x, "identifier", max_length=100) for x in ids]
    if any("," in x for x in values):
        raise ConfigurationError("Each identifier must be a single ID.")
    return ",".join(values)


def ncbi_fetch_dna_records(
    database: str,
    ids: Sequence[str],
    *,
    rettype: str = "fasta",
    retmode: str = "text",
    email: str | None = None,
    api_key: str | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch nucleotide FASTA/gbc/ft/native/CDS, Gene XML/table/ASN.1, SNP, ClinVar or SRA records.

    Explicit ``rettype`` and ``retmode`` reproduce upstream EFetch options.
    CDS peptide output is deliberately excluded from this DNA-only adapter.
    """
    if database not in _FETCH or rettype not in _FETCH[database] or retmode not in {"text", "xml"}:
        raise ConfigurationError("Unsupported DNA EFetch database/rettype/retmode combination.")
    return _call(
        "efetch",
        {"db": database, "id": _ids(ids), "rettype": rettype, "retmode": retmode},
        email,
        api_key,
        timeout,
        max_bytes,
        False,
    )


def ncbi_search_dna_database(
    database: str,
    term: str,
    *,
    retmax: int = 20,
    retstart: int = 0,
    email: str | None = None,
    api_key: str | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """ESearch assembly, annotinfo, snp, dbvar, grasp, geoprofiles, gtr or seqannot; return IDs
    and warnings.
    """
    if database not in _SEARCH:
        raise ConfigurationError("Unsupported DNA ESearch database.")
    if (
        any(isinstance(v, bool) or not isinstance(v, int) for v in (retmax, retstart))
        or not 0 <= retmax <= 10000
        or retstart < 0
    ):
        raise ConfigurationError("retmax must be 0-10000 and retstart non-negative.")
    result = _call(
        "esearch",
        {
            "db": database,
            "term": require_text(term, "term"),
            "retmode": "json",
            "retmax": retmax,
            "retstart": retstart,
        },
        email,
        api_key,
        timeout,
        max_bytes,
        True,
    )
    payload = result["data"].get("esearchresult", {})
    if payload.get("error"):
        raise QueryError(str(payload["error"]))
    # Preserve warninglist/errorlist rather than silently dropping unmatched terms.
    return result


def ncbi_summarize_dna_records(
    database: str,
    ids: Sequence[str],
    *,
    email: str | None = None,
    api_key: str | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return unmodified ESummary JSON records from snp, sra or omim."""
    if database not in _SUMMARY:
        raise ConfigurationError("database must be snp, sra or omim.")
    return _call(
        "esummary",
        {"db": database, "id": _ids(ids), "retmode": "json"},
        email,
        api_key,
        timeout,
        max_bytes,
        True,
    )


def ncbi_link_dna_records(
    linkname: str,
    ids: Sequence[str],
    *,
    email: str | None = None,
    api_key: str | None = None,
    timeout: float = 60,
    max_bytes: int = 20_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """ELink one of the 28 curated nuccore/assembly/gene/variant/SRA/GEO/OMIM relationships.

    ``linkname`` selects an exact upstream link name from ``_LINKS``; results
    preserve NCBI's link sets, ID ordering and empty-link responses.
    """
    if linkname not in _LINKS:
        raise ConfigurationError(f"Unsupported linkname; choose from {sorted(_LINKS)}.")
    source, target = _LINKS[linkname]
    return _call(
        "elink",
        {"dbfrom": source, "db": target, "id": _ids(ids), "linkname": linkname, "retmode": "json"},
        email,
        api_key,
        timeout,
        max_bytes,
        True,
    )


__all__ = [
    "ncbi_fetch_dna_records",
    "ncbi_link_dna_records",
    "ncbi_search_dna_database",
    "ncbi_summarize_dna_records",
]
