"""ENA browser, portal, cross-reference and CRAM reference operations.

The browser routes and query names follow its official OpenAPI v3 schema.
Downloads stream to atomic files and inherit DNAKit byte limits and progress.
"""

from __future__ import annotations

import re
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote

from dnakit.databases._network_p44 import _request
from dnakit.download.files import download_file
from dnakit.download.models import DownloadConfig, RemoteFile
from dnakit.exceptions import ConfigurationError
from dnakit.search._http import build_url, require_text

_BROWSER = "https://www.ebi.ac.uk/ena/browser/api"
_FORMATS = {"embl", "text", "fasta", "xml"}
_RESULTS = {
    "sequence",
    "coding",
    "noncoding",
    "wgs_set",
    "tsa_set",
    "tls_set",
    "read_run",
    "read_experiment",
    "analysis",
    "assembly",
}


def _format(value: str) -> str:
    if value not in _FORMATS:
        raise ConfigurationError("format must be embl, text, fasta or xml.")
    return value


def _download(url: str, output_path: str, overwrite: bool, verbose: bool) -> dict[str, Any]:
    record = download_file(
        RemoteFile(url, filename=Path(output_path).name),
        output_path,
        config=DownloadConfig(overwrite=overwrite, show_progress=verbose),
        verbose=verbose,
    )
    return {
        "status": "completed",
        "file": asdict(record),
        "provenance": {"service": "ENA", "retrieved_at": datetime.now(timezone.utc).isoformat()},
    }


def ena_download_record(
    accession: str,
    output_path: str,
    *,
    format: str = "embl",
    full_contig_set: bool = False,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download EMBL contig-set records or XML assembly records by accession.

    Use ``full_contig_set=True`` for all WGS/TSA/TLS sequences; the default
    retains ENA's master-record behaviour. RefSeq identifiers require NCBI.
    """
    fmt = _format(format)
    value = quote(require_text(accession, "accession", max_length=100), safe="")
    if not isinstance(full_contig_set, bool):
        raise ConfigurationError("full_contig_set must be boolean.")
    if full_contig_set and fmt not in {"embl", "text"}:
        raise ConfigurationError("full_contig_set applies to EMBL/text output.")
    params = [("set", True)] if full_contig_set else []
    return _download(
        build_url(_BROWSER, f"/{fmt}/{value}", params), output_path, overwrite, verbose
    )


def ena_download_search(
    query: str,
    result: str,
    output_path: str,
    *,
    format: str = "embl",
    limit: int = 100000,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download advanced-search EMBL/FASTA sequences or XML run/experiment/analysis/assembly
    records.

    ``result`` selects sequence, coding, noncoding, wgs_set, tsa_set, tls_set, read_run,
    read_experiment, analysis or assembly. ENA default limit is 100000; limit=0 requests all
    matches subject to the streaming byte limit.
    """
    fmt = _format(format)
    if result not in _RESULTS:
        raise ConfigurationError("Unsupported DNA result type.")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ConfigurationError("limit must be a non-negative integer.")
    if result in {"read_run", "read_experiment", "analysis", "assembly"} and fmt != "xml":
        raise ConfigurationError("Sequencing/assembly metadata search requires XML format.")
    if result.endswith("_set") and fmt not in {"embl", "text"}:
        raise ConfigurationError("Contig-set metadata search requires EMBL format.")
    url = build_url(
        _BROWSER,
        f"/{fmt}/search",
        [("result", result), ("query", require_text(query, "query")), ("limit", limit)],
    )
    return _download(url, output_path, overwrite, verbose)


def ena_download_linked_records(
    link_type: str,
    accession: str,
    result: str,
    output_path: str,
    *,
    format: str = "embl",
    subtree: bool = False,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download EMBL/FASTA/XML records linked to a study, sample or taxon accession."""
    fmt = _format(format)
    if link_type not in {"study", "sample", "taxon"} or result not in _RESULTS:
        raise ConfigurationError("Unsupported ENA link type or DNA result type.")
    if not isinstance(subtree, bool):
        raise ConfigurationError("subtree must be boolean.")
    url = build_url(
        _BROWSER,
        f"/{fmt}/links/{link_type}",
        [
            ("accession", require_text(accession, "accession")),
            ("result", result),
            ("subtree", subtree),
        ],
    )
    return _download(url, output_path, overwrite, verbose)


def ena_download_keyword_results(
    query: str,
    domain: str,
    result: str,
    output_path: str,
    *,
    format: str = "embl",
    limit: int = 100000,
    offset: int = 0,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download keyword-search EMBL/FASTA/XML records, including WGS/TSA/TLS master sets.

    ``domain`` is the upstream EBI Search domain (for example ``embl`` or
    ``sra-run``); ``result`` selects its ENA dataset. Defaults match upstream.
    """
    fmt = _format(format)
    if result not in _RESULTS:
        raise ConfigurationError("Unsupported DNA result type.")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in (limit, offset)):
        raise ConfigurationError("limit/offset must be non-negative integers.")
    url = build_url(
        _BROWSER,
        f"/{fmt}/textsearch",
        [
            ("query", require_text(query, "query")),
            ("domain", require_text(domain, "domain")),
            ("result", result),
            ("limit", limit),
            ("offset", offset),
        ],
    )
    return _download(url, output_path, overwrite, verbose)


def ena_count_keyword_results(
    query: str, domain: str, result: str, *, format: str = "embl", verbose: bool = True
) -> dict[str, Any]:
    """Count ENA keyword-search hits for a record format and domain."""
    fmt = _format(format)
    if result not in _RESULTS:
        raise ConfigurationError("Unsupported DNA result type.")
    return _request(
        build_url(
            _BROWSER,
            f"/{fmt}/textsearch/count",
            [
                ("query", require_text(query, "query")),
                ("domain", require_text(domain, "domain")),
                ("result", result),
            ],
        ),
        json_response=False,
    )


def ena_linked_accessions(
    link_type: str, accession: str, result: str, *, verbose: bool = True
) -> dict[str, Any]:
    """List portal DNA record accessions linked to a study, sample or taxon."""
    if link_type not in {"study", "sample", "taxon"} or result not in _RESULTS:
        raise ConfigurationError("Unsupported ENA link type or DNA result type.")
    return _request(
        build_url(
            "https://www.ebi.ac.uk/ena/portal/api",
            "/links/" + link_type,
            [
                ("accession", require_text(accession, "accession")),
                ("result", result),
                ("format", "tsv"),
            ],
        ),
        json_response=False,
    )


def ena_search_cross_references(
    source: str,
    *,
    source_accession: str | None = None,
    target: str | None = None,
    target_accession: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search ENA cross-references to an external database and return original TSV."""
    params = [("source", require_text(source, "source"))]
    for name, value in [
        ("source_accession", source_accession),
        ("target", target),
        ("accession", target_accession),
    ]:
        if value is not None:
            params.append((name, require_text(value, name)))
    return _request(
        build_url("https://www.ebi.ac.uk/ena/xref/rest", "/tsv/search", params), json_response=False
    )


def ena_download_cram_reference(
    checksum: str,
    output_path: str,
    *,
    algorithm: str = "md5",
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download exact reference sequence bytes by CRAM registry MD5 or SHA1 checksum."""
    sizes = {"md5": 32, "sha1": 40}
    if (
        algorithm not in sizes
        or not isinstance(checksum, str)
        or not re.fullmatch("[a-fA-F0-9]{" + str(sizes.get(algorithm, 0)) + "}", checksum)
    ):
        raise ConfigurationError("checksum must be a hexadecimal MD5 (32) or SHA1 (40).")
    return _download(
        f"https://www.ebi.ac.uk/ena/cram/{algorithm}/{checksum.lower()}",
        output_path,
        overwrite,
        verbose,
    )


__all__ = [
    "ena_count_keyword_results",
    "ena_download_cram_reference",
    "ena_download_keyword_results",
    "ena_download_linked_records",
    "ena_download_record",
    "ena_download_search",
    "ena_linked_accessions",
    "ena_search_cross_references",
]
