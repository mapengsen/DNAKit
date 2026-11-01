"""UCSC, ENCODE and RCSB nucleotide-specific public network operations."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any
from urllib.parse import quote

from dnakit.databases._network_p44 import _params, _request
from dnakit.exceptions import ConfigurationError
from dnakit.search._http import build_url, require_text


def ucsc_track_schema(genome: str, track: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return original UCSC annotation-track column schema."""
    return _request(
        build_url(
            "https://api.genome.ucsc.edu",
            "/list/schema",
            [("genome", require_text(genome, "genome")), ("track", require_text(track, "track"))],
        )
    )


def ucsc_search_assembly(
    genome: str, term: str, *, categories: str | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Search UCSC track descriptions, names and genomic feature data in one assembly."""
    if categories is not None and categories not in {"helpDocs", "publicHubs", "trackDb"}:
        raise ConfigurationError(
            "categories must be helpDocs, publicHubs or trackDb, or None for all."
        )
    return _request(
        build_url(
            "https://api.genome.ucsc.edu",
            "/search",
            [
                ("genome", require_text(genome, "genome")),
                ("search", require_text(term, "term")),
                ("categories", categories),
            ],
        )
    )


_ENCODE_FILTERS = {
    "type",
    "assay_title",
    "assay_term_name",
    "biosample_ontology.term_name",
    "biosample_ontology.classification",
    "organ_slims",
    "system_slims",
    "replicates.library.biosample.donor.organism.scientific_name",
    "status",
    "file_format",
    "file_format_type",
    "output_type",
    "assembly",
    "target.label",
    "accession",
    "searchTerm",
    "limit",
    "field",
    "frame",
    "cart",
    "dataset",
}


def _encode(path: str, filters: Mapping[str, Any] | None, *, text: bool = False) -> dict[str, Any]:
    query = _params(filters, _ENCODE_FILTERS)
    if not text:
        query["format"] = "json"
    return _request(
        build_url("https://www.encodeproject.org", path, query.items()), json_response=not text
    )


def encode_object_metadata(accession: str, *, verbose: bool = True) -> dict[str, Any]:
    """Retrieve ENCODE object metadata by accession or @id, using frame=object."""
    value = require_text(accession, "accession", max_length=256).strip("/")
    if not re.fullmatch(r"(?:[A-Za-z0-9_-]+/)?[A-Za-z0-9_-]+", value):
        raise ConfigurationError(
            "accession must be an ENCODE accession or /collection/accession/ @id."
        )
    return _encode("/" + quote(value, safe="/") + "/", {"frame": "object"})


def encode_tabular_report(
    fields: Sequence[str], *, filters: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Return ENCODE report.tsv for selected fields and a bounded or explicit all-results search."""
    if (
        isinstance(fields, (str, bytes))
        or not isinstance(fields, Sequence)
        or not 1 <= len(fields) <= 100
    ):
        raise ConfigurationError("fields must contain 1-100 field names.")
    query = _params(filters, _ENCODE_FILTERS - {"field"})
    query["field"] = [require_text(x, "field", max_length=256) for x in fields]
    return _encode("/report.tsv", query, text=True)


def encode_experiment_matrix(
    *, filters: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Aggregate ENCODE Experiment counts in the assay x biosample matrix."""
    query = _params(filters, _ENCODE_FILTERS)
    query.setdefault("type", "Experiment")
    return _encode("/matrix/", query)


def encode_body_system_summary(
    *, filters: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Aggregate ENCODE Experiment counts by body system and organ."""
    query = _params(filters, _ENCODE_FILTERS)
    query.setdefault("type", "Experiment")
    return _encode("/summary/", query)


def encode_download_urls(
    *, filters: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Retrieve ENCODE batch-download files.txt URL list for a search or cart.

    This operation intentionally returns a URL list. For actual file bytes use
    dnakit.download.encode.files on the discovered files.
    """
    result = _encode("/batch_download/", filters, text=True)
    result["urls"] = [
        line.strip() for line in result["data"].splitlines() if line.strip().startswith("https://")
    ]
    return result


def _rcsb_options(start: int, rows: int) -> dict[str, Any]:
    if (
        any(isinstance(v, bool) or not isinstance(v, int) for v in (start, rows))
        or start < 0
        or not 1 <= rows <= 10000
    ):
        raise ConfigurationError("start must be non-negative and rows must be 1-10000.")
    return {"paginate": {"start": start, "rows": rows}}


def rcsb_nucleotide_sequence_search(
    sequence: str,
    *,
    sequence_type: str = "dna",
    identity_cutoff: float = 0.0,
    evalue_cutoff: float = 0.1,
    start: int = 0,
    rows: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search PDB DNA or RNA chains by sequence similarity, returning polymer-entity hits.

    RCSB defaults identity=0 and evalue=0.1 are retained; its sequence search
    requires at least 25 residues. Paging is explicit through start/rows.
    """
    seq = require_text(sequence, "sequence", max_length=1_000_000)
    if (
        sequence_type not in {"dna", "rna"}
        or len(seq) < 25
        or not re.fullmatch("[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+", seq)
    ):
        raise ConfigurationError(
            "Provide a DNA/RNA sequence of at least 25 residues and sequence_type dna or rna."
        )
    for name, value in [("identity_cutoff", identity_cutoff), ("evalue_cutoff", evalue_cutoff)]:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise ConfigurationError(f"{name} must be finite.")
    if not 0 <= identity_cutoff <= 1 or evalue_cutoff <= 0:
        raise ConfigurationError("identity_cutoff must be 0-1 and evalue_cutoff must be positive.")
    body = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": seq,
                "sequence_type": sequence_type,
                "identity_cutoff": identity_cutoff,
                "evalue_cutoff": evalue_cutoff,
            },
        },
        "return_type": "polymer_entity",
        "request_options": _rcsb_options(start, rows),
    }
    return _request("https://search.rcsb.org/rcsbsearch/v2/query", body=body)


def rcsb_nucleotide_motif_search(
    pattern: str,
    *,
    sequence_type: str = "dna",
    pattern_type: str = "simple",
    start: int = 0,
    rows: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search PDB DNA/RNA chain motifs in simple, prosite or regex notation."""
    if sequence_type not in {"dna", "rna"} or pattern_type not in {"simple", "prosite", "regex"}:
        raise ConfigurationError("Unsupported sequence_type or pattern_type.")
    body = {
        "query": {
            "type": "terminal",
            "service": "seqmotif",
            "parameters": {
                "value": require_text(pattern, "pattern", max_length=4096),
                "sequence_type": sequence_type,
                "pattern_type": pattern_type,
            },
        },
        "return_type": "polymer_entity",
        "request_options": _rcsb_options(start, rows),
    }
    return _request("https://search.rcsb.org/rcsbsearch/v2/query", body=body)


__all__ = [
    "encode_body_system_summary",
    "encode_download_urls",
    "encode_experiment_matrix",
    "encode_object_metadata",
    "encode_tabular_report",
    "rcsb_nucleotide_motif_search",
    "rcsb_nucleotide_sequence_search",
    "ucsc_search_assembly",
    "ucsc_track_schema",
]
