"""Average nucleotide identity through PyFastANI's in-process FastANI bindings.

This adapter does not reimplement FastANI or invoke its command-line program.
PyFastANI is installed separately under MIT; its bundled FastANI code is Apache-2.0.
API reference: https://pyfastani.readthedocs.io/en/stable/api/ani.html
FastANI method: https://doi.org/10.1038/s41467-018-07641-9
Initial version: directional ANI for sequences/FASTA, with stage progress and
validated JSON output, without a separately installed FastANI executable.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

from dnakit._progress import Progress
from dnakit.bioinformatics._common import positive, require, sequence_text
from dnakit.exceptions import BackendExecutionError, ConfigurationError, InputFormatError
from dnakit.io import read

_IUPAC = frozenset("ACGTRYSWKMBDHVN")
_FRAGMENT_LENGTH = 3000


def _genome(value: str | Path, label: str) -> tuple[list[str], dict[str, Any]]:
    """Resolve explicitly distinguishable DNA and FASTA paths, retaining contigs."""
    if not isinstance(value, (str, Path)) or (isinstance(value, str) and not value):
        raise ConfigurationError(f"{label} must be non-empty DNA or a FASTA path.")
    # Interpret alphabet-only strings as sequences before touching the filesystem:
    # chromosome-sized inputs exceed filesystem filename limits.
    if isinstance(value, str) and set(value.upper()) <= _IUPAC:
        contigs = [sequence_text(value)]
        source: str | None = None
    else:
        try:
            path = Path(value).expanduser().resolve()
            if not path.is_file():
                raise ConfigurationError(f"{label} must be IUPAC DNA or an existing FASTA file.")
        except (OSError, ValueError) as exc:
            raise ConfigurationError(f"{label} must be IUPAC DNA or a valid FASTA path.") from exc
        with read(path, format="fasta", verbose=False) as records:
            contigs = [sequence_text(record) for record in records]
        if not contigs:
            raise InputFormatError(f"{label} FASTA contains no sequences.")
        source = str(path)

    digest = hashlib.sha256()
    ambiguous = 0
    normalized = []
    for contig in contigs:
        # Delimit lengths as well as contents, so different contig boundaries
        # produce different fingerprints even when concatenated DNA is identical.
        digest.update(f"{len(contig)}:".encode("ascii"))
        digest.update(contig.encode("ascii"))
        ambiguous += sum(base not in "ACGT" for base in contig)
        normalized.append("".join(base if base in "ACGT" else "N" for base in contig))
    return normalized, {
        "source": source,
        "contigs": len(contigs),
        "length": sum(map(len, contigs)),
        "input_sha256": digest.hexdigest(),
        "ambiguous_bases_to_n": ambiguous,
    }


def average_nucleotide_identity(
    query: str | Path,
    reference: str | Path,
    *,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate directional genome ANI with the optional PyFastANI library.

    Install with ``pip install 'dnakit[genomics]'``. Each input is a non-empty,
    ungapped IUPAC DNA string or a local FASTA path (gzip is supported by DNAKit
    I/O). A string consisting only of IUPAC letters is always DNA; use a
    ``Path`` object or an explicit ``./`` prefix for such a filename. FASTA
    text itself is not an input format. Every FASTA record is a contig of ONE
    genome; contigs are never concatenated. Inputs are loaded into memory.
    Ambiguous bases are represented as N without removing positions.

    FastANI's standard settings are retained: 3,000-base query fragments,
    16-mers, 80-percent identity tuning, and a 0.2 shared-genome reporting
    threshold. Contigs shorter than a query fragment and incomplete terminal
    fragments cannot contribute query fragments. Reference contigs shorter
    than the backend's sketching window are skipped by PyFastANI.

    ANI is a percentage for the reported homologous fragments, not an identity
    across every input base. ``query_aligned_fraction`` is matching fragments
    divided by all complete query fragments, not base-level genome coverage.
    The comparison runs query -> reference only and can be asymmetric; call
    again with swapped inputs if both directions are needed. This API targets
    related prokaryotic genomes, not short-read alignment or distant genomes.

    ``status='no_match'`` means the backend reported no estimate (or no query
    fragment exists); ANI, matching counts and aligned fraction are then None.
    It does not mean zero-percent identity. Invalid or non-finite backend
    results raise ``BackendExecutionError``. Output is JSON serializable.
    Progress reports loading, indexing and mapping stages, not base-level work.
    """
    positive(threads, "threads")
    if threads > 1024:
        raise ConfigurationError("threads must not exceed 1024.")
    if not isinstance(verbose, bool):
        raise ConfigurationError("verbose must be a boolean.")

    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Compute genome ANI", total=3)
        query_contigs, query_info = _genome(query, "query")
        reference_contigs, reference_info = _genome(reference, "reference")
        progress.advance(task)
        pyfastani = require("pyfastani", "genomics")
        fragments = sum(len(contig) // _FRAGMENT_LENGTH for contig in query_contigs)
        try:
            sketch = pyfastani.Sketch()
            sketch.add_draft("reference", reference_contigs)
            mapper = sketch.index()
            progress.advance(task)
            reference_has_minimizers = bool(len(mapper.minimizers))
            hits = (
                mapper.query_draft(query_contigs, threads=threads)
                if fragments and reference_has_minimizers
                else []
            )
            progress.advance(task)
        except (TypeError, ValueError, OverflowError, RuntimeError) as exc:
            raise BackendExecutionError(f"PyFastANI genome comparison failed: {exc}") from exc

    result: dict[str, Any] = {
        "status": "no_match",
        "ani_percent": None,
        "matching_fragments": None,
        "query_fragments": fragments,
        "query_aligned_fraction": None,
        "direction": "query_to_reference",
        "query": query_info,
        "reference": reference_info,
        "backend": "pyfastani",
        "backend_version": str(pyfastani.__version__),
        "parameters": {
            "threads": threads,
            "fragment_length": int(mapper.fragment_length),
            "k": int(mapper.k),
            "minimum_fraction": float(mapper.minimum_fraction),
            "percentage_identity": float(mapper.percentage_identity),
        },
        "citations": [
            "https://doi.org/10.1038/s41467-018-07641-9",
            "https://doi.org/10.1093/nargab/lqaf095",
        ],
    }
    if not hits:
        result["reason"] = (
            "no_complete_query_fragments"
            if not fragments
            else "no_reference_minimizers"
            if not reference_has_minimizers
            else "no_reportable_ani"
        )
        return result
    if len(hits) != 1 or hits[0].name != "reference":
        raise BackendExecutionError("PyFastANI returned unexpected reference hits.")
    hit = hits[0]
    identity = float(hit.identity)
    if (
        not math.isfinite(identity)
        or not 0 <= identity <= 100
        or type(hit.matches) is not int
        or type(hit.fragments) is not int
        or hit.fragments != fragments
        or not 0 < hit.matches <= fragments
    ):
        raise BackendExecutionError("PyFastANI returned an invalid ANI or fragment count.")
    result.update(
        status="matched",
        ani_percent=identity,
        matching_fragments=hit.matches,
        query_aligned_fraction=hit.matches / fragments,
    )
    return result


__all__ = ["average_nucleotide_identity"]
