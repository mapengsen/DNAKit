"""Offline genomic-benchmarks helper: find FASTA sequences on reference contigs.

``locate_sequences_in_reference`` calls ``genomic_benchmarks.seq2loc.fasta2loc``, which
matches exact sequences and their reverse complements against reference contigs. Dataset
downloads, the cloud cache and torch dataset classes need the network and are not wrapped.
Local dataset summaries and cache checks are in ``dnakit.predictions.dna_model_backends``.
"""

from __future__ import annotations

import contextlib
import io
from collections.abc import Mapping
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import _integer, _provenance, _rejected, _text

_GENOMIC_BENCHMARKS = "genomic-benchmarks"
_EXTRAS = "bio-genomic-benchmarks"
_DNA = frozenset("ACGTN")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_MAX_CONTIG_NAME = 128
_MAX_REFERENCE_CONTIGS = 1_000
_MAX_REFERENCE_LENGTH = 10_000_000


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _contigs(reference: object, *, limit: int) -> dict[str, str]:
    """Validate a mapping of contig name to uppercase DNA and return it with sorted keys."""

    if not isinstance(reference, Mapping) or not reference:
        raise ConfigurationError("reference must be a non-empty mapping of contig name to DNA.")
    if len(reference) > _MAX_REFERENCE_CONTIGS:
        raise ConfigurationError(f"reference may contain at most {_MAX_REFERENCE_CONTIGS} contigs.")
    total = 0
    contigs: dict[str, str] = {}
    for name in sorted(reference, key=str):
        contig = _text(name, "reference contig name", maximum=_MAX_CONTIG_NAME)
        sequence = reference[name]
        if not isinstance(sequence, str) or not sequence:
            raise ConfigurationError(f"reference[{contig!r}] must be a non-empty DNA string.")
        text = sequence.upper()
        if set(text) - _DNA:
            raise ConfigurationError(f"reference[{contig!r}] may contain only A, C, G, T and N.")
        total += len(text)
        if total > limit:
            raise ConfigurationError(f"reference is longer than max_reference_length={limit}.")
        contigs[contig] = text
    return contigs


def _overlapping_count(haystack: str, needle: str) -> int:
    count = 0
    start = haystack.find(needle)
    while start != -1:
        count += 1
        start = haystack.find(needle, start + 1)
    return count


def locate_sequences_in_reference(
    fasta_path: str,
    reference: Mapping[str, str],
    *,
    max_reference_length: int = 2_000_000,
    max_sequences: int = 100,
    max_query_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate exact FASTA sequences (or their reverse complements) on reference contigs.

    Each sequence is reported once, at the last match found while scanning contigs in name
    order from left to right. ``occurrences`` counts every exact match, so a value above 1
    marks a location that is not unique.
    """

    path = input_path(fasta_path)
    limit = _integer(
        max_reference_length, "max_reference_length", minimum=1, maximum=_MAX_REFERENCE_LENGTH
    )
    contigs = _contigs(reference, limit=limit)
    count_limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=1_000)
    query_limit = _integer(max_query_length, "max_query_length", minimum=1, maximum=1_000_000)
    seqio = require("Bio.SeqIO", _EXTRAS)
    queries: dict[str, str] = {}
    total = 0
    for record in seqio.parse(str(path), "fasta"):
        name = str(record.id)
        text = str(record.seq)
        if name in queries:
            raise ConfigurationError(f"FASTA sequence name {name!r} is used more than once.")
        if not text or set(text) - _DNA:
            raise ConfigurationError(f"FASTA sequence {name!r} must use uppercase A, C, G, T or N.")
        queries[name] = text
        total += len(text)
        if len(queries) > count_limit:
            raise ConfigurationError(f"FASTA has more than max_sequences={count_limit} sequences.")
        if total > query_limit:
            raise ConfigurationError(f"FASTA is longer than max_query_length={query_limit}.")
    if not queries:
        raise ConfigurationError(f"no FASTA sequences found in {path.name}.")
    # upstream keys its results by name and stores one terminal per trie node, so two queries
    # that are equal or reverse complements of each other would overwrite one another.
    canonical: dict[str, str] = {}
    for name in sorted(queries):
        text = queries[name]
        key = min(text, _reverse_complement(text))
        if key in canonical:
            raise ConfigurationError(
                f"FASTA sequences {canonical[key]!r} and {name!r} are equal or reverse "
                "complements; locate_sequences_in_reference cannot report both."
            )
        canonical[key] = name
    upstream = require("genomic_benchmarks.seq2loc", _EXTRAS)
    printed = io.StringIO()
    try:
        with contextlib.redirect_stdout(printed), contextlib.redirect_stderr(io.StringIO()):
            located = upstream.fasta2loc(str(path), contigs, use_seq_ids=True)
    except (ValueError, KeyError, OSError) as exc:
        raise _rejected(_GENOMIC_BENCHMARKS, exc) from exc
    matches = []
    for name in sorted(located):
        chrom, start, end, strand = located[name]
        occurrences = sum(_overlapping_count(seq, queries[name]) for seq in contigs.values())
        reverse = _reverse_complement(queries[name])
        if reverse != queries[name]:
            occurrences += sum(_overlapping_count(seq, reverse) for seq in contigs.values())
        matches.append(
            {
                "sequence": name,
                "chromosome": str(chrom),
                "start": int(start),
                "end": int(end),
                "strand": str(strand),
                "occurrences": occurrences,
                "ambiguous": occurrences > 1,
            }
        )
    unmatched = sorted(name for name in queries if name not in located)
    messages = [line.strip() for line in printed.getvalue().splitlines() if line.strip()]
    return {
        "status": "completed",
        "operation": "locate_sequences_in_reference",
        "query_count": len(queries),
        "contig_count": len(contigs),
        "matched_count": len(matches),
        "matches": matches,
        "unmatched": unmatched,
        "messages": messages if verbose else [],
        "provenance": _provenance(
            _GENOMIC_BENCHMARKS,
            coordinates="0-based, half-open",
            strand="'-' means the reverse complement of the query matched",
            network="not used",
        ),
    }


__all__ = ["locate_sequences_in_reference"]

_register_public_verbose(globals(), tuple(__all__))
