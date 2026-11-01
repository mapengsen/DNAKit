"""gtars transcript-store operations (integration pack P24, part 3).

A transcript store is one binary file built once from transcript records, cdot JSON files and a MANE
summary. These functions build such a store, look transcripts up and convert transcript (c. and n.)
positions to genomic positions with ``gtars.reftx``. Chromosomes of a store are refget sequence
digests (base64url, 32 characters, optionally written with the ``SQ.`` prefix). gtars is imported
only when a function runs; a missing package raises ``BackendUnavailableError`` with the install
hint."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.gtars_functions import (
    _choice,
    _integer,
    _new_output,
    _plain,
    _provenance,
    _text,
    _upstream,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-gtars"
_COORDINATE_TYPES = ("c", "n")
_MAX_TRANSCRIPTS = 100_000


def _reftx() -> Any:
    return require("gtars.reftx", _EXTRA)


def _open_store(store: str | Path) -> Any:
    reftx = _reftx()
    path = input_path(store)
    return _upstream(lambda: reftx.ReadonlyTxStore.open(str(path)), "open the transcript store")


def _prov(note: str) -> dict[str, str]:
    return _provenance("c. and n. numbering as defined by gtars; " + note)


def _record(transcript: Any) -> dict[str, Any]:
    return _plain(transcript.to_dict())


def transcript_coordinate_to_genomic(
    transcript_store: str | Path,
    position: int,
    *,
    accession: str | None = None,
    gene: str | None = None,
    coordinate_type: str = "c",
    datum: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a transcript position to a genomic position with gtars ``CoordinateMapper``.

    ``coordinate_type="c"`` converts a coding (c.) position of ``accession`` (``c_to_g_full``) or,
    when ``gene`` is given instead, of that gene's MANE Select transcript (``c_to_g_by_gene``);
    ``coordinate_type="n"`` converts a non-coding (n.) position of ``accession`` (``n_to_g_full``).
    ``datum`` is the optional gtars datum argument for c. positions. The result holds the genomic
    position, the chromosome accession and the strand."""

    kind = _choice(coordinate_type, "coordinate_type", _COORDINATE_TYPES)
    pos = _integer(position, "position", -(10**9), 10**10)
    ref = None if datum is None else _integer(datum, "datum", -(10**9), 10**10)
    if (accession is None) == (gene is None):
        raise ConfigurationError("give exactly one of accession and gene.")
    if gene is not None and kind != "c":
        raise ConfigurationError("a gene can only be used with coordinate_type='c'.")
    reftx = _reftx()
    store = _open_store(transcript_store)
    mapper = reftx.CoordinateMapper(store)
    if gene is not None:
        name = _text(gene, "gene", maximum=128)
        raw = _upstream(lambda: mapper.c_to_g_by_gene(name, pos, ref), "map the gene position")
    else:
        acc = _text(accession, "accession", maximum=128)
        if kind == "c":
            raw = _upstream(lambda: mapper.c_to_g_full(acc, pos, ref), "map the c. position")
        else:
            raw = _upstream(lambda: mapper.n_to_g_full(acc, pos), "map the n. position")
    row = {str(key): _plain(value) for key, value in raw.items()}
    row["strand"] = (
        str(raw["strand"]).replace("Strand.", "").replace("Plus", "+").replace("Minus", "-")
    )
    return {
        "status": "completed",
        "operation": "transcript_coordinate_to_genomic",
        "coordinate_type": kind,
        "position": pos,
        **row,
        "provenance": _prov(
            "gene mode uses the MANE Select transcript and needs a store built with MANE flags"
        ),
    }


def transcript_record_lookup(
    transcript_store: str | Path,
    accession: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up a transcript record (gene, chromosome, strand, CDS span, exons, MANE flags) in a
    transcript store.

    Also returns the chromosome accession (``get_chrom_accession``) and strand (``get_strand``) of
    the transcript. ``found`` is false when the accession is not in the store."""

    acc = _text(accession, "accession", maximum=128)
    store = _open_store(transcript_store)
    transcript = store.lookup(acc)
    if transcript is None:
        return {
            "status": "completed",
            "operation": "transcript_record_lookup",
            "accession": acc,
            "found": False,
            "transcript_count": len(store),
            "provenance": _prov("accession is exact, including its version"),
        }
    strand = store.get_strand(acc)
    return {
        "status": "completed",
        "operation": "transcript_record_lookup",
        "accession": acc,
        "found": True,
        "transcript": _record(transcript),
        "chromosome_accession": store.get_chrom_accession(acc),
        "strand": None if strand is None else str(strand),
        "transcript_count": len(store),
        "provenance": _prov("accession is exact, including its version"),
    }


def transcript_store_build(
    output: str | Path,
    *,
    transcripts: Sequence[Mapping[str, Any]] | None = None,
    cdot_files: Sequence[str | Path] | None = None,
    mane_summary: str | Path | None = None,
    chrom_mappings: Mapping[str, str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a gtars transcript store file from transcript records, cdot JSON files and a MANE
    summary.

    ``transcripts`` are dicts with accession, gene, chrom (refget digest), strand (``+`` or
    ``-``), cds_start, cds_end, exons ([{start, end}]) and mane (null or {select,
    plus_clinical}); ``TxStoreBuilder.add_transcript``. ``cdot_files`` are cdot JSON files
    (``load_cdot``, gzip detected; strand 1 or -1). ``chrom_mappings`` maps contig names such as
    ``NC_000001.11`` to refget digests (``add_chrom_mapping``) and must come with cdot files.
    ``mane_summary`` is an NCBI MANE summary TSV (``load_mane_summary``). The output file must
    not exist."""

    target = _new_output(output)
    if transcripts is None and not cdot_files:
        raise ConfigurationError("give transcripts and/or cdot_files.")
    if transcripts is not None and (
        isinstance(transcripts, (str, bytes)) or len(transcripts) > _MAX_TRANSCRIPTS
    ):
        raise ConfigurationError(f"transcripts must be a list of at most {_MAX_TRANSCRIPTS} dicts.")
    cdot = [input_path(item) for item in (cdot_files or [])]
    mane = None if mane_summary is None else input_path(mane_summary)
    mappings = {}
    for name, digest in (chrom_mappings or {}).items():
        mappings[_text(name, "contig name", maximum=128)] = _text(
            digest, "chromosome digest", maximum=128
        )
    reftx = _reftx()
    builder = reftx.TxStoreBuilder()
    for name, digest in mappings.items():
        _upstream(
            lambda n=name, d=digest: builder.add_chrom_mapping(n, d),
            "register a chromosome mapping",
        )
    added = 0
    for record in transcripts or []:
        if not isinstance(record, Mapping):
            raise ConfigurationError("each transcript must be a dict.")
        _upstream(lambda r=dict(record): builder.add_transcript(r), "add a transcript")
        added += 1
    loaded = 0
    for path in cdot:
        loaded += int(_upstream(lambda p=path: builder.load_cdot(str(p)), f"read {path.name}"))
    mane_rows = 0
    if mane is not None:
        mane_rows = int(
            _upstream(lambda: builder.load_mane_summary(str(mane)), "read the MANE summary")
        )
    total = len(builder)
    _upstream(lambda: builder.build(str(target)), "write the transcript store")
    return {
        "status": "completed",
        "operation": "transcript_store_build",
        "output": str(target),
        "transcripts_added": added,
        "cdot_transcripts_loaded": loaded,
        "mane_rows_loaded": mane_rows,
        "chromosome_mappings": len(mappings),
        "transcript_count": total,
        "provenance": _provenance("store is a single binary file; chromosomes are refget digests"),
    }


__all__ = [
    "transcript_coordinate_to_genomic",
    "transcript_record_lookup",
    "transcript_store_build",
]

_register_public_verbose(globals(), tuple(__all__))
