"""pyensembl annotation queries (pack P36): features at a locus, ID look-ups, relations, sequences, loci.

Each function builds a ``pyensembl.Genome`` from local files (a GTF and optional genome, transcript or protein
FASTA) and calls the pyensembl method that carries the operation. Nothing is downloaded, except the explicit
``download=True`` of ``load_hgnc_gene_aliases``. The GTF is indexed once; the index is cached under the
``cache_dir`` folder (pyensembl's default cache when ``None``) with a name that carries the GTF's content
digest, so an edited GTF gets a fresh index. Coordinates are 1-based inclusive, as in pyensembl.
Missing pyensembl raises ``BackendUnavailableError`` with the ``bio-pyensembl`` hint.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-pyensembl"
_REFERENCE = "dnakit_custom"
_MAX_RESULTS = 1_000_000
_FEATURES = ("gene", "transcript", "exon", "protein")
_FROM_KINDS = (
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_name",
    "exon_id",
    "protein_id",
)
_TO_KINDS = ("gene", "transcript", "exon")
_ERRORS = (ValueError, KeyError, OSError, TypeError, AttributeError, LookupError)


def _provenance(operation: str) -> dict[str, Any]:
    return {
        "package": "pyensembl",
        "version": importlib.metadata.version("pyensembl"),
        "operation": operation,
        "coordinates": "1-based inclusive",
    }


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > 256
        or any(c in value for c in "\x00\r\n\t")
    ):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value.strip()


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return str(value)


def _strand(value: object) -> str | None:
    return None if value is None else _choice(value, "strand", ("+", "-"))


def _digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()[:16]


def _cache(cache_dir: object) -> str | None:
    if cache_dir is None:
        return None
    if not isinstance(cache_dir, (str, Path)):
        raise ConfigurationError("cache_dir must be a folder path.")
    return str(Path(cache_dir).expanduser().resolve())


def _genome(
    gtf: object = None,
    *,
    genome_fasta: object = None,
    transcript_fasta: object = None,
    protein_fasta: object = None,
    cache_dir: object = None,
    index: bool = True,
) -> Any:
    pyensembl = require("pyensembl", _EXTRA)
    gtf_path = None if gtf is None else input_path(gtf)  # type: ignore[arg-type]
    kwargs: dict[str, Any] = {
        "reference_name": _REFERENCE,
        "annotation_name": "gtf_" + (_digest(gtf_path) if gtf_path else "sequence_only"),
        "cache_directory_path": _cache(cache_dir),
    }
    if gtf_path is not None:
        kwargs["gtf_path_or_url"] = str(gtf_path)
    if genome_fasta is not None:
        kwargs["genome_fasta_path_or_url"] = str(input_path(genome_fasta))  # type: ignore[arg-type]
    if transcript_fasta is not None:
        kwargs["transcript_fasta_paths_or_urls"] = [str(input_path(transcript_fasta))]  # type: ignore[arg-type]
    if protein_fasta is not None:
        kwargs["protein_fasta_paths_or_urls"] = [str(input_path(protein_fasta))]  # type: ignore[arg-type]
    try:
        genome = pyensembl.Genome(**kwargs)
        if index and gtf_path is not None:
            genome.index()
    except _ERRORS as exc:
        raise BackendExecutionError(f"pyensembl could not load the annotation: {exc}") from exc
    return genome


def _call(operation: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except _ERRORS as exc:
        raise BackendExecutionError(f"pyensembl {operation} failed: {exc}") from exc


def _locus(item: Any) -> dict[str, Any]:
    return {
        "contig": item.contig,
        "start": int(item.start),
        "end": int(item.end),
        "strand": item.strand,
    }


def _record(kind: str, item: Any) -> dict[str, Any]:
    if kind == "gene":
        return {
            "gene_id": item.gene_id,
            "gene_name": item.gene_name,
            "biotype": item.biotype,
            **_locus(item),
        }
    if kind == "transcript":
        return {
            "transcript_id": item.transcript_id,
            "transcript_name": item.transcript_name,
            "gene_id": item.gene_id,
            "biotype": item.biotype,
            **_locus(item),
        }
    if kind == "exon":
        return {
            "exon_id": item.exon_id,
            "gene_id": item.gene_id,
            "gene_name": item.gene_name,
            **_locus(item),
        }
    raise ConfigurationError(f"Unsupported feature {kind}.")


def _limit(value: object) -> int:
    return _int(value, "max_results", 1, _MAX_RESULTS)


def _cut(items: Sequence[Any], limit: int) -> tuple[list[Any], bool]:
    return list(items[:limit]), len(items) > limit


def load_hgnc_gene_aliases(
    hgnc_tsv: str | Path | None = None,
    *,
    download: bool = False,
    source_url: str | None = None,
    cache_directory: str | Path | None = None,
    names: Sequence[str] | None = None,
    max_names: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load the HGNC gene-name alias table (symbols, aliases, previous symbols) mapped to Ensembl gene IDs.

    With ``hgnc_tsv`` the file is read offline (``GeneNameAliases.from_hgnc``; plain or gzip complete-set
    TSV). With ``download=True`` the current HGNC complete set (or ``source_url``) is fetched and cached
    (``GeneNameAliases.download_hgnc``; network). ``names`` restricts the result to those symbols; otherwise
    at most ``max_names`` entries are returned and ``total_names`` gives the full size.
    """

    if (hgnc_tsv is None) == (not download):
        raise ConfigurationError("Give hgnc_tsv or set download=True, not both and not neither.")
    limit = _int(max_names, "max_names", 1, _MAX_RESULTS)
    wanted = None if names is None else [_text(n, "names item") for n in names]
    pyensembl = require("pyensembl", _EXTRA)
    try:
        if hgnc_tsv is not None:
            table = pyensembl.GeneNameAliases.from_hgnc(str(input_path(hgnc_tsv)))
            origin = "local file"
        else:
            kwargs: dict[str, Any] = {}
            if source_url is not None:
                kwargs["source_url"] = _text(source_url, "source_url")
            if cache_directory is not None:
                kwargs["cache_directory_path"] = _cache(cache_directory)
            table = pyensembl.GeneNameAliases.download_hgnc(**kwargs)
            origin = "download (cached)"
    except _ERRORS as exc:
        raise BackendExecutionError(f"pyensembl HGNC alias loading failed: {exc}") from exc
    keys = sorted(table)[:limit] if wanted is None else [n for n in wanted if n in table]
    return {
        "status": "completed",
        "operation": "load_hgnc_gene_aliases",
        "origin": origin,
        "source": table.source,
        "species": table.species,
        "total_names": len(table),
        "aliases": {k: list(table[k]) for k in keys},
        "missing": [] if wanted is None else [n for n in wanted if n not in table],
        "provenance": _provenance("GeneNameAliases"),
    }


def fetch_genomic_sequence(
    genome_fasta: str | Path,
    contig: str,
    start: int,
    end: int,
    *,
    strand: str = "+",
    mask: str = "upper",
    verbose: bool = True,
) -> dict[str, Any]:
    """Read genomic DNA for a 1-based inclusive interval from a local genome FASTA (``Genome.sequence``).

    ``strand="-"`` returns the reverse complement. ``mask="upper"`` upper-cases the result, ``mask="raw"``
    keeps soft masking (pyensembl's default is ``upper``). The contig may be named as in the FASTA.
    """

    chrom = _text(contig, "contig")
    first = _int(start, "start", 1, 2**31 - 1)
    last = _int(end, "end", first, 2**31 - 1)
    if last - first + 1 > 50_000_000:
        raise ConfigurationError("The interval is longer than 50 Mb.")
    genome = _genome(genome_fasta=genome_fasta, index=False)
    sequence = _call(
        "sequence",
        genome.sequence,
        chrom,
        first,
        last,
        mask=_choice(mask, "mask", ("upper", "raw")),
        strand=_choice(strand, "strand", ("+", "-")),
    )
    return {
        "status": "completed",
        "operation": "fetch_genomic_sequence",
        "contig": chrom,
        "start": first,
        "end": last,
        "strand": strand,
        "length": len(sequence),
        "sequence": str(sequence),
        "provenance": _provenance("Genome.sequence"),
    }


def fetch_transcript_sequence(
    transcript_fasta: str | Path,
    transcript_ids: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up the cDNA sequences of transcript IDs in a local transcript FASTA (``Genome.transcript_sequence``).

    IDs absent from the FASTA map to ``null`` in ``sequences`` and are listed in ``missing``. A versioned ID
    such as ``ENST00000269305.8`` is accepted when the FASTA records that version.
    """

    if (
        isinstance(transcript_ids, (str, bytes))
        or not isinstance(transcript_ids, Sequence)
        or not transcript_ids
    ):
        raise ConfigurationError("transcript_ids must be a non-empty list.")
    ids = [_text(i, "transcript_ids item") for i in transcript_ids]
    genome = _genome(transcript_fasta=transcript_fasta, index=False)
    sequences: dict[str, str | None] = {}
    for transcript_id in ids:
        value = _call("transcript_sequence", genome.transcript_sequence, transcript_id)
        sequences[transcript_id] = None if value is None else str(value)
    return {
        "status": "completed",
        "operation": "fetch_transcript_sequence",
        "sequences": sequences,
        "lengths": {k: (None if v is None else len(v)) for k, v in sequences.items()},
        "missing": [k for k, v in sequences.items() if v is None],
        "provenance": _provenance("Genome.transcript_sequence"),
    }


def query_annotation_at_locus(
    gtf: str | Path,
    contig: str,
    start: int,
    end: int | None = None,
    *,
    feature: str = "transcript",
    strand: str | None = None,
    output: str = "records",
    cache_dir: str | Path | None = None,
    max_results: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the transcripts, exons or protein IDs that overlap a locus (``Genome.transcripts_at_locus``,
    ``exons_at_locus``, ``protein_ids_at_locus`` and their ``*_ids_at_locus`` / ``*_names_at_locus`` forms).

    ``feature`` is ``transcript``, ``exon`` or ``protein``; ``output`` is ``records`` (objects as
    dictionaries), ``ids``, or ``names`` (transcripts only). ``end`` defaults to ``start``. For genes use the
    existing ``query_genes_at_locus``.
    """

    kind = _choice(feature, "feature", ("transcript", "exon", "protein"))
    form = _choice(output, "output", ("records", "ids", "names"))
    if form == "names" and kind != "transcript":
        raise ConfigurationError("output='names' is only available for transcripts.")
    if form == "records" and kind == "protein":
        raise ConfigurationError("protein has no record form; use output='ids'.")
    chrom = _text(contig, "contig")
    first = _int(start, "start", 1, 2**31 - 1)
    last = first if end is None else _int(end, "end", first, 2**31 - 1)
    side = _strand(strand)
    limit = _limit(max_results)
    genome = _genome(gtf, cache_dir=cache_dir)
    if chrom not in genome.contigs():
        raise ConfigurationError(f"contig {chrom!r} is not in the annotation.")
    method = {
        ("transcript", "records"): "transcripts_at_locus",
        ("transcript", "ids"): "transcript_ids_at_locus",
        ("transcript", "names"): "transcript_names_at_locus",
        ("exon", "records"): "exons_at_locus",
        ("exon", "ids"): "exon_ids_at_locus",
        ("protein", "ids"): "protein_ids_at_locus",
    }[(kind, form)]
    found = _call(method, getattr(genome, method), chrom, first, last, strand=side)
    rows, truncated = _cut(list(found), limit)
    if form == "records":
        rows = [_record(kind, item) for item in rows]
    return {
        "status": "completed",
        "operation": "query_annotation_at_locus",
        "feature": kind,
        "output": form,
        "locus": {"contig": chrom, "start": first, "end": last, "strand": side},
        "count": len(found),
        "truncated": truncated,
        "results": rows,
        "provenance": _provenance(f"Genome.{method}"),
    }


def merge_gene_intervals(
    gtf: str | Path,
    contig: str,
    *,
    strand: str | None = None,
    cache_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Union of all gene loci on a contig as sorted non-overlapping intervals (``Genome.merged_gene_intervals``).

    Adjacent intervals (``end + 1 == next start``) are merged; ``strand`` restricts to one strand.
    """

    chrom = _text(contig, "contig")
    side = _strand(strand)
    genome = _genome(gtf, cache_dir=cache_dir)
    if chrom not in genome.contigs():
        raise ConfigurationError(f"contig {chrom!r} is not in the annotation.")
    merged = _call("merged_gene_intervals", genome.merged_gene_intervals, chrom, strand=side)
    intervals = [[int(a), int(b)] for a, b in merged]
    return {
        "status": "completed",
        "operation": "merge_gene_intervals",
        "contig": chrom,
        "strand": side,
        "interval_count": len(intervals),
        "covered_bases": sum(b - a + 1 for a, b in intervals),
        "intervals": intervals,
        "provenance": _provenance("Genome.merged_gene_intervals"),
    }


def find_nearest_transcript(
    gtf: str | Path,
    contig: str,
    position: int,
    *,
    end: int | None = None,
    strand: str | None = None,
    cache_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the transcript whose locus is nearest to a position or interval, even if none overlaps
    (``Genome.nearest_transcript``).

    ``distance`` is 0 when the transcript overlaps the position (``position`` to ``end``).
    """

    chrom = _text(contig, "contig")
    first = _int(position, "position", 1, 2**31 - 1)
    last = None if end is None else _int(end, "end", first, 2**31 - 1)
    side = _strand(strand)
    genome = _genome(gtf, cache_dir=cache_dir)
    if chrom not in genome.contigs():
        raise ConfigurationError(f"contig {chrom!r} is not in the annotation.")
    found = _call(
        "nearest_transcript", genome.nearest_transcript, chrom, first, end=last, strand=side
    )
    distance, transcript = found if found else (None, None)
    return {
        "status": "completed",
        "operation": "find_nearest_transcript",
        "query": {"contig": chrom, "position": first, "end": last, "strand": side},
        "distance": None if distance is None else int(distance),
        "transcript": None if transcript is None else _record("transcript", transcript),
        "provenance": _provenance("Genome.nearest_transcript"),
    }


def get_annotation_by_id(
    gtf: str | Path,
    feature: str,
    identifier: str,
    *,
    cache_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Get the annotation of one gene, transcript or exon by ID (``Genome.gene_by_id``, ``transcript_by_id``,
    ``exon_by_id``; the name and locus look-ups ``gene_name_of_gene_id``, ``locus_of_*_id`` and
    ``transcript_name_of_transcript_id`` are included in the record).

    ``feature`` is ``gene``, ``transcript`` or ``exon``.
    """

    kind = _choice(feature, "feature", _TO_KINDS)
    ident = _text(identifier, "identifier")
    genome = _genome(gtf, cache_dir=cache_dir)
    getter = {
        "gene": genome.gene_by_id,
        "transcript": genome.transcript_by_id,
        "exon": genome.exon_by_id,
    }[kind]
    item = _call(f"{kind}_by_id", getter, ident)
    record = _record(kind, item)
    if kind == "transcript":
        record["exon_ids"] = list(
            _call("exon_ids_of_transcript_id", genome.exon_ids_of_transcript_id, ident)
        )
    elif kind == "gene":
        record["transcript_ids"] = list(
            _call("transcript_ids_of_gene_id", genome.transcript_ids_of_gene_id, ident)
        )
    return {
        "status": "completed",
        "operation": "get_annotation_by_id",
        "feature": kind,
        "record": record,
        "provenance": _provenance(f"Genome.{kind}_by_id"),
    }


def list_annotation_contigs(
    gtf: str | Path, *, cache_dir: str | Path | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List the contig names that carry genes in an annotation (``Genome.contigs``)."""

    genome = _genome(gtf, cache_dir=cache_dir)
    contigs = [str(c) for c in _call("contigs", genome.contigs)]
    return {
        "status": "completed",
        "operation": "list_annotation_contigs",
        "contig_count": len(contigs),
        "contigs": contigs,
        "provenance": _provenance("Genome.contigs"),
    }


def list_annotation_features(
    gtf: str | Path,
    feature: str = "gene",
    *,
    contig: str | None = None,
    strand: str | None = None,
    biotype: str | None = None,
    output: str = "ids",
    cache_dir: str | Path | None = None,
    max_results: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List genes, transcripts, exons or protein IDs, filtered by contig, strand or biotype
    (``Genome.gene_ids`` / ``gene_names`` / ``genes``, ``transcript_ids`` / ``transcript_names`` /
    ``transcripts``, ``exon_ids`` / ``exons``, ``protein_ids``).

    ``feature`` is ``gene``, ``transcript``, ``exon`` or ``protein``; ``output`` is ``ids``, ``names``
    (genes and transcripts) or ``records`` (genes, transcripts, exons). ``biotype`` applies to genes and
    transcripts only (pyensembl has no exon or protein biotype).
    """

    kind = _choice(feature, "feature", _FEATURES)
    form = _choice(output, "output", ("ids", "names", "records"))
    if form == "names" and kind not in ("gene", "transcript"):
        raise ConfigurationError("output='names' is only available for genes and transcripts.")
    if form == "records" and kind == "protein":
        raise ConfigurationError("protein has no record form; use output='ids'.")
    if biotype is not None and kind not in ("gene", "transcript"):
        raise ConfigurationError("biotype filters genes and transcripts only.")
    chrom = None if contig is None else _text(contig, "contig")
    side = _strand(strand)
    limit = _limit(max_results)
    genome = _genome(gtf, cache_dir=cache_dir)
    filters: dict[str, Any] = {"contig": chrom, "strand": side}
    if biotype is not None:
        filters["biotype"] = _text(biotype, "biotype")
    method = {
        ("gene", "ids"): "gene_ids",
        ("gene", "names"): "gene_names",
        ("gene", "records"): "genes",
        ("transcript", "ids"): "transcript_ids",
        ("transcript", "names"): "transcript_names",
        ("transcript", "records"): "transcripts",
        ("exon", "ids"): "exon_ids",
        ("exon", "records"): "exons",
        ("protein", "ids"): "protein_ids",
    }[(kind, form)]
    if method in ("gene_names", "transcript_names", "exon_ids", "exons", "protein_ids"):
        filters.pop("biotype", None)
        if biotype is not None:
            raise ConfigurationError(
                f"biotype is not supported together with output={form!r} for {kind}."
            )
    found = _call(method, getattr(genome, method), **filters)
    rows, truncated = _cut(list(found), limit)
    if form == "records":
        rows = [_record(kind, item) for item in rows]
    return {
        "status": "completed",
        "operation": "list_annotation_features",
        "feature": kind,
        "output": form,
        "filters": {"contig": chrom, "strand": side, "biotype": biotype},
        "count": len(found),
        "truncated": truncated,
        "results": rows,
        "provenance": _provenance(f"Genome.{method}"),
    }


def related_annotation_features(
    gtf: str | Path,
    from_kind: str,
    identifier: str,
    to_kind: str,
    *,
    cache_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Follow the relations between genes, transcripts, exons and proteins in an annotation.

    ``from_kind`` is ``gene_id``, ``gene_name``, ``transcript_id``, ``transcript_name``, ``exon_id`` or
    ``protein_id``; ``to_kind`` is ``gene``, ``transcript`` or ``exon``. Supported pairs: protein to gene
    (``gene_by_protein_id``) and to transcript (``transcript_by_protein_id``); transcript to gene
    (``gene_name_of_transcript_id`` / ``_name``) and to exons (``exon_ids_of_transcript_id`` / ``_name``);
    transcript name to transcripts (``transcripts_by_name``); exon to gene (``gene_name_of_exon_id``) and to
    transcripts (``transcript_ids_of_exon_id``); gene to transcripts (``transcript_ids_of_gene_id`` /
    ``_gene_name``) and to exons (``exon_ids_of_gene_id`` / ``_gene_name``).
    """

    source = _choice(from_kind, "from_kind", _FROM_KINDS)
    target = _choice(to_kind, "to_kind", _TO_KINDS)
    ident = _text(identifier, "identifier")
    genome = _genome(gtf, cache_dir=cache_dir)
    kind_of_source = source.split("_")[0]
    if kind_of_source == target and source.endswith("_id"):
        raise ConfigurationError(
            "from_kind and to_kind name the same feature; use get_annotation_by_id."
        )
    plan = {
        ("protein_id", "gene"): ("gene_by_protein_id", "gene"),
        ("protein_id", "transcript"): ("transcript_by_protein_id", "transcript"),
        ("transcript_id", "gene"): ("gene_name_of_transcript_id", "gene_from_name"),
        ("transcript_name", "gene"): ("gene_name_of_transcript_name", "gene_from_name"),
        ("transcript_id", "exon"): ("exon_ids_of_transcript_id", "exon_ids"),
        ("transcript_name", "exon"): ("exon_ids_of_transcript_name", "exon_ids"),
        ("transcript_name", "transcript"): ("transcripts_by_name", "transcripts"),
        ("exon_id", "gene"): ("gene_name_of_exon_id", "gene_from_name"),
        ("exon_id", "transcript"): ("transcript_ids_of_exon_id", "transcript_ids"),
        ("gene_id", "transcript"): ("transcript_ids_of_gene_id", "transcript_ids"),
        ("gene_name", "transcript"): ("transcript_ids_of_gene_name", "transcript_ids"),
        ("gene_id", "exon"): ("exon_ids_of_gene_id", "exon_ids"),
        ("gene_name", "exon"): ("exon_ids_of_gene_name", "exon_ids"),
    }
    if (source, target) not in plan:
        raise ConfigurationError(f"No relation from {source} to {target} is available.")
    method, shape = plan[(source, target)]
    raw = _call(method, getattr(genome, method), ident)
    results: list[dict[str, Any]]
    if shape == "gene":
        results = [_record("gene", raw)]
    elif shape == "transcript":
        results = [_record("transcript", raw)]
    elif shape == "transcripts":
        results = [_record("transcript", t) for t in raw]
    elif shape == "gene_from_name":
        names = [raw] if isinstance(raw, str) else list(raw)
        results = []
        for name in names:
            results.extend(
                _record("gene", g) for g in _call("genes_by_name", genome.genes_by_name, name)
            )
    elif shape == "transcript_ids":
        results = [
            _record("transcript", _call("transcript_by_id", genome.transcript_by_id, t))
            for t in raw
        ]
    else:
        results = [_record("exon", _call("exon_by_id", genome.exon_by_id, e)) for e in raw]
    return {
        "status": "completed",
        "operation": "related_annotation_features",
        "from_kind": source,
        "identifier": ident,
        "to_kind": target,
        "count": len(results),
        "results": results,
        "provenance": _provenance(f"Genome.{method}"),
    }


def map_position_to_transcript_offset(
    gtf: str | Path,
    transcript_id: str,
    position: int,
    *,
    cache_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a genomic position to the 0-based offset in the spliced transcript sequence
    (``Transcript.spliced_offset``).

    The offset counts exon bases only, from the transcript's 5' end (strand-aware). A position that lies
    in an intron or outside every exon is reported with ``offset: null`` and the reason.
    """

    genome = _genome(gtf, cache_dir=cache_dir)
    ident = _text(transcript_id, "transcript_id")
    pos = _int(position, "position", 1, 2**31 - 1)
    transcript = _call("transcript_by_id", genome.transcript_by_id, ident)
    try:
        offset: int | None = int(transcript.spliced_offset(pos))
        reason = None
    except ValueError as exc:
        offset, reason = None, str(exc)
    return {
        "status": "completed",
        "operation": "map_position_to_transcript_offset",
        "transcript": _record("transcript", transcript),
        "position": pos,
        "offset": offset,
        "reason": reason,
        "provenance": _provenance("Transcript.spliced_offset"),
    }


def genomic_locus_offset(
    contig: str,
    start: int,
    end: int,
    strand: str,
    *,
    position: int | None = None,
    range_start: int | None = None,
    range_end: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Offset of a position (or a position range) from the stranded start of a locus
    (``Locus.offset``, ``Locus.offset_range``).

    On the plus strand the offset of position 13 in 10..20 is 3, on the minus strand 7. Give ``position``
    or both ``range_start`` and ``range_end``; positions must lie inside the locus.
    """

    pyensembl = require("pyensembl", _EXTRA)
    first = _int(start, "start", 1, 2**31 - 1)
    last = _int(end, "end", first, 2**31 - 1)
    locus = _call(
        "Locus",
        pyensembl.Locus,
        _text(contig, "contig"),
        first,
        last,
        _choice(strand, "strand", ("+", "-")),
    )
    if position is not None and (range_start is not None or range_end is not None):
        raise ConfigurationError("Give position or range_start/range_end, not both.")
    if position is not None:
        offset: Any = _call("offset", locus.offset, _int(position, "position", 1, 2**31 - 1))
        method = "Locus.offset"
    elif range_start is not None and range_end is not None:
        offset = list(
            _call(
                "offset_range",
                locus.offset_range,
                _int(range_start, "range_start", 1, 2**31 - 1),
                _int(range_end, "range_end", 1, 2**31 - 1),
            )
        )
        method = "Locus.offset_range"
    else:
        raise ConfigurationError("Give position, or both range_start and range_end.")
    return {
        "status": "completed",
        "operation": "genomic_locus_offset",
        "locus": {"contig": locus.contig, "start": first, "end": last, "strand": strand},
        "offset": offset,
        "provenance": _provenance(method),
    }


def compare_genomic_loci(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    *,
    ignore_strand: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distance and overlap length of two loci (``Locus.distance_to_locus`` / ``distance_to_interval`` and
    ``Locus.overlap_length``).

    Each locus is ``{"contig", "start", "end", "strand"}`` (1-based inclusive). ``distance`` is 0 for
    overlapping loci and ``null`` for loci on different contigs (or strands, unless ``ignore_strand``);
    ``overlap_length`` is the number of shared bases (0 when they cannot overlap).
    """

    pyensembl = require("pyensembl", _EXTRA)
    flag = ignore_strand
    if not isinstance(flag, bool):
        raise ConfigurationError("ignore_strand must be a boolean.")

    def build(item: Mapping[str, Any], name: str) -> Any:
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"{name} must be a mapping with contig, start, end, strand.")
        begin = _int(item.get("start"), f"{name}.start", 1, 2**31 - 1)
        stop = _int(item.get("end"), f"{name}.end", begin, 2**31 - 1)
        return _call(
            "Locus",
            pyensembl.Locus,
            _text(item.get("contig"), f"{name}.contig"),
            begin,
            stop,
            _choice(item.get("strand", "+"), f"{name}.strand", ("+", "-")),
        )

    one, two = build(first, "first"), build(second, "second")
    distance = _call("distance_to_locus", one.distance_to_locus, two, ignore_strand=flag)
    overlap = _call("overlap_length", one.overlap_length, two, ignore_strand=flag)
    return {
        "status": "completed",
        "operation": "compare_genomic_loci",
        "first": _locus(one),
        "second": _locus(two),
        "ignore_strand": flag,
        "distance": None if distance is None or distance == float("inf") else int(distance),
        "overlap_length": int(overlap),
        "provenance": _provenance("Locus.distance_to_locus / Locus.overlap_length"),
    }


__all__ = [
    "compare_genomic_loci",
    "fetch_genomic_sequence",
    "fetch_transcript_sequence",
    "find_nearest_transcript",
    "genomic_locus_offset",
    "get_annotation_by_id",
    "list_annotation_contigs",
    "list_annotation_features",
    "load_hgnc_gene_aliases",
    "map_position_to_transcript_offset",
    "merge_gene_intervals",
    "query_annotation_at_locus",
    "related_annotation_features",
]
_register_public_verbose(globals())
