"""Biopython alignment, search-result and alignment-matrix functions of integration pack P18.

The functions cover the formats of ``Bio.Align`` (reading and writing A2M, bigBed, bigMaf,
bigPsl, chain, Clustal, EMBOSS, Exonerate, MAF, Mauve, MSF, NEXUS, PHYLIP, PSL, Stockholm and
tabular alignments), MAF index queries, ``Bio.SearchIO`` readers for HMMER and Infernal
output (with indexing, filtering and sorting), the FOGSAA and codon-aware pairwise aligners,
the analysis methods of ``Bio.Align.Alignment`` (sorting, substitution counts, column
frequencies, identity counts, coordinate maps, merging and mapping) and the matrix operations
of ``Bio.Nexus`` (constant and gap-only columns, site states, subsetting, gap insertion).

Alignments are passed as mappings of sequence name to gapped sequence (``-`` for a gap) and
returned the same way, so every result is JSON-serialisable. Biopython is imported only when a
function runs. File writers never overwrite; the output path must be new.
"""

from __future__ import annotations

import contextlib
import re
import shutil
import tempfile
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _biopython_p18_common as _c
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = _c.EXTRA
_READ_FORMATS = (
    "a2m",
    "bigmaf",
    "bigpsl",
    "chain",
    "emboss",
    "exonerate",
    "mauve",
    "msf",
    "nexus",
    "stockholm",
    "tabular",
)
_ROW_WRITE_FORMATS = ("a2m", "clustal", "fasta", "maf", "mauve", "nexus", "phylip", "stockholm")
_CONVERT_SOURCES = (
    "a2m",
    "bigbed",
    "bigmaf",
    "bigpsl",
    "chain",
    "clustal",
    "emboss",
    "exonerate",
    "fasta",
    "maf",
    "mauve",
    "msf",
    "nexus",
    "phylip",
    "psl",
    "stockholm",
    "tabular",
)
_CONVERT_TARGETS = (
    "a2m",
    "bigbed",
    "bigmaf",
    "bigpsl",
    "chain",
    "clustal",
    "exonerate",
    "maf",
    "mauve",
    "phylip",
    "psl",
    "stockholm",
)
_SEQUENCE_TARGETS = frozenset({"a2m", "clustal", "maf", "mauve", "phylip", "stockholm"})
_BINARY_TARGETS = frozenset({"bigbed", "bigmaf", "bigpsl"})
_SEARCH_FORMATS = (
    "hmmer2-text",
    "hmmer3-text",
    "hmmer3-tab",
    "hmmscan3-domtab",
    "hmmsearch3-domtab",
    "phmmer3-domtab",
    "infernal-tab",
    "infernal-text",
    "blast-tab",
    "blast-xml",
    "blat-psl",
    "exonerate-text",
    "exonerate-cigar",
    "exonerate-vulgar",
    "fasta-m10",
)
_UNINDEXABLE = ("hmmer3-tab", "hmmscan3-domtab", "hmmsearch3-domtab", "phmmer3-domtab")
_SORT_KEYS = ("evalue", "bitscore", "id", "length", "hsp_count")
_ROW_SORT_KEYS = ("id", "sequence", "gc_fraction", "ungapped_length", "gap_count")
_NAME = re.compile(r"^[A-Za-z0-9_.|:-]+$")
_ROW_SYMBOLS = re.compile(r"^[A-Za-z*.\-?]+$")
_NEXUS_SYMBOLS = frozenset("ACGTURYSWKMBDHVN?-")
_MAX_ROWS = 10_000
_MAX_COLUMNS = 1_000_000
_MAX_ALIGNMENTS = 10_000


# ---------------------------------------------------------------- shared helpers


def _rows(
    sequences: object, name: str = "aligned_sequences", *, minimum: int = 1, nexus: bool = False
) -> dict[str, str]:
    """Validate a mapping of sequence name to gapped sequence (equal lengths)."""

    if not isinstance(sequences, Mapping) or len(sequences) < minimum:
        raise ConfigurationError(f"{name} must be a mapping with at least {minimum} sequence(s).")
    if len(sequences) > _MAX_ROWS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_ROWS} sequences.")
    cleaned: dict[str, str] = {}
    for key, value in sequences.items():
        if not isinstance(key, str) or not _NAME.match(key):
            raise ConfigurationError(f"{name} keys must be names without spaces; got {key!r}.")
        if not isinstance(value, str) or not value or not _ROW_SYMBOLS.match(value):
            raise ConfigurationError(f"{name}[{key!r}] must be a gapped sequence of letters.")
        if len(value) > _MAX_COLUMNS:
            raise ConfigurationError(f"{name}[{key!r}] is longer than {_MAX_COLUMNS}.")
        cleaned[key] = value.upper() if nexus else value
    if len({len(value) for value in cleaned.values()}) != 1:
        raise ConfigurationError(f"{name} must all have the same aligned length.")
    if nexus:
        for key, value in cleaned.items():
            if set(value) - _NEXUS_SYMBOLS:
                raise ConfigurationError(f"{name}[{key!r}] has symbols outside DNA, '?' and '-'.")
    return cleaned


def _build_alignment(rows: Mapping[str, str], *, molecule_type: str | None = None) -> Any:
    """Return a ``Bio.Align.Alignment`` for gapped rows (sequence records carry the names)."""

    align = require("Bio.Align", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    try:
        sequences, coordinates = align.Alignment.parse_printed_alignment(
            [value.encode("ascii") for value in rows.values()]
        )
    except Exception as exc:
        raise ConfigurationError(f"The rows are not a valid printed alignment: {exc}") from exc
    records = []
    for name, sequence in zip(rows, sequences, strict=True):
        record = record_module.SeqRecord(
            seq_module.Seq(sequence.decode("ascii")), id=name, name=name, description=""
        )
        if molecule_type:
            record.annotations["molecule_type"] = molecule_type
        records.append(record)
    return align.Alignment(records, coordinates)


def _gapped_rows(alignment: Any) -> dict[str, str] | None:
    try:
        return {
            str(sequence.id): str(alignment[index])
            for index, sequence in enumerate(alignment.sequences)
        }
    except Exception:
        return None


def _sequence_length(sequence: object) -> int | None:
    try:
        length = len(sequence)  # type: ignore[arg-type]
    except Exception:
        return None
    return length if 0 < length < 2**40 else None


def _alignment_summary(item: Any, index: int, names: Sequence[str] | None = None) -> dict[str, Any]:
    ids = [str(getattr(sequence, "id", "")) for sequence in item.sequences]
    labels = list(names) if names is not None else ids
    rows: dict[str, str] | None
    try:
        if len(labels) == len(set(labels)) and all(labels):
            rows = {label: str(item[i]) for i, label in enumerate(labels)}
        else:
            rows = {f"{label}#{i}": str(item[i]) for i, label in enumerate(labels)}
    except Exception:
        rows = None
    try:
        columns: int | None = int(item.length)
    except Exception:  # alignments with unequal step sizes (codon layout) have no column count
        columns = None
    return {
        "index": index,
        "sequence_ids": ids,
        "sequence_lengths": [_sequence_length(sequence) for sequence in item.sequences],
        "aligned_columns": columns,
        "score": _c.optional_float(getattr(item, "score", None)),
        "coordinates": [[int(v) for v in row] for row in item.coordinates.tolist()],
        "aligned_rows": rows,
        "annotations": _c.plain(getattr(item, "annotations", {})),
        "column_annotations": _c.plain(getattr(item, "column_annotations", {})),
    }


def _ignore_warnings() -> Any:
    manager = warnings.catch_warnings()
    manager.__enter__()
    warnings.simplefilter("ignore")
    return manager


# ---------------------------------------------------------------- alignment files


def read_alignment_records(
    path: str | Path,
    *,
    file_format: str = "nexus",
    max_alignments: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read alignments of a format that the other alignment readers do not cover.

    ``file_format`` is one of ``a2m``, ``bigmaf`` (bigMaf, a binary bigBed file), ``bigpsl``
    (bigPsl), ``chain`` (UCSC chain), ``emboss`` (EMBOSS needle/water/matcher output),
    ``exonerate`` (text, cigar and vulgar output are detected), ``mauve`` (XMFA), ``msf``,
    ``nexus``, ``stockholm`` and ``tabular`` (BLAST or FASTA tabular output with ``-outfmt 7``
    comments). Each alignment has its sequence ids and lengths, the number of columns, the
    score when the format stores one, the block coordinates (0-based, half-open, one row per
    sequence), the gapped rows when the sequences are stored in the file, and the annotations.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _READ_FORMATS, "file_format")
    limit = _c.integer(max_alignments, "max_alignments", minimum=1, maximum=_MAX_ALIGNMENTS)
    align = require("Bio.Align", _EXTRA)
    found: list[dict[str, Any]] = []
    guard = _ignore_warnings()
    try:
        for item in align.parse(str(source), fmt):
            if len(found) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_alignments={limit}.")
            found.append(_alignment_summary(item, len(found)))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(
            f"{source.name} is not a valid {fmt} alignment file: {exc}"
        ) from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "read_alignment_records",
        "source": str(source),
        "file_format": fmt,
        "alignment_count": len(found),
        "alignments": found,
        "provenance": _c.provenance(
            parser="Bio.Align.parse",
            coordinates="block boundaries, 0-based and half-open, one row per sequence",
        ),
    }


def write_aligned_sequences(
    aligned_sequences: Mapping[str, str],
    output_path: str | Path,
    *,
    file_format: str = "clustal",
    molecule_type: str = "DNA",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one multiple alignment (name to gapped sequence) to a new file.

    ``file_format`` is ``a2m``, ``clustal``, ``fasta``, ``maf``, ``mauve`` (XMFA), ``nexus``,
    ``phylip`` or ``stockholm``. ``molecule_type`` (``DNA`` or ``RNA``) is used by the NEXUS
    writer. PHYLIP writes names of at most ten characters (Biopython rejects longer names).
    """

    fmt = _c.choice(file_format, _ROW_WRITE_FORMATS, "file_format")
    kind = _c.choice(molecule_type, ("DNA", "RNA"), "molecule_type")
    rows = _rows(
        aligned_sequences, minimum=1 if fmt in ("fasta", "maf") else 2, nexus=fmt == "nexus"
    )
    destination = _c.new_output(output_path)
    align = require("Bio.Align", _EXTRA)
    alignment = _build_alignment(rows, molecule_type=kind)
    options: dict[str, Any] = {}
    if fmt == "a2m":
        alignment.column_annotations = {"state": "D" * alignment.length}
    elif fmt == "stockholm":
        alignment.column_annotations = {}
    elif fmt == "mauve":
        options = {"metadata": {"FormatVersion": "Mauve1"}, "identifiers": list(rows)}
    try:
        with destination.open("x", encoding="utf-8") as handle:
            written = align.write([alignment], handle, fmt, **options)
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(f"Biopython could not write {fmt}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_aligned_sequences",
        "output_path": str(destination),
        "file_format": fmt,
        "sequence_count": len(rows),
        "aligned_columns": len(next(iter(rows.values()))),
        "alignment_count": int(written),
        "provenance": _c.provenance(writer="Bio.Align.write"),
    }


def _prepare_for_target(alignments: Any, target: str) -> Any:
    """Add the attributes that a Biopython writer reads from alignments of another format."""

    align = require("Bio.Align", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    items = list(alignments)
    if target in _BINARY_TARGETS:
        items.sort(key=lambda a: (str(getattr(a.sequences[0], "id", "")), int(a.coordinates[0, 0])))
    prepared = align.Alignments(items)
    for item in items:
        if not hasattr(item, "annotations"):
            item.annotations = {}
        if target in ("stockholm", "a2m") and not hasattr(item, "column_annotations"):
            item.column_annotations = {}
        if target == "a2m":
            item.column_annotations = {"state": "D" * item.length}
    if target in _BINARY_TARGETS:
        if hasattr(alignments, "targets"):
            prepared.targets = alignments.targets
        else:
            sizes: dict[str, int] = {}
            for item in items:
                first = item.sequences[0]
                sizes.setdefault(
                    str(getattr(first, "id", "")),
                    int(getattr(first, "annotations", {}).get("srcSize", len(first))),
                )
            prepared.targets = [
                record_module.SeqRecord(seq_module.Seq(None, length=size), id=name)
                for name, size in sizes.items()
            ]
    if target == "mauve":
        identifiers: list[str] = []
        for item in items:
            for sequence in item.sequences:
                name = str(sequence.id)
                if name not in identifiers:
                    identifiers.append(name)
        prepared.metadata = getattr(alignments, "metadata", None) or {"FormatVersion": "Mauve1"}
        prepared.identifiers = getattr(alignments, "identifiers", None) or identifiers
    return prepared


def convert_alignment_file(
    source: str | Path,
    output_path: str | Path,
    *,
    source_format: str = "psl",
    target_format: str = "chain",
    exonerate_variant: str = "vulgar",
    max_alignments: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an alignment file with ``Bio.Align`` to A2M, bigBed, bigMaf, bigPsl, chain, Clustal,
    Exonerate, MAF, Mauve, PHYLIP, PSL or Stockholm.

    ``source_format`` may be a2m, bigbed, bigmaf, bigpsl, chain, clustal, emboss, exonerate,
    fasta, maf, mauve, msf, nexus, phylip, psl, stockholm or tabular. ``exonerate_variant`` is
    ``vulgar`` or ``cigar`` for the Exonerate target. The text targets that print the sequences
    (a2m, clustal, maf, mauve, phylip, stockholm) need a source that stores sequences; the
    coordinate targets (bigBed, bigMaf, bigPsl, chain, Exonerate, PSL) need a pairwise source
    that Biopython can write in that format, otherwise the function stops with Biopython's
    message. Alignments are sorted by target name and start for the bigBed family. The output
    file must not exist.
    """

    path = _c.source_file(source, "source")
    in_format = _c.choice(source_format, _CONVERT_SOURCES, "source_format")
    out_format = _c.choice(target_format, _CONVERT_TARGETS, "target_format")
    variant = _c.choice(exonerate_variant, ("vulgar", "cigar"), "exonerate_variant")
    limit = _c.integer(max_alignments, "max_alignments", minimum=1, maximum=_MAX_ALIGNMENTS)
    destination = _c.new_output(output_path)
    align = require("Bio.Align", _EXTRA)
    guard = _ignore_warnings()
    try:
        parsed = align.parse(str(path), in_format)
        items: list[Any] = []
        for item in parsed:
            if len(items) >= limit:
                raise ConfigurationError(f"{path.name} exceeds max_alignments={limit}.")
            items.append(item)
        prepared = _prepare_for_target(_Carrier(items, parsed), out_format)
        binary = out_format in _BINARY_TARGETS
        with _c.open_new(destination, binary=binary) as handle:
            if out_format == "exonerate":
                # Align.write has its own ``fmt`` argument, so the writer class is used directly.
                writer_module = require("Bio.Align.exonerate", _EXTRA)
                written = writer_module.AlignmentWriter(handle, fmt=variant).write(prepared)
            else:
                written = align.write(prepared, handle, out_format)
    except ConfigurationError:
        destination.unlink(missing_ok=True)
        raise
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(
            f"Biopython could not convert {path.name} from {in_format} to {out_format}: "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "convert_alignment_file",
        "source": str(path),
        "output_path": str(destination),
        "source_format": in_format,
        "target_format": out_format,
        "alignment_count": len(items) if written is None else int(written),
        "bytes_written": destination.stat().st_size,
        "provenance": _c.provenance(reader="Bio.Align.parse", writer="Bio.Align.write"),
    }


class _Carrier(list[Any]):
    """List of alignments that keeps the attributes (targets, metadata) of the reader."""

    def __init__(self, items: list[Any], reader: Any) -> None:
        super().__init__(items)
        for name in ("targets", "metadata", "identifiers"):
            if hasattr(reader, name):
                setattr(self, name, getattr(reader, name))


# ---------------------------------------------------------------- MAF index


def _maf_index(
    maf_path: object, target: object, index_path: object
) -> tuple[Any, Path, Path | None]:
    source = _c.source_file(maf_path, "maf_path")
    name = _c.text(target, "target_seqname", limit=500).strip()
    maf_io = require("Bio.AlignIO.MafIO", _EXTRA)
    temporary: Path | None = None
    if index_path is None:
        temporary = Path(tempfile.mkdtemp(prefix="dnakit_mafindex_"))
        index_file = temporary / "maf.index"
    else:
        if not isinstance(index_path, (str, Path)):
            raise ConfigurationError("index_path must be a file path or None.")
        index_file = Path(index_path).expanduser().resolve()
        if not index_file.parent.is_dir():
            raise ConfigurationError(
                f"index_path parent directory does not exist: {index_file.parent}"
            )
    try:
        index = maf_io.MafIndex(str(index_file), str(source), name)
    except Exception as exc:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)
        raise ConfigurationError(f"Biopython could not index {source.name}: {exc}") from exc
    return index, source, temporary


def _intervals(starts: object, ends: object) -> tuple[list[int], list[int]]:
    for label, value in (("starts", starts), ("ends", ends)):
        if not isinstance(value, Sequence) or isinstance(value, str) or not value:
            raise ConfigurationError(f"{label} must be a non-empty list of integers.")
    assert isinstance(starts, Sequence) and isinstance(ends, Sequence)
    if len(starts) != len(ends) or len(starts) > 10_000:
        raise ConfigurationError("starts and ends must have the same length (at most 10000).")
    first = [_c.integer(v, "starts", minimum=0, maximum=2**40) for v in starts]
    last = [_c.integer(v, "ends", minimum=0, maximum=2**40) for v in ends]
    if any(e <= s for s, e in zip(first, last, strict=True)):
        raise ConfigurationError("every end must be greater than its start.")
    return first, last


def _block_rows(block: Any) -> list[dict[str, Any]]:
    return [
        {
            "id": str(record.id),
            "sequence": str(record.seq),
            "annotations": _c.plain(record.annotations),
        }
        for record in block
    ]


def search_maf_index(
    maf_path: str | Path,
    target_seqname: str,
    starts: Sequence[int],
    ends: Sequence[int],
    *,
    index_path: str | Path | None = None,
    max_blocks: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the alignment blocks that overlap regions of the reference sequence of a MAF file.

    A ``Bio.AlignIO.MafIO.MafIndex`` is built over ``target_seqname`` (for example
    ``mm9.chr10``), kept in the SQLite file ``index_path`` and reused when it exists; without
    ``index_path`` a temporary index is used and removed. ``starts`` are 0-based, ``ends`` are
    half-open, as in the UCSC convention. Each returned block lists its rows (id, gapped
    sequence, MAF annotations such as start, size, strand and srcSize).
    """

    first, last = _intervals(starts, ends)
    limit = _c.integer(max_blocks, "max_blocks", minimum=1, maximum=_MAX_ALIGNMENTS)
    index, source, temporary = _maf_index(maf_path, target_seqname, index_path)
    blocks: list[dict[str, Any]] = []
    try:
        for block in index.search(first, last):
            if len(blocks) >= limit:
                raise ConfigurationError(f"The search returns more than max_blocks={limit} blocks.")
            blocks.append(
                {
                    "index": len(blocks),
                    "row_count": len(block),
                    "alignment_length": int(block.get_alignment_length()),
                    "rows": _block_rows(block),
                }
            )
        indexed = len(index)
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"MAF search failed: {exc}") from exc
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)
    return {
        "status": "completed",
        "operation": "search_maf_index",
        "source": str(source),
        "target_seqname": target_seqname,
        "indexed_blocks": int(indexed),
        "block_count": len(blocks),
        "blocks": blocks,
        "provenance": _c.provenance(method="Bio.AlignIO.MafIO.MafIndex.search"),
    }


def splice_maf_index(
    maf_path: str | Path,
    target_seqname: str,
    starts: Sequence[int],
    ends: Sequence[int],
    *,
    strand: int = 1,
    index_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the in-silico spliced alignment of several reference segments (for example exons).

    The segments ``starts``/``ends`` (0-based, half-open) on ``target_seqname`` are cut out of
    the MAF alignment and joined, giving one multiple alignment with a row per species
    (``MafIndex.get_spliced``). ``strand`` is 1 or -1; -1 reverse-complements the result.
    """

    first, last = _intervals(starts, ends)
    if isinstance(strand, bool) or strand not in (1, -1):
        raise ConfigurationError("strand must be 1 or -1.")
    direction = int(strand)
    index, source, temporary = _maf_index(maf_path, target_seqname, index_path)
    try:
        spliced = index.get_spliced(first, last, direction)
    except Exception as exc:
        raise ConfigurationError(f"MAF splicing failed: {exc}") from exc
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)
    return {
        "status": "completed",
        "operation": "splice_maf_index",
        "source": str(source),
        "target_seqname": target_seqname,
        "strand": direction,
        "segment_count": len(first),
        "alignment_length": int(spliced.get_alignment_length()),
        "rows": {str(r.id): str(r.seq) for r in sorted(spliced, key=lambda r: str(r.id))},
        "provenance": _c.provenance(method="Bio.AlignIO.MafIO.MafIndex.get_spliced"),
    }


# ---------------------------------------------------------------- search results


_HSP_POSITION_NAMES = (
    "query_start",
    "query_end",
    "hit_start",
    "hit_end",
    "query_strand",
    "hit_strand",
    "query_frame",
    "hit_frame",
    "aln_span",
)


def _public_attributes(obj: object, skip: tuple[str, ...] = ()) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in vars(obj).items():
        if key.startswith("_") or key in skip or key in ("fragments",):
            continue
        result[key] = _c.plain(value)
    return result


def _hsp_dict(hsp: Any) -> dict[str, Any]:
    data = _public_attributes(hsp, skip=("aln_annotation",))
    for name in _HSP_POSITION_NAMES:
        try:
            value = getattr(hsp, name)
        except Exception:
            continue
        if value is not None:
            data[name] = _c.plain(value)
    for name in ("query_range", "hit_range"):
        with contextlib.suppress(Exception):
            data[name] = [int(v) for v in getattr(hsp, name)]
    data["fragment_count"] = len(hsp.fragments)
    return data


def _hit_dict(hit: Any) -> dict[str, Any]:
    data = _public_attributes(hit)
    data.update({"id": str(hit.id), "description": str(hit.description)})
    data["hsp_count"] = len(hit.hsps)
    data["hsps"] = [_hsp_dict(hsp) for hsp in hit.hsps]
    return data


def _query_dict(query: Any) -> dict[str, Any]:
    data = _public_attributes(query)
    data.update({"id": str(query.id), "description": str(query.description)})
    data["hit_count"] = len(query)
    data["hits"] = [_hit_dict(hit) for hit in query.hits]
    return data


def _search_parse_kwargs(fmt: str, comments: bool) -> dict[str, Any]:
    return {"comments": True} if fmt == "blast-tab" and comments else {}


def read_search_results(
    path: str | Path,
    *,
    file_format: str = "hmmer3-text",
    comments: bool = False,
    max_queries: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a search-program report (HMMER, Infernal and other formats) with ``Bio.SearchIO``.

    ``file_format`` is one of hmmer2-text (hmmpfam/hmmsearch 2.x), hmmer3-text, hmmer3-tab
    (``--tblout``), hmmscan3-domtab, hmmsearch3-domtab, phmmer3-domtab (``--domtblout``),
    infernal-tab, infernal-text, blast-tab, blast-xml, blat-psl, exonerate-text,
    exonerate-cigar, exonerate-vulgar and fasta-m10. ``comments=True`` reads BLAST tabular
    files written with ``-outfmt 7``. Every query is returned with its program, target and
    version, its hits (identifier, description, E-value, bit score and the other attributes
    the format provides) and their HSPs with the ranges, scores and strands.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _SEARCH_FORMATS, "file_format")
    use_comments = _c.flag(comments, "comments")
    limit = _c.integer(max_queries, "max_queries", minimum=1, maximum=_c.MAX_RECORDS)
    search_io = require("Bio.SearchIO", _EXTRA)
    queries: list[dict[str, Any]] = []
    guard = _ignore_warnings()
    try:
        for query in search_io.parse(str(source), fmt, **_search_parse_kwargs(fmt, use_comments)):
            if len(queries) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_queries={limit}.")
            queries.append(_query_dict(query))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} report: {exc}") from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "read_search_results",
        "source": str(source),
        "file_format": fmt,
        "query_count": len(queries),
        "hit_count": sum(item["hit_count"] for item in queries),
        "queries": queries,
        "provenance": _c.provenance(
            parser="Bio.SearchIO.parse",
            coordinates="HSP ranges are 0-based and half-open, as Bio.SearchIO stores them",
        ),
    }


def index_search_results(
    path: str | Path,
    *,
    file_format: str = "hmmer3-text",
    index_path: str | Path | None = None,
    query_ids: Sequence[str] | None = None,
    comments: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index a search-output file for random access by query id and return chosen queries.

    Without ``index_path`` an in-memory offset index is built (``SearchIO.index``); with it the
    index is stored in a SQLite file (``SearchIO.index_db``), created when missing and reused
    otherwise. The result lists all query ids; the queries named in ``query_ids`` are read
    through the index and returned in full (see ``read_search_results`` for the fields).
    The HMMER tabular formats (``hmmer3-tab`` and the three ``domtab`` formats) cannot be
    indexed by Biopython 1.88 and are refused with a clear message.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _SEARCH_FORMATS, "file_format")
    if fmt in _UNINDEXABLE:
        raise ConfigurationError(
            f"Bio.SearchIO cannot index {fmt} files (its tabular indexer fails with an IndexError "
            "in Biopython 1.88); read the file with read_search_results instead."
        )
    use_comments = _c.flag(comments, "comments")
    wanted: list[str] = []
    if query_ids is not None:
        if (
            isinstance(query_ids, str)
            or not isinstance(query_ids, Sequence)
            or len(query_ids) > 1000
        ):
            raise ConfigurationError("query_ids must be a list of at most 1000 strings.")
        wanted = [_c.text(item, "query id", limit=1000) for item in query_ids]
    database: Path | None = None
    if index_path is not None:
        if not isinstance(index_path, (str, Path)):
            raise ConfigurationError("index_path must be a file path or None.")
        database = Path(index_path).expanduser().resolve()
        if not database.parent.is_dir():
            raise ConfigurationError(
                f"index_path parent directory does not exist: {database.parent}"
            )
    search_io = require("Bio.SearchIO", _EXTRA)
    kwargs = _search_parse_kwargs(fmt, use_comments)
    guard = _ignore_warnings()
    try:
        if database is None:
            indexed = search_io.index(str(source), fmt, **kwargs)
        else:
            indexed = search_io.index_db(str(database), [str(source)], fmt, **kwargs)
        keys = sorted(str(key) for key in indexed)
        missing = [item for item in wanted if item not in keys]
        if missing:
            raise ConfigurationError(f"query ids not found in the index: {missing}")
        selected = [_query_dict(indexed[item]) for item in wanted]
        indexed.close()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"Biopython could not index {source.name}: {exc}") from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "index_search_results",
        "source": str(source),
        "file_format": fmt,
        "index_path": None if database is None else str(database),
        "query_count": len(keys),
        "query_ids": keys,
        "queries": selected,
        "provenance": _c.provenance(
            method="Bio.SearchIO.index" if database is None else "Bio.SearchIO.index_db"
        ),
    }


def filter_search_results(
    path: str | Path,
    *,
    file_format: str = "hmmer3-text",
    max_evalue: float | None = None,
    min_bitscore: float | None = None,
    min_alignment_span: int | None = None,
    min_hsps_per_hit: int | None = None,
    hit_ids: Sequence[str] | None = None,
    comments: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter the hits and HSPs of a search report (``QueryResult.hsp_filter`` / ``hit_filter``).

    HSPs are kept when their E-value is at most ``max_evalue``, their bit score is at least
    ``min_bitscore`` and their alignment span is at least ``min_alignment_span``; hits that lose
    all HSPs disappear. Hits are then kept when they have at least ``min_hsps_per_hit`` HSPs
    and, if ``hit_ids`` is given, when their id is in that list. Unset criteria are ignored.
    Returns the remaining queries (see ``read_search_results``) and the counts before and after.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _SEARCH_FORMATS, "file_format")
    use_comments = _c.flag(comments, "comments")
    evalue = None if max_evalue is None else _c.number(max_evalue, "max_evalue")
    bitscore = None if min_bitscore is None else _c.number(min_bitscore, "min_bitscore")
    span = _c.optional_integer(min_alignment_span, "min_alignment_span", minimum=0, maximum=2**40)
    hsp_floor = _c.optional_integer(min_hsps_per_hit, "min_hsps_per_hit", minimum=1, maximum=10**6)
    allowed: set[str] | None = None
    if hit_ids is not None:
        if isinstance(hit_ids, str) or not isinstance(hit_ids, Sequence) or len(hit_ids) > 100_000:
            raise ConfigurationError("hit_ids must be a list of at most 100000 strings.")
        allowed = {_c.text(item, "hit id", limit=1000) for item in hit_ids}
    if (
        evalue is None
        and bitscore is None
        and span is None
        and hsp_floor is None
        and allowed is None
    ):
        raise ConfigurationError("Give at least one filter criterion.")
    search_io = require("Bio.SearchIO", _EXTRA)

    def hsp_ok(hsp: Any) -> bool:
        if evalue is not None and not (
            getattr(hsp, "evalue", None) is not None and hsp.evalue <= evalue
        ):
            return False
        if bitscore is not None and not (
            getattr(hsp, "bitscore", None) is not None and hsp.bitscore >= bitscore
        ):
            return False
        return span is None or hsp.aln_span >= span

    def hit_ok(hit: Any) -> bool:
        if hsp_floor is not None and len(hit.hsps) < hsp_floor:
            return False
        return allowed is None or str(hit.id) in allowed

    kept: list[dict[str, Any]] = []
    before = after = 0
    guard = _ignore_warnings()
    try:
        for query in search_io.parse(str(source), fmt, **_search_parse_kwargs(fmt, use_comments)):
            before += len(query)
            result = query
            if evalue is not None or bitscore is not None or span is not None:
                result = result.hsp_filter(hsp_ok)
            if hsp_floor is not None or allowed is not None:
                result = result.hit_filter(hit_ok)
            after += len(result)
            kept.append(_query_dict(result))
    except Exception as exc:
        raise ConfigurationError(f"Could not filter {source.name} as {fmt}: {exc}") from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "filter_search_results",
        "source": str(source),
        "file_format": fmt,
        "hits_before": before,
        "hits_after": after,
        "queries": kept,
        "provenance": _c.provenance(
            method="Bio.SearchIO QueryResult.hsp_filter and QueryResult.hit_filter"
        ),
    }


def _hit_sort_key(name: str) -> Any:
    if name == "evalue":
        return lambda hit: (getattr(hit, "evalue", None) is None, getattr(hit, "evalue", 0.0))
    if name == "bitscore":
        return lambda hit: getattr(hit, "bitscore", float("-inf")) or float("-inf")
    if name == "id":
        return lambda hit: str(hit.id)
    if name == "length":
        return lambda hit: getattr(hit, "seq_len", None) or 0
    return lambda hit: len(hit.hsps)


def sort_search_results(
    path: str | Path,
    *,
    file_format: str = "hmmer3-text",
    sort_by: str = "evalue",
    reverse: bool = False,
    sort_hsps: bool = True,
    comments: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort the hits of every query, and the HSPs of every hit, of a search report.

    ``sort_by`` is ``evalue`` (ascending; hits without an E-value last), ``bitscore``, ``id``,
    ``length`` (hit sequence length when the format stores it) or ``hsp_count``; ``reverse``
    flips the order. HSPs are sorted by the same key when ``sort_hsps`` is true and the key
    applies to HSPs (evalue, bitscore). This is ``QueryResult.sort`` and ``Hit.sort``.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _SEARCH_FORMATS, "file_format")
    key_name = _c.choice(sort_by, _SORT_KEYS, "sort_by")
    flip = _c.flag(reverse, "reverse")
    with_hsps = _c.flag(sort_hsps, "sort_hsps")
    use_comments = _c.flag(comments, "comments")
    search_io = require("Bio.SearchIO", _EXTRA)
    hit_key = _hit_sort_key(key_name)
    queries: list[dict[str, Any]] = []
    guard = _ignore_warnings()
    try:
        for query in search_io.parse(str(source), fmt, **_search_parse_kwargs(fmt, use_comments)):
            query.sort(key=hit_key, reverse=flip, in_place=True)
            if with_hsps and key_name in ("evalue", "bitscore"):
                for hit in query.hits:
                    hit.sort(key=hit_key, reverse=flip, in_place=True)
            queries.append(_query_dict(query))
    except Exception as exc:
        raise ConfigurationError(f"Could not sort {source.name} as {fmt}: {exc}") from exc
    finally:
        guard.__exit__(None, None, None)
    return {
        "status": "completed",
        "operation": "sort_search_results",
        "source": str(source),
        "file_format": fmt,
        "sort_by": key_name,
        "reverse": flip,
        "queries": queries,
        "provenance": _c.provenance(method="Bio.SearchIO QueryResult.sort and Hit.sort"),
    }


# ---------------------------------------------------------------- pairwise aligners


def _alignment_list(
    alignments: Any, limit: int, names: Sequence[str] | None = None
) -> tuple[list[dict[str, Any]], int]:
    found: list[dict[str, Any]] = []
    for item in alignments:
        if len(found) >= limit:
            break
        found.append(_alignment_summary(item, len(found), names))
    try:
        total = len(alignments)
    except OverflowError:
        total = -1
    return found, int(total)


def fogsaa_global_alignment(
    sequence_a: str,
    sequence_b: str,
    *,
    match_score: float = 1.0,
    mismatch_score: float = 0.0,
    open_gap_score: float = 0.0,
    extend_gap_score: float = 0.0,
    max_alignments: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two sequences globally with the FOGSAA algorithm (``PairwiseAligner(mode='fogsaa')``).

    FOGSAA is a branch-and-bound search for the optimal global alignment. The scores default
    to Biopython's (match 1, mismatch 0, no gap penalty); pass explicit scores for a useful
    DNA alignment. Returns the optimal score, the alignment count, and up to ``max_alignments``
    alignments with their gapped rows and block coordinates.
    """

    first = _c.text(sequence_a, "sequence_a", limit=50_000)
    second = _c.text(sequence_b, "sequence_b", limit=50_000)
    scores = {
        "match_score": _c.number(match_score, "match_score"),
        "mismatch_score": _c.number(mismatch_score, "mismatch_score"),
        "open_gap_score": _c.number(open_gap_score, "open_gap_score"),
        "extend_gap_score": _c.number(extend_gap_score, "extend_gap_score"),
    }
    limit = _c.integer(max_alignments, "max_alignments", minimum=1, maximum=1000)
    align = require("Bio.Align", _EXTRA)
    try:
        aligner = align.PairwiseAligner(mode="fogsaa", **scores)
        alignments = aligner.align(first, second)
        best = float(aligner.score(first, second))
        found, total = _alignment_list(alignments, limit, ("sequence_a", "sequence_b"))
    except Exception as exc:
        raise ConfigurationError(f"FOGSAA alignment failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "fogsaa_global_alignment",
        "score": best,
        "alignment_count": total,
        "alignments": found,
        "parameters": scores,
        "provenance": _c.provenance(method="Bio.Align.PairwiseAligner(mode='fogsaa')"),
    }


def align_nucleotide_to_protein(
    protein: str,
    nucleotides: str,
    *,
    codon_table: int = 1,
    anchor_length: int = 10,
    max_alignments: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align a coding nucleotide sequence to its protein with the codon-aware aligner.

    ``Bio.Align.CodonAligner`` pairs each amino acid with its codon and tolerates frame shifts
    and introns of the nucleotide sequence (``anchor_length`` is the seed length in codons
    used to find the next codon after a gap). DNA or RNA is accepted. ``codon_table`` is an
    NCBI genetic-code id. The result has the number of aligned codons (the score), and up to
    ``max_alignments`` alignments with block coordinates (protein in residues, nucleotide
    in bases).
    """

    amino = _c.text(protein, "protein", limit=100_000)
    nucleotide = _c.text(nucleotides, "nucleotides", limit=300_000)
    table_id = _c.integer(codon_table, "codon_table", minimum=1, maximum=33)
    anchor = _c.integer(anchor_length, "anchor_length", minimum=1, maximum=1000)
    limit = _c.integer(max_alignments, "max_alignments", minimum=1, maximum=1000)
    if "U" in nucleotide.upper() and "T" not in nucleotide.upper():
        nucleotide = nucleotide.replace("U", "T").replace("u", "t")  # RNA codons are read as DNA
    align = require("Bio.Align", _EXTRA)
    codon_tables = require("Bio.Data.CodonTable", _EXTRA)
    try:
        table = codon_tables.unambiguous_dna_by_id[table_id]
    except KeyError as exc:
        raise ConfigurationError(
            f"codon_table {table_id} is not an unambiguous DNA table."
        ) from exc
    try:
        aligner = align.CodonAligner(codon_table=table, anchor_len=anchor)
        score = float(aligner.score(amino, nucleotide))
        found, total = _alignment_list(
            aligner.align(amino, nucleotide), limit, ("protein", "nucleotides")
        )
    except Exception as exc:
        raise ConfigurationError(f"Codon alignment failed: {exc}") from exc
    for item in found:
        item["aligned_rows"] = None  # the codon layout has no single-character gapped rows
    return {
        "status": "completed",
        "operation": "align_nucleotide_to_protein",
        "codon_table": table_id,
        "score": score,
        "alignment_count": total,
        "alignments": found,
        "provenance": _c.provenance(
            method="Bio.Align.CodonAligner.align",
            coordinates="protein coordinates count residues, nucleotide coordinates count bases",
        ),
    }


# ---------------------------------------------------------------- alignment analysis


def sort_alignment_rows(
    aligned_sequences: Mapping[str, str],
    *,
    sort_by: str = "id",
    reverse: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort the rows of an alignment (``Alignment.sort``).

    ``sort_by`` is ``id`` (the Biopython default), ``sequence`` (the ungapped letters),
    ``gc_fraction`` (of the ungapped sequence), ``ungapped_length`` or ``gap_count``. Ties keep
    the input order. Returns the new row order and the rows.
    """

    rows = _rows(aligned_sequences, minimum=2)
    key_name = _c.choice(sort_by, _ROW_SORT_KEYS, "sort_by")
    flip = _c.flag(reverse, "reverse")
    seq_utils = require("Bio.SeqUtils", _EXTRA)
    alignment = _build_alignment(rows)
    keys: dict[str, Any] = {
        "id": None,
        "sequence": lambda record: str(record.seq),
        "gc_fraction": lambda record: seq_utils.gc_fraction(record.seq),
        "ungapped_length": lambda record: len(record.seq),
        "gap_count": lambda record: alignment.length - len(record.seq),
    }
    alignment.sort(key=keys[key_name], reverse=flip)
    ordered = _gapped_rows(alignment) or {}
    return {
        "status": "completed",
        "operation": "sort_alignment_rows",
        "sort_by": key_name,
        "reverse": flip,
        "order": list(ordered),
        "aligned_sequences": ordered,
        "provenance": _c.provenance(method="Bio.Align.Alignment.sort"),
    }


def alignment_substitution_counts(
    aligned_sequences: Mapping[str, str], *, verbose: bool = True
) -> dict[str, Any]:
    """Count letter substitutions between the rows of an alignment (``Alignment.substitutions``).

    Returns the alphabet and a square matrix whose entry (row letter, column letter) is the
    number of columns where the first letter is aligned to the second one (symmetrised for
    alignments of several rows as Biopython does). Gap columns are not counted.
    """

    rows = _rows(aligned_sequences, minimum=2)
    alignment = _build_alignment(rows)
    try:
        counts = alignment.substitutions
        alphabet = str(counts.alphabet)
        matrix = [[float(counts[a, b]) for b in alphabet] for a in alphabet]
    except Exception as exc:
        raise ConfigurationError(f"Could not count substitutions: {exc}") from exc
    return {
        "status": "completed",
        "operation": "alignment_substitution_counts",
        "alphabet": alphabet,
        "matrix": matrix,
        "provenance": _c.provenance(method="Bio.Align.Alignment.substitutions"),
    }


def alignment_column_frequencies(
    aligned_sequences: Mapping[str, str], *, verbose: bool = True
) -> dict[str, Any]:
    """Count each letter per alignment column (``Alignment.frequencies``).

    Returns a mapping letter to a list of counts, one per column; the gap character is ``-``.
    """

    rows = _rows(aligned_sequences, minimum=2)
    alignment = _build_alignment(rows)
    frequencies = {
        str(letter): [float(v) for v in values] for letter, values in alignment.frequencies.items()
    }
    return {
        "status": "completed",
        "operation": "alignment_column_frequencies",
        "column_count": int(alignment.length),
        "row_count": len(rows),
        "frequencies": dict(sorted(frequencies.items())),
        "provenance": _c.provenance(method="Bio.Align.Alignment.frequencies"),
    }


_COUNT_FIELDS = (
    "aligned",
    "identities",
    "mismatches",
    "positives",
    "substitution_score",
    "gap_score",
    "gaps",
    "left_gaps",
    "internal_gaps",
    "right_gaps",
    "insertions",
    "deletions",
    "left_insertions",
    "left_deletions",
    "internal_insertions",
    "internal_deletions",
    "right_insertions",
    "right_deletions",
    "open_internal_insertions",
    "open_internal_deletions",
    "extend_internal_insertions",
    "extend_internal_deletions",
)


def alignment_identity_counts(
    aligned_sequences: Mapping[str, str],
    *,
    wildcard: str | None = None,
    scoring_matrix: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count identities, mismatches and gaps of an alignment (``Alignment.counts``).

    The result lists the aligned columns, identities, mismatches, gaps and their split into
    left, internal and right gaps. ``wildcard`` (a single character such as ``N``) excludes
    that letter from the counts; ``scoring_matrix`` (a bundled substitution matrix name such as
    ``NUC.4.4``) adds the ``substitution_score`` and the number of ``positives``.
    """

    rows = _rows(aligned_sequences, minimum=2)
    alignment = _build_alignment(rows)
    scoring: Any = None
    if wildcard is not None:
        scoring = _c.text(wildcard, "wildcard", limit=1)
    if scoring_matrix is not None:
        if wildcard is not None:
            raise ConfigurationError("Give wildcard or scoring_matrix, not both.")
        matrices = require("Bio.Align.substitution_matrices", _EXTRA)
        name = _c.text(scoring_matrix, "scoring_matrix", limit=100)
        if name not in matrices.load():
            raise ConfigurationError(f"unknown substitution matrix {name!r}.")
        scoring = matrices.load(name)
    try:
        counts = alignment.counts(scoring)
    except Exception as exc:
        raise ConfigurationError(f"Could not count the alignment: {exc}") from exc
    values = {
        name: _c.plain(getattr(counts, name)) for name in _COUNT_FIELDS if hasattr(counts, name)
    }
    values = {key: value for key, value in values.items() if value is not None}
    return {
        "status": "completed",
        "operation": "alignment_identity_counts",
        "row_count": len(rows),
        "column_count": int(alignment.length),
        "counts": values,
        "provenance": _c.provenance(method="Bio.Align.Alignment.counts"),
    }


def alignment_coordinate_maps(
    aligned_sequences: Mapping[str, str], *, verbose: bool = True
) -> dict[str, Any]:
    """Map between alignment columns and sequence positions.

    This covers ``Alignment.aligned``, ``indices`` and ``inverse_indices``.

    ``indices`` has one list per row with, for each column, the 0-based position in the
    ungapped sequence, or -1 at a gap. ``inverse_indices`` has one list per row with, for each
    ungapped position, its column. ``aligned_blocks`` (two rows only) lists the aligned
    block intervals ``[[start, end], ...]`` per sequence; it is null for more than two rows.
    """

    rows = _rows(aligned_sequences, minimum=2)
    alignment = _build_alignment(rows)
    try:
        blocks = alignment.aligned.tolist() if len(rows) == 2 else None
    except NotImplementedError:
        blocks = None
    return {
        "status": "completed",
        "operation": "alignment_coordinate_maps",
        "names": list(rows),
        "column_count": int(alignment.length),
        "indices": {
            name: [int(v) for v in line]
            for name, line in zip(rows, alignment.indices.tolist(), strict=True)
        },
        "inverse_indices": {
            name: [int(v) for v in line]
            for name, line in zip(rows, alignment.inverse_indices, strict=True)
        },
        "aligned_blocks": blocks,
        "provenance": _c.provenance(
            method="Bio.Align.Alignment.aligned, .indices and .inverse_indices"
        ),
    }


def parse_printed_alignment(
    aligned_lines: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Infer the ungapped sequences and block coordinates from printed alignment lines.

    ``aligned_lines`` holds one gapped line per sequence, all with the same number of
    columns, ``-`` for a gap (``Alignment.parse_printed_alignment``). Returned are the
    ungapped sequences and the coordinate matrix (rows are sequences, columns are block
    boundaries, 0-based).
    """

    if isinstance(aligned_lines, str) or not isinstance(aligned_lines, Sequence):
        raise ConfigurationError("aligned_lines must be a list of strings.")
    if not 1 <= len(aligned_lines) <= _MAX_ROWS:
        raise ConfigurationError(f"aligned_lines must hold 1 to {_MAX_ROWS} lines.")
    cleaned = []
    for index, line in enumerate(aligned_lines):
        if not isinstance(line, str) or not line or not _ROW_SYMBOLS.match(line):
            raise ConfigurationError(f"aligned_lines[{index}] must be a gapped sequence.")
        cleaned.append(line)
    if len({len(line) for line in cleaned}) != 1:
        raise ConfigurationError("aligned_lines must all have the same number of columns.")
    align = require("Bio.Align", _EXTRA)
    try:
        sequences, coordinates = align.Alignment.parse_printed_alignment(
            [line.encode("ascii") for line in cleaned]
        )
    except Exception as exc:
        raise ConfigurationError(f"Could not parse the printed alignment: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_printed_alignment",
        "sequences": [sequence.decode("ascii") for sequence in sequences],
        "coordinates": [[int(v) for v in row] for row in coordinates.tolist()],
        "provenance": _c.provenance(method="Bio.Align.Alignment.parse_printed_alignment"),
    }


def _pair(item: object, label: str) -> tuple[str, str]:
    if not isinstance(item, Mapping) or "target" not in item or "query" not in item:
        raise ConfigurationError(
            f"{label} must be a mapping with 'target' and 'query' gapped rows."
        )
    rows = _rows({"target": item["target"], "query": item["query"]}, label, minimum=2)
    return rows["target"], rows["query"]


def merge_alignments_by_reference(
    reference_alignments: Mapping[str, Mapping[str, str]], *, verbose: bool = True
) -> dict[str, Any]:
    """Combine pairwise alignments that share one reference sequence into a multiple alignment.

    ``reference_alignments`` maps a query name to ``{"target": gapped reference, "query":
    gapped query}``. After removing gaps, every reference must be the same sequence
    (``Alignment.from_alignments_with_same_reference``). The result has the reference first and
    the queries after it, with columns added where a query inserts letters.
    """

    if not isinstance(reference_alignments, Mapping) or not reference_alignments:
        raise ConfigurationError("reference_alignments must be a non-empty mapping.")
    if len(reference_alignments) > 1000:
        raise ConfigurationError("reference_alignments may hold at most 1000 alignments.")
    align = require("Bio.Align", _EXTRA)
    pairs = []
    for name, item in reference_alignments.items():
        if not isinstance(name, str) or not _NAME.match(name) or name == "reference":
            raise ConfigurationError(f"invalid query name {name!r}.")
        target_row, query_row = _pair(item, f"reference_alignments[{name!r}]")
        built = _build_alignment({"reference": target_row, name: query_row})
        pairs.append(built)
    try:
        merged = align.Alignment.from_alignments_with_same_reference(pairs)
    except Exception as exc:
        raise ConfigurationError(f"The alignments cannot be merged: {exc}") from exc
    merged_rows = {str(sequence.id): str(merged[i]) for i, sequence in enumerate(merged.sequences)}
    return {
        "status": "completed",
        "operation": "merge_alignments_by_reference",
        "row_count": len(merged_rows),
        "column_count": int(merged.length),
        "aligned_sequences": merged_rows,
        "provenance": _c.provenance(
            method="Bio.Align.Alignment.from_alignments_with_same_reference"
        ),
    }


def map_alignment_through(
    alignment: Mapping[str, str],
    through: Sequence[Mapping[str, str]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map (lift over) an alignment through one or more alignments of its query sequence.

    ``alignment`` is ``{"target": gapped row, "query": gapped row}``, for example a chromosome
    against a transcript. ``through`` lists alignments of the same form whose target is the
    first alignment's query (ungapped), for example reads aligned to the transcript. One
    alignment uses ``Alignment.map``; several use ``Alignment.mapall``. The result is the
    alignment of the mapped sequences to the first target, with block coordinates and rows.
    """

    first_target, first_query = _pair(alignment, "alignment")
    candidate: object = through
    if isinstance(candidate, str) or not isinstance(candidate, Sequence) or not candidate:
        raise ConfigurationError("through must be a non-empty list of alignments.")
    if len(through) > 1000:
        raise ConfigurationError("through may hold at most 1000 alignments.")
    base = _build_alignment({"target": first_target, "query": first_query})
    others = [
        _build_alignment(dict(zip(("target", "query"), _pair(item, f"through[{i}]"), strict=True)))
        for i, item in enumerate(through)
    ]
    try:
        mapped = base.map(others[0]) if len(others) == 1 else base.mapall(others)
    except Exception as exc:
        raise ConfigurationError(f"The alignment cannot be mapped: {exc}") from exc
    summary = _alignment_summary(mapped, 0)
    return {
        "status": "completed",
        "operation": "map_alignment_through",
        "method": "map" if len(others) == 1 else "mapall",
        "coordinates": summary["coordinates"],
        "aligned_rows": summary["aligned_rows"],
        "provenance": _c.provenance(
            method="Bio.Align.Alignment.map" if len(others) == 1 else "Bio.Align.Alignment.mapall"
        ),
    }


# ---------------------------------------------------------------- Nexus matrix operations


def _nexus(rows: Mapping[str, str]) -> Any:
    nexus = require("Bio.Nexus.Nexus", _EXTRA)
    length = len(next(iter(rows.values())))
    lines = "\n".join(f"{name} {value}" for name, value in rows.items())
    text = (
        "#NEXUS\nbegin data;\n"
        f"dimensions ntax={len(rows)} nchar={length};\n"
        "format datatype=dna missing=? gap=-;\nmatrix\n"
        f"{lines}\n;\nend;\n"
    )
    try:
        return nexus.Nexus(text)
    except Exception as exc:
        raise ConfigurationError(f"The rows are not a valid DNA matrix: {exc}") from exc


def _taxa(value: object, rows: Mapping[str, str], name: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of sequence names.")
    unknown = [item for item in value if item not in rows]
    if unknown:
        raise ConfigurationError(f"{name} names unknown sequences: {unknown}")
    return list(value)


def _sites(value: object, length: int, name: str) -> list[int]:
    if value is None:
        return []
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of 0-based column indices.")
    return [_c.integer(item, name, minimum=0, maximum=length - 1) for item in value]


def find_alignment_columns(
    aligned_sequences: Mapping[str, str],
    *,
    kind: str = "constant",
    delete_taxa: Sequence[str] | None = None,
    exclude_sites: Sequence[int] | None = None,
    include_missing: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the constant or the gap-only columns of a DNA alignment (``Bio.Nexus``).

    ``kind="constant"`` returns the columns where all (undeleted) sequences agree, treating
    ``?`` as compatible and ambiguity codes as sets of bases (``Nexus.constant``);
    ``kind="gap_only"`` returns the columns that hold only gaps (``Nexus.gaponly``), counting
    ``?`` as a gap when ``include_missing`` is true. ``delete_taxa`` leaves sequences out and
    ``exclude_sites`` leaves columns out of the constant search. Columns are 0-based.
    """

    rows = _rows(aligned_sequences, minimum=2, nexus=True)
    selected = _c.choice(kind, ("constant", "gap_only"), "kind")
    removed = _taxa(delete_taxa, rows, "delete_taxa")
    length = len(next(iter(rows.values())))
    skipped = _sites(exclude_sites, length, "exclude_sites")
    missing = _c.flag(include_missing, "include_missing")
    matrix = _nexus(rows)
    if selected == "constant":
        found = matrix.constant(delete=removed, exclude=skipped)
        method = "Bio.Nexus.Nexus.Nexus.constant"
    else:
        found = matrix.gaponly(include_missing=missing)
        method = "Bio.Nexus.Nexus.Nexus.gaponly"
    columns = sorted(int(v) for v in (found or []))
    return {
        "status": "completed",
        "operation": "find_alignment_columns",
        "kind": selected,
        "column_count": length,
        "columns": columns,
        "found": len(columns),
        "provenance": _c.provenance(method=method),
    }


def site_character_states(
    aligned_sequences: Mapping[str, str],
    *,
    sites: Sequence[int] | None = None,
    delete_taxa: Sequence[str] | None = None,
    narrow: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the character states of alignment columns (``Nexus.cstatus``).

    For each column in ``sites`` (all columns when null; 0-based) returns the sorted list of
    bases seen in the undeleted sequences, with ambiguity codes expanded. ``narrow=True``
    follows PAUP (``A C ?`` gives ``A C``; all ``?`` gives ``?``); ``narrow=False`` expands ``?``
    and gaps to all states.
    """

    rows = _rows(aligned_sequences, minimum=2, nexus=True)
    length = len(next(iter(rows.values())))
    removed = _taxa(delete_taxa, rows, "delete_taxa")
    chosen = _sites(sites, length, "sites") if sites is not None else list(range(length))
    narrowed = _c.flag(narrow, "narrow")
    matrix = _nexus(rows)
    states = [
        {"site": site, "states": matrix.cstatus(site, delete=removed, narrow=narrowed)}
        for site in chosen
    ]
    return {
        "status": "completed",
        "operation": "site_character_states",
        "column_count": length,
        "sites": states,
        "provenance": _c.provenance(method="Bio.Nexus.Nexus.Nexus.cstatus"),
    }


def subset_alignment(
    aligned_sequences: Mapping[str, str],
    *,
    delete_taxa: Sequence[str] | None = None,
    exclude_sites: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete sequences and exclude columns from a DNA alignment (``Nexus.crop_matrix``).

    ``delete_taxa`` names the sequences to drop and ``exclude_sites`` the 0-based columns to
    remove. At least one of them must be given and at least one sequence must remain.
    """

    rows = _rows(aligned_sequences, minimum=2, nexus=True)
    length = len(next(iter(rows.values())))
    removed = _taxa(delete_taxa, rows, "delete_taxa")
    skipped = _sites(exclude_sites, length, "exclude_sites")
    if not removed and not skipped:
        raise ConfigurationError("Give delete_taxa or exclude_sites.")
    if len(set(removed)) >= len(rows):
        raise ConfigurationError("At least one sequence must remain.")
    cropped = _nexus(rows).crop_matrix(delete=removed, exclude=skipped)
    result = {str(name): str(value) for name, value in cropped.items()}
    return {
        "status": "completed",
        "operation": "subset_alignment",
        "row_count": len(result),
        "column_count": len(next(iter(result.values()))) if result else 0,
        "aligned_sequences": result,
        "provenance": _c.provenance(method="Bio.Nexus.Nexus.Nexus.crop_matrix"),
    }


def insert_alignment_gaps(
    aligned_sequences: Mapping[str, str],
    position: int,
    *,
    count: int = 1,
    left_greedy: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert gap columns into a DNA alignment before a 0-based column (``Nexus.insert_gap``).

    ``count`` columns of ``-`` are added to every sequence. When the insertion point lies
    inside a run of gaps, the new gaps are placed at the right end of the run, or at the left
    end with ``left_greedy=True``.
    """

    rows = _rows(aligned_sequences, minimum=2, nexus=True)
    length = len(next(iter(rows.values())))
    where = _c.integer(position, "position", minimum=0, maximum=length)
    number = _c.integer(count, "count", minimum=1, maximum=100_000)
    greedy = _c.flag(left_greedy, "left_greedy")
    matrix = _nexus(rows)
    matrix.insert_gap(where, number, greedy)
    result = {str(name): str(value) for name, value in matrix.matrix.items()}
    return {
        "status": "completed",
        "operation": "insert_alignment_gaps",
        "column_count": int(matrix.nchar),
        "aligned_sequences": result,
        "provenance": _c.provenance(method="Bio.Nexus.Nexus.Nexus.insert_gap"),
    }


__all__ = [
    "align_nucleotide_to_protein",
    "alignment_column_frequencies",
    "alignment_coordinate_maps",
    "alignment_identity_counts",
    "alignment_substitution_counts",
    "convert_alignment_file",
    "filter_search_results",
    "find_alignment_columns",
    "fogsaa_global_alignment",
    "index_search_results",
    "insert_alignment_gaps",
    "map_alignment_through",
    "merge_alignments_by_reference",
    "parse_printed_alignment",
    "read_alignment_records",
    "read_search_results",
    "search_maf_index",
    "site_character_states",
    "sort_alignment_rows",
    "sort_search_results",
    "splice_maf_index",
    "subset_alignment",
    "write_aligned_sequences",
]

_register_public_verbose(globals(), tuple(__all__))
