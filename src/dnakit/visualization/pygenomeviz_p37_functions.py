"""Genome file parsers, alignment-coordinate tools, BLAST/LAST aligners and feature drawing after pyGenomeViz (pack P37).

pyGenomeViz 1.x ships parsers for GenBank, GFF and FASTA files and an ``AlignCoord`` class for genome
alignments. The pinned ``pygenomeviz<1`` extra does not have them, so the parsing, filtering and
file-format code here is written with Biopython and the standard library and follows the pyGenomeViz
1.7.0 behaviour; it needs no pyGenomeViz. Coordinates of features and alignment blocks are 0-based
half-open, as in pyGenomeViz. Alignment coordinates are plain dictionaries with the keys of
``AlignCoord`` (``query_id``, ``query_name``, ``query_start``, ``query_end``, ``ref_id``, ``ref_name``,
``ref_start``, ``ref_end``, ``identity``, ``evalue``). One deliberate difference: ``ref_length`` is
``abs(ref_end - ref_start)``; pyGenomeViz 1.7.0 computes it from the query columns, which is a bug.

``align_genomes_blast`` and ``align_genomes_last`` run BLAST+ or LAST with fixed argument lists and no
shell. ``draw_genome_segment_features`` and ``export_genome_html_viewer`` need pyGenomeViz 0.4.x.
"""

from __future__ import annotations

import bz2
import csv
import gzip
import io
import re
import zipfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.visualization.pygenomeviz_gap_functions import (
    _check_new_folder,
    _executable,
    _require_isolated,
    _Steps,
)

_EXTRA = "bio-pygenomeviz"
_BIOPYTHON_EXTRA = "bio-biopython"
_MAX_TAIL = 4000
_LARGE_OUTPUT = 100_000_000
_FILE_FORMATS = ("auto", "genbank", "gff", "fasta")
_ALIGN_FORMATS = ("blast", "mummer", "progressive_mauve")
_SEQTYPES = ("nucleotide", "protein")
_BLAST_HINT = (
    "Install BLAST+ (conda install -c bioconda blast, or apt install ncbi-blast+) and put makeblastdb, "
    "blastn and tblastx on PATH or pass blast_dir."
)
_LAST_HINT = (
    "Install LAST (conda install -c bioconda last, or https://gitlab.com/mcfrith/last) and put lastdb, "
    "lastal, last-split and maf-convert on PATH or pass last_dir."
)
_TSV_HEADER = (
    "QUERY_ID", "QUERY_NAME", "QUERY_START", "QUERY_END", "QUERY_LENGTH",
    "REF_ID", "REF_NAME", "REF_START", "REF_END", "REF_LENGTH", "IDENTITY", "EVALUE",
)  # fmt: skip
_COORD_KEYS = (
    "query_id", "query_name", "query_start", "query_end",
    "ref_id", "ref_name", "ref_start", "ref_end", "identity", "evalue",
)  # fmt: skip


def _provenance(**notes: object) -> dict[str, Any]:
    return {"reference": "pyGenomeViz 1.7.0 semantics, implemented with Biopython", **notes}


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _strand(value: object) -> int | None:
    if value is None:
        return None
    if value not in (1, -1, 0) or isinstance(value, bool):
        raise ConfigurationError("strand must be None, 1, -1 or 0.")
    return int(value)


def _types(value: str | Sequence[str] | None, name: str) -> list[str] | None:
    if value is None:
        return None
    items = [value] if isinstance(value, str) else list(value)
    if not items or not all(isinstance(item, str) and item for item in items):
        raise ConfigurationError(f"{name} must be None, a feature type or a list of feature types.")
    return items


def _range(value: Sequence[int] | None, name: str) -> tuple[int, int] | None:
    if value is None:
        return None
    if (
        isinstance(value, (str, bytes))
        or len(value) != 2
        or not all(isinstance(item, int) and not isinstance(item, bool) for item in value)
    ):
        raise ConfigurationError(f"{name} must be None or two integers (start, end).")
    return (min(value), max(value))


def _positive(value: object, name: str, *, allow_none: bool = False) -> int | None:
    if value is None and allow_none:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _new_file(value: str | Path, name: str, suffixes: tuple[str, ...] | None = None) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists: {path.name}; choose a new file.")
    if suffixes is not None and path.suffix.lower() not in suffixes:
        raise ConfigurationError(f"{name} must end with one of {list(suffixes)}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder of {name} does not exist: {path.parent}.")
    return path


def _open_text(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    if path.suffix == ".bz2":
        return bz2.open(path, mode="rt", encoding="utf-8")
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            return io.StringIO(archive.read(archive.namelist()[0]).decode("utf-8"))
    return open(path, encoding="utf-8")


def _stem(path: Path) -> str:
    return (
        path.with_suffix("").with_suffix("").name
        if path.suffix in (".gz", ".bz2", ".zip")
        else path.with_suffix("").name
    )


def _detect(path: Path, file_format: str) -> str:
    _choice(file_format, _FILE_FORMATS, "file_format")
    if file_format != "auto":
        return file_format
    name = path.name.lower()
    for suffix in (".gz", ".bz2", ".zip"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith((".gb", ".gbk", ".gbff", ".genbank")):
        return "genbank"
    if name.endswith((".gff", ".gff3")):
        return "gff"
    if name.endswith((".fa", ".fna", ".fasta", ".fas")):
        return "fasta"
    raise ConfigurationError(f"cannot tell the format of {path.name}; pass file_format.")


def _biopython() -> tuple[Any, Any]:
    seqio = require("Bio.SeqIO", _BIOPYTHON_EXTRA)
    feature = require("Bio.SeqFeature", _BIOPYTHON_EXTRA)
    return seqio, feature


def _records(path: Path, file_format: str) -> list[Any]:
    seqio, _ = _biopython()
    try:
        with _open_text(path) as handle:
            records = list(seqio.parse(handle, file_format))
    except Exception as exc:
        raise ConfigurationError(
            f"{path.name} could not be parsed as {file_format}: {exc}"
        ) from exc
    if not records:
        raise ConfigurationError(f"{path.name} has no {file_format} records.")
    return records


def _feature_row(
    seqid: str,
    kind: str,
    start: int,
    end: int,
    strand: int | None,
    identifier: str,
    qualifiers: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "seqid": seqid,
        "type": kind,
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "id": identifier,
        "qualifiers": {
            str(k): [str(v) for v in (values if isinstance(values, list) else [values])]
            for k, values in qualifiers.items()
        },
    }


def _straddles(feature: Any) -> bool:
    location = feature.location
    if location.strand == -1:
        start, end = int(location.parts[-1].start), int(location.parts[0].end)
    else:
        start, end = int(location.parts[0].start), int(location.parts[-1].end)
    return start > end


# GFF parsing (pyGenomeViz GffRecord semantics) ------------------------------------------------


def _gff_records(path: Path) -> tuple[list[dict[str, Any]], dict[str, int], list[str]]:
    with _open_text(path) as handle:
        lines = handle.read().splitlines()
    records: list[dict[str, Any]] = []
    for line in lines:
        if line.startswith("#") or len(line.split("\t")) < 9:
            continue
        parts = line.split("\t")[:9]
        try:
            start, end = int(parts[3]), int(parts[4])
            score = None if parts[5] in (".", "") else float(parts[5])
            phase = None if parts[7] in (".", "") else int(parts[7])
        except ValueError as exc:
            raise ConfigurationError(
                f"{path.name} has an invalid GFF line: {line[:80]!r}."
            ) from exc
        attrs: dict[str, list[str]] = {}
        for attribute in parts[8].split(";"):
            if attribute and "=" in attribute:
                key, value = attribute.split("=", 1)
                attrs[key] = value.split(",")
        records.append(
            {
                "seqid": parts[0], "source": parts[1], "type": parts[2], "start": start, "end": end,
                "score": score, "strand": 1 if parts[6] == "+" else -1 if parts[6] == "-" else 0,
                "phase": phase, "attrs": attrs,
            }
        )  # fmt: skip
    if not records:
        raise ConfigurationError(f"{path.name} has no GFF records.")
    seqids = list(dict.fromkeys(record["seqid"] for record in records))
    sizes: dict[str, int] = {}
    for seqid in seqids:
        start = end = None
        for line in lines:
            fields = line.split()
            if line.startswith("##sequence-region") and len(fields) == 4 and fields[1] == seqid:
                start, end = int(fields[2]) - 1, int(fields[3])
                break
        if start is None or end is None:
            start, end = 0, max(record["end"] for record in records if record["seqid"] == seqid)
        sizes[seqid] = end - start
    return records, sizes, seqids


def _gff_feature(record: dict[str, Any]) -> dict[str, Any]:
    identifier = record["attrs"].get("ID", [""])[0]
    return _feature_row(
        record["seqid"],
        record["type"],
        record["start"] - 1,
        record["end"],
        record["strand"],
        identifier,
        record["attrs"],
    )


def _gff_select(
    records: Sequence[dict[str, Any]],
    types: list[str] | None,
    strand: int | None,
    target_range: tuple[int, int] | None,
) -> list[dict[str, Any]]:
    chosen = []
    for record in records:
        if types is not None and record["type"] not in types:
            continue
        if strand is not None and record["strand"] != strand:
            continue
        if (
            target_range is not None
            and not target_range[0] <= record["start"] <= record["end"] <= target_range[1]
        ):
            continue
        chosen.append(record)
    return chosen


# Sequence and feature readers -------------------------------------------------------------------


def read_contig_sequences(
    source: str | Path,
    *,
    file_format: str = "auto",
    include_sequences: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the contig (record) sequences of a GenBank or FASTA file as an ID to sequence listing.

    Covers ``Genbank.get_seqid2seq`` (``file_format='genbank'``) and ``Fasta.get_seqid2seq`` /
    ``get_seqid2record`` (``'fasta'``). FASTA files with duplicate IDs are rejected, as pyGenomeViz does.
    """

    path = input_path(source)
    kind = _detect(path, file_format)
    if kind not in ("genbank", "fasta"):
        raise ConfigurationError("file_format must be genbank or fasta (GFF has no sequences).")
    records = _records(path, kind)
    ids = [str(record.id) for record in records]
    if len(set(ids)) != len(ids):
        raise ConfigurationError("Duplicate IDs are contained in the file.")
    contigs = [
        {"id": str(record.id), "length": len(record.seq), "description": str(record.description)}
        | ({"sequence": str(record.seq)} if include_sequences else {})
        for record in records
    ]
    return {
        "status": "completed",
        "operation": "read_contig_sequences",
        "file_format": kind,
        "name": _stem(path),
        "contig_count": len(contigs),
        "total_length": sum(item["length"] for item in contigs),
        "contigs": contigs,
        "provenance": _provenance(),
    }


def read_gff_contig_sizes(source: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Contig lengths from GFF ``##sequence-region`` headers (``Gff.get_seqid2size``).

    A contig without a header gets the largest end coordinate of its features, which can be shorter than the contig.
    """

    path = input_path(source)
    records, sizes, seqids = _gff_records(path)
    headers = {
        line.split()[1]
        for line in _open_text(path).read().splitlines()
        if line.startswith("##sequence-region")
    }
    return {
        "status": "completed",
        "operation": "read_gff_contig_sizes",
        "contigs": [
            {"id": seqid, "size": sizes[seqid], "from_header": seqid in headers} for seqid in seqids
        ],
        "record_count": len(records),
        "provenance": _provenance(),
    }


def list_contig_features(
    source: str | Path,
    *,
    file_format: str = "auto",
    feature_type: str | Sequence[str] | None = "CDS",
    strand: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the features of a type and strand for every contig of a GenBank or GFF file.

    Covers ``Genbank.get_seqid2features`` and ``Gff.get_seqid2features``. ``feature_type=None`` keeps every type.
    GenBank features that wrap over the genome origin are dropped. ``strand`` is 1, -1 or None (any).
    Coordinates are 0-based half-open.
    """

    path = input_path(source)
    kind = _detect(path, file_format)
    types, wanted = _types(feature_type, "feature_type"), _strand(strand)
    per_contig: dict[str, list[dict[str, Any]]] = {}
    if kind == "genbank":
        for record in _records(path, "genbank"):
            for feature in record.features:
                if feature.location is None or (types is not None and feature.type not in types):
                    continue
                if wanted is not None and feature.location.strand != wanted:
                    continue
                if _straddles(feature):
                    continue
                per_contig.setdefault(str(record.id), []).append(
                    _feature_row(
                        str(record.id),
                        feature.type,
                        int(feature.location.start),
                        int(feature.location.end),
                        feature.location.strand,
                        "",
                        feature.qualifiers,
                    )
                )
    elif kind == "gff":
        records, _, seqids = _gff_records(path)
        per_contig = {seqid: [] for seqid in seqids}
        for record in _gff_select(records, types, wanted, None):
            per_contig[record["seqid"]].append(_gff_feature(record))
    else:
        raise ConfigurationError("file_format must be genbank or gff.")
    return {
        "status": "completed",
        "operation": "list_contig_features",
        "file_format": kind,
        "feature_type": types,
        "strand": wanted,
        "contigs": {seqid: rows for seqid, rows in per_contig.items()},
        "feature_count": sum(len(rows) for rows in per_contig.values()),
        "provenance": _provenance(),
    }


def extract_genome_features(
    source: str | Path,
    *,
    file_format: str = "auto",
    feature_type: str | Sequence[str] | None = "CDS",
    strand: int | None = None,
    target_range: Sequence[int] | None = None,
    seqid: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract features of one contig, optionally inside a range (``Genbank.extract_features``, ``Gff.extract_features``).

    GenBank uses the first record that has a matching feature and GFF the first seqid, unless ``seqid`` is given.
    ``target_range`` is ``(start, end)``: GenBank compares 0-based half-open feature coordinates, GFF the 1-based
    record coordinates, as pyGenomeViz does.
    """

    path = input_path(source)
    kind = _detect(path, file_format)
    window = _range(target_range, "target_range")
    listing = list_contig_features(
        path, file_format=kind, feature_type=feature_type, strand=strand, verbose=False
    )
    contigs = listing["contigs"]
    if kind == "gff":
        records, _, seqids = _gff_records(path)
        chosen = seqid if seqid is not None else seqids[0]
        if chosen not in seqids:
            raise ConfigurationError(f"seqid {chosen!r} was not found in {path.name}.")
        rows = [
            _gff_feature(item)
            for item in _gff_select(
                [r for r in records if r["seqid"] == chosen],
                _types(feature_type, "feature_type"),
                _strand(strand),
                window,
            )
        ]
    else:
        nonempty = [name for name, items in contigs.items() if items]
        chosen = seqid if seqid is not None else (nonempty[0] if nonempty else None)
        if chosen is None:
            raise ConfigurationError("the GenBank file has no matching feature.")
        rows = [
            row
            for row in contigs.get(chosen, [])
            if window is None or window[0] <= row["start"] <= row["end"] <= window[1]
        ]
    return {
        "status": "completed",
        "operation": "extract_genome_features",
        "file_format": kind,
        "seqid": chosen,
        "target_range": list(window) if window else None,
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(),
    }


def extract_gff_exon_features(
    source: str | Path,
    *,
    feature_type: str = "mRNA",
    strand: int | None = None,
    target_range: Sequence[int] | None = None,
    seqid: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build exon-structured features from the ``ID``/``Parent`` relation of a feature type and its exons (``Gff.extract_exon_features``).

    A parent with one exon keeps its own span; with two or more exons it gets a ``parts`` list of exon
    ``[start, end]`` pairs (0-based half-open, sorted); a parent without exons is skipped. Exons are matched to the
    most recent parent of ``feature_type`` in file order, as pyGenomeViz does.
    """

    path = input_path(source)
    if not isinstance(feature_type, str) or not feature_type:
        raise ConfigurationError("feature_type must be a feature type name such as mRNA.")
    wanted, window = _strand(strand), _range(target_range, "target_range")
    records, _, seqids = _gff_records(path)
    chosen = seqid if seqid is not None else seqids[0]
    if chosen not in seqids:
        raise ConfigurationError(f"seqid {chosen!r} was not found in {path.name}.")
    parent_id = None
    parents: dict[str, dict[str, Any]] = {}
    exons: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in (r for r in records if r["seqid"] == chosen):
        if record["type"] == feature_type:
            parent_id = record["attrs"].get("ID", [None])[0]
            if parent_id is None:
                continue
            parents[parent_id] = record
        if (
            record["type"] == "exon"
            and parent_id is not None
            and parent_id == record["attrs"].get("Parent", [None])[0]
        ):
            exons[parent_id].append(record)
    rows: list[dict[str, Any]] = []
    for name, parent in parents.items():
        found = sorted(exons.get(name, []), key=lambda item: item["start"])
        if not found:
            continue
        row = _gff_feature(parent)
        if len(found) >= 2:
            row["parts"] = [[item["start"] - 1, item["end"]] for item in found]
            row["start"], row["end"] = row["parts"][0][0], row["parts"][-1][1]
        else:
            row["parts"] = [[row["start"], row["end"]]]
        if wanted is not None and row["strand"] != wanted:
            continue
        if window is not None and not window[0] <= row["start"] <= row["end"] <= window[1]:
            continue
        rows.append(row)
    return {
        "status": "completed",
        "operation": "extract_gff_exon_features",
        "seqid": chosen,
        "feature_type": feature_type,
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(),
    }


def select_features_within_range(
    features: Sequence[Mapping[str, Any]],
    target_range: Sequence[int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the features that lie completely inside ``(start, end)`` (``pygenomeviz.utils.extract_features_within_range``).

    Each feature is a mapping with integer ``start`` and ``end`` (for example rows from ``extract_genome_features``).
    """

    window = _range(target_range, "target_range")
    if window is None or isinstance(features, (str, bytes)):
        raise ConfigurationError("target_range and a list of features are required.")
    kept = []
    for feature in features:
        if not isinstance(feature, Mapping) or not all(
            isinstance(feature.get(k), int) for k in ("start", "end")
        ):
            raise ConfigurationError("every feature needs integer start and end.")
        if window[0] <= feature["start"] <= feature["end"] <= window[1]:
            kept.append(dict(feature))
    return {
        "status": "completed",
        "operation": "select_features_within_range",
        "target_range": list(window),
        "input_count": len(features),
        "feature_count": len(kept),
        "features": kept,
    }


def write_genome_fasta(
    source: str | Path,
    output_path: str | Path,
    *,
    file_format: str = "auto",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the genome sequences of a GenBank or FASTA file as one-line-per-sequence FASTA (``write_genome_fasta``)."""

    path = input_path(source)
    target = _new_file(output_path, "output_path")
    contigs = read_contig_sequences(path, file_format=file_format, verbose=False)["contigs"]
    target.write_text(
        "".join(f">{item['id']}\n{item['sequence']}\n" for item in contigs), encoding="utf-8"
    )
    return {
        "status": "completed",
        "operation": "write_genome_fasta",
        "output_path": str(target),
        "sequence_count": len(contigs),
        "total_length": sum(item["length"] for item in contigs),
    }


def write_genbank_cds_fasta(
    source: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the CDS features of a GenBank file with a translation to protein FASTA (``Genbank.write_cds_fasta``).

    Each record is named ``GENE<n>_<protein_id>|<seqid>|<start>_<end>_<strand>|`` with the product as description.
    CDS features without a ``/translation`` and features that wrap the origin are skipped. For the nucleotide
    sequence of every feature use ``extract_gffutils_feature_sequence`` or Biopython.
    """

    path = input_path(source)
    target = _new_file(output_path, "output_path")
    listing = list_contig_features(path, file_format="genbank", feature_type="CDS", verbose=False)
    lines: list[str] = []
    written = counter = 0
    for seqid, rows in listing["contigs"].items():
        for row in rows:
            translation = row["qualifiers"].get("translation", [None])[0]
            if translation is None:
                continue
            strand = -1 if row["strand"] == -1 else 1
            location_id = f"|{seqid}|{row['start']}_{row['end']}_{strand}|"
            protein_id = row["qualifiers"].get("protein_id", [None])[0]
            name = (
                f"GENE{counter:06d}{location_id}"
                if protein_id is None
                else f"GENE{counter:06d}_{protein_id}{location_id}"
            )
            counter += 1
            product = row["qualifiers"].get("product", [""])[0]
            lines.append(f">{name} {product}".rstrip() if product else f">{name}")
            lines.append(translation)
            written += 1
    target.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return {
        "status": "completed",
        "operation": "write_genbank_cds_fasta",
        "output_path": str(target),
        "cds_written": written,
    }


def calc_genbank_gc_content(
    source: str | Path,
    *,
    window_size: int | None = None,
    step_size: int | None = None,
    seqid: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """GC content (percent) in a sliding window along a GenBank genome (``Genbank.calc_gc_content``).

    Positions run from 0 to the sequence length in steps of ``step_size`` (default length/1000); each value is the GC
    percentage of the window centred on the position (default window length/500). The first record is used unless ``seqid`` is given.
    """

    path = input_path(source)
    records = _records(path, "genbank")
    record = records[0] if seqid is None else next((r for r in records if r.id == seqid), None)
    if record is None:
        raise ConfigurationError(f"seqid {seqid!r} was not found in {path.name}.")
    window = _positive(window_size, "window_size", allow_none=True)
    step = _positive(step_size, "step_size", allow_none=True)
    sequence = str(record.seq)
    size = len(sequence)
    window = int(size / 500) if window is None else window
    step = int(size / 1000) if step is None else step
    if window == 0 or step == 0:
        window, step = size, int(size / 2)
    if step == 0:
        raise ConfigurationError("the sequence is too short for a sliding window.")
    seq_utils = require("Bio.SeqUtils", _BIOPYTHON_EXTRA)
    positions = [*range(0, size, step), size]
    values = []
    for position in positions:
        subsequence = sequence[
            max(position - int(window / 2), 0) : min(position + int(window / 2), size)
        ]
        values.append(float(seq_utils.gc_fraction(subsequence) * 100) if subsequence else 0.0)
    return {
        "status": "completed",
        "operation": "calc_genbank_gc_content",
        "seqid": str(record.id),
        "sequence_length": size,
        "window_size": window,
        "step_size": step,
        "positions": positions,
        "gc_content": values,
        "genome_gc_content": float(seq_utils.gc_fraction(sequence) * 100) if sequence else 0.0,
        "provenance": _provenance(),
    }


# Alignment coordinates -----------------------------------------------------------------------------


def _coord(query_id: str, query_name: str, query_start: int, query_end: int, ref_id: str, ref_name: str,
           ref_start: int, ref_end: int, identity: float | None, evalue: float | None) -> dict[str, Any]:  # fmt: skip
    return {
        "query_id": query_id, "query_name": query_name, "query_start": int(query_start), "query_end": int(query_end),
        "ref_id": ref_id, "ref_name": ref_name, "ref_start": int(ref_start), "ref_end": int(ref_end),
        "identity": identity, "evalue": evalue,
    }  # fmt: skip


def _checked_coords(coords: object) -> list[dict[str, Any]]:
    if isinstance(coords, (str, bytes)) or not isinstance(coords, Sequence):
        raise ConfigurationError("coords must be a list of alignment coordinate mappings.")
    checked = []
    for item in coords:
        if not isinstance(item, Mapping) or not set(_COORD_KEYS[:8]) <= set(item):
            raise ConfigurationError(f"every coordinate needs the keys {list(_COORD_KEYS[:8])}.")
        try:
            row = _coord(
                str(item["query_id"]), str(item["query_name"]), item["query_start"], item["query_end"],
                str(item["ref_id"]), str(item["ref_name"]), item["ref_start"], item["ref_end"],
                None if item.get("identity") is None else float(item["identity"]),
                None if item.get("evalue") is None else float(item["evalue"]),
            )  # fmt: skip
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("coordinates and scores must be numbers.") from exc
        checked.append(row)
    return checked


def _lengths(row: Mapping[str, Any]) -> tuple[int, int]:
    return abs(row["query_end"] - row["query_start"]), abs(row["ref_end"] - row["ref_start"])


def _tsv_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    qlen, rlen = _lengths(row)
    return (row["query_id"], row["query_name"], row["query_start"], row["query_end"], qlen, row["ref_id"], row["ref_name"],
            row["ref_start"], row["ref_end"], rlen, row["identity"], row["evalue"])  # fmt: skip


def _contained(inner: Mapping[str, Any], outer: Mapping[str, Any]) -> bool:
    if any(inner[key] != outer[key] for key in ("query_id", "query_name", "ref_id", "ref_name")):
        return False
    qmin, qmax = sorted((inner["query_start"], inner["query_end"]))
    rmin, rmax = sorted((inner["ref_start"], inner["ref_end"]))
    oqmin, oqmax = sorted((outer["query_start"], outer["query_end"]))
    ormin, ormax = sorted((outer["ref_start"], outer["ref_end"]))
    return oqmin <= qmin <= qmax <= oqmax and ormin <= rmin <= rmax <= ormax


def _drop_overlaps(coords: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in coords:
        groups[f"{row['query_id']}{row['query_name']}-{row['ref_id']}{row['ref_name']}"].append(row)
    kept: list[dict[str, Any]] = []
    for rows in groups.values():
        for first in rows:
            if not any(
                _contained(first, second) and _tsv_key(first) != _tsv_key(second) for second in rows
            ):
                kept.append(first)
    return kept


def _parse_blast(path: Path, query_id: str, ref_id: str) -> list[dict[str, Any]]:
    coords = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            try:
                identity = float(row[2])
                if 0 <= identity <= 1.0:
                    identity = int(identity * 10000) / 100
                qstart, qend, sstart, send = (int(value) for value in row[6:10])
                evalue = float(row[10]) if len(row) >= 11 else None
            except (ValueError, IndexError) as exc:
                raise ConfigurationError(
                    f"{path.name} is not a BLAST tabular (outfmt 6) file."
                ) from exc
            coords.append(
                _coord(
                    query_id,
                    row[0],
                    qstart - 1,
                    qend,
                    ref_id,
                    row[1],
                    sstart - 1,
                    send,
                    identity,
                    evalue,
                )
            )
    return _drop_overlaps(coords)


def _parse_mummer(path: Path, query_id: str, ref_id: str, seqtype: str) -> list[dict[str, Any]]:
    coords = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            if seqtype == "nucleotide":
                if len(row) != 9:
                    raise ConfigurationError(
                        f"{path.name} is not a nucmer show-coords file (9 columns)."
                    )
            elif len(row) != 13:
                raise ConfigurationError(
                    f"{path.name} is not a promer show-coords file (13 columns)."
                )
            else:
                row = row[0:7] + row[11:13]
            try:
                rstart, rend, qstart, qend, identity = (
                    int(row[0]) - 1,
                    int(row[1]),
                    int(row[2]) - 1,
                    int(row[3]),
                    float(row[6]),
                )
            except ValueError as exc:
                raise ConfigurationError(f"{path.name} has non-numeric coordinates.") from exc
            coords.append(
                _coord(
                    query_id,
                    str(row[8]),
                    qstart,
                    qend,
                    ref_id,
                    str(row[7]),
                    rstart,
                    rend,
                    identity,
                    None,
                )
            )
    return _drop_overlaps(coords)


def _parse_pmauve(path: Path, names: Sequence[str], ref_index: int) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ConfigurationError(f"{path.name} is empty.")
        genome_count = len(header) // 2
        if len(names) < genome_count:
            raise ConfigurationError(
                f"names must list {genome_count} sequence names for this backbone file."
            )
        if not 0 <= ref_index < genome_count:
            raise ConfigurationError("ref_index is outside the genomes of the backbone file.")
        rows = []
        for raw in reader:
            try:
                row = [int(value) for value in raw]
            except ValueError as exc:
                raise ConfigurationError(f"{path.name} has non-numeric backbone values.") from exc
            position = ref_index * 2
            if row[position] < 0:
                row = [-value for value in row]
            if row.count(0) >= 2 or abs(row[position] - row[position + 1]) < 20:
                continue
            rows.append(row)
        rows.sort(key=lambda row: row[ref_index * 2])
    coords = []
    for row in rows:
        for i in range(genome_count - 1):
            idx = i * 2
            qstart, qend = row[idx], row[idx + 1]
            if qstart < 0 and qend < 0:
                qstart, qend = abs(qend), abs(qstart)
            rstart, rend = row[idx + 2], row[idx + 3]
            if rstart < 0 and rend < 0:
                rstart, rend = abs(rend), abs(rstart)
            coords.append(
                _coord(
                    names[i],
                    names[i],
                    qstart,
                    qend,
                    names[i + 1],
                    names[i + 1],
                    rstart,
                    rend,
                    None,
                    None,
                )
            )
    return _drop_overlaps(coords)


def parse_alignment_coordinates(
    alignment_file: str | Path,
    *,
    file_format: str = "blast",
    query_id: str = "query",
    ref_id: str = "ref",
    seqtype: str = "nucleotide",
    names: Sequence[str] | None = None,
    ref_index: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a genome alignment result file into alignment coordinates (``AlignCoord.parse_*_file``).

    ``file_format`` is ``blast`` (BLAST or LAST tabular, ``-outfmt 6``), ``mummer`` (``show-coords -H -T -r -k`` of nucmer,
    ``seqtype='nucleotide'``, or promer, ``'protein'``) or ``progressive_mauve`` (backbone file; ``names`` lists the genomes,
    blocks shorter than 20 bp or missing in two genomes are dropped). Completely overlapping blocks are removed.
    Coordinates are converted to 0-based starts.
    """

    path = input_path(alignment_file)
    kind = _choice(file_format, _ALIGN_FORMATS, "file_format")
    _choice(seqtype, _SEQTYPES, "seqtype")
    if kind == "blast":
        coords = _parse_blast(path, str(query_id), str(ref_id))
    elif kind == "mummer":
        coords = _parse_mummer(path, str(query_id), str(ref_id), seqtype)
    else:
        if names is None or isinstance(names, (str, bytes)):
            raise ConfigurationError(
                "names (the sequence names in file order) is required for progressive_mauve."
            )
        coords = _parse_pmauve(path, [str(item) for item in names], int(ref_index))
    return {
        "status": "completed",
        "operation": "parse_alignment_coordinates",
        "file_format": kind,
        "alignment_count": len(coords),
        "coords": coords,
        "provenance": _provenance(),
    }


def write_alignment_coordinates_tsv(
    coords: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write alignment coordinates as the 12-column pyGenomeViz TSV (``AlignCoord.write``); unknown scores are ``na``."""

    rows = _checked_coords(coords)
    target = _new_file(output_path, "output_path")
    lines = ["\t".join(_TSV_HEADER)]
    for row in rows:
        lines.append("\t".join("na" if value is None else str(value) for value in _tsv_key(row)))
    target.write_text("\n".join(lines), encoding="utf-8")
    return {
        "status": "completed",
        "operation": "write_alignment_coordinates_tsv",
        "output_path": str(target),
        "alignment_count": len(rows),
    }


def read_alignment_coordinates_tsv(
    coords_file: str | Path, *, verbose: bool = True
) -> dict[str, Any]:
    """Read alignment coordinates from the pyGenomeViz TSV written by ``write_alignment_coordinates_tsv`` (``AlignCoord.read``)."""

    path = input_path(coords_file)
    coords = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(header) != _TSV_HEADER:
            raise ConfigurationError(
                f"{path.name} does not have the alignment coordinate TSV header."
            )
        for row in reader:
            if not row:
                continue
            try:
                coords.append(
                    _coord(
                        row[0],
                        row[1],
                        int(row[2]),
                        int(row[3]),
                        row[5],
                        row[6],
                        int(row[7]),
                        int(row[8]),
                        None if row[10] == "na" else float(row[10]),
                        None if row[11] == "na" else float(row[11]),
                    )
                )
            except (ValueError, IndexError) as exc:
                raise ConfigurationError(f"{path.name} has an invalid row: {row[:3]}.") from exc
    return {
        "status": "completed",
        "operation": "read_alignment_coordinates_tsv",
        "alignment_count": len(coords),
        "coords": coords,
    }


def filter_alignment_coordinates(
    coords: Sequence[Mapping[str, Any]],
    *,
    length_thr: int | None = None,
    identity_thr: float | None = None,
    evalue_thr: float | None = None,
    remove_overlapping: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter alignment coordinates by length, identity and e-value (``AlignCoord.filter``) and optionally drop contained blocks.

    A block is dropped when either of its lengths is below ``length_thr``, its identity is below ``identity_thr`` or its
    e-value is above ``evalue_thr``; a threshold of ``None`` or 0 disables the test, and blocks without a score pass the
    score tests. ``remove_overlapping=True`` also removes blocks completely inside another block of the same
    query/reference pair (``AlignCoord.filter_overlap``).
    """

    rows = _checked_coords(coords)
    kept = []
    for row in rows:
        qlen, rlen = _lengths(row)
        if length_thr and (qlen < length_thr or rlen < length_thr):
            continue
        if identity_thr and row["identity"] and row["identity"] < identity_thr:
            continue
        if evalue_thr and row["evalue"] and row["evalue"] > evalue_thr:
            continue
        kept.append(row)
    if remove_overlapping:
        kept = _drop_overlaps(kept)
    return {
        "status": "completed",
        "operation": "filter_alignment_coordinates",
        "input_count": len(rows),
        "alignment_count": len(kept),
        "coords": kept,
    }


# Aligners -------------------------------------------------------------------------------------------


def _fasta_inputs(values: Sequence[str | Path], minimum: int = 2) -> list[Path]:
    if isinstance(values, (str, bytes)) or len(values) < minimum:
        raise ConfigurationError(f"seq_files must list at least {minimum} FASTA files.")
    paths = [input_path(value) for value in values]
    names = [_stem(path) for path in paths]
    if len(set(names)) != len(names):
        raise ConfigurationError("FASTA file names (without suffix) must be unique.")
    for name in names:
        if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
            raise ConfigurationError(
                f"file name {name!r} must use only letters, digits, '.', '_' and '-'."
            )
    return paths


def _located(directory: str | Path | None, name: str, hint: str) -> str:
    explicit = None if directory is None else str(Path(directory).expanduser() / name)
    return _executable(name, explicit, hint)


def _failure(steps: Any, backend: str) -> dict[str, Any]:
    return {
        "status": "failed", "backend": backend, "commands": steps.commands, "return_codes": steps.return_codes,
        "outputs": [], "alignment_count": 0, "coords": [], "log_tail": steps.tail(),
    }  # fmt: skip


def align_genomes_blast(
    seq_files: Sequence[str | Path],
    output_dir: str | Path,
    *,
    seqtype: str = "nucleotide",
    evalue: float = 1e-3,
    threads: int = 1,
    timeout_seconds: float = 3600.0,
    blast_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align each genome FASTA with the next one by BLAST (``Blast.run``, ``pgv-blast``).

    ``seqtype='nucleotide'`` runs ``blastn`` and ``'protein'`` runs ``tblastx`` (translated comparison). For every
    neighbouring pair the second file becomes a ``makeblastdb`` database, the first is the query, and the tabular
    result (``-outfmt 6 -evalue <evalue>``) is parsed into alignment coordinates with contained hits removed.
    """

    files = _fasta_inputs(seq_files)
    _choice(seqtype, _SEQTYPES, "seqtype")
    if isinstance(evalue, bool) or not isinstance(evalue, (int, float)) or evalue <= 0:
        raise ConfigurationError("evalue must be a positive number.")
    _positive(threads, "threads")
    folder = _check_new_folder(str(output_dir))
    tool = "blastn" if seqtype == "nucleotide" else "tblastx"
    programs = {name: _located(blast_dir, name, _BLAST_HINT) for name in ("makeblastdb", tool)}
    folder.mkdir(parents=True, exist_ok=True)
    steps = _Steps("blast")
    coords: list[dict[str, Any]] = []
    outputs: list[str] = []
    for index in range(len(files) - 1):
        query, ref = files[index], files[index + 1]
        database = folder / f"{_stem(ref)}_blastdb"
        code, _ = steps.run(
            programs["makeblastdb"], ["-in", str(ref), "-dbtype", "nucl", "-out", str(database)],
            cwd=folder, timeout_seconds=timeout_seconds,
        )  # fmt: skip
        if code != 0:
            return _failure(steps, "blast")
        result = folder / f"{index + 1:02d}_{_stem(query)}_vs_{_stem(ref)}.tsv"
        code, _ = steps.run(
            programs[tool],
            ["-query", str(query), "-db", str(database), "-out", str(result), "-outfmt", "6", "-evalue", str(evalue), "-num_threads", str(threads)],
            cwd=folder, timeout_seconds=timeout_seconds,
        )  # fmt: skip
        if code != 0 or not result.is_file():
            return _failure(steps, "blast")
        outputs.append(str(result))
        coords.extend(_parse_blast(result, _stem(query), _stem(ref)))
    return {
        "status": "completed", "backend": "blast", "tool": tool, "commands": steps.commands,
        "return_codes": steps.return_codes, "outputs": outputs, "alignment_count": len(coords),
        "coords": coords, "log_tail": steps.tail(),
    }  # fmt: skip


def align_genomes_last(
    seq_files: Sequence[str | Path],
    output_dir: str | Path,
    *,
    threads: int = 1,
    timeout_seconds: float = 3600.0,
    last_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align each genome FASTA with the next one by LAST (``Last.run``), following the LAST cookbook genome recipe.

    For every neighbouring pair: ``lastdb`` on the second file, ``lastal -D 1e9 --split-f=MAF+`` for the first, ``last-split -r``
    to make the alignment one-to-one, and ``maf-convert -n blasttab``; the tabular result is parsed into alignment coordinates.
    """

    files = _fasta_inputs(seq_files)
    _positive(threads, "threads")
    folder = _check_new_folder(str(output_dir))
    programs = {
        name: _located(last_dir, name, _LAST_HINT)
        for name in ("lastdb", "lastal", "last-split", "maf-convert")
    }
    folder.mkdir(parents=True, exist_ok=True)
    steps = _Steps("last")
    coords: list[dict[str, Any]] = []
    outputs: list[str] = []
    for index in range(len(files) - 1):
        query, ref = files[index], files[index + 1]
        prefix = f"{index + 1:02d}_{_stem(query)}_vs_{_stem(ref)}"
        database = folder / f"{_stem(ref)}_lastdb"
        code, _ = steps.run(
            programs["lastdb"],
            [str(database), str(ref), "-P", str(threads)],
            cwd=folder,
            timeout_seconds=timeout_seconds,
        )
        if code != 0:
            return _failure(steps, "last")
        stages = [
            (
                "lastal",
                [str(database), str(query), "-P", str(threads), "-D", "1e9", "--split-f=MAF+"],
                f"{prefix}_many-to-one.maf",
            ),
            (
                "last-split",
                ["-r", str(folder / f"{prefix}_many-to-one.maf")],
                f"{prefix}_one-to-one.maf",
            ),
            (
                "maf-convert",
                ["-n", "blasttab", str(folder / f"{prefix}_one-to-one.maf")],
                f"{prefix}.tsv",
            ),
        ]
        for program, arguments, filename in stages:
            code, text = steps.run(
                programs[program],
                arguments,
                cwd=folder,
                timeout_seconds=timeout_seconds,
                max_output_bytes=_LARGE_OUTPUT,
            )
            if code != 0:
                return _failure(steps, "last")
            (folder / filename).write_text(text, encoding="utf-8")
            outputs.append(str(folder / filename))
        coords.extend(_parse_blast(folder / f"{prefix}.tsv", _stem(query), _stem(ref)))
    return {
        "status": "completed", "backend": "last", "commands": steps.commands, "return_codes": steps.return_codes,
        "outputs": outputs, "alignment_count": len(coords), "coords": coords, "log_tail": steps.tail(),
    }  # fmt: skip


# Drawing --------------------------------------------------------------------------------------------


def _feature_kws(item: Mapping[str, Any]) -> dict[str, Any]:
    kws: dict[str, Any] = {}
    for key, kind in (
        ("label", str),
        ("facecolor", str),
        ("edgecolor", str),
        ("plotstyle", str),
        ("linewidth", float),
    ):
        if key in item:
            kws[key] = kind(item[key])
    return kws


def draw_genome_segment_features(
    features: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    segment_name: str = "genome",
    segment_size: int,
    fig_width: float = 15.0,
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw gene glyphs on one genome segment with pyGenomeViz (``add_feature(s)``, ``add_exon_feature(s)``) and save the figure.

    Each feature is a mapping. A gene arrow needs ``start``, ``end`` and ``strand`` (1 or -1); a gene with exons
    gives ``parts`` (a list of ``[start, end]``) instead of ``start``/``end``, and is drawn as exon boxes joined by
    intron lines. Optional keys are ``label``, ``facecolor``, ``edgecolor``, ``plotstyle`` (``bigarrow``, ``arrow``,
    ``box``, ``bigbox``, ``rbox``, ``bigrbox``) and ``linewidth``. Needs pyGenomeViz 0.4.x and matplotlib.
    The promoter glyph of pyGenomeViz 1.x is not available.
    """

    target = _new_file(output_path, "output_path", (".png", ".svg", ".pdf"))
    size = _positive(segment_size, "segment_size")
    if isinstance(features, (str, bytes)) or not features:
        raise ConfigurationError("features must be a non-empty list of mappings.")
    if isinstance(fig_width, bool) or not isinstance(fig_width, (int, float)) or fig_width <= 0:
        raise ConfigurationError("fig_width must be a positive number.")
    _positive(dpi, "dpi")
    pgv = _require_isolated("pygenomeviz", _EXTRA)
    plt = require("matplotlib.pyplot", _EXTRA)
    graph = pgv.GenomeViz(
        fig_width=fig_width, fig_track_height=1.0, feature_track_ratio=1.0, plot_size_thr=0
    )
    track = graph.add_feature_track(str(segment_name), size=size)
    drawn = []
    for item in features:
        if not isinstance(item, Mapping):
            raise ConfigurationError("every feature must be a mapping.")
        strand = int(item.get("strand", 1))
        if strand not in (1, -1):
            raise ConfigurationError("feature strand must be 1 or -1.")
        try:
            if "parts" in item:
                parts = [(int(a), int(b)) for a, b in item["parts"]]
                track.add_exon_feature(parts, strand, **_feature_kws(item))
                drawn.append({"parts": [list(p) for p in parts], "strand": strand})
            else:
                start, end = int(item["start"]), int(item["end"])
                track.add_feature(start, end, strand, **_feature_kws(item))
                drawn.append({"start": start, "end": end, "strand": strand})
        except (KeyError, TypeError, ValueError) as exc:
            raise ConfigurationError(f"a feature is invalid or outside the segment: {exc}") from exc
    figure = graph.plotfig(dpi=dpi)
    figure.savefig(target, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return {
        "status": "completed",
        "operation": "draw_genome_segment_features",
        "image_path": str(target),
        "segment": {"name": segment_name, "size": size},
        "features": drawn,
    }


def export_genome_html_viewer(
    gbk_resources: Sequence[str | Path],
    html_path: str | Path,
    *,
    fig_width: float = 15.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export GenBank genomes (CDS features as tracks) as one interactive HTML viewer with hover tooltips (``GenomeViz.savefig_html``).

    Needs pyGenomeViz 0.4.x. The output is a single HTML file that embeds the SVG figure and the viewer script.
    """

    from dnakit.visualization.pygenomeviz_gap_functions import (
        _draw_cds,
        _genome_summary,
        _genome_viz,
        _genomes,
    )

    genomes = _genomes([str(item) for item in gbk_resources], minimum=1)
    target = _new_file(html_path, "html_path", (".html", ".htm"))
    pgv = _require_isolated("pygenomeviz", _EXTRA)
    graph = _genome_viz(
        pgv,
        fig_width=fig_width,
        fig_track_height=1.0,
        align_type="left",
        tick_style=None,
        link_track_ratio=5.0,
    )
    for gbk in genomes:
        track = graph.add_feature_track(
            gbk.name, size=gbk.range_size, start_pos=gbk.min_range, labelsize=20
        )
        _draw_cds(track, gbk, plotstyle="bigarrow", facecolor="orange", linewidth=0.0)
    graph.savefig_html(target)
    return {
        "status": "completed",
        "operation": "export_genome_html_viewer",
        "html_path": str(target),
        "bytes": target.stat().st_size,
        "genomes": [_genome_summary(gbk) for gbk in genomes],
    }


__all__ = [
    "align_genomes_blast",
    "align_genomes_last",
    "calc_genbank_gc_content",
    "draw_genome_segment_features",
    "export_genome_html_viewer",
    "extract_genome_features",
    "extract_gff_exon_features",
    "filter_alignment_coordinates",
    "list_contig_features",
    "parse_alignment_coordinates",
    "read_alignment_coordinates_tsv",
    "read_contig_sequences",
    "read_gff_contig_sizes",
    "select_features_within_range",
    "write_alignment_coordinates_tsv",
    "write_genbank_cds_fasta",
    "write_genome_fasta",
]

_register_public_verbose(globals(), tuple(__all__))
