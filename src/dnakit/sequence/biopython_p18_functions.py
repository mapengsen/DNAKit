"""Biopython sequence-file and sequence-string functions of integration pack P18 (sequence area).

The functions cover the sequence formats of ``Bio.SeqIO`` that DNAKit did not read or write
before (QUAL, IMGT, GenBank and EMBL CDS records, GCK, GFA1/GFA2, IntelliGenetics, nib, PIR,
SeqXML, tab, Xdna, PHD and SFF), exact-match operations of ``Bio.Seq.Seq``, INSDC location
strings, EMBOSS ``primersearch`` input, GO-annotation filtering, substitution-matrix text,
NCBI Entrez XML files and BGZF virtual offsets.

Biopython is imported only when a function runs. Readers read local files; writers never
overwrite (the output path must be new). Sequence-string operations keep the case and the
characters of the input, as ``Bio.Seq`` does. Nothing here uses the network.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _biopython_p18_common as _c
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = _c.EXTRA
_READ_FORMATS = (
    "abi",
    "abi-trim",
    "qual",
    "imgt",
    "genbank-cds",
    "embl-cds",
    "gck",
    "gfa1",
    "gfa2",
    "ig",
    "nib",
    "pir",
    "seqxml",
    "tab",
    "xdna",
)
_WRITE_FORMATS = ("imgt", "nib", "phd", "pir", "seqxml", "tab", "xdna")
_CONVERT_SOURCE_FORMATS = (
    "abi",
    "abi-trim",
    "embl",
    "embl-cds",
    "fasta",
    "fastq",
    "gck",
    "genbank",
    "genbank-cds",
    "gfa1",
    "gfa2",
    "ig",
    "imgt",
    "nib",
    "pir",
    "qual",
    "seqxml",
    "sff",
    "sff-trim",
    "tab",
    "xdna",
)
_CONVERT_TARGET_FORMATS = (
    "embl",
    "fasta",
    "fastq",
    "genbank",
    "imgt",
    "nib",
    "phd",
    "pir",
    "qual",
    "seqxml",
    "sff",
    "tab",
    "xdna",
)
_BINARY_TARGETS = frozenset({"nib", "seqxml", "sff", "xdna"})
_SINGLE_RECORD_TARGETS = frozenset({"nib", "xdna"})
_TRIM_MODES = ("strip", "lstrip", "rstrip", "removeprefix", "removesuffix")
_SEQUENCE_SYMBOLS = re.compile(r"^[A-Za-z*.\-]+$")
_PRIMER = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$")
_DOCTYPE = re.compile(r'<!DOCTYPE\s+\S+\s+(?:PUBLIC\s+"[^"]*"\s+|SYSTEM\s+)"([^"]+)"')


def _record_summary(record: Any, include_sequence: bool, include_raw: bool) -> dict[str, Any]:
    try:
        sequence: str | None = str(record.seq) if include_sequence else None
    except Exception:  # undefined sequence content (for example QUAL files)
        sequence = None
    return {
        "id": str(record.id),
        "name": str(record.name),
        "description": str(record.description),
        "length": len(record),
        "sequence": sequence,
        "annotations": _c.plain(
            {k: v for k, v in record.annotations.items() if include_raw or k != "abif_raw"}
        ),
        "letter_annotations": {
            str(key): [_c.plain(item) for item in value]
            for key, value in record.letter_annotations.items()
        },
        "feature_count": len(record.features),
        "features": [
            {
                "type": str(feature.type),
                "location": str(feature.location),
                "qualifiers": _c.plain(feature.qualifiers),
            }
            for feature in record.features[:1000]
        ],
        "dbxrefs": [str(item) for item in record.dbxrefs],
    }


def read_sequence_file(
    path: str | Path,
    *,
    file_format: str = "tab",
    max_records: int = 10_000,
    include_sequence: bool = True,
    include_raw_trace: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read records of a less common ``Bio.SeqIO`` format and return them as plain dictionaries.

    ``file_format`` is one of ``qual`` (PHRED qualities only, the sequence is undefined and
    returned as null), ``imgt`` (IMGT/HLA EMBL-style flat file), ``genbank-cds`` and
    ``embl-cds`` (one record per CDS feature, holding the translation), ``gck`` (Gene
    Construction Kit), ``gfa1`` and ``gfa2`` (assembly-graph segments), ``ig``
    (IntelliGenetics), ``nib`` (UCSC four-bit), ``pir`` (NBRF/PIR), ``seqxml`` and ``xdna``,
    or ``tab`` (two tab-separated columns, identifier and sequence). Each record has its id,
    name, description, length, sequence, annotations, per-letter annotations (for QUAL the
    ``phred_quality`` list), features with locations and qualifiers, and cross-references.
    ``abi`` reads an ABI Sanger trace (called bases and PHRED qualities; the raw ABIF data
    block is left out unless ``include_raw_trace``), and ``abi-trim`` the same record after
    Biopython's Mott-style trimming with error cutoff 0.05; the chromatogram channels are
    available from ``dnakit.bioinformatics.cloning.read_sanger``. The result stops with an
    error when the file holds more than ``max_records`` records.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _READ_FORMATS, "file_format")
    limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    with_sequence = _c.flag(include_sequence, "include_sequence")
    raw_trace = _c.flag(include_raw_trace, "include_raw_trace")
    seq_io = require("Bio.SeqIO", _EXTRA)
    records: list[dict[str, Any]] = []
    try:
        for record in seq_io.parse(str(source), fmt):
            if len(records) >= limit:
                raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
            records.append(_record_summary(record, with_sequence, raw_trace))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_sequence_file",
        "source": str(source),
        "file_format": fmt,
        "record_count": len(records),
        "records": records,
        "provenance": _c.provenance(
            parser="Bio.SeqIO.parse",
            coordinates="feature locations are written as Biopython prints them (0-based start)",
        ),
    }


def _build_record(item: object, index: int, fmt: str) -> Any:
    seq_module = require("Bio.Seq", _EXTRA)
    record_module = require("Bio.SeqRecord", _EXTRA)
    if not isinstance(item, Mapping):
        raise ConfigurationError(f"records[{index}] must be a mapping.")
    identifier = _c.identifier(item.get("id"), f"records[{index}].id")
    sequence = _c.text(item.get("sequence"), f"records[{index}].sequence", limit=_c.MAX_SEQUENCE)
    if not _SEQUENCE_SYMBOLS.match(sequence):
        raise ConfigurationError(f"records[{index}].sequence has symbols outside letters.")
    description = item.get("description", "")
    name = item.get("name", identifier)
    if not isinstance(description, str) or not isinstance(name, str):
        raise ConfigurationError(f"records[{index}] description and name must be strings.")
    record = record_module.SeqRecord(
        seq_module.Seq(sequence), id=identifier, name=name or identifier, description=description
    )
    record.annotations["molecule_type"] = _c.choice(
        item.get("molecule_type", "DNA"), ("DNA", "RNA"), f"records[{index}].molecule_type"
    )
    topology = item.get("topology")
    if topology is not None:
        record.annotations["topology"] = _c.choice(
            topology, ("linear", "circular"), f"records[{index}].topology"
        )
    qualities = item.get("qualities")
    if qualities is not None:
        if not isinstance(qualities, Sequence) or isinstance(qualities, str):
            raise ConfigurationError(f"records[{index}].qualities must be a list of integers.")
        if len(qualities) != len(sequence):
            raise ConfigurationError(f"records[{index}] needs one quality per base.")
        record.letter_annotations["phred_quality"] = [
            _c.integer(value, f"records[{index}].qualities", minimum=0, maximum=93)
            for value in qualities
        ]
    elif fmt == "phd":
        raise ConfigurationError(f"records[{index}] needs qualities for the phd format.")
    return record


def write_sequence_file(
    records: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    file_format: str = "tab",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write sequence records to a new file in an IMGT, nib, PHD, PIR, SeqXML, tab or Xdna format.

    Each record is a mapping with ``id`` and ``sequence`` and, optionally, ``name``,
    ``description``, ``molecule_type`` (``DNA`` or ``RNA``), ``topology`` (``linear`` or
    ``circular``, used by Xdna and IMGT) and ``qualities`` (PHRED integers, required for
    ``phd``). The ``nib`` and ``xdna`` formats hold exactly one record. SFF files carry flow
    data, so they are written by ``convert_sequence_format`` from an SFF source instead.
    """

    fmt = _c.choice(file_format, _WRITE_FORMATS, "file_format")
    candidate: object = records
    if not isinstance(candidate, Sequence) or isinstance(candidate, str) or not candidate:
        raise ConfigurationError("records must be a non-empty list of mappings.")
    if len(records) > _c.MAX_RECORDS:
        raise ConfigurationError(f"records may hold at most {_c.MAX_RECORDS} entries.")
    if fmt in _SINGLE_RECORD_TARGETS and len(records) != 1:
        raise ConfigurationError(f"The {fmt} format holds exactly one record.")
    destination = _c.new_output(output_path)
    seq_io = require("Bio.SeqIO", _EXTRA)
    built = [_build_record(item, index, fmt) for index, item in enumerate(records)]
    binary = fmt in _BINARY_TARGETS
    try:
        with _c.open_new(destination, binary=binary) as handle:
            written = seq_io.write(built, handle, fmt)
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(f"Biopython could not write {fmt}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_sequence_file",
        "output_path": str(destination),
        "file_format": fmt,
        "record_count": int(written),
        "bytes_written": destination.stat().st_size,
        "provenance": _c.provenance(writer="Bio.SeqIO.write"),
    }


def convert_sequence_format(
    source: str | Path,
    output_path: str | Path,
    *,
    source_format: str = "genbank",
    target_format: str = "tab",
    max_records: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sequence file between ``Bio.SeqIO`` formats, including SFF output.

    ``source_format`` is any of the formats listed in ``read_sequence_file`` plus abi,
    abi-trim, embl, fasta, fastq, genbank, sff and sff-trim. ``target_format`` is one of embl,
    fasta, fastq, genbank, imgt, nib, phd, pir, qual, seqxml, sff, tab and xdna. Information
    that a target format cannot hold (features in FASTA, qualities in tab) is dropped, as
    Biopython does; the ``nib`` and ``xdna`` targets accept exactly one record. When the source
    is an SFF file with a Roche XML manifest, the manifest is copied into an SFF target. The
    output file must not exist.
    """

    path = _c.source_file(source, "source")
    in_format = _c.choice(source_format, _CONVERT_SOURCE_FORMATS, "source_format")
    out_format = _c.choice(target_format, _CONVERT_TARGET_FORMATS, "target_format")
    limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    destination = _c.new_output(output_path)
    seq_io = require("Bio.SeqIO", _EXTRA)
    sff_io = require("Bio.SeqIO.SffIO", _EXTRA)
    manifest: str | None = None
    if out_format == "sff" and in_format in ("sff", "sff-trim"):
        try:
            with path.open("rb") as handle:
                manifest = sff_io.ReadRocheXmlManifest(handle)
        except ValueError:
            manifest = None
    binary = out_format in _BINARY_TARGETS
    count = 0
    try:
        records: list[Any] = []
        for record in seq_io.parse(str(path), in_format):
            if len(records) >= limit:
                raise ConfigurationError(f"{path.name} exceeds max_records={limit}.")
            records.append(record)
        with _c.open_new(destination, binary=binary) as handle:
            if out_format == "sff":
                count = sff_io.SffWriter(handle, xml=manifest).write_file(records)
            else:
                count = seq_io.write(records, handle, out_format)
    except ConfigurationError:
        destination.unlink(missing_ok=True)
        raise
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(
            f"Biopython could not convert {path.name} from {in_format} to {out_format}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "convert_sequence_format",
        "source": str(path),
        "output_path": str(destination),
        "source_format": in_format,
        "target_format": out_format,
        "record_count": int(count),
        "roche_manifest_copied": manifest is not None,
        "provenance": _c.provenance(
            reader="Bio.SeqIO.parse",
            writer="Bio.SeqIO.SffIO.SffWriter" if out_format == "sff" else "Bio.SeqIO.write",
        ),
    }


def read_sff_manifest(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Return the Roche XML manifest stored in an SFF file (``ReadRocheXmlManifest``).

    454 SFF files may carry an XML manifest describing the run. The result holds the XML text
    and its length, or ``has_manifest = False`` with a null ``xml`` when the file has none.
    """

    source = _c.source_file(path)
    sff_io = require("Bio.SeqIO.SffIO", _EXTRA)
    try:
        with source.open("rb") as handle:
            xml: str | None = sff_io.ReadRocheXmlManifest(handle)
    except ValueError as exc:
        if "No XML manifest" in str(exc):
            xml = None
        else:
            raise ConfigurationError(f"{source.name} is not a valid SFF file: {exc}") from exc
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid SFF file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_sff_manifest",
        "source": str(source),
        "has_manifest": xml is not None,
        "xml": xml,
        "xml_length": 0 if xml is None else len(xml),
        "provenance": _c.provenance(parser="Bio.SeqIO.SffIO.ReadRocheXmlManifest"),
    }


def _plain_sequence(value: object, name: str) -> str:
    return _c.text(value, name, limit=_c.MAX_SEQUENCE)


def count_subsequence(
    sequence: str,
    subsequence: str,
    *,
    overlap: bool = False,
    start: int | None = None,
    end: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count occurrences of an exact subsequence, non-overlapping or overlapping.

    ``overlap=False`` is ``Seq.count`` (non-overlapping, scanning left to right) and
    ``overlap=True`` is ``Seq.count_overlap``. The match is case sensitive. ``start`` and
    ``end`` are optional 0-based slice limits, as in Python slicing.
    """

    seq = _plain_sequence(sequence, "sequence")
    sub = _plain_sequence(subsequence, "subsequence")
    use_overlap = _c.flag(overlap, "overlap")
    first = _c.optional_integer(start, "start", minimum=-_c.MAX_SEQUENCE, maximum=_c.MAX_SEQUENCE)
    last = _c.optional_integer(end, "end", minimum=-_c.MAX_SEQUENCE, maximum=_c.MAX_SEQUENCE)
    seq_module = require("Bio.Seq", _EXTRA)
    biopython_sequence = seq_module.Seq(seq)
    begin = 0 if first is None else first
    finish = len(seq) if last is None else last
    method = biopython_sequence.count_overlap if use_overlap else biopython_sequence.count
    return {
        "status": "completed",
        "operation": "count_subsequence",
        "subsequence": sub,
        "overlap": use_overlap,
        "sequence_length": len(seq),
        "count": int(method(sub, begin, finish)),
        "provenance": _c.provenance(
            method="Bio.Seq.Seq.count_overlap" if use_overlap else "Bio.Seq.Seq.count"
        ),
    }


def search_subsequences(
    sequence: str, subsequences: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Find every position of several exact subsequences (``Seq.search``).

    Matches may overlap each other; positions are 0-based and listed in increasing order.
    Each hit has the position and the matched subsequence. The match is case sensitive.
    """

    seq = _plain_sequence(sequence, "sequence")
    if not isinstance(subsequences, Sequence) or isinstance(subsequences, str) or not subsequences:
        raise ConfigurationError("subsequences must be a non-empty list of strings.")
    if len(subsequences) > 1000:
        raise ConfigurationError("subsequences may hold at most 1000 entries.")
    words = [_plain_sequence(item, "subsequence") for item in subsequences]
    seq_module = require("Bio.Seq", _EXTRA)
    hits = [
        {"position": int(position), "subsequence": str(word)}
        for position, word in seq_module.Seq(seq).search(words)
    ]
    return {
        "status": "completed",
        "operation": "search_subsequences",
        "sequence_length": len(seq),
        "hit_count": len(hits),
        "hits": hits,
        "provenance": _c.provenance(method="Bio.Seq.Seq.search"),
    }


def split_sequence(
    sequence: str,
    separator: str | None = None,
    *,
    maxsplit: int = -1,
    from_right: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a sequence at a separator (``Seq.split``, or ``Seq.rsplit`` with ``from_right``).

    Without a separator the split is on whitespace, as for Python strings. ``maxsplit`` limits
    the number of splits; -1 means no limit. The pieces are returned in order.
    """

    seq = _plain_sequence(sequence, "sequence")
    sep = None if separator is None else _plain_sequence(separator, "separator")
    limit = _c.integer(maxsplit, "maxsplit", minimum=-1, maximum=_c.MAX_SEQUENCE)
    reverse = _c.flag(from_right, "from_right")
    seq_module = require("Bio.Seq", _EXTRA)
    biopython_sequence = seq_module.Seq(seq)
    parts = (
        biopython_sequence.rsplit(sep, limit) if reverse else biopython_sequence.split(sep, limit)
    )
    pieces = [str(part) for part in parts]
    return {
        "status": "completed",
        "operation": "split_sequence",
        "separator": sep,
        "from_right": reverse,
        "piece_count": len(pieces),
        "pieces": pieces,
        "provenance": _c.provenance(
            method="Bio.Seq.Seq.rsplit" if reverse else "Bio.Seq.Seq.split"
        ),
    }


def trim_sequence(
    sequence: str,
    *,
    mode: str = "strip",
    pattern: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim the ends of a sequence the way ``Bio.Seq.Seq`` does.

    ``mode`` is ``strip`` (both ends), ``lstrip`` or ``rstrip``: ``pattern`` is then a set of
    characters to remove (whitespace when null). ``removeprefix`` and ``removesuffix`` remove
    the exact ``pattern`` string once when the sequence starts or ends with it. Case is kept.
    """

    seq = _plain_sequence(sequence, "sequence")
    selected = _c.choice(mode, _TRIM_MODES, "mode")
    chars = None if pattern is None else _plain_sequence(pattern, "pattern")
    if selected in ("removeprefix", "removesuffix") and chars is None:
        raise ConfigurationError(f"{selected} needs a pattern.")
    seq_module = require("Bio.Seq", _EXTRA)
    trimmed = str(getattr(seq_module.Seq(seq), selected)(chars))
    return {
        "status": "completed",
        "operation": "trim_sequence",
        "mode": selected,
        "pattern": chars,
        "original_length": len(seq),
        "trimmed": trimmed,
        "trimmed_length": len(trimmed),
        "provenance": _c.provenance(method=f"Bio.Seq.Seq.{selected}"),
    }


def _position_summary(position: Any) -> dict[str, Any]:
    return {"type": type(position).__name__, "text": str(position), "value": int(position)}


def parse_feature_location(location: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse an INSDC feature-location string such as ``complement(join(1..10,20..30))``.

    Uses ``Bio.SeqFeature.Location.fromstring``. Returned are the overall start (0-based) and
    end, the strand (1, -1, 0 or null), the operator (``join``/``order``/``single``), the
    reference accession when the location names one, and every part with its positions. Fuzzy
    positions (``<1``, ``>200``, ``(5.9)``, ``12^13``) are reported with their Biopython type.
    """

    value = _c.text(location, "location", limit=100_000).strip()
    feature = require("Bio.SeqFeature", _EXTRA)
    try:
        parsed = feature.Location.fromstring(value)
    except Exception as exc:
        raise ConfigurationError(f"{value!r} is not a valid INSDC location: {exc}") from exc
    parts = []
    for part in parsed.parts:
        parts.append(
            {
                "start": _position_summary(part.start),
                "end": _position_summary(part.end),
                "strand": part.strand,
                "ref": part.ref,
            }
        )
    return {
        "status": "completed",
        "operation": "parse_feature_location",
        "location": value,
        "location_type": type(parsed).__name__,
        "start": _position_summary(parsed.start),
        "end": _position_summary(parsed.end),
        "strand": parsed.strand,
        "operator": getattr(parsed, "operator", "single"),
        "ref": parsed.ref,
        "part_count": len(parts),
        "parts": parts,
        "normalized": str(parsed),
        "provenance": _c.provenance(
            method="Bio.SeqFeature.Location.fromstring",
            coordinates="start is 0-based and end is exclusive; the input string is 1-based",
        ),
    }


def format_primersearch_input(
    primer_sets: Mapping[str, Sequence[str]],
    *,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the input file of EMBOSS ``primersearch`` from named forward and reverse primers.

    ``primer_sets`` maps a set name to ``[forward, reverse]``. The text has one line per set:
    name, forward primer, reverse primer (``InputRecord.add_primer_set``). When ``output_path``
    is given the text is also written to that new file.
    """

    if not isinstance(primer_sets, Mapping) or not primer_sets:
        raise ConfigurationError("primer_sets must be a non-empty mapping.")
    if len(primer_sets) > 10_000:
        raise ConfigurationError("primer_sets may hold at most 10000 entries.")
    destination = None if output_path is None else _c.new_output(output_path)
    primer_search = require("Bio.Emboss.PrimerSearch", _EXTRA)
    record = primer_search.InputRecord()
    for name, pair in primer_sets.items():
        _c.identifier(name, "primer set name")
        if not isinstance(pair, Sequence) or isinstance(pair, str) or len(pair) != 2:
            raise ConfigurationError(f"primer set {name!r} needs [forward, reverse].")
        for primer in pair:
            if not isinstance(primer, str) or not _PRIMER.match(primer):
                raise ConfigurationError(f"primer set {name!r} has a non-IUPAC primer.")
        record.add_primer_set(name, pair[0], pair[1])
    content = str(record)
    if destination is not None:
        with destination.open("x", encoding="utf-8") as handle:
            handle.write(content)
    return {
        "status": "completed",
        "operation": "format_primersearch_input",
        "primer_set_count": len(primer_sets),
        "text": content,
        "output_path": None if destination is None else str(destination),
        "provenance": _c.provenance(method="Bio.Emboss.PrimerSearch.InputRecord"),
    }


def filter_go_annotations(
    path: str | Path,
    field_values: Mapping[str, Sequence[str]],
    *,
    max_records: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the records of a GAF annotation file whose fields match given values.

    ``field_values`` maps GAF field names (for example ``DB_Object_Symbol``, ``GO_ID``,
    ``Evidence``, ``Aspect``, ``Taxon_ID``) to accepted values. A record is kept when any one of
    the named fields holds any accepted value (``Bio.UniProt.GOA.record_has``). List-valued
    fields such as ``DB:Reference`` match when any element does. The count of kept records is
    bounded by ``max_records``.
    """

    source = _c.source_file(path)
    limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    if not isinstance(field_values, Mapping) or not field_values:
        raise ConfigurationError("field_values must be a non-empty mapping.")
    wanted: dict[str, set[str]] = {}
    for key, values in field_values.items():
        _c.identifier(key, "field name")
        if isinstance(values, str) or not isinstance(values, Sequence) or not values:
            raise ConfigurationError(f"field_values[{key!r}] must be a non-empty list of strings.")
        wanted[key] = {_c.text(item, f"field_values[{key!r}] value") for item in values}
    goa = require("Bio.UniProt.GOA", _EXTRA)
    kept: list[dict[str, Any]] = []
    total = 0
    try:
        with source.open(encoding="utf-8") as handle:
            for record in goa.gafiterator(handle):
                total += 1
                for key in wanted:
                    if key not in record:
                        raise ConfigurationError(f"{key!r} is not a field of this GAF file.")
                if goa.record_has(record, wanted):
                    if len(kept) >= limit:
                        raise ConfigurationError(f"{source.name} keeps more than {limit} records.")
                    kept.append(_c.plain(dict(record)))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid GAF file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "filter_go_annotations",
        "source": str(source),
        "records_read": total,
        "kept_count": len(kept),
        "records": kept,
        "provenance": _c.provenance(
            parser="Bio.UniProt.GOA.gafiterator", filter="Bio.UniProt.GOA.record_has"
        ),
    }


def format_substitution_matrix(
    name: str = "NUC.4.4", *, fmt: str = "%i", verbose: bool = True
) -> dict[str, Any]:
    """Return a bundled substitution matrix (for example ``NUC.4.4``) as formatted text.

    ``fmt`` is the printf-style number format used by ``Array.format`` (``%i``, ``%.1f`` ...).
    The text starts with the matrix header comments, then the alphabet and one row per letter.
    """

    selected = _c.text(name, "name", limit=100).strip()
    number_format = _c.text(fmt, "fmt", limit=20)
    if not re.fullmatch(r"%[-+ 0#]*\d*(?:\.\d+)?[diouxXeEfFgG]", number_format):
        raise ConfigurationError("fmt must be one printf-style number format such as %i or %.1f.")
    matrices = require("Bio.Align.substitution_matrices", _EXTRA)
    available = sorted(matrices.load())
    if selected not in available:
        raise ConfigurationError(
            f"unknown substitution matrix {selected!r}; choose from {available}."
        )
    matrix = matrices.load(selected)
    content = str(matrix.format(number_format))
    return {
        "status": "completed",
        "operation": "format_substitution_matrix",
        "name": selected,
        "alphabet": str(matrix.alphabet),
        "fmt": number_format,
        "text": content,
        "provenance": _c.provenance(method="Bio.Align.substitution_matrices.Array.format"),
    }


def _dtd_available(head: str, handler: Any) -> tuple[bool, str | None]:
    match = _DOCTYPE.search(head)
    if match is None:
        return True, None
    name = match.group(1).rsplit("/", 1)[-1]
    for folder in (
        handler.local_dtd_dir,
        handler.global_dtd_dir,
        handler.local_xsd_dir,
        handler.global_xsd_dir,
    ):
        if folder and (Path(folder) / name).is_file():
            return True, name
    return False, name


def parse_entrez_xml(
    path: str | Path,
    *,
    multiple: bool = False,
    max_records: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a saved NCBI Entrez (E-utilities) XML result file into plain Python data.

    ``multiple=False`` is ``Bio.Entrez.read`` (one result, for example ESearch, EInfo, ELink);
    ``multiple=True`` is ``Bio.Entrez.parse`` (a list of records, for example EFetch or
    ESummary sets). Biopython needs the DTD named in the file; this function does not download
    it, so it stops with an error when the DTD is neither bundled with Biopython nor in
    the user DTD directory. Validation of unknown tags is on, as in Biopython.
    """

    source = _c.source_file(path)
    many = _c.flag(multiple, "multiple")
    limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    entrez = require("Bio.Entrez", _EXTRA)
    parser = require("Bio.Entrez.Parser", _EXTRA)
    with source.open("rb") as handle:
        head = handle.read(4096).decode("utf-8", errors="replace")
    present, dtd_name = _dtd_available(head, parser.DataHandler)
    if not present:
        raise ConfigurationError(
            f"The DTD {dtd_name} named in {source.name} is not available locally; this function "
            "does not download it. Place it in the Biopython DTD directory first."
        )
    try:
        with source.open("rb") as handle:
            if many:
                data: Any = []
                for item in entrez.parse(handle):
                    if len(data) >= limit:
                        raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
                    data.append(_c.plain(item))
            else:
                data = _c.plain(entrez.read(handle))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not valid Entrez XML: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parse_entrez_xml",
        "source": str(source),
        "multiple": many,
        "dtd": dtd_name,
        "record_count": len(data) if many else 1,
        "data": data,
        "provenance": _c.provenance(
            parser="Bio.Entrez.parse" if many else "Bio.Entrez.read",
            network="none; the DTD must already be installed",
        ),
    }


def bgzf_virtual_offset(
    block_start: int | None = None,
    within_block: int | None = None,
    *,
    virtual_offset: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine or split a BGZF virtual file offset.

    Give ``block_start`` (compressed offset of the BGZF block, below 2**48) and
    ``within_block`` (offset in the uncompressed block, below 2**16) to compute the 64-bit
    virtual offset (``make_virtual_offset``); or give ``virtual_offset`` to split it into the
    two parts (``split_virtual_offset``). Exactly one of the two forms must be used.
    """

    making = block_start is not None or within_block is not None
    if making == (virtual_offset is not None):
        raise ConfigurationError(
            "Give block_start and within_block, or virtual_offset, but not both or neither."
        )
    bgzf = require("Bio.bgzf", _EXTRA)
    if making:
        start = _c.integer(block_start, "block_start", minimum=0, maximum=2**48 - 1)
        within = _c.integer(within_block, "within_block", minimum=0, maximum=2**16 - 1)
        combined = int(bgzf.make_virtual_offset(start, within))
        return {
            "status": "completed",
            "operation": "bgzf_virtual_offset",
            "direction": "make",
            "block_start": start,
            "within_block": within,
            "virtual_offset": combined,
            "provenance": _c.provenance(method="Bio.bgzf.make_virtual_offset"),
        }
    value = _c.integer(virtual_offset, "virtual_offset", minimum=0, maximum=2**64 - 1)
    start, within = bgzf.split_virtual_offset(value)
    return {
        "status": "completed",
        "operation": "bgzf_virtual_offset",
        "direction": "split",
        "block_start": int(start),
        "within_block": int(within),
        "virtual_offset": value,
        "provenance": _c.provenance(method="Bio.bgzf.split_virtual_offset"),
    }


__all__ = [
    "bgzf_virtual_offset",
    "convert_sequence_format",
    "count_subsequence",
    "filter_go_annotations",
    "format_primersearch_input",
    "format_substitution_matrix",
    "parse_entrez_xml",
    "parse_feature_location",
    "read_sequence_file",
    "read_sff_manifest",
    "search_subsequences",
    "split_sequence",
    "trim_sequence",
    "write_sequence_file",
]

_register_public_verbose(globals(), tuple(__all__))
