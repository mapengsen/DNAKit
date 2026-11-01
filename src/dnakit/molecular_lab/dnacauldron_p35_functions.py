"""DNA Cauldron wrappers, pack P35: records, sticky ends, assembly mixes, graphs, plans and BASIC assembly.

Every function imports DNA Cauldron (and Biopython, networkx, matplotlib) only when it runs. A missing
package raises ``BackendUnavailableError`` with the ``dnakit[bio-dnacauldron]`` install hint.

Conventions shared by the functions in this module:

- DNA is passed as plain strings of A, C, G and T. Input is upper-cased (DNA Cauldron itself treats sticky
  ends as upper case).
- A *record* is a dictionary ``{"id", "sequence", "topology", "features"}`` where each feature is
  ``{"type", "start", "end", "strand", "label"}`` (zero-based, end exclusive, strand 1 or -1).
- A *sticky-end fragment* is ``{"sequence", "left_end", "right_end"}`` where an end is ``None`` or
  ``{"sequence", "strand"}``. As in DNA Cauldron, the strand is the strand that carries the protruding bases.
- A *mix* is built from ``parts`` (name to DNA string). ``mix_type`` selects ``restriction_ligation``
  (Type-2S or other enzymes; the parts are circular by default), ``homology`` (Gibson-style end homologies;
  linear parts) or ``ligase_cycling`` (bridging oligos; linear parts). Fragment labels are the part name,
  with ``_r`` appended for the reverse-complement fragment.
- Plot functions write one PNG to a new path and return the drawn data.

Upstream quirks that are handled explicitly: the Ligase Cycling Reaction mix of DNA Cauldron matches
bridging oligos as given (they must be written as they appear in the concatenated fragments, although the
class names them "direct sense"); the BASIC assembly error path raises a bare ``TypeError`` in DNA Cauldron
2.0.12 and is reported as a ``ConfigurationError``.
"""

from __future__ import annotations

import contextlib
import importlib
import importlib.metadata
import itertools
import logging
import math
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnacauldron"
_MIX_TYPES = ("restriction_ligation", "homology", "ligase_cycling")
_MAX_DNA = 200_000
_MAX_PARTS = 200
_MAX_RESULTS = 1000
_MAX_NODES = 2000
_PLOT_DPI = 150
_NAME = re.compile(r"^[^\x00\n\r]{1,80}$")
_STRANDS = (1, -1)


# ---------------------------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------------------------


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not minimum <= float(value) <= maximum:
        raise ConfigurationError(f"{name} must be in [{minimum}, {maximum}].")
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of: {', '.join(allowed)}.")
    return value


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not _NAME.match(value.strip() or "\n"):
        raise ConfigurationError(f"{name} must be a non-empty name of at most 80 characters.")
    return value.strip()


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA, minimum: int = 1) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a DNA string.")
    text = "".join(value.split()).upper()
    if not minimum <= len(text) <= maximum:
        raise ConfigurationError(f"{name} must have between {minimum} and {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _named_dna(value: object, name: str, *, minimum_items: int = 1) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping of name to DNA string.")
    if not minimum_items <= len(value) <= _MAX_PARTS:
        raise ConfigurationError(f"{name} must contain between {minimum_items} and {_MAX_PARTS} items.")
    result: dict[str, str] = {}
    total = 0
    for key, text in value.items():
        label = _name(key, f"{name} key")
        result[label] = _dna(text, f"{name}[{label!r}]")
        total += len(result[label])
    if total > _MAX_DNA:
        raise ConfigurationError(f"{name} exceeds {_MAX_DNA} bases in total.")
    return result


def _new_file(value: object, name: str, suffixes: Sequence[str] | None = None) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(value).expanduser().resolve()
    if suffixes is not None and target.suffix.lower() not in suffixes:
        raise ConfigurationError(f"{name} must end with {' or '.join(suffixes)}.")
    if target.exists():
        raise ConfigurationError(f"Output file already exists and will not be overwritten: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def _existing_dir(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing folder: {path}")
    return path


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version("dnacauldron")
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        installed = None
    return {"package": "dnacauldron", "version": installed, **notes}


def _run(action: Any, message: str) -> Any:
    mix_error = require("dnacauldron.AssemblyMix", _EXTRA).AssemblyMixError
    try:
        return action()
    except ConfigurationError:
        raise
    except (
        AssertionError,
        AttributeError,
        IndexError,
        KeyError,
        RuntimeError,
        TypeError,
        ValueError,
        mix_error,
    ) as exc:
        node_error = type(exc).__name__ == "NodeNotFound"
        detail = str(exc) or type(exc).__name__
        raise ConfigurationError(f"{message}: {'no connection in the mix: ' if node_error else ''}{detail}") from exc
    except Exception as exc:  # networkx and Biopython raise their own classes
        if type(exc).__module__.split(".")[0] in {"networkx", "Bio"}:
            raise ConfigurationError(f"{message}: {exc}") from exc
        raise


# ---------------------------------------------------------------------------------------------
# Records and sticky ends
# ---------------------------------------------------------------------------------------------


def _enzyme_objects(names: object, name: str = "enzymes") -> list[Any]:
    restriction = require("Bio.Restriction", _EXTRA)
    if isinstance(names, str):
        names = [names]
    if not isinstance(names, Sequence) or not 1 <= len(names) <= 8:
        raise ConfigurationError(f"{name} must list 1 to 8 restriction enzyme names.")
    objects = []
    for item in names:
        try:
            enzyme = restriction.AllEnzymes.get(item, add=False) if isinstance(item, str) else None
        except ValueError:
            enzyme = None
        if enzyme is None:
            raise ConfigurationError(f"{name} holds an unknown restriction enzyme: {item!r}.")
        objects.append(enzyme)
    return objects


def _enzyme_names(names: object, name: str = "enzymes") -> list[str]:
    objects = _enzyme_objects(names, name)
    return [str(obj) for obj in objects]


def _make_record(record_id: str, sequence: str, topology: str) -> Any:
    biotools = require("dnacauldron.biotools", _EXTRA)
    record = biotools.sequence_to_biopython_record(sequence, id=record_id)
    biotools.set_record_topology(record, topology)
    return record


def _feature_rows(record: Any) -> list[dict[str, Any]]:
    rows = []
    for feature in record.features:
        if feature.location is None:
            continue
        label = feature.qualifiers.get("label", feature.qualifiers.get("note", ""))
        if isinstance(label, list):
            label = " ".join(str(item) for item in label)
        rows.append(
            {
                "type": str(feature.type),
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "strand": feature.location.strand,
                "label": str(label),
            }
        )
    return rows


def _record_dict(record: Any) -> dict[str, Any]:
    sequence = str(record.seq).upper()
    return {
        "id": str(record.id),
        "sequence": sequence,
        "length": len(sequence),
        "topology": str(record.annotations.get("topology", "linear")),
        "features": _feature_rows(record),
    }


def _apply_features(record: Any, features: object) -> int:
    biotools = require("dnacauldron.biotools", _EXTRA)
    if features is None:
        return 0
    if isinstance(features, (str, bytes)) or not isinstance(features, Sequence):
        raise ConfigurationError("features must be a list of feature dictionaries.")
    if len(features) > 1000:
        raise ConfigurationError("At most 1000 features can be added at once.")
    length = len(record.seq)
    for index, item in enumerate(features):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"features[{index}] must be a dictionary.")
        unknown = set(item) - {"type", "start", "end", "strand", "label", "qualifiers"}
        if unknown:
            raise ConfigurationError(f"features[{index}] has unknown keys: {sorted(unknown)}.")
        start = item.get("start", 0)
        end = item.get("end", length)
        _integer(start, f"features[{index}].start", minimum=0, maximum=length)
        _integer(end, f"features[{index}].end", minimum=1, maximum=length)
        if end <= start:
            raise ConfigurationError(f"features[{index}] must have start < end.")
        strand = item.get("strand", 1)
        if strand not in _STRANDS:
            raise ConfigurationError(f"features[{index}].strand must be 1 or -1.")
        feature_type = _name(item.get("type", "misc_feature"), f"features[{index}].type")
        qualifiers = dict(item.get("qualifiers") or {})
        if "label" in item:
            qualifiers["label"] = str(item["label"])
        biotools.annotate_record(
            record, location=(start, end, strand), feature_type=feature_type, **{str(k): v for k, v in qualifiers.items()}
        )
    return len(features)


def _record_from_spec(spec: object, name: str = "record") -> Any:
    if not isinstance(spec, Mapping):
        raise ConfigurationError(f"{name} must be a record dictionary with id and sequence.")
    unknown = set(spec) - {"id", "sequence", "topology", "features", "length"}
    if unknown:
        raise ConfigurationError(f"{name} has unknown keys: {sorted(unknown)}.")
    topology = _choice(spec.get("topology", "linear"), f"{name}.topology", ("linear", "circular"))
    record = _make_record(
        _name(spec.get("id"), f"{name}.id"), _dna(spec.get("sequence"), f"{name}.sequence"), topology
    )
    _apply_features(record, spec.get("features"))
    return record


def build_dnacauldron_record(
    sequence: str,
    *,
    record_id: str = "record",
    topology: str = "linear",
    features: list[dict[str, Any]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an annotated DNA record from a raw DNA string with DNA Cauldron.

    ``topology`` is ``linear`` or ``circular`` (DNA Cauldron's ``set_record_topology``). Each item of
    ``features`` is ``{"type", "start", "end", "strand", "label"}`` and is added with DNA Cauldron's
    ``annotate_record`` (type default ``misc_feature``, strand default 1, span default the whole record).
    """

    _choice(topology, "topology", ("linear", "circular"))
    record = _make_record(_name(record_id, "record_id"), _dna(sequence, "sequence"), topology)
    added = _apply_features(record, features)
    return {
        "status": "completed",
        "operation": "build_dnacauldron_record",
        "record": _record_dict(record),
        "features_added": added,
        "provenance": _provenance(method="dnacauldron.biotools.sequence_to_biopython_record"),
    }


def annotate_dnacauldron_record(
    record: Mapping[str, Any],
    *,
    features: list[dict[str, Any]] | None = None,
    topology: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add feature annotations to a DNA record and/or set its topology (circular or linear).

    ``record`` is a record dictionary (see the module description). Each item of ``features`` is
    ``{"type", "start", "end", "strand", "label"}`` and is added with DNA Cauldron's ``annotate_record``;
    ``topology`` (``linear`` or ``circular``) is applied with ``set_record_topology``. The record you pass is
    not changed; the annotated record is returned.
    """

    biotools = require("dnacauldron.biotools", _EXTRA)
    if topology is not None:
        _choice(topology, "topology", ("linear", "circular"))
    rec = _record_from_spec(record)
    added = _apply_features(rec, features)
    if topology is not None:
        biotools.set_record_topology(rec, topology)
    return {
        "status": "completed",
        "operation": "annotate_dnacauldron_record",
        "record": _record_dict(rec),
        "features_added": added,
        "topology_set": topology,
        "provenance": _provenance(method="dnacauldron.biotools.annotate_record / set_record_topology"),
    }


def write_dnacauldron_records(
    records: list[dict[str, Any]],
    output_dir: str,
    *,
    output_format: str = "genbank",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write DNA records to GenBank or FASTA files, one file per record, with DNA Cauldron's ``write_record``.

    ``output_format`` is ``genbank`` (files ``<id>.gb``, as ``BackboneChoice.write_final_records`` writes) or
    ``fasta`` (files ``<id>.fa``). Files are created in ``output_dir`` and are never overwritten. As in DNA
    Cauldron, GenBank locus names are cut to 20 characters; the file name keeps the full id.
    """

    biotools = require("dnacauldron.biotools", _EXTRA)
    fmt = _choice(output_format, "output_format", ("genbank", "fasta"))
    if isinstance(records, (str, bytes)) or not isinstance(records, Sequence) or not 1 <= len(records) <= 500:
        raise ConfigurationError("records must be a list of 1 to 500 record dictionaries.")
    prepared = [_record_from_spec(item, f"records[{i}]") for i, item in enumerate(records)]
    ids = [rec.id for rec in prepared]
    if len(set(ids)) != len(ids):
        raise ConfigurationError("record ids must be unique because they name the files.")
    if any("/" in item or "\\" in item or item in {".", ".."} for item in ids):
        raise ConfigurationError("record ids must not contain path separators.")
    folder = Path(output_dir).expanduser().resolve()
    suffix = ".gb" if fmt == "genbank" else ".fa"
    targets = [folder / f"{rec_id}{suffix}" for rec_id in ids]
    for target in targets:
        if target.exists():
            raise ConfigurationError(f"Output file already exists and will not be overwritten: {target}")
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for rec, target in zip(prepared, targets):
        with open(target, "w", encoding="utf-8") as handle:
            biotools.write_record(rec, handle, fmt)
        written.append({"id": rec.id, "path": str(target), "bytes": target.stat().st_size})
    return {
        "status": "completed",
        "operation": "write_dnacauldron_records",
        "format": fmt,
        "files": written,
        "provenance": _provenance(method="dnacauldron.biotools.sequence_io.write_record"),
    }


def get_dnacauldron_repository_records(
    record_ids: list[str],
    *,
    files: list[str] | None = None,
    folder: str | None = None,
    collection: str = "parts",
    use_file_names_as_ids: bool = True,
    topology: str = "default_to_linear",
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve sequence records by ID from a DNA Cauldron sequence repository built from local files.

    The repository is filled from ``files`` (FASTA, GenBank, SnapGene or zip) and/or a ``folder`` through
    ``SequenceRepository.import_records`` (``collection`` names the collection; ``topology`` is ``linear``,
    ``circular``, ``default_to_linear`` or ``default_to_circular``). ``record_ids`` lists the records to
    return (``get_record`` for one ID, ``get_records`` for several); an ID that is missing raises
    ``ConfigurationError`` with DNA Cauldron's name suggestions.
    """

    dc = require("dnacauldron", _EXTRA)
    _choice(topology, "topology", ("linear", "circular", "default_to_linear", "default_to_circular"))
    if isinstance(record_ids, (str, bytes)) or not isinstance(record_ids, Sequence) or not 1 <= len(record_ids) <= 500:
        raise ConfigurationError("record_ids must list 1 to 500 record IDs.")
    ids = [_name(item, "record id") for item in record_ids]
    paths = [str(_existing_file(item, "files item")) for item in (files or [])]
    folder_path = None if folder is None else str(_existing_dir(folder, "folder"))
    if not paths and folder_path is None:
        raise ConfigurationError("Provide files or a folder to fill the repository.")
    repository = dc.SequenceRepository()
    _run(
        lambda: repository.import_records(
            files=paths or None,
            folder=None if paths else folder_path,
            collection=_name(collection, "collection"),
            use_file_names_as_ids=_flag(use_file_names_as_ids, "use_file_names_as_ids"),
            topology=topology,
        ),
        "DNA Cauldron could not import the records",
    )
    with _clean_ids(dc):
        found = _run(
            lambda: [repository.get_record(ids[0])] if len(ids) == 1 else repository.get_records(ids),
            "DNA Cauldron could not find the records",
        )
    return {
        "status": "completed",
        "operation": "get_dnacauldron_repository_records",
        "repository_ids": repository.get_all_part_names(),
        "records": [_record_dict(item) for item in found],
        "provenance": _provenance(method="dnacauldron.SequenceRepository.get_record / get_records"),
    }


@contextlib.contextmanager
def _clean_ids(dc: Any) -> Any:
    """Turn the repository's NotInRepositoryError into a plain message."""

    try:
        yield
    except dc.NotInRepositoryError as exc:
        raise ConfigurationError(str(exc)) from exc


def _end_to_object(spec: object, name: str) -> Any:
    fragment_module = require("dnacauldron.Fragment", _EXTRA)
    if spec is None:
        return None
    if not isinstance(spec, Mapping) or set(spec) - {"sequence", "strand"}:
        raise ConfigurationError(f"{name} must be null or {{\"sequence\", \"strand\"}}.")
    strand = spec.get("strand")
    if strand not in _STRANDS:
        raise ConfigurationError(f"{name}.strand must be 1 or -1.")
    return fragment_module.StickyEnd(_dna(spec.get("sequence"), f"{name}.sequence", maximum=2000), strand=strand)


def _end_dict(end: Any) -> dict[str, Any] | None:
    if end is None:
        return None
    return {"sequence": str(end), "strand": int(end.strand)}


def _sticky_from_spec(spec: object, name: str = "fragment") -> Any:
    fragment_module = require("dnacauldron.Fragment", _EXTRA)
    if not isinstance(spec, Mapping) or set(spec) - {"sequence", "left_end", "right_end", "length"}:
        raise ConfigurationError(f"{name} must be {{\"sequence\", \"left_end\", \"right_end\"}}.")
    return fragment_module.StickyEndSeq(
        _dna(spec.get("sequence"), f"{name}.sequence", minimum=0 if not spec.get("sequence") else 1),
        left_end=_end_to_object(spec.get("left_end"), f"{name}.left_end"),
        right_end=_end_to_object(spec.get("right_end"), f"{name}.right_end"),
    )


def _sticky_dict(seq: Any) -> dict[str, Any]:
    plain = str(seq)
    return {
        "sequence": plain,
        "length": len(plain),
        "left_end": _end_dict(getattr(seq, "left_end", None)),
        "right_end": _end_dict(getattr(seq, "right_end", None)),
    }


def digest_dnacauldron_sticky_ends(
    sequence: str,
    enzymes: list[str],
    *,
    circular: bool = False,
    method: str = "sequence",
    verbose: bool = True,
) -> dict[str, Any]:
    """Digest DNA with one or several restriction enzymes into sticky-end fragments (DNA Cauldron).

    ``method`` is ``sequence`` (``StickyEndSeq.list_from_sequence_digestion``) or ``record``
    (``StickyEndFragment.list_from_record_digestion``, the route used by assembly mixes). Each fragment
    lists its middle ``sequence`` and the left and right protruding ends with their strand. A sequence with
    no site returns one fragment without ends. ``enzymes`` are Biopython restriction enzyme names.
    """

    dc = require("dnacauldron", _EXTRA)
    seq_class = require("Bio.Seq", _EXTRA).Seq
    chosen = _enzyme_objects(enzymes)
    text = _dna(sequence, "sequence")
    kind = _choice(method, "method", ("sequence", "record"))
    is_circular = _flag(circular, "circular")
    if kind == "sequence":
        fragments = _run(
            lambda: dc.StickyEndSeq.list_from_sequence_digestion(
                seq_class(text), chosen if len(chosen) > 1 else chosen[0], linear=not is_circular
            ),
            "DNA Cauldron could not digest the sequence",
        )
        rows = [_sticky_dict(item) for item in fragments]
    else:
        record = _make_record("digested", text, "circular" if is_circular else "linear")
        fragments = _run(
            lambda: dc.StickyEndFragment.list_from_record_digestion(record, chosen if len(chosen) > 1 else chosen[0]),
            "DNA Cauldron could not digest the record",
        )
        rows = [_sticky_dict(item.seq) for item in fragments]
    return {
        "status": "completed",
        "operation": "digest_dnacauldron_sticky_ends",
        "enzymes": [str(item) for item in chosen],
        "circular": is_circular,
        "fragment_count": len(rows),
        "fragments": rows,
        "provenance": _provenance(method=f"dnacauldron StickyEnd {kind} digestion"),
    }


def transform_dnacauldron_sticky_sequence(
    fragment: dict[str, Any],
    *,
    operation: str = "reverse_complement",
    start: int | None = None,
    end: int | None = None,
    discard_sticky_ends: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement, flatten or slice a sticky-end sequence.

    ``operation`` is ``reverse_complement`` (the ends are swapped and reverse-complemented), ``to_standard``
    (the plain DNA string left end + middle + right end; with ``discard_sticky_ends`` only the middle) or
    ``slice`` (``start`` and ``end`` are zero-based coordinates on the middle sequence; the overhangs are
    dropped, as ``StickyEndSeq.slice_seq`` does).
    """

    kind = _choice(operation, "operation", ("reverse_complement", "to_standard", "slice"))
    seq = _sticky_from_spec(fragment)
    if kind == "reverse_complement":
        result = _run(seq.reverse_complement, "DNA Cauldron could not reverse-complement the fragment")
        return {
            "status": "completed",
            "operation": "transform_dnacauldron_sticky_sequence",
            "transform": kind,
            "fragment": _sticky_dict(result),
            "provenance": _provenance(method="StickyEndSeq.reverse_complement"),
        }
    if kind == "to_standard":
        text = str(seq.to_standard_sequence(discard_sticky_ends=_flag(discard_sticky_ends, "discard_sticky_ends")))
        return {
            "status": "completed",
            "operation": "transform_dnacauldron_sticky_sequence",
            "transform": kind,
            "sequence": text,
            "length": len(text),
            "provenance": _provenance(method="StickyEndSeq.to_standard_sequence"),
        }
    for label, value in (("start", start), ("end", end)):
        if value is not None:
            _integer(value, label, minimum=-_MAX_DNA, maximum=_MAX_DNA)
    sliced = seq.slice_seq(start=start, end=end)
    return {
        "status": "completed",
        "operation": "transform_dnacauldron_sticky_sequence",
        "transform": kind,
        "fragment": _sticky_dict(sliced),
        "provenance": _provenance(method="StickyEndSeq.slice_seq"),
    }


def check_dnacauldron_sticky_end_compatibility(
    fragments: list[dict[str, Any]],
    *,
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check whether sticky-end fragments clip together in the given order.

    Fragment ``i`` clips with fragment ``i+1`` when its right end has the same bases as the next left end on
    the opposite strand (``StickyEnd.will_clip_directly_with``; empty ends never clip). With
    ``circular=True`` the last fragment is also checked against the first.
    """

    if isinstance(fragments, (str, bytes)) or not isinstance(fragments, Sequence) or not 2 <= len(fragments) <= 200:
        raise ConfigurationError("fragments must list 2 to 200 sticky-end fragments.")
    objects = [_sticky_from_spec(item, f"fragments[{i}]") for i, item in enumerate(fragments)]
    pairs = list(zip(range(len(objects)), range(1, len(objects))))
    if _flag(circular, "circular"):
        pairs.append((len(objects) - 1, 0))
    rows = []
    for first, second in pairs:
        right, left = objects[first].right_end, objects[second].left_end
        rows.append(
            {
                "from": first,
                "to": second,
                "right_end": None if right is None else str(right),
                "left_end": None if left is None else str(left),
                "compatible": bool(right is not None and right.will_clip_directly_with(left)),
            }
        )
    return {
        "status": "completed",
        "operation": "check_dnacauldron_sticky_end_compatibility",
        "pairs": rows,
        "all_compatible": all(row["compatible"] for row in rows),
        "provenance": _provenance(method="StickyEnd.will_clip_directly_with"),
    }


def sticky_fragment_to_dnacauldron_record(
    fragment: dict[str, Any],
    *,
    part: str = "fragment",
    record_id: str = "fragment",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a sticky-end fragment or one of its overhangs to an annotated record.

    ``part`` is ``fragment`` (left overhang + middle + right overhang, the overhangs annotated with their
    strand, as ``StickyEndFragment.as_biopython_record``), ``left_end`` or ``right_end``
    (``StickyEnd.as_biopython_record``, one record labelled ``(+) strand`` or ``(-) strand``).
    """

    dc = require("dnacauldron", _EXTRA)
    kind = _choice(part, "part", ("fragment", "left_end", "right_end"))
    seq = _sticky_from_spec(fragment)
    if kind == "fragment":
        if seq.left_end is None or seq.right_end is None:
            raise ConfigurationError("fragment needs both a left_end and a right_end.")
        sticky = dc.StickyEndFragment(seq=seq, id=_name(record_id, "record_id"))
        record = _run(sticky.as_biopython_record, "DNA Cauldron could not convert the fragment")
    else:
        end = seq.left_end if kind == "left_end" else seq.right_end
        if end is None:
            raise ConfigurationError(f"fragment has no {kind}.")
        record = end.as_biopython_record()
        record.id = _name(record_id, "record_id")
    return {
        "status": "completed",
        "operation": "sticky_fragment_to_dnacauldron_record",
        "part": kind,
        "record": _record_dict(record),
        "provenance": _provenance(method="as_biopython_record"),
    }


def check_dnacauldron_homology(
    sequence: str,
    other_sequence: str | None = None,
    *,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a homology between two sequences by length, melting temperature and edit distance.

    This is ``HomologyChecker.check_homology`` with DNA Cauldron's defaults (15 to 80 bases, no Tm bounds,
    edit distance 0). ``other_sequence`` is optional: without it only the length and Tm are checked.
    The melting temperature is DNA Cauldron's rule 2 degrees per A/T and 4 per G/C. With
    ``max_distance`` above 0 small edits are accepted (Levenshtein distance).
    """

    dc = require("dnacauldron", _EXTRA)
    first = _dna(sequence, "sequence", maximum=10_000)
    second = None if other_sequence is None else _dna(other_sequence, "other_sequence", maximum=10_000)
    checker = _homology_checker(dc, min_size, max_size, min_tm, max_tm, max_distance)
    distance = None
    if second is not None:
        distance = int(require("Levenshtein", _EXTRA).distance(first, second))
    return {
        "status": "completed",
        "operation": "check_dnacauldron_homology",
        "acceptable": bool(checker.check_homology(first, second)),
        "length": len(first),
        "tm": checker.compute_tm(first),
        "edit_distance": distance,
        "parameters": checker.parameters_as_string(),
        "provenance": _provenance(method="HomologyChecker.check_homology"),
    }


def _homology_checker(
    dc: Any, min_size: int, max_size: int | None, min_tm: float, max_tm: float | None, max_distance: int
) -> Any:
    _integer(min_size, "min_size", minimum=1, maximum=10_000)
    if max_size is not None:
        _integer(max_size, "max_size", minimum=min_size, maximum=100_000)
    _number(min_tm, "min_tm", minimum=0, maximum=1000)
    if max_tm is not None:
        _number(max_tm, "max_tm", minimum=min_tm, maximum=100_000)
    _integer(max_distance, "max_distance", minimum=0, maximum=100)
    return dc.HomologyChecker(
        min_size=min_size, max_size=max_size, min_tm=min_tm, max_tm=max_tm, max_distance=max_distance
    )


# ---------------------------------------------------------------------------------------------
# Assembly mixes
# ---------------------------------------------------------------------------------------------


def _build_mix(
    mix_type: str,
    parts: Mapping[str, str],
    *,
    enzymes: Sequence[str],
    topology: str,
    filter_site_fragments: bool,
    bridging_oligos: Mapping[str, str] | None,
    min_size: int,
    max_size: int | None,
    min_tm: float,
    max_tm: float | None,
    max_distance: int,
    mix_name: str = "mix",
) -> Any:
    dc = require("dnacauldron", _EXTRA)
    kind = _choice(mix_type, "mix_type", _MIX_TYPES)
    named = _named_dna(parts, "parts")
    topology = _choice(topology, "topology", ("auto", "linear", "circular"))
    if topology == "auto":
        topology = "circular" if kind == "restriction_ligation" else "linear"
    records = [_make_record(key, text, topology) for key, text in sorted(named.items())]
    if kind == "restriction_ligation":
        names = _enzyme_names(enzymes)
        filters = [dc.NoRestrictionSiteFilter(item) for item in names] if _flag(filter_site_fragments, "filter_site_fragments") else []
        return _run(
            lambda: dc.RestrictionLigationMix(
                parts=records, enzymes=names, fragment_filters=filters, name=mix_name
            ),
            "DNA Cauldron could not build the restriction-ligation mix",
        )
    checker = _homology_checker(dc, min_size, max_size, min_tm, max_tm, max_distance)
    if kind == "homology":
        return _run(
            lambda: dc.HomologousAssemblyMix(records, checker, name=mix_name),
            "DNA Cauldron could not build the homology mix",
        )
    oligos = _named_dna(bridging_oligos, "bridging_oligos")
    oligo_records = [_make_record(key, text, "linear") for key, text in sorted(oligos.items())]
    lcr = importlib.import_module("dnacauldron.AssemblyMix.LigaseCyclingReactionMix")
    return _run(
        lambda: lcr.LigaseCyclingReactionMix(records, oligo_records, homology_checker=checker, name=mix_name),
        "DNA Cauldron could not build the ligase cycling mix",
    )


def _label(fragment: Any) -> str:
    return f"{fragment.original_part.id}{'_r' if getattr(fragment, 'is_reversed', False) else ''}"


def _fragment_row(index: int, fragment: Any) -> dict[str, Any]:
    left = getattr(fragment.seq, "left_end", None)
    right = getattr(fragment.seq, "right_end", None)
    return {
        "index": int(index),
        "label": _label(fragment),
        "part": str(fragment.original_part.id),
        "reversed": bool(getattr(fragment, "is_reversed", False)),
        "length": len(fragment),
        "left_end": None if left is None else str(left),
        "right_end": None if right is None else str(right),
    }


def _construct_row(record: Any) -> dict[str, Any]:
    seq = record.seq
    sequence = str(seq.to_standard_sequence()).upper() if hasattr(seq, "to_standard_sequence") else str(seq).upper()
    fragments = getattr(record, "fragments", None)
    return {
        "id": str(record.id),
        "sequence": sequence,
        "length": len(sequence),
        "topology": str(record.annotations.get("topology", "linear")),
        "parts": [] if fragments is None else [_label(item) if hasattr(item, "original_part") else str(item.id) for item in fragments],
    }


def _cap(value: object, name: str = "max_results") -> int:
    return _integer(value, name, minimum=1, maximum=_MAX_RESULTS)


def compute_dnacauldron_slots(
    parts: dict[str, str],
    *,
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the slot of each fragment of a restriction-ligation mix (standardized overhang pairs).

    A fragment's slot is the alphabetically smaller of its (left, right) overhangs and the reverse
    complement of the reversed pair. Parts that share a slot are interchangeable in an assembly. ``parts``
    maps names to DNA strings (circular by default); ``enzymes`` default to ``["BsmBI"]``.
    """

    mix = _build_mix(
        "restriction_ligation",
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=None,
        min_size=15,
        max_size=80,
        min_tm=0,
        max_tm=None,
        max_distance=0,
    )
    slots = _run(mix.compute_slots, "DNA Cauldron could not compute the slots")
    rows = [
        {"left_overhang": left, "right_overhang": right, "parts": sorted(names)}
        for (left, right), names in sorted(slots.items())
    ]
    return {
        "status": "completed",
        "operation": "compute_dnacauldron_slots",
        "slot_count": len(rows),
        "slots": rows,
        "provenance": _provenance(method="StickyEndAssemblyMix.compute_slots"),
    }


def _node_name(node: Any) -> str:
    return "-".join(node) if isinstance(node, tuple) else str(node)


def compute_dnacauldron_slots_graph(
    parts: dict[str, str],
    *,
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    with_overhangs: bool = True,
    directed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the slots graph of a restriction-ligation mix as nodes and edges.

    A slot node is written ``LEFT-RIGHT`` (its two overhangs). With ``with_overhangs`` the overhangs are
    nodes too, and with ``directed`` the edges carry DNA Cauldron's arrow directions; without overhangs two
    slots are linked when they share an overhang.
    """

    mix = _build_mix(
        "restriction_ligation",
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=None,
        min_size=15,
        max_size=80,
        min_tm=0,
        max_tm=None,
        max_distance=0,
    )
    graph = _run(
        lambda: mix.slots_graph(with_overhangs=_flag(with_overhangs, "with_overhangs"), directed=_flag(directed, "directed")),
        "DNA Cauldron could not compute the slots graph",
    )
    slot_nodes = {_node_name(slot) for slot in mix.compute_slots()}
    nodes = [{"id": _node_name(n), "kind": "slot" if _node_name(n) in slot_nodes else "overhang"} for n in graph.nodes()]
    edges = [[_node_name(a), _node_name(b)] for a, b in graph.edges()]
    return {
        "status": "completed",
        "operation": "compute_dnacauldron_slots_graph",
        "directed": bool(graph.is_directed()),
        "node_count": len(nodes),
        "edge_count": len(edges),
        "nodes": sorted(nodes, key=lambda item: item["id"]),
        "edges": sorted(edges),
        "provenance": _provenance(method="StickyEndAssemblyMix.slots_graph"),
    }


_MIX_ARGUMENTS = (
    "mix_type",
    "enzymes",
    "topology",
    "filter_site_fragments",
    "bridging_oligos",
    "min_size",
    "max_size",
    "min_tm",
    "max_tm",
    "max_distance",
)


def compute_dnacauldron_connections_graph(
    parts: dict[str, str],
    *,
    mix_type: str = "restriction_ligation",
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    bridging_oligos: dict[str, str] | None = None,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute which fragments of an assembly mix clip together and in which order.

    ``mix_type`` is ``restriction_ligation`` (digestion by ``enzymes``, default ``["BsmBI"]``),
    ``homology`` (end homologies checked with ``min_size``/``max_size``/``min_tm``/``max_tm``/
    ``max_distance``, Gibson style) or ``ligase_cycling`` (needs ``bridging_oligos``). The result lists the
    fragments (index, label, overhangs) and the directed edges ``[i, j]``: fragment ``i`` clips before ``j``.
    Fragments rejected by the site filter are left out of the edges.
    """

    mix = _build_mix(
        mix_type,
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )
    graph = mix.filtered_connections_graph
    if len(mix.fragments_dict) > _MAX_NODES:
        raise ConfigurationError(f"The mix has more than {_MAX_NODES} fragments.")
    fragments = [_fragment_row(index, mix.fragments_dict[index]) for index in sorted(mix.fragments_dict)]
    edges = sorted([int(a), int(b)] for a, b in graph.edges())
    return {
        "status": "completed",
        "operation": "compute_dnacauldron_connections_graph",
        "mix_type": mix_type,
        "fragment_count": len(fragments),
        "fragments": fragments,
        "edge_count": len(edges),
        "edges": edges,
        "provenance": _provenance(method="AssemblyMix.compute_connections_graph / filtered_connections_graph"),
    }


def _with_pyplot() -> Any:
    matplotlib = require("matplotlib", "bio-plots")
    try:
        matplotlib.use("Agg", force=False)
    except Exception:  # pragma: no cover - backend already chosen
        pass
    return require("matplotlib.pyplot", "bio-plots")


@contextlib.contextmanager
def _quiet_fonts() -> Any:
    logger = logging.getLogger("matplotlib.font_manager")
    previous = logger.level
    logger.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            yield
    finally:
        logger.setLevel(previous)


def _save_png(figure: Any, target: Path, pyplot: Any) -> None:
    try:
        figure.savefig(str(target), format="png", dpi=_PLOT_DPI, bbox_inches="tight")
    finally:
        pyplot.close(figure)


def plot_dnacauldron_connections_graph(
    parts: dict[str, str],
    output_path: str,
    *,
    mix_type: str = "restriction_ligation",
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    bridging_oligos: dict[str, str] | None = None,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    figure_size: float = 12.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the fragment connection graph of an assembly mix to a new PNG file.

    The mix is built as in ``compute_dnacauldron_connections_graph``. DNA Cauldron draws the
    de-duplicated graph (``plot_connections_graph``): one node per fragment, labelled with the overhangs and
    the part name, one arrow per allowed junction. ``figure_size`` is the figure edge in inches.
    The result gives the PNG path and the drawn nodes and edges.
    """

    target = _new_file(output_path, "output_path", (".png",))
    size = _number(figure_size, "figure_size", minimum=2, maximum=40)
    mix = _build_mix(
        mix_type,
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
        mix_name=Path(target).stem,
    )
    pyplot = _with_pyplot()
    graph = mix.uniquified_connection_graph
    with _quiet_fonts():
        figure, axis = pyplot.subplots(1, 1, figsize=(size, size))
        try:
            _run(lambda: mix.plot_connections_graph(ax=axis), "DNA Cauldron could not plot the graph")
        except Exception:
            pyplot.close(figure)
            raise
        _save_png(figure, target, pyplot)
    return {
        "status": "completed",
        "operation": "plot_dnacauldron_connections_graph",
        "path": str(target),
        "data": {
            "node_count": graph.number_of_nodes(),
            "edge_count": graph.number_of_edges(),
            "nodes": [_fragment_row(int(n), mix.fragments_dict[n]) for n in sorted(graph.nodes())],
            "edges": sorted([int(a), int(b)] for a, b in graph.edges()),
        },
        "provenance": _provenance(method="AssemblyMix.plot_connections_graph"),
    }


def plot_dnacauldron_slots_graph(
    parts: dict[str, str],
    output_path: str,
    *,
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    with_overhangs: bool = False,
    directed: bool = True,
    highlighted_parts: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the slots graph of a restriction-ligation mix to a new PNG file.

    Slots are numbered nodes with a legend of the parts that fill each slot; ``highlighted_parts`` marks
    parts in blue; ``with_overhangs`` adds the overhangs as nodes (``StickyEndAssemblyMix.plot_slots_graph``).
    The result gives the PNG path and the slots graph data.
    """

    target = _new_file(output_path, "output_path", (".png",))
    mix = _build_mix(
        "restriction_ligation",
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=None,
        min_size=15,
        max_size=80,
        min_tm=0,
        max_tm=None,
        max_distance=0,
        mix_name=Path(target).stem,
    )
    highlighted = [_name(item, "highlighted_parts item") for item in (highlighted_parts or [])]
    pyplot = _with_pyplot()
    with _quiet_fonts():
        axis = _run(
            lambda: mix.plot_slots_graph(
                with_overhangs=_flag(with_overhangs, "with_overhangs"),
                directed=_flag(directed, "directed"),
                highlighted_parts=highlighted,
            ),
            "DNA Cauldron could not plot the slots graph",
        )
        _save_png(axis.figure, target, pyplot)
    slots = mix.compute_slots()
    return {
        "status": "completed",
        "operation": "plot_dnacauldron_slots_graph",
        "path": str(target),
        "data": {
            "slots": [
                {"left_overhang": left, "right_overhang": right, "parts": sorted(names)}
                for (left, right), names in sorted(slots.items())
            ],
            "highlighted_parts": highlighted,
        },
        "provenance": _provenance(method="StickyEndAssemblyMix.plot_slots_graph"),
    }


def _seeded(seed: int | None) -> Any:
    @contextlib.contextmanager
    def manager() -> Any:
        numpy = require("numpy", _EXTRA)
        state = numpy.random.get_state()
        try:
            if seed is not None:
                numpy.random.seed(_integer(seed, "random_seed", minimum=0, maximum=2**32 - 1))
            yield
        finally:
            numpy.random.set_state(state)

    return manager()


def enumerate_dnacauldron_circular_assemblies(
    parts: dict[str, str],
    *,
    mix_type: str = "restriction_ligation",
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    bridging_oligos: dict[str, str] | None = None,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    output: str = "records",
    sampling: str = "exhaustive",
    max_results: int = 20,
    annotate_homologies: bool = False,
    random_seed: int | None = None,
    staling_cutoff: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate or sample the circular constructs an assembly mix can form.

    ``output`` is ``fragments_sets`` (the ordered fragment labels of each construct,
    ``compute_circular_fragments_sets``) or ``records`` (the assembled circular sequences,
    ``compute_circular_assemblies``). ``sampling`` is ``exhaustive`` (networkx cycle order, at most
    ``max_results``) or ``random`` (``compute_random_circular_fragments_sets``; ``random_seed`` makes it
    repeatable; ``staling_cutoff`` stops when only repeats are found, reported as ``exhausted``). Each set
    is standardized, so rotations and reverse complements of one construct are reported once.
    """

    kind_output = _choice(output, "output", ("fragments_sets", "records"))
    kind_sampling = _choice(sampling, "sampling", ("exhaustive", "random"))
    limit = _cap(max_results)
    _integer(staling_cutoff, "staling_cutoff", minimum=1, maximum=100_000)
    mix = _build_mix(
        mix_type,
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )
    random_mode = kind_sampling == "random"
    exhausted = False
    results: list[Any] = []

    def collect() -> None:
        nonlocal exhausted
        if kind_output == "fragments_sets":
            iterator = (
                mix.compute_random_circular_fragments_sets(staling_cutoff=staling_cutoff)
                if random_mode
                else mix.compute_circular_fragments_sets()
            )
        else:
            iterator = mix.compute_circular_assemblies(
                annotate_parts_homologies=_flag(annotate_homologies, "annotate_homologies"),
                randomize=random_mode,
                randomization_staling_cutoff=staling_cutoff,
            )
        try:
            for item in itertools.islice(iterator, limit):
                results.append(item)
        except ValueError as exc:
            if "staled" not in str(exc):
                raise
            exhausted = True

    with _seeded(random_seed):
        _run(collect, "DNA Cauldron could not enumerate the circular assemblies")
    if kind_output == "fragments_sets":
        rows: list[Any] = [[_label(f) for f in item] for item in results]
    else:
        rows = [_construct_row(item) for item in results]
    return {
        "status": "completed",
        "operation": "enumerate_dnacauldron_circular_assemblies",
        "output": kind_output,
        "sampling": kind_sampling,
        "returned": len(rows),
        "limit_reached": len(rows) >= limit,
        "exhausted": exhausted,
        "assemblies": rows,
        "provenance": _provenance(method="ConstructsMixin circular enumeration"),
    }


def _linear_assemblies(mix: Any, min_parts: int, limit: int, annotate: bool) -> list[Any]:
    """The algorithm of ``ConstructsMixin.compute_linear_assemblies`` (shortest paths between fragments).

    The upstream method iterates ``nx.shortest_path(graph).items()``, which fails with networkx 3.x
    because a generator is returned; here the all-pairs paths are taken through
    ``dict(nx.shortest_path(graph))``.
    """

    nx = require("networkx", _EXTRA)
    chain_class = require("dnacauldron.Fragment.FragmentChain", _EXTRA).FragmentChain
    seen: set[int] = set()
    found_assemblies: list[Any] = []
    for _source, targets in dict(nx.shortest_path(mix.filtered_connections_graph)).items():
        for _target, path in targets.items():
            if len(path) < min_parts:
                continue
            fragments = [mix.fragments_dict[i] for i in path]
            chain_hash = hash(chain_class(fragments).standardized())
            if chain_hash in seen:
                continue
            seen.add(chain_hash)
            assembled = mix.assemble(fragments, annotate_homologies=annotate)
            if not hasattr(assembled, "fragments"):
                assembled.fragments = fragments
            found_assemblies.append(assembled)
            if len(found_assemblies) >= limit:
                return found_assemblies
    return found_assemblies


def enumerate_dnacauldron_linear_assemblies(
    parts: dict[str, str],
    *,
    mix_type: str = "restriction_ligation",
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    bridging_oligos: dict[str, str] | None = None,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    min_parts: int = 2,
    max_results: int = 20,
    annotate_homologies: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate the linear constructs of an assembly mix (shortest path between fragments).

    This is ``compute_linear_assemblies``: for every pair of connected fragments the shortest path through
    the connection graph is assembled, assemblies with fewer than ``min_parts`` fragments are dropped and
    equivalent chains are reported once. For sticky-end constructs ``left_end`` and ``right_end`` give the
    free overhangs and ``sequence`` includes them. DNA Cauldron calls this feature "a bit undertested".
    """

    limit = _cap(max_results)
    _integer(min_parts, "min_parts", minimum=1, maximum=1000)
    mix = _build_mix(
        mix_type,
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )

    annotate = _flag(annotate_homologies, "annotate_homologies")

    def collect() -> list[Any]:
        return _linear_assemblies(mix, min_parts, limit, annotate)

    found = _run(collect, "DNA Cauldron could not enumerate the linear assemblies")
    rows = []
    for item in found:
        row = _construct_row(item)
        row["left_end"] = None if getattr(item.seq, "left_end", None) is None else str(item.seq.left_end)
        row["right_end"] = None if getattr(item.seq, "right_end", None) is None else str(item.seq.right_end)
        rows.append(row)
    return {
        "status": "completed",
        "operation": "enumerate_dnacauldron_linear_assemblies",
        "returned": len(rows),
        "limit_reached": len(rows) >= limit,
        "assemblies": rows,
        "provenance": _provenance(method="ConstructsMixin.compute_linear_assemblies"),
    }


def transform_dnacauldron_fragment_chain(
    parts: dict[str, str],
    *,
    operation: str = "standardized",
    mix_type: str = "restriction_ligation",
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    bridging_oligos: dict[str, str] | None = None,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    max_results: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement or standardize the circular fragment chains of an assembly mix.

    The raw cycles of the connection graph are wrapped in DNA Cauldron ``FragmentChain`` objects. With
    ``operation="reverse_complement"`` each chain is replaced by the inverted list of reversed fragments;
    with ``"standardized"`` it is given its canonical rotation and orientation (backbone first when one part
    is marked as backbone, otherwise at most half reversed nucleotides, rotated to start at the largest
    fragment). The result lists, per cycle, the fragment labels before and after.
    """

    kind = _choice(operation, "operation", ("reverse_complement", "standardized"))
    limit = _cap(max_results)
    dc = require("dnacauldron", _EXTRA)
    nx = require("networkx", _EXTRA)
    chain_class = require("dnacauldron.Fragment.FragmentChain", _EXTRA).FragmentChain
    mix = _build_mix(
        mix_type,
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )
    del dc

    def collect() -> list[dict[str, Any]]:
        rows = []
        for cycle in itertools.islice(nx.simple_cycles(mix.filtered_connections_graph), limit):
            chain = chain_class([mix.fragments_dict[i] for i in cycle], is_cycle=True)
            result = chain.reverse_complement() if kind == "reverse_complement" else chain.standardized()
            rows.append(
                {
                    "before": [_label(f) for f in chain.fragments],
                    "after": [_label(f) for f in result.fragments],
                    "is_standardized": bool(result.is_standardized),
                }
            )
        return rows

    rows = _run(collect, "DNA Cauldron could not transform the fragment chains")
    return {
        "status": "completed",
        "operation": "transform_dnacauldron_fragment_chain",
        "transform": kind,
        "chain_count": len(rows),
        "chains": rows,
        "provenance": _provenance(method=f"FragmentChain.{kind}"),
    }


def select_dnacauldron_connectors(
    parts: dict[str, str],
    connectors: dict[str, str],
    *,
    enzymes: list[str] | None = None,
    topology: str = "auto",
    filter_site_fragments: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select the bridging connector parts needed to close a circular assembly (``autoselect_connectors``).

    ``parts`` are the main parts of the assembly (one per slot) and ``connectors`` is a pool of candidate
    bridging parts; both are circular DNA strings with the enzyme sites by default. When a cycle that
    contains every main part exists, the connectors it uses are returned (with the circular constructs of
    the completed mix); otherwise ``ConfigurationError`` says that no construct was found.
    """

    names = _named_dna(connectors, "connectors")
    mix = _build_mix(
        "restriction_ligation",
        parts,
        enzymes=enzymes or ["BsmBI"],
        topology=topology,
        filter_site_fragments=filter_site_fragments,
        bridging_oligos=None,
        min_size=15,
        max_size=80,
        min_tm=0,
        max_tm=None,
        max_distance=0,
    )
    topology_used = "circular" if topology == "auto" else topology
    records = [_make_record(key, text, topology_used) for key, text in sorted(names.items())]
    selected = _run(lambda: mix.autoselect_connectors(records), "DNA Cauldron could not select connectors")
    assemblies = _run(lambda: list(itertools.islice(mix.compute_circular_assemblies(), 20)), "DNA Cauldron could not assemble")
    return {
        "status": "completed",
        "operation": "select_dnacauldron_connectors",
        "selected_connectors": sorted({item.id for item in selected}),
        "connector_count": len(selected),
        "constructs": [_construct_row(item) for item in assemblies],
        "provenance": _provenance(method="ConnectorsMixin.autoselect_connectors"),
    }


def simulate_dnacauldron_lcr_assembly(
    parts: dict[str, str],
    bridging_oligos: dict[str, str],
    *,
    circular: bool = True,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    min_parts: int = 2,
    max_results: int = 20,
    annotate_homologies: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate a ligase cycling reaction assembly of linear parts joined by bridging oligos.

    ``parts`` are linear DNA strings; each bridging oligo must span a junction (its two halves are
    checked with DNA Cauldron's ``HomologyChecker`` defaults: 15 to 80 bases) and is written as it appears
    in the concatenation of the two neighbouring parts. With ``circular=True`` the circular constructs are
    returned, otherwise the linear ones (at least ``min_parts`` parts). Each construct lists its parts and
    the bridging oligos it used.
    """

    limit = _cap(max_results)
    mix = _build_mix(
        "ligase_cycling",
        parts,
        enzymes=[],
        topology="linear",
        filter_site_fragments=False,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )
    flag = _flag(annotate_homologies, "annotate_homologies")
    if _flag(circular, "circular"):
        generator = mix.compute_circular_assemblies(annotate_parts_homologies=flag)
    else:
        _integer(min_parts, "min_parts", minimum=1, maximum=1000)
        generator = iter(_linear_assemblies(mix, min_parts, limit, flag))
    found = _run(lambda: list(itertools.islice(generator, limit)), "DNA Cauldron could not simulate the assembly")
    oligo_ids = {oligo.id for oligo in mix.bridging_oligos}
    rows = []
    for item in found:
        row = _construct_row(item)
        fragments = [f for f in getattr(item, "fragments", []) if str(f.id) not in oligo_ids or hasattr(f, "original_part")]
        part_fragments = [f for f in fragments if hasattr(f, "original_part")]
        pairs = list(zip(part_fragments, part_fragments[1:]))
        if circular and len(part_fragments) > 1:
            pairs.append((part_fragments[-1], part_fragments[0]))
        used = []
        for first, second in pairs:
            hit = mix.find_oligo_homology(first, second)
            if hit is not None:
                used.append(str(hit[0].id))
        row["parts"] = [_label(f) for f in part_fragments]
        row["bridging_oligos"] = used
        rows.append(row)
    return {
        "status": "completed",
        "operation": "simulate_dnacauldron_lcr_assembly",
        "circular": bool(circular),
        "returned": len(rows),
        "constructs": rows,
        "provenance": _provenance(method="LigaseCyclingReactionMix.assemble"),
    }


def find_dnacauldron_bridging_oligo(
    left_sequence: str,
    right_sequence: str,
    bridging_oligos: dict[str, str],
    *,
    min_size: int = 15,
    max_size: int | None = 80,
    min_tm: float = 0,
    max_tm: float | None = None,
    max_distance: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the bridging oligo that joins two fragments in a ligase cycling reaction.

    The oligo must occur in ``left_sequence + right_sequence`` across the junction, with both halves
    accepted by the homology checker. Returns ``found`` and, when found, the oligo name and the span
    ``[start, end)`` of the oligo in the concatenation (``find_oligo_homology``). The first matching oligo in
    name order is returned, as in DNA Cauldron.
    """

    left = _named_dna({"left": left_sequence}, "left_sequence")["left"]
    right = _named_dna({"right": right_sequence}, "right_sequence")["right"]
    mix = _build_mix(
        "ligase_cycling",
        {"left": left, "right": right},
        enzymes=[],
        topology="linear",
        filter_site_fragments=False,
        bridging_oligos=bridging_oligos,
        min_size=min_size,
        max_size=max_size,
        min_tm=min_tm,
        max_tm=max_tm,
        max_distance=max_distance,
    )
    by_id = {fragment.original_part.id: fragment for fragment in mix.fragments}
    hit = _run(lambda: mix.find_oligo_homology(by_id["left"], by_id["right"]), "DNA Cauldron could not search the oligos")
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "find_dnacauldron_bridging_oligo",
        "found": hit is not None,
        "provenance": _provenance(method="LigaseCyclingReactionMix.find_oligo_homology"),
    }
    if hit is not None:
        oligo, (start, end) = hit
        result.update({"oligo": str(oligo.id), "start": int(start), "end": int(end)})
    return result


# ---------------------------------------------------------------------------------------------
# Record-level helpers, cut plot
# ---------------------------------------------------------------------------------------------


def plot_dnacauldron_cuts(
    record: dict[str, Any],
    enzymes: list[str],
    output_path: str,
    *,
    linear: str = "auto",
    figure_width: float = 5.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a construct map that marks where restriction enzymes cut, to a new PNG file.

    ``record`` is a record dictionary; ``linear`` is ``auto`` (from the record topology), ``true`` or
    ``false``. Each cut position of each enzyme (Biopython ``search``, one-based) is added as a feature and
    drawn on a linear or circular map (``plot_cuts``). The result gives the PNG path and the cut positions.
    """

    target = _new_file(output_path, "output_path", (".png",))
    module = require("dnacauldron.Assembly.AssemblyReportWriter.plot_cuts", _EXTRA)
    names = _enzyme_names(enzymes)
    mode = _choice(linear, "linear", ("auto", "true", "false"))
    width = _number(figure_width, "figure_width", minimum=1, maximum=40)
    rec = _record_from_spec(record)
    is_linear = rec.annotations.get("topology", "linear") == "linear" if mode == "auto" else mode == "true"
    restriction = require("Bio.Restriction", _EXTRA)
    cuts = {item: [int(c) for c in getattr(restriction, item).search(rec.seq, linear=is_linear)] for item in names}
    pyplot = _with_pyplot()
    with _quiet_fonts():
        result = _run(
            lambda: module.plot_cuts(rec, names, linear=is_linear, figure_width=width),
            "DNA Cauldron could not plot the cuts",
        )
        axis = result[0] if isinstance(result, tuple) else result
        _save_png(axis.figure, target, pyplot)
    return {
        "status": "completed",
        "operation": "plot_dnacauldron_cuts",
        "path": str(target),
        "data": {"linear": is_linear, "length": len(rec.seq), "cuts": cuts},
        "provenance": _provenance(method="dnacauldron.Assembly.AssemblyReportWriter.plot_cuts.plot_cuts"),
    }


def list_dnacauldron_record_overhangs(
    record: dict[str, Any] | None = None,
    *,
    genbank_path: str | None = None,
    with_locations: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the overhang sequences recorded as annotations in a construct record.

    Give a ``record`` dictionary or a GenBank file (``genbank_path``). Four-base ``homology`` features
    (labels written when a construct is assembled with homologies annotated) are returned in position order;
    when there are none, the 4-base gaps between ``From <part>`` features are used
    (``list_overhangs_from_record_annotations``). ``with_locations`` also returns each start position.
    """

    utils = require("dnacauldron.utils", _EXTRA)
    if (record is None) == (genbank_path is None):
        raise ConfigurationError("Provide exactly one of record and genbank_path.")
    if record is not None:
        rec = _record_from_spec(record)
    else:
        rec = _run(
            lambda: require("dnacauldron.biotools", _EXTRA).load_record(str(_existing_file(genbank_path, "genbank_path"))),
            "DNA Cauldron could not read the file",
        )
    found = _run(
        lambda: utils.list_overhangs_from_record_annotations(rec, with_locations=_flag(with_locations, "with_locations")),
        "DNA Cauldron could not read the overhangs",
    )
    rows: list[Any] = [[int(start), str(o)] for start, o in found] if with_locations else [str(o) for o in found]
    return {
        "status": "completed",
        "operation": "list_dnacauldron_record_overhangs",
        "overhang_count": len(rows),
        "overhangs": rows,
        "provenance": _provenance(method="dnacauldron.utils.list_overhangs_from_record_annotations"),
    }


def get_dnacauldron_insert_overhangs(
    sequence: str,
    *,
    enzyme: str = "BsmBI",
    standardize: bool = True,
    topology: str = "circular",
    verbose: bool = True,
) -> dict[str, Any]:
    """Get the overhang pair of the insert of a Type-2S part by digesting it with an enzyme.

    The part is digested into sticky-end fragments; the one fragment without a site is the insert and its
    ``(left, right)`` overhangs are returned. With ``standardize`` the pair is the smaller of ``(o1, o2)`` and
    ``(revcomp(o2), revcomp(o1))`` so that a part and its reverse-complement give the same pair
    (``get_overhangs_from_record``).
    """

    module = require("dnacauldron.utils.insert_parts_on_backbones", _EXTRA)
    name = _enzyme_names(enzyme)[0]
    rec = _make_record("part", _dna(sequence, "sequence"), _choice(topology, "topology", ("linear", "circular")))
    pair = _run(
        lambda: module.get_overhangs_from_record(rec, enzyme=name, standardize=_flag(standardize, "standardize")),
        "DNA Cauldron could not find the insert overhangs (the part needs exactly one site-free fragment)",
    )
    return {
        "status": "completed",
        "operation": "get_dnacauldron_insert_overhangs",
        "enzyme": name,
        "standardized": bool(standardize),
        "left_overhang": str(pair[0]),
        "right_overhang": str(pair[1]),
        "provenance": _provenance(method="dnacauldron.utils.insert_parts_on_backbones.get_overhangs_from_record"),
    }


# ---------------------------------------------------------------------------------------------
# Oligo annealing, BASIC assembly, plans
# ---------------------------------------------------------------------------------------------


def _simulate_plan_in_memory(records: Mapping[str, str], assemblies: list[Any], dc: Any) -> Any:
    repository = dc.SequenceRepository()
    for key, text in records.items():
        repository.add_record((key, text))
    plan = dc.AssemblyPlan(assemblies, logger=None)
    return _run(lambda: plan.simulate(repository), "DNA Cauldron could not simulate the plan")


def _flaw_rows(simulation: Any) -> list[dict[str, Any]]:
    return [{"message": str(flaw.message), "suggestion": str(getattr(flaw, "suggestion", "")), "data": _plain_data(flaw.data)} for flaw in simulation.errors]


def _plain_data(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _plain_data(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_plain_data(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def anneal_dnacauldron_oligo_pair(
    oligo: str,
    other_oligo: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate the annealing of an oligo pair into a sticky-end fragment (``OligoPairAnnealing``).

    Both oligos are written 5' to 3'. One oligo must contain the reverse complement of the other; the
    annealed part is the double-stranded ``sequence`` and the unpaired bases on each side are the sticky
    ends (strand 1). Both orders are tried, so up to two products are returned; none raises
    ``ConfigurationError``. This is the adapter-making step of BASIC assembly.
    """

    dc = require("dnacauldron", _EXTRA)
    records = {"oligo_1": _dna(oligo, "oligo", maximum=2000), "oligo_2": _dna(other_oligo, "other_oligo", maximum=2000)}
    assembly = dc.OligoPairAnnealing(parts=["oligo_1", "oligo_2"], name="annealed")
    repository = dc.SequenceRepository()
    for key, text in records.items():
        repository.add_record((key, text))
    simulation = _run(lambda: assembly.simulate(repository), "DNA Cauldron could not anneal the oligos")
    if simulation.errors:
        raise ConfigurationError("; ".join(str(flaw.message) for flaw in simulation.errors))
    return {
        "status": "completed",
        "operation": "anneal_dnacauldron_oligo_pair",
        "product_count": len(simulation.construct_records),
        "products": [_sticky_dict(item.seq) for item in simulation.construct_records],
        "provenance": _provenance(method="OligoPairAnnealing.simulate"),
    }


def simulate_dnacauldron_basic_assembly(
    parts: dict[str, str],
    adapter_oligos: dict[str, list[str]],
    assembly_parts: list[str],
    *,
    max_constructs: int = 40,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate a BASIC assembly: BsaI adapters on each part, then sticky-end circular assembly.

    ``parts`` are the linear DNA parts (with BsaI sites); ``adapter_oligos`` maps each adapter name to its
    two oligos, which are annealed first (``OligoPairAnnealing``); ``assembly_parts`` is the ordered list of
    triplets ``adapter, part, adapter`` (a multiple of 3 names). Each triplet is ligated through BsaI
    overhangs and the adapted parts are assembled into circular constructs through their long linker ends.
    Up to ``max_constructs`` constructs are returned; assembly flaws are listed in ``errors``.
    """

    dc = require("dnacauldron", _EXTRA)
    records = _named_dna(parts, "parts")
    if not isinstance(adapter_oligos, Mapping) or not 1 <= len(adapter_oligos) <= 100:
        raise ConfigurationError("adapter_oligos must map 1 to 100 adapter names to a pair of oligos.")
    assemblies: list[Any] = []
    for adapter, pair in adapter_oligos.items():
        label = _name(adapter, "adapter name")
        if isinstance(pair, (str, bytes)) or not isinstance(pair, Sequence) or len(pair) != 2:
            raise ConfigurationError(f"adapter_oligos[{label!r}] must be a list of two oligos.")
        for index, text in enumerate(pair, start=1):
            records[f"{label}__oligo_{index}"] = _dna(text, f"adapter_oligos[{label!r}][{index - 1}]", maximum=2000)
        assemblies.append(dc.OligoPairAnnealing(parts=[f"{label}__oligo_1", f"{label}__oligo_2"], name=label))
    if isinstance(assembly_parts, (str, bytes)) or not isinstance(assembly_parts, Sequence) or len(assembly_parts) % 3 or not assembly_parts:
        raise ConfigurationError("assembly_parts must list adapter, part, adapter triplets (a multiple of 3 names).")
    ordered = [_name(item, "assembly_parts item") for item in assembly_parts]
    unknown = [item for item in ordered if item not in records and item not in adapter_oligos]
    if unknown:
        raise ConfigurationError(f"assembly_parts names are neither parts nor adapters: {unknown}")
    assemblies.append(
        dc.BASICAssembly(
            parts=ordered,
            name="basic_assembly",
            max_constructs=_integer(max_constructs, "max_constructs", minimum=1, maximum=1000),
            expected_constructs="any_number",
        )
    )
    simulation = _simulate_plan_in_memory(records, assemblies, dc)
    basic = simulation.assembly_simulations[-1]
    constructs = [
        {**_construct_row(item), "parts": [_label_part(f) for f in item.fragments]}
        for item in basic.construct_records
    ]
    return {
        "status": "completed",
        "operation": "simulate_dnacauldron_basic_assembly",
        "adapters_annealed": sorted(adapter_oligos),
        "construct_count": len(constructs),
        "constructs": constructs,
        "errors": _flaw_rows(basic),
        "outcomes": {k: int(v) for k, v in sorted(simulation.compute_stats().items())},
        "provenance": _provenance(method="BASICAssembly.simulate (adapters via OligoPairAnnealing)"),
    }


def _label_part(fragment: Any) -> str:
    part = fragment.original_part.id if hasattr(fragment, "original_part") else fragment.id
    return f"{part}{'_r' if getattr(fragment, 'is_reversed', False) else ''}"


def read_dnacauldron_assembly_plan(
    plan_path: str,
    *,
    sheet_name: str = "all",
    header: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import an assembly plan from a CSV or Excel spreadsheet (``AssemblyPlan.from_spreadsheet``).

    Each row is ``assembly name, class: <type>, part, part, ...`` with the class one of
    ``type2s_assembly``, ``gibson_assembly``, ``BASIC_assembly``, ``biobrick_assembly``, ``oligo_annealing``
    or ``lcr_assembly``; other cells may be ``parameter: value``. The rows of a CSV must have the same
    number of columns (pad with empty cells). The result lists each assembly with its class, parts and
    the dependency level that DNA Cauldron computes from the shared names.
    """

    dc = require("dnacauldron", _EXTRA)
    path = _existing_file(plan_path, "plan_path")
    if path.stat().st_size > 20_000_000:
        raise ConfigurationError("The plan file is larger than 20 MB.")
    plan = _run(
        lambda: dc.AssemblyPlan.from_spreadsheet(
            path=str(path),
            sheet_name=_name(sheet_name, "sheet_name"),
            header=0 if _flag(header, "header") else None,
            logger=None,
        ),
        "DNA Cauldron could not read the plan",
    )
    return {
        "status": "completed",
        "operation": "read_dnacauldron_assembly_plan",
        "plan_name": str(plan.name),
        "assembly_count": len(plan.assemblies),
        "assemblies": [
            {
                "name": str(item.name),
                "class": type(item).__name__,
                "parts": [str(p) for p in item.parts],
                "level": int(item.dependencies["level"]),
                "depends_on": sorted(str(x) for x in item.dependencies["depends_on"]),
                "used_in": sorted(str(x) for x in item.dependencies["used_in"]),
            }
            for item in plan.assemblies
        ],
        "provenance": _provenance(method="AssemblyPlan.from_spreadsheet"),
    }


def write_dnacauldron_assembly_plan(
    assemblies: list[dict[str, Any]],
    *,
    output_path: str | None = None,
    assembly_class: str = "type2s_assembly",
    verbose: bool = True,
) -> dict[str, Any]:
    """Export an assembly plan as a table, and optionally as a CSV file.

    ``assemblies`` lists ``{"name", "parts"}`` (parts are part names or the names of other assemblies).
    The table has one row per assembly, sorted by dependency level then name, with the columns
    ``assembly`` and ``parts`` (``AssemblyPlan.to_dataframe``). With ``output_path`` (a new ``.csv``
    file) the plan is also written in DNA Cauldron's ``construct,parts`` format (``to_spreadsheet``).
    """

    dc = require("dnacauldron", _EXTRA)
    classes = require("dnacauldron.Assembly", _EXTRA).ASSEMBLY_CLASS_DICT
    kind = _choice(assembly_class, "assembly_class", tuple(classes))
    if isinstance(assemblies, (str, bytes)) or not isinstance(assemblies, Sequence) or not 1 <= len(assemblies) <= 1000:
        raise ConfigurationError("assemblies must list 1 to 1000 assemblies.")
    built = []
    for index, item in enumerate(assemblies):
        if not isinstance(item, Mapping) or set(item) - {"name", "parts"}:
            raise ConfigurationError(f"assemblies[{index}] must be {{\"name\", \"parts\"}}.")
        part_names = item.get("parts")
        if isinstance(part_names, (str, bytes)) or not isinstance(part_names, Sequence) or not part_names:
            raise ConfigurationError(f"assemblies[{index}].parts must be a non-empty list.")
        built.append(classes[kind](name=_name(item.get("name"), f"assemblies[{index}].name"), parts=[_name(p, "part") for p in part_names]))
    plan = _run(lambda: dc.AssemblyPlan(built, logger=None), "DNA Cauldron could not build the plan")
    frame = plan.to_dataframe()
    rows = [{"assembly": str(index), "parts": str(row["parts"])} for index, row in frame.iterrows()]
    written = None
    if output_path is not None:
        target = _new_file(output_path, "output_path", (".csv",))
        plan.to_spreadsheet(str(target))
        written = str(target)
    return {
        "status": "completed",
        "operation": "write_dnacauldron_assembly_plan",
        "rows": rows,
        "path": written,
        "provenance": _provenance(method="AssemblyPlan.to_dataframe / to_spreadsheet"),
    }


def tabulate_dnacauldron_plan_constructs(
    plan_path: str,
    sequences_path: str,
    *,
    connectors_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate an assembly plan and tabulate the constructs of every assembly.

    The plan is a CSV or spreadsheet (see ``read_dnacauldron_assembly_plan``); ``sequences_path`` is a folder
    of part files (FASTA, GenBank or SnapGene, the file name is the part name); ``connectors_path`` is an
    optional folder whose sub-folders are connector collections. One row is returned per construct with the
    columns of DNA Cauldron's construct table: construct id, parts, number of parts, size, assembly name,
    dependencies and level. ``outcomes`` counts valid, errored and cancelled assemblies.
    """

    dc = require("dnacauldron", _EXTRA)
    plan_file = _existing_file(plan_path, "plan_path")
    sequences = _existing_dir(sequences_path, "sequences_path")
    repository = dc.SequenceRepository()
    _run(lambda: repository.import_records(folder=str(sequences)), "DNA Cauldron could not read the sequences")
    if connectors_path is not None:
        root = _existing_dir(connectors_path, "connectors_path")
        for sub in sorted(root.iterdir()):
            if sub.is_dir():
                _run(lambda s=sub: repository.import_records(folder=str(s), collection=s.name), "DNA Cauldron could not read the connectors")
    plan = _run(lambda: dc.AssemblyPlan.from_spreadsheet(str(plan_file), logger=None), "DNA Cauldron could not read the plan")
    simulation = _run(lambda: plan.simulate(sequence_repository=repository), "DNA Cauldron could not simulate the plan")
    rows = [_plain_data(item) for item in simulation.compute_all_construct_data_dicts()]
    return {
        "status": "completed",
        "operation": "tabulate_dnacauldron_plan_constructs",
        "construct_count": len(rows),
        "constructs": rows,
        "outcomes": {k: int(v) for k, v in sorted(simulation.compute_stats().items())},
        "provenance": _provenance(method="AssemblyPlanSimulation.compute_all_construct_data_dicts"),
    }


__all__ = [
    "anneal_dnacauldron_oligo_pair",
    "annotate_dnacauldron_record",
    "build_dnacauldron_record",
    "check_dnacauldron_homology",
    "check_dnacauldron_sticky_end_compatibility",
    "compute_dnacauldron_connections_graph",
    "compute_dnacauldron_slots",
    "compute_dnacauldron_slots_graph",
    "digest_dnacauldron_sticky_ends",
    "enumerate_dnacauldron_circular_assemblies",
    "enumerate_dnacauldron_linear_assemblies",
    "find_dnacauldron_bridging_oligo",
    "get_dnacauldron_insert_overhangs",
    "get_dnacauldron_repository_records",
    "list_dnacauldron_record_overhangs",
    "plot_dnacauldron_connections_graph",
    "plot_dnacauldron_cuts",
    "plot_dnacauldron_slots_graph",
    "read_dnacauldron_assembly_plan",
    "select_dnacauldron_connectors",
    "simulate_dnacauldron_basic_assembly",
    "simulate_dnacauldron_lcr_assembly",
    "sticky_fragment_to_dnacauldron_record",
    "tabulate_dnacauldron_plan_constructs",
    "transform_dnacauldron_fragment_chain",
    "transform_dnacauldron_sticky_sequence",
    "write_dnacauldron_assembly_plan",
    "write_dnacauldron_records",
]

_register_public_verbose(globals(), tuple(__all__))
