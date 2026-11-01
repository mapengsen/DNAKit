"""AlphaGenome genome, transcript, track, junction, FASTA, metadata and plot functions (pack P34).

The data classes of the alphagenome package (``Interval``, ``Variant``, ``Junction``,
``Transcript``, ``TrackData``, ``JunctionData``) are built from JSON values and their own methods do
the work, so results agree with alphagenome. Values:

* an interval is ``"chr1:100-200:+"`` or ``{"chromosome", "start", "end", "strand"}`` (0-based,
  half-open; strand ``+``, ``-`` or ``.``);
* a variant is ``"chr1:1024:A>C"`` or ``{"chromosome", "position", "reference_bases",
  "alternate_bases"}`` (1-based position);
* track data is ``{"values": bins x tracks, "names": [...], "strands": [...], "resolution": int,
  "interval": interval}`` (3-D values give a contact map);
* a transcript is ``{"chromosome", "strand", "exons": [[start, end]], "cds": [[start, end(, frame)]],
  "start_codon", "stop_codon", "transcript_id", "gene_id"}``.

alphagenome, numpy, pandas and matplotlib are imported when a function runs. Only
``alphagenome_output_metadata`` and ``alphagenome_scorer_metadata`` use the network.
"""

from __future__ import annotations

import gzip
import importlib.metadata
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-alphagenome"
_PLOT_EXTRA = "bio-plots"
_NUMPY_EXTRA = "bio-analysis"
_PACKAGE = "alphagenome"
_GENOME = "alphagenome.data.genome"
_TRACK = "alphagenome.data.track_data"
_ERRORS = (ValueError, TypeError, KeyError, IndexError, AssertionError, AttributeError)
_STRANDS = ("+", "-", ".")
_ORGANISMS = ("homo_sapiens", "mus_musculus")
_MAX_ITEMS = 100_000
_MAX_CELLS = 50_000_000
_MAX_ROWS = 5_000
_GTF_FEATURES = frozenset({"exon", "CDS", "start_codon", "stop_codon", "Selenocysteine"})
_GTF_ATTRIBUTE = re.compile(r'(\w+)\s+"([^"]*)"')


# Helpers.


def _provenance() -> dict[str, Any]:
    try:
        version: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": _PACKAGE, "package_version": version}


def _done(operation: str, **fields: Any) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": _provenance()}


def _guard(operation: str, call: Any) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except _ERRORS as exc:
        raise ConfigurationError(f"alphagenome could not run {operation}: {exc}") from exc


def _integer(value: object, name: str, minimum: int = -(10**12), maximum: int = 10**12) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return value


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _text(value: object, name: str, maximum: int = 4096) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum} characters.")
    return value.strip()


def _interval(value: object, name: str) -> Any:
    genome = require(_GENOME, _EXTRA)
    if isinstance(value, str):
        return _guard("interval parsing", lambda: genome.Interval.from_str(_text(value, name)))
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be 'chr:start-end:strand' or an object.")
    unknown = set(value) - {"chromosome", "start", "end", "strand", "name"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    return _guard(
        "interval parsing",
        lambda: genome.Interval(
            _text(value.get("chromosome"), f"{name}.chromosome", 256),
            _integer(value.get("start"), f"{name}.start"),
            _integer(value.get("end"), f"{name}.end"),
            _choice(value.get("strand", "."), _STRANDS, f"{name}.strand"),
            str(value.get("name", "")),
        ),
    )


def _intervals(value: object, name: str) -> list[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of intervals.")
    if len(value) > _MAX_ITEMS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_ITEMS} intervals.")
    return [_interval(item, f"{name}[{index}]") for index, item in enumerate(value)]


def _interval_dict(interval: Any) -> dict[str, Any]:
    return {
        "chromosome": str(interval.chromosome),
        "start": int(interval.start),
        "end": int(interval.end),
        "strand": str(interval.strand),
        "width": int(interval.width),
    }


def _variant(value: object, name: str) -> Any:
    genome = require(_GENOME, _EXTRA)
    if isinstance(value, str):
        return _guard("variant parsing", lambda: genome.Variant.from_str(_text(value, name)))
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be 'chr:pos:REF>ALT' or an object.")
    unknown = set(value) - {"chromosome", "position", "reference_bases", "alternate_bases", "name"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    return _guard(
        "variant parsing",
        lambda: genome.Variant(
            _text(value.get("chromosome"), f"{name}.chromosome", 256),
            _integer(value.get("position"), f"{name}.position", 1),
            str(value.get("reference_bases", "")).upper(),
            str(value.get("alternate_bases", "")).upper(),
            str(value.get("name", "")),
        ),
    )


def _variant_dict(variant: Any | None) -> dict[str, Any] | None:
    if variant is None:
        return None
    return {
        "chromosome": str(variant.chromosome),
        "position": int(variant.position),
        "reference_bases": str(variant.reference_bases),
        "alternate_bases": str(variant.alternate_bases),
        "text": str(variant),
    }


def _plain(value: Any) -> Any:
    """JSON-safe copy of numpy, pandas and nested values (NaN becomes None)."""

    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (bool, int, str)):
        return value
    return str(value)


def _pairs(value: object, name: str, *, frame: bool = False) -> list[tuple[int, ...]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of [start, end] pairs.")
    rows: list[tuple[int, ...]] = []
    width = (2, 3) if frame else (2,)
    for index, item in enumerate(value):
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) not in width:
            raise ConfigurationError(f"{name}[{index}] must be {'[start, end(, frame)]' if frame else '[start, end]'}.")
        rows.append(tuple(_integer(x, f"{name}[{index}]") for x in item))
    return rows


# Interval and variant operations.


def alphagenome_adjust_interval(
    interval: str | dict[str, Any],
    *,
    operation: str = "pad",
    start: int = 0,
    end: int = 0,
    use_strand: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pad or shift the boundaries of a genomic interval (alphagenome ``Interval`` methods).

    ``operation`` is ``pad`` (``Interval.pad`` / ``pad_inplace``: ``start`` and ``end`` bases are
    ADDED on each side, negative values trim) or ``boundary_shift`` (``Interval.boundary_shift``:
    the start boundary moves by ``start`` bases and the end boundary by ``end`` bases, so a
    positive ``start`` moves the start to the right). With ``use_strand`` the two sides swap for a
    negative-strand interval, so ``start`` always means the 5' side.
    """

    item = _interval(interval, "interval")
    mode = _choice(operation, ("pad", "boundary_shift"), "operation")
    first, second = _integer(start, "start"), _integer(end, "end")
    flag = _flag(use_strand, "use_strand")
    if mode == "pad":
        new = _guard("Interval.pad", lambda: item.pad(first, second, use_strand=flag))
    else:
        new = _guard("Interval.boundary_shift", lambda: item.boundary_shift(first, second, flag))
    return _done("adjust_interval", method=mode, input=_interval_dict(item), result=_interval_dict(new))


def alphagenome_interval_center(
    interval: str | dict[str, Any], *, use_strand: bool = True, verbose: bool = True
) -> dict[str, Any]:
    """Centre coordinate of an interval (alphagenome ``Interval.center``).

    With ``use_strand`` the rounding of an even width follows the strand, as alphagenome does.
    """

    item = _interval(interval, "interval")
    return _done("interval_center", input=_interval_dict(item),
                 center=int(item.center(_flag(use_strand, "use_strand"))))


def alphagenome_interval_within_reference(
    interval: str | dict[str, Any], *, reference_length: int, verbose: bool = True
) -> dict[str, Any]:
    """Test whether an interval lies inside ``[0, reference_length]`` (``Interval.within_reference``)."""

    item = _interval(interval, "interval")
    length = _integer(reference_length, "reference_length", 1)
    return _done("interval_within_reference", input=_interval_dict(item),
                 reference_length=length, within_reference=bool(item.within_reference(length)))


def alphagenome_interval_strand(
    interval: str | dict[str, Any], *, mode: str = "swap", verbose: bool = True
) -> dict[str, Any]:
    """Flip or remove the strand of an interval.

    ``mode`` is ``swap`` (``Interval.swap_strand``: ``+`` becomes ``-`` and the reverse) or
    ``drop`` (``Interval.as_unstranded``: strand becomes ``.``).
    """

    item = _interval(interval, "interval")
    chosen = _choice(mode, ("swap", "drop"), "mode")
    new = _guard("strand change", lambda: item.swap_strand() if chosen == "swap" else item.as_unstranded())
    return _done("interval_strand", mode=chosen, input=_interval_dict(item), result=_interval_dict(new))


def alphagenome_interval_relation(
    interval: str | dict[str, Any],
    other: str | dict[str, Any],
    *,
    relation: str = "overlaps",
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether an interval overlaps or contains another one (same chromosome required).

    ``relation`` is ``overlaps`` (``Interval.overlaps``: they share at least one base) or
    ``contains`` (``Interval.contains``: ``interval`` completely covers ``other``).
    """

    first, second = _interval(interval, "interval"), _interval(other, "other")
    chosen = _choice(relation, ("overlaps", "contains"), "relation")
    answer = first.overlaps(second) if chosen == "overlaps" else first.contains(second)
    return _done("interval_relation", relation=chosen, interval=_interval_dict(first),
                 other=_interval_dict(second), result=bool(answer))


def alphagenome_interval_coverage(
    interval: str | dict[str, Any],
    intervals: list[Any],
    *,
    bin_size: int = 1,
    stranded: bool = False,
    mask: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage of an interval by other intervals, as a count track or a boolean mask.

    ``stranded=False, mask=False`` is ``Interval.coverage`` (number of intervals per bin);
    ``mask=True`` is ``Interval.binary_mask`` (coverage above zero); ``stranded=True`` is
    ``Interval.coverage_stranded`` / ``binary_mask_stranded``: two columns, same-strand and
    opposite-strand coverage relative to ``interval``. ``bin_size`` must divide the interval width;
    counts are summed within a bin.
    """

    item = _interval(interval, "interval")
    others = _intervals(intervals, "intervals")
    size = _integer(bin_size, "bin_size", 1)
    by_strand, as_mask = _flag(stranded, "stranded"), _flag(mask, "mask")
    if item.width > _MAX_CELLS:
        raise ConfigurationError("interval is too wide for a coverage track.")
    if by_strand:
        method = item.binary_mask_stranded if as_mask else item.coverage_stranded
        track = _guard("coverage_stranded", lambda: method(others, bin_size=size))
    else:
        method = item.binary_mask if as_mask else item.coverage
        track = _guard("coverage", lambda: method(others, bin_size=size))
    return _done("interval_coverage", interval=_interval_dict(item), bin_size=size,
                 stranded=by_strand, mask=as_mask, shape=list(track.shape), values=_plain(track))


def alphagenome_overlap_ranges(
    interval: str | dict[str, Any], intervals: list[Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Relative sub-ranges of an interval that other intervals overlap (``Interval.overlap_ranges``).

    Each range is ``[start, end]`` relative to the start of ``interval``, clipped to it.
    """

    item = _interval(interval, "interval")
    others = _intervals(intervals, "intervals")
    ranges = _guard("overlap_ranges", lambda: item.overlap_ranges(others))
    return _done("overlap_ranges", interval=_interval_dict(item), ranges=_plain(ranges))


def alphagenome_variant_overlap(
    variant: str | dict[str, Any],
    interval: str | dict[str, Any],
    *,
    allele: str = "reference",
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether the reference or alternate allele of a variant overlaps an interval.

    ``allele`` is ``reference`` (``Variant.reference_overlaps``) or ``alternate``
    (``Variant.alternate_overlaps``: the span the ALT bases would occupy from the variant start).
    """

    var, item = _variant(variant, "variant"), _interval(interval, "interval")
    chosen = _choice(allele, ("reference", "alternate"), "allele")
    answer = var.reference_overlaps(item) if chosen == "reference" else var.alternate_overlaps(item)
    return _done("variant_overlap", allele=chosen, variant=_variant_dict(var),
                 interval=_interval_dict(item), overlaps=bool(answer))


def alphagenome_split_variant(
    variant: str | dict[str, Any], *, anchor: int, verbose: bool = True
) -> dict[str, Any]:
    """Split a variant into an upstream and a downstream part at a 0-based anchor (``Variant.split``).

    If the anchor lies outside the reference span the whole variant is on one side and the other
    part is None.
    """

    var = _variant(variant, "variant")
    point = _integer(anchor, "anchor", 0)
    upstream, downstream = _guard("Variant.split", lambda: var.split(point))
    return _done("split_variant", anchor=point, variant=_variant_dict(var),
                 upstream=_variant_dict(upstream), downstream=_variant_dict(downstream))


def alphagenome_junction_regions(
    junction: str | dict[str, Any],
    *,
    region: str = "dinucleotide",
    overhang: list[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Regions around the splice sites of a junction.

    ``region`` is ``dinucleotide`` (``Junction.dinucleotide_region``: the first and last two bases
    of the junction), ``acceptor`` (``Junction.acceptor_region``) or ``donor``
    (``Junction.donor_region``); the last two are the site padded by ``overhang`` = [upstream,
    downstream] bases (default [250, 250]). The junction must be stranded (``+`` or ``-``).
    """

    genome = require(_GENOME, _EXTRA)
    item = _interval(junction, "junction")
    chosen = _choice(region, ("dinucleotide", "acceptor", "donor"), "region")
    pad = (250, 250) if overhang is None else tuple(_pairs([overhang], "overhang")[0])

    def run() -> Any:
        junc = genome.Junction(item.chromosome, item.start, item.end, item.strand)
        if chosen == "dinucleotide":
            return [_interval_dict(x) for x in junc.dinucleotide_region()], junc
        method = junc.acceptor_region if chosen == "acceptor" else junc.donor_region
        return [_interval_dict(method(pad))], junc

    regions, junc = _guard("Junction regions", run)
    return _done("junction_regions", region=chosen, junction=_interval_dict(junc),
                 acceptor=int(junc.acceptor), donor=int(junc.donor), regions=regions)


# Transcripts.


def _transcript(value: object, name: str) -> Any:
    genome = require(_GENOME, _EXTRA)
    transcript_module = require("alphagenome.data.transcript", _EXTRA)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be an object.")
    unknown = set(value) - {"chromosome", "strand", "exons", "cds", "start_codon", "stop_codon",
                            "transcript_id", "gene_id", "protein_id"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    chromosome = _text(value.get("chromosome"), f"{name}.chromosome", 256)
    strand = _choice(value.get("strand"), ("+", "-"), f"{name}.strand")

    def build(key: str, frame: bool) -> list[Any] | None:
        if value.get(key) is None:
            return None
        out = []
        for row in _pairs(value[key], f"{name}.{key}", frame=frame):
            item = genome.Interval(chromosome, row[0], row[1], strand)
            if len(row) == 3:
                item.info["frame"] = row[2]
            out.append(item)
        return out

    exons = build("exons", False)
    if not exons:
        raise ConfigurationError(f"{name}.exons must hold at least one exon.")
    return _guard("transcript building", lambda: transcript_module.Transcript(
        exons=exons,
        cds=build("cds", True),
        start_codon=build("start_codon", False),
        stop_codon=build("stop_codon", True),
        transcript_id=value.get("transcript_id"),
        gene_id=value.get("gene_id"),
        protein_id=value.get("protein_id"),
    ))


def _transcript_dict(transcript: Any) -> dict[str, Any]:
    def spans(items: Any) -> list[list[int]] | None:
        if items is None:
            return None
        return [
            [int(x.start), int(x.end)] + ([int(x.info["frame"])] if "frame" in x.info else [])
            for x in items
        ]

    return {
        "transcript_id": transcript.transcript_id,
        "gene_id": transcript.gene_id,
        "protein_id": transcript.protein_id,
        "chromosome": transcript.chromosome,
        "strand": transcript.strand,
        "exons": spans(transcript.exons),
        "cds": spans(transcript.cds),
        "start_codon": spans(transcript.start_codon),
        "stop_codon": spans(transcript.stop_codon),
        "is_coding": bool(transcript.is_coding),
        "truncation_fixed": bool(transcript.info.get("truncation_fixed", False)),
    }


def alphagenome_transcript_cds_offset(
    transcript: dict[str, Any], *, position: int, verbose: bool = True
) -> dict[str, Any]:
    """Offset of a 0-based genome position within the coding sequence (``Transcript.offset_in_cds``).

    The offset counts coding bases (CDS plus stop codon) from the start of the CDS along the
    transcript direction, so it accounts for strand and introns. None when the position is not in
    the CDS.
    """

    item = _transcript(transcript, "transcript")
    point = _integer(position, "position", 0)
    offset = _guard("offset_in_cds", lambda: item.offset_in_cds(point))
    return _done("transcript_cds_offset", position=point, offset=None if offset is None else int(offset))


def alphagenome_fix_transcript_truncation(
    transcript: dict[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Repair the CDS frame of a transcript whose coding sequence is cut (``Transcript.fix_truncation``).

    Each CDS and stop-codon entry needs a frame: ``[start, end, frame]``. The result has the CDS
    trimmed to whole codons and new start and stop codon intervals; a transcript with 6 or fewer
    coding bases is returned unchanged.
    """

    item = _transcript(transcript, "transcript")
    transcript_module = require("alphagenome.data.transcript", _EXTRA)
    fixed = _guard("fix_truncation", lambda: transcript_module.Transcript.fix_truncation(item))
    return _done("fix_transcript_truncation", input=_transcript_dict(item), result=_transcript_dict(fixed))


def alphagenome_transcripts_from_gtf(
    gtf_path: str,
    *,
    transcript_ids: list[str],
    fix_truncation: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build transcript models from the rows of a GENCODE GTF (``Transcript.from_gtf_df``).

    ``gtf_path`` is a GTF or GTF.gz file; ``transcript_ids`` lists the transcripts to build (their
    exon, CDS, start_codon, stop_codon and Selenocysteine rows are used). With ``fix_truncation`` the
    CDS frame is repaired as ``alphagenome_fix_transcript_truncation`` does. Coordinates are
    0-based half-open (GTF start minus one).
    """

    pandas = require("pandas", _EXTRA)
    transcript_module = require("alphagenome.data.transcript", _EXTRA)
    path = input_path(gtf_path)
    if isinstance(transcript_ids, (str, bytes)) or not isinstance(transcript_ids, Sequence) or not transcript_ids:
        raise ConfigurationError("transcript_ids must be a non-empty list.")
    wanted = {_text(item, "transcript_ids[]", 256) for item in transcript_ids}
    repair = _flag(fix_truncation, "fix_truncation")
    opener = gzip.open if path.suffix == ".gz" else open
    rows: dict[str, list[dict[str, Any]]] = {name: [] for name in wanted}
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ConfigurationError(f"GTF line {number} has fewer than 9 columns.")
            if fields[2] not in _GTF_FEATURES:
                continue
            attributes = dict(_GTF_ATTRIBUTE.findall(fields[8]))
            if attributes.get("transcript_id") not in wanted:
                continue
            try:
                row = {
                    "Chromosome": fields[0], "Start": int(fields[3]) - 1, "End": int(fields[4]),
                    "Strand": fields[6], "Feature": fields[2],
                    "Frame": int(fields[7]) if fields[7] != "." else 0,
                    **{key: attributes.get(key) for key in (
                        "gene_id", "transcript_id", "protein_id", "gene_type", "transcript_type")},
                }
            except ValueError as exc:
                raise ConfigurationError(f"GTF line {number} has invalid coordinates.") from exc
            rows[attributes["transcript_id"]].append(row)
    built: list[dict[str, Any]] = []
    missing: list[str] = []
    for name in sorted(rows):
        if not rows[name]:
            missing.append(name)
            continue
        frame = pandas.DataFrame(rows[name])
        transcript = _guard("from_gtf_df", lambda frame=frame: transcript_module.Transcript.from_gtf_df(
            frame, fix_truncation=repair))
        built.append(_transcript_dict(transcript))
    return _done("transcripts_from_gtf", fix_truncation=repair, transcript_count=len(built),
                 missing=missing, transcripts=built)


# Track data.


def _track_data(value: object, name: str) -> Any:
    np = require("numpy", _NUMPY_EXTRA)
    pandas = require("pandas", _EXTRA)
    module = require(_TRACK, _EXTRA)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be an object with values and names.")
    unknown = set(value) - {"values", "names", "strands", "resolution", "interval"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    try:
        values = np.asarray(value.get("values"), dtype=np.float32)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"{name}.values must be a numeric array.") from exc
    if values.ndim < 1 or values.size == 0 or values.size > _MAX_CELLS:
        raise ConfigurationError(f"{name}.values must be a non-empty numeric array.")
    names = value.get("names")
    if isinstance(names, (str, bytes)) or not isinstance(names, Sequence):
        raise ConfigurationError(f"{name}.names must be a list of track names.")
    strands = value.get("strands") or ["."] * len(names)
    if len(strands) != len(names):
        raise ConfigurationError(f"{name}.strands must have one entry per name.")
    resolution = _integer(value.get("resolution", 1), f"{name}.resolution", 1)
    interval = None if value.get("interval") is None else _interval(value["interval"], f"{name}.interval")
    metadata = pandas.DataFrame({"name": [str(x) for x in names],
                                 "strand": [_choice(s, _STRANDS, f"{name}.strands[]") for s in strands]})
    return _guard("TrackData", lambda: module.TrackData(
        values=values, metadata=metadata, resolution=resolution, interval=interval))


def _track_dict(track: Any) -> dict[str, Any]:
    return {
        "values": _plain(track.values),
        "shape": list(track.values.shape),
        "names": [str(x) for x in track.names],
        "strands": [str(x) for x in track.strands],
        "resolution": int(track.resolution),
        "interval": None if track.interval is None else _interval_dict(track.interval),
    }


def alphagenome_slice_track_data(
    track: dict[str, Any],
    *,
    interval: str | dict[str, Any] | None = None,
    start: int | None = None,
    end: int | None = None,
    match_resolution: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut track data down to a genomic region.

    Give either ``interval`` (``TrackData.slice_by_interval``: genome coordinates, which the track's
    own interval must fully contain; ``match_resolution`` widens it to whole bins) or ``start`` and
    ``end`` (``TrackData.slice_by_positions``: 0-based base positions relative to the track start).
    """

    data = _track_data(track, "track")
    if interval is not None:
        if start is not None or end is not None:
            raise ConfigurationError("give either interval or start and end, not both.")
        region = _interval(interval, "interval")
        new = _guard("slice_by_interval", lambda: data.slice_by_interval(region, _flag(match_resolution, "match_resolution")))
    else:
        if start is None or end is None:
            raise ConfigurationError("give interval, or both start and end.")
        low, high = _integer(start, "start", 0), _integer(end, "end", 0)
        new = _guard("slice_by_positions", lambda: data.slice_by_positions(low, high))
    return _done("slice_track_data", input=_track_dict(data), result=_track_dict(new))


def alphagenome_track_bin_index(
    track: dict[str, Any], *, relative_position: int, verbose: bool = True
) -> dict[str, Any]:
    """Bin index of a position relative to the track start (``TrackData.bin_index``): position // resolution."""

    data = _track_data(track, "track")
    position = _integer(relative_position, "relative_position", 0)
    return _done("track_bin_index", resolution=int(data.resolution), relative_position=position,
                 bin_index=int(data.bin_index(position)))


def alphagenome_resize_track_data(
    track: dict[str, Any], *, width: int, verbose: bool = True
) -> dict[str, Any]:
    """Crop or zero-pad track data to ``width`` bases around its centre.

    This is ``TrackData.resize`` (the model-output ``Output.resize`` does the same on its data).
    ``width`` must be divisible by the resolution when growing.
    """

    data = _track_data(track, "track")
    size = _integer(width, "width", 1, _MAX_CELLS)
    new = _guard("resize", lambda: data.resize(size))
    return _done("resize_track_data", width=size, input=_track_dict(data), result=_track_dict(new))


def alphagenome_change_track_resolution(
    track: dict[str, Any],
    *,
    resolution: int,
    aggregation: str = "sum",
    verbose: bool = True,
) -> dict[str, Any]:
    """Change the resolution of track data, finer (upsample) or coarser (downsample).

    A target below the current resolution is ``TrackData.upsample`` (values repeated; with ``sum``
    they are divided by the repeat so the total is kept), a target above it is
    ``TrackData.downsample`` (values pooled with ``sum`` or ``max``); the target must divide or be
    divisible by the current resolution. The result reports which direction ran.
    """

    data = _track_data(track, "track")
    target = _integer(resolution, "resolution", 1)
    rule = _choice(aggregation, ("sum", "max"), "aggregation")
    module = require(_TRACK, _EXTRA)
    kind = module.AggregationType.SUM if rule == "sum" else module.AggregationType.MAX
    direction = "same" if target == data.resolution else ("downsample" if target > data.resolution else "upsample")
    method = data.downsample if direction == "downsample" else data.upsample
    new = _guard(direction, lambda: method(target, kind))
    return _done("change_track_resolution", direction=direction, aggregation=rule,
                 input=_track_dict(data), result=_track_dict(new))


def alphagenome_reverse_complement_track_data(
    track: dict[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Reverse-complement track data (``TrackData.reverse_complement``).

    Positions are reversed and the strands of stranded tracks swap; the interval, if present, must
    be stranded (it is flipped too).
    """

    data = _track_data(track, "track")
    new = _guard("reverse_complement", lambda: data.reverse_complement())
    return _done("reverse_complement_track_data", input=_track_dict(data), result=_track_dict(new))


# Junction data.


def _junction_data(value: object, name: str) -> Any:
    np = require("numpy", _NUMPY_EXTRA)
    pandas = require("pandas", _EXTRA)
    genome = require(_GENOME, _EXTRA)
    module = require("alphagenome.data.junction_data", _EXTRA)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be an object with junctions, values and names.")
    unknown = set(value) - {"junctions", "values", "names", "interval"}
    if unknown:
        raise ConfigurationError(f"{name} has unsupported keys: {sorted(unknown)}.")
    raw = value.get("junctions")
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence):
        raise ConfigurationError(f"{name}.junctions must be a list of junctions.")
    junctions = [_interval(item, f"{name}.junctions[{i}]") for i, item in enumerate(raw)]
    try:
        values = np.asarray(value.get("values"), dtype=np.float32)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"{name}.values must be a numeric matrix.") from exc
    names = value.get("names")
    if isinstance(names, (str, bytes)) or not isinstance(names, Sequence):
        raise ConfigurationError(f"{name}.names must be a list of track names.")
    if values.ndim != 2 or values.shape != (len(junctions), len(names)):
        raise ConfigurationError(f"{name}.values must have one row per junction and one column per name.")
    array = np.empty(len(junctions), dtype=object)
    for index, item in enumerate(junctions):
        array[index] = _guard(
            "Junction", lambda item=item: genome.Junction(item.chromosome, item.start, item.end, item.strand)
        )
    interval = None if value.get("interval") is None else _interval(value["interval"], f"{name}.interval")
    return _guard("JunctionData", lambda: module.JunctionData(
        junctions=array, values=values,
        metadata=pandas.DataFrame({"name": [str(x) for x in names]}), interval=interval))


def _junction_dict(data: Any) -> dict[str, Any]:
    return {
        "junctions": [_interval_dict(j) for j in data.junctions],
        "names": [str(x) for x in data.names],
        "values": _plain(data.values),
    }


def alphagenome_normalize_junction_counts(
    junction_data: dict[str, Any], *, total_k: float = 10.0, verbose: bool = True
) -> dict[str, Any]:
    """Scale splice-junction read counts so that all values sum to ``total_k`` (``JunctionData.normalize_values``).

    ``junction_data`` is ``{"junctions": [intervals], "values": junctions x tracks, "names": [...]}``.
    """

    data = _junction_data(junction_data, "junction_data")
    scale = _number(total_k, "total_k")
    new = _guard("normalize_values", lambda: data.normalize_values(scale))
    return _done("normalize_junction_counts", total_k=scale, result=_junction_dict(new))


def alphagenome_junctions_in_interval(
    junction_data: dict[str, Any], *, interval: str | dict[str, Any], verbose: bool = True
) -> dict[str, Any]:
    """Keep the splice junctions that overlap an interval (``JunctionData.intersect_with_interval``)."""

    data = _junction_data(junction_data, "junction_data")
    region = _interval(interval, "interval")
    new = _guard("intersect_with_interval", lambda: data.intersect_with_interval(region))
    return _done("junctions_in_interval", interval=_interval_dict(region),
                 input_count=len(data), result_count=len(new), result=_junction_dict(new))


# FASTA and metadata.


def alphagenome_extract_fasta_interval(
    fasta_path: str, *, interval: str | dict[str, Any], verbose: bool = True
) -> dict[str, Any]:
    """Read the reference sequence of an interval from an indexed FASTA (alphagenome ``FastaExtractor``).

    The FASTA needs a ``.fai`` index. Bases outside the chromosome are padded with ``N``, the
    sequence is upper-case, and a negative-strand interval returns the reverse complement.
    """

    fasta = require("alphagenome.io.fasta", _EXTRA)
    path = input_path(fasta_path)
    region = _interval(interval, "interval")
    sequence = _guard("FastaExtractor.extract", lambda: fasta.FastaExtractor(str(path)).extract(region))
    return _done("extract_fasta_interval", interval=_interval_dict(region), sequence=sequence,
                 length=len(sequence))


def _metadata_rows(frame: Any, limit: int) -> dict[str, Any]:
    records = frame.reset_index().head(limit).to_dict("records")
    return {"row_count": int(len(frame)), "columns": [str(c) for c in frame.columns],
            "rows": [_plain(r) for r in records]}


def alphagenome_output_metadata(
    api_key: str,
    *,
    organism: str = "homo_sapiens",
    max_rows: int = 50,
    timeout: float | None = None,
    address: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the output tracks of the hosted AlphaGenome DNA model (``DnaClient.output_metadata``).

    Network call with an API key (sent only to the service). For each output type (atac, cage,
    dnase, rna_seq, chip_histone, chip_tf, splice_sites, splice_site_usage, splice_junctions,
    contact_maps, procap) the result gives the number of tracks, the columns and the first
    ``max_rows`` rows. ``organism`` is ``homo_sapiens`` or ``mus_musculus``.
    """

    client_module = require("alphagenome.models.dna_client", _EXTRA)
    key = _text(api_key, "api_key", 512)
    species = _choice(organism, _ORGANISMS, "organism")
    limit = _integer(max_rows, "max_rows", 0, _MAX_ROWS)
    wait = None if timeout is None else _number(timeout, "timeout")

    def run() -> dict[str, Any]:
        client = client_module.create(key, timeout=wait, address=address)
        meta = client.output_metadata(client_module.Organism[species.upper()])
        outputs: dict[str, Any] = {}
        for field in meta.__dataclass_fields__:
            frame = getattr(meta, field)
            if frame is not None:
                outputs[field] = _metadata_rows(frame, limit)
        return outputs

    outputs = _guard("DnaClient.output_metadata", run)
    return _done("output_metadata", organism=species, source="AlphaGenome service", outputs=outputs)


def alphagenome_scorer_metadata(
    api_key: str,
    *,
    max_rows: int = 50,
    timeout: float | None = None,
    address: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the variant scorers of the AlphaGenome Atlas service (``AtlasClient.scorer_metadata``).

    Network call with an API key. For each scorer: its name, whether scores are signed, and the
    track metadata (row count, columns, first ``max_rows`` rows).
    """

    atlas = require("alphagenome.atlas.atlas", _EXTRA)
    key = _text(api_key, "api_key", 512)
    limit = _integer(max_rows, "max_rows", 0, _MAX_ROWS)
    wait = None if timeout is None else _number(timeout, "timeout")

    def run() -> dict[str, Any]:
        client = atlas.create(key, timeout=wait, address=address)
        return {
            name: {"is_signed": bool(item.is_signed), "track_metadata": _metadata_rows(item.track_metadata, limit)}
            for name, item in client.scorer_metadata().items()
        }

    scorers = _guard("AtlasClient.scorer_metadata", run)
    return _done("scorer_metadata", source="AlphaGenome Atlas service", scorer_count=len(scorers), scorers=scorers)


# Plots.


def _png(output_path: object) -> Path:
    if not isinstance(output_path, str) or not output_path.strip():
        raise ConfigurationError("output_path must be a non-empty string.")
    target = Path(output_path).expanduser().resolve()
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _save(target: Path, draw: Any) -> None:
    matplotlib = require("matplotlib", _PLOT_EXTRA)
    matplotlib.use("Agg")
    plt = require("matplotlib.pyplot", _PLOT_EXTRA)
    before = set(plt.get_fignums())
    try:
        figure = _guard("plot", draw)
        figure.savefig(target, format="png", dpi=120)
    finally:
        for number in plt.get_fignums():
            if number not in before:
                plt.close(number)


def alphagenome_plot_contact_map_difference(
    contact_maps: dict[str, Any],
    output_path: str,
    *,
    vmin: float = -1.0,
    vmax: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the difference between two contact maps (ALT minus REF) with a red-blue scale centred on 0.

    This is alphagenome's ``ContactMapsDiff`` component. ``contact_maps`` is track data whose values
    are bins x bins x tracks (one square map per track) and that has an ``interval``. The scale is
    clipped to ``[vmin, vmax]``. Writes a PNG and returns its path with the data.
    """

    target = _png(output_path)
    data = _track_data(contact_maps, "contact_maps")
    low, high = _number(vmin, "vmin"), _number(vmax, "vmax")
    if not low < 0 < high:
        raise ConfigurationError("vmin must be below 0 and vmax above 0.")
    if data.interval is None or data.values.ndim != 3:
        raise ConfigurationError("contact_maps needs an interval and 3-D values (bins x bins x tracks).")
    components = require("alphagenome.visualization.plot_components", _EXTRA)
    _save(target, lambda: components.plot(
        [components.ContactMapsDiff(data, vmin=low, vmax=high)], interval=data.interval))
    return _done("plot_contact_map_difference", image_path=str(target), vmin=low, vmax=high,
                 data=_track_dict(data))


def alphagenome_plot_overlaid_tracks(
    tracks: dict[str, Any],
    output_path: str,
    *,
    interval: str | dict[str, Any] | None = None,
    shared_y_scale: bool = False,
    alpha: float = 0.8,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw several track sets on top of each other, such as REF and ALT predictions.

    This is alphagenome's ``OverlaidTracks`` component. ``tracks`` maps a label (``"REF"``,
    ``"ALT"``) to track data; the sets need the same track names and an ``interval``. ``interval``
    is the plotted window (default: the interval of the first set). Writes a PNG.
    """

    target = _png(output_path)
    if not isinstance(tracks, Mapping) or len(tracks) < 1:
        raise ConfigurationError("tracks must map labels to track data.")
    mapping = {str(label): _track_data(item, f"tracks[{label}]") for label, item in tracks.items()}
    first = next(iter(mapping.values()))
    window = _interval(interval, "interval") if interval is not None else first.interval
    if window is None or any(d.interval is None for d in mapping.values()):
        raise ConfigurationError("every track set needs an interval.")
    share, opacity = _flag(shared_y_scale, "shared_y_scale"), _number(alpha, "alpha")
    components = require("alphagenome.visualization.plot_components", _EXTRA)
    _save(target, lambda: components.plot(
        [components.OverlaidTracks(mapping, shared_y_scale=share, alpha=opacity)], interval=window))
    return _done("plot_overlaid_tracks", image_path=str(target), labels=list(mapping),
                 interval=_interval_dict(window), data={k: _track_dict(v) for k, v in mapping.items()})


def alphagenome_plot_annotations(
    interval: str | dict[str, Any],
    output_path: str,
    *,
    intervals: list[Any] | None = None,
    variants: list[Any] | None = None,
    track: dict[str, Any] | None = None,
    labels: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mark genomic intervals and/or variants across a plot window.

    ``intervals`` are drawn with alphagenome's ``IntervalAnnotation`` (shaded rectangles) and
    ``variants`` with ``VariantAnnotation`` (thin lines with labels). ``track`` (optional) is track
    data drawn underneath; without it an empty panel is used. ``labels`` replaces the default
    labels of the intervals (when given) or of the variants. Writes a PNG.
    """

    target = _png(output_path)
    window = _interval(interval, "interval")
    spans = _intervals(intervals, "intervals") if intervals else []
    snvs = [_variant(v, f"variants[{i}]") for i, v in enumerate(variants or [])]
    if not spans and not snvs:
        raise ConfigurationError("give at least one interval or variant to annotate.")
    names = None
    if labels is not None:
        if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence):
            raise ConfigurationError("labels must be a list of text.")
        names = [str(x) for x in labels]
    data = _track_data(track, "track") if track is not None else None
    components = require("alphagenome.visualization.plot_components", _EXTRA)

    def draw() -> Any:
        notes: list[Any] = []
        if spans:
            notes.append(components.IntervalAnnotation(spans, labels=names if not snvs else None))
        if snvs:
            notes.append(components.VariantAnnotation(snvs, labels=names if not spans else None))
        body = [components.Tracks(data)] if data is not None else [components.EmptyComponent(2.0)]
        return components.plot(body, interval=window, annotations=notes)

    _save(target, draw)
    return _done("plot_annotations", image_path=str(target), interval=_interval_dict(window),
                 intervals=[_interval_dict(x) for x in spans], variants=[_variant_dict(v) for v in snvs])


__all__ = [
    "alphagenome_adjust_interval",
    "alphagenome_change_track_resolution",
    "alphagenome_extract_fasta_interval",
    "alphagenome_fix_transcript_truncation",
    "alphagenome_interval_center",
    "alphagenome_interval_coverage",
    "alphagenome_interval_relation",
    "alphagenome_interval_strand",
    "alphagenome_interval_within_reference",
    "alphagenome_junction_regions",
    "alphagenome_junctions_in_interval",
    "alphagenome_normalize_junction_counts",
    "alphagenome_output_metadata",
    "alphagenome_overlap_ranges",
    "alphagenome_plot_annotations",
    "alphagenome_plot_contact_map_difference",
    "alphagenome_plot_overlaid_tracks",
    "alphagenome_resize_track_data",
    "alphagenome_reverse_complement_track_data",
    "alphagenome_scorer_metadata",
    "alphagenome_slice_track_data",
    "alphagenome_split_variant",
    "alphagenome_track_bin_index",
    "alphagenome_transcript_cds_offset",
    "alphagenome_transcripts_from_gtf",
    "alphagenome_variant_overlap",
]

_register_public_verbose(globals(), tuple(__all__))
