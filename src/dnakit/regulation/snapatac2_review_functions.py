"""Fragment export, coverage, BAM conversion, contact and base-pair imports, 10x QC, MACS3 peaks and
motif binding edges over snapATAC2.

Fragment inputs use the five-column format of ``read_fragment_qc()``. Cell groups are a mapping from
barcode to group name and must name every barcode in the fragment file. Gene annotations are local
GTF files, genomes are local FASTA files and no data is downloaded. Output files are created only
when they are absent; an existing output path raises ``ConfigurationError``.
"""

from __future__ import annotations

import functools
import math
import re
import tempfile
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _DEFAULT_MAX_PEAKS,
    _SNAPATAC2_EXTRA,
    _as_list,
    _choose,
    _chrom_sizes,
    _flag,
    _fraction,
    _imported_fragments,
    _integer,
    _label,
    _peak_file,
    _snap_provenance,
    _upstream,
)

_MAX_CELLS = 100_000
_MAX_FRAGMENTS = 5_000_000
_MAX_CONTACTS = 5_000_000
_MAX_BIN_SIZE = 10_000_000
_MAX_REGIONS = 100_000
_MAX_MOTIFS = 1_000
_MAX_NAME = 64
_GROUP = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")
_TAG = re.compile(r"[A-Za-z][A-Za-z0-9]")
_REGION = re.compile(r"([^:\s]+):(\d+)-(\d+)")
_VALUE_HEADER = "chrom\tpos\tmethyl\tunmethyl"
_SCRATCH = "dnakit-snapatac2-"


def _keep_tempdir(function: Callable[..., Any]) -> Callable[..., Any]:
    """Restore ``tempfile.tempdir`` after the call, because snapATAC2 MACS3 sets it globally."""

    @functools.wraps(function)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        saved = tempfile.tempdir
        try:
            return function(*args, **kwargs)
        finally:
            tempfile.tempdir = saved

    return wrapper


def _file_part(value: object, name: str) -> str:
    if not isinstance(value, str) or len(value) > _MAX_NAME or set(value) & set("/\\\x00"):
        raise ConfigurationError(
            f"{name} must be a string of at most {_MAX_NAME} characters without path separators."
        )
    return value


def _suffix(value: object, name: str) -> str:
    text = _file_part(value, name)
    if not text.startswith(".") or len(text) < 2:
        raise ConfigurationError(f"{name} must start with a dot, for example '.bed'.")
    return text


def _new_output(value: str | Path, name: str, *, suffix: str | None = None) -> Path:
    path = Path(value).expanduser().resolve()
    if suffix is not None and path.suffix != suffix:
        raise ConfigurationError(f"{name} must end with {suffix}; got {path.name}.")
    if path.exists():
        raise ConfigurationError(f"Output file already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _output_dir(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path}")
    return path


def _group_labels(barcodes: Sequence[str], groups: object) -> list[str]:
    if not isinstance(groups, Mapping) or not groups:
        raise ConfigurationError("groups must be a non-empty mapping of barcode to group name.")
    missing = [barcode for barcode in barcodes if barcode not in groups]
    if missing:
        raise ConfigurationError(
            f"groups must name every barcode in the fragment file; {len(missing)} missing, "
            f"first {missing[0]!r}."
        )
    return [_label(groups[barcode], f"group of {barcode!r}") for barcode in barcodes]


def _group_targets(
    directory: Path, names: Sequence[str], prefix: str, suffix: str
) -> dict[str, Path]:
    targets: dict[str, Path] = {}
    for name in names:
        if not _GROUP.fullmatch(name):
            raise ConfigurationError(
                f"group name {name!r} must start with a letter or digit and use only letters, "
                "digits, '.', '_' or '-' (at most 64 characters)."
            )
        path = directory / f"{prefix}{name}{suffix}"
        if path.exists():
            raise ConfigurationError(f"Output file already exists: {path}")
        targets[name] = path
    return targets


def _line_count(path: Path) -> int:
    if not path.is_file():
        raise BackendExecutionError(f"snapatac2 did not write {path.name}.")
    with path.open("rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def _int_field(text: str, name: str, source: str, number: int) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise ConfigurationError(f"{source} line {number} has a non-integer {name}.") from exc


def _barcode_tag(value: object) -> str:
    if not isinstance(value, str) or not _TAG.fullmatch(value):
        raise ConfigurationError("barcode_tag must be a two-character SAM tag such as CB.")
    return value


def _plain(value: Any) -> Any:
    """Convert numpy and polars scalars into JSON types and sort mapping keys."""

    if isinstance(value, Mapping):
        return {str(key): _plain(value[key]) for key in sorted(value, key=str)}
    if isinstance(value, list | tuple):
        return [_plain(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, str | bytes):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _contact_pairs(path: Path, sizes: Mapping[str, int]) -> tuple[int, int]:
    rows = 0
    total = 0
    with path.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 6 or not fields[0]:
                raise ConfigurationError(
                    f"{path.name} line {number} must have six tab-separated fields: barcode, "
                    "chrom1, pos1, chrom2, pos2, count."
                )
            for chrom, text in ((fields[1], fields[2]), (fields[3], fields[4])):
                if chrom not in sizes:
                    raise ConfigurationError(
                        f"{path.name} line {number} uses chromosome {chrom!r}, not in chrom_sizes."
                    )
                position = _int_field(text, "position", path.name, number)
                if not 0 <= position < sizes[chrom]:
                    raise ConfigurationError(
                        f"{path.name} line {number} must satisfy 0 <= pos < {sizes[chrom]} "
                        f"on {chrom}."
                    )
            count = _int_field(fields[5], "count", path.name, number)
            if count < 1:
                raise ConfigurationError(
                    f"{path.name} line {number} must have a count of at least 1."
                )
            rows += 1
            total += count
            if rows > _MAX_CONTACTS:
                raise ConfigurationError(
                    f"{path.name} has more than {_MAX_CONTACTS} contact pairs."
                )
    if rows == 0:
        raise ConfigurationError(f"{path.name} contains no contact pairs.")
    return rows, total


def _value_rows(path: Path, sizes: Mapping[str, int]) -> int:
    rows = 0
    seen: set[tuple[str, int]] = set()
    with path.open("rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n")
        if header != _VALUE_HEADER:
            raise ConfigurationError(
                f"{path.name} must start with the header chrom, pos, methyl, unmethyl."
            )
        for number, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 4:
                raise ConfigurationError(
                    f"{path.name} line {number} must have four tab-separated fields."
                )
            chrom = fields[0]
            if chrom not in sizes:
                raise ConfigurationError(
                    f"{path.name} line {number} uses chromosome {chrom!r}, not in chrom_sizes."
                )
            position = _int_field(fields[1], "position", path.name, number)
            if not 0 <= position < sizes[chrom]:
                raise ConfigurationError(
                    f"{path.name} line {number} must satisfy 0 <= pos < {sizes[chrom]} on {chrom}."
                )
            if (chrom, position) in seen:
                raise ConfigurationError(f"{path.name} line {number} repeats {chrom}:{position}.")
            seen.add((chrom, position))
            methyl = _int_field(fields[2], "methyl count", path.name, number)
            unmethyl = _int_field(fields[3], "unmethyl count", path.name, number)
            if methyl < 0 or unmethyl < 0 or methyl + unmethyl < 1:
                raise ConfigurationError(
                    f"{path.name} line {number} needs non-negative counts, total at least 1."
                )
            rows += 1
    if rows == 0:
        raise ConfigurationError(f"{path.name} contains no values.")
    return rows


def _regions(values: object) -> list[str]:
    if isinstance(values, str | bytes) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("regions must be a non-empty list of chrom:start-end strings.")
    if len(values) > _MAX_REGIONS:
        raise ConfigurationError(f"regions has more than {_MAX_REGIONS} entries.")
    checked: list[str] = []
    for value in values:
        match = _REGION.fullmatch(value) if isinstance(value, str) else None
        if match is None or int(match.group(2)) >= int(match.group(3)):
            raise ConfigurationError(
                f"region {value!r} must look like chrom:start-end with start < end."
            )
        checked.append(value)
    return checked


@_keep_tempdir
def export_fragment_files(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    groups: Mapping[str, str],
    out_dir: str | Path,
    *,
    prefix: str = "",
    suffix: str = ".bed",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one five-column fragment file per cell group, named ``{prefix}{group}{suffix}``.

    The fragment file is imported with snapATAC2 and every barcode must appear in ``groups``.
    Each exported line is one stored fragment with read count 1.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    directory = _output_dir(out_dir)
    text_prefix = _file_part(prefix, "prefix")
    text_suffix = _suffix(suffix, "suffix")
    with _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=1,
        max_fragments=_MAX_FRAGMENTS,
        max_cells=_MAX_CELLS,
    ) as (adata, _work):
        labels = _group_labels(list(adata.obs_names), groups)
        names = sorted(set(labels))
        targets = _group_targets(directory, names, text_prefix, text_suffix)
        _upstream(
            lambda: snap.ex.export_fragments(
                adata,
                groupby=labels,
                out_dir=str(directory),
                prefix=text_prefix,
                suffix=text_suffix,
                compression=None,
            )
        )
    records = [
        {
            "group": name,
            "file": targets[name].name,
            "fragment_count": _line_count(targets[name]),
        }
        for name in names
    ]
    return {
        "status": "completed",
        "operation": "export_fragment_files",
        "group_count": len(records),
        "groups": records,
        "provenance": _snap_provenance(group_source="mapping of barcode to group name"),
    }


@_keep_tempdir
def export_coverage_tracks(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    groups: Mapping[str, str],
    out_dir: str | Path,
    *,
    bin_size: int = 10,
    normalization: str | None = None,
    output_format: str = "bedgraph",
    max_fragment_length: int | None = None,
    prefix: str = "",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one genome-wide coverage track per cell group as bedGraph or bigWig.

    Tracks are named ``{prefix}{group}.bedgraph`` or ``{prefix}{group}.bw``. Without normalization
    the values are fragment counts per bin.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    directory = _output_dir(out_dir)
    width = _integer(bin_size, "bin_size", 1, _MAX_BIN_SIZE)
    if normalization is not None:
        normalization = _choose(normalization, ("RPKM", "CPM", "BPM"), "normalization")
    output = _choose(output_format, ("bedgraph", "bigwig"), "output_format")
    suffix = ".bedgraph" if output == "bedgraph" else ".bw"
    max_length = (
        None
        if max_fragment_length is None
        else _integer(max_fragment_length, "max_fragment_length", 1, 1_000_000)
    )
    text_prefix = _file_part(prefix, "prefix")
    with _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=1,
        max_fragments=_MAX_FRAGMENTS,
        max_cells=_MAX_CELLS,
    ) as (adata, work):
        labels = _group_labels(list(adata.obs_names), groups)
        names = sorted(set(labels))
        targets = _group_targets(directory, names, text_prefix, suffix)
        _upstream(
            lambda: snap.ex.export_coverage(
                adata,
                groupby=labels,
                bin_size=width,
                normalization=normalization,
                max_frag_length=max_length,
                counting_strategy="fragment",
                out_dir=str(directory),
                prefix=text_prefix,
                suffix=suffix,
                output_format=output,
                compression=None,
                n_jobs=1,
                tempdir=str(work),
            )
        )
    records = []
    for name in names:
        if not targets[name].is_file():
            raise BackendExecutionError(f"snapatac2 did not write {targets[name].name}.")
        records.append({"group": name, "file": targets[name].name})
    return {
        "status": "completed",
        "operation": "export_coverage_tracks",
        "output_format": output,
        "bin_size": width,
        "normalization": normalization,
        "track_count": len(records),
        "tracks": records,
        "provenance": _snap_provenance(group_source="mapping of barcode to group name"),
    }


@_keep_tempdir
def fragments_from_bam(
    bam_file: str | Path,
    output_fragment_file: str | Path,
    *,
    barcode_tag: str = "CB",
    is_paired: bool = True,
    min_mapq: int = 30,
    mitochondrial_chroms: Sequence[str] = ("chrM", "M"),
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a barcode-tagged BAM file into a sorted five-column fragment file.

    Fragments use the snapATAC2 Tn5 shift (+4 at the start, -5 at the end). The returned QC values
    are the snapATAC2 sequencing, mapping and duplicate fractions for the same BAM file.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    bam = input_path(bam_file)
    output = _new_output(output_fragment_file, "output_fragment_file")
    tag = _barcode_tag(barcode_tag)
    paired = _flag(is_paired, "is_paired")
    mapq = _integer(min_mapq, "min_mapq", 0, 255)
    if isinstance(mitochondrial_chroms, str | bytes) or not isinstance(
        mitochondrial_chroms, Sequence
    ):
        raise ConfigurationError("mitochondrial_chroms must be a list of chromosome names.")
    chroms = [_label(name, "mitochondrial chromosome") for name in mitochondrial_chroms]
    with tempfile.TemporaryDirectory(prefix=_SCRATCH) as scratch:
        qc = _upstream(
            lambda: snap.pp.make_fragment_file(
                str(bam),
                str(output),
                is_paired=paired,
                barcode_tag=tag,
                min_mapq=mapq,
                chrM=chroms,
                source=None,
                compression=None,
                tempdir=scratch,
            )
        )
    fragments = _line_count(output)
    with output.open("rt", encoding="utf-8") as handle:
        barcodes = {line.split("\t")[3] for line in handle if line.count("\t") >= 4}
    return {
        "status": "completed",
        "operation": "fragments_from_bam",
        "fragment_count": fragments,
        "barcode_count": len(barcodes),
        "qc": _plain(qc),
        "output_file": str(output),
        "provenance": _snap_provenance(
            input_format="BAM alignment file with a cell barcode tag",
            barcode_tag=tag,
            min_mapq=mapq,
        ),
    }


@_keep_tempdir
def import_contact_pairs(
    contact_file: str | Path,
    chrom_sizes: Mapping[str, int],
    output_h5ad: str | Path,
    *,
    bin_size: int = 500_000,
    sorted_by_barcode: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin single-cell contact pairs into a backed snapATAC2 h5ad file.

    Each input line is ``barcode, chrom1, pos1, chrom2, pos2, count``. Positions are 0-based and
    must lie inside ``chrom_sizes``. Rows of the matrix are cells; columns are bin pairs.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    sizes = _chrom_sizes(chrom_sizes)
    path = input_path(contact_file)
    output = _new_output(output_h5ad, "output_h5ad", suffix=".h5ad")
    width = _integer(bin_size, "bin_size", 1, _MAX_BIN_SIZE)
    sorted_flag = _flag(sorted_by_barcode, "sorted_by_barcode")
    pair_count, contact_count = _contact_pairs(path, sizes)
    with tempfile.TemporaryDirectory(prefix=_SCRATCH) as scratch:
        adata = _upstream(
            lambda: snap.pp.import_contacts(
                str(path),
                sizes,
                sorted_by_barcode=sorted_flag,
                bin_size=width,
                file=str(output),
                tempdir=scratch,
            )
        )
        try:
            cell_count = int(adata.n_obs)
            feature_count = int(adata.n_vars)
        finally:
            adata.close()
    return {
        "status": "completed",
        "operation": "import_contact_pairs",
        "cell_count": cell_count,
        "feature_count": feature_count,
        "bin_size": width,
        "contact_pair_count": pair_count,
        "contact_count": contact_count,
        "chromosomes": [{"chrom": name, "length": sizes[name]} for name in sorted(sizes)],
        "output_file": str(output),
        "provenance": _snap_provenance(
            input_format="contact pairs: barcode, chrom1, pos1, chrom2, pos2, count",
            coordinates="0-based positions binned into fixed-width bins",
        ),
    }


@_keep_tempdir
def import_base_pair_values(
    input_dir: str | Path,
    chrom_sizes: Mapping[str, int],
    output_h5ad: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import one per-cell base-pair count table into a snapATAC2 h5ad file.

    ``input_dir`` must hold only ``.tsv`` files, one per cell. The cell name is the file stem and
    each file has the header ``chrom, pos, methyl, unmethyl`` followed by one row per base.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    sizes = _chrom_sizes(chrom_sizes)
    directory = Path(input_dir).expanduser().resolve()
    if not directory.is_dir():
        raise ConfigurationError(f"Input directory does not exist: {directory}")
    entries = sorted(directory.iterdir())
    if any(not (entry.is_file() and entry.suffix == ".tsv") for entry in entries):
        raise ConfigurationError("input_dir must hold only .tsv files, one per cell.")
    if not entries:
        raise ConfigurationError("input_dir contains no cell files.")
    if len(entries) > _MAX_CELLS:
        raise ConfigurationError(f"input_dir has more than {_MAX_CELLS} cell files.")
    output = _new_output(output_h5ad, "output_h5ad", suffix=".h5ad")
    expected = {entry.stem: _value_rows(entry, sizes) for entry in entries}
    adata = _upstream(lambda: snap.pp.import_values(str(directory), sizes, file=str(output)))
    try:
        cells = [
            {"cell": str(name), "value_count": int(count)}
            for name, count in zip(adata.obs_names, _as_list(adata.obs["num_values"]), strict=True)
        ]
    finally:
        adata.close()
    cells.sort(key=lambda item: item["cell"])
    if {item["cell"]: item["value_count"] for item in cells} != expected:
        raise BackendExecutionError("snapatac2 stored different value counts than the input files.")
    return {
        "status": "completed",
        "operation": "import_base_pair_values",
        "cell_count": len(cells),
        "reference_length": sum(sizes.values()),
        "total_value_count": sum(item["value_count"] for item in cells),
        "cells": cells,
        "output_file": str(output),
        "provenance": _snap_provenance(
            input_format="one tsv per cell: chrom, pos, methyl, unmethyl",
            coordinates="0-based base positions on each chromosome",
        ),
    }


@_keep_tempdir
def tenx_atac_metrics(
    bam_file: str | Path,
    output_fragment_file: str | Path,
    output_h5ad: str | Path,
    chrom_sizes: Mapping[str, int],
    gene_gtf: str | Path,
    peaks_bed: str | Path,
    *,
    barcode_tag: str = "CB",
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the snapATAC2 10x ATAC QC recipe on a barcode-tagged BAM file with given peaks.

    The recipe converts the BAM to fragments, imports them, computes TSS enrichment from the local
    GTF, and summarizes read, cell, complexity, mapping and peak metrics. Peaks must be a BED file;
    MACS3 peak calling is not run here.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    bam = input_path(bam_file)
    gtf = input_path(gene_gtf)
    sizes = _chrom_sizes(chrom_sizes)
    peaks = _peak_file(peaks_bed, sizes, _DEFAULT_MAX_PEAKS)
    fragments = _new_output(output_fragment_file, "output_fragment_file")
    h5ad = _new_output(output_h5ad, "output_h5ad", suffix=".h5ad")
    tag = _barcode_tag(barcode_tag)
    with tempfile.TemporaryDirectory(prefix=_SCRATCH) as scratch:
        metrics = _upstream(
            lambda: snap.pp.recipe_10x_metrics(
                str(bam),
                str(fragments),
                str(h5ad),
                peaks=str(peaks),
                barcode_tag=tag,
                source=None,
                chrom_sizes=sizes,
                gene_anno=str(gtf),
                sorted_by_barcode=False,
                n_jobs=1,
                tempdir=scratch,
            )
        )
    if not fragments.is_file() or not h5ad.is_file():
        raise BackendExecutionError("snapatac2 did not write the fragment and h5ad files.")
    return {
        "status": "completed",
        "operation": "tenx_atac_metrics",
        "metrics": _plain(metrics),
        "fragment_file": str(fragments),
        "h5ad_file": str(h5ad),
        "provenance": _snap_provenance(
            input_format="BAM alignment file with a cell barcode tag, GTF genes and BED peaks",
            barcode_tag=tag,
        ),
    }


@_keep_tempdir
def call_peaks_macs3(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    groups: Mapping[str, str],
    *,
    qvalue: float = 0.05,
    shift: int = -100,
    extsize: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call open chromatin peaks with MACS3 separately for each cell group.

    Fragments are imported into a temporary snapATAC2 file; every barcode must appear in ``groups``.
    Peaks are sorted by chromosome and position within each group.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    threshold = _fraction(qvalue, "qvalue")
    if threshold == 0.0:
        raise ConfigurationError("qvalue must be greater than 0 and at most 1.")
    shift_value = _integer(shift, "shift", -10_000, 10_000)
    extension = _integer(extsize, "extsize", 1, 10_000)
    with _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=1,
        max_fragments=_MAX_FRAGMENTS,
        max_cells=_MAX_CELLS,
    ) as (adata, work):
        labels = _group_labels(list(adata.obs_names), groups)
        frames = _upstream(
            lambda: snap.tl.macs3(
                adata,
                groupby=labels,
                qvalue=threshold,
                shift=shift_value,
                extsize=extension,
                inplace=False,
                n_jobs=1,
                tempdir=str(work),
            )
        )
    records = []
    for name in sorted(frames):
        rows = frames[name].sort(["chrom", "start", "end"]).to_dicts()
        peaks = [_plain(row) for row in rows]
        records.append({"group": name, "peak_count": len(peaks), "peaks": peaks})
    return {
        "status": "completed",
        "operation": "call_peaks_macs3",
        "group_count": len(records),
        "qvalue": threshold,
        "shift": shift_value,
        "extsize": extension,
        "groups": records,
        "provenance": _snap_provenance(
            peak_caller="MACS3 on snapATAC2 fragment pileups, one call per group",
            group_source="mapping of barcode to group name",
        ),
    }


@_keep_tempdir
def add_tf_binding_edges(
    regions: Sequence[str],
    gtf_file: str | Path,
    motif_file: str | Path,
    genome_fasta: str | Path,
    *,
    pvalue: float = 1e-5,
    upstream: int = 250_000,
    downstream: int = 250_000,
    coding_gene_only: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a region-to-gene network from a GTF and add motif-to-region edges by genome scanning.

    ``regions`` are ``chrom:start-end`` strings. Motifs come from a MEME file and are scanned in the
    region sequences of the local FASTA genome. The result lists the motif-to-region edges only.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    region_list = _regions(regions)
    gtf = input_path(gtf_file)
    motifs_path = input_path(motif_file)
    fasta = input_path(genome_fasta)
    threshold = _fraction(pvalue, "pvalue")
    if threshold == 0.0:
        raise ConfigurationError("pvalue must be greater than 0 and at most 1.")
    up = _integer(upstream, "upstream", 0, 1_000_000)
    down = _integer(downstream, "downstream", 0, 1_000_000)
    coding = _flag(coding_gene_only, "coding_gene_only")
    motifs = _upstream(lambda: snap.read_motifs(str(motifs_path)))
    if not motifs:
        raise ConfigurationError(f"{motifs_path.name} contains no motifs.")
    if len(motifs) > _MAX_MOTIFS:
        raise ConfigurationError(f"{motifs_path.name} has more than {_MAX_MOTIFS} motifs.")
    network = _upstream(
        lambda: snap.tl.init_network_from_annotation(
            region_list,
            str(gtf),
            upstream=up,
            downstream=down,
            id_type="gene_name",
            coding_gene_only=coding,
        )
    )
    _upstream(
        lambda: snap.tl.add_tf_binding(
            network, motifs=motifs, genome_fasta=str(fasta), pvalue=threshold
        )
    )
    node_types = Counter(node.type for node in network.nodes())
    edges = []
    for src, dst, _data in network.edge_index_map().values():
        if network[src].type == "motif" and network[dst].type == "region":
            edges.append({"motif": str(network[src].id), "region": str(network[dst].id)})
    edges.sort(key=lambda item: (item["motif"], item["region"]))
    return {
        "status": "completed",
        "operation": "add_tf_binding_edges",
        "motif_count": len(motifs),
        "region_count": node_types.get("region", 0),
        "gene_count": node_types.get("gene", 0),
        "binding_edge_count": len(edges),
        "binding_edges": edges,
        "provenance": _snap_provenance(
            input_format="MEME motif file, GTF annotation and local FASTA genome",
            coordinates="chrom:start-end region strings as read by snapATAC2",
            pvalue=threshold,
        ),
    }


__all__ = [
    "add_tf_binding_edges",
    "call_peaks_macs3",
    "export_coverage_tracks",
    "export_fragment_files",
    "fragments_from_bam",
    "import_base_pair_values",
    "import_contact_pairs",
    "tenx_atac_metrics",
]

_register_public_verbose(globals(), tuple(__all__))
