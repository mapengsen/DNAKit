"""deepTools wrappers for GC-bias tables, the synthetic fingerprint JSD and bigWig tile bedGraph.

These functions complete the deepTools coverage of ``deeptools_functions``. They call deepTools'
Python helpers on local files and run every call in one process, so results do not depend on
process scheduling.

Two upstream details shape the wrappers. ``computeGCBias.tabulateGCcontent`` and
``countReadsPerGC`` read a module-level dict ``global_vars`` that only the command line sets; the
wrappers build that dict as the command line does, install it for the call and restore the
previous state afterwards. ``writeBedGraph_bam_and_bw.writeBedGraph`` calls
``getCoverageFromBam`` and ``getSmoothRange``, which deepTools 4.0.0 no longer defines, so its BAM
branch raises ``NameError``; only the bigWig branch works and only that branch is wrapped.

deepTools calls ``sys.exit`` on some bad inputs; those exits are converted into
``ConfigurationError``. Packages are imported only when a function runs.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.deeptools_functions import (
    _flag,
    _integer,
    _new_path,
    _number,
    _path,
    _plain,
    _provenance,
    _sequence,
    _upstream_call,
)
from dnakit.exceptions import ConfigurationError

_DEEPTOOLS_EXTRA = "bio-deeptools"
_TWOBIT_EXTRA = "bio-py2bit"
_BIGWIG_EXTRA = "bio-pybigwig"
_MAX_FILES = 32
_MAX_BINS = 10_000_000
_MAX_FRAGMENT = 1_000_000
_MAX_GENOME = 100_000_000_000
_MAX_SAMPLE_SIZE = 2_000_000_000
_MAX_POSITIONS = 200_000_000
_MAX_WINDOWS = 1_000_000
_MAX_PREVIEW = 20
_DEFAULT_SAMPLE_SIZE = 50_000_000
_REGION = re.compile(r"([^:\s]+):([0-9]+):([0-9]+)")
_REDUCTIONS = ("mean", "sum")


def _import_deeptools(module: str) -> Any:
    """Import a deepTools module and undo the numpy error settings it changes on import."""

    numpy = require("numpy", _DEEPTOOLS_EXTRA)
    saved = numpy.geterr()
    try:
        return require(module, _DEEPTOOLS_EXTRA)
    finally:
        numpy.seterr(**saved)


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a deepTools function with the numpy error handling that deepTools sets for itself."""

    numpy = require("numpy", _DEEPTOOLS_EXTRA)
    with numpy.errstate(all="ignore"):
        return _upstream_call(function, *args, **kwargs)


def _region_text(region: object) -> str | None:
    if region is None:
        return None
    if not isinstance(region, str):
        raise ConfigurationError("region must be a string such as chr1:1000:5000.")
    match = _REGION.fullmatch(region)
    if match is None or int(match.group(2)) >= int(match.group(3)):
        raise ConfigurationError(
            "region must have the form chrom:start:end with start below end, for example "
            "chr1:1000:5000."
        )
    return region


def _check_region(region: str | None, sizes: Sequence[tuple[str, int]]) -> None:
    if region is None:
        return
    mapreduce = _import_deeptools("deeptools.mapReduce")
    try:
        mapreduce.getUserRegion(dict(sizes), region)
    except (NameError, ValueError) as exc:
        raise ConfigurationError(str(exc).splitlines()[0]) from exc


class _GcBiasState:
    """Install ``global_vars`` on ``deeptools.computeGCBias`` and restore the old state."""

    def __init__(self, module: Any, values: dict[str, Any]) -> None:
        self._module = module
        self._values = values
        self._missing = object()
        self._previous: Any = self._missing

    def __enter__(self) -> None:
        self._previous = getattr(self._module, "global_vars", self._missing)
        self._module.global_vars = self._values

    def __exit__(self, *exc_info: object) -> None:
        if self._previous is self._missing:
            if hasattr(self._module, "global_vars"):
                del self._module.global_vars
        else:
            self._module.global_vars = self._previous


def _prepare_gc_inputs(
    twobit_path: object,
    bam_path: object,
    effective_genome_size: object,
    sample_size: object,
    step_size: object,
    step_factor: int,
    region: object,
) -> dict[str, Any]:
    """Check the inputs of the GC-bias functions and derive what the command line derives."""

    twobit = _path(twobit_path, "twobit_path")
    bam = _path(bam_path, "bam_path")
    effective = _integer(effective_genome_size, "effective_genome_size", 1, _MAX_GENOME)
    sample = _integer(sample_size, "sample_size", 2, _MAX_SAMPLE_SIZE)
    if step_size is not None:
        _integer(step_size, "step_size", 1, _MAX_GENOME)
    region_text = _region_text(region)
    py2bit = require("py2bit", _TWOBIT_EXTRA)
    bam_handler = _import_deeptools("deeptools.bamHandler")
    utilities = _import_deeptools("deeptools.utilities")
    try:
        handle = py2bit.open(str(twobit))
    except RuntimeError as exc:
        raise ConfigurationError(f"{twobit.name} is not a readable 2bit file.") from exc
    try:
        lengths = {str(name): int(length) for name, length in handle.chroms().items()}
    finally:
        handle.close()
    reader, mapped, _unmapped, _stats = _call(
        bam_handler.openBam, str(bam), returnStats=True, nThreads=1
    )
    try:
        references = [str(name) for name in reader.references]
        bam_lengths = [int(length) for length in reader.lengths]
    finally:
        reader.close()
    if int(mapped) < 1:
        raise ConfigurationError(f"{bam.name} has no mapped reads in its index.")
    try:
        mapping = _call(utilities.tbitToBamChrName, list(lengths), references)
    except (ConfigurationError, IndexError, KeyError, NameError) as exc:
        raise ConfigurationError(
            "The chromosome names of the 2bit and BAM files cannot be matched "
            f"(2bit: {sorted(lengths)[:5]}, BAM: {references[:5]})."
        ) from exc
    bam_to_bit = {str(bam_name): str(bit_name) for bit_name, bam_name in mapping.items()}
    sizes = [
        (name, length)
        for name, length in zip(references, bam_lengths, strict=True)
        if bam_to_bit.get(name) in lengths
    ]
    if not sizes:
        raise ConfigurationError("The 2bit and BAM files share no chromosome.")
    genome_size = sum(lengths.values())
    step = step_size if step_size is not None else max(int(genome_size / sample), 1) * step_factor
    positions = sum(-(-length // step) for _, length in sizes)
    _check_region(region_text, sizes)
    return {
        "twobit": twobit,
        "bam": bam,
        "effective": effective,
        "sample": sample,
        "mapped": int(mapped),
        "reads_per_bp": float(mapped) / effective,
        "genome_size": genome_size,
        "step": step,
        "positions": positions,
        "sizes": sizes,
        "mapping": {str(bit): str(bam_name) for bit, bam_name in mapping.items()},
        "region": region_text,
    }


def tabulate_gc_fragment_counts(
    twobit_path: str | Path,
    bam_path: str | Path,
    fragment_length: int,
    effective_genome_size: int,
    *,
    sample_size: int = _DEFAULT_SAMPLE_SIZE,
    step_size: int | None = None,
    region: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tabulate GC content of genome fragments against the forward read starts, per GC count."""

    length = _integer(fragment_length, "fragment_length", 1, _MAX_FRAGMENT)
    context = _prepare_gc_inputs(
        twobit_path, bam_path, effective_genome_size, sample_size, step_size, 1, region
    )
    if context["positions"] > _MAX_POSITIONS:
        raise ConfigurationError(
            f"step_size {context['step']} would sample {context['positions']} positions; "
            f"the limit is {_MAX_POSITIONS}."
        )
    stats = require("scipy.stats", _DEEPTOOLS_EXTRA)
    module = _import_deeptools("deeptools.computeGCBias")
    sample = context["sample"]
    rate = context["reads_per_bp"]
    max_reads = float(stats.poisson(4 * rate * length).isf(1 / sample))
    min_reads = float(stats.poisson(0.25 * rate * length).ppf(1 / sample))
    values = {
        "2bit": str(context["twobit"]),
        "bam": str(context["bam"]),
        "filter_out": None,
        "extra_sampling_file": None,
        "genome_size": context["genome_size"],
        "total_reads": context["mapped"],
        "reads_per_bp": rate,
        "max_reads": max_reads,
        "min_reads": min_reads,
    }
    with _GcBiasState(module, values):
        table = _call(
            module.tabulateGCcontent,
            {"median": length},
            context["mapping"],
            context["step"],
            context["sizes"],
            numberOfProcessors=1,
            verbose=False,
            region=context["region"],
        )
    rows = _plain(table)
    return {
        "status": "completed",
        "operation": "tabulate_gc_fragment_counts",
        "twobit": str(context["twobit"]),
        "bam": str(context["bam"]),
        "fragment_length": length,
        "effective_genome_size": context["effective"],
        "sample_size": sample,
        "step_size": context["step"],
        "region": context["region"],
        "mapped_reads": context["mapped"],
        "reads_per_bp": rate,
        "max_reads": _plain(max_reads),
        "min_reads": _plain(min_reads),
        "positions_sampled": sum(int(row[1]) for row in rows),
        "read_starts_counted": sum(int(row[0]) for row in rows),
        "gc_table": [
            {
                "gc_bases": index,
                "read_starts": int(row[0]),
                "genome_positions": int(row[1]),
                "read_ratio": row[2],
            }
            for index, row in enumerate(rows)
        ],
        "provenance": _provenance(
            "deeptools",
            "fragments of fragment_length bases are sampled every step_size positions; the ratio "
            "is the read starts per sampled position of that GC count, relative to the same "
            "density over all sampled positions",
            note="sampled positions with at least max_reads read starts are skipped as peaks; "
            "fragments with more than 5% N bases are not counted; forward reads only; "
            "a BAM that lacks the chr prefix of the 2bit names is matched, unlike in computeGCBias",
            coordinate_system="0-based half-open",
        ),
    }


def count_reads_per_gc_bins(
    twobit_path: str | Path,
    bam_path: str | Path,
    effective_genome_size: int,
    *,
    region_size: int = 300,
    sample_size: int = _DEFAULT_SAMPLE_SIZE,
    step_size: int | None = None,
    region: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """GC fraction and overlapping read count of genome windows of ``region_size`` bases."""

    size = _integer(region_size, "region_size", 1, _MAX_FRAGMENT)
    context = _prepare_gc_inputs(
        twobit_path, bam_path, effective_genome_size, sample_size, step_size, 10, region
    )
    if context["positions"] > _MAX_WINDOWS:
        raise ConfigurationError(
            f"step_size {context['step']} would give {context['positions']} windows; "
            f"the limit is {_MAX_WINDOWS}."
        )
    module = _import_deeptools("deeptools.computeGCBias")
    values = {
        "2bit": str(context["twobit"]),
        "bam": str(context["bam"]),
        "filter_out": None,
        "extra_sampling_file": None,
        "genome_size": context["genome_size"],
        "total_reads": context["mapped"],
        "reads_per_bp": context["reads_per_bp"],
    }
    with _GcBiasState(module, values):
        table = _call(
            module.countReadsPerGC,
            size,
            context["mapping"],
            context["step"],
            context["sizes"],
            numberOfProcessors=1,
            verbose=False,
            region=context["region"],
        )
    rows = _plain(table)
    windows = [{"reads": int(row[0]), "gc_fraction": row[1]} for row in rows]
    return {
        "status": "completed",
        "operation": "count_reads_per_gc_bins",
        "twobit": str(context["twobit"]),
        "bam": str(context["bam"]),
        "region_size": size,
        "effective_genome_size": context["effective"],
        "sample_size": context["sample"],
        "step_size": context["step"],
        "region": context["region"],
        "mapped_reads": context["mapped"],
        "window_count": len(windows),
        "mean_reads": sum(window["reads"] for window in windows) / len(windows)
        if windows
        else None,
        "windows": windows,
        "provenance": _provenance(
            "deeptools",
            "windows of region_size bases start every step_size positions; reads of both strands "
            "that overlap the window are counted",
            note="windows in genome order; windows with more than 5% N bases are left out; "
            "the default step is ten times the one of tabulate_gc_fragment_counts",
            coordinate_system="0-based half-open",
        ),
    }


def fingerprint_synthetic_jsd(coverage: Sequence[float], *, verbose: bool = True) -> dict[str, Any]:
    """Jensen-Shannon distance between per-bin coverage and a Poisson null with the same mean."""

    items = _sequence(coverage, "coverage", _MAX_BINS)
    values = [_number(item, "coverage entry", 0.0, 1e12) for item in items]
    if max(values) < 1:
        raise ConfigurationError("coverage needs at least one bin with a value of 1 or more.")
    plot = _import_deeptools("deeptools.plotFingerprint")
    ceiling = int(plot.MAXLEN)
    result = _call(plot.getSyntheticJSD, values)
    jsd = _plain(result)
    if not isinstance(jsd, float) or not math.isfinite(jsd):
        raise ConfigurationError("The synthetic JSD could not be computed for this coverage.")
    return {
        "status": "completed",
        "operation": "fingerprint_synthetic_jsd",
        "bins": len(values),
        "mean_coverage": sum(values) / len(values),
        "bins_above_maximum": sum(1 for value in values if value >= ceiling),
        "jsd": jsd,
        "provenance": _provenance(
            "deeptools",
            "Jensen-Shannon distance between the cumulative coverage curve of the bins and that "
            "of a Poisson sample with the mean coverage of the bins",
            note="bins with no coverage are dropped; coverage fractions are cut to integers",
        ),
    }


def _scaled_tile_mean(tile_coverage: Sequence[float], args: dict[str, Any]) -> float:
    """Tile value for writeBedGraph: mean over the files of coverage times its scale factor."""

    scales = args["scale_factors"]
    total = sum(float(value) * scale for value, scale in zip(tile_coverage, scales, strict=True))
    return total / len(tile_coverage)


def _scaled_tile_sum(tile_coverage: Sequence[float], args: dict[str, Any]) -> float:
    """Tile value for writeBedGraph: sum over the files of coverage times its scale factor."""

    scales = args["scale_factors"]
    return sum(float(value) * scale for value, scale in zip(tile_coverage, scales, strict=True))


_TILE_FUNCTIONS = {"mean": _scaled_tile_mean, "sum": _scaled_tile_sum}


def _bigwig_chromosomes(sources: Sequence[Path]) -> list[tuple[str, int]]:
    """Chromosomes shared by all bigWig files, in the order of the first file."""

    pybigwig = require("pyBigWig", _BIGWIG_EXTRA)
    tables: list[dict[str, int]] = []
    for source in sources:
        try:
            handle = pybigwig.open(str(source))
        except RuntimeError as exc:
            raise ConfigurationError(f"{source.name} is not a bigWig file.") from exc
        try:
            if not handle.isBigWig():
                raise ConfigurationError(f"{source.name} is not a bigWig file.")
            tables.append({str(name): int(length) for name, length in handle.chroms().items()})
        finally:
            handle.close()
    shared = [name for name in tables[0] if all(name in table for table in tables)]
    if not shared:
        raise ConfigurationError("The bigWig files share no chromosome.")
    for name in shared:
        if len({table[name] for table in tables}) != 1:
            raise ConfigurationError(f"{name} has different lengths in the bigWig files.")
    return [(name, tables[0][name]) for name in shared]


def _read_bedgraph(path: Path) -> tuple[int, list[dict[str, Any]], list[str]]:
    count = 0
    preview: list[dict[str, Any]] = []
    chromosomes: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            count += 1
            if fields[0] not in chromosomes:
                chromosomes.append(fields[0])
            if len(preview) < _MAX_PREVIEW:
                preview.append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[1]),
                        "end": int(fields[2]),
                        "value": _plain(float(fields[3])),
                    }
                )
    return count, preview, chromosomes


def bigwig_bedgraph_tiles(
    bigwig_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    tile_size: int = 50,
    reduction: str = "mean",
    scale_factors: Sequence[float] | None = None,
    region: str | None = None,
    missing_data_as_zero: bool = True,
    fixed_step: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a bedGraph of per-tile scaled mean or sum of coverage from one or more bigWig files."""

    files = _sequence(bigwig_paths, "bigwig_paths", _MAX_FILES)
    sources = [_path(item, "bigwig_paths") for item in files]
    target = _new_path(output_path, "output_path")
    tile = _integer(tile_size, "tile_size", 1, 10_000_000)
    if reduction not in _REDUCTIONS:
        raise ConfigurationError("reduction must be 'mean' or 'sum'.")
    if scale_factors is None:
        factors = [1.0] * len(sources)
    else:
        given = _sequence(scale_factors, "scale_factors", _MAX_FILES)
        if len(given) != len(sources):
            raise ConfigurationError("scale_factors needs one value per bigWig file.")
        factors = [_number(item, "scale_factors entry", -1e12, 1e12) for item in given]
    region_text = _region_text(region)
    as_zero = _flag(missing_data_as_zero, "missing_data_as_zero")
    fixed = _flag(fixed_step, "fixed_step")
    sizes = _bigwig_chromosomes(sources)
    _check_region(region_text, sizes)
    writer = _import_deeptools("deeptools.writeBedGraph_bam_and_bw")
    _call(
        writer.writeBedGraph,
        [(str(source), "bigwig") for source in sources],
        str(target),
        0,
        _TILE_FUNCTIONS[reduction],
        {"scale_factors": factors},
        tileSize=tile,
        region=region_text,
        numberOfProcessors=1,
        format="bedgraph",
        extendPairedEnds=False,
        missingDataAsZero=as_zero,
        smoothLength=0,
        fixedStep=fixed,
        verbose=False,
    )
    count, preview, chromosomes = _read_bedgraph(target)
    return {
        "status": "completed",
        "operation": "bigwig_bedgraph_tiles",
        "bigwigs": [str(source) for source in sources],
        "output": str(target),
        "tile_size": tile,
        "reduction": reduction,
        "scale_factors": factors,
        "region": region_text,
        "missing_data_as_zero": as_zero,
        "fixed_step": fixed,
        "intervals": count,
        "chromosomes": chromosomes,
        "first_intervals": preview,
        "provenance": _provenance(
            "deeptools",
            "each file is averaged over tiles of tile_size bases; the per-tile values are scaled "
            "and reduced over the files, and neighbouring tiles with equal values are merged",
            note="a region start is moved down to a multiple of tile_size; a final tile with "
            "value 0 is not written; BAM input does not work in deepTools 4.0.0",
            coordinate_system="0-based half-open",
        ),
    }


__all__ = [
    "bigwig_bedgraph_tiles",
    "count_reads_per_gc_bins",
    "fingerprint_synthetic_jsd",
    "tabulate_gc_fragment_counts",
]

_register_public_verbose(globals(), tuple(__all__))
