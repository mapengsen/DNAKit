"""Summit, coverage, motif-scan, motif-activity and score-statistic functions over gimmemotifs.

Each function runs one gimmemotifs routine on local inputs: FASTA genomes, peak or summit BED files,
bigWig tracks, PFM files, region lists and score tables. No genome, background or motif database is
downloaded, and no external motif tool is started.
"""

from __future__ import annotations

import math
import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _choose,
    _flag,
    _fraction,
    _gimme_provenance,
    _integer,
    _pfm,
)

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_MAX_SCORES = 1_000_000
_MAX_REGIONS = 100_000
_MAX_SUMMIT_FILES = 50
_MAX_SUMMITS = 1_000_000
_MAX_COVERAGE_TRACKS = 50
_MAX_REGION_LENGTH = 10_000_000
_REGION = re.compile(r"^([A-Za-z0-9_.\-]+):(\d+)-(\d+)$")
_OUTPUT_KEYS = ("pvalue", "significance")
_METRICS = (
    "roc_auc",
    "pr_auc",
    "max_fmeasure",
    "mncp",
    "max_enrichment",
    "recall_at_fdr",
    "matches_at_fpr",
    "phyper_at_fpr",
    "fraction_fpr",
    "score_at_fpr",
    "enr_at_fpr",
    "roc_auc_xlim",
    "roc_values",
)


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except SystemExit as exc:  # gimmemotifs calls sys.exit on some input errors
        raise BackendExecutionError(f"gimmemotifs exited with code {exc.code}.") from exc
    except Exception as exc:  # gimmemotifs raises several non-DNAKit exception types
        raise BackendExecutionError(f"gimmemotifs failed: {exc}") from exc


def _scores(values: object, name: str) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of finite numbers.")
    if not 1 <= len(values) <= _MAX_SCORES:
        raise ConfigurationError(f"{name} must hold 1 to {_MAX_SCORES} values.")
    numbers: list[float] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigurationError(f"{name} must contain only finite numbers.")
        number = float(value)
        if not math.isfinite(number):
            raise ConfigurationError(f"{name} must contain only finite numbers.")
        numbers.append(number)
    return numbers


def _probabilities(counts: list[list[float]]) -> list[list[float]]:
    return [[value / sum(row) for value in row] for row in counts]


def _bed_hits(path: Path) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                raise BackendExecutionError(f"{path.name} row has {len(fields)} columns.")
            hits.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "motif": fields[3],
                    "score": float(fields[4]),
                    "strand": fields[5],
                }
            )
    hits.sort(key=lambda hit: (hit["chrom"], hit["start"], hit["end"], hit["motif"], hit["strand"]))
    return hits


def fasta_record_count(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Count the sequence records in a FASTA file with gimmemotifs ``number_of_seqs_in_file``."""

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    fasta = input_path(path)
    count = _call(utils.number_of_seqs_in_file, str(fasta))
    return {
        "status": "completed",
        "operation": "fasta_record_count",
        "file": fasta.name,
        "record_count": int(count),
        "provenance": _gimme_provenance(function="gimmemotifs.utils.number_of_seqs_in_file"),
    }


def equal_length_motif_pair(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    pos: int = 0,
    truncate: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shift and pad two PFMs to equal length with gimmemotifs ``make_equal_length``.

    Counts are normalised to probabilities first. ``pos`` shifts ``pfm_b`` against ``pfm_a``;
    missing positions are filled with the uniform background (0.25 per base). ``truncate`` is
    ``None``, ``"first"``, ``"second"`` or ``"both"`` and removes the padded ends of one or both
    matrices instead of keeping them.
    """

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    rows_a = _probabilities(_pfm(pfm_a, "pfm_a"))
    rows_b = _probabilities(_pfm(pfm_b, "pfm_b"))
    shift = _integer(pos, "pos", -100, 100)
    if truncate is not None:
        _choose(truncate, ("first", "second", "both"), "truncate")
    aligned_a, aligned_b = _call(
        utils.make_equal_length,
        np.asarray(rows_a, dtype=float),
        np.asarray(rows_b, dtype=float),
        shift,
        truncate=truncate,
    )
    left = [[float(value) for value in row] for row in aligned_a]
    right = [[float(value) for value in row] for row in aligned_b]
    if len(left) != len(right):
        raise BackendExecutionError("gimmemotifs returned matrices of different lengths.")
    return {
        "status": "completed",
        "operation": "equal_length_motif_pair",
        "length": len(left),
        "pos": shift,
        "truncate": truncate,
        "pfm_a": left,
        "pfm_b": right,
        "provenance": _gimme_provenance(
            function="gimmemotifs.utils.make_equal_length",
            normalisation="rows divided by their sum before alignment",
        ),
    }


def combine_summit_files(
    summit_files: Sequence[str | Path],
    genome_fasta: str | Path,
    *,
    window: int = 200,
    scale_value: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge MACS2-style summit files and keep the strongest summit per window with gimmemotifs.

    ``genome_fasta`` is a FASTA genome; ``window`` is the width used to merge nearby summits.
    With ``scale_value`` the summit values are rescaled per file before the maximum is taken.
    """

    preprocessing = require("gimmemotifs.preprocessing", _GIMMEMOTIFS_EXTRA)
    if isinstance(summit_files, (str, bytes, Path)) or not isinstance(summit_files, Sequence):
        raise ConfigurationError("summit_files must be a list of summit file paths.")
    if not 1 <= len(summit_files) <= _MAX_SUMMIT_FILES:
        raise ConfigurationError(f"summit_files must hold 1 to {_MAX_SUMMIT_FILES} files.")
    paths = [str(input_path(value)) for value in summit_files]
    genome = input_path(genome_fasta)
    width = _integer(window, "window", 1, 100_000)
    scaled = _flag(scale_value, "scale_value")
    frame = _call(
        preprocessing.combine_peaks,
        paths,
        str(genome),
        window=width,
        scale_value=scaled,
    )
    summits: list[dict[str, Any]] = []
    for chrom, start, end in zip(
        frame["chrom"].tolist(), frame["start"].tolist(), frame["end"].tolist(), strict=True
    ):
        summits.append({"chrom": str(chrom), "start": int(start), "end": int(end)})
    if len(summits) > _MAX_SUMMITS:
        raise BackendExecutionError("gimmemotifs returned more summits than the limit.")
    summits.sort(key=lambda item: (item["chrom"], item["start"], item["end"]))
    return {
        "status": "completed",
        "operation": "combine_summit_files",
        "file_count": len(paths),
        "window": width,
        "scale_value": scaled,
        "summit_count": len(summits),
        "summits": summits,
        "provenance": _gimme_provenance(
            function="gimmemotifs.preprocessing.combine_peaks",
            genome="FASTA chromosome sizes read from the given genome file",
        ),
    }


def peak_coverage_table(
    peak_file: str | Path,
    coverage_files: Sequence[str | Path],
    *,
    window: int = 200,
    log_transform: bool = True,
    top: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read coverage around each peak summit from bigWig tracks with gimmemotifs ``coverage_table``.

    Each peak is widened to ``window`` bases around its summit and the mean coverage of every
    track is reported. ``log_transform`` applies log(1 + x). ``top`` keeps the highest-variance
    peaks when it is positive; ``-1`` keeps all peaks. Peaks are listed in input order.
    """

    preprocessing = require("gimmemotifs.preprocessing", _GIMMEMOTIFS_EXTRA)
    peaks = input_path(peak_file)
    if isinstance(coverage_files, (str, bytes, Path)) or not isinstance(coverage_files, Sequence):
        raise ConfigurationError("coverage_files must be a list of bigWig paths.")
    if not 1 <= len(coverage_files) <= _MAX_COVERAGE_TRACKS:
        raise ConfigurationError(f"coverage_files must hold 1 to {_MAX_COVERAGE_TRACKS} tracks.")
    tracks = [str(input_path(value)) for value in coverage_files]
    width = _integer(window, "window", 1, 100_000)
    logged = _flag(log_transform, "log_transform")
    if top != -1:
        _integer(top, "top", 1, 1_000_000)
    frame = _call(
        preprocessing.coverage_table,
        str(peaks),
        tracks,
        window=width,
        log_transform=logged,
        normalization=None,
        top=top,
        topmethod="var",
        rmdup=True,
        rmrepeats=True,
        ncpus=1,
    )
    columns = [str(name) for name in frame.columns]
    order = sorted(range(len(columns)), key=lambda index: columns[index])
    values = [[float(row[index]) for index in order] for row in frame.to_numpy().tolist()]
    return {
        "status": "completed",
        "operation": "peak_coverage_table",
        "window": width,
        "log_transform": logged,
        "top": top,
        "peak_count": len(values),
        "tracks": [columns[index] for index in order],
        "peaks": [str(name) for name in frame.index.tolist()],
        "values": values,
        "provenance": _gimme_provenance(
            function="gimmemotifs.preprocessing.coverage_table",
            normalization="none; log(1 + x) when log_transform is true",
        ),
    }


def scan_motifs_to_file(
    sequences_fasta: str | Path,
    motifs_pfm: str | Path,
    output_bed: str | Path,
    *,
    cutoff: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write every motif hit above ``cutoff`` to a BED file with gimmemotifs ``scan_to_file``.

    ``cutoff`` is a fraction from 0 to 1 of the gimmemotifs log-odds threshold; ``0`` keeps all
    hits. Scores are not z-normalised and no GC correction is applied. The output file must not
    exist. Hits are returned in coordinate order.
    """

    scanner = require("gimmemotifs.scanner", _GIMMEMOTIFS_EXTRA)
    fasta = input_path(sequences_fasta)
    pfm = input_path(motifs_pfm)
    output = Path(output_bed).expanduser().resolve()
    if output.exists():
        raise ConfigurationError(f"Output file already exists: {output}")
    if not output.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {output.parent}")
    threshold = _fraction(cutoff, "cutoff")
    _call(
        scanner.scan_to_file,
        str(fasta),
        pfmfile=str(pfm),
        filepath_or_buffer=str(output),
        cutoff=threshold,
        bed=True,
        ncpus=1,
        zscore=False,
        gcnorm=False,
        progress=False,
    )
    if not output.is_file():
        raise BackendExecutionError("gimmemotifs did not write the BED output file.")
    hits = _bed_hits(output)
    return {
        "status": "completed",
        "operation": "scan_motifs_to_file",
        "cutoff": threshold,
        "output": output.name,
        "hit_count": len(hits),
        "hits": hits,
        "provenance": _gimme_provenance(
            function="gimmemotifs.scanner.scan_to_file",
            score="gimmemotifs log-odds score; no z-score and no GC correction",
            coordinates="0-based start, exclusive end, strand of the hit",
        ),
    }


def scan_region_scores(
    regions_file: str | Path,
    genome_fasta: str | Path,
    motifs_pfm: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score every motif in every region of a list with gimmemotifs ``scan_regionfile_to_table``.

    ``regions_file`` holds one ``chrom:start-end`` region per line; sequences come from the
    FASTA genome whose contigs must be longer than 1,280 bases (gimmemotifs samples 128-base
    background windows from them). The table holds the best log-odds score of each motif in each
    region. Z-scores and GC normalisation are not applied, so no GC background is needed.
    """

    scanner = require("gimmemotifs.scanner", _GIMMEMOTIFS_EXTRA)
    regions = input_path(regions_file)
    genome = input_path(genome_fasta)
    pfm = input_path(motifs_pfm)
    names: list[str] = []
    with regions.open("rt", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            text = raw.strip()
            if not text:
                continue
            match = _REGION.match(text)
            if match is None or int(match.group(3)) <= int(match.group(2)):
                raise ConfigurationError(f"{regions.name} line {number} is not chrom:start-end.")
            if int(match.group(3)) - int(match.group(2)) > _MAX_REGION_LENGTH:
                raise ConfigurationError(f"{regions.name} line {number} is too long.")
            names.append(text)
            if len(names) > _MAX_REGIONS:
                raise ConfigurationError(f"{regions.name} has more than {_MAX_REGIONS} regions.")
    if not names:
        raise ConfigurationError(f"{regions.name} contains no regions.")
    with tempfile.TemporaryDirectory(prefix="gimme_regions_") as scratch:
        # gimmemotifs reads the first line of the table as a header, so one is added here.
        table = Path(scratch) / "regions.tsv"
        table.write_text("region\n" + "\n".join(names) + "\n", encoding="utf-8")
        frame = _call(
            scanner.scan_regionfile_to_table,
            str(table),
            str(genome),
            "score",
            pfmfile=str(pfm),
            ncpus=1,
            zscore=False,
            gc=False,
            progress=False,
        )
    motifs = sorted(str(column) for column in frame.columns)
    order = sorted(range(len(frame.columns)), key=lambda index: str(frame.columns[index]))
    rows: list[tuple[str, list[float]]] = []
    for region, values in zip(frame.index.tolist(), frame.to_numpy().tolist(), strict=True):
        scored = [float(values[index]) for index in order]
        if not all(math.isfinite(value) for value in scored):
            raise BackendExecutionError("gimmemotifs returned a non-finite region score.")
        rows.append((str(region), scored))
    rows.sort(key=lambda item: item[0])
    return {
        "status": "completed",
        "operation": "scan_region_scores",
        "region_count": len(rows),
        "motif_count": len(motifs),
        "motifs": motifs,
        "regions": [region for region, _ in rows],
        "scores": [scored for _, scored in rows],
        "provenance": _gimme_provenance(
            function="gimmemotifs.scanner.scan_regionfile_to_table",
            scoring="best log-odds score per motif and region",
            normalisation="no z-score and no GC correction",
        ),
    }


def motif_activity_scores(
    labels_file: str | Path,
    motif_counts_file: str | Path,
    *,
    method: str = "hypergeom",
    verbose: bool = True,
) -> dict[str, Any]:
    """Motif activity per motif from a precomputed count table with gimmemotifs ``moap``.

    ``labels_file`` is a tab-separated table with a header line, indexed by region name, with one
    value column (1 for the foreground, 0 for the background). ``motif_counts_file`` is a
    region-by-motif table of non-negative integer hit counts with a header line and the same index.
    Only the ``hypergeom`` method is exposed here. The result has one row per motif and one column
    per label class.
    """

    moap = require("gimmemotifs.maelstrom.moap", _GIMMEMOTIFS_EXTRA)
    chosen = _choose(method, ("hypergeom",), "method")
    labels = input_path(labels_file)
    counts = input_path(motif_counts_file)
    frame = _call(
        moap.moap,
        str(labels),
        method=chosen,
        scoring="count",
        outfile=None,
        motiffile=str(counts),
        ncpus=1,
        zscore=False,
        gc=False,
        progress=False,
    )
    classes = [str(name) for name in frame.columns]
    class_order = sorted(range(len(classes)), key=lambda index: classes[index])
    motif_rows = sorted(
        zip([str(name) for name in frame.index.tolist()], frame.to_numpy().tolist(), strict=True),
        key=lambda item: item[0],
    )
    values = [[float(row[index]) for index in class_order] for _, row in motif_rows]
    if not all(math.isfinite(value) for row in values for value in row):
        raise BackendExecutionError("gimmemotifs returned a non-finite activity value.")
    return {
        "status": "completed",
        "operation": "motif_activity_scores",
        "method": chosen,
        "classes": [classes[index] for index in class_order],
        "motifs": [name for name, _ in motif_rows],
        "values": values,
        "provenance": _gimme_provenance(
            function="gimmemotifs.maelstrom.moap.moap",
            method=chosen,
            input="precomputed motif count table (motiffile); no genome scan",
        ),
    }


def motif_score_metric(
    foreground: Sequence[float],
    background: Sequence[float],
    *,
    metric: str = "roc_auc",
    fpr: float = 0.01,
    fdr: float = 0.1,
    minbg: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """One ROC, PR or enrichment statistic for foreground and background motif scores.

    ``metric`` selects a gimmemotifs ``rocmetrics`` function: ``roc_auc``, ``pr_auc``,
    ``max_fmeasure``, ``mncp``, ``max_enrichment`` (uses ``minbg``), ``recall_at_fdr`` (uses
    ``fdr``), ``matches_at_fpr``, ``phyper_at_fpr``, ``fraction_fpr``, ``score_at_fpr``,
    ``enr_at_fpr`` (use ``fpr``), ``roc_auc_xlim`` (area up to ``fpr``) and ``roc_values``
    (returns the ROC curve in ``curve``). Non-finite results are returned as None.
    """

    rocmetrics = require("gimmemotifs.rocmetrics", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    chosen = _choose(metric, _METRICS, "metric")
    fg = np.asarray(_scores(foreground, "foreground"), dtype=float)
    bg = np.asarray(_scores(background, "background"), dtype=float)
    fpr_value = _fraction(fpr, "fpr")
    if fpr_value <= 0.0:
        raise ConfigurationError("fpr must be greater than 0.")
    fdr_value = _fraction(fdr, "fdr")
    if fdr_value <= 0.0:
        raise ConfigurationError("fdr must be greater than 0.")
    minimum = _integer(minbg, "minbg", 1, 1_000)
    provenance = _gimme_provenance(
        function=f"gimmemotifs.rocmetrics.{chosen}",
        input="foreground and background score vectors",
    )
    base: dict[str, Any] = {
        "status": "completed",
        "operation": "motif_score_metric",
        "metric": chosen,
        "foreground_count": int(fg.size),
        "background_count": int(bg.size),
        "provenance": provenance,
    }
    if chosen == "roc_values":
        x_values, y_values = _call(rocmetrics.roc_values, fg, bg)
        base["curve"] = {
            "fpr": [float(value) for value in x_values],
            "tpr": [float(value) for value in y_values],
        }
        return base
    if chosen == "roc_auc_xlim":
        x_values, y_values = _call(rocmetrics.roc_values, fg, bg)
        value = _call(rocmetrics.roc_auc_xlim, x_values, y_values, xlim=fpr_value)
    elif chosen == "recall_at_fdr":
        value = _call(rocmetrics.recall_at_fdr, fg, bg, fdr_cutoff=fdr_value)
    elif chosen == "max_enrichment":
        value = _call(rocmetrics.max_enrichment, fg, bg, minbg=minimum)
    elif chosen in (
        "matches_at_fpr",
        "phyper_at_fpr",
        "fraction_fpr",
        "score_at_fpr",
        "enr_at_fpr",
    ):
        value = _call(getattr(rocmetrics, chosen), fg, bg, fpr=fpr_value)
    else:
        value = _call(getattr(rocmetrics, chosen), fg, bg)
    number = float(value)
    base["value"] = number if math.isfinite(number) else None
    return base


def motif_position_ks(
    foreground_positions: Sequence[float],
    background_positions: Sequence[float] | None = None,
    *,
    output: str = "pvalue",
    verbose: bool = True,
) -> dict[str, Any]:
    """Kolmogorov-Smirnov test of motif hit positions with gimmemotifs ``ks_pvalue``.

    ``output`` is ``"pvalue"`` for the p-value or ``"significance"`` for -log10 of the p-value.
    Positions are offsets within the sequence; without ``background_positions`` the upstream
    default comparison is used. Non-finite values are returned as None.
    """

    rocmetrics = require("gimmemotifs.rocmetrics", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    chosen = _choose(output, _OUTPUT_KEYS, "output")
    fg = np.asarray(_scores(foreground_positions, "foreground_positions"), dtype=float)
    bg = None
    if background_positions is not None:
        bg = np.asarray(_scores(background_positions, "background_positions"), dtype=float)
    function = rocmetrics.ks_pvalue if chosen == "pvalue" else rocmetrics.ks_significance
    value = float(_call(function, fg, bg))
    return {
        "status": "completed",
        "operation": "motif_position_ks",
        "output": chosen,
        "foreground_count": int(fg.size),
        "background_count": None if bg is None else int(bg.size),
        "value": value if math.isfinite(value) else None,
        "provenance": _gimme_provenance(
            function=f"gimmemotifs.rocmetrics.{function.__name__}",
        ),
    }


__all__ = [
    "combine_summit_files",
    "equal_length_motif_pair",
    "fasta_record_count",
    "motif_activity_scores",
    "motif_position_ks",
    "motif_score_metric",
    "peak_coverage_table",
    "scan_motifs_to_file",
    "scan_region_scores",
]

_register_public_verbose(globals(), tuple(__all__))
