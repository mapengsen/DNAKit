"""Motif-matrix, motif-scan and BED-interval functions over gimmemotifs.

Each function calls one gimmemotifs routine on local files or explicit lists. Genome-dependent
routines take a FASTA file as the genome; no genome, background or motif database is downloaded.
Motif matrices are PFM counts with one row per position and columns A, C, G and T.
"""

from __future__ import annotations

import math
import random
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _choose,
    _dna,
    _flag,
    _gimme_provenance,
    _integer,
    _nested,
    _pfm,
)

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_IUPAC = "ACGTRYSWKMBDHVN"
_MAX_CONSENSUS = 100
_MAX_ALIGNMENT_ROWS = 10_000
_MAX_MOTIF_FILE_ROWS = 1_000
_MAX_SEQUENCES = 100_000
_MAX_INTERVAL_ROWS = 1_000_000
_METRICS = ("akl", "chisq", "ssd")
_CLUSTER_METRICS = ("wic", "pcc", "ed", "distance", "chisq")
_CLUSTER_MATCH = ("total", "subtotal")


def _call(package: str, function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:  # gimmemotifs raises several non-DNAKit exception types
        raise BackendExecutionError(f"{package} failed: {exc}") from exc


def _count(value: object) -> int | float:
    number = float(value)  # type: ignore[arg-type]
    return int(number) if number.is_integer() else number


def _matrix_rows(motif: Any) -> list[list[int | float]]:
    return [[_count(value) for value in row] for row in motif.pfm]


def _consensus(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("consensus must be a non-empty IUPAC string.")
    text = value.strip().upper()
    if len(text) > _MAX_CONSENSUS:
        raise ConfigurationError(f"consensus is longer than {_MAX_CONSENSUS} symbols.")
    if set(text) - set(_IUPAC):
        raise ConfigurationError("consensus contains symbols outside IUPAC nucleotide codes.")
    return text


def _ppm(counts: Sequence[Sequence[float]]) -> list[list[float]]:
    """Normalise PFM counts to probabilities and add gimmemotifs' small pseudocount."""

    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    array = np.asarray(counts, dtype=float)
    probabilities = array / array.sum(axis=1, keepdims=True)
    return [[float(v) for v in row] for row in utils.ppm_pseudocount(probabilities)]


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _fasta_names(path: Path) -> list[str]:
    """Sequence names in file order, using the first token after '>' as gimmemotifs does."""

    names: list[str] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                tokens = line[1:].split()
                if not tokens:
                    raise ConfigurationError(f"{path.name} has a FASTA header without a name.")
                names.append(tokens[0])
    if not names:
        raise ConfigurationError(f"{path.name} contains no FASTA records.")
    return names


def _bed_rows(path: Path, minimum: int, source: str) -> list[list[str]]:
    """Read BED-like rows, skipping track, browser and blank lines."""

    rows: list[list[str]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("track", "browser")):
                continue
            if line.startswith("#"):
                raise ConfigurationError(f"{source} line {number} is a comment; remove it.")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < minimum:
                raise ConfigurationError(
                    f"{source} line {number} has {len(fields)} columns; at least {minimum} needed."
                )
            rows.append(fields)
            if len(rows) > _MAX_INTERVAL_ROWS:
                raise ConfigurationError(f"{source} has more than {_MAX_INTERVAL_ROWS} rows.")
    if not rows:
        raise ConfigurationError(f"{source} contains no intervals.")
    return rows


def _interval(fields: list[str], source: str, number: int) -> dict[str, Any]:
    try:
        start, end = int(fields[1]), int(fields[2])
    except ValueError as exc:
        raise BackendExecutionError(
            f"{source} row {number} has non-integer coordinates: {exc}"
        ) from exc
    return {"chrom": fields[0], "start": start, "end": end, "columns": fields[3:]}


def _intervals_from(path: Path, source: str) -> list[dict[str, Any]]:
    rows = _bed_rows(path, 3, source)
    return [_interval(fields, source, number) for number, fields in enumerate(rows, start=1)]


def _motif_length_by_id(motif_module: Any, path: Path) -> dict[str, int]:
    parsed = _call(
        "gimmemotifs",
        lambda: motif_module.read_motifs(str(path), fmt="pfm", as_dict=False),
    )
    return {str(motif.id): len(motif.pfm) for motif in parsed}


def motif_from_consensus(
    consensus: str, *, n_sites: int = 1200, verbose: bool = True
) -> dict[str, Any]:
    """Build a PFM from an IUPAC consensus string with gimmemotifs ``motif_from_consensus``.

    ``n_sites`` is the total count per position; a base named by the consensus gets that count,
    and an ambiguous code spreads it over its bases.
    """

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    text = _consensus(consensus)
    total = _integer(n_sites, "n_sites", 1, 1_000_000)
    motif = _call("gimmemotifs", lambda: motif_module.motif_from_consensus(text, n=total))
    return {
        "status": "completed",
        "operation": "motif_from_consensus",
        "consensus_input": text,
        "consensus": str(motif.consensus),
        "length": len(motif.pfm),
        "pfm": _matrix_rows(motif),
        "provenance": _gimme_provenance(
            function="gimmemotifs.motif.motif_from_consensus",
            counts=f"n_sites={total} per position, as gimmemotifs stores it",
        ),
    }


def motif_from_alignment(sequences: Sequence[str], *, verbose: bool = True) -> dict[str, Any]:
    """Build a PFM from equal-length aligned DNA strings with gimmemotifs ``motif_from_align``.

    The result has the base counts of the alignment (``pfm``) and the probabilities that
    gimmemotifs derives from them (``ppm``).
    """

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    if isinstance(sequences, (str, bytes)) or not isinstance(sequences, Sequence):
        raise ConfigurationError("sequences must be a list of aligned DNA strings.")
    if len(sequences) < 2 or len(sequences) > _MAX_ALIGNMENT_ROWS:
        raise ConfigurationError(
            f"an alignment needs 2 to {_MAX_ALIGNMENT_ROWS} sequences; got {len(sequences)}."
        )
    texts = [_dna(value, f"sequence {index}", 1_000) for index, value in enumerate(sequences)]
    if len({len(text) for text in texts}) != 1:
        raise ConfigurationError("aligned sequences must all have the same length.")
    motif = _call("gimmemotifs", lambda: motif_module.motif_from_align(texts))
    return {
        "status": "completed",
        "operation": "motif_from_alignment",
        "sequence_count": len(texts),
        "consensus": str(motif.consensus),
        "length": len(motif.pfm),
        "pfm": _matrix_rows(motif),
        "ppm": _nested([[round(float(v), 6) for v in row] for row in motif.ppm]),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.motif_from_align"),
    }


def compare_motif_positions(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    metric: str = "akl",
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two equal-length PFMs position by position with a gimmemotifs position metric.

    ``metric`` is ``"akl"`` (average Kullback-Leibler similarity), ``"chisq"`` or ``"ssd"``
    (sum of squared distances). Both matrices are normalised to probabilities with gimmemotifs'
    pseudocount first, so no position has a zero probability.
    """

    comparison = require("gimmemotifs.comparison", _GIMMEMOTIFS_EXTRA)
    selected = _choose(metric, _METRICS, "metric")
    rows_a = _ppm(_pfm(pfm_a, "pfm_a"))
    rows_b = _ppm(_pfm(pfm_b, "pfm_b"))
    if len(rows_a) != len(rows_b):
        raise ConfigurationError("pfm_a and pfm_b must have the same length.")
    scorer = getattr(comparison, selected)
    scores = []
    for left, right in zip(rows_a, rows_b, strict=True):
        value = float(_call("gimmemotifs", scorer, left, right))
        if not math.isfinite(value):
            raise BackendExecutionError(f"gimmemotifs {selected} returned a non-finite score.")
        scores.append(value)
    return {
        "status": "completed",
        "operation": "compare_motif_positions",
        "metric": selected,
        "length": len(scores),
        "position_scores": scores,
        "mean_score": sum(scores) / len(scores),
        "provenance": _gimme_provenance(
            function=f"gimmemotifs.comparison.{selected}",
            normalisation="rows divided by their sum, then gimmemotifs ppm_pseudocount",
        ),
    }


def motif_score_correlation(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pearson correlation of two PFMs' scores over a de Bruijn scan, with gimmemotifs ``seqcor``.

    Returns the best correlation, the offset of the second motif against the first, and the
    strand (1 or -1) at which it is found.
    """

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    comparison = require("gimmemotifs.comparison", _GIMMEMOTIFS_EXTRA)
    matrix_a = _pfm(pfm_a, "pfm_a")
    matrix_b = _pfm(pfm_b, "pfm_b")
    motif_a = motif_module.Motif(matrix_a)
    motif_b = motif_module.Motif(matrix_b)
    score, position, strand = _call("gimmemotifs", lambda: comparison.seqcor(motif_a, motif_b))
    if not math.isfinite(float(score)):
        raise BackendExecutionError("gimmemotifs seqcor returned a non-finite score.")
    return {
        "status": "completed",
        "operation": "motif_score_correlation",
        "length_a": len(matrix_a),
        "length_b": len(matrix_b),
        "score": float(score),
        "position": int(position),
        "strand": int(strand),
        "provenance": _gimme_provenance(function="gimmemotifs.comparison.seqcor"),
    }


def scan_best_motif_hits(
    sequences: str | Path,
    motifs: str | Path,
    *,
    max_sequences: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Best log-odds hit per PFM motif and FASTA sequence (gimmemotifs ``scan_to_best_match``).

    ``motifs`` is a PFM file. Each sequence and motif gives one hit: the highest score on either
    strand, its 0-based forward start, its exclusive end and its strand. Hits are sorted by
    sequence order in the FASTA file, then by motif ID.
    """

    scanner = require("gimmemotifs.scanner", _GIMMEMOTIFS_EXTRA)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    fasta = _existing_file(sequences, "sequences")
    motif_file = _existing_file(motifs, "motifs")
    limit = _integer(max_sequences, "max_sequences", 1, _MAX_SEQUENCES)
    names = _fasta_names(fasta)
    if len(names) > limit:
        raise ConfigurationError(f"{fasta.name} has more than max_sequences={limit} records.")
    lengths = _motif_length_by_id(motif_module, motif_file)
    if not lengths:
        raise ConfigurationError(f"{motif_file.name} contains no PFM motifs.")
    found = _call(
        "gimmemotifs",
        lambda: scanner.scan_to_best_match(
            str(fasta),
            pfmfile=str(motif_file),
            ncpus=1,
            score=False,
            zscore=False,
            gc=False,
            progress=False,
        ),
    )
    hits: list[dict[str, Any]] = []
    for motif_id, rows in found.items():
        if motif_id not in lengths or len(rows) != len(names):
            raise BackendExecutionError(
                f"gimmemotifs returned an unexpected result for {motif_id}."
            )
        for index, (score, start, strand) in enumerate(rows):
            if not math.isfinite(float(score)):
                raise BackendExecutionError(
                    f"gimmemotifs returned a non-finite score for {motif_id}."
                )
            length = lengths[str(motif_id)]
            hits.append(
                {
                    "sequence_index": index,
                    "sequence": names[index],
                    "motif": str(motif_id),
                    "score": float(score),
                    "start": int(start),
                    "end": int(start) + length,
                    "strand": "+" if int(strand) > 0 else "-",
                }
            )
    hits.sort(key=lambda hit: (hit["sequence_index"], hit["motif"]))
    for hit in hits:
        del hit["sequence_index"]
    return {
        "status": "completed",
        "operation": "scan_best_motif_hits",
        "sequence_count": len(names),
        "motif_count": len(lengths),
        "hit_count": len(hits),
        "hits": hits,
        "provenance": _gimme_provenance(
            function="gimmemotifs.scanner.scan_to_best_match",
            score="gimmemotifs log-odds: sum of ln(ppm / 0.25 + 0.01) over motif positions",
            coordinates="0-based forward start; minus-strand hits keep the forward start",
        ),
    }


def motif_enrichment_stats(
    foreground: str | Path,
    background: str | Path,
    motifs: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """ROC, PR, recall and enrichment statistics per motif, with gimmemotifs ``calc_stats``.

    ``foreground`` and ``background`` are FASTA files; ``motifs`` is a PFM file. Motif scores
    are not z-normalised and no GC correction is applied, so no genome is required.
    """

    stats_module = require("gimmemotifs.stats", _GIMMEMOTIFS_EXTRA)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    fg_file = _existing_file(foreground, "foreground")
    bg_file = _existing_file(background, "background")
    motif_file = _existing_file(motifs, "motifs")
    _fasta_names(fg_file)
    _fasta_names(bg_file)
    if not _motif_length_by_id(motif_module, motif_file):
        raise ConfigurationError(f"{motif_file.name} contains no PFM motifs.")
    raw = _call(
        "gimmemotifs",
        lambda: stats_module.calc_stats(
            fg_file=str(fg_file),
            bg_file=str(bg_file),
            motifs=str(motif_file),
            zscore=False,
            gc=False,
            ncpus=1,
        ),
    )
    motif_rows = []
    for motif_id in sorted(raw):
        values = raw[motif_id]
        stats = {}
        for name in sorted(values):
            number = float(values[name])
            stats[name] = number if math.isfinite(number) else None
        motif_rows.append({"motif": str(motif_id), "stats": stats})
    return {
        "status": "completed",
        "operation": "motif_enrichment_stats",
        "motif_count": len(motif_rows),
        "motifs": motif_rows,
        "provenance": _gimme_provenance(
            function="gimmemotifs.stats.calc_stats",
            scores="motif log-odds scores of each FASTA sequence",
            normalisation="no z-score and no GC correction",
        ),
    }


def cluster_motif_set(
    motifs: str | Path,
    *,
    metric: str = "wic",
    match: str = "total",
    threshold: float = 0.95,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hierarchically cluster PFM motifs by similarity, with gimmemotifs ``cluster_motifs``.

    At least two motifs are needed. ``metric`` is one of wic, pcc, ed, distance or chisq;
    ``match`` is ``"total"`` or ``"subtotal"``; ``threshold`` (0 to 1) is the merge cut-off.
    Each cluster lists its member motif IDs and the ID of its averaged representative.
    """

    cluster_module = require("gimmemotifs.motif.cluster", _GIMMEMOTIFS_EXTRA)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    motif_file = _existing_file(motifs, "motifs")
    chosen_metric = _choose(metric, _CLUSTER_METRICS, "metric")
    chosen_match = _choose(match, _CLUSTER_MATCH, "match")
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
        raise ConfigurationError("threshold must be a number in [0, 1].")
    cut = float(threshold)
    if not math.isfinite(cut) or not 0.0 <= cut <= 1.0:
        raise ConfigurationError("threshold must be a number in [0, 1].")
    count = len(_motif_length_by_id(motif_module, motif_file))
    if count < 2:
        raise ConfigurationError("clustering needs at least two PFM motifs.")
    root = _call(
        "gimmemotifs",
        lambda: cluster_module.cluster_motifs(
            str(motif_file),
            match=chosen_match,
            metric=chosen_metric,
            combine="mean",
            pval=True,
            threshold=cut,
            progress=False,
            ncpus=1,
        ),
    )
    clusters = []
    for representative, members in root.getResult():
        ids = sorted(str(motif.id) for motif in members)
        clusters.append({"representative": str(representative.id), "members": ids})
    clusters.sort(key=lambda item: item["members"][0])
    return {
        "status": "completed",
        "operation": "cluster_motif_set",
        "motif_count": count,
        "cluster_count": len(clusters),
        "clusters": clusters,
        "provenance": _gimme_provenance(
            function="gimmemotifs.motif.cluster.cluster_motifs",
            metric=chosen_metric,
            match=chosen_match,
            threshold=cut,
        ),
    }


def read_peak_summits(
    path: str | Path, *, max_peaks: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """Read MACS2 summit or narrowPeak rows with gimmemotifs ``read_peak_file_to_df``.

    Summit files have five columns (chromosome, start, end, name, value). narrowPeak files
    keep their ten columns. Other BED files are rejected by gimmemotifs.
    """

    preprocessing = require("gimmemotifs.preprocessing", _GIMMEMOTIFS_EXTRA)
    source = _existing_file(path, "path")
    limit = _integer(max_peaks, "max_peaks", 1, _MAX_INTERVAL_ROWS)
    try:
        frame = preprocessing.read_peak_file_to_df(str(source))
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(
            f"{source.name} is not a MACS2 summit or narrowPeak file: {exc}"
        ) from exc
    except BaseException as exc:
        raise BackendExecutionError(f"gimmemotifs failed: {exc}") from exc
    if len(frame) > limit:
        raise ConfigurationError(f"{source.name} has more than max_peaks={limit} peaks.")
    peaks = []
    for row in frame.to_dict(orient="records"):
        clean = {}
        for key, value in row.items():
            if hasattr(value, "item"):
                value = value.item()
            clean[str(key)] = (
                value if not isinstance(value, float) or math.isfinite(value) else None
            )
        peaks.append(clean)
    return {
        "status": "completed",
        "operation": "read_peak_summits",
        "peak_count": len(peaks),
        "columns": [str(column) for column in frame.columns],
        "peaks": peaks,
        "provenance": _gimme_provenance(function="gimmemotifs.preprocessing.read_peak_file_to_df"),
    }


def summit_intervals_from_narrowpeak(
    path: str | Path,
    *,
    size: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fixed-width intervals centred on each narrowPeak summit (gimmemotifs ``narrowpeak_to_bed``).

    The summit is column 10 of the narrowPeak row. Each output interval is ``size`` bases long
    and starts at ``summit - size // 2`` in the peak; ``columns`` keeps the upstream BED columns.
    """

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    source = _existing_file(path, "path")
    width = _integer(size, "size", 1, 1_000_000)
    rows = _bed_rows(source, 10, source.name)
    for number, fields in enumerate(rows, start=1):
        try:
            summit = int(fields[9])
        except ValueError as exc:
            raise ConfigurationError(
                f"{source.name} row {number} has a non-integer summit."
            ) from exc
        if summit < -1:
            raise ConfigurationError(f"{source.name} row {number} has a negative summit.")
    with tempfile.TemporaryDirectory(prefix="dnakit-gimmemotifs-") as scratch:
        output = Path(scratch) / "summits.bed"
        _call("gimmemotifs", lambda: utils.narrowpeak_to_bed(str(source), str(output), size=width))
        intervals = _intervals_from(output, "gimmemotifs output")
    return {
        "status": "completed",
        "operation": "summit_intervals_from_narrowpeak",
        "size": width,
        "interval_count": len(intervals),
        "intervals": intervals,
        "provenance": _gimme_provenance(
            function="gimmemotifs.utils.narrowpeak_to_bed",
            anchor="summit (narrowPeak column 10); -1 uses the peak centre",
        ),
    }


def equal_size_intervals(
    path: str | Path,
    *,
    size: int,
    keep_columns: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resize BED intervals to ``size`` bases around the midpoint (``write_equalsize_bedfile``).

    With ``keep_columns=False`` only the first six BED columns are returned; the upstream
    ``peakfile2bedfile`` writes the same midpoint intervals in that six-column form.
    """

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    source = _existing_file(path, "path")
    width = _integer(size, "size", 1, 1_000_000)
    flag = _flag(keep_columns, "keep_columns")
    rows = _bed_rows(source, 3, source.name)
    for number, fields in enumerate(rows, start=1):
        try:
            if int(fields[1]) >= int(fields[2]):
                raise ConfigurationError(f"{source.name} row {number} must have start < end.")
        except ValueError as exc:
            raise ConfigurationError(
                f"{source.name} row {number} has non-integer coordinates."
            ) from exc
    with tempfile.TemporaryDirectory(prefix="dnakit-gimmemotifs-") as scratch:
        output = Path(scratch) / "equal.bed"
        _call("gimmemotifs", lambda: utils.write_equalsize_bedfile(str(source), width, str(output)))
        intervals = _intervals_from(output, "gimmemotifs output")
    if not flag:
        for item in intervals:
            item["columns"] = item["columns"][:3]
    return {
        "status": "completed",
        "operation": "equal_size_intervals",
        "size": width,
        "keep_columns": flag,
        "interval_count": len(intervals),
        "intervals": intervals,
        "provenance": _gimme_provenance(
            function="gimmemotifs.utils.write_equalsize_bedfile",
            anchor="midpoint of each interval; start = midpoint - size // 2",
        ),
    }


def random_genomic_regions(
    genome: str | Path,
    size: int,
    n: int,
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Random fixed-width genome intervals (gimmemotifs ``create_random_genomic_bedfile``).

    Intervals are sorted by chromosome and start. ``seed`` is passed to gimmemotifs; with
    ``seed=None`` the intervals differ between calls.
    """

    background = require("gimmemotifs.background", _GIMMEMOTIFS_EXTRA)
    fasta = _existing_file(genome, "genome")
    width = _integer(size, "size", 1, 10_000_000)
    count = _integer(n, "n", 1, 100_000)
    chosen = None if seed is None else _integer(seed, "seed", 0, 2**32 - 1)
    with tempfile.TemporaryDirectory(prefix="dnakit-gimmemotifs-") as scratch:
        output = Path(scratch) / "random.bed"
        _call(
            "gimmemotifs",
            lambda: background.create_random_genomic_bedfile(
                str(output), str(fasta), width, count, seed=chosen
            ),
        )
        intervals = sorted(
            _intervals_from(output, "gimmemotifs output"),
            key=lambda item: (item["chrom"], item["start"], item["end"]),
        )
    return {
        "status": "completed",
        "operation": "random_genomic_regions",
        "size": width,
        "interval_count": len(intervals),
        "intervals": intervals,
        "provenance": _gimme_provenance(
            function="gimmemotifs.background.create_random_genomic_bedfile",
            genome=fasta.name,
            seed="none" if chosen is None else chosen,
        ),
    }


def promoter_regions(
    genes: str | Path,
    size: int,
    n: int,
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample promoter windows upstream of gene starts (gimmemotifs ``create_promoter_bedfile``).

    ``genes`` is a six-column BED file (chromosome, start, end, name, score, strand). A window
    of ``size`` bases lies before the start on the plus strand and after the end on the minus
    strand. Windows are sampled without replacement when there are more than ``n``. gimmemotifs
    draws from Python's global random module, so the state is saved and restored around the
    call; ``seed=None`` leaves that state untouched and the sample is not reproducible.
    """

    background = require("gimmemotifs.background", _GIMMEMOTIFS_EXTRA)
    source = _existing_file(genes, "genes")
    width = _integer(size, "size", 1, 10_000_000)
    count = _integer(n, "n", 1, 100_000)
    chosen = None if seed is None else _integer(seed, "seed", 0, 2**32 - 1)
    rows = _bed_rows(source, 6, source.name)
    for number, fields in enumerate(rows, start=1):
        if fields[5] not in {"+", "-", "1", "-1"}:
            raise ConfigurationError(f"{source.name} row {number} needs strand +, - , 1 or -1.")
        try:
            int(fields[1])
            int(fields[2])
        except ValueError as exc:
            raise ConfigurationError(
                f"{source.name} row {number} has non-integer coordinates."
            ) from exc
    state = random.getstate()
    try:
        if chosen is not None:
            random.seed(chosen)
        with tempfile.TemporaryDirectory(prefix="dnakit-gimmemotifs-") as scratch:
            output = Path(scratch) / "promoters.bed"
            _call(
                "gimmemotifs",
                lambda: background.create_promoter_bedfile(str(output), str(source), width, count),
            )
            intervals = sorted(
                _intervals_from(output, "gimmemotifs output"),
                key=lambda item: (item["chrom"], item["start"], item["end"]),
            )
    finally:
        random.setstate(state)
    promoters = [
        {
            "chrom": item["chrom"],
            "start": item["start"],
            "end": item["end"],
            "strand": item["columns"][-1],
        }
        for item in intervals
    ]
    return {
        "status": "completed",
        "operation": "promoter_regions",
        "size": width,
        "requested": count,
        "interval_count": len(promoters),
        "intervals": promoters,
        "provenance": _gimme_provenance(
            function="gimmemotifs.background.create_promoter_bedfile",
            seed="none" if chosen is None else chosen,
            sampling="global random module, state restored after the call",
        ),
    }


__all__ = [
    "cluster_motif_set",
    "compare_motif_positions",
    "equal_size_intervals",
    "motif_enrichment_stats",
    "motif_from_alignment",
    "motif_from_consensus",
    "motif_score_correlation",
    "promoter_regions",
    "random_genomic_regions",
    "read_peak_summits",
    "scan_best_motif_hits",
    "summit_intervals_from_narrowpeak",
]

_register_public_verbose(globals(), tuple(__all__))
