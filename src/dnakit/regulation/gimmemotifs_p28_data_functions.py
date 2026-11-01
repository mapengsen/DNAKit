"""FASTA, peak-input, statistics, threshold, match-count and motif-database functions (pack P28).

Each function runs one gimmemotifs routine on local files or explicit lists and returns a
JSON-serialisable dictionary. Genome arguments are FASTA paths (no genomepy download). The motif
database functions use the network; ordinary tests mock the requests.
"""

from __future__ import annotations

import contextlib
import math
import random
import re
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import _choose, _flag, _gimme_provenance, _integer

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_MAX_SEQUENCES = 1_000_000
_MAX_MOTIF_FILE_ROWS = 10_000
_MAX_STATS_MOTIFS = 100_000
_INPUT_TYPES = ("fasta", "bed", "narrowpeak")
_DATABASES = ("jaspar", "homer", "hocomoco", "encode", "swissregulon", "image")
_HOCOMOCO_GROUPS = ("HUMAN", "MOUSE")


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures (also sys.exit) into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except BaseException as exc:  # gimmemotifs raises several non-DNAKit exception types
        raise BackendExecutionError(f"gimmemotifs failed: {type(exc).__name__}: {exc}") from exc


@contextlib.contextmanager
def _seeded(seed: int | None) -> Iterator[None]:
    """Seed Python's global random generator for one call and restore it afterwards."""

    state = random.getstate()
    try:
        if seed is not None:
            random.seed(seed)
        yield
    finally:
        random.setstate(state)


def _seed(seed: object) -> int | None:
    if seed is None:
        return None
    return _integer(seed, "seed", 0, 2**32 - 1)


def _new_output(path: str | Path, name: str) -> Path:
    destination = Path(path).expanduser()
    if destination.exists():
        raise ConfigurationError(f"{name} already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination


def _fasta(path: str | Path, name: str = "FASTA") -> Any:
    fasta_module = require("gimmemotifs.fasta", _GIMMEMOTIFS_EXTRA)
    source = input_path(path)
    fa = _call(fasta_module.Fasta, str(source))
    if not 1 <= len(fa) <= _MAX_SEQUENCES:
        raise ConfigurationError(f"{name} must hold 1 to {_MAX_SEQUENCES} sequences.")
    return fa


def gimmemotifs_hardmask_fasta(
    input_fasta: str | Path,
    output_fasta: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the lower-case (soft-masked) bases of a FASTA file by N (``Fasta.hardmask``).

    Lower-case a, c, g, t and n become ``N``; upper-case bases are kept. The output file is
    written without line wrapping and must not exist. Returns the number of masked bases.
    """

    fa = _fasta(input_fasta)
    destination = _new_output(output_fasta, "output_fasta")
    masked = sum(sum(1 for base in seq if base in "acgtn") for seq in fa.seqs)
    _call(fa.hardmask)
    _call(fa.writefasta, str(destination))
    return {
        "status": "completed",
        "operation": "hardmask_fasta",
        "n_sequences": len(fa),
        "masked_bases": masked,
        "outputs": {"fasta": str(destination)},
        "provenance": _gimme_provenance(function="gimmemotifs.fasta.Fasta.hardmask"),
    }


def gimmemotifs_random_fasta_subset(
    input_fasta: str | Path,
    n_sequences: int,
    *,
    length: int | None = None,
    seed: int | None = None,
    output_fasta: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw random sequences from a FASTA file (``Fasta.get_random``).

    Without ``length`` the function samples ``n_sequences`` whole records without replacement and
    keeps their names. With ``length`` it draws, from randomly chosen records at least that long,
    one random window of ``length`` bases each and names them ``random1``, ``random2`` and so on.
    ``seed`` seeds Python's random generator for the call. The subset can be written to
    ``output_fasta`` (must not exist).
    """

    fa = _fasta(input_fasta)
    count = _integer(n_sequences, "n_sequences", 1, _MAX_SEQUENCES)
    window = None if length is None else _integer(length, "length", 1, 100_000_000)
    if window is None and count > len(fa):
        raise ConfigurationError("n_sequences is larger than the number of records.")
    chosen_seed = _seed(seed)
    destination = None if output_fasta is None else _new_output(output_fasta, "output_fasta")
    with _seeded(chosen_seed):
        subset = _call(fa.get_random, count, length=window)
    if subset is None:
        raise ConfigurationError("Not enough sequences of the requested length.")
    outputs: dict[str, str] = {}
    if destination is not None:
        _call(subset.writefasta, str(destination))
        outputs["fasta"] = str(destination)
    return {
        "status": "completed",
        "operation": "random_fasta_subset",
        "seed": chosen_seed,
        "length": window,
        "n_sequences": len(subset),
        "sequences": {name: seq for name, seq in subset.items()},
        "outputs": outputs,
        "provenance": _gimme_provenance(function="gimmemotifs.fasta.Fasta.get_random"),
    }


def gimmemotifs_sequence_dict(
    source: str | Path | Sequence[str],
    *,
    genome_fasta: str | Path | None = None,
    minsize: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read sequences into a name-keyed dictionary (``gimmemotifs.utils.as_seqdict``).

    ``source`` is a FASTA file (names are the complete header lines), a BED or region
    (``chrom:start-end``) file or list that is cut out of ``genome_fasta``, or a list of DNA
    strings (named ``1``, ``2``, ...). With ``minsize`` every sequence must have at least that many
    bases.
    """

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    limit = None if minsize is None else _integer(minsize, "minsize", 0, 100_000_000)
    genome = None if genome_fasta is None else str(input_path(genome_fasta))
    if isinstance(source, (str, Path)):
        converted: Any = str(input_path(source))
    else:
        if not isinstance(source, Sequence) or not source or len(source) > _MAX_SEQUENCES:
            raise ConfigurationError("source must be a path or a non-empty list of strings.")
        if not all(isinstance(item, str) and item.strip() for item in source):
            raise ConfigurationError("source list items must be non-empty strings.")
        converted = [str(item).strip() for item in source]
    sequences = _call(utils.as_seqdict, converted, genome=genome, minsize=limit)
    return {
        "status": "completed",
        "operation": "sequence_dict",
        "n_sequences": len(sequences),
        "sequences": {str(name): str(seq) for name, seq in sequences.items()},
        "provenance": _gimme_provenance(function="gimmemotifs.utils.as_seqdict"),
    }


def _stats_table(
    stats: Mapping[str, Mapping[str, Mapping[str, float]]],
) -> dict[str, dict[str, dict[str, float]]]:
    if not isinstance(stats, Mapping) or not 1 <= len(stats) <= _MAX_STATS_MOTIFS:
        raise ConfigurationError("stats must map motif names to {background: {metric: value}}.")
    table: dict[str, dict[str, dict[str, float]]] = {}
    for motif, per_background in stats.items():
        if not isinstance(per_background, Mapping) or not per_background:
            raise ConfigurationError(
                f"stats[{motif!r}] must map background names to metric dictionaries."
            )
        table[str(motif)] = {}
        for background, metrics in per_background.items():
            if not isinstance(metrics, Mapping) or not metrics:
                raise ConfigurationError(
                    f"stats[{motif!r}][{background!r}] must be a non-empty metric dictionary."
                )
            table[str(motif)][str(background)] = dict(metrics)
    return table


def gimmemotifs_rank_motif_stats(
    stats: Mapping[str, Mapping[str, Mapping[str, float]]],
    *,
    metrics: Sequence[str] = ("roc_auc", "recall_at_fdr"),
    verbose: bool = True,
) -> dict[str, Any]:
    """Rank motifs by their mean rank over several metrics (``gimmemotifs.stats.rank_motifs``).

    ``stats`` maps motif name to background name to metric name to value, as produced by the
    motif-enrichment statistics. For each metric the motifs' values are averaged over the
    backgrounds and ranked (rank 1 is the lowest value, ties get the average rank); the result is
    the mean of those ranks per motif, so a higher value means a better motif.
    """

    stats_module = require("gimmemotifs.stats", _GIMMEMOTIFS_EXTRA)
    table = _stats_table(stats)
    names = [str(metric) for metric in metrics]
    if not names or isinstance(metrics, str):
        raise ConfigurationError("metrics must be a non-empty list of metric names.")
    backgrounds = list(next(iter(table.values())).keys())
    for motif, per_background in table.items():
        for background in backgrounds:
            if background not in per_background:
                raise ConfigurationError(f"stats[{motif!r}] lacks background {background!r}.")
            for metric in names:
                if metric not in per_background[background]:
                    raise ConfigurationError(
                        f"stats[{motif!r}][{background!r}] lacks metric {metric!r}."
                    )
    ranks = _call(stats_module.rank_motifs, table, metrics=tuple(names))
    ordered = sorted(ranks.items(), key=lambda item: (-float(item[1]), item[0]))
    return {
        "status": "completed",
        "operation": "rank_motif_stats",
        "metrics": names,
        "backgrounds": backgrounds,
        "mean_rank": {str(name): float(value) for name, value in ranks.items()},
        "order": [str(name) for name, _ in ordered],
        "provenance": _gimme_provenance(function="gimmemotifs.stats.rank_motifs"),
    }


def gimmemotifs_write_motif_stats(
    stats: Mapping[str, Mapping[str, Mapping[str, float]]],
    output_template: str,
    *,
    header: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the motif statistics table, one tab-separated file per background.

    Calls ``gimmemotifs.stats.write_stats``.

    ``output_template`` is a path containing one ``{}`` that is replaced by each background name,
    for example ``stats.{}.txt``. A file has a ``Motif`` column and one sorted column per metric.
    As upstream does, the motif name is written without its last ``_`` part (``MA0001_1_run``
    becomes ``MA0001_1``). ``header`` is text written first. Existing files are not overwritten.
    """

    stats_module = require("gimmemotifs.stats", _GIMMEMOTIFS_EXTRA)
    table = _stats_table(stats)
    if not isinstance(output_template, str) or output_template.count("{}") != 1:
        raise ConfigurationError("output_template must be a string with exactly one {}.")
    backgrounds = list(next(iter(table.values())).keys())
    paths = {bg: _new_output(output_template.format(bg), "output file") for bg in backgrounds}
    text = None if header is None else str(header)
    _call(stats_module.write_stats, table, output_template, header=text)
    return {
        "status": "completed",
        "operation": "write_motif_stats",
        "backgrounds": backgrounds,
        "n_motifs": len(table),
        "outputs": {bg: str(path) for bg, path in paths.items()},
        "provenance": _gimme_provenance(function="gimmemotifs.stats.write_stats"),
    }


def _scanner(seed: int | None) -> Any:
    scanner_module = require("gimmemotifs.scanner", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    state = None if seed is None else np.random.RandomState(seed)
    return _call(scanner_module.Scanner, ncpus=1, random_state=state, progress=False)


def _fpr(value: object) -> float:
    number = float(value)  # type: ignore[arg-type]
    if not math.isfinite(number) or not 0.0 < number < 1.0:
        raise ConfigurationError("fpr must be a number strictly between 0 and 1.")
    return number


def _background_thresholds(
    scanner: Any, motifs_path: Path, background: Any, fpr: float
) -> dict[str, Any]:
    """Scan a background FASTA and set the FPR-based per-motif thresholds on a scanner."""

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    motifs = _call(motif_module.read_motifs, str(motifs_path))
    if not 1 <= len(motifs) <= _MAX_MOTIF_FILE_ROWS:
        raise ConfigurationError(f"motif file must hold 1 to {_MAX_MOTIF_FILE_ROWS} motifs.")
    scanner.set_motifs(str(motifs_path))
    scanner.background = background
    scanner.fpr = fpr
    # The upstream set_threshold() would also cache the thresholds under Motif.hash, which fails
    # with xxhash >= 4; the scan of the background sequences itself is the same call it makes.
    scanner._threshold = _call(scanner._threshold_from_seqs, motifs, background)
    return _call(scanner.get_gc_thresholds, background.seqs, motifs=motifs)


def gimmemotifs_motif_thresholds(
    motifs_pfm: str | Path,
    background_fasta: str | Path,
    *,
    fpr: float = 0.01,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Motif score thresholds at a target false-positive rate on background sequences.

    Equivalent of ``Scanner.set_threshold`` and the ``gimme threshold`` command.

    Every motif of the PFM file is scanned over the background FASTA (best hit per sequence, both
    strands); the threshold of a motif is the ``1 - fpr`` quantile of those best scores, so a
    fraction ``fpr`` of background sequences scores above it. A threshold of ``None`` means the
    quantile reaches the motif's maximum score (no hit can pass). gimmemotifs resamples the
    background scores; pass ``seed`` for reproducible values. No GC correction is applied.
    """

    chosen_fpr = _fpr(fpr)
    motifs_path = input_path(motifs_pfm)
    background = _fasta(background_fasta, "background_fasta")
    scanner = _scanner(_seed(seed))
    thresholds = _background_thresholds(scanner, motifs_path, background, chosen_fpr)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    motifs = {m.id: m for m in motif_module.read_motifs(str(motifs_path))}
    rows = {}
    for motif_id, value in thresholds.items():
        motif = motifs[motif_id]
        score = None if value is None else float(value)
        span = float(motif.max_score) - float(motif.min_score)
        fraction = None if score is None or span == 0 else (score - float(motif.min_score)) / span
        rows[str(motif_id)] = {
            "threshold": score,
            "fraction": fraction,
            "max_score": float(motif.max_score),
        }
    return {
        "status": "completed",
        "operation": "motif_thresholds",
        "fpr": chosen_fpr,
        "n_background_sequences": len(background),
        "thresholds": rows,
        "provenance": _gimme_provenance(function="gimmemotifs.scanner.base.Scanner.set_threshold"),
    }


def gimmemotifs_count_motif_matches(
    sequences_fasta: str | Path,
    motifs_pfm: str | Path,
    *,
    threshold: float = 0.95,
    background_fasta: str | Path | None = None,
    fpr: float | None = None,
    nreport: int = 100,
    scan_rc: bool = True,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the motif matches per sequence and in total.

    Calls ``Scanner.count`` and ``Scanner.total_count``.

    A match is a hit of a PFM-file motif above the motif's cutoff, on either strand when
    ``scan_rc``. By default the cutoff is ``threshold``, the fraction (0 to 1) of the distance
    from the minimum to the maximum motif score (0.95, the gimmemotifs default without a genome).
    With ``fpr`` and ``background_fasta`` the cutoff is the score reached by a fraction ``fpr`` of
    the background sequences (see ``gimmemotifs_motif_thresholds``). At most ``nreport`` matches
    are counted per motif and sequence. Returns the counts per sequence and the total per motif.
    """

    scan_threshold = float(threshold)
    if not math.isfinite(scan_threshold) or not 0.0 <= scan_threshold <= 1.0:
        raise ConfigurationError("threshold must be a number in [0, 1].")
    limit = _integer(nreport, "nreport", 1, 1_000_000)
    strand = _flag(scan_rc, "scan_rc")
    motifs_path = input_path(motifs_pfm)
    fa = _fasta(sequences_fasta, "sequences_fasta")
    scanner = _scanner(_seed(seed))
    if (fpr is None) != (background_fasta is None):
        raise ConfigurationError("fpr and background_fasta must be given together.")
    mode = "threshold"
    if fpr is not None and background_fasta is not None:
        mode = "fpr"
        background = _fasta(background_fasta, "background_fasta")
        _background_thresholds(scanner, motifs_path, background, _fpr(fpr))
    else:
        scanner.set_motifs(str(motifs_path))
        _call(scanner.set_threshold, threshold=scan_threshold)
    counted = _call(lambda: list(scanner.count(fa, nreport=limit, scan_rc=strand)))
    totals = _call(scanner.total_count, fa, nreport=limit, scan_rc=strand)
    motif_ids = [str(name) for name in scanner.motif_ids]
    return {
        "status": "completed",
        "operation": "count_motif_matches",
        "cutoff_mode": mode,
        "threshold": scan_threshold if mode == "threshold" else None,
        "fpr": None if fpr is None else float(fpr),
        "motifs": motif_ids,
        "n_sequences": len(counted),
        "counts": {
            str(name): {m: int(c) for m, c in zip(motif_ids, row, strict=True)}
            for name, row in counted
        },
        "total": {m: int(c) for m, c in zip(motif_ids, totals, strict=True)},
        "provenance": _gimme_provenance(function="gimmemotifs.scanner.base.Scanner.count"),
    }


def gimmemotifs_narrowpeak_to_bed(
    input_file: str | Path,
    output_bed: str | Path,
    *,
    size: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a narrowPeak file to a four-column BED file.

    Calls ``gimmemotifs.utils.narrowpeak_to_bed``.

    With ``size=0`` every peak keeps its coordinates. With a positive ``size`` the interval becomes
    ``size`` bases long, starting at ``peak start + summit - size // 2`` (the summit is column 10;
    a summit of -1 means no summit and the peak centre offset ``(end - start) // 2`` is used). The
    BED columns are chromosome, start, end and the signal value (narrowPeak column 7). Comment,
    ``track`` and ``browser`` lines are skipped. The output file must not exist.
    """

    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    source = input_path(input_file)
    width = _integer(size, "size", 0, 1_000_000)
    destination = _new_output(output_bed, "output_bed")
    rows = 0
    for line in source.read_text().splitlines():
        if line.startswith(("#", "track", "browser")) or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            raise ConfigurationError("narrowPeak rows need 10 tab-separated columns.")
        rows += 1
    _call(utils.narrowpeak_to_bed, str(source), str(destination), size=width)
    intervals = [line.split("\t") for line in destination.read_text().splitlines()]
    return {
        "status": "completed",
        "operation": "narrowpeak_to_bed",
        "size": width,
        "n_intervals": len(intervals),
        "intervals": [[c, int(a), int(b), v] for c, a, b, v in intervals],
        "outputs": {"bed": str(destination)},
        "provenance": _gimme_provenance(function="gimmemotifs.utils.narrowpeak_to_bed"),
    }


def gimmemotifs_prepare_denovo_input(
    input_file: str | Path,
    output_dir: str | Path,
    *,
    input_type: str = "fasta",
    genome_fasta: str | Path | None = None,
    size: int | None = None,
    fraction: float | None = None,
    abs_max: int | None = None,
    lsize: int | None = None,
    use_strand: bool = False,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prepare the input files of a de novo motif run.

    Calls ``motif.denovo.prepare_denovo_input_fa``, ``_bed`` and ``_narrowpeak``.

    ``input_type`` is ``fasta`` (the records are split into ``prediction.fa`` and
    ``validation.fa``; ``localization.fa`` is a copy of the validation set), ``bed`` (regions are
    resized to ``size`` bases around their midpoint into ``input.bed``, split into
    ``prediction.bed`` and ``validation.bed`` and cut out of ``genome_fasta``; the localization
    FASTA is extended to ``lsize`` bases) or ``narrowpeak`` (first converted to a BED of ``size``
    bases around the peak summit, then as ``bed``). The prediction set takes ``fraction`` of the
    regions, at most ``abs_max``. Unset values use the gimmemotifs defaults (size 200, fraction 0.2,
    abs_max 1000, lsize 500). ``seed`` seeds the random split. ``output_dir`` must be empty or new.
    """

    denovo = require("gimmemotifs.motif.denovo", _GIMMEMOTIFS_EXTRA)
    config = require("gimmemotifs.config", _GIMMEMOTIFS_EXTRA).MotifConfig().get_default_params()
    chosen = _choose(input_type, _INPUT_TYPES, "input_type")
    source = input_path(input_file)
    window = int(config["size"]) if size is None else _integer(size, "size", 1, 1_000_000)
    portion = float(config["fraction"]) if fraction is None else float(fraction)
    if not math.isfinite(portion) or not 0.0 < portion <= 1.0:
        raise ConfigurationError("fraction must be a number in (0, 1].")
    maximum = (
        int(config["abs_max"]) if abs_max is None else _integer(abs_max, "abs_max", 1, 10_000_000)
    )
    wide = int(config["lsize"]) if lsize is None else _integer(lsize, "lsize", 1, 10_000_000)
    params = {
        "size": window,
        "fraction": portion,
        "abs_max": maximum,
        "lsize": wide,
        "use_strand": _flag(use_strand, "use_strand"),
    }
    if chosen != "fasta":
        if genome_fasta is None:
            raise ConfigurationError(f"input_type={chosen!r} needs genome_fasta.")
        params["genome"] = str(input_path(genome_fasta))
    destination = Path(output_dir).expanduser()
    if destination.exists() and any(destination.iterdir()):
        raise ConfigurationError(f"output_dir is not empty: {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    function = {
        "fasta": denovo.prepare_denovo_input_fa,
        "bed": denovo.prepare_denovo_input_bed,
        "narrowpeak": denovo.prepare_denovo_input_narrowpeak,
    }[chosen]
    chosen_seed = _seed(seed)
    with _seeded(chosen_seed):
        _call(function, str(source), params, str(destination))
    files = sorted(path.name for path in destination.iterdir())
    counts: dict[str, int] = {}
    for name in files:
        if name.endswith((".fa", ".bed")):
            text = (destination / name).read_text()
            counts[name] = (
                text.count(">")
                if name.endswith(".fa")
                else len([ln for ln in text.splitlines() if ln])
            )
    return {
        "status": "completed",
        "operation": "prepare_denovo_input",
        "input_type": chosen,
        "parameters": {key: value for key, value in params.items() if key != "genome"},
        "seed": chosen_seed,
        "records_per_file": counts,
        "outputs": {name: str(destination / name) for name in files},
        "provenance": _gimme_provenance(
            function=f"gimmemotifs.motif.denovo.prepare_denovo_input_{chosen}"
        ),
    }


def _db(name: str) -> Any:
    db_module = require("gimmemotifs.db", _GIMMEMOTIFS_EXTRA)
    return db_module, _call(db_module.MotifDb.create, name)


def gimmemotifs_download_motif_database(
    database: str,
    output_dir: str | Path,
    *,
    jaspar_version: str = "2022",
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a public motif database and write it in gimmemotifs PFM format (network).

    Calls ``MotifDb.download``.

    ``database`` is ``jaspar`` (all JASPAR CORE groups of ``jaspar_version``), ``homer``,
    ``hocomoco`` (v11 human and mouse), ``encode``, ``swissregulon`` or ``image``. The upstream
    downloads of Factorbook, CIS-BP and RSAT are not offered: their source URLs answer 404 (checked
    2026-10-20). For each PFM file the matching
    ``.motif2factors.txt`` annotation is written next to it. The output directory is created if
    needed; files are overwritten as upstream does. Returns the files that were written, with the
    database date.
    """

    chosen = _choose(database, _DATABASES, "database")
    version = str(jaspar_version)
    if not re.fullmatch(r"\d{4}", version):
        raise ConfigurationError("jaspar_version must be a four-digit year such as 2022.")
    _, instance = _db(chosen)
    destination = Path(output_dir).expanduser()
    destination.mkdir(parents=True, exist_ok=True)
    before = {p.name: p.stat().st_mtime_ns for p in destination.iterdir()}
    if chosen == "jaspar":
        _call(instance.download, version=version, outdir=str(destination))
    else:
        _call(instance.download, outdir=str(destination))
    written = sorted(
        p.name
        for p in destination.iterdir()
        if p.name not in before or p.stat().st_mtime_ns != before[p.name]
    )
    return {
        "status": "completed",
        "operation": "download_motif_database",
        "database": chosen,
        "date": str(instance.date),
        "outputs": {name: str(destination / name) for name in written},
        "provenance": _gimme_provenance(
            function=f"gimmemotifs.db.{type(instance).__name__}.download"
        ),
    }


def gimmemotifs_annotate_motif_factors(
    database: str,
    motifs_pfm: str | Path,
    *,
    hocomoco_group: str = "HUMAN",
    output_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate the motifs of a PFM file with transcription factors from a motif database (network).

    Calls ``MotifDb.annotate_factors``.

    ``database`` is ``jaspar`` (looks up every ``MA`` matrix at JASPAR), ``homer``,
    ``hocomoco`` (``hocomoco_group`` ``HUMAN`` or ``MOUSE``), ``encode``,
    ``swissregulon`` or ``image``. The motif IDs must follow the naming of that database's own PFM
    file (the file written by ``gimmemotifs_download_motif_database``). Returns,
    per motif, the list of ``[factor, evidence, curated]``; with ``output_dir`` the table is also
    written as ``<pfm name>.motif2factors.txt`` (columns Motif, Factor, Evidence, Curated).
    """

    chosen = _choose(database, _DATABASES, "database")
    group = _choose(hocomoco_group, _HOCOMOCO_GROUPS, "hocomoco_group")
    path = input_path(motifs_pfm)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    motifs = _call(motif_module.read_motifs, str(path))
    if not 1 <= len(motifs) <= _MAX_MOTIF_FILE_ROWS:
        raise ConfigurationError(f"motifs_pfm must hold 1 to {_MAX_MOTIF_FILE_ROWS} motifs.")
    _, instance = _db(chosen)
    if chosen == "hocomoco":
        annotation = _call(instance.annotate_factors, motifs, instance.ANNO_URL.format(group))
    else:
        annotation = _call(instance.annotate_factors, motifs)
    table: dict[str, list[list[str]]] = {}
    for motif_id, rows in annotation.items():
        table[str(motif_id)] = [[str(v) for v in list(row)] for row in rows]
    outputs: dict[str, str] = {}
    if output_dir is not None:
        directory = Path(output_dir).expanduser()
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / (path.stem + ".motif2factors.txt")
        if target.exists():
            raise ConfigurationError(f"output file already exists: {target}")
        _call(instance.create_annotation, str(directory / path.name), annotation)
        outputs["motif2factors"] = str(target)
    return {
        "status": "completed",
        "operation": "annotate_motif_factors",
        "database": chosen,
        "n_motifs": len(motifs),
        "n_annotated": len(table),
        "annotation": table,
        "outputs": outputs,
        "provenance": _gimme_provenance(
            function=f"gimmemotifs.db.{type(instance).__name__}.annotate_factors"
        ),
    }


__all__ = [
    "gimmemotifs_annotate_motif_factors",
    "gimmemotifs_count_motif_matches",
    "gimmemotifs_download_motif_database",
    "gimmemotifs_hardmask_fasta",
    "gimmemotifs_motif_thresholds",
    "gimmemotifs_narrowpeak_to_bed",
    "gimmemotifs_prepare_denovo_input",
    "gimmemotifs_random_fasta_subset",
    "gimmemotifs_rank_motif_stats",
    "gimmemotifs_sequence_dict",
    "gimmemotifs_write_motif_stats",
]

_register_public_verbose(globals(), tuple(__all__))
