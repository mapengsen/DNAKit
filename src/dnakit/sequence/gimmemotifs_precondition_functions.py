"""GimmeMotifs wrappers for operations that need local inputs (Fig. 2 coverage gap, batch F).

Every public function takes explicit local files: FASTA sequences, PFM motif files, score or count
tables, a gene annotation GTF, an ortholog SQLite database or an OrthoFinder result folder. None of
them looks up a genomepy genome, downloads a database or queries an online service. The
command-line programs used here (gffread, orthofinder) are resolved on PATH; a missing program
raises BackendUnavailableError with an install hint. Upstream gimmemotifs functions are called
with their own default parameters unless a parameter below says otherwise.
"""

from __future__ import annotations

import importlib.metadata
import re
import shutil
import sqlite3
import subprocess
from collections.abc import Mapping, Sequence
from contextlib import closing
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import artifacts, require, workspace
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "gimmemotifs"
_OFFLINE_STRATEGY = "strict"
_BACKGROUND_TYPES = ("random", "genomic", "gc", "promoter", "custom")
_VALIDATION_TYPES = ("random", "genomic", "gc", "custom")
_LOG_TAIL = 2000


def _require(module: str) -> Any:
    return require(module, _EXTRA)


def _version() -> str | None:
    try:
        return importlib.metadata.version("gimmemotifs")
    except importlib.metadata.PackageNotFoundError:
        return None


def _file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} does not exist: {path}")
    return path


def _fresh(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _names(value: object, name: str, *, minimum: int = 1) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of names.")
    items = list(value)
    if len(items) < minimum:
        raise ConfigurationError(
            f"{name} needs at least {minimum} entr{'y' if minimum == 1 else 'ies'}."
        )
    for item in items:
        if not isinstance(item, str) or not item.strip() or "'" in item:
            raise ConfigurationError(f"{name} must hold non-empty strings without quotes.")
    return items


def _number(value: object, name: str, *, low: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    if low is not None and value < low:
        raise ConfigurationError(f"{name} must be at least {low}.")
    return float(value)


def _positive(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _vector(value: object, length: int, name: str) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != length:
        raise ConfigurationError(f"{name} must hold {length} values, one per motif.")
    return [_number(item, name) for item in value]


def _matrix(value: object, rows: int, columns: int, name: str) -> list[list[float]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != rows:
        raise ConfigurationError(f"{name} must hold {rows} rows, one per motif.")
    table: list[list[float]] = []
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != columns:
            raise ConfigurationError(f"each row of {name} must hold {columns} values.")
        table.append([_number(item, name) for item in row])
    return table


def _motifs(pfm_file: Path) -> dict[str, Any]:
    motif_module = _require("gimmemotifs.motif")
    motifs = motif_module.read_motifs(str(pfm_file), fmt="pfm", as_dict=True)
    if not motifs:
        raise ConfigurationError(f"No motifs were read from {pfm_file}.")
    return dict(motifs)


def _program(name: str, hint: str) -> str:
    found = shutil.which(name)
    if found is None:
        raise BackendUnavailableError(f"The {name} program was not found on PATH.", hint=hint)
    return str(Path(found).resolve())


def _written(path: Path, what: str) -> Path:
    if not path.exists():
        raise BackendExecutionError(f"gimmemotifs did not write {what}: {path}")
    return path


def _copy(source: Path, target: Path) -> None:
    shutil.copyfile(source, target)


def _tail(text: str) -> str:
    return text[-_LOG_TAIL:]


def _result(operation: str, outputs: Mapping[str, object], **extra: Any) -> dict[str, Any]:
    return {
        "status": "complete",
        "operation": operation,
        "backend": "gimmemotifs",
        "version": _version(),
        "outputs": {key: str(value) for key, value in outputs.items()},
        **extra,
    }


def _backgrounds(value: object) -> dict[str, str]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError("background_fastas must map background names to FASTA paths.")
    return {str(name): str(_file(path, f"background {name}")) for name, path in value.items()}


def gimmemotifs_diff_plot(
    motif_ids: Sequence[str],
    pfm_file: str,
    sample_names: Sequence[str],
    frequencies: Sequence[Sequence[float]],
    counts: Sequence[Sequence[float]],
    background_frequencies: Sequence[float],
    background_counts: Sequence[float],
    output_path: str,
    *,
    min_difference: float = 0,
    min_enrichment: float = 3,
    min_frequency: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot motif frequencies of input sets against a background (``gimme diff`` plot).

    ``frequencies`` and ``counts`` hold one row per motif in ``motif_ids`` and one column per entry
    of ``sample_names``; ``background_frequencies`` and ``background_counts`` hold one value per
    motif. The PNG is written to ``output_path`` and the tables are returned as given.
    Upstream keeps only motifs whose enrichment differs between sample sets, so at least two
    sample sets are needed for a plot; with one set ``plot_written`` is False.
    """
    motifs = _motifs(_file(pfm_file, "pfm_file"))
    ids = _names(motif_ids, "motif_ids")
    names = _names(sample_names, "sample_names")
    missing = [motif for motif in ids if motif not in motifs]
    if missing:
        raise ConfigurationError(f"Motifs are not in {pfm_file}: {missing[:5]}")
    freq = _matrix(frequencies, len(ids), len(names), "frequencies")
    count = _matrix(counts, len(ids), len(names), "counts")
    bg_freq = _vector(background_frequencies, len(ids), "background_frequencies")
    bg_count = _vector(background_counts, len(ids), "background_counts")
    out = _fresh(output_path, "output_path")
    np = _require("numpy")
    plot = _require("gimmemotifs.plot")
    plot.diff_plot(
        ids,
        motifs,
        names,
        np.array(freq, dtype=float),
        np.array(count, dtype=float),
        np.array(bg_freq, dtype=float),
        np.array(bg_count, dtype=float),
        str(out),
        mindiff=_number(min_difference, "min_difference"),
        minenr=_number(min_enrichment, "min_enrichment"),
        minfreq=_number(min_frequency, "min_frequency", low=0),
    )
    # Upstream draws nothing when no motif passes the enrichment and frequency thresholds.
    drawn = out.exists()
    return _result(
        "diff_plot",
        {"plot": out} if drawn else {},
        plot_written=drawn,
        motif_ids=ids,
        sample_names=names,
        frequencies=freq,
        counts=count,
        background_frequencies=bg_freq,
        background_counts=bg_count,
    )


def gimmemotifs_diff_motifs(
    foreground_fastas: Sequence[str],
    background_fasta: str,
    pfm_file: str,
    output_path: str,
    *,
    cutoff: float = 0.9,
    min_enrichment: float = 3,
    min_frequency: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare motif frequencies of FASTA sets with a background FASTA (``gimme diff``).

    The steps follow the ``gimme diff`` command: each motif is scanned with ``cutoff`` and the
    number of sequences with a hit is counted; frequencies are ``(count + 0.01) / n`` for the set
    of ``n`` sequences. Set names are the FASTA base names. The plot is written to
    ``output_path``; the frequency and count tables are returned.
    """
    fg_files = [
        _file(item, "foreground_fastas item")
        for item in _names(foreground_fastas, "foreground_fastas")
    ]
    names = [path.stem for path in fg_files]
    if len(set(names)) != len(names):
        raise ConfigurationError("Foreground FASTA base names must be unique.")
    background = _file(background_fasta, "background_fasta")
    pfm = _file(pfm_file, "pfm_file")
    motifs = _motifs(pfm)
    motif_ids = list(motifs)
    threshold = _number(cutoff, "cutoff", low=0)
    fasta_module = _require("gimmemotifs.fasta")
    scanner_module = _require("gimmemotifs.scanner")
    scanner = scanner_module.Scanner()
    scanner.set_motifs(str(pfm))
    scanner.set_threshold(threshold=threshold)
    background_size = float(len(fasta_module.Fasta(str(background)).seqs))
    background_counts = [float(c) for c in scanner.total_count(str(background), nreport=1)]
    background_frequencies = [(c + 0.01) / background_size for c in background_counts]
    frequency_columns: list[list[float]] = []
    count_columns: list[list[float]] = []
    for fasta in fg_files:
        counts = [float(c) for c in scanner.total_count(str(fasta), nreport=1)]
        size = float(len(fasta_module.Fasta(str(fasta)).seqs))
        count_columns.append(counts)
        frequency_columns.append([(c + 0.01) / size for c in counts])
    frequencies = [list(row) for row in zip(*frequency_columns, strict=True)]
    count_table = [list(row) for row in zip(*count_columns, strict=True)]
    plotted = gimmemotifs_diff_plot(
        motif_ids,
        str(pfm),
        names,
        frequencies,
        count_table,
        background_frequencies,
        background_counts,
        output_path,
        min_enrichment=min_enrichment,
        min_frequency=min_frequency,
        verbose=False,
    )
    return {
        **plotted,
        "operation": "diff_motifs",
        "cutoff": threshold,
        "background_size": int(background_size),
        "foreground_sizes": {
            name: _size(path, fasta_module) for name, path in zip(names, fg_files, strict=True)
        },
    }


def _size(path: Path, fasta_module: Any) -> int:
    return len(fasta_module.Fasta(str(path)).seqs)


def gimmemotifs_select_nonredundant_motifs(
    roc_report: str,
    pfm_file: str,
    foreground_table: str,
    background_table: str,
    *,
    tolerance: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select non-redundant motifs from a ROC report and two motif score tables.

    ``roc_report`` is a tab-separated table with motif ids in the first column and the columns
    "ROC AUC", "PR AUC", "Enr. at 1% FPR" and "Recall at 10% FDR". ``foreground_table`` and
    ``background_table`` are tab-separated score tables with one column per motif id. Motifs are
    ranked on the report, redundant ones (seqcor at least 0.7) are removed, and a recursive
    feature elimination keeps the motifs that improve the average precision.
    """
    comparison = _require("gimmemotifs.comparison")
    selected = comparison.select_nonredundant_motifs(
        str(_file(roc_report, "roc_report")),
        str(_file(pfm_file, "pfm_file")),
        str(_file(foreground_table, "foreground_table")),
        str(_file(background_table, "background_table")),
        tolerance=_number(tolerance, "tolerance", low=0),
    )
    return _result(
        "select_nonredundant_motifs",
        {},
        selected_motifs=[str(motif) for motif in selected],
    )


def gimmemotifs_denovo_create_background(
    background_type: str,
    fasta_file: str,
    output_path: str,
    *,
    genome_fasta: str | None = None,
    annotation_bed: str | None = None,
    size: int = 200,
    nr_times: int = 10,
    custom_background: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create one background FASTA for de novo motif analysis (``create_background``).

    ``background_type`` is ``random`` (first-order Markov model of ``fasta_file``), ``genomic``
    (random regions of ``genome_fasta``), ``gc`` (regions of ``genome_fasta`` matched to the GC
    content of ``fasta_file``), ``promoter`` (random promoters from ``annotation_bed``, a BED of
    genes, on ``genome_fasta``) or ``custom`` (``custom_background`` copied as it is). ``nr_times``
    gives the number of background sequences per foreground sequence. The genome is a FASTA path,
    never a genomepy name.
    """
    if background_type not in _BACKGROUND_TYPES:
        raise ConfigurationError(f"background_type must be one of {list(_BACKGROUND_TYPES)}.")
    source = _file(fasta_file, "fasta_file")
    target = _fresh(output_path, "output_path")
    width = _positive(size, "size")
    times = _positive(nr_times, "nr_times")
    genome = (
        _file(genome_fasta, "genome_fasta")
        if background_type in ("genomic", "gc", "promoter")
        else None
    )
    custom = _file(custom_background, "custom_background") if background_type == "custom" else None
    if background_type == "promoter":
        genes = _file(annotation_bed, "annotation_bed")
        background = _require("gimmemotifs.background")
        foreground = _require("gimmemotifs.fasta").Fasta(str(source))
        sequences = background.PromoterFasta(
            str(genes), str(genome), width, times * len(foreground)
        )
        sequences.writefasta(str(target))
        count = len(sequences)
    else:
        denovo = _require("gimmemotifs.motif.denovo")
        count = denovo.create_background(
            background_type,
            str(source),
            str(target),
            genome=None if genome is None else str(genome),
            size=width,
            nr_times=times,
            custom_background=None if custom is None else str(custom),
        )
    return _result(
        "denovo_create_background",
        {"background": _written(target, "the background FASTA")},
        background_type=background_type,
        sequence_count=int(count),
    )


def gimmemotifs_denovo_create_backgrounds(
    output_dir: str,
    *,
    prediction_fasta: str,
    validation_fasta: str,
    background_types: Sequence[str] = ("random",),
    genome_fasta: str | None = None,
    size: int = 200,
    custom_background: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the prediction and validation backgrounds of a de novo run (``create_backgrounds``).

    ``prediction_fasta`` and ``validation_fasta`` are copied into ``output_dir`` as the upstream
    layout expects. The prediction background uses the first type (or ``gc`` when listed); one
    background file per type is written for validation. Types are ``random``, ``genomic``, ``gc``
    and ``custom``; ``promoter`` needs an annotation and is available through
    gimmemotifs_denovo_create_background.
    """
    types = _names(background_types, "background_types")
    unknown = [item for item in types if item not in _VALIDATION_TYPES]
    if unknown:
        raise ConfigurationError(
            f"background_types must be among {list(_VALIDATION_TYPES)}; got {unknown}."
        )
    prediction = _file(prediction_fasta, "prediction_fasta")
    validation = _file(validation_fasta, "validation_fasta")
    needs_genome = any(item in ("genomic", "gc") for item in types)
    genome = _file(genome_fasta, "genome_fasta") if needs_genome else None
    custom = _file(custom_background, "custom_background") if "custom" in types else None
    width = _positive(size, "size")
    folder = workspace(output_dir)
    _copy(prediction, folder / "prediction.fa")
    _copy(validation, folder / "validation.fa")
    denovo = _require("gimmemotifs.motif.denovo")
    info = denovo.create_backgrounds(
        str(folder),
        background=types,
        genome=None if genome is None else str(genome),
        size=width,
        custom_background=None if custom is None else str(custom),
    )
    outputs: dict[str, Path] = {
        "prediction_background": _written(folder / "prediction.bg.fa", "the prediction background")
    }
    for background_name, path in info.items():
        outputs[f"validation_{background_name}"] = _written(
            Path(path), f"the {background_name} background"
        )
    return _result(
        "denovo_create_backgrounds",
        outputs,
        background_types=types,
        background_files={name: str(path) for name, path in info.items()},
    )


def gimmemotifs_filter_significant_motifs(
    motif_file: str,
    output_path: str,
    *,
    background_type: str,
    stats: Mapping[str, Mapping[str, Mapping[str, float]]],
    metrics: Mapping[str, float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the de novo motifs whose statistics pass minimum values (``filter_significant_motifs``).

    ``stats`` maps each motif id to background type to statistics (``max_enrichment``,
    ``roc_auc``, ``enr_at_fpr``, ...). ``metrics`` maps statistic names to minimum values; the
    upstream defaults (max_enrichment 3, roc_auc 0.55, enr_at_fpr 0.55) apply when it is None.
    The significant motifs are written as PFM text to ``output_path``.
    """
    motifs = list(_motifs(_file(motif_file, "motif_file")).values())
    target = _fresh(output_path, "output_path")
    if not isinstance(background_type, str) or not background_type:
        raise ConfigurationError("background_type must be a non-empty string.")
    if not isinstance(stats, Mapping):
        raise ConfigurationError("stats must map motif ids to background types to statistics.")
    keyed: dict[str, dict[str, dict[str, float]]] = {}
    for motif in motifs:
        per_background = stats.get(motif.id)
        if per_background is not None:
            key = f"{motif.id}_{motif.to_consensus()}"
            keyed[key] = {str(name): dict(values) for name, values in per_background.items()}
    rule = None
    if metrics is not None:
        if not isinstance(metrics, Mapping) or not metrics:
            raise ConfigurationError("metrics must map statistic names to minimum values.")
        rule = tuple((str(name), _number(value, name)) for name, value in metrics.items())
    denovo = _require("gimmemotifs.motif.denovo")
    significant = denovo.filter_significant_motifs(
        str(target),
        SimpleNamespace(motifs=motifs, stats=keyed),
        background_type,
        metrics=rule,
    )
    return _result(
        "filter_significant_motifs",
        {"motifs": _written(target, "the significant motif file")},
        significant_motifs=[motif.id for motif in significant],
        background_type=background_type,
    )


def gimmemotifs_best_motif_in_cluster(
    single_motif_file: str,
    cluster_motif_file: str,
    clusters: Sequence[Mapping[str, Any]],
    foreground_fasta: str,
    background_fastas: Mapping[str, str],
    genome_fasta: str,
    *,
    metrics: Sequence[str] = ("roc_auc", "recall_at_fdr"),
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the best motif of each de novo cluster (``best_motif_in_cluster``).

    Each cluster is a mapping with ``representative`` (the id of the cluster motif, read from
    ``cluster_motif_file``) and ``members`` (ids of the single motifs in ``single_motif_file``).
    Statistics for every motif of the clusters are computed against each background FASTA in
    ``background_fastas`` (name to path) before the upstream ranking runs, so the upstream step
    never computes statistics itself. Those statistics use the z-score and GC options off (see
    ``_roc_values``); the ranking metrics (ROC AUC, recall at FDR) are rank statistics and do not
    change under a per-motif z-score. Best motif ids are returned from best to worst.
    """
    singles = _motifs(_file(single_motif_file, "single_motif_file"))
    clustered = _motifs(_file(cluster_motif_file, "cluster_motif_file"))
    known = {**singles, **clustered}
    if isinstance(clusters, (str, bytes)) or not isinstance(clusters, Sequence) or not clusters:
        raise ConfigurationError("clusters must be a non-empty list of cluster mappings.")
    cluster_objects: list[tuple[Any, list[Any]]] = []
    for index, cluster in enumerate(clusters):
        if not isinstance(cluster, Mapping):
            raise ConfigurationError(f"clusters[{index}] must be a mapping.")
        representative = cluster.get("representative")
        members = _names(cluster.get("members"), f"clusters[{index}].members")
        for name in [representative, *members]:
            if not isinstance(name, str) or name not in known:
                raise ConfigurationError(
                    f"clusters[{index}] names a motif that is not in the files: {name!r}"
                )
        cluster_objects.append((known[representative], [known[name] for name in members]))
    foreground = _file(foreground_fasta, "foreground_fasta")
    backgrounds = _backgrounds(background_fastas)
    genome = _file(genome_fasta, "genome_fasta")
    metric_names = tuple(_names(metrics, "metrics"))
    involved = list(
        {id(motif): motif for pair in cluster_objects for motif in [pair[0], *pair[1]]}.values()
    )
    stats_module = _require("gimmemotifs.stats")
    stats: dict[str, dict[str, Any]] = {}
    for background_name, background_path in backgrounds.items():
        computed = stats_module.calc_stats(
            fg_file=str(foreground),
            bg_file=background_path,
            motifs=involved,
            genome=str(genome),
            zscore=False,
            gc=False,
        )
        for motif_key, values in computed.items():
            stats.setdefault(motif_key, {})[background_name] = values
    denovo = _require("gimmemotifs.motif.denovo")
    best = denovo.best_motif_in_cluster(
        str(_file(single_motif_file, "single_motif_file")),
        str(_file(cluster_motif_file, "cluster_motif_file")),
        cluster_objects,
        str(foreground),
        backgrounds,
        str(genome),
        stats=stats,
        metrics=metric_names,
    )
    return _result("best_motif_in_cluster", {}, best_motifs=[motif.id for motif in best])


def gimmemotifs_denovo_motif_report(
    input_file: str,
    pfm_file: str,
    foreground_fasta: str,
    background_fastas: Mapping[str, str],
    location_fasta: str,
    output_dir: str,
    *,
    genome_fasta: str,
    cutoff_fpr: float = 0.9,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the text and HTML report of de novo motifs (``create_denovo_motif_report``).

    The steps are those of the upstream report: ROC plots, closest database match, statistics,
    localization histograms, then the text and HTML reports. ``input_file`` names the de novo input
    in the report; ``background_fastas`` maps background names to FASTA paths; ``location_fasta``
    holds the sequences of the localization histograms (centred on the peaks). The ROC curves and
    statistics use the z-score and GC options off (see ``_roc_values``); ``cutoff_fpr`` is the
    upstream default 0.9.
    """
    source = _file(input_file, "input_file")
    pfm = _file(pfm_file, "pfm_file")
    foreground = _file(foreground_fasta, "foreground_fasta")
    backgrounds = _backgrounds(background_fastas)
    location = _file(location_fasta, "location_fasta")
    genome = _file(genome_fasta, "genome_fasta")
    cutoff = _number(cutoff_fpr, "cutoff_fpr", low=0)
    folder = workspace(output_dir)
    images = folder / "images"
    images.mkdir()
    report = _require("gimmemotifs.report")
    motifs = _require("gimmemotifs.motif").read_motifs(str(pfm), fmt="pwm")
    for background_name, background_path in backgrounds.items():
        for motif in motifs:
            x_values, y_values = _roc_values(motif, str(foreground), background_path, str(genome))
            report.roc_plot(
                str(images / f"{motif.id}_roc.{background_name}.png"), x_values, y_values
            )
    closest_match = _require("gimmemotifs.comparison").MotifComparer().get_closest_match(motifs)
    stats_module = _require("gimmemotifs.stats")
    stats: dict[str, dict[str, Any]] = {}
    for background_name, background_path in backgrounds.items():
        computed = stats_module.calc_stats(
            fg_file=str(foreground),
            bg_file=background_path,
            motifs=motifs,
            zscore=False,
            gc=False,
        )
        for motif_key, values in computed.items():
            stats.setdefault(motif_key, {})[background_name] = values
    stats = report.add_star(stats)
    lsize = _require("numpy").median(
        [len(seq) for seq in _require("gimmemotifs.fasta").Fasta(str(location)).seqs]
    )
    for motif in motifs:
        report.motif_localization(
            str(location),
            motif,
            lsize,
            str(images / f"{motif.id}_histogram.svg"),
            cutoff=cutoff,
        )
    report._create_text_report(str(source), motifs, closest_match, stats, str(folder))
    report._create_graphical_report(
        str(source), str(pfm), backgrounds, closest_match, str(folder), stats
    )
    files = artifacts(folder)
    if not any(name.endswith(".html") for name in files):
        raise BackendExecutionError("The de novo motif report wrote no HTML report.")
    return _result("denovo_motif_report", dict(sorted(files.items())), report_folder=str(folder))


def _roc_values(
    motif: Any, foreground: str, background: str, genome: str
) -> tuple[list[float], list[float]]:
    """ROC curve of one motif (the upstream ``get_roc_values`` step).

    Upstream calls calc_stats with the default z-score and GC options, which build a genome
    background and hash each motif with xxhash; xxhash 3 and later reject str input, so that path
    raises TypeError with gimmemotifs 0.18.4. A ROC curve is unchanged by the per-motif z-score,
    so the curve is computed with both options off.
    """
    stats = _require("gimmemotifs.stats").calc_stats(
        fg_file=foreground,
        bg_file=background,
        motifs=[motif],
        genome=genome,
        stats=["roc_values"],
        ncpus=1,
        zscore=False,
        gc=False,
    )
    x_values, y_values = next(iter(stats.values()))["roc_values"]
    return list(x_values), list(y_values)


def gimmemotifs_roc_plots(
    pfm_file: str,
    foreground_fasta: str,
    background_fastas: Mapping[str, str],
    output_dir: str,
    *,
    genome_fasta: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw one ROC curve per motif and background (``create_roc_plots``).

    The PNG files are written to ``output_dir/images`` as ``<motif>_roc.<background>.png``, the
    names used upstream. The curves come from scanning ``foreground_fasta`` and each background
    FASTA; the z-score and GC options are off (see ``_roc_values``).
    """
    motifs = _motifs(_file(pfm_file, "pfm_file"))
    foreground = _file(foreground_fasta, "foreground_fasta")
    backgrounds = _backgrounds(background_fastas)
    genome = _file(genome_fasta, "genome_fasta")
    folder = workspace(output_dir)
    images = folder / "images"
    images.mkdir()
    report = _require("gimmemotifs.report")
    for background_name, background_path in backgrounds.items():
        for motif_id, motif in motifs.items():
            x_values, y_values = _roc_values(motif, str(foreground), background_path, str(genome))
            report.roc_plot(
                str(images / f"{motif_id}_roc.{background_name}.png"), x_values, y_values
            )
    written = artifacts(images)
    expected = len(motifs) * len(backgrounds)
    if len(written) != expected:
        raise BackendExecutionError(f"Expected {expected} ROC plots, found {len(written)}.")
    return _result(
        "roc_plots",
        {name: path for name, path in sorted(written.items())},
        motif_count=len(motifs),
        background_count=len(backgrounds),
    )


def gimmemotifs_roc_html_report(
    output_dir: str,
    input_table: str,
    pfm_file: str,
    output_name: str,
    *,
    threshold: float = 0.01,
    use_motifs: Sequence[str] | None = None,
    link_matches: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the HTML report of a motif enrichment table (``roc_html_report``).

    ``input_table`` is the tab-separated enrichment table of a finished run (motif ids in the
    first column; the P-value, match, ROC and recall columns of gimmemotifs). Motifs with a
    corrected P-value above ``threshold`` are left out. Logos are written to ``output_dir/logos``.
    """
    if (
        not isinstance(output_name, str)
        or not output_name
        or "/" in output_name
        or output_name.startswith(".")
    ):
        raise ConfigurationError("output_name must be a plain file name, for example report.html.")
    table = _file(input_table, "input_table")
    pfm = _file(pfm_file, "pfm_file")
    cutoff = _number(threshold, "threshold", low=0)
    if cutoff > 1:
        raise ConfigurationError("threshold must be at most 1.")
    chosen = None if use_motifs is None else _names(use_motifs, "use_motifs")
    folder = workspace(output_dir)
    report = _require("gimmemotifs.report")
    report.roc_html_report(
        str(folder),
        str(table),
        str(pfm),
        output_name,
        threshold=cutoff,
        use_motifs=chosen,
        link_matches=bool(link_matches),
    )
    html = _written(folder / output_name, "the ROC HTML report")
    return _result("roc_html_report", {"report": html, **artifacts(folder)}, threshold=cutoff)


def gimmemotifs_motif_tree_plot(
    output_path: str,
    newick_tree: str,
    leaf_values: Mapping[str, Sequence[float]],
    *,
    circle: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    dpi: int = 300,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a motif clustering tree coloured by values per motif (``motif_tree_plot``).

    ``newick_tree`` is the Newick text of the tree whose leaves are motif ids. ``leaf_values``
    maps each motif id to its values (one per row of the data, for example per cell type); the
    colour of a branch is the mean over the leaves below it. Drawing uses ete3, which must be
    installed; otherwise BackendUnavailableError is raised before anything is written.
    """
    _require("ete3")
    target = _fresh(output_path, "output_path")
    if not isinstance(newick_tree, str) or not newick_tree.strip().endswith(";"):
        raise ConfigurationError("newick_tree must be Newick text ending with ';'.")
    if not isinstance(leaf_values, Mapping) or not leaf_values:
        raise ConfigurationError("leaf_values must map motif ids to values.")
    columns = {str(name): _vector_of_any(values, str(name)) for name, values in leaf_values.items()}
    lengths = {len(values) for values in columns.values()}
    if len(lengths) != 1:
        raise ConfigurationError("every motif in leaf_values needs the same number of values.")
    frame = _require("pandas").DataFrame(columns)
    plot = _require("gimmemotifs.plot")
    plot.motif_tree_plot(
        str(target),
        newick_tree,
        frame,
        circle=bool(circle),
        vmin=None if vmin is None else _number(vmin, "vmin"),
        vmax=None if vmax is None else _number(vmax, "vmax"),
        dpi=_positive(dpi, "dpi"),
    )
    return _result(
        "motif_tree_plot",
        {"plot": _written(target, "the motif tree plot")},
        leaf_count=len(columns),
    )


def _vector_of_any(values: object, name: str) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"values for {name} must be a non-empty list of numbers.")
    return [_number(item, f"values for {name}") for item in values]


def _run_orthofinder(program: str, peptides: Path, threads: int) -> Path:
    completed = subprocess.run(
        [program, "-f", str(peptides), "-t", str(threads)],
        capture_output=True,
        text=True,
        check=False,
    )
    match = re.search(r"Results directory:\s*(\S+)", completed.stdout)
    if completed.returncode != 0 or match is None:
        raise BackendExecutionError(
            "orthofinder did not report a results directory.",
            context={"stderr_tail": _tail(completed.stderr)},
        )
    return Path(match.group(1)).resolve()


def _strict_only(strategy: object) -> str:
    if strategy != _OFFLINE_STRATEGY:
        raise ConfigurationError(
            "strategy must be 'strict': the medium and lenient strategies query mygene.info "
            "over the network, which these offline wrappers do not do."
        )
    return _OFFLINE_STRATEGY


def gimmemotifs_annot2primpep(
    genome_fasta: str,
    annotation_gtf: str,
    output_dir: str,
    *,
    name: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the primary peptide FASTA of an assembly (``annot2primpep``).

    gffread turns the genome and the GTF annotation into peptides and only the longest protein per
    gene is kept. Records are named ``gene_name|gene_id`` (``.`` when the name is missing). The
    genome and annotation are copied into ``output_dir/<name>``; the result is
    ``output_dir/prim_genes/<name>.pep.fa``. gffread must be on PATH.
    """
    _program("gffread", "Install gffread with: conda install -c bioconda gffread")
    fasta = _file(genome_fasta, "genome_fasta")
    gtf = _file(annotation_gtf, "annotation_gtf")
    label = fasta.stem if name is None else name
    if not isinstance(label, str) or not label or "/" in label or "'" in label:
        raise ConfigurationError("name must be a plain identifier.")
    folder = workspace(output_dir)
    assembly = folder / label
    assembly.mkdir()
    _copy(fasta, assembly / f"{label}.fa")
    _copy(gtf, assembly / f"{label}.annotation.gtf")
    orthologs = _require("gimmemotifs.orthologs")
    orthologs.annot2primpep(label, str(folder))
    peptides = folder / "prim_genes" / f"{label}.pep.fa"
    records = sum(
        1
        for line in _written(peptides, "the primary peptide FASTA")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.startswith(">")
    )
    return _result("annot2primpep", {"peptides": peptides}, assembly=label, protein_count=records)


def gimmemotifs_load_orthogroups_in_db(
    database_path: str,
    assemblies: Sequence[str],
    orthofinder_results_dir: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Store OrthoFinder orthogroups in a new SQLite database (``load_orthogroups_in_db``).

    ``orthofinder_results_dir`` must contain ``Orthogroups/Orthogroups.tsv`` and
    ``Orthogroups/Orthogroups_UnassignedGenes.tsv`` with one ``<assembly>.pep`` column per
    assembly. The database path must not exist yet (the upstream function would replace it).
    """
    names = _names(assemblies, "assemblies")
    results = Path(orthofinder_results_dir).expanduser().resolve()
    groups = _file(str(results / "Orthogroups" / "Orthogroups.tsv"), "Orthogroups.tsv")
    unassigned = _file(
        str(results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"),
        "Orthogroups_UnassignedGenes.tsv",
    )
    pandas = _require("pandas")
    for table in (groups, unassigned):
        columns = set(pandas.read_table(table, nrows=0).columns)
        missing = [f"{name}.pep" for name in names if f"{name}.pep" not in columns]
        if missing:
            raise ConfigurationError(f"{table.name} has no column for {missing}.")
    database = _fresh(database_path, "database_path")
    orthologs = _require("gimmemotifs.orthologs")
    orthologs.load_orthogroups_in_db(str(database), names, str(results))
    with closing(sqlite3.connect(str(database))) as connection:
        orthogroup_count = connection.execute("SELECT COUNT(*) FROM orthogroups").fetchone()[0]
        gene_counts = {
            name: connection.execute(
                "SELECT COUNT(*) FROM genes WHERE assembly = ?", (name,)
            ).fetchone()[0]
            for name in names
        }
    return _result(
        "load_orthogroups_in_db",
        {"database": _written(database, "the ortholog database")},
        orthogroup_count=int(orthogroup_count),
        gene_counts={key: int(value) for key, value in gene_counts.items()},
    )


def gimmemotifs_factor2orthogroups(
    factor: str,
    references: Sequence[str],
    database_path: str,
    *,
    strategy: str = "strict",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the orthogroups of a transcription factor in an ortholog database.

    The factor is matched against gene names and ids of the ``references`` assemblies in the
    database. Only the ``strict`` strategy is offered here: the upstream medium and lenient
    strategies look names up on mygene.info.
    """
    _strict_only(strategy)
    name = _names([factor], "factor")[0]
    assemblies = _names(references, "references")
    database = _file(database_path, "database_path")
    orthologs = _require("gimmemotifs.orthologs")
    groups = orthologs.factor2orthogroups(name, tuple(assemblies), str(database), _OFFLINE_STRATEGY)
    return _result(
        "factor2orthogroups", {}, factor=name, orthogroups=sorted({str(group) for group in groups})
    )


def gimmemotifs_make_motif2factors(
    output_prefix: str,
    *,
    motif_file: str,
    motif_factors: Mapping[str, Sequence[str]],
    database_path: str,
    new_reference: str,
    database_references: Sequence[str],
    strategy: str = "strict",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a motif-to-factor table for a new assembly from an ortholog database.

    ``motif_factors`` maps each motif id to its factor names. For every factor the orthogroups are
    looked up in the ``database_references`` assemblies, and every gene of ``new_reference`` in
    those orthogroups is written as a factor of the motif. Output: ``<prefix>.motif2factors.txt``
    and ``<prefix>.pfm``. Only the ``strict`` strategy is available offline.
    """
    _strict_only(strategy)
    motifs = list(_motifs(_file(motif_file, "motif_file")).values())
    if not isinstance(motif_factors, Mapping) or not motif_factors:
        raise ConfigurationError("motif_factors must map motif ids to factor names.")
    factors = {
        str(motif): _names(list(names), f"motif_factors[{motif}]", minimum=0)
        for motif, names in motif_factors.items()
    }
    references = _names(database_references, "database_references")
    reference = _names([new_reference], "new_reference")[0]
    database = _file(database_path, "database_path")
    table = _fresh(f"{output_prefix}.motif2factors.txt", "output_prefix")
    pfm = _fresh(f"{output_prefix}.pfm", "output_prefix")
    orthologs = _require("gimmemotifs.orthologs")
    orthologs.make_motif2factors(
        str(Path(output_prefix).expanduser().resolve()),
        new_reference=reference,
        database_references=references,
        motifsandfactors=factors,
        database=str(database),
        strategy=_OFFLINE_STRATEGY,
        motifs=motifs,
    )
    return _result(
        "make_motif2factors",
        {
            "motif2factors": _written(table, "the motif2factors table"),
            "pfm": _written(pfm, "the motif PFM file"),
        },
        motif_count=len(motifs),
    )


def gimmemotifs_motif2factor_from_orthologs(
    genomes: Mapping[str, Mapping[str, str]],
    database_references: Sequence[str],
    new_references: Sequence[str],
    motif_file: str,
    motif_factors: Mapping[str, Sequence[str]],
    output_dir: str,
    *,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build motif-to-factor tables for new assemblies from local genomes.

    ``genomes`` maps each assembly name to ``{"fasta": path, "annotation": path}`` (the genome and
    its GTF annotation). The pipeline runs annot2primpep on every assembly, OrthoFinder on the
    peptides (``orthofinder`` on PATH), stores the orthogroups in a SQLite database and writes
    ``<new reference>.<motif file name>.motif2factors.txt`` for each name in ``new_references``.
    ``database_references`` are the assemblies on which the factor names of ``motif_factors`` are
    based. Nothing is downloaded. gffread and orthofinder are required.
    """
    _program("gffread", "Install gffread with: conda install -c bioconda gffread")
    orthofinder = _program(
        "orthofinder", "Install OrthoFinder with: conda install -c bioconda orthofinder"
    )
    if not isinstance(genomes, Mapping) or not genomes:
        raise ConfigurationError("genomes must map assembly names to fasta and annotation paths.")
    names = _names(list(genomes), "genomes")
    references = _names(database_references, "database_references")
    targets = _names(new_references, "new_references")
    for name in [*references, *targets]:
        if name not in names:
            raise ConfigurationError(f"{name!r} is not in genomes.")
    sources: dict[str, tuple[Path, Path]] = {}
    for name in names:
        entry = genomes[name]
        if not isinstance(entry, Mapping) or "fasta" not in entry or "annotation" not in entry:
            raise ConfigurationError(f"genomes[{name!r}] needs 'fasta' and 'annotation' paths.")
        sources[name] = (
            _file(entry["fasta"], f"{name} fasta"),
            _file(entry["annotation"], f"{name} annotation"),
        )
    motif_path = _file(motif_file, "motif_file")
    motifs = list(_motifs(motif_path).values())
    factors = {
        str(motif): _names(list(values), f"motif_factors[{motif}]", minimum=0)
        for motif, values in motif_factors.items()
    }
    folder = workspace(output_dir)
    work = folder / "intermediate"
    annotated = work / "annotated"
    annotated.mkdir(parents=True)
    orthologs = _require("gimmemotifs.orthologs")
    for name, (fasta, gtf) in sources.items():
        assembly = annotated / name
        assembly.mkdir()
        _copy(fasta, assembly / f"{name}.fa")
        _copy(gtf, assembly / f"{name}.annotation.gtf")
        orthologs.annot2primpep(name, str(annotated))
    results = _run_orthofinder(
        orthofinder, annotated / "prim_genes", int(_positive(threads, "threads"))
    )
    database = _fresh(str(work / "orthologs.sqlite"), "orthologs database")
    orthologs.load_orthogroups_in_db(str(database), names, str(results))
    outputs: dict[str, Path] = {}
    for reference in targets:
        prefix = folder / f"{reference}.{motif_path.name}"
        orthologs.make_motif2factors(
            str(prefix),
            new_reference=reference,
            database_references=references,
            motifsandfactors=factors,
            database=str(database),
            strategy=_OFFLINE_STRATEGY,
            motifs=motifs,
        )
        outputs[f"{reference}_motif2factors"] = _written(
            Path(f"{prefix}.motif2factors.txt"), "a motif2factors table"
        )
    return _result("motif2factor_from_orthologs", outputs, orthofinder_results=str(results))


__all__ = [
    "gimmemotifs_annot2primpep",
    "gimmemotifs_best_motif_in_cluster",
    "gimmemotifs_denovo_create_background",
    "gimmemotifs_denovo_create_backgrounds",
    "gimmemotifs_denovo_motif_report",
    "gimmemotifs_diff_motifs",
    "gimmemotifs_diff_plot",
    "gimmemotifs_factor2orthogroups",
    "gimmemotifs_filter_significant_motifs",
    "gimmemotifs_load_orthogroups_in_db",
    "gimmemotifs_make_motif2factors",
    "gimmemotifs_motif2factor_from_orthologs",
    "gimmemotifs_motif_tree_plot",
    "gimmemotifs_roc_html_report",
    "gimmemotifs_roc_plots",
    "gimmemotifs_select_nonredundant_motifs",
]

_register_public_verbose(globals())
