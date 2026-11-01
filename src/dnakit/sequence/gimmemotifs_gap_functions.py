"""GimmeMotifs operations not covered by the existing DNAKit motif helpers (Fig. 2 coverage gap).

Two kinds of wrapper live here:

- Command-line operations (``scan``, ``logo``, ``location``, ``cluster``) run the installed
  ``gimme``
  program with a fixed argument list and no shell. The program is resolved from ``gimme_path``,
  from PATH, or from the directory of the running Python interpreter.
- Library operations call the gimmemotifs Python API directly (peak combination, coverage tables,
  plots, background sequences, reports, JASPAR lookup, de novo motif tools).

Operations whose upstream code fails in this environment are not wrapped (see the disposition file):
``gimme threshold`` (Motif.hash passes a str to hashlib), ``gimme match`` (numpy broadcasting error
for equal-length motifs) and the automatic BAM indexing in ``coverage_table`` (undefined name).
"""

from __future__ import annotations

import importlib.metadata
import os
import shutil
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import (
    ExecutionConfig,
    artifacts,
    require,
    run,
    run_to_file,
    workspace,
)
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "gimmemotifs"
_SCAN_FORMATS = {"table": ("-t", "tsv"), "bed": ("-b", "bed"), "score_table": ("-T", "tsv")}
_LOGO_KINDS = ("information", "frequency", "energy", "ensembl")
_BACKGROUND_TYPES = ("random", "genomic", "gc")
_BACKGROUND_FORMATS = ("fasta", "bed")
_DENOVO_TOOLS = (
    "meme",
    "homer",
    "dreme",
    "weeder",
    "trawler",
    "improbizer",
    "gadem",
    "xxmotif",
    "amd",
    "posmo",
    "yamda",
    "motifsampler",
    "chipmunk",
    "hms",
    "dinamo",
    "prosampler",
    "rpmcmc",
    "bioprospector",
    "memew",
    "mdmodule",
)
_LOG_TAIL = 2000


def _require(module: str) -> Any:
    return require(module, _EXTRA)


def _version() -> str | None:
    try:
        return importlib.metadata.version("gimmemotifs")
    except importlib.metadata.PackageNotFoundError:
        return None


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} does not exist: {path}")
    return path


def _fresh_path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _positive(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _unit_interval(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value <= 1:
        raise ConfigurationError(f"{name} must be a number in (0, 1].")
    return float(value)


def _plain(value: Any) -> Any:
    """Convert numpy and pandas scalars to JSON-compatible Python values."""
    item = getattr(value, "item", None)
    if callable(item) and not isinstance(value, (str, bytes)):
        try:
            return item()
        except (TypeError, ValueError):
            return value
    return value


def _genome_fasta(value: object) -> str:
    """Genome for gimmemotifs: a FASTA path only (no registered genome lookup or download)."""
    return str(_existing_file(value, "genome"))


def _centred_peaks(peak_file: Path, target: Path) -> Path:
    """Copy peaks with a '+' strand column, so upstream windows are centred on each peak.

    Upstream reverses the window for any other strand value (start > end) and then fails, so the
    strand is normalised here; only the peak centre is used, so the result is unchanged.
    """
    lines: list[str] = []
    for number, raw in enumerate(peak_file.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip() or raw.startswith("#") or raw.startswith("track"):
            continue
        fields = raw.split("\t")
        if len(fields) < 3:
            raise ConfigurationError(f"Peak line {number} needs chrom, start and end columns.")
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError as exc:
            raise ConfigurationError(f"Peak line {number} has non-integer coordinates.") from exc
        if end <= start:
            raise ConfigurationError(f"Peak line {number} has end <= start.")
        name = fields[3] if len(fields) > 3 and fields[3] else "peak"
        score = fields[4] if len(fields) > 4 and fields[4] else "0"
        lines.append(f"{fields[0]}\t{start}\t{end}\t{name}\t{score}\t+")
    if not lines:
        raise ConfigurationError("The peak file has no peak records.")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _tail(text: str) -> str:
    return text[-_LOG_TAIL:]


def _gimme_program(gimme_path: str | Path | None) -> str:
    if gimme_path is not None:
        candidate = Path(gimme_path).expanduser()
        if not (candidate.is_file() and os.access(candidate, os.X_OK)):
            raise BackendUnavailableError(
                f"gimme executable is not available at {candidate}.",
                hint="Pass the path of an installed gimme program, or install GimmeMotifs.",
            )
        return str(candidate.resolve())
    found = shutil.which("gimme")
    if found is not None:
        return str(Path(found).resolve())
    sibling = Path(sys.executable).with_name("gimme")
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling.resolve())
    raise BackendUnavailableError(
        "The gimme command-line program was not found.",
        hint='Install GimmeMotifs with: python -m pip install "dnakit[gimmemotifs]", '
        "or pass gimme_path.",
    )


def _parse_table(text: str) -> list[dict[str, str]]:
    """Parse a tab-separated table with optional '#' comment lines and a header row."""
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    if not header[0]:
        header[0] = "sequence"
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        values = line.split("\t")
        rows.append({key: values[i] if i < len(values) else "" for i, key in enumerate(header)})
    return rows


def _numbered(paths: Sequence[Path], prefix: str) -> dict[str, Path]:
    return {f"{prefix}_{index}": path for index, path in enumerate(paths)}


def _result(operation: str, outputs: Mapping[str, object], log: str = "", **extra: Any) -> dict:
    return {
        "status": "complete",
        "operation": operation,
        "backend": "gimmemotifs",
        "version": _version(),
        "outputs": {key: str(value) for key, value in outputs.items()},
        "log_tail": _tail(log),
        **extra,
    }


def _motifs_from_file(pfm_file: Path) -> dict[str, Any]:
    motif_module = _require("gimmemotifs.motif")
    return {motif.id: motif for motif in motif_module.read_motifs(str(pfm_file), fmt="pfm")}


def _motif(pfm_file: Path, motif_id: str) -> Any:
    motifs = _motifs_from_file(pfm_file)
    if motif_id not in motifs:
        raise ConfigurationError(f"Motif {motif_id!r} is not in {pfm_file}.")
    return motifs[motif_id]


def gimmemotifs_pfm_scan(
    input_file: str,
    *,
    pfm_file: str,
    output_dir: str,
    cutoff: float,
    nreport: int = 1,
    reverse_complement: bool = True,
    output_format: str = "table",
    threads: int = 1,
    gimme_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan FASTA or BED sequences with PFM motifs using a fixed score cutoff (``gimme scan``).

    The cutoff is required: the FPR mode of ``gimme scan`` needs a background set (-B or a
    genome), which this wrapper does not take. ``output_format`` is ``table``, ``bed`` or
    ``score_table``. The scan output is written to ``output_dir``; for tabular formats the parsed
    rows are returned.
    """
    source = _existing_file(input_file, "input_file")
    pfm = _existing_file(pfm_file, "pfm_file")
    if output_format not in _SCAN_FORMATS:
        raise ConfigurationError(f"output_format must be one of {sorted(_SCAN_FORMATS)}.")
    flag, extension = _SCAN_FORMATS[output_format]
    score = _unit_interval(cutoff, "cutoff")
    reports = _positive(nreport, "nreport")
    workers = _positive(threads, "threads")
    directory = workspace(output_dir)
    output = directory / f"scan.{extension}"
    arguments = [
        "scan",
        "-p",
        str(pfm),
        "-c",
        repr(score),
        "-n",
        str(reports),
        "-N",
        str(workers),
        flag,
    ]
    if not reverse_complement:
        arguments.append("-r")
    arguments.append(str(source))
    log = run_to_file(
        _gimme_program(gimme_path), arguments, output, ExecutionConfig(threads=workers)
    )
    rows = _parse_table(output.read_text(encoding="utf-8")) if output_format != "bed" else []
    return _result(
        "pfm_scan",
        {"scan": output},
        log,
        rows=rows,
        row_count=len(rows),
    )


def gimmemotifs_match_motifs(
    pfm_file: str,
    *,
    database_pfm: str,
    output_path: str,
    nreport: int = 1,
    gimme_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the best matching database motifs for each query motif (``gimme match``, table output).

    The graphical report (``gimme match -o``) is not offered: it fails in this gimmemotifs release
    with a numpy broadcasting error. The table is written to ``output_path`` (must not exist) and
    the parsed rows are returned.
    """
    query = _existing_file(pfm_file, "pfm_file")
    database = _existing_file(database_pfm, "database_pfm")
    reports = _positive(nreport, "nreport")
    target = _fresh_path(output_path, "output_path")
    arguments = ["match", "-d", str(database), "-n", str(reports), str(query)]
    log = run_to_file(_gimme_program(gimme_path), arguments, target, ExecutionConfig())
    rows = _parse_table(target.read_text(encoding="utf-8"))
    return _result("motif_match", {"table": target}, log, rows=rows, row_count=len(rows))


def gimmemotifs_logo_images(
    pfm_file: str,
    *,
    output_dir: str,
    motif_ids: Sequence[str] | None = None,
    kind: str | None = None,
    show_title: bool = True,
    gimme_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw sequence logos for the motifs of a PFM file (``gimme logo``); returns the PNG paths."""
    pfm = _existing_file(pfm_file, "pfm_file")
    if kind is not None and kind not in _LOGO_KINDS:
        raise ConfigurationError(f"kind must be one of {list(_LOGO_KINDS)}.")
    arguments = ["logo", "-p", str(pfm)]
    if motif_ids:
        arguments += ["-i", ",".join(str(item) for item in motif_ids)]
    if kind is not None:
        arguments += ["-k", kind]
    if not show_title:
        arguments.append("--notitle")
    directory = workspace(output_dir)
    log = run(_gimme_program(gimme_path), arguments, directory, ExecutionConfig())
    images = sorted(path for path in artifacts(directory).values() if path.suffix.lower() == ".png")
    if not images:
        raise BackendExecutionError("gimme logo wrote no PNG files.", context={"log": _tail(log)})
    return _result("logo_images", _numbered(images, "image"), log)


def gimmemotifs_motif_location_plots(
    fasta_file: str,
    *,
    pfm_file: str,
    output_dir: str,
    size: int | None = None,
    motif_ids: Sequence[str] | None = None,
    cutoff: float = 0.95,
    gimme_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot where each motif occurs relative to the sequence centre (``gimme location``)."""
    fasta = _existing_file(fasta_file, "fasta_file")
    pfm = _existing_file(pfm_file, "pfm_file")
    score = _unit_interval(cutoff, "cutoff")
    arguments = ["location"]
    if size is not None:
        arguments += ["-s", str(_positive(size, "size"))]
    if motif_ids:
        arguments += ["-i", ",".join(str(item) for item in motif_ids)]
    arguments += ["-c", repr(score), str(pfm), str(fasta)]
    directory = workspace(output_dir)
    log = run(_gimme_program(gimme_path), arguments, directory, ExecutionConfig())
    files = sorted(artifacts(directory).values())
    return _result("motif_location_plots", _numbered(files, "file"), log)


def gimmemotifs_cluster_pfm(
    input_pfm: str,
    *,
    output_dir: str,
    threshold: float | None = None,
    compare_reverse_complement: bool = True,
    threads: int = 1,
    gimme_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster PFM motifs by similarity; write the clustered PFM and report (``gimme cluster``)."""
    pfm = _existing_file(input_pfm, "input_pfm")
    workers = _positive(threads, "threads")
    target = _fresh_path(output_dir, "output_dir")
    arguments = ["cluster"]
    if threshold is not None:
        arguments += ["-t", repr(_unit_interval(threshold, "threshold"))]
    if not compare_reverse_complement:
        arguments.append("-s")
    arguments += ["-N", str(workers), str(pfm), str(target)]
    log = run(
        _gimme_program(gimme_path), arguments, target.parent, ExecutionConfig(threads=workers)
    )
    outputs = {
        "clustered_pfm": target / "clustered_motifs.pfm",
        "report": target / "gimme.clustered.html",
    }
    missing = [name for name, path in outputs.items() if not path.is_file()]
    if missing:
        raise BackendExecutionError(
            f"gimme cluster did not write {missing}.", context={"log": _tail(log)}
        )
    return _result("cluster_pfm", outputs, log)


def gimmemotifs_combine_peak_summits(
    peak_files: Sequence[str],
    *,
    genome_fasta: str,
    window: int = 200,
    scale_value: bool = False,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine MACS2 summit or narrowPeak files, keeping the top summit per window.

    Wraps ``combine_peaks``.

    ``genome_fasta`` is the FASTA of the genome (combine_peaks indexes it with pyfaidx). Returns the
    combined table rows and
    optionally writes them as a headerless tab-separated file (``output_path``, must not exist).
    """
    if isinstance(peak_files, (str, Path)) or not peak_files:
        raise ConfigurationError("peak_files must be a non-empty list of file paths.")
    files = [str(_existing_file(item, "peak_files item")) for item in peak_files]
    genome = _existing_file(genome_fasta, "genome_fasta")
    span = _positive(window, "window")
    preprocessing = _require("gimmemotifs.preprocessing")
    frame = preprocessing.combine_peaks(files, str(genome), span, bool(scale_value))
    rows = [
        {str(key): _plain(value) for key, value in record.items()}
        for record in frame.to_dict(orient="records")
    ]
    outputs: dict[str, Path] = {}
    if output_path is not None:
        target = _fresh_path(output_path, "output_path")
        frame.to_csv(target, sep="\t", index=False, header=False)
        outputs["table"] = target
    return _result(
        "combine_peak_summits",
        outputs,
        columns=[str(column) for column in frame.columns],
        rows=rows,
        row_count=len(rows),
    )


def gimmemotifs_coverage_table(
    peak_file: str,
    data_files: Sequence[str],
    *,
    output_path: str | None = None,
    window: int = 200,
    log_transform: bool = True,
    normalization: str | None = None,
    top: int = -1,
    top_method: str = "var",
    remove_duplicates: bool = True,
    remove_repeats: bool = True,
    threads: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count reads from indexed BAM or BED data files in each peak region (``coverage_table``).

    BAM inputs must already have a ``.bai`` index: the upstream automatic indexing refers to an
    undefined name and is not used. Returns the table rows; ``output_path`` writes a new TSV file.
    """
    peaks = _existing_file(peak_file, "peak_file")
    if isinstance(data_files, (str, Path)) or not data_files:
        raise ConfigurationError("data_files must be a non-empty list of file paths.")
    files = [str(_existing_file(item, "data_files item")) for item in data_files]
    for item in files:
        if item.endswith(".bam") and not Path(item + ".bai").is_file():
            raise ConfigurationError(f"BAM file needs an index first (samtools index): {item}")
    if normalization not in (None, "quantile", "scale"):
        raise ConfigurationError("normalization must be None, 'quantile' or 'scale'.")
    if top_method not in ("var", "std", "mean", "random"):
        raise ConfigurationError("top_method must be 'var', 'std', 'mean' or 'random'.")
    preprocessing = _require("gimmemotifs.preprocessing")
    with tempfile.TemporaryDirectory(prefix="gimmemotifs-peaks-") as scratch:
        centred = _centred_peaks(peaks, Path(scratch) / "peaks.bed")
        frame = preprocessing.coverage_table(
            str(centred),
            files,
            window=_positive(window, "window"),
            log_transform=bool(log_transform),
            normalization=normalization,
            top=int(top),
            topmethod=top_method,
            rmdup=bool(remove_duplicates),
            rmrepeats=bool(remove_repeats),
            ncpus=_positive(threads, "threads"),
        )
    outputs: dict[str, Path] = {}
    if output_path is not None:
        target = _fresh_path(output_path, "output_path")
        frame.to_csv(target, sep="\t")
        outputs["table"] = target
    rows = [
        {str(key): _plain(value) for key, value in record.items()}
        for record in frame.reset_index().to_dict(orient="records")
    ]
    return _result(
        "coverage_table",
        outputs,
        columns=[str(column) for column in frame.columns],
        rows=rows,
        row_count=len(rows),
    )


def gimmemotifs_histogram_plot(
    values: Sequence[float],
    output_path: str,
    *,
    x_range: tuple[float, float] | None = None,
    bins: int = 10,
    title: str | None = None,
    x_label: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a histogram of numeric values as SVG (``gimmemotifs.plot.plot_histogram``).

    gimmemotifs appends ``.svg`` to the output name when it is missing; the returned path is the
    file that was written.
    """
    data = [float(value) for value in values]
    if not data:
        raise ConfigurationError("values must not be empty.")
    if x_range is not None and (len(x_range) != 2 or not x_range[0] < x_range[1]):
        raise ConfigurationError("x_range must be (low, high) with low < high.")
    target = output_path if str(output_path).endswith(".svg") else f"{output_path}.svg"
    final = _fresh_path(target, "output_path")
    plot = _require("gimmemotifs.plot")
    plot.plot_histogram(
        data,
        str(final)[: -len(".svg")],
        xrange=x_range,
        breaks=_positive(bins, "bins"),
        title=title,
        xlabel=x_label,
    )
    return _result(
        "histogram_plot",
        {"image": final},
        count=len(data),
        minimum=min(data),
        maximum=max(data),
    )


def gimmemotifs_roc_curve_plot(
    output_path: str,
    *,
    plot_x: Sequence[Any],
    plot_y: Sequence[Any],
    curve_ids: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot one or several ROC curves (``gimmemotifs.plot.roc_plot``) with 1 - specificity on x.

    ``plot_x`` and ``plot_y`` are either one curve (a list of numbers) or a list of curves. A PNG
    extension is added when the path has no image extension; the written path is returned.
    """
    if not plot_x or not plot_y or len(plot_x) != len(plot_y):
        raise ConfigurationError("plot_x and plot_y must be non-empty and of equal length.")
    plot = _require("gimmemotifs.plot")
    import numpy as np

    multiple = isinstance(plot_x[0], (list, tuple))
    if multiple:
        xs = [np.asarray(curve, dtype=float) for curve in plot_x]
        ys = [np.asarray(curve, dtype=float) for curve in plot_y]
    else:
        xs = [float(value) for value in plot_x]
        ys = [float(value) for value in plot_y]
    valid = list(getattr(plot, "VALID_EXTENSIONS", [".png", ".pdf", ".svg", ".ps"]))
    name = str(output_path)
    if os.path.splitext(name)[-1] not in valid:
        name += ".png"
    final = _fresh_path(name, "output_path")
    plot.roc_plot(str(final), xs, ys, ids=list(curve_ids or []))
    return _result("roc_curve_plot", {"image": final}, curves=len(xs) if multiple else 1)


def gimmemotifs_match_plot(
    matches: Sequence[Mapping[str, Any]],
    output_path: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot each query motif with its best database match and p-value (``match_plot``).

    Each item needs ``query_pfm``, ``query_id``, ``database_pfm``, ``database_id`` and ``p_value``.
    """
    if not matches:
        raise ConfigurationError("matches must not be empty.")
    plotdata = []
    for item in matches:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Each match must be a mapping.")
        query = _motif(
            _existing_file(item.get("query_pfm"), "query_pfm"), str(item.get("query_id"))
        )
        target = _motif(
            _existing_file(item.get("database_pfm"), "database_pfm"),
            str(item.get("database_id")),
        )
        p_value = item.get("p_value")
        if isinstance(p_value, bool) or not isinstance(p_value, (int, float)):
            raise ConfigurationError("p_value must be a number.")
        plotdata.append((query, target, float(p_value)))
    final = _fresh_path(output_path, "output_path")
    _require("gimmemotifs.plot").match_plot(plotdata, str(final))
    return _result("match_plot", {"image": final}, pairs=len(plotdata))


def gimmemotifs_motif_localization(
    fasta_file: str,
    *,
    pfm_file: str,
    motif_id: str,
    output_path: str,
    size: int,
    cutoff: float = 0.9,
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of motif positions in a FASTA file with a Kolmogorov-Smirnov p-value.

    Returns the motif id and p-value. The SVG is written only when the motif has matches (upstream
    behaviour); ``image`` is absent otherwise.
    """
    fasta = _existing_file(fasta_file, "fasta_file")
    pfm = _existing_file(pfm_file, "pfm_file")
    motif = _motif(pfm, str(motif_id))
    span = _positive(size, "size")
    if span <= len(motif):
        raise ConfigurationError("size must be larger than the motif length.")
    score = _unit_interval(cutoff, "cutoff")
    final = _fresh_path(output_path, "output_path")
    name = str(final)[: -len(".svg")] if str(final).endswith(".svg") else str(final)
    utils = _require("gimmemotifs.utils")
    found_id, p_value = utils.motif_localization(str(fasta), motif, span, name, cutoff=score)
    svg = Path(name + ".svg")
    outputs = {"image": svg} if svg.is_file() else {}
    return _result("motif_localization", outputs, motif_id=str(found_id), p_value=float(p_value))


def gimmemotifs_cluster_motif_report(
    pfm_file: str,
    *,
    output_dir: str,
    threshold: float,
    title: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster motifs by similarity and write an HTML report with logos.

    Wraps ``cluster_motifs_with_report``.

    The clustered PFM is written to ``output_dir/clustered_motifs.pfm``. Returns each cluster's
    representative and members.
    """
    pfm = _existing_file(pfm_file, "pfm_file")
    score = _unit_interval(threshold, "threshold")
    directory = workspace(output_dir)
    (directory / "images").mkdir()
    outfile = directory / "clustered_motifs.pfm"
    cluster_module = _require("gimmemotifs.motif.cluster")
    clusters = cluster_module.cluster_motifs_with_report(
        str(pfm), str(outfile), str(directory), score, title=title
    )
    described = [
        {
            "representative": str(getattr(group[0], "id", group[0])),
            "members": [str(getattr(member, "id", member)) for member in group[1]],
        }
        for group in clusters
    ]
    outputs: dict[str, Path] = {"clustered_pfm": outfile}
    report = directory / "gimme.clustereds.html"
    if report.is_file():
        outputs["report"] = report
    return _result(
        "cluster_motif_report", outputs, clusters=described, cluster_count=len(described)
    )


def gimmemotifs_motif_image_series(
    motif_ids: Sequence[str],
    *,
    pfm_file: str,
    output_dir: str,
    subdir: str = "logos",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one logo PNG per listed motif id (``gimmemotifs.report.motif_to_img_series``)."""
    if isinstance(motif_ids, (str, Path)) or not motif_ids:
        raise ConfigurationError("motif_ids must be a non-empty list of motif ids.")
    pfm = _existing_file(pfm_file, "pfm_file")
    motifs = _motifs_from_file(pfm)
    missing = [item for item in motif_ids if item not in motifs]
    if missing:
        raise ConfigurationError(f"Motifs not found in {pfm}: {missing}")
    pandas = _require("pandas")
    report = _require("gimmemotifs.report")
    directory = workspace(output_dir)
    index = pandas.Index([str(item) for item in motif_ids])
    series = report.motif_to_img_series(index, motifs=motifs, outdir=str(directory), subdir=subdir)
    images = {str(key): str(directory / str(value)) for key, value in series.items()}
    for path in images.values():
        if not Path(path).is_file():
            raise BackendExecutionError(f"Logo image was not written: {path}")
    return _result("motif_image_series", {}, images=images, image_count=len(images))


def gimmemotifs_create_background(
    output_path: str,
    *,
    background_type: str,
    fmt: str = "fasta",
    size: int | None = None,
    genome: str | None = None,
    input_file: str | None = None,
    number: int = 10000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create background sequences for motif analysis (``create_background_file``).

    ``background_type`` is ``random`` (Markov model of ``input_file``), ``genomic`` (random genome
    regions) or ``gc`` (GC-matched regions of ``input_file``). ``genome`` is a genomepy name or a
    FASTA path. The ``promoter`` type needs a genomepy annotation registered under the genome name
    and is rejected here. Upstream exits the process on invalid arguments, so inputs are
    checked first.
    """
    if background_type == "promoter":
        raise ConfigurationError(
            "The promoter background needs a genomepy annotation registered under the genome name."
        )
    if background_type not in _BACKGROUND_TYPES:
        raise ConfigurationError(f"background_type must be one of {list(_BACKGROUND_TYPES)}.")
    if fmt not in _BACKGROUND_FORMATS:
        raise ConfigurationError(f"fmt must be one of {list(_BACKGROUND_FORMATS)}.")
    if background_type == "random" and fmt == "bed":
        raise ConfigurationError("A random background can only be written as FASTA.")
    target = _fresh_path(output_path, "output_path")
    count = _positive(number, "number")
    if background_type == "random":
        if input_file is None:
            raise ConfigurationError("A random background needs input_file (its Markov model).")
        input_file = str(_existing_file(input_file, "input_file"))
    if background_type in ("genomic", "gc"):
        if not genome:
            raise ConfigurationError(f"A {background_type} background needs genome (a FASTA path).")
        genome = _genome_fasta(genome)
    if background_type == "gc":
        if input_file is None:
            raise ConfigurationError("A gc background needs input_file (sequences to match).")
        input_file = str(_existing_file(input_file, "input_file"))
    if background_type == "genomic" and size is None and input_file is None:
        raise ConfigurationError("A genomic background needs size or input_file.")
    if size is not None:
        _positive(size, "size")
    background = _require("gimmemotifs.background")
    background.create_background_file(
        str(target),
        background_type,
        fmt=fmt,
        size=size,
        genome=genome,
        inputfile=input_file,
        number=count,
    )
    if not target.is_file():
        raise BackendExecutionError("gimmemotifs did not write the background file.")
    return _result("create_background", {"background": target}, sequence_count=count)


def gimmemotifs_gc_bin_index(
    genome: str,
    output_path: str,
    *,
    min_bin_size: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Index the GC content of genome windows for GC-matched sampling (``create_gc_bin_index``)."""
    genome_ref = _genome_fasta(genome)
    target = _fresh_path(output_path, "output_path")
    background = _require("gimmemotifs.background")
    background.create_gc_bin_index(
        genome_ref, str(target), min_bin_size=_positive(min_bin_size, "min_bin_size")
    )
    if not target.exists():
        raise BackendExecutionError("gimmemotifs did not write the GC index.")
    return _result("gc_bin_index", {"index": target})


def gimmemotifs_gc_bin_bed(
    output_path: str,
    *,
    genome: str,
    number: int,
    length: int = 200,
    bins: Sequence[float] | None = None,
    random_state: int | None = None,
    min_bin_size: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample BED regions spread evenly over GC bins of the genome (``gc_bin_bedfile``)."""
    genome_ref = _genome_fasta(genome)
    target = _fresh_path(output_path, "output_path")
    background = _require("gimmemotifs.background")
    background.gc_bin_bedfile(
        str(target),
        genome_ref,
        _positive(number, "number"),
        length=_positive(length, "length"),
        bins=list(bins) if bins is not None else None,
        random_state=random_state,
        min_bin_size=_positive(min_bin_size, "min_bin_size"),
    )
    if not target.is_file():
        raise BackendExecutionError("gimmemotifs did not write the GC-bin BED file.")
    return _result("gc_bin_bed", {"bed": target}, region_count=_positive(number, "number"))


def gimmemotifs_matched_gc_bed(
    output_path: str,
    *,
    match_bed: str,
    genome: str,
    number: int,
    size: int | None = None,
    min_bin_size: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample BED regions whose GC content follows the regions of ``match_bed``."""
    match = _existing_file(match_bed, "match_bed")
    genome_ref = _genome_fasta(genome)
    target = _fresh_path(output_path, "output_path")
    background = _require("gimmemotifs.background")
    background.matched_gc_bedfile(
        str(target),
        str(match),
        genome_ref,
        _positive(number, "number"),
        size=size,
        min_bin_size=_positive(min_bin_size, "min_bin_size"),
    )
    if not target.is_file():
        raise BackendExecutionError("gimmemotifs did not write the GC-matched BED file.")
    return _result("matched_gc_bed", {"bed": target})


def gimmemotifs_jaspar_motif_info(motif_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Fetch the metadata of a JASPAR matrix from the JASPAR web service (network).

    Calls ``gimmemotifs.utils.get_jaspar_motif_info``; the response is returned as the JSON object
    JASPAR sends. Needs network access; tests mock the request.
    """
    if not isinstance(motif_id, str) or not motif_id.strip():
        raise ConfigurationError("motif_id must be a non-empty string.")
    utils = _require("gimmemotifs.utils")
    return {"motif_id": motif_id, "info": utils.get_jaspar_motif_info(motif_id.strip())}


def _denovo_tools_installed(tools_module: Any) -> list[str]:
    installed: list[str] = []
    for name in _DENOVO_TOOLS:
        try:
            program = tools_module.get_tool(name)
        except ValueError:
            continue
        if program.is_installed():
            installed.append(name)
    return installed


def _require_denovo_tool() -> list[str]:
    tools = _require("gimmemotifs.tools")
    installed = _denovo_tools_installed(tools)
    if not installed:
        raise BackendUnavailableError(
            "No de novo motif program (MEME, HOMER, DREME, ...) is installed.",
            hint="Install at least one of MEME or HOMER (for example conda install -c bioconda "
            "meme homer), put it on PATH and retry.",
        )
    return installed


def gimmemotifs_denovo_motifs(
    input_file: str,
    *,
    output_dir: str,
    params: Mapping[str, Any] | None = None,
    filter_significant: bool = True,
    cluster: bool = True,
    create_report: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict de novo motifs with an ensemble of external tools (``gimme_motifs``).

    Requires at least one installed de novo program; otherwise BackendUnavailableError is raised
    before gimmemotifs is called. ``params`` follows gimmemotifs' parameter mapping.
    """
    source = _existing_file(input_file, "input_file")
    _require_denovo_tool()
    denovo = _require("gimmemotifs.motif.denovo")
    directory = workspace(output_dir)
    denovo.gimme_motifs(
        str(source),
        str(directory),
        params=dict(params) if params is not None else None,
        filter_significant=bool(filter_significant),
        cluster=bool(cluster),
        create_report=bool(create_report),
    )
    files = sorted(artifacts(directory).values())
    if not files:
        raise BackendExecutionError("gimme_motifs wrote no output files.")
    return _result("denovo_motifs", _numbered(files, "file"))


def gimmemotifs_denovo_predict(
    input_fasta: str,
    background_fasta: str,
    output_path: str,
    *,
    params: Mapping[str, Any] | None = None,
    stats_fg: Mapping[str, Any] | None = None,
    stats_bg: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the configured de novo tools on a FASTA file against a background."""
    source = _existing_file(input_fasta, "input_fasta")
    background = _existing_file(background_fasta, "background_fasta")
    target = _fresh_path(output_path, "output_path")
    _require_denovo_tool()
    prediction = _require("gimmemotifs.motif.prediction")
    prediction.predict_motifs(
        str(source),
        str(background),
        str(target),
        params=dict(params) if params is not None else None,
        stats_fg=dict(stats_fg) if stats_fg is not None else None,
        stats_bg=dict(stats_bg) if stats_bg is not None else None,
    )
    if not target.is_file():
        raise BackendExecutionError("predict_motifs wrote no motif file.")
    return _result("denovo_predict", {"motifs": target})


def gimmemotifs_parallel_denovo_predict(
    fasta_file: str,
    output_path: str,
    *,
    analysis: str = "small",
    organism: str = "hg19",
    single: bool = False,
    background: str = "",
    tools: str | None = None,
    threads: int = 8,
    max_time: int = -1,
    stats_fg: Mapping[str, Any] | None = None,
    stats_bg: Mapping[str, Any] | None = None,
    gc: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parallel de novo motif prediction over several external tools (``pp_predict_motifs``)."""
    source = _existing_file(fasta_file, "fasta_file")
    target = _fresh_path(output_path, "output_path")
    _require_denovo_tool()
    prediction = _require("gimmemotifs.motif.prediction")
    prediction.pp_predict_motifs(
        str(source),
        str(target),
        analysis=analysis,
        organism=organism,
        single=bool(single),
        background=background,
        tools=tools,
        ncpus=_positive(threads, "threads"),
        max_time=int(max_time),
        stats_fg=dict(stats_fg) if stats_fg is not None else None,
        stats_bg=dict(stats_bg) if stats_bg is not None else None,
        gc=bool(gc),
    )
    if not target.is_file():
        raise BackendExecutionError("pp_predict_motifs wrote no motif file.")
    return _result("parallel_denovo_predict", {"motifs": target})


__all__ = [
    "gimmemotifs_cluster_motif_report",
    "gimmemotifs_cluster_pfm",
    "gimmemotifs_combine_peak_summits",
    "gimmemotifs_coverage_table",
    "gimmemotifs_create_background",
    "gimmemotifs_denovo_motifs",
    "gimmemotifs_denovo_predict",
    "gimmemotifs_gc_bin_bed",
    "gimmemotifs_gc_bin_index",
    "gimmemotifs_histogram_plot",
    "gimmemotifs_jaspar_motif_info",
    "gimmemotifs_logo_images",
    "gimmemotifs_match_motifs",
    "gimmemotifs_match_plot",
    "gimmemotifs_matched_gc_bed",
    "gimmemotifs_motif_image_series",
    "gimmemotifs_motif_localization",
    "gimmemotifs_motif_location_plots",
    "gimmemotifs_parallel_denovo_predict",
    "gimmemotifs_pfm_scan",
    "gimmemotifs_roc_curve_plot",
]

_register_public_verbose(globals())
