"""NanoPlot command-line facade: the long-read plots and the per-read table of NanoPlot.

NanoPlot (https://github.com/wdecoster/NanoPlot; 1.48.0 is the tested version) draws, in one run, every plot that the
input allows: read length histograms, yield by length, read length versus quality, the alignment plots of BAM/CRAM
input and the time and channel plots of sequencing summaries. NanoPlot cannot be told to draw only one of them, so the
facade runs the program once in a new output directory and reports each plot under a stable name. The user installs
NanoPlot (``conda install -c bioconda nanoplot`` or ``pip install NanoPlot``) and puts ``NanoPlot`` on PATH or passes
``nanoplot_path``. The program is started with an explicit, validated argument list and no shell.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p12_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="NanoPlot",
    executable="NanoPlot",
    hint=(
        "Install NanoPlot 1.40 or later (for example `conda install -c bioconda nanoplot` or `pip install NanoPlot`) "
        "and put NanoPlot on PATH or pass nanoplot_path."
    ),
    tested_version="1.48.0",
    citation="https://doi.org/10.1093/bioinformatics/bty149",
    version_arguments=("--version",),
    version_pattern=r"NanoPlot\s+(\S+)",
)
_INPUT_FLAGS = {
    "fastq": "--fastq",
    "fasta": "--fasta",
    "fastq_rich": "--fastq_rich",
    "fastq_minimal": "--fastq_minimal",
    "summary": "--summary",
    "bam": "--bam",
    "ubam": "--ubam",
    "cram": "--cram",
    "feather": "--feather",
}
_BIVARIATE = ("kde", "hex", "dot")
_FORMATS = ("png", "jpg", "jpeg", "webp", "svg", "pdf", "eps", "json")
_COLOR = re.compile(r"^#?[A-Za-z0-9]{1,30}$")
_TITLE = re.compile(r"^[^\x00-\x1f]{1,200}$")

# Plot identifiers of this facade and the file-name prefixes NanoPlot gives them (see nanoplotter).
PLOT_FILES: dict[str, tuple[str, ...]] = {
    "read_length_histogram": (
        "Non_weightedHistogramReadlength",
        "WeightedHistogramReadlength",
        "Non_weightedLogTransformed_HistogramReadlength",
        "WeightedLogTransformed_HistogramReadlength",
    ),
    "yield_by_minimal_length": ("Yield_By_Length",),
    "read_length_vs_quality_scatter": ("LengthvsQualityScatterPlot",),
    "aligned_vs_sequenced_length_scatter": ("AlignedReadlengthvsSequencedReadLength",),
    "mapq_vs_base_quality_scatter": ("MappingQualityvsAverageBaseQuality",),
    "mapq_vs_read_length_scatter": ("MappingQualityvsReadLength",),
    "percent_identity_vs_base_quality_scatter": ("PercentIdentityvsAverageBaseQuality",),
    "percent_identity_vs_aligned_length_scatter": ("PercentIdentityvsAlignedReadLength",),
    "percent_identity_histogram": ("PercentIdentityHistogramDynamic_Histogram_percent_identity",),
    "reads_per_channel_activity_map": ("ActivityMap_ReadsPerChannel",),
    "reads_over_time": ("NumberOfReads_Over_Time",),
    "active_pores_over_time": ("ActivePores_Over_Time",),
    "cumulative_yield_gigabases": ("CumulativeYieldPlot_Gigabases",),
    "cumulative_yield_reads": ("CumulativeYieldPlot_NumberOfReads",),
    "read_length_over_time_violin": ("TimeLengthViolinPlot", "TimeLogLengthViolinPlot"),
    "quality_over_time_violin": ("TimeQualityViolinPlot",),
    "sequencing_speed_over_time_violin": ("TimeSequencingSpeed_ViolinPlot",),
}


def _parse_stats(path: Path) -> dict[str, Any]:
    """Read ``NanoStats.txt`` written with ``--tsv_stats`` (metric, value) into a dictionary."""
    stats: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0]:
            stats[fields[0]] = _cli.convert_number(fields[1])
    return stats


def _raw_table(path: Path) -> tuple[list[str], int]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        columns = handle.readline().rstrip("\n").split("\t")
        return columns, sum(1 for _ in handle)


def nanoplot_plots(
    input_paths: Sequence[str],
    output_dir: str,
    *,
    input_type: str = "fastq",
    plot_names: Sequence[str] | None = None,
    bivariate_plots: Sequence[str] | None = None,
    static_images: bool = True,
    image_formats: Sequence[str] | None = None,
    log_length: bool = False,
    percent_quality: bool = False,
    aligned_length: bool = False,
    min_length: int | None = None,
    max_length: int | None = None,
    min_quality: float | None = None,
    downsample: int | None = None,
    drop_outliers: bool = False,
    no_supplementary: bool = False,
    show_n50: bool | None = None,
    color: str | None = None,
    colormap: str | None = None,
    title: str | None = None,
    dpi: int | None = None,
    raw_table: bool = False,
    threads: int = 1,
    timeout_seconds: float = 86_400.0,
    nanoplot_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the NanoPlot long-read plots and statistics of one dataset and return the plot files by name (``NanoPlot``).

    ``input_paths`` are one or more files of ``input_type``: ``fastq``, ``fasta``, ``fastq_rich`` and ``fastq_minimal``
    (reads), ``summary`` (a Guppy/MinKNOW sequencing summary, which also gives the time and channel plots), ``bam``,
    ``ubam``, ``cram`` (alignments, which also give the alignment plots) or ``feather``. ``output_dir`` must not exist.
    NanoPlot writes every plot the input allows, as interactive HTML and, when ``static_images`` is true (the NanoPlot
    default), as static images in ``image_formats`` (default PNG); with ``static_images=False`` only HTML is written, which
    is much faster. The result groups the files under these plot names, which cover the plots of the work list:
    ``read_length_histogram`` (weighted and non-weighted, plus log-transformed with ``log_length``),
    ``yield_by_minimal_length``, ``read_length_vs_quality_scatter``, ``aligned_vs_sequenced_length_scatter``,
    ``mapq_vs_base_quality_scatter``, ``mapq_vs_read_length_scatter``, ``percent_identity_vs_base_quality_scatter``,
    ``percent_identity_vs_aligned_length_scatter``, ``percent_identity_histogram``, ``reads_per_channel_activity_map``,
    ``reads_over_time``, ``active_pores_over_time``, ``cumulative_yield_gigabases``, ``cumulative_yield_reads``,
    ``read_length_over_time_violin``, ``quality_over_time_violin`` and ``sequencing_speed_over_time_violin``.
    ``plot_names`` limits what is reported (NanoPlot still writes all of them). ``bivariate_plots`` selects the bivariate
    plot kinds among ``kde``, ``hex`` and ``dot`` (``--plots``). Filters and transforms: ``min_length``/``max_length``,
    ``min_quality``, ``downsample``, ``drop_outliers``, ``no_supplementary``, ``log_length`` (``--loglength``),
    ``percent_quality`` (``--percentqual``), ``aligned_length`` (``--alength``) and ``show_n50`` (True ``--N50``, False
    ``--no-N50``). Plot style: ``color``, ``colormap``, ``title`` and ``dpi``. ``raw_table=True`` also writes the per-read
    table (``--raw``: read id, quality, length, and for alignments aligned length, mapping quality) as
    ``NanoPlot-data.tsv.gz``.

    Result keys: ``output_dir``, ``plots`` (plot name to the list of files written), ``missing_plots`` (requested names
    for which NanoPlot wrote nothing, for example the time plots of a FASTQ input or the mapping-quality plots when the
    mapping qualities do not vary), ``statistics`` (the ``NanoStats`` table as a dictionary: ``number_of_reads``,
    ``number_of_bases``, ``median_read_length``, ``mean_read_length``, ``n50``, ``mean_qual``, ``median_qual``, ...),
    ``stats_file``, ``report_file`` and, with ``raw_table``, ``raw_table_file``, ``raw_table_columns`` and
    ``raw_table_rows``.
    """
    kind = _cli.choice(input_type, tuple(_INPUT_FLAGS), "input_type")
    inputs = _cli.existing_files(list(input_paths), "input_paths")
    directory = Path(output_dir).expanduser().resolve()
    if directory.exists():
        raise ConfigurationError(f"output_dir must not exist yet: {directory}")
    if not directory.parent.is_dir():
        raise ConfigurationError(f"The parent folder of output_dir does not exist: {directory.parent}")
    wanted = list(PLOT_FILES) if plot_names is None else list(plot_names)
    unknown = [name for name in wanted if name not in PLOT_FILES]
    if unknown or not wanted:
        raise ConfigurationError(
            f"plot_names must be a non-empty list of {sorted(PLOT_FILES)}; got {unknown}."
        )
    arguments = ["-t", str(_cli.threads_value(threads)), "-o", str(directory), "--tsv_stats"]
    if bivariate_plots is not None:
        kinds = list(bivariate_plots)
        if not kinds or any(item not in _BIVARIATE for item in kinds):
            raise ConfigurationError(
                f"bivariate_plots must be a non-empty list of {list(_BIVARIATE)}."
            )
        arguments += ["--plots", *kinds]
    if _cli.flag(static_images, "static_images"):
        formats = list(image_formats) if image_formats is not None else ["png"]
        if not formats or any(item not in _FORMATS for item in formats):
            raise ConfigurationError(f"image_formats must be a non-empty list of {list(_FORMATS)}.")
        arguments += ["-f", *formats]
    else:
        arguments.append("--no_static")
    for switch, switch_name, option in (
        (log_length, "log_length", "--loglength"),
        (percent_quality, "percent_quality", "--percentqual"),
        (aligned_length, "aligned_length", "--alength"),
        (drop_outliers, "drop_outliers", "--drop_outliers"),
        (no_supplementary, "no_supplementary", "--no_supplementary"),
        (raw_table, "raw_table", "--raw"),
    ):
        if _cli.flag(switch, switch_name):
            arguments.append(option)
    for count, count_name, option in (
        (min_length, "min_length", "--minlength"),
        (max_length, "max_length", "--maxlength"),
        (downsample, "downsample", "--downsample"),
        (dpi, "dpi", "--dpi"),
    ):
        if count is not None:
            arguments += [option, str(_cli.integer(count, count_name, minimum=1))]
    if min_quality is not None:
        arguments += [
            "--minqual",
            format(_cli.number(min_quality, "min_quality", minimum=0.0), ".12g"),
        ]
    if show_n50 is not None:
        arguments.append("--N50" if _cli.flag(show_n50, "show_n50") else "--no-N50")
    for colour, colour_name, option in ((color, "color", "-c"), (colormap, "colormap", "-cm")):
        if colour is not None:
            if not isinstance(colour, str) or not _COLOR.match(colour):
                raise ConfigurationError(
                    f"{colour_name} must be a matplotlib name or a hex colour."
                )
            arguments += [option, colour]
    if title is not None:
        if not isinstance(title, str) or not _TITLE.match(title):
            raise ConfigurationError(
                "title must be text without control characters (at most 200 characters)."
            )
        arguments += ["--title", title]
    arguments += [_INPUT_FLAGS[kind], *(str(path) for path in inputs)]
    directory.mkdir()
    stats_file = directory / "NanoStats.txt"
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=directory,
        expected=[stats_file],
        explicit=nanoplot_path,
        timeout_seconds=timeout_seconds,
    )
    files = sorted(path.name for path in directory.iterdir() if path.is_file())
    plots: dict[str, list[str]] = {}
    for name in wanted:
        prefixes = PLOT_FILES[name]
        found = [
            str(directory / file)
            for file in files
            if file.startswith(prefixes) and not file.endswith((".log", ".tsv.gz", ".txt"))
        ]
        if found:
            plots[name] = found
    extra: dict[str, Any] = {
        "input_type": kind,
        "output_dir": str(directory),
        "plots": plots,
        "missing_plots": [name for name in wanted if name not in plots],
        "statistics": _parse_stats(stats_file),
        "stats_file": str(stats_file),
        "report_file": str(directory / "NanoPlot-report.html")
        if (directory / "NanoPlot-report.html").is_file()
        else None,
    }
    if raw_table:
        raw = directory / "NanoPlot-data.tsv.gz"
        columns, rows = _raw_table(raw)
        extra.update(raw_table_file=str(raw), raw_table_columns=columns, raw_table_rows=rows)
    record["outputs"] = [str(directory / name) for name in files]
    return _cli.finish(record, "nanoplot_plots", **extra)


__all__ = ["nanoplot_plots"]
_register_public_verbose(globals())
