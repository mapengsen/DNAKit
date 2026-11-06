"""More deepTools command-line tools as typed DNAKit facades (second batch).

Covers computeMatrix (scale-regions and reference-point), multiBamSummary bins, multiBigwigSummary bins,
bigwigAverage, plotHeatmap, plotProfile, plotCorrelation and plotPCA. Same contract as
``deeptools_cli_facades``: the user installs deepTools, each function has an explicit typed signature
mapped to fixed flags, no shell, outputs are new explicit paths, results are JSON-serialisable.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.epigenome.deeptools_cli_facades import _existing, _new, _run
from dnakit.exceptions import ConfigurationError


def _many(paths: Sequence[str], name: str) -> list[str]:
    if not paths:
        raise ConfigurationError(f"{name} must contain at least one file.")
    return [str(_existing(path, name)) for path in paths]


def compute_matrix_scale_regions(
    score_files: Sequence[str],
    regions: Sequence[str],
    output_matrix: str,
    *,
    region_body_length: int = 1000,
    before_region_length: int = 500,
    after_region_length: int = 500,
    bin_size: int = 10,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a scores-per-region matrix with scaled regions (deepTools computeMatrix scale-regions)."""
    scores = _many(score_files, "score_files")
    beds = _many(regions, "regions")
    out = _new(output_matrix, "output_matrix")
    args = [
        "scale-regions",
        "-S",
        *scores,
        "-R",
        *beds,
        "-o",
        str(out),
        "-m",
        str(region_body_length),
        "-b",
        str(before_region_length),
        "-a",
        str(after_region_length),
        "--binSize",
        str(bin_size),
        "-p",
        str(number_of_processors),
    ]
    return _run("computeMatrix", args, [out], deeptools_path)


def compute_matrix_reference_point(
    score_files: Sequence[str],
    regions: Sequence[str],
    output_matrix: str,
    *,
    reference_point: Literal["TSS", "TES", "center"] = "TSS",
    before_region_length: int = 500,
    after_region_length: int = 500,
    bin_size: int = 10,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a scores-per-region matrix around a reference point (deepTools computeMatrix reference-point)."""
    scores = _many(score_files, "score_files")
    beds = _many(regions, "regions")
    out = _new(output_matrix, "output_matrix")
    args = [
        "reference-point",
        "--referencePoint",
        reference_point,
        "-S",
        *scores,
        "-R",
        *beds,
        "-o",
        str(out),
        "-b",
        str(before_region_length),
        "-a",
        str(after_region_length),
        "--binSize",
        str(bin_size),
        "-p",
        str(number_of_processors),
    ]
    return _run("computeMatrix", args, [out], deeptools_path)


def multi_bam_summary(
    bams: Sequence[str],
    output: str,
    *,
    bin_size: int = 10000,
    distance_between_bins: int = 0,
    labels: Sequence[str] | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise read counts of several BAM files in genomic bins (deepTools multiBamSummary bins)."""
    paths = _many(bams, "bams")
    out = _new(output, "output")
    args = [
        "bins",
        "-b",
        *paths,
        "-o",
        str(out),
        "--binSize",
        str(bin_size),
        "-n",
        str(distance_between_bins),
        "-p",
        str(number_of_processors),
    ]
    if labels is not None:
        if len(labels) != len(paths):
            raise ConfigurationError("labels must have one entry per BAM file.")
        args += ["--labels", *list(labels)]
    return _run("multiBamSummary", args, [out], deeptools_path)


def multi_bigwig_summary(
    bigwigs: Sequence[str],
    output: str,
    *,
    bin_size: int = 10000,
    labels: Sequence[str] | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise several bigWig tracks in genomic bins (deepTools multiBigwigSummary bins)."""
    paths = _many(bigwigs, "bigwigs")
    out = _new(output, "output")
    args = [
        "bins",
        "-b",
        *paths,
        "-o",
        str(out),
        "--binSize",
        str(bin_size),
        "-p",
        str(number_of_processors),
    ]
    if labels is not None:
        if len(labels) != len(paths):
            raise ConfigurationError("labels must have one entry per bigWig file.")
        args += ["--labels", *list(labels)]
    return _run("multiBigwigSummary", args, [out], deeptools_path)


def bigwig_average(
    bigwigs: Sequence[str],
    output: str,
    *,
    bin_size: int = 50,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average several bigWig tracks bin by bin into one bigWig (deepTools bigwigAverage)."""
    paths = _many(bigwigs, "bigwigs")
    out = _new(output, "output")
    args = [
        "-b",
        *paths,
        "-o",
        str(out),
        "--binSize",
        str(bin_size),
        "-p",
        str(number_of_processors),
    ]
    return _run("bigwigAverage", args, [out], deeptools_path)


def plot_heatmap(
    matrix: str,
    output_figure: str,
    *,
    kmeans: int | None = None,
    dpi: int = 200,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a heatmap from a computeMatrix output (deepTools plotHeatmap)."""
    matrix_path = _existing(matrix, "matrix")
    out = _new(output_figure, "output_figure")
    args = ["-m", str(matrix_path), "-out", str(out), "--dpi", str(dpi)]
    if kmeans is not None:
        if kmeans < 1:
            raise ConfigurationError("kmeans must be a positive integer.")
        args += ["--kmeans", str(kmeans)]
    return _run("plotHeatmap", args, [out], deeptools_path)


def plot_profile(
    matrix: str,
    output_figure: str,
    *,
    kmeans: int | None = None,
    dpi: int = 200,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot average profiles from a computeMatrix output (deepTools plotProfile)."""
    matrix_path = _existing(matrix, "matrix")
    out = _new(output_figure, "output_figure")
    args = ["-m", str(matrix_path), "-out", str(out), "--dpi", str(dpi)]
    if kmeans is not None:
        if kmeans < 1:
            raise ConfigurationError("kmeans must be a positive integer.")
        args += ["--kmeans", str(kmeans)]
    return _run("plotProfile", args, [out], deeptools_path)


def plot_correlation(
    summary_matrix: str,
    output_figure: str,
    *,
    method: Literal["spearman", "pearson"] = "spearman",
    plot_type: Literal["heatmap", "scatterplot"] = "heatmap",
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the correlation between samples from a multiBamSummary or multiBigwigSummary matrix (plotCorrelation)."""
    data = _existing(summary_matrix, "summary_matrix")
    out = _new(output_figure, "output_figure")
    args = ["-in", str(data), "-c", method, "-p", plot_type, "-o", str(out)]
    return _run("plotCorrelation", args, [out], deeptools_path)


def plot_pca(
    summary_matrix: str,
    output_figure: str,
    *,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot a principal component analysis of samples from a summary matrix (deepTools plotPCA)."""
    data = _existing(summary_matrix, "summary_matrix")
    out = _new(output_figure, "output_figure")
    args = ["-in", str(data), "-o", str(out)]
    return _run("plotPCA", args, [out], deeptools_path)


__all__ = [
    "bigwig_average",
    "compute_matrix_reference_point",
    "compute_matrix_scale_regions",
    "multi_bam_summary",
    "multi_bigwig_summary",
    "plot_correlation",
    "plot_heatmap",
    "plot_pca",
    "plot_profile",
]
_register_public_verbose(globals())
