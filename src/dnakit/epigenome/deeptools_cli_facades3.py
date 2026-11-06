"""deepTools plotEnrichment and estimateReadFiltering as typed DNAKit facades (third batch)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.epigenome.deeptools_cli_facades import _existing, _new, _run
from dnakit.exceptions import ConfigurationError


def _many(paths: Sequence[str], name: str) -> list[str]:
    if not paths:
        raise ConfigurationError(f"{name} must contain at least one file.")
    return [str(_existing(path, name)) for path in paths]


def plot_enrichment(
    bams: Sequence[str],
    regions: Sequence[str],
    output_figure: str,
    *,
    labels: Sequence[str] | None = None,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the fraction of reads falling in each set of regions (deepTools plotEnrichment)."""
    bam_paths = _many(bams, "bams")
    bed_paths = _many(regions, "regions")
    out = _new(output_figure, "output_figure")
    args = ["-b", *bam_paths, "--BED", *bed_paths, "-o", str(out)]
    if labels is not None:
        if len(labels) != len(bam_paths):
            raise ConfigurationError("labels must have one entry per BAM file.")
        args += ["--labels", *list(labels)]
    return _run("plotEnrichment", args, [out], deeptools_path)


def estimate_read_filtering(
    bams: Sequence[str],
    output: str,
    *,
    bin_size: int = 1000,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate how much of each BAM file is removed by filtering (deepTools estimateReadFiltering)."""
    paths = _many(bams, "bams")
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
    return _run("estimateReadFiltering", args, [out], deeptools_path)


__all__ = ["estimate_read_filtering", "plot_enrichment"]
_register_public_verbose(globals())
