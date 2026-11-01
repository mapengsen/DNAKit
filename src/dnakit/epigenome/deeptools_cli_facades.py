"""deepTools command-line tools as typed DNAKit facades (bamCoverage, bamCompare, bigwigCompare,
bamPEFragmentSize, alignmentSieve, estimateScaleFactor, plotCoverage, plotFingerprint).

The user installs deepTools (``pip install deeptools`` or conda). Each function runs one tool with an
explicit, typed keyword list that maps to fixed flags; no other option can be passed, and no shell is used.
Outputs are written to explicit new paths. Results are JSON-serialisable dictionaries with status,
backend, command, outputs and a bounded log tail. Verified against deepTools 4.0.0.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Literal, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_TOOL_HINT = "Install deepTools 4.0.0 or later (python -m pip install deeptools, or conda install -c bioconda deeptools)."
_MAX_TAIL = 4000


def _program(tool: str, deeptools_path: str | None) -> str:
    if deeptools_path is not None:
        candidate = Path(deeptools_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"deeptools_path does not exist: {candidate}", hint=_TOOL_HINT)
        return str(candidate)
    found = shutil.which(tool)
    if found is None:
        raise BackendUnavailableError(f"{tool} was not found on PATH.", hint=_TOOL_HINT)
    return found


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _run(tool: str, arguments: list[str], outputs: Sequence[Path], deeptools_path: str | None) -> dict[str, Any]:
    program = _program(tool, deeptools_path)
    cwd = outputs[0].parent
    result = execute_bounded_command(program, arguments, backend_id=f"deeptools:{tool}", cwd=cwd, verbose=False)
    status = "completed" if result.return_code == 0 and all(path.is_file() for path in outputs) else "failed"
    return {
        "status": status,
        "backend": "deeptools",
        "tool": tool,
        "command": [tool, *arguments],
        "return_code": result.return_code,
        "outputs": [str(path) for path in outputs if path.is_file()],
        "log_tail": (result.output or "")[-_MAX_TAIL:],
    }


_OPERATIONS = Literal["log2", "ratio", "subtract", "add", "mean", "reciprocal_ratio", "first", "second"]
_NORMALIZATION = Literal["RPGC", "CPM", "BPM", "RPKM", "None"]


def bam_coverage(
    bam: str,
    output: str,
    *,
    output_format: Literal["bigwig", "bedgraph"] = "bigwig",
    bin_size: int = 50,
    normalize_using: _NORMALIZATION = "None",
    effective_genome_size: int | None = None,
    extend_reads: int | None = None,
    scale_factor: float | None = None,
    region: str | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a coverage track from one BAM file (deepTools bamCoverage)."""
    bam_path = _existing(bam, "bam")
    out = _new(output, "output")
    if bin_size < 1:
        raise ConfigurationError("bin_size must be a positive integer.")
    args = ["-b", str(bam_path), "-o", str(out), "-of", output_format, "--binSize", str(bin_size),
            "--normalizeUsing", normalize_using, "-p", str(number_of_processors)]
    if effective_genome_size is not None:
        args += ["--effectiveGenomeSize", str(effective_genome_size)]
    if extend_reads is not None:
        args += ["--extendReads", str(extend_reads)]
    if scale_factor is not None:
        args += ["--scaleFactor", str(scale_factor)]
    if region is not None:
        args += ["--region", region]
    return _run("bamCoverage", args, [out], deeptools_path)


def bam_compare(
    bam1: str,
    bam2: str,
    output: str,
    *,
    operation: _OPERATIONS = "log2",
    output_format: Literal["bigwig", "bedgraph"] = "bigwig",
    bin_size: int = 50,
    pseudocount: float = 1.0,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two BAM files bin by bin (deepTools bamCompare)."""
    b1, b2 = _existing(bam1, "bam1"), _existing(bam2, "bam2")
    out = _new(output, "output")
    args = ["-b1", str(b1), "-b2", str(b2), "-o", str(out), "-of", output_format, "--binSize", str(bin_size),
            "--operation", operation, "--pseudocount", str(pseudocount), "-p", str(number_of_processors)]
    return _run("bamCompare", args, [out], deeptools_path)


def bigwig_compare(
    bigwig1: str,
    bigwig2: str,
    output: str,
    *,
    operation: _OPERATIONS = "log2",
    bin_size: int = 50,
    pseudocount: float = 1.0,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two bigWig tracks bin by bin (deepTools bigwigCompare)."""
    b1, b2 = _existing(bigwig1, "bigwig1"), _existing(bigwig2, "bigwig2")
    out = _new(output, "output")
    args = ["-b1", str(b1), "-b2", str(b2), "-o", str(out), "--binSize", str(bin_size),
            "--operation", operation, "--pseudocount", str(pseudocount), "-p", str(number_of_processors)]
    return _run("bigwigCompare", args, [out], deeptools_path)


def bam_pe_fragment_size(
    bams: Sequence[str],
    histogram_output: str,
    *,
    table_output: str | None = None,
    max_fragment_length: int | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate paired-end fragment length distributions (deepTools bamPEFragmentSize)."""
    if not bams:
        raise ConfigurationError("bams must contain at least one BAM file.")
    paths = [str(_existing(bam, "bam")) for bam in bams]
    hist = _new(histogram_output, "histogram_output")
    args = ["-b", *paths, "-hist", str(hist), "-p", str(number_of_processors)]
    outputs = [hist]
    if table_output is not None:
        table = _new(table_output, "table_output")
        args += ["--table", str(table)]
        outputs.append(table)
    if max_fragment_length is not None:
        args += ["--maxFragmentLength", str(max_fragment_length)]
    return _run("bamPEFragmentSize", args, outputs, deeptools_path)


def alignment_sieve(
    bam: str,
    output: str,
    *,
    shift: Sequence[int] | None = None,
    bed_output: bool = False,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter and optionally shift alignments (deepTools alignmentSieve); output is BAM or BEDPE."""
    bam_path = _existing(bam, "bam")
    out = _new(output, "output")
    args = ["--bam", str(bam_path), "-o", str(out), "-p", str(number_of_processors)]
    if shift is not None:
        if len(shift) != 4 or not all(isinstance(value, int) for value in shift):
            raise ConfigurationError("shift must be four integers: forward start, reverse start, forward end, reverse end.")
        args += ["--shift", *[str(value) for value in shift]]
    if bed_output:
        args.append("--BED")
    return _run("alignmentSieve", args, [out], deeptools_path)


def estimate_scale_factor(
    bams: Sequence[str],
    *,
    sample_window_length: int = 1000,
    number_of_samples: int = 1_000_000,
    normalization_length: int = 1000,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate scale factors from sampled read counts (deepTools estimateScaleFactor). Prints the factors."""
    if len(bams) < 2:
        raise ConfigurationError("estimateScaleFactor needs at least two BAM files (scale factors are defined per pair).")
    paths = [str(_existing(bam, "bam")) for bam in bams]
    args = ["-b", *paths, "-l", str(sample_window_length), "-n", str(number_of_samples),
            "-nl", str(normalization_length), "-p", str(number_of_processors)]
    program = _program("estimateScaleFactor", deeptools_path)
    result = execute_bounded_command(program, args, backend_id="deeptools:estimateScaleFactor",
                                     cwd=Path(paths[0]).parent, verbose=False)
    return {
        "status": "completed" if result.return_code == 0 else "failed",
        "backend": "deeptools",
        "tool": "estimateScaleFactor",
        "command": ["estimateScaleFactor", *args],
        "return_code": result.return_code,
        "stdout_tail": (result.output or "")[-_MAX_TAIL:],
    }


def plot_coverage(
    bams: Sequence[str],
    plot_output: str,
    *,
    labels: Sequence[str] | None = None,
    number_of_samples: int = 1_000_000,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the coverage distribution of one or more BAM files (deepTools plotCoverage)."""
    if not bams:
        raise ConfigurationError("bams must contain at least one BAM file.")
    paths = [str(_existing(bam, "bam")) for bam in bams]
    out = _new(plot_output, "plot_output")
    args = ["-b", *paths, "-o", str(out), "-n", str(number_of_samples)]
    if labels is not None:
        if len(labels) != len(paths):
            raise ConfigurationError("labels must have one entry per BAM file.")
        args += ["--labels", *list(labels)]
    return _run("plotCoverage", args, [out], deeptools_path)


def plot_fingerprint(
    bams: Sequence[str],
    plot_output: str,
    *,
    extend_reads: int | None = None,
    min_mapping_quality: int | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the cumulative read-coverage fingerprint of ChIP-like BAM files (deepTools plotFingerprint)."""
    if not bams:
        raise ConfigurationError("bams must contain at least one BAM file.")
    paths = [str(_existing(bam, "bam")) for bam in bams]
    out = _new(plot_output, "plot_output")
    args = ["-b", *paths, "-plot", str(out), "-p", str(number_of_processors)]
    if extend_reads is not None:
        args += ["--extendReads", str(extend_reads)]
    if min_mapping_quality is not None:
        args += ["--minMappingQuality", str(min_mapping_quality)]
    return _run("plotFingerprint", args, [out], deeptools_path)


__all__ = [
    "alignment_sieve",
    "bam_compare",
    "bam_coverage",
    "bam_pe_fragment_size",
    "bigwig_compare",
    "estimate_scale_factor",
    "plot_coverage",
    "plot_fingerprint",
]
_register_public_verbose(globals())
