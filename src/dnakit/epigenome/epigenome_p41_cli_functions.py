"""Fixed pairtools, MACS3 and deepTools operations missing from existing facades.

Commands use the bounded shared runner, preserve upstream defaults, and never
overwrite outputs. Inputs are explicit local files; no shell or arbitrary flags.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError, DNAKitError


def _paths(values: Sequence[str | Path], name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not 1 <= len(values) <= 100:
        raise ConfigurationError(f"{name} must list between 1 and 100 local files.")
    return [str(cli.existing_file(value, name)) for value in values]


def _preview(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        result = []
        for row in reader:
            result.append(dict(row))
            if len(result) == 20:
                break
        return result


def pairtools_tile_duplication_stats(pairs_path: str | Path, output_path: str | Path, *,
                                     tile_output_path: str | Path, executable_path: str | None = None,
                                     timeout_seconds: float = 3600, verbose: bool = True) -> dict[str, Any]:
    """pairtools stats --bytile-dups, returning normal stats and a by-tile duplicate table.

    The pairs header must include parent_readID (saved by pairtools dedup) and
    Illumina read IDs must retain tile information. The two outputs are new files.
    """
    source = cli.existing_file(pairs_path, "pairs_path")
    opener = gzip.open if source.suffix in (".gz", ".bgz") else open
    columns: list[str] = []
    with opener(source, "rt") as handle:
        for line in handle:
            if line.startswith("#columns:"):
                columns = line.partition(":")[2].split()
                break
            if not line.startswith("#"):
                break
    if "parent_readID" not in columns:
        raise ConfigurationError("By-tile stats require a parent_readID column from pairtools dedup.")
    target = cli.new_file(output_path, "output_path")
    tiles = cli.new_file(tile_output_path, "tile_output_path")
    if target == tiles:
        raise ConfigurationError("Stats and by-tile outputs must be different files.")
    result = cli.run("pairtools", ["stats", "--bytile-dups", "--output", str(target),
                      "--output-bytile-stats", str(tiles), str(source)], cwd=target.parent,
                      expected=[target, tiles], explicit=executable_path,
                      hint="Install dnakit[bio-pairtools] or pairtools.", timeout_seconds=timeout_seconds)
    if result["status"] == "completed":
        result["tile_duplicate_preview"] = _preview(tiles)
    return result


def macs3_bedgraph_cutoff_report(bedgraph: str | Path, output_path: str | Path, *,
                                min_length: int | None = None, max_gap: int | None = None,
                                cutoff_max: int | None = None, cutoff_steps: int | None = None,
                                executable_path: str | None = None, timeout_seconds: float = 3600,
                                verbose: bool = True) -> dict[str, Any]:
    """MACS3 bdgpeakcall --cutoff-analysis: peak counts/lengths versus score cutoff.

    Omitted min_length, max_gap, cutoff_max and cutoff_steps retain MACS3 defaults.
    Returns the real report path and up to 20 report rows.
    """
    source = cli.existing_file(bedgraph, "bedgraph")
    target = cli.new_file(output_path, "output_path")
    args = ["bdgpeakcall", "-i", str(source), "--cutoff-analysis", "--outdir", str(target.parent), "-o", target.name]
    for value, flag, name, minimum in ((min_length, "-l", "min_length", 1), (max_gap, "-g", "max_gap", 0),
                                       (cutoff_max, "--cutoff-analysis-max", "cutoff_max", 1),
                                       (cutoff_steps, "--cutoff-analysis-steps", "cutoff_steps", 1)):
        if value is not None:
            args += [flag, str(cli.integer(value, name, minimum=minimum))]
    result = cli.run("macs3", args, cwd=target.parent, expected=[target], explicit=executable_path,
                      hint="Install dnakit[bio-macs3] or MACS3.", timeout_seconds=timeout_seconds)
    if result["status"] == "completed":
        result["cutoff_preview"] = _preview(target)
    return result


def macs3_duplicate_statistics(alignment: str | Path, *, file_format: str = "AUTO",
                               genome_size: str | int = "hs", keep_duplicates: str | int = "auto",
                               executable_path: str | None = None, timeout_seconds: float = 3600,
                               verbose: bool = True) -> dict[str, Any]:
    """MACS3 filterdup --dry-run: allowed duplicates, reads before/after and redundancy.

    No filtered alignment file is written. The upstream default duplicate policy
    is auto; all and positive integers are also accepted.
    """
    source = cli.existing_file(alignment, "alignment")
    kind = cli.choice(file_format, ("AUTO", "BAM", "SAM", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "BOWTIE", "BAMPE", "BEDPE"), "file_format")
    if isinstance(genome_size, int) and not isinstance(genome_size, bool):
        genome = str(cli.integer(genome_size, "genome_size", minimum=1))
    else:
        genome = cli.choice(genome_size, ("hs", "mm", "ce", "dm"), "genome_size")
    if isinstance(keep_duplicates, int) and not isinstance(keep_duplicates, bool):
        duplicates = str(cli.integer(keep_duplicates, "keep_duplicates", minimum=1))
    else:
        duplicates = cli.choice(keep_duplicates, ("auto", "all"), "keep_duplicates")
    result = cli.run("macs3", ["filterdup", "-i", str(source), "-f", kind, "-g", genome,
                      "--keep-dup", duplicates, "--dry-run"], cwd=source.parent, expected=[],
                      explicit=executable_path, hint="Install dnakit[bio-macs3] or MACS3.",
                      timeout_seconds=timeout_seconds)
    if result["status"] != "completed":
        return result
    log = result["log_tail"]
    patterns = {"reads_before": r"total (?:tags|fragments/pairs) in alignment file:\s*(\d+)",
                "reads_after": r"tags after filtering in alignment file:\s*(\d+)",
                "max_duplicates": r"(?:at most |max_dup_tags based on binomal = )(\d+)",
                "redundancy_rate": r"Redundant rate of alignment file:\s*([\d.]+)"}
    statistics: dict[str, Any] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, log)
        if match:
            statistics[key] = float(match.group(1)) if key == "redundancy_rate" else int(match.group(1))
    if "reads_before" not in statistics:
        raise DNAKitError("MACS3 finished without reporting the number of input reads.")
    if duplicates == "all":
        statistics.update(reads_after=statistics["reads_before"], max_duplicates=None, redundancy_rate=0.0)
    elif "reads_after" not in statistics:
        raise DNAKitError("MACS3 finished without reporting duplicate-filtered read counts.")
    result.update(statistics=statistics, duplicate_policy=duplicates)
    return result


def _region_matrix(program: str, inputs: Sequence[str | Path], regions: Sequence[str | Path],
                   output_path: str | Path, *, raw_counts_path: str | Path | None,
                   labels: Sequence[str] | None, processors: int, executable_path: str | None,
                   timeout_seconds: float) -> dict[str, Any]:
    files = _paths(inputs, "inputs")
    beds = _paths(regions, "regions")
    count = cli.integer(processors, "processors", minimum=1, maximum=256)
    target = cli.new_file(output_path, "output_path")
    if target.suffix != ".npz":
        raise ConfigurationError("output_path must end in .npz (deepTools NumPy archive).")
    raw = None if raw_counts_path is None else cli.new_file(raw_counts_path, "raw_counts_path")
    if raw == target:
        raise ConfigurationError("Raw counts and NumPy archive must use different paths.")
    args = ["BED-file", "--bamfiles" if program == "multiBamSummary" else "--bwfiles", *files,
            "--BED", *beds, "-o", str(target), "--numberOfProcessors", str(count)]
    if labels is not None:
        if isinstance(labels, (str, bytes)) or len(labels) != len(files):
            raise ConfigurationError("labels must contain one label per input file.")
        args += ["--labels", *[cli.text(label, "label") for label in labels]]
    if raw is not None:
        args += ["--outRawCounts", str(raw)]
    result = cli.run(program, args, cwd=target.parent, expected=[target] + ([] if raw is None else [raw]),
                      explicit=executable_path, hint="Install dnakit[bio-deeptools] or deepTools.",
                      timeout_seconds=timeout_seconds)
    if result["status"] == "completed":
        np = require("numpy", "bio-deeptools")
        with np.load(target, allow_pickle=False) as data:
            matrix = data["matrix"]
            result["matrix_shape"] = list(matrix.shape)
            result["labels"] = [bytes(value).decode("utf-8") if getattr(value, "ndim", 0) == 1 and getattr(value.dtype, "kind", None) in ("u", "i") else (value.decode("utf-8") if isinstance(value, bytes) else str(value)) for value in data["labels"]]
            result["matrix_preview"] = [[float(v) if math.isfinite(v) else None for v in row]
                                        for row in matrix[:20, :100]]
    return result


def bam_region_count_matrix(bams: Sequence[str | Path], regions: Sequence[str | Path],
                             output_path: str | Path, *, raw_counts_path: str | Path | None = None,
                             labels: Sequence[str] | None = None, processors: int = 1,
                             executable_path: str | None = None, timeout_seconds: float = 3600,
                             verbose: bool = True) -> dict[str, Any]:
    """multiBamSummary BED-file: read-count matrix over caller-supplied BED/GTF regions."""
    return _region_matrix("multiBamSummary", bams, regions, output_path, raw_counts_path=raw_counts_path,
                          labels=labels, processors=processors, executable_path=executable_path,
                          timeout_seconds=timeout_seconds)


def bigwig_region_score_matrix(bigwigs: Sequence[str | Path], regions: Sequence[str | Path],
                                output_path: str | Path, *, raw_counts_path: str | Path | None = None,
                                labels: Sequence[str] | None = None, processors: int = 1,
                                executable_path: str | None = None, timeout_seconds: float = 3600,
                                verbose: bool = True) -> dict[str, Any]:
    """multiBigwigSummary BED-file: mean bigWig scores over BED/GTF regions."""
    return _region_matrix("multiBigwigSummary", bigwigs, regions, output_path, raw_counts_path=raw_counts_path,
                          labels=labels, processors=processors, executable_path=executable_path,
                          timeout_seconds=timeout_seconds)


__all__ = ["pairtools_tile_duplication_stats", "macs3_bedgraph_cutoff_report", "macs3_duplicate_statistics",
           "bam_region_count_matrix", "bigwig_region_score_matrix"]
_register_public_verbose(globals(), tuple(__all__))
