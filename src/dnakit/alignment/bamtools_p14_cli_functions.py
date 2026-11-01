"""BamTools command-line facades (tested with BamTools 2.5.3): convert BAM to JSON/YAML and other text formats, count,
filter, random-sample, resolve paired-end flags, revert, split and index (BamTools ``.bti``) BAM files.

The user installs BamTools (``conda install -c bioconda bamtools``, or from https://github.com/pezmaster31/bamtools)
and puts it on PATH, or passes ``bamtools_path``. Each public function runs one fixed ``bamtools`` sub-command with an
argument list built from typed, validated parameters: no shell, no generic runner, outputs are new files. Only the
options that the caller sets are passed, so every other option keeps BamTools' own default. Results are
JSON-serialisable dictionaries (status, backend, version, command, return code, outputs, log tail).
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = (
    "Install BamTools 2.5 or later (for example `conda install -c bioconda bamtools`, or build it from "
    "https://github.com/pezmaster31/bamtools), then put it on PATH or pass bamtools_path."
)
_CONVERT_FORMATS = ("json", "yaml", "bed", "fasta", "fastq", "pileup", "sam")
_SPLIT_MODES = ("mapped", "paired", "reference", "tag")
_RESOLVE_MODES = ("markPairs", "twoPass", "makeStats")
_PATTERN = re.compile(r"^(?:>=|<=|>|<|!)?-?\d+$")
_TAG_FILTER = re.compile(r"^[A-Za-z][A-Za-z0-9]:(?:>=|<=|>|<|!)?\S+$")
_TAG_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]$")
_COUNT = re.compile(r"^\s*(\d+)\s*$")


def _bamtools(
    sub: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    bamtools_path: str | None,
    stdout_file: Path | None = None,
) -> dict[str, Any]:
    return _p14_cli.run(
        "bamtools",
        [sub, *arguments],
        cwd=cwd,
        expected=expected,
        explicit=bamtools_path,
        hint=_HINT,
        stdout_file=stdout_file,
    )


def _inputs(bams: Sequence[str]) -> list[str]:
    if isinstance(bams, str) or not bams:
        raise ConfigurationError("bams must be a non-empty list of BAM file paths.")
    arguments: list[str] = []
    for item in bams:
        arguments += ["-in", str(_gap_cli.existing_file(item, "bams"))]
    return arguments


def _region(region: str | None) -> list[str]:
    return [] if region is None else ["-region", _gap_cli.text(region, "region")]


def _pattern(value: int | str, name: str) -> str:
    """A number, optionally with a comparison prefix: ``50``, ``'>=50'``, ``'<100'``, ``'!0'``."""
    text = str(value)
    if isinstance(value, bool) or not _PATTERN.fullmatch(text):
        raise ConfigurationError(f"{name} must be a number, optionally prefixed by >, >=, <, <= or !.")
    return text


def bamtools_convert_alignments(
    bams: Sequence[str],
    output: str,
    *,
    output_format: str = "json",
    region: str | None = None,
    no_header: bool = False,
    reference_fasta: str | None = None,
    print_mapping_quality: bool = False,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert BAM alignments to JSON, YAML or another text format (``bamtools convert``).

    ``output_format`` is 'json' (one JSON object per alignment: name, flag, reference, position, mapping quality,
    CIGAR, bases, qualities, tags), 'yaml', 'bed', 'fasta', 'fastq', 'pileup' or 'sam'. ``region`` restricts the
    alignments (for example ``'chr1:1000..2000'``; an index next to the BAM is used when it exists). ``no_header``
    omits the SAM header (sam only); ``reference_fasta`` and ``print_mapping_quality`` are the pileup options.
    Several BAM files are merged. Writes a new file.
    """
    arguments = _inputs(bams)
    output_path = _gap_cli.new_file(output, "output")
    arguments += ["-format", _gap_cli.choice(output_format, _CONVERT_FORMATS, "output_format")]
    arguments += _region(region)
    if _p14_cli.flag(no_header, "no_header"):
        arguments.append("-noheader")
    if reference_fasta is not None:
        arguments += ["-fasta", str(_gap_cli.existing_file(reference_fasta, "reference_fasta"))]
    if _p14_cli.flag(print_mapping_quality, "print_mapping_quality"):
        arguments.append("-mapqual")
    arguments += ["-out", str(output_path)]
    return _bamtools("convert", arguments, cwd=output_path.parent, expected=[output_path], bamtools_path=bamtools_path)


def bamtools_count_alignments(
    bams: Sequence[str],
    *,
    region: str | None = None,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the alignments in one or more BAM files (``bamtools count``), optionally within a region.

    The BAM files are merged before counting. ``count`` in the result is the number printed by BamTools.
    """
    arguments = _inputs(bams) + _region(region)
    with tempfile.TemporaryDirectory(prefix="dnakit_bamtools_") as folder:
        work = Path(folder)
        result = _bamtools(
            "count", arguments, cwd=work, expected=[], bamtools_path=bamtools_path, stdout_file=work / "count.txt"
        )
        match = _COUNT.match((work / "count.txt").read_text(encoding="utf-8").splitlines()[0]) if result["status"] == "completed" else None
    result["outputs"] = []
    result["count"] = int(match.group(1)) if match is not None else None
    if result["status"] == "completed" and match is None:
        result["status"] = "failed"
    return result


def bamtools_filter_alignments(
    bams: Sequence[str],
    output: str,
    *,
    region: str | None = None,
    script: str | None = None,
    alignment_flag: int | None = None,
    insert_size: int | str | None = None,
    length: int | str | None = None,
    map_quality: int | str | None = None,
    name: str | None = None,
    query_bases: str | None = None,
    tag: str | None = None,
    is_duplicate: bool | None = None,
    is_failed_qc: bool | None = None,
    is_first_mate: bool | None = None,
    is_mapped: bool | None = None,
    is_mate_mapped: bool | None = None,
    is_mate_reverse_strand: bool | None = None,
    is_paired: bool | None = None,
    is_primary_alignment: bool | None = None,
    is_proper_pair: bool | None = None,
    is_reverse_strand: bool | None = None,
    is_second_mate: bool | None = None,
    is_singleton: bool | None = None,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the BAM alignments that match all given criteria (``bamtools filter``).

    ``alignment_flag`` keeps alignments with exactly that flag; ``insert_size``, ``length`` and ``map_quality`` are a
    number with an optional comparison prefix (``50``, ``'>=50'``, ``'<100'``, ``'!0'``); ``name`` matches the read
    name (``*`` wildcard); ``query_bases`` matches a base motif; ``tag`` is ``'NM:<4'``. The ``is_*`` parameters keep
    (True) or drop (False) alignments with that flag, for example ``is_mapped=True`` or ``is_reverse_strand=False``.
    ``script`` is a BamTools filter script (JSON) and cannot be combined with the other criteria; ``region``
    restricts the input. Writes a new BAM file.
    """
    arguments = _inputs(bams)
    output_path = _gap_cli.new_file(output, "output")
    arguments += _region(region)
    criteria: list[str] = []
    if alignment_flag is not None:
        criteria += ["-alignmentFlag", str(_gap_cli.integer(alignment_flag, "alignment_flag", minimum=0, maximum=65535))]
    for flag_name, value in (("insertSize", insert_size), ("length", length), ("mapQuality", map_quality)):
        if value is not None:
            criteria += [f"-{flag_name}", _pattern(value, flag_name)]
    if name is not None:
        criteria += ["-name", _gap_cli.text(name, "name")]
    if query_bases is not None:
        criteria += ["-queryBases", _gap_cli.text(query_bases, "query_bases")]
    if tag is not None:
        if not _TAG_FILTER.fullmatch(tag):
            raise ConfigurationError("tag must look like 'NM:<4' or 'RG:group1'.")
        criteria += ["-tag", tag]
    for flag_name, value in (
        ("isDuplicate", is_duplicate), ("isFailedQC", is_failed_qc), ("isFirstMate", is_first_mate),
        ("isMapped", is_mapped), ("isMateMapped", is_mate_mapped), ("isMateReverseStrand", is_mate_reverse_strand),
        ("isPaired", is_paired), ("isPrimaryAlignment", is_primary_alignment), ("isProperPair", is_proper_pair),
        ("isReverseStrand", is_reverse_strand), ("isSecondMate", is_second_mate), ("isSingleton", is_singleton),
    ):
        if value is not None:
            criteria += [f"-{flag_name}", "true" if _p14_cli.flag(value, flag_name) else "false"]
    if script is not None:
        if criteria:
            raise ConfigurationError("script cannot be combined with the other filter criteria.")
        arguments += ["-script", str(_gap_cli.existing_file(script, "script"))]
    elif not criteria:
        raise ConfigurationError("Give at least one filter criterion or a script.")
    arguments += criteria + ["-out", str(output_path)]
    return _bamtools("filter", arguments, cwd=output_path.parent, expected=[output_path], bamtools_path=bamtools_path)


def bamtools_index_bam_bti(
    bam: str,
    *,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the BamTools-specific index ``<bam>.bti`` of a coordinate-sorted BAM file (``bamtools index -bti``).

    The standard BAI index is not made. The ``.bti`` file is written next to the BAM and must not exist yet.
    """
    source = _gap_cli.existing_file(bam, "bam")
    target = Path(str(source) + ".bti")
    if target.exists():
        raise ConfigurationError(f"The index already exists: {target}")
    return _bamtools(
        "index", ["-in", str(source), "-bti"], cwd=source.parent, expected=[target], bamtools_path=bamtools_path
    )


def bamtools_sample_alignments(
    bams: Sequence[str],
    output: str,
    *,
    count: int | None = None,
    seed: int | None = None,
    region: str | None = None,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a random subset of the alignments of one or more BAM files (``bamtools random``).

    ``count`` is the number of alignments to draw (BamTools default 10000; no duplicate check is made, so a draw
    can contain the same alignment twice) and ``seed`` makes the draw repeatable. ``region`` limits the sampled
    alignments. Writes a new BAM file.
    """
    arguments = _inputs(bams)
    output_path = _gap_cli.new_file(output, "output")
    arguments += _region(region)
    if count is not None:
        arguments += ["-n", str(_gap_cli.integer(count, "count", minimum=1))]
    if seed is not None:
        arguments += ["-seed", str(_gap_cli.integer(seed, "seed", minimum=0, maximum=2**32 - 1))]
    arguments += ["-out", str(output_path)]
    return _bamtools("random", arguments, cwd=output_path.parent, expected=[output_path], bamtools_path=bamtools_path)


def bamtools_resolve_pairs(
    bam: str,
    stats: str,
    *,
    output: str | None = None,
    mode: str = "twoPass",
    min_mapping_quality: int | None = None,
    confidence_interval: float | None = None,
    unused_model_threshold: float | None = None,
    force: bool = False,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mark the proper-pair flag of paired-end alignments from the fragment-length model (``bamtools resolve``).

    ``mode='twoPass'`` (the default) learns the fragment-length statistics from the BAM, writes them to ``stats`` (a
    new file) and then marks the alignments into ``output`` (a new BAM). ``mode='makeStats'`` only writes the
    statistics file ``stats`` (no ``output``). ``mode='markPairs'`` marks ``bam`` into ``output`` from an existing
    ``stats`` file made earlier by makeStats or twoPass (BamTools also needs its ``<stats>.uniq_names.txt`` file next
    to it). ``min_mapping_quality`` is the minimum mapping quality; ``confidence_interval`` and
    ``unused_model_threshold`` tune the model (makeStats/twoPass); ``force`` marks all read groups even when their
    orientation model is ambiguous (markPairs/twoPass).
    """
    source = _gap_cli.existing_file(bam, "bam")
    chosen = _gap_cli.choice(mode, _RESOLVE_MODES, "mode")
    arguments = [f"-{chosen}", "-in", str(source)]
    expected: list[Path] = []
    if chosen == "markPairs":
        arguments += ["-stats", str(_gap_cli.existing_file(stats, "stats"))]
    else:
        stats_path = _gap_cli.new_file(stats, "stats")
        arguments += ["-stats", str(stats_path)]
        expected.append(stats_path)
    if chosen == "makeStats":
        if output is not None:
            raise ConfigurationError("output is not used with mode='makeStats'.")
        if force:
            raise ConfigurationError("force applies to markPairs and twoPass only.")
        cwd = Path(expected[0]).parent
    else:
        if output is None:
            raise ConfigurationError("output is required with mode='markPairs' or 'twoPass'.")
        output_path = _gap_cli.new_file(output, "output")
        arguments += ["-out", str(output_path)]
        expected.append(output_path)
        cwd = output_path.parent
        if chosen == "twoPass" and expected[0].parent != cwd:
            raise ConfigurationError("stats and output must be in the same folder.")
    if min_mapping_quality is not None:
        arguments += ["-minMQ", str(_gap_cli.integer(min_mapping_quality, "min_mapping_quality", minimum=0, maximum=255))]
    if chosen != "markPairs":
        if confidence_interval is not None:
            arguments += ["-ci", str(_gap_cli.number(confidence_interval, "confidence_interval", minimum=0.0, maximum=1.0))]
        if unused_model_threshold is not None:
            arguments += ["-umt", str(_gap_cli.number(unused_model_threshold, "unused_model_threshold", minimum=0.0, maximum=1.0))]
    elif confidence_interval is not None or unused_model_threshold is not None:
        raise ConfigurationError("confidence_interval and unused_model_threshold do not apply to mode='markPairs'.")
    if _p14_cli.flag(force, "force"):
        arguments.append("-force")
    result = _bamtools("resolve", arguments, cwd=cwd, expected=expected, bamtools_path=bamtools_path)
    result["mode"] = chosen
    return result


def bamtools_revert_alignments(
    bams: Sequence[str],
    output: str,
    *,
    keep_duplicate_marks: bool = False,
    keep_qualities: bool = False,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove duplicate marks and restore the original base qualities (OQ tag) of BAM alignments (``bamtools revert``).

    By default both are done; ``keep_duplicate_marks=True`` keeps the duplicate flag and ``keep_qualities=True``
    keeps the current base qualities instead of replacing them with the OQ tag. Writes a new BAM file.
    """
    arguments = _inputs(bams)
    output_path = _gap_cli.new_file(output, "output")
    if _p14_cli.flag(keep_duplicate_marks, "keep_duplicate_marks"):
        arguments.append("-keepDuplicate")
    if _p14_cli.flag(keep_qualities, "keep_qualities"):
        arguments.append("-keepQualities")
    arguments += ["-out", str(output_path)]
    return _bamtools("revert", arguments, cwd=output_path.parent, expected=[output_path], bamtools_path=bamtools_path)


def bamtools_split_alignments(
    bam: str,
    output_stub: str,
    *,
    by: str = "mapped",
    tag: str | None = None,
    reference_prefix: str | None = None,
    tag_prefix: str | None = None,
    tag_list_delimiter: str | None = None,
    bamtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a BAM file into one BAM file per value of a property (``bamtools split``).

    ``by`` is 'mapped' (``<stub>.MAPPED.bam`` and ``<stub>.UNMAPPED.bam``), 'paired' (``<stub>.PAIRED_END.bam`` and
    ``<stub>.SINGLE_END.bam``), 'reference' (``<stub>.REF_<name>.bam``, one per reference; ``reference_prefix``
    replaces ``REF_``) or 'tag' (``<stub>.TAG_<tag>_<value>.bam`` for every value of the tag named by ``tag``;
    ``tag_prefix`` replaces ``TAG_`` and ``tag_list_delimiter`` joins the values of list-type tags). BamTools only
    writes the files for values that occur. ``output_stub`` is the path prefix of the outputs; existing files with
    that prefix are refused. ``outputs`` lists the files written.
    """
    source = _gap_cli.existing_file(bam, "bam")
    stub = Path(output_stub).expanduser().resolve()
    if not stub.parent.is_dir():
        raise ConfigurationError(f"The folder for output_stub does not exist: {stub.parent}")
    if any(stub.parent.glob(f"{stub.name}.*")):
        raise ConfigurationError(f"Files with the prefix {stub.name}. already exist in {stub.parent}")
    mode = _gap_cli.choice(by, _SPLIT_MODES, "by")
    arguments = ["-in", str(source), "-stub", str(stub)]
    if mode == "tag":
        if tag is None or not _TAG_NAME.fullmatch(tag):
            raise ConfigurationError("tag must be a two-character BAM tag name such as 'RG' when by='tag'.")
        arguments += ["-tag", tag]
    else:
        if tag is not None:
            raise ConfigurationError("tag is only used with by='tag'.")
        arguments.append(f"-{mode}")
    if reference_prefix is not None:
        if mode != "reference":
            raise ConfigurationError("reference_prefix is only used with by='reference'.")
        arguments += ["-refPrefix", _gap_cli.text(reference_prefix, "reference_prefix")]
    if tag_prefix is not None:
        if mode != "tag":
            raise ConfigurationError("tag_prefix is only used with by='tag'.")
        arguments += ["-tagPrefix", _gap_cli.text(tag_prefix, "tag_prefix")]
    if tag_list_delimiter is not None:
        if mode != "tag":
            raise ConfigurationError("tag_list_delimiter is only used with by='tag'.")
        arguments += ["-tagListDelim", _gap_cli.text(tag_list_delimiter, "tag_list_delimiter")]
    result = _bamtools("split", arguments, cwd=stub.parent, expected=[], bamtools_path=bamtools_path)
    produced = sorted(str(path) for path in stub.parent.glob(f"{stub.name}.*.bam") if path.is_file())
    result["outputs"] = produced
    if result["status"] == "completed" and not produced:
        result["status"] = "failed"
    return result


__all__ = [
    "bamtools_convert_alignments",
    "bamtools_count_alignments",
    "bamtools_filter_alignments",
    "bamtools_index_bam_bti",
    "bamtools_resolve_pairs",
    "bamtools_revert_alignments",
    "bamtools_sample_alignments",
    "bamtools_split_alignments",
]
_register_public_verbose(globals())
