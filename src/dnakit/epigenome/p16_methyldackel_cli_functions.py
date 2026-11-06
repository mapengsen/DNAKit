"""MethylDackel facades (pack P16): cytosine report and context merging, methylation bias, per-read summary, BBM mappability.

Each function runs one MethylDackel sub-command with a fixed argument whitelist (no shell) through
``dnakit.bioinformatics._p16_cli``; the tested version is 0.6.1. The reference FASTA must be indexed (``samtools faidx``) and the
alignment sorted (and indexed for region queries). Results are JSON-serialisable dictionaries.
"""

from __future__ import annotations

import re
import shutil
import tempfile
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p16_cli as cli
from dnakit.exceptions import ConfigurationError

_TOOL = cli.Tool(
    "methyldackel",
    "MethylDackel",
    "Install MethylDackel (conda install -c bioconda methyldackel) and put it on PATH or pass methyldackel_path.",
    "0.6.1",
    "https://github.com/dpryan79/MethylDackel",
    ("--version",),
    None,
)
_BOUNDS = re.compile(r"^\d+,\d+,\d+,\d+$")
_SUGGESTION = re.compile(r"--(OT|OB|CTOT|CTOB)\s+(\d+,\d+,\d+,\d+)")


def _bounds(value: str | None, name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, str) or not _BOUNDS.match(value):
        raise ConfigurationError(
            f"{name} must be four comma-separated integers such as '1,100,1,100'."
        )
    return [f"--{name}", value]


def _prefix(value: str) -> Path:
    prefix = Path(value).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder of the output prefix does not exist: {prefix.parent}")
    if any(prefix.parent.glob(prefix.name + "[._]*")):
        raise ConfigurationError(
            f"Files with the prefix {prefix} already exist and will not be overwritten."
        )
    return prefix


def _filters(
    min_mapq: int, min_base_quality: int, region: str | None, regions_bed: str | None, threads: int
) -> list[str]:
    arguments = [
        "-q",
        str(cli.integer(min_mapq, "min_mapq", minimum=0, maximum=255)),
        "-p",
        str(cli.integer(min_base_quality, "min_base_quality", minimum=1, maximum=93)),
        "-@",
        str(cli.threads_value(threads)),
    ]
    if region is not None:
        arguments += ["-r", cli.text(region, "region")]
    if regions_bed is not None:
        arguments += ["-l", str(cli.existing_file(regions_bed, "regions_bed"))]
    return arguments


def extract_methylation_report_methyldackel(
    reference: str,
    alignment: str,
    output_prefix: str,
    *,
    report: str = "cytosine_report",
    include_chg: bool = False,
    include_chh: bool = False,
    include_cpg: bool = True,
    min_depth: int = 1,
    min_mapq: int = 10,
    min_base_quality: int = 5,
    region: str | None = None,
    regions_bed: str | None = None,
    inclusion_bounds_ot: str | None = None,
    inclusion_bounds_ob: str | None = None,
    threads: int = 1,
    methyldackel_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-cytosine methylation report, or per-CpG/CHG merged metrics, from a bisulfite BAM (``MethylDackel extract``).

    ``report`` selects the output: ``cytosine_report`` (``--cytosine_report``) writes the exhaustive per-cytosine table with
    chromosome, 1-based position, strand, methylated count, unmethylated count, CG/CHG/CHH and trinucleotide
    context for every cytosine including uncovered ones, to ``<prefix>.cytosine_report.txt``; ``merge_context``
    (``--mergeContext``) merges the two strands of each CpG/CHG into ``<prefix>_CpG.bedGraph`` (and ``_CHG``/``_CHH``);
    ``fraction`` (``--fraction``) writes ``<prefix>_CpG.meth.bedGraph``; ``per_cytosine`` is the plain extract
    (``<prefix>_CpG.bedGraph`` with one line per cytosine, the input of :func:`merge_cytosine_context_methyldackel`). ``include_chg``/``include_chh`` add contexts and
    ``include_cpg=False`` is ``--noCpG``. ``inclusion_bounds_ot``/``_ob`` are ``--OT``/``--OB`` bounds from
    :func:`methylation_bias_methyldackel`. The first 20 lines of the main result are returned in ``head``.
    """
    ref = cli.existing_file(reference, "reference")
    bam = cli.existing_file(alignment, "alignment")
    prefix = _prefix(output_prefix)
    mode = cli.choice(
        report, ("cytosine_report", "merge_context", "fraction", "per_cytosine"), "report"
    )
    arguments = [
        "extract",
        *_filters(min_mapq, min_base_quality, region, regions_bed, threads),
        "-d",
        str(cli.integer(min_depth, "min_depth", minimum=1)),
        "-o",
        str(prefix),
    ]
    if mode == "cytosine_report":
        arguments.append("--cytosine_report")
    elif mode == "merge_context":
        arguments.append("--mergeContext")
    elif mode == "fraction":
        arguments.append("--fraction")
    for enabled, name, switch in (
        (include_chg, "include_chg", "--CHG"),
        (include_chh, "include_chh", "--CHH"),
    ):
        if cli.flag(enabled, name):
            arguments.append(switch)
    if not cli.flag(include_cpg, "include_cpg"):
        arguments.append("--noCpG")
    arguments += _bounds(inclusion_bounds_ot, "OT") + _bounds(inclusion_bounds_ob, "OB")
    arguments += [str(ref), str(bam)]
    main = Path(
        str(prefix)
        + (
            ".cytosine_report.txt"
            if mode == "cytosine_report"
            else (
                "_CpG.bedGraph"
                if mode == "merge_context" and include_cpg
                else "_CpG.meth.bedGraph"
                if mode == "fraction" and include_cpg
                else "_CpG.bedGraph"
                if mode == "per_cytosine" and include_cpg
                else "_CHG.bedGraph"
                if include_chg
                else "_CHH.bedGraph"
            )
        )
    )
    record = cli.invoke(
        _TOOL, arguments, cwd=prefix.parent, expected=[main], explicit=methyldackel_path
    )
    outputs = sorted(str(p) for p in prefix.parent.glob(prefix.name + "[._]*"))
    head = [line.rstrip("\n") for _, line in zip(range(20), cli.open_text(main), strict=False)]
    return cli.finish(
        record, "extract_methylation_report_methyldackel", report=mode, files=outputs, head=head
    )


def merge_cytosine_context_methyldackel(
    reference: str,
    per_cytosine_metrics: str,
    output_path: str,
    *,
    methyldackel_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge per-cytosine metrics from ``MethylDackel extract`` into per-CpG/CHG metrics (``MethylDackel mergeContext``).

    The input must be coordinate sorted (for example ``<prefix>_CpG.bedGraph`` of a plain extract run).
    """
    ref = cli.existing_file(reference, "reference")
    source = cli.existing_file(per_cytosine_metrics, "per_cytosine_metrics")
    target = cli.new_file(output_path, "output_path")
    record = cli.invoke(
        _TOOL,
        ["mergeContext", "-o", str(target), str(ref), str(source)],
        cwd=target.parent,
        expected=[target],
        explicit=methyldackel_path,
    )
    return cli.finish(
        record, "merge_cytosine_context_methyldackel", sites=max(cli.count_lines(target) - 1, 0)
    )


def methylation_bias_methyldackel(
    reference: str,
    alignment: str,
    output_prefix: str,
    *,
    include_chg: bool = False,
    include_chh: bool = False,
    min_mapq: int = 10,
    min_base_quality: int = 5,
    region: str | None = None,
    regions_bed: str | None = None,
    max_depth: int = 2000,
    threads: int = 1,
    methyldackel_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Methylation bias by read position: SVG plots and suggested trim bounds (``MethylDackel mbias``).

    SVG files ``<prefix>_<strand>_<context>.svg`` are written; the suggested inclusion bounds the program prints
    (for use as ``--OT``, ``--OB``, ``--CTOT``, ``--CTOB`` in :func:`extract_methylation_report_methyldackel`) are parsed
    into ``suggested_bounds``.
    """
    ref = cli.existing_file(reference, "reference")
    bam = cli.existing_file(alignment, "alignment")
    prefix = _prefix(output_prefix)
    arguments = [
        "mbias",
        *_filters(min_mapq, min_base_quality, region, regions_bed, threads),
        "-D",
        str(cli.integer(max_depth, "max_depth", minimum=1)),
    ]
    for enabled, name, switch in (
        (include_chg, "include_chg", "--CHG"),
        (include_chh, "include_chh", "--CHH"),
    ):
        if cli.flag(enabled, name):
            arguments.append(switch)
    arguments += [str(ref), str(bam), str(prefix)]
    record = cli.invoke(_TOOL, arguments, cwd=prefix.parent, explicit=methyldackel_path)
    svgs = sorted(str(p) for p in prefix.parent.glob(prefix.name + "_*.svg"))
    if not svgs:
        raise cli.BackendExecutionError(
            "MethylDackel mbias wrote no SVG files.", context={"log_tail": record["log_tail"]}
        )
    record["outputs"] = svgs
    bounds = {name: value for name, value in _SUGGESTION.findall(record["log_tail"])}
    return cli.finish(record, "methylation_bias_methyldackel", suggested_bounds=bounds)


def summarize_per_read_methylation_methyldackel(
    reference: str,
    alignment: str,
    output_path: str,
    *,
    min_mapq: int = 10,
    min_base_quality: int = 5,
    region: str | None = None,
    regions_bed: str | None = None,
    threads: int = 1,
    methyldackel_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average CpG methylation per read (``MethylDackel perRead``): read, chromosome, position, methylation %, informative bases."""
    ref = cli.existing_file(reference, "reference")
    bam = cli.existing_file(alignment, "alignment")
    target = cli.new_file(output_path, "output_path")
    arguments = [
        "perRead",
        *_filters(min_mapq, min_base_quality, region, regions_bed, threads),
        "-o",
        str(target),
        str(ref),
        str(bam),
    ]
    record = cli.invoke(
        _TOOL, arguments, cwd=target.parent, expected=[target], explicit=methyldackel_path
    )
    rows = cli.read_table(
        target,
        header=False,
        names=["read", "chrom", "position", "methylation_percent", "informative_bases"],
        text_columns=frozenset({"read", "chrom"}),
        limit=1000,
    )
    return cli.finish(
        record,
        "summarize_per_read_methylation_methyldackel",
        reads=cli.count_lines(target),
        rows=rows,
    )


def convert_mappability_bigwig_to_bbm_methyldackel(
    reference: str,
    alignment: str,
    bigwig: str,
    bbm_path: str,
    output_prefix: str,
    *,
    mappability_threshold: float = 0.01,
    min_mappable_bases: int = 15,
    region: str | None = None,
    methyldackel_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a mappability bigWig to a Binary BisMap (``.bbm``) file (``MethylDackel extract -M ... -N``).

    MethylDackel writes the BBM as a by-product of an ``extract`` run that filters reads by the bigWig; this function runs
    that command with ``-O`` (the CpG bedGraph goes to ``output_prefix``) and returns the BBM path. ``mappability_threshold``
    is ``-t`` and ``min_mappable_bases`` is ``-b``. MethylDackel 0.6.1 does not accept ``-N <file>``, so the file is produced
    next to a symlink of the bigWig in a temporary folder and moved to ``bbm_path``.
    """
    ref = cli.existing_file(reference, "reference")
    bam = cli.existing_file(alignment, "alignment")
    wig = cli.existing_file(bigwig, "bigwig")
    bbm = cli.new_file(bbm_path, "bbm_path")
    prefix = _prefix(output_prefix)
    # MethylDackel 0.6.1 rejects -N <file>; -O writes <bigwig basename>.bbm next to the bigWig, so a symlink in a
    # private folder is used and the result is moved to bbm_path.
    with tempfile.TemporaryDirectory(dir=prefix.parent) as scratch:
        link = Path(scratch) / (wig.stem + ".bw")
        link.symlink_to(wig)
        arguments = [
            "extract",
            "-M",
            str(link),
            "-O",
            "-t",
            str(cli.number(mappability_threshold, "mappability_threshold", minimum=0.0)),
            "-b",
            str(cli.integer(min_mappable_bases, "min_mappable_bases", minimum=1)),
            "-o",
            str(prefix),
        ]
        if region is not None:
            arguments += ["-r", cli.text(region, "region")]
        arguments += [str(ref), str(bam)]
        produced = Path(scratch) / (wig.stem + ".bbm")
        record = cli.invoke(
            _TOOL, arguments, cwd=prefix.parent, expected=[produced], explicit=methyldackel_path
        )
        shutil.move(str(produced), str(bbm))
    record["outputs"] = [str(bbm)]
    return cli.finish(
        record,
        "convert_mappability_bigwig_to_bbm_methyldackel",
        bbm=str(bbm),
        bbm_bytes=bbm.stat().st_size,
    )


__all__ = [
    "convert_mappability_bigwig_to_bbm_methyldackel",
    "extract_methylation_report_methyldackel",
    "merge_cytosine_context_methyldackel",
    "methylation_bias_methyldackel",
    "summarize_per_read_methylation_methyldackel",
]

_register_public_verbose(globals(), tuple(__all__))
