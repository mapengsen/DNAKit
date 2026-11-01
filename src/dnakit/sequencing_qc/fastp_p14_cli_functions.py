"""fastp command-line facade with the options that ``run_fastp`` does not expose (tested with fastp 1.4.0):
adapter FASTA files, low-complexity filtering, index-barcode filtering, overlap-based base correction, UMI
preprocessing and over-represented sequence analysis, together with the common trimming and filtering options.

The user installs fastp (``conda install -c bioconda fastp``, or a release binary from
https://github.com/OpenGene/fastp) and puts it on PATH, or passes ``fastp_path``. The function runs ``fastp`` with an
argument list built from typed, validated parameters: no shell, outputs are new files. Only the options that the caller
sets are passed, so every other option keeps fastp's own default. Results are JSON-serialisable dictionaries (status,
backend, version, command, return code, outputs, log tail, parsed JSON report summary).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.exceptions import ConfigurationError

_HINT = "Install fastp 0.23 or later (for example `conda install -c bioconda fastp`) and put it on PATH or pass fastp_path."
_UMI_LOCATIONS = ("index1", "index2", "read1", "read2", "per_index", "per_read")


def _optional_file(value: str | None, flag_name: str, name: str) -> list[str]:
    return [] if value is None else [flag_name, str(_gap_cli.existing_file(value, name))]


def fastp_process_reads(
    read1: str,
    output_read1: str,
    json_report: str,
    *,
    read2: str | None = None,
    output_read2: str | None = None,
    html_report: str | None = None,
    adapter_sequence: str | None = None,
    adapter_sequence_read2: str | None = None,
    adapter_fasta: str | None = None,
    detect_adapter_for_pe: bool = False,
    disable_adapter_trimming: bool = False,
    low_complexity_filter: bool = False,
    complexity_threshold: int | None = None,
    filter_by_index1: str | None = None,
    filter_by_index2: str | None = None,
    index_filter_threshold: int | None = None,
    correction: bool = False,
    overlap_length_required: int | None = None,
    overlap_diff_limit: int | None = None,
    overlap_diff_percent_limit: int | None = None,
    umi: bool = False,
    umi_location: str | None = None,
    umi_length: int | None = None,
    umi_prefix: str | None = None,
    umi_skip: int | None = None,
    umi_delimiter: str | None = None,
    overrepresentation_analysis: bool = False,
    overrepresentation_sampling: int | None = None,
    trim_poly_x: bool = False,
    cut_front: bool = False,
    cut_tail: bool = False,
    cut_right: bool = False,
    cut_window_size: int | None = None,
    cut_mean_quality: int | None = None,
    qualified_quality_phred: int | None = None,
    unqualified_percent_limit: int | None = None,
    n_base_limit: int | None = None,
    average_quality: int | None = None,
    length_required: int | None = None,
    threads: int | None = None,
    fastp_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim and filter FASTQ reads with fastp, including the options that ``run_fastp`` does not offer.

    ``read1`` (and ``read2`` for pairs; then ``output_read2`` is required) are FASTQ files (gzip is read; an output
    name ending in ``.gz`` is written compressed). Adapters: ``adapter_sequence`` / ``adapter_sequence_read2`` give the
    adapters, ``adapter_fasta`` a FASTA of adapters trimmed from both reads, ``detect_adapter_for_pe`` detects
    paired-end adapters, ``disable_adapter_trimming`` turns trimming off. Filters: ``low_complexity_filter`` keeps
    reads whose complexity (percent of bases that differ from the next base) is at least ``complexity_threshold``
    (default 30); ``filter_by_index1`` / ``filter_by_index2`` are files of index barcodes (one per line) whose reads
    are removed, allowing ``index_filter_threshold`` differences; ``qualified_quality_phred``,
    ``unqualified_percent_limit``, ``n_base_limit``, ``average_quality`` and ``length_required`` are the standard
    quality filters. ``correction`` corrects mismatched bases inside the overlap of read pairs
    (``overlap_length_required``, ``overlap_diff_limit``, ``overlap_diff_percent_limit`` define the overlap).
    UMI: ``umi=True`` moves the UMI into the read name; ``umi_location`` is 'index1', 'index2', 'read1', 'read2',
    'per_index' or 'per_read', ``umi_length`` is needed for read1/read2, ``umi_prefix``, ``umi_skip`` and
    ``umi_delimiter`` shape the name. ``overrepresentation_analysis`` reports over-represented sequences (every
    ``overrepresentation_sampling``-th read is examined). Trimming: ``trim_poly_x``, ``cut_front`` / ``cut_tail`` /
    ``cut_right`` (sliding windows of ``cut_window_size`` with mean quality ``cut_mean_quality``). ``json_report`` (and
    optionally ``html_report``) is the quality profile report; a summary of the JSON is returned in ``summary``.
    Writes new files.
    """
    first = _gap_cli.existing_file(read1, "read1")
    out1 = _gap_cli.new_file(output_read1, "output_read1")
    json_path = _gap_cli.new_file(json_report, "json_report")
    arguments = ["-i", str(first), "-o", str(out1), "--json", str(json_path)]
    expected = [out1, json_path]
    if (read2 is None) != (output_read2 is None):
        raise ConfigurationError("read2 and output_read2 must be given together.")
    if read2 is not None and output_read2 is not None:
        out2 = _gap_cli.new_file(output_read2, "output_read2")
        arguments += ["-I", str(_gap_cli.existing_file(read2, "read2")), "-O", str(out2)]
        expected.append(out2)
    if html_report is not None:
        html = _gap_cli.new_file(html_report, "html_report")
        arguments += ["--html", str(html)]
        expected.append(html)
    folders = {path.parent for path in expected}
    if len(folders) != 1:
        raise ConfigurationError("All output files must be in the same folder.")
    for value, flag_name, name in (
        (adapter_sequence, "--adapter_sequence", "adapter_sequence"),
        (adapter_sequence_read2, "--adapter_sequence_r2", "adapter_sequence_read2"),
    ):
        if value is not None:
            arguments += [flag_name, _p14_cli.nucleotides(value, name, alphabet="ACGTNacgtn")]
    arguments += _optional_file(adapter_fasta, "--adapter_fasta", "adapter_fasta")
    arguments += _optional_file(filter_by_index1, "--filter_by_index1", "filter_by_index1")
    arguments += _optional_file(filter_by_index2, "--filter_by_index2", "filter_by_index2")
    if umi_location is not None and not umi:
        raise ConfigurationError("umi_location needs umi=True.")
    if umi:
        arguments.append("--umi")
        if umi_location is None:
            raise ConfigurationError("umi=True needs umi_location.")
        arguments += ["--umi_loc", _gap_cli.choice(umi_location, _UMI_LOCATIONS, "umi_location")]
        if umi_location in ("read1", "read2"):
            if umi_length is None:
                raise ConfigurationError("umi_length is required when umi_location is read1 or read2.")
        if umi_length is not None:
            arguments += ["--umi_len", str(_gap_cli.integer(umi_length, "umi_length", minimum=1))]
        if umi_prefix is not None:
            arguments += ["--umi_prefix", _gap_cli.text(umi_prefix, "umi_prefix")]
        if umi_skip is not None:
            arguments += ["--umi_skip", str(_gap_cli.integer(umi_skip, "umi_skip", minimum=0))]
        if umi_delimiter is not None:
            arguments += ["--umi_delim", _gap_cli.text(umi_delimiter, "umi_delimiter")]
    elif any(item is not None for item in (umi_length, umi_prefix, umi_skip, umi_delimiter)):
        raise ConfigurationError("umi_length, umi_prefix, umi_skip and umi_delimiter need umi=True.")
    numbers = (
        ("complexity_threshold", complexity_threshold, "--complexity_threshold", 0, 100),
        ("index_filter_threshold", index_filter_threshold, "--filter_by_index_threshold", 0, None),
        ("overlap_length_required", overlap_length_required, "--overlap_len_require", 1, None),
        ("overlap_diff_limit", overlap_diff_limit, "--overlap_diff_limit", 0, None),
        ("overlap_diff_percent_limit", overlap_diff_percent_limit, "--overlap_diff_percent_limit", 0, 100),
        ("overrepresentation_sampling", overrepresentation_sampling, "--overrepresentation_sampling", 1, 10000),
        ("cut_window_size", cut_window_size, "--cut_window_size", 1, 1000),
        ("cut_mean_quality", cut_mean_quality, "--cut_mean_quality", 1, 36),
        ("qualified_quality_phred", qualified_quality_phred, "--qualified_quality_phred", 0, 93),
        ("unqualified_percent_limit", unqualified_percent_limit, "--unqualified_percent_limit", 0, 100),
        ("n_base_limit", n_base_limit, "--n_base_limit", 0, None),
        ("average_quality", average_quality, "--average_qual", 0, 93),
        ("length_required", length_required, "--length_required", 0, None),
        ("threads", threads, "--thread", 1, 256),
    )
    for name, value, flag_name, low, high in numbers:
        if value is not None:
            arguments += [flag_name, str(_gap_cli.integer(value, name, minimum=low, maximum=high))]
    for value, flag_name, name in (
        (detect_adapter_for_pe, "--detect_adapter_for_pe", "detect_adapter_for_pe"),
        (disable_adapter_trimming, "--disable_adapter_trimming", "disable_adapter_trimming"),
        (low_complexity_filter, "--low_complexity_filter", "low_complexity_filter"),
        (correction, "--correction", "correction"),
        (overrepresentation_analysis, "--overrepresentation_analysis", "overrepresentation_analysis"),
        (trim_poly_x, "--trim_poly_x", "trim_poly_x"),
        (cut_front, "--cut_front", "cut_front"),
        (cut_tail, "--cut_tail", "cut_tail"),
        (cut_right, "--cut_right", "cut_right"),
    ):
        if _p14_cli.flag(value, name):
            arguments.append(flag_name)
    if (correction or detect_adapter_for_pe) and read2 is None:
        raise ConfigurationError("correction and detect_adapter_for_pe need paired reads (read2).")
    result = _p14_cli.run(
        "fastp", arguments, cwd=next(iter(folders)), expected=expected, explicit=fastp_path, hint=_HINT,
        timeout_seconds=_gap_cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0),
    )
    if result["status"] == "completed":
        report = json.loads(Path(json_path).read_text(encoding="utf-8"))
        summary = report.get("summary", {})
        result["summary"] = {
            "before_filtering": summary.get("before_filtering"),
            "after_filtering": summary.get("after_filtering"),
            "filtering_result": report.get("filtering_result"),
            "duplication_rate": (report.get("duplication") or {}).get("rate"),
            "adapter_cutting": report.get("adapter_cutting"),
            "overrepresented": {key: value for key, value in report.items() if "overrepresented" in key} or None,
        }
    return result


__all__ = ["fastp_process_reads"]
_register_public_verbose(globals())
