"""More SeqKit command-line facades: FASTA index, FASTX sanitising, streaming, common records, header replace, BAM toolbox.

Each public function runs one SeqKit subcommand (2.14.0 is the tested version) with a fixed argument whitelist and no shell.
The user installs SeqKit (``conda install -c bioconda seqkit`` or a release binary) and puts it on PATH or passes
``seqkit_path``. Results are JSON-serialisable dictionaries (status, backend, version, command, outputs, parsed rows, log
tail). Output files must not exist yet; nothing is ever overwritten.

Complements ``dnakit.fasta_fastq.seqkit_cli_functions`` (which holds the first 21 SeqKit operations).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="seqkit",
    executable="seqkit",
    hint=(
        "Install SeqKit 2.14.0 or later (for example `conda install -c bioconda seqkit`, or a release binary from "
        "https://github.com/shenwei356/seqkit/releases), then put it on PATH or pass seqkit_path."
    ),
    tested_version="2.14.0",
    citation="https://doi.org/10.1002/imt2.191",
    version_arguments=("--help",),
    version_pattern=r"^Version:\s*(\S+)",
)
_FORMATS = ("fastq", "fasta")
_BAM_DUMP_FIELDS = (
    "Read",
    "Ref",
    "Pos",
    "EndPos",
    "MapQual",
    "Acc",
    "Match",
    "Mismatch",
    "Ins",
    "Del",
    "AlnLen",
    "ReadLen",
    "RefLen",
    "RefAln",
    "RefCov",
    "ReadAln",
    "ReadCov",
    "Strand",
    "MeanQual",
    "LeftClip",
    "RightClip",
    "Flags",
    "IsSec",
    "IsSup",
    "ReadSeq",
    "ReadAlnSeq",
    "LeftSoftClipSeq",
    "RightSoftClipSeq",
    "RightSoftClip",
    "LeftHardClip",
    "RightHardClip",
)
_PASS = re.compile(r"Pass records:\s*(\d+)\s+Discarded lines:\s*(\d+)")
_COMMON = re.compile(r"(\d+) common/shared sequences saved")
_DURATION = re.compile(r"^\d+(?:\.\d+)?(?:ns|us|ms|s|m|h)$")


def _inputs(values: Sequence[str | Path], name: str, minimum: int = 1) -> list[Path]:
    return _cli.existing_files(values, name, minimum=minimum)


def _common_flags(threads: int) -> list[str]:
    return ["-j", str(_cli.threads_value(threads))]


def _run(
    subcommand: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    seqkit_path: str | None,
    timeout_seconds: float,
    keep_info: bool = False,
) -> dict[str, Any]:
    arguments = [subcommand, *([] if keep_info else ["--quiet"]), *arguments]
    return _cli.invoke(
        _TOOL,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=seqkit_path,
        timeout_seconds=timeout_seconds,
    )


def _count_records(path: Path) -> int:
    return _cli.count_fastx_records(path)


def seqkit_faidx(
    fasta_path: str,
    *,
    regions: Sequence[str] | None = None,
    output_path: str | None = None,
    full_head: bool = False,
    use_regexp: bool = False,
    ignore_case: bool = False,
    update_index: bool = True,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the FASTA index (``.fai``) of a FASTA file and optionally extract regions (``seqkit faidx``).

    Without ``regions`` the function only builds the index and returns it as rows (``name``, ``length``, ``offset``,
    ``line_bases``, ``line_width``); this is the SeqKit counterpart of ``samtools faidx file.fa``. With ``regions`` (IDs or
    regions such as ``chr1:100-200``, 1-based, negative positions count from the end) the sequences are written to the new
    file ``output_path`` (required then). ``full_head`` prints the full header line and builds ``<fasta>.seqkit.fai`` instead
    of ``<fasta>.fai``; ``use_regexp`` reads the IDs as regular expressions; ``update_index`` (``-U``, the default) refreshes an
    index that no longer matches the file.
    """
    source = _cli.existing_file(fasta_path, "fasta_path")
    flags = ["-U"] if _cli.flag(update_index, "update_index") else []
    for enabled, option in ((full_head, "-f"), (use_regexp, "-r"), (ignore_case, "-i")):
        if _cli.flag(enabled, option):
            flags.append(option)
    expected: list[Path] = []
    positional: list[str] = [str(source)]
    target: Path | None = None
    if regions is not None:
        if isinstance(regions, str) or not regions:
            raise ConfigurationError("regions must be a non-empty list of IDs or regions.")
        positional += [_cli.free_text(item, "regions") for item in regions]
        if output_path is None:
            raise ConfigurationError("output_path is required when regions are given.")
    if output_path is not None:
        if regions is None:
            raise ConfigurationError(
                "output_path needs regions; the index itself is written next to the FASTA file."
            )
        target = _cli.new_file(output_path, "output_path")
        flags += ["-o", str(target)]
        expected.append(target)
    index = source.with_name(source.name + (".seqkit.fai" if full_head else ".fai"))
    record = _run(
        "faidx",
        [*flags, *positional],
        cwd=source.parent,
        expected=expected + [index],
        seqkit_path=seqkit_path,
        timeout_seconds=timeout_seconds,
    )
    rows = _cli.read_table(
        index,
        header=False,
        names=("name", "length", "offset", "line_bases", "line_width"),
        text_columns=frozenset({"name"}),
    )
    return _cli.finish(
        record,
        "build_fasta_index",
        index_file=str(index),
        sequence_count=len(rows),
        index=rows[:10000],
        output=str(target) if target else None,
        extracted_records=_count_records(target) if target else None,
    )


def seqkit_sanitize_fastx(
    input_path: str,
    output_path: str,
    *,
    format: str = "fastq",
    input_format: str | None = None,
    output_format: str | None = None,
    allow_gaps: bool = False,
    quality_ascii_base: int = 33,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Repair a broken single-line FASTQ or FASTA file by skipping malformed records (``seqkit sana``).

    SeqKit's resilient parser does not stop at the first error: records whose sequence and quality lengths differ, truncated
    final records and stray lines are dropped and the rest is written to the new file ``output_path``. ``format`` sets the
    input and output format (``fastq`` or ``fasta``); ``input_format`` and ``output_format`` override one side (for example
    ``input_format='fastq', output_format='fasta'`` converts while sanitising). ``allow_gaps`` keeps ``-`` in sequences;
    ``quality_ascii_base`` is 33 for Phred+33. Result keys: ``passed_records``, ``discarded_lines`` (parsed from SeqKit's
    summary) and ``output_record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = ["-i", _cli.choice(format, _FORMATS, "format")]
    if input_format is not None:
        arguments += ["-I", _cli.choice(input_format, _FORMATS, "input_format")]
    if output_format is not None:
        arguments += ["-O", _cli.choice(output_format, _FORMATS, "output_format")]
    if _cli.flag(allow_gaps, "allow_gaps"):
        arguments.append("-A")
    arguments += [
        "-b",
        str(_cli.integer(quality_ascii_base, "quality_ascii_base", minimum=0, maximum=126)),
    ]
    arguments += ["-o", str(target), str(source)]
    record = _run(
        "sana",
        arguments,
        cwd=target.parent,
        expected=[target],
        seqkit_path=seqkit_path,
        timeout_seconds=timeout_seconds,
        keep_info=True,
    )
    match = _PASS.search(record["log_tail"])
    return _cli.finish(
        record,
        "sanitize_broken_fastx",
        output=str(target),
        passed_records=int(match.group(1)) if match else None,
        discarded_lines=int(match.group(2)) if match else None,
        output_record_count=_count_records(target),
    )


def seqkit_stream_directory(
    directory: str,
    output_path: str,
    *,
    format: str = "fastq",
    find_only: bool = False,
    time_limit: str | None = None,
    gz_only: bool = False,
    file_pattern: str | None = None,
    minimum_size_increase_kb: int = 5,
    quality_ascii_base: int = 33,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate and stream the FASTX files of a directory as they appear (``seqkit scat``).

    SeqKit watches ``directory`` recursively (for example a sequencer output folder) and writes the records of every matching
    file, including data appended later, to the new file ``output_path``. ``find_only`` concatenates the files that exist now
    and quits. Otherwise ``time_limit`` (a duration such as ``30s`` or ``5m``) is required so the call ends: SeqKit stops when
    no file has grown for that long. ``format`` is ``fastq`` or ``fasta``; ``gz_only`` looks for ``.gz`` files only;
    ``file_pattern`` replaces the default file-name regular expression; ``minimum_size_increase_kb`` is the growth that triggers
    parsing. Result keys: ``passed_records``, ``discarded_lines`` and ``output_record_count``.
    """
    folder = _cli.existing_directory(directory, "directory")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        "-i",
        _cli.choice(format, _FORMATS, "format"),
        "-d",
        str(_cli.integer(minimum_size_increase_kb, "minimum_size_increase_kb", minimum=1)),
    ]
    arguments += [
        "-b",
        str(_cli.integer(quality_ascii_base, "quality_ascii_base", minimum=0, maximum=126)),
    ]
    if _cli.flag(find_only, "find_only"):
        arguments.append("-f")
    elif time_limit is None:
        raise ConfigurationError(
            "time_limit is required unless find_only is True, otherwise scat would never stop."
        )
    if time_limit is not None:
        if not isinstance(time_limit, str) or not _DURATION.match(time_limit):
            raise ConfigurationError("time_limit must be a duration such as '30s', '5m' or '1h'.")
        arguments += ["-T", time_limit]
    if _cli.flag(gz_only, "gz_only"):
        arguments.append("-g")
    if file_pattern is not None:
        arguments += ["-r", _cli.free_text(file_pattern, "file_pattern")]
    arguments += ["-o", str(target)]
    record = _cli.invoke(
        _TOOL,
        ["scat", *arguments, str(folder)],
        cwd=target.parent,
        expected=[target],
        explicit=seqkit_path,
        timeout_seconds=timeout_seconds,
    )
    match = _PASS.search(record["log_tail"])
    return _cli.finish(
        record,
        "realtime_fastx_directory_streaming",
        output=str(target),
        passed_records=int(match.group(1)) if match else None,
        discarded_lines=int(match.group(2)) if match else None,
        output_record_count=_count_records(target),
    )


def seqkit_common_records(
    input_paths: Sequence[str],
    output_path: str,
    *,
    by: str = "id",
    ignore_case: bool = False,
    only_positive_strand: bool = False,
    check_embedded_seqs: bool = False,
    threads: int = 1,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Records shared by two or more FASTA/FASTQ files, matched by ID, full name or sequence (``seqkit common``).

    ``by`` selects the key: ``'id'`` (first word of the header, the default), ``'name'`` (the full header, ``-n``) or
    ``'sequence'`` (``-s``). Sequence matching compares both strands unless ``only_positive_strand`` is True;
    ``check_embedded_seqs`` also reports sequences that are embedded in sequences of the other files (sequence mode only).
    ``ignore_case`` ignores case. The records of the first file that are shared by all files are written to the new file
    ``output_path`` in the format of the input; every record of the first file with a shared key is kept, so the count can be
    larger than the number of keys. Result keys: ``shared_record_count`` and ``input_file_count``.
    """
    files = _inputs(input_paths, "input_paths", minimum=2)
    target = _cli.new_file(output_path, "output_path")
    chosen = _cli.choice(by, ("id", "name", "sequence"), "by")
    arguments = _common_flags(threads)
    if chosen == "name":
        arguments.append("-n")
    elif chosen == "sequence":
        arguments.append("-s")
    if _cli.flag(ignore_case, "ignore_case"):
        arguments.append("-i")
    if _cli.flag(only_positive_strand, "only_positive_strand"):
        if chosen != "sequence":
            raise ConfigurationError("only_positive_strand applies to by='sequence'.")
        arguments.append("-P")
    if _cli.flag(check_embedded_seqs, "check_embedded_seqs"):
        if chosen != "sequence":
            raise ConfigurationError("check_embedded_seqs applies to by='sequence'.")
        arguments.append("-e")
    arguments += ["-o", str(target), *map(str, files)]
    record = _run(
        "common",
        arguments,
        cwd=target.parent,
        expected=[target],
        seqkit_path=seqkit_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "common_records_by_id_name_or_sequence",
        output=str(target),
        key=chosen,
        input_file_count=len(files),
        shared_record_count=_count_records(target),
    )


def seqkit_replace_header(
    input_path: str,
    output_path: str,
    pattern: str,
    replacement: str,
    *,
    ignore_case: bool = False,
    key_value_file: str | None = None,
    key_capture_index: int = 1,
    keep_key: bool = False,
    keep_untouched: bool = False,
    missing_key_replacement: str | None = None,
    number_width: int = 1,
    only_records_matching: str | None = None,
    only_records_matching_regexp: bool = False,
    threads: int = 1,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace sequence names (headers) with a regular expression (``seqkit replace`` without ``--by-seq``).

    ``pattern`` is a regular expression searched in the ID and ``replacement`` may use capture variables (``$1``) and the
    special symbols ``{nr}`` (record number), ``{fn}``, ``{fbn}``, ``{fbne}`` (file name parts), ``{uuid}`` and ``{kv}``
    (value looked up in ``key_value_file`` by the captured key number ``key_capture_index``). ``number_width`` zero-pads
    ``{nr}``; ``keep_key``, ``keep_untouched`` and ``missing_key_replacement`` choose what happens to keys without a value.
    ``only_records_matching`` restricts the edit to records whose ID matches that pattern (``--f-pattern``; a regular
    expression when ``only_records_matching_regexp`` is True). The sequence-replacement mode of the
    same subcommand is ``seqkit_replace_sequence_regex`` in ``dnakit.fasta_fastq.seqkit_cli_functions``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _cli.new_file(output_path, "output_path")
    arguments = _common_flags(threads)
    arguments += [
        "-p",
        _cli.free_text(pattern, "pattern"),
        "-r",
        _cli.free_text(replacement, "replacement"),
    ]
    if _cli.flag(ignore_case, "ignore_case"):
        arguments.append("-i")
    if key_value_file is not None:
        arguments += ["-k", str(_cli.existing_file(key_value_file, "key_value_file"))]
    arguments += ["-I", str(_cli.integer(key_capture_index, "key_capture_index", minimum=1))]
    if _cli.flag(keep_key, "keep_key"):
        arguments.append("-K")
    if _cli.flag(keep_untouched, "keep_untouched"):
        arguments.append("-U")
    if missing_key_replacement is not None:
        arguments += ["-m", _cli.free_text(missing_key_replacement, "missing_key_replacement")]
    arguments += ["--nr-width", str(_cli.integer(number_width, "number_width", minimum=1))]
    if only_records_matching is not None:
        arguments += ["--f-pattern", _cli.free_text(only_records_matching, "only_records_matching")]
        if _cli.flag(only_records_matching_regexp, "only_records_matching_regexp"):
            arguments.append("--f-use-regexp")
    arguments += ["-o", str(target), str(source)]
    record = _run(
        "replace",
        arguments,
        cwd=target.parent,
        expected=[target],
        seqkit_path=seqkit_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(
        record,
        "replace_header_regex",
        output=str(target),
        output_record_count=_count_records(target),
    )


def _toolbox(
    yaml: str, bam: Path, output: Path | None, seqkit_path: str | None, timeout: float
) -> dict[str, Any]:
    arguments = ["bam", "--quiet", "-Z", "-T", yaml]
    expected: list[Path] = []
    if output is not None:
        arguments += ["-o", str(output)]
        expected.append(output)
    arguments.append(str(bam))
    return _cli.invoke(
        _TOOL,
        arguments,
        cwd=(output or bam).parent,
        expected=expected,
        explicit=seqkit_path,
        timeout_seconds=timeout,
    )


def _yaml_text(value: str, name: str) -> str:
    """A YAML flow scalar for a path or a regular expression (JSON strings are valid YAML)."""
    import json

    return json.dumps(_cli.free_text(value, name))


def seqkit_bam_dump_records(
    bam_path: str,
    tsv_path: str,
    *,
    fields: Sequence[str] | None = None,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dump per-record alignment features of a BAM file as a TSV table (``seqkit bam -T Dump``).

    ``fields`` chooses the columns, in order, from ``Read, Ref, Pos, EndPos, MapQual, Acc, Match, Mismatch, Ins, Del, AlnLen,
    ReadLen, RefLen, RefAln, RefCov, ReadAln, ReadCov, Strand, MeanQual, LeftClip, RightClip, Flags, IsSec, IsSup, ReadSeq,
    ReadAlnSeq, LeftSoftClipSeq, RightSoftClipSeq, RightSoftClip, LeftHardClip, RightHardClip``; the default is all of them.
    Features that need the ``NM`` tag (``Acc``, ``Match``, ``Mismatch``) require a BAM written with it (``samtools calmd``).
    The new file ``tsv_path`` holds the table; the first ``rows`` (up to 10000) are returned as dictionaries.
    """
    bam = _cli.existing_file(bam_path, "bam_path")
    tsv = _cli.new_file(tsv_path, "tsv_path")
    chosen = list(_BAM_DUMP_FIELDS if fields is None else fields)
    if not chosen:
        raise ConfigurationError("fields must not be empty.")
    for name in chosen:
        _cli.choice(name, _BAM_DUMP_FIELDS, "fields")
    yaml = f"Dump: {{Fields: [{', '.join(chosen)}], Tsv: {_yaml_text(str(tsv), 'tsv_path')}}}\nSink: true"
    record = _toolbox(yaml, bam, None, seqkit_path, timeout_seconds)
    if not tsv.is_file():
        raise ConfigurationError("SeqKit did not write the dump table.")
    rows = _cli.read_table(
        tsv,
        text_columns=frozenset(
            name
            for name in chosen
            if name
            in {"Read", "Ref", "ReadSeq", "ReadAlnSeq", "LeftSoftClipSeq", "RightSoftClipSeq"}
        ),
    )
    record["outputs"] = [str(tsv)]
    return _cli.finish(
        record,
        "bam_record_feature_dump",
        tsv=str(tsv),
        fields=chosen,
        row_count=len(rows),
        rows=rows[:10000],
    )


def seqkit_bam_filter_by_alignment_context(
    bam_path: str,
    reference_fasta: str,
    output_bam: str,
    *,
    regex_start: str | None = None,
    regex_end: str | None = None,
    left_shift: int = -10,
    right_shift: int = 10,
    stranded: bool = False,
    invert: bool = False,
    tsv_path: str | None = None,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep BAM records by the reference sequence context at their start and end (``seqkit bam -T AlnContext``).

    For every mapped record SeqKit reads the reference window ``[position + left_shift, position + right_shift)`` around the
    alignment start and around the alignment end, and tests ``regex_start`` / ``regex_end`` (Go regular expressions) on
    them. A record is kept when either window matches (with ``invert=True`` when neither does). ``stranded`` reverse-
    complements the windows of reverse-strand records so the context is read in the direction of the read. The kept records
    go to the new BAM ``output_bam``; unmapped records pass through. ``reference_fasta`` must be indexed (SeqKit builds the
    ``.fai`` next to it). ``tsv_path`` (optional, new file) receives the table of examined windows with columns ``Read, Ref,
    Strand, StartSeq, StartMatch, EndSeq, EndMatch``.
    """
    bam = _cli.existing_file(bam_path, "bam_path")
    reference = _cli.existing_file(reference_fasta, "reference_fasta")
    target = _cli.new_file(output_bam, "output_bam")
    if regex_start is None and regex_end is None:
        raise ConfigurationError("give regex_start, regex_end or both.")
    parts = [
        f"Ref: {_yaml_text(str(reference), 'reference_fasta')}",
        f"LeftShift: {_cli.integer(left_shift, 'left_shift')}",
        f"RightShift: {_cli.integer(right_shift, 'right_shift')}",
        f"Stranded: {'true' if _cli.flag(stranded, 'stranded') else 'false'}",
        f"Invert: {'true' if _cli.flag(invert, 'invert') else 'false'}",
    ]
    if regex_start is not None:
        parts.append(f"RegexStart: {_yaml_text(regex_start, 'regex_start')}")
    if regex_end is not None:
        parts.append(f"RegexEnd: {_yaml_text(regex_end, 'regex_end')}")
    table: Path | None = None
    if tsv_path is not None:
        table = _cli.new_file(tsv_path, "tsv_path")
        parts.append(f"Tsv: {_yaml_text(str(table), 'tsv_path')}")
    record = _toolbox(
        "AlnContext: {" + ", ".join(parts) + "}", bam, target, seqkit_path, timeout_seconds
    )
    rows = (
        _cli.read_table(table, text_columns=frozenset({"Read", "Ref", "StartSeq", "EndSeq"}))
        if table and table.is_file()
        else []
    )
    if table is not None and table.is_file():
        record["outputs"].append(str(table))
    return _cli.finish(
        record,
        "bam_filter_by_alignment_context",
        output_bam=str(target),
        tsv=str(table) if table else None,
        context_rows=rows[:10000],
    )


def seqkit_bam_accuracy_stats(
    bam_path: str,
    *,
    tsv_path: str,
    seqkit_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean alignment accuracy of a BAM file, plain and weighted by alignment length (``seqkit bam -T AccStats``).

    Accuracy of a record is ``(1 - NM / alignment length) * 100`` computed from the ``NM`` tag, so the BAM must carry it
    (``samtools calmd -b``). The two values (``AccMean`` and ``WeightedAccMean``) are written to the new file ``tsv_path`` and
    returned as ``accuracy_mean`` and ``accuracy_weighted_mean``.
    """
    bam = _cli.existing_file(bam_path, "bam_path")
    tsv = _cli.new_file(tsv_path, "tsv_path")
    record = _toolbox(
        f"AccStats: {{Tsv: {_yaml_text(str(tsv), 'tsv_path')}}}\nSink: true",
        bam,
        None,
        seqkit_path,
        timeout_seconds,
    )
    if not tsv.is_file():
        raise ConfigurationError("SeqKit did not write the accuracy table.")
    rows = _cli.read_table(tsv)
    record["outputs"] = [str(tsv)]
    row = rows[0] if rows else {}
    return _cli.finish(
        record,
        "bam_weighted_mean_accuracy",
        tsv=str(tsv),
        accuracy_mean=row.get("AccMean"),
        accuracy_weighted_mean=row.get("WeightedAccMean"),
    )


__all__ = [
    "seqkit_bam_accuracy_stats",
    "seqkit_bam_dump_records",
    "seqkit_bam_filter_by_alignment_context",
    "seqkit_common_records",
    "seqkit_faidx",
    "seqkit_replace_header",
    "seqkit_sanitize_fastx",
    "seqkit_stream_directory",
]
_register_public_verbose(globals())
