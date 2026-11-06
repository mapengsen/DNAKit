"""seqmagick command-line facades: residue slicing, alignment clean-up, record selection, quality filter and file summary.

Each public function runs one ``seqmagick`` subcommand (0.8.6 is the tested version) with a fixed argument whitelist and no
shell. The user installs seqmagick (``pip install seqmagick`` or ``conda install -c bioconda seqmagick``) and puts it on PATH or
passes ``seqmagick_path``. Output files must not exist yet. Results are JSON-serialisable dictionaries (status, backend,
version, command, outputs, parsed tables, log tail).

Complements ``dnakit.fasta_fastq.external.convert_seqmagick`` (format conversion, ungap, reverse, upper/lower, length and
head/tail filters, sequence de-duplication).
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
    name="seqmagick",
    executable="seqmagick",
    hint="Install seqmagick 0.8 or later (`pip install seqmagick` or `conda install -c bioconda seqmagick`) and put it on PATH or pass seqmagick_path.",
    tested_version="0.8.6",
    citation="https://github.com/fhcrc/seqmagick",
    version_arguments=("--version",),
    version_pattern=r"v?(\d+\.\d+\.\d+)",
)
_RANGES = re.compile(r"^(?:-?\d*:-?\d*)(?:,(?:-?\d*:-?\d*))*$")
_ALPHABETS = ("dna", "dna-ambiguous", "rna", "rna-ambiguous", "protein")
_FORMAT = re.compile(r"^[A-Za-z0-9-]{1,32}$")


def _ranges(value: str, name: str) -> str:
    if not isinstance(value, str) or not _RANGES.match(value) or value.strip(":,") == "":
        raise ConfigurationError(
            f"{name} must be 'start:end' ranges separated by commas, for example '3:10' or '1:5,20:25'."
        )
    return value


def _format(value: str | None, name: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, str) or not _FORMAT.match(value):
        raise ConfigurationError(
            f"{name} must be a seqmagick format name such as fasta, fastq, nexus or stockholm."
        )
    return [f"--{name.replace('_', '-')}", value]


def _convert(
    arguments: list[str],
    source: Path,
    destination: Path,
    *,
    input_format: str | None,
    output_format: str | None,
    alphabet: str | None,
    line_wrap: int | None,
    seqmagick_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    command = [
        "convert",
        *arguments,
        *_format(input_format, "input_format"),
        *_format(output_format, "output_format"),
    ]
    if alphabet is not None:
        command += ["--alphabet", _cli.choice(alphabet, _ALPHABETS, "alphabet")]
    if line_wrap is not None:
        command += ["--line-wrap", str(_cli.integer(line_wrap, "line_wrap", minimum=0))]
    command += [str(source), str(destination)]
    record = _cli.invoke(
        _TOOL,
        command,
        cwd=destination.parent,
        expected=[destination],
        explicit=seqmagick_path,
        timeout_seconds=timeout_seconds,
    )
    record["output"] = str(destination)
    record["output_record_count"] = (
        _cli.count_fastx_records(destination)
        if destination.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq"}
        else None
    )
    return record


def seqmagick_slice_residues(
    source: str,
    destination: str,
    *,
    cut: str | None = None,
    drop: str | None = None,
    mask: str | None = None,
    relative_to: str | None = None,
    input_format: str | None = None,
    output_format: str | None = None,
    alphabet: str | None = None,
    line_wrap: int | None = None,
    seqmagick_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut, delete or mask residues of every sequence by position range (``seqmagick convert --cut/--drop/--mask``).

    Ranges are 1-based, inclusive and written ``start:end`` (several separated by commas, either end may be empty; a negative
    start counts from the end of the sequence, for example ``-10:``). ``cut`` keeps only the residues in the ranges, ``drop``
    deletes them and ``mask`` replaces them with gap characters (``-``). ``relative_to`` is the ID of a sequence in an alignment
    whose non-gap residue numbers the ranges refer to. At least one of the three must be given; seqmagick applies them in the
    order cut, drop, mask. The result is written to the new file ``destination`` (format from its extension, or
    ``output_format``). Result keys: ``output`` and ``output_record_count``.
    """
    src = _cli.existing_file(source, "source")
    dst = _cli.new_file(destination, "destination")
    arguments: list[str] = []
    for flag, value in (("--cut", cut), ("--drop", drop), ("--mask", mask)):
        if value is not None:
            arguments.append(f"{flag}={_ranges(value, flag.lstrip('-'))}")
    if not arguments:
        raise ConfigurationError("give at least one of cut, drop and mask.")
    if relative_to is not None:
        arguments += ["--relative-to", _cli.text(relative_to, "relative_to")]
    record = _convert(
        arguments,
        src,
        dst,
        input_format=input_format,
        output_format=output_format,
        alphabet=alphabet,
        line_wrap=line_wrap,
        seqmagick_path=seqmagick_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(record, "slice_residues_by_position", cut=cut, drop=drop, mask=mask)


def seqmagick_clean_alignment(
    source: str,
    destination: str,
    *,
    dash_gap: bool = False,
    squeeze: bool = False,
    squeeze_threshold: float | None = None,
    prune_empty: bool = False,
    input_format: str | None = None,
    output_format: str | None = None,
    alphabet: str | None = None,
    line_wrap: int | None = None,
    seqmagick_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalise gap characters and remove gap-dominated columns or all-gap sequences of an alignment (``seqmagick convert``).

    ``dash_gap`` replaces the characters ``? . : ~`` with ``-`` in all sequences (``--dash-gap``). ``squeeze`` removes the
    columns that hold a gap in every sequence (``--squeeze``, the same as ``squeeze_threshold=1.0``); ``squeeze_threshold``
    (0 to 1) removes the columns whose gap proportion is at least that value (``--squeeze-threshold``). ``prune_empty`` removes
    sequences that consist only of gaps (``--prune-empty``). At least one option is required. The result goes to the new file
    ``destination``.
    """
    src = _cli.existing_file(source, "source")
    dst = _cli.new_file(destination, "destination")
    arguments: list[str] = []
    if _cli.flag(dash_gap, "dash_gap"):
        arguments.append("--dash-gap")
    if _cli.flag(squeeze, "squeeze"):
        arguments.append("--squeeze")
    if squeeze_threshold is not None:
        if squeeze:
            raise ConfigurationError("give squeeze or squeeze_threshold, not both.")
        arguments += [
            "--squeeze-threshold",
            format(
                _cli.number(squeeze_threshold, "squeeze_threshold", minimum=0.0, maximum=1.0),
                ".12g",
            ),
        ]
    if _cli.flag(prune_empty, "prune_empty"):
        arguments.append("--prune-empty")
    if not arguments:
        raise ConfigurationError(
            "give at least one of dash_gap, squeeze, squeeze_threshold and prune_empty."
        )
    record = _convert(
        arguments,
        src,
        dst,
        input_format=input_format,
        output_format=output_format,
        alphabet=alphabet,
        line_wrap=line_wrap,
        seqmagick_path=seqmagick_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(record, "clean_alignment_gaps")


def seqmagick_select_sequences(
    source: str,
    destination: str,
    *,
    include_ids_file: str | None = None,
    exclude_ids_file: str | None = None,
    seq_pattern_include: str | None = None,
    seq_pattern_exclude: str | None = None,
    input_format: str | None = None,
    output_format: str | None = None,
    alphabet: str | None = None,
    line_wrap: int | None = None,
    seqmagick_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep or drop sequences by an ID list or by a regular expression on the sequence (``seqmagick convert``).

    ``include_ids_file`` keeps only the sequences whose ID is listed (one per line, ``--include-from-file``) and
    ``exclude_ids_file`` removes the listed ones (``--exclude-from-file``). ``seq_pattern_include`` keeps the sequences whose
    residues match the regular expression (``--seq-pattern-include``) and ``seq_pattern_exclude`` removes them
    (``--seq-pattern-exclude``). At least one selector is required; several are applied together. The selected records go to the
    new file ``destination``.
    """
    src = _cli.existing_file(source, "source")
    dst = _cli.new_file(destination, "destination")
    arguments: list[str] = []
    if include_ids_file is not None:
        arguments += [
            "--include-from-file",
            str(_cli.existing_file(include_ids_file, "include_ids_file")),
        ]
    if exclude_ids_file is not None:
        arguments += [
            "--exclude-from-file",
            str(_cli.existing_file(exclude_ids_file, "exclude_ids_file")),
        ]
    if seq_pattern_include is not None:
        arguments += [
            "--seq-pattern-include",
            _cli.free_text(seq_pattern_include, "seq_pattern_include"),
        ]
    if seq_pattern_exclude is not None:
        arguments += [
            "--seq-pattern-exclude",
            _cli.free_text(seq_pattern_exclude, "seq_pattern_exclude"),
        ]
    if not arguments:
        raise ConfigurationError("give at least one selector.")
    record = _convert(
        arguments,
        src,
        dst,
        input_format=input_format,
        output_format=output_format,
        alphabet=alphabet,
        line_wrap=line_wrap,
        seqmagick_path=seqmagick_path,
        timeout_seconds=timeout_seconds,
    )
    return _cli.finish(record, "select_sequences")


def seqmagick_quality_filter(
    sequence_file: str,
    output_file: str,
    *,
    input_qual: str | None = None,
    min_mean_quality: float = 25.0,
    min_length: int = 200,
    max_length: int = 1000,
    ambiguous_action: str | None = None,
    max_ambiguous: int | None = None,
    pct_ambiguous: float | None = None,
    primer: str | None = None,
    barcode_file: str | None = None,
    barcode_header: bool = False,
    sample_map_out: str | None = None,
    quality_window: int = 0,
    quality_window_mean_qual: float | None = None,
    quality_window_prop: float = 1.0,
    details_out: str | None = None,
    seqmagick_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter reads by quality, length, ambiguous bases and a leading primer (``seqmagick quality-filter``).

    The defaults are seqmagick's: mean quality 25, length 200 to 1000 (longer reads are truncated). ``ambiguous_action``
    ``'truncate'`` cuts a read at the first ambiguous base (``N``) and ``'drop'`` removes the read; ``max_ambiguous`` and
    ``pct_ambiguous`` remove reads with too many ambiguous bases. ``primer`` is an IUPAC primer that a read must start with
    (the read is kept when its first bases match the primer, ambiguity codes included). ``quality_window`` (non-zero) truncates
    a read where the mean quality in that window falls below ``quality_window_mean_qual``. ``input_qual`` gives the quality file
    when ``sequence_file`` is FASTA. The filtered reads go to the new file ``output_file`` (format from its extension);
    ``details_out`` (optional, new file) receives the fate of every read. Result keys: ``report`` (rows of seqmagick's report
    table: filter name and the counts passed unchanged, passed changed, failed, total and proportion passed).
    """
    src = _cli.existing_file(sequence_file, "sequence_file")
    dst = _cli.new_file(output_file, "output_file")
    arguments = [
        "--min-mean-quality",
        format(_cli.number(min_mean_quality, "min_mean_quality", minimum=0.0), ".12g"),
        "--min-length",
        str(_cli.integer(min_length, "min_length", minimum=0)),
        "--max-length",
        str(_cli.integer(max_length, "max_length", minimum=1)),
    ]
    if input_qual is not None:
        arguments += ["--input-qual", str(_cli.existing_file(input_qual, "input_qual"))]
    if ambiguous_action is not None:
        arguments += [
            "--ambiguous-action",
            _cli.choice(ambiguous_action, ("truncate", "drop"), "ambiguous_action"),
        ]
    if max_ambiguous is not None:
        arguments += [
            "--max-ambiguous",
            str(_cli.integer(max_ambiguous, "max_ambiguous", minimum=0)),
        ]
    if pct_ambiguous is not None:
        arguments += [
            "--pct-ambiguous",
            format(_cli.number(pct_ambiguous, "pct_ambiguous", minimum=0.0, maximum=100.0), ".12g"),
        ]
    if primer is not None:
        if not isinstance(primer, str) or not re.fullmatch(
            r"[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+", primer
        ):
            raise ConfigurationError("primer must be an IUPAC nucleotide string.")
        arguments += ["--primer", primer]
    temporary: Path | None = None
    if barcode_file is not None:
        arguments += ["--barcode-file", str(_cli.existing_file(barcode_file, "barcode_file"))]
        if _cli.flag(barcode_header, "barcode_header"):
            arguments.append("--barcode-header")
    elif primer is not None:
        temporary = dst.with_name(dst.name + ".primer_only.csv")
        if temporary.exists():
            raise ConfigurationError(
                f"A helper file already exists next to output_file: {temporary}"
            )
        temporary.write_text("primer_only,\n", encoding="utf-8")
        arguments += ["--barcode-file", str(temporary)]
    if sample_map_out is not None:
        arguments += ["--map-out", str(_cli.new_file(sample_map_out, "sample_map_out"))]
    if quality_window:
        arguments += [
            "--quality-window",
            str(_cli.integer(quality_window, "quality_window", minimum=1)),
        ]
        if quality_window_mean_qual is not None:
            arguments += [
                "--quality-window-mean-qual",
                format(
                    _cli.number(quality_window_mean_qual, "quality_window_mean_qual", minimum=0.0),
                    ".12g",
                ),
            ]
        arguments += [
            "--quality-window-prop",
            format(
                _cli.number(quality_window_prop, "quality_window_prop", minimum=0.0, maximum=1.0),
                ".12g",
            ),
        ]
    expected: list[Path] = [dst]
    report = dst.with_name(dst.name + ".report.tsv")
    if report.exists():
        raise ConfigurationError(f"A report file already exists next to output_file: {report}")
    arguments += ["--report-out", str(report)]
    expected.append(report)
    if details_out is not None:
        details = _cli.new_file(details_out, "details_out")
        arguments += ["--details-out", str(details)]
        expected.append(details)
    try:
        record = _cli.invoke(
            _TOOL,
            ["quality-filter", *arguments, str(src), str(dst)],
            cwd=dst.parent,
            expected=expected,
            explicit=seqmagick_path,
            timeout_seconds=timeout_seconds,
        )
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    rows = _cli.read_table(report, text_columns=frozenset({"name"}))
    return _cli.finish(
        record,
        "filter_reads_by_quality_ambiguity_primer",
        output=str(dst),
        report_file=str(report),
        report=rows,
        output_record_count=_cli.count_fastx_records(dst),
    )


def seqmagick_info(
    sequence_files: Sequence[str],
    output_file: str,
    *,
    input_format: str | None = None,
    table_format: str = "tab",
    threads: int = 1,
    seqmagick_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise sequence files: record count, length range and mean length, and whether they are alignments (``seqmagick info``).

    One row per file with the columns ``name``, ``alignment`` (``TRUE`` when all sequences have the same length), ``min_len``,
    ``max_len``, ``avg_len`` and ``num_seqs``. ``input_format`` overrides the format guessed from the extensions;
    ``table_format`` is ``'tab'``, ``'csv'`` or ``'align'``. The table is written to the new file ``output_file`` and the rows
    are returned as ``files`` (only for ``tab`` and ``csv``).
    """
    sources = _cli.existing_files(sequence_files, "sequence_files")
    target = _cli.new_file(output_file, "output_file")
    chosen = _cli.choice(table_format, ("tab", "csv", "align"), "table_format")
    arguments = [
        "info",
        "--format",
        chosen,
        "--threads",
        str(_cli.threads_value(threads)),
        "--out-file",
        str(target),
    ]
    arguments += _format(input_format, "input_format")
    arguments += [str(path) for path in sources]
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=target.parent,
        expected=[target],
        explicit=seqmagick_path,
        timeout_seconds=timeout_seconds,
    )
    rows = (
        _cli.read_table(
            target,
            delimiter="," if chosen == "csv" else "\t",
            text_columns=frozenset({"name", "alignment"}),
        )
        if chosen != "align"
        else []
    )
    return _cli.finish(
        record, "summarize_sequence_files", output=str(target), file_count=len(sources), files=rows
    )


__all__ = [
    "seqmagick_clean_alignment",
    "seqmagick_info",
    "seqmagick_quality_filter",
    "seqmagick_select_sequences",
    "seqmagick_slice_residues",
]
_register_public_verbose(globals())
