"""seqtk command-line facades: fixed argument whitelists, no shell, bounded execution.

Each public function runs one seqtk sub-command (https://github.com/lh3/seqtk; 1.5-r133 is the tested version)
through ``dnakit.backends.external``. The user installs seqtk (``conda install -c bioconda seqtk``, a distribution
package or a build from source) and puts it on PATH or passes ``seqtk_path``. seqtk writes its result to standard
output; the facade streams it into a new explicit output file (an existing file is never overwritten) and returns a
JSON-serialisable dictionary with the status, backend, version, command, outputs, the parsed table where the result is
small, and a bounded log tail.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p12_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="seqtk",
    executable="seqtk",
    hint=(
        "Install seqtk 1.4 or later (for example `conda install -c bioconda seqtk`, or build it from "
        "https://github.com/lh3/seqtk), then put it on PATH or pass seqtk_path."
    ),
    tested_version="1.5-r133",
    citation="https://github.com/lh3/seqtk",
    version_arguments=(),
    version_pattern=r"Version:\s*(\S+)",
)
_MAX_RECORDS = 10_000
_BASES = re.compile(r"^[ACGTacgt]+$")
_MASK_CHARACTER = re.compile(r"^[!-~]$")


def _output(path: str, name: str = "output_path") -> Path:
    return _cli.new_file(path, name)


def _result(
    operation: str,
    arguments: list[str],
    output: Path,
    seqtk_path: str | None,
    **extra: Any,
) -> dict[str, Any]:
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=output.parent,
        stdout_to=output,
        explicit=seqtk_path,
    )
    return _cli.finish(record, operation, output_file=str(output), **extra)


def _fastx_records(path: Path) -> int:
    """Number of FASTA or FASTQ records in a plain-text seqtk output file."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if not first:
            return 0
        if first.startswith(">"):
            return 1 + sum(1 for line in handle if line.startswith(">"))
        if first.startswith("@"):
            return (1 + sum(1 for _ in handle) + 3) // 4
    return 0


def _rows(
    path: Path, names: list[str], text_columns: frozenset[str], limit: int
) -> tuple[list[dict[str, Any]], int]:
    """Parse a headerless tab-separated seqtk table: at most ``limit`` rows and the total row count."""
    rows: list[dict[str, Any]] = []
    total = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            total += 1
            if len(rows) < limit:
                rows.append(
                    {
                        name: (value if name in text_columns else _cli.convert_number(value))
                        for name, value in zip(names, fields, strict=False)
                    }
                )
    return rows, total


def _limit(value: int) -> int:
    return _cli.integer(value, "max_records", minimum=0, maximum=10_000_000)


def _flag(value: object, name: str) -> bool:
    return _cli.flag(value, name)


def seqtk_seq(
    input_path: str,
    output_path: str,
    *,
    drop_comments: bool = False,
    mask_quality_below: int | None = None,
    mask_character: str | None = None,
    quality_shift: int = 33,
    uppercase: bool = False,
    mate: int | None = None,
    line_length: int | None = None,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transform FASTA or FASTQ records with ``seqtk seq`` (comment removal, quality masking, upper-casing, mate selection).

    This single function covers four operations of the work list through its parameters: ``drop_comments=True`` removes
    the comment after the first white space of every header (``seqtk seq -C``); ``mask_quality_below=Q`` masks bases
    whose Phred quality is lower than Q (``-q Q``, the bases become lower case unless ``mask_character`` gives the
    replacement character, ``-n``); ``uppercase=True`` converts all bases to upper case (``-U``); ``mate=1`` keeps the
    odd records (reads 1, ``-1``) and ``mate=2`` the even records (reads 2, ``-2``) of an interleaved paired-end file.
    ``quality_shift`` is the ASCII offset of the quality characters (``-Q``, default 33). ``line_length`` re-wraps
    sequences (``-l``; None keeps every sequence on one line). With no option the records are copied unchanged.
    ``input_path`` may be gzip-compressed. Result keys: ``output_file`` and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    arguments = ["seq"]
    if _flag(drop_comments, "drop_comments"):
        arguments.append("-C")
    if mask_quality_below is not None:
        arguments += [
            "-q",
            str(_cli.integer(mask_quality_below, "mask_quality_below", minimum=1, maximum=93)),
        ]
    if mask_character is not None:
        if not isinstance(mask_character, str) or not _MASK_CHARACTER.match(mask_character):
            raise ConfigurationError("mask_character must be one visible ASCII character.")
        arguments += ["-n", mask_character]
    if quality_shift != 33:
        arguments += [
            "-Q",
            str(_cli.integer(quality_shift, "quality_shift", minimum=0, maximum=126)),
        ]
    if _flag(uppercase, "uppercase"):
        arguments.append("-U")
    if mate is not None:
        if isinstance(mate, bool) or mate not in (1, 2):
            raise ConfigurationError("mate must be 1, 2 or None.")
        arguments.append(f"-{mate}")
    if line_length is not None:
        arguments += ["-l", str(_cli.integer(line_length, "line_length", minimum=1))]
    arguments.append(str(source))
    result = _result("seqtk_seq", arguments, target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_size(
    input_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the sequences and the bases of a FASTA or FASTQ file with ``seqtk size``.

    ``input_path`` may be gzip-compressed. The tab-separated line written by seqtk is stored in ``output_path`` and
    parsed into ``sequence_count`` and ``base_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_size", ["size", str(source)], target, seqtk_path)
    rows, _ = _rows(target, ["sequence_count", "base_count"], frozenset(), 1)
    result["sequence_count"] = rows[0]["sequence_count"] if rows else 0
    result["base_count"] = rows[0]["base_count"] if rows else 0
    return result


def seqtk_fqchk(
    input_path: str,
    output_path: str,
    *,
    quality_threshold: int = 20,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise base composition and base quality per read position with ``seqtk fqchk``.

    ``quality_threshold`` is ``-q`` (default 20): the table then ends with the percentage of bases below and at or above
    that quality (``%low`` and ``%high``); 0 gives the percentage of every observed quality value instead. The table is
    saved in ``output_path``. Result keys: ``summary`` (``min_len``, ``max_len``, ``avg_len``, ``distinct_quality_values``),
    ``columns``, ``rows`` (the row ``ALL`` first, then one row per read position, at most ``max_records``) and
    ``row_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    threshold = _cli.integer(quality_threshold, "quality_threshold", minimum=0, maximum=93)
    limit = _limit(max_records)
    result = _result(
        "seqtk_fqchk", ["fqchk", "-q", str(threshold), str(source)], target, seqtk_path
    )
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    summary: dict[str, Any] = {}
    match = re.match(
        r"min_len: (\d+); max_len: (\d+); avg_len: ([\d.]+); (\d+) distinct quality values",
        lines[0] if lines else "",
    )
    if match:
        summary = {
            "min_len": int(match.group(1)),
            "max_len": int(match.group(2)),
            "avg_len": float(match.group(3)),
            "distinct_quality_values": int(match.group(4)),
        }
    columns = lines[1].split("\t") if len(lines) > 1 else []
    rows: list[dict[str, Any]] = []
    for line in lines[2 : 2 + limit]:
        fields = line.split("\t")
        rows.append(
            {
                name: (value if name == "POS" else _cli.convert_number(value))
                for name, value in zip(columns, fields, strict=False)
            }
        )
    result.update(summary=summary, columns=columns, rows=rows, row_count=max(len(lines) - 2, 0))
    return result


def seqtk_dropse(
    input_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop the unpaired reads from an interleaved paired-end FASTA or FASTQ file with ``seqtk dropse``.

    Two consecutive records are kept only when their names are the same except for a trailing ``/1`` ``/2``.
    ``output_path`` receives the paired reads. Result keys: ``output_file`` and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_dropse", ["dropse", str(source)], target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_hety(
    input_path: str,
    output_path: str,
    *,
    window_size: int = 50000,
    start_positions: int = 5,
    mask_lowercase: bool = False,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Measure the regional heterozygosity of a FASTA with IUPAC ambiguity codes with ``seqtk hety``.

    ``window_size`` is ``-w`` (default 50000), ``start_positions`` the number of window start positions per window
    (``-t``, default 5; the window slides by ``window_size / start_positions``) and ``mask_lowercase`` treats lower-case
    bases as masked (``-m``). Result keys: ``records`` (``seqid``, ``start``, ``end``, ``hets_per_window``,
    ``called_sites``, ``het_sites``; at most ``max_records``) and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    window = _cli.integer(window_size, "window_size", minimum=1)
    starts = _cli.integer(start_positions, "start_positions", minimum=1, maximum=window)
    arguments = ["hety", "-w", str(window), "-t", str(starts)]
    if _flag(mask_lowercase, "mask_lowercase"):
        arguments.append("-m")
    arguments.append(str(source))
    result = _result("seqtk_hety", arguments, target, seqtk_path)
    rows, total = _rows(
        target,
        ["seqid", "start", "end", "hets_per_window", "called_sites", "het_sites"],
        frozenset({"seqid"}),
        _limit(max_records),
    )
    result.update(records=rows, record_count=total)
    return result


def seqtk_gc(
    input_path: str,
    output_path: str,
    *,
    at_rich: bool = False,
    min_fraction: float = 0.6,
    min_length: int = 20,
    x_dropoff: float = 10.0,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Identify GC-rich (default) or AT-rich regions of a FASTA with the X-dropoff algorithm of ``seqtk gc``.

    ``at_rich=True`` looks for high-AT regions (``-w``) and then ``min_fraction`` is an AT fraction. ``min_fraction`` is
    ``-f`` (default 0.60), ``min_length`` the shortest region reported (``-l``, default 20) and ``x_dropoff`` the
    X-dropoff (``-x``, default 10.0). Result keys: ``records`` (``seqid``, ``start`` 0-based, ``end``, ``hit_bases``; at
    most ``max_records``) and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    fraction = _cli.number(min_fraction, "min_fraction", minimum=0.0, maximum=1.0)
    if not fraction > 0:
        raise ConfigurationError("min_fraction must be greater than 0.")
    arguments = ["gc"]
    if _flag(at_rich, "at_rich"):
        arguments.append("-w")
    arguments += [
        "-f",
        format(fraction, ".12g"),
        "-l",
        str(_cli.integer(min_length, "min_length", minimum=1)),
        "-x",
        format(_cli.number(x_dropoff, "x_dropoff", minimum=0.0), ".12g"),
        str(source),
    ]
    result = _result("seqtk_gc", arguments, target, seqtk_path)
    rows, total = _rows(
        target, ["seqid", "start", "end", "hit_bases"], frozenset({"seqid"}), _limit(max_records)
    )
    result.update(records=rows, record_count=total)
    return result


def seqtk_mutfa(
    input_path: str,
    snp_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Introduce point mutations into a FASTA at listed positions with ``seqtk mutfa``.

    ``snp_path`` is a text file whose lines have at least four columns: sequence name, 1-based position, any value
    (ignored) and the base to change to. The mutated FASTA is written to ``output_path``. Result keys: ``output_file``
    and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    snps = _cli.existing_file(snp_path, "snp_path")
    target = _output(output_path)
    result = _result("seqtk_mutfa", ["mutfa", str(source), str(snps)], target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_mergefa(
    first_path: str,
    second_path: str,
    output_path: str,
    *,
    quality_threshold: int = 0,
    intersection: bool = False,
    lowercase_on_n: bool = False,
    random_allele: bool = False,
    suppress_hets: bool = False,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge two FASTA/FASTQ allele sets (for example the two haplotypes of one individual) with ``seqtk mergefa``.

    Bases that differ become IUPAC heterozygous codes. ``quality_threshold`` is ``-q`` (default 0), ``intersection``
    takes the intersection of the inputs (``-i``), ``lowercase_on_n`` converts to lower case where one input base is N
    (``-m``), ``random_allele`` picks a random allele from heterozygous input bases (``-r``) and ``suppress_hets``
    removes heterozygous bases from the input (``-h``). Result keys: ``output_file`` and ``record_count``.
    """
    first = _cli.existing_file(first_path, "first_path")
    second = _cli.existing_file(second_path, "second_path")
    target = _output(output_path)
    arguments = [
        "mergefa",
        "-q",
        str(_cli.integer(quality_threshold, "quality_threshold", minimum=0, maximum=93)),
    ]
    for value, flag_name, option in (
        (intersection, "intersection", "-i"),
        (lowercase_on_n, "lowercase_on_n", "-m"),
        (random_allele, "random_allele", "-r"),
        (suppress_hets, "suppress_hets", "-h"),
    ):
        if _flag(value, flag_name):
            arguments.append(option)
    arguments += [str(first), str(second)]
    result = _result("seqtk_mergefa", arguments, target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_famask(
    source_path: str,
    mask_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply an X-coded mask FASTA to a source FASTA with ``seqtk famask``.

    Positions that are lower case or ``X`` in ``mask_path`` become lower case in the copy of ``source_path`` written to
    ``output_path``. Result keys: ``output_file`` and ``record_count``.
    """
    source = _cli.existing_file(source_path, "source_path")
    mask = _cli.existing_file(mask_path, "mask_path")
    target = _output(output_path)
    result = _result("seqtk_famask", ["famask", str(source), str(mask)], target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_randbase(
    input_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resolve every heterozygous IUPAC base of a FASTA to one of its two alleles at random with ``seqtk randbase``.

    seqtk draws with the C ``drand48`` generator without seeding, so a given file always gives the same output.
    Sequences are written 60 bases per line. Result keys: ``output_file`` and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_randbase", ["randbase", str(source)], target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_cutn(
    input_path: str,
    output_path: str,
    *,
    min_n_tract: int = 1000,
    non_n_penalty: int = 10,
    gaps_only: bool = False,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut sequences at long N tracts with ``seqtk cutN``, or list the tracts only.

    ``min_n_tract`` is ``-n`` (default 1000) and ``non_n_penalty`` the penalty for a non-N base inside a tract (``-p``,
    default 10). With ``gaps_only=True`` (``-g``) no sequence is written; ``output_path`` then holds the gaps as
    ``seqid``, ``start``, ``end`` rows, which are also returned in ``records``. Otherwise the pieces are written as FASTA
    records named ``seqid:start-end`` and ``record_count`` is their number.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    arguments = [
        "cutN",
        "-n",
        str(_cli.integer(min_n_tract, "min_n_tract", minimum=1)),
        "-p",
        str(_cli.integer(non_n_penalty, "non_n_penalty", minimum=0)),
    ]
    only_gaps = _flag(gaps_only, "gaps_only")
    if only_gaps:
        arguments.append("-g")
    arguments.append(str(source))
    result = _result("seqtk_cutn", arguments, target, seqtk_path)
    if only_gaps:
        rows, total = _rows(
            target, ["seqid", "start", "end"], frozenset({"seqid"}), _limit(max_records)
        )
        result.update(records=rows, record_count=total)
    else:
        result["record_count"] = _fastx_records(target)
    return result


def seqtk_gap(
    input_path: str,
    output_path: str,
    *,
    min_length: int = 50,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the N-gap locations (runs of non-ACGT bases) of at least ``min_length`` with ``seqtk gap``.

    ``min_length`` is ``-l`` (default 50). Result keys: ``records`` (``seqid``, 0-based ``start``, ``end``; at most
    ``max_records``) and ``record_count``; the same table is saved in ``output_path``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    arguments = ["gap", "-l", str(_cli.integer(min_length, "min_length", minimum=1)), str(source)]
    result = _result("seqtk_gap", arguments, target, seqtk_path)
    rows, total = _rows(
        target, ["seqid", "start", "end"], frozenset({"seqid"}), _limit(max_records)
    )
    result.update(records=rows, record_count=total)
    return result


def seqtk_listhet(
    input_path: str,
    output_path: str,
    *,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the position of every heterozygous (two-allele IUPAC) base of a FASTA with ``seqtk listhet``.

    Result keys: ``records`` (``seqid``, 1-based ``position``, ``base``; at most ``max_records``) and ``record_count``;
    the same table is saved in ``output_path``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_listhet", ["listhet", str(source)], target, seqtk_path)
    rows, total = _rows(
        target, ["seqid", "position", "base"], frozenset({"seqid", "base"}), _limit(max_records)
    )
    result.update(records=rows, record_count=total)
    return result


def seqtk_hpc(
    input_path: str,
    output_path: str,
    *,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Homopolymer-compress every sequence of a FASTA or FASTQ file with ``seqtk hpc``.

    Each run of identical bases is replaced by a single base. Result keys: ``output_file`` and ``record_count``.
    """
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_hpc", ["hpc", str(source)], target, seqtk_path)
    result["record_count"] = _fastx_records(target)
    return result


def seqtk_kfreq(
    kmer: str,
    input_path: str,
    output_path: str,
    *,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the occurrences of a k-mer and of its single-mismatch neighbours in every sequence with ``seqtk kfreq``.

    ``kmer`` is a DNA word of 1 to 12 bases (A, C, G, T). Both strands are scanned; the reported strand is the one with
    more neighbour hits. Result keys: ``records`` (``seqid``, ``length``, ``strand``, ``neighbour_count`` and
    ``exact_count``; at most ``max_records``) and ``record_count``; the same table is saved in ``output_path``.
    """
    if not isinstance(kmer, str) or not _BASES.match(kmer) or len(kmer) > 12:
        raise ConfigurationError("kmer must be 1 to 12 bases from A, C, G and T.")
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    result = _result("seqtk_kfreq", ["kfreq", kmer.upper(), str(source)], target, seqtk_path)
    rows, total = _rows(
        target,
        ["seqid", "length", "strand", "neighbour_count", "exact_count"],
        frozenset({"seqid", "strand"}),
        _limit(max_records),
    )
    result.update(records=rows, record_count=total, kmer=kmer.upper())
    return result


def seqtk_telo(
    input_path: str,
    output_path: str,
    *,
    motif: str = "CCCTAA",
    penalty: int = 1,
    max_drop: int = 2000,
    min_score: int = 300,
    max_records: int = _MAX_RECORDS,
    seqtk_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Identify telomere repeats at the ends of assembly sequences or long reads with ``seqtk telo``.

    ``motif`` is the repeat unit as written on the 5' end (``-m``, default ``CCCTAA``; the reverse complement is searched
    on the 3' end), ``penalty`` the score subtracted for a non-motif base (``-p``, default 1), ``max_drop`` the maximum
    score drop before extension stops (``-d``, default 2000) and ``min_score`` the minimum score to report (``-s``, default
    300). Result keys: ``records`` (``seqid``, ``start``, ``end``, ``length``; at most ``max_records``),
    ``record_count``, ``telomere_bases`` and ``total_bases`` (the totals seqtk prints to standard error).
    """
    if not isinstance(motif, str) or not _BASES.match(motif) or len(motif) > 30:
        raise ConfigurationError("motif must be 1 to 30 bases from A, C, G and T.")
    source = _cli.existing_file(input_path, "input_path")
    target = _output(output_path)
    arguments = [
        "telo",
        "-m",
        motif.upper(),
        "-p",
        str(_cli.integer(penalty, "penalty", minimum=0)),
        "-d",
        str(_cli.integer(max_drop, "max_drop", minimum=1)),
        "-s",
        str(_cli.integer(min_score, "min_score", minimum=1)),
        str(source),
    ]
    result = _result("seqtk_telo", arguments, target, seqtk_path)
    rows, total = _rows(
        target, ["seqid", "start", "end", "length"], frozenset({"seqid"}), _limit(max_records)
    )
    totals = re.findall(r"^(\d+)\t(\d+)$", result["log_tail"], re.MULTILINE)
    result.update(
        records=rows,
        record_count=total,
        telomere_bases=int(totals[-1][0]) if totals else None,
        total_bases=int(totals[-1][1]) if totals else None,
    )
    return result


__all__ = [
    "seqtk_cutn",
    "seqtk_dropse",
    "seqtk_famask",
    "seqtk_fqchk",
    "seqtk_gap",
    "seqtk_gc",
    "seqtk_hety",
    "seqtk_hpc",
    "seqtk_kfreq",
    "seqtk_listhet",
    "seqtk_mergefa",
    "seqtk_mutfa",
    "seqtk_randbase",
    "seqtk_seq",
    "seqtk_size",
    "seqtk_telo",
]
_register_public_verbose(globals())
