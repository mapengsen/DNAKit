"""minimap2 and paftools.js command-line facades (pack P09): fixed argument lists, no shell, bounded execution.

``minimap2_map_reads`` runs the ``minimap2`` aligner (2.31 is the tested version); every other function runs one
sub-command of ``paftools.js`` (shipped with minimap2, executed by the ``k8`` JavaScript shell) or the companion script
``pafcluster.js``. The user installs the tools (``conda install -c bioconda minimap2 k8``, or the release tarball from
https://github.com/lh3/minimap2/releases, which contains ``misc/paftools.js`` and ``misc/pafcluster.js``) and either puts
them on PATH or passes ``minimap2_path``, ``paftools_path``, ``k8_path`` and ``pafcluster_path``. Programs that print to
standard output are streamed into a file inside a fresh ``output_directory`` (never overwritten); standard error goes to
a log file next to it. Results are JSON-serialisable dictionaries with the status, version, command, output files, small
parsed results and a bounded log tail. Optional parameters default to ``None`` and then are not passed, so the tool's own
defaults apply.
"""

from __future__ import annotations

import contextlib
import os
import re
import signal
import subprocess
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    ConfigurationError,
)

_TESTED_MINIMAP2 = "2.31"
_CITATION = "https://doi.org/10.1093/bioinformatics/bty191"
_MINIMAP2_HINT = (
    "Install minimap2 2.31 or later (for example `conda install -c bioconda minimap2`, or the release tarball from "
    "https://github.com/lh3/minimap2/releases), then put it on PATH or pass minimap2_path."
)
_PAFTOOLS_HINT = (
    "paftools.js is part of minimap2 (`conda install -c bioconda minimap2 k8` installs both; the release tarball has "
    "misc/paftools.js). Put paftools.js and the k8 JavaScript shell on PATH or pass paftools_path and k8_path."
)
_PAFCLUSTER_HINT = (
    "pafcluster.js is in the misc/ folder of the minimap2 source tree (https://github.com/lh3/minimap2/tree/master/misc). "
    "Download it, make it available and pass pafcluster_path (k8 must be installed)."
)
_PRESETS = (
    "map-ont",
    "map-hifi",
    "map-pb",
    "lr:hq",
    "splice",
    "splice:hq",
    "splice:sr",
    "sr",
    "asm5",
    "asm10",
    "asm20",
    "ava-pb",
    "ava-ont",
)
_PAF_COLUMNS = (
    "query",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target",
    "target_length",
    "target_start",
    "target_end",
    "matches",
    "alignment_length",
    "mapping_quality",
)
_SIZE = re.compile(r"^\d+(\.\d+)?[kKmMgG]?$")
_MAX_ARTIFACT_BYTES = 2_000_000_000
_MAX_LOG_BYTES = 5_000_000
_ROW_LIMIT = 1000
_REPORT_LIMIT = 200_000


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _switch(arguments: list[str], flag: str, value: object, name: str) -> None:
    if _flag(value, name):
        arguments.append(flag)


def _opt_int(
    arguments: list[str],
    flag: str,
    value: int | None,
    name: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> None:
    if value is not None:
        arguments.extend((flag, str(_cli.integer(value, name, minimum=minimum, maximum=maximum))))


def _opt_num(
    arguments: list[str],
    flag: str,
    value: float | None,
    name: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> None:
    if value is not None:
        arguments.extend(
            (flag, format(_cli.number(value, name, minimum=minimum, maximum=maximum), ".12g"))
        )


def _opt_text(arguments: list[str], flag: str, value: str | None, name: str) -> None:
    if value is not None:
        arguments.extend((flag, _cli.text(value, name)))


def _opt_file(arguments: list[str], flag: str, value: str | Path | None, name: str) -> None:
    if value is not None:
        arguments.extend((flag, str(_cli.existing_file(value, name))))


def _opt_size(arguments: list[str], flag: str, value: int | str | None, name: str) -> None:
    """A length such as ``1000000`` or ``1m`` (paftools.js accepts the k, m and g suffixes)."""
    if value is None:
        return
    text = str(value)
    if isinstance(value, bool) or _SIZE.fullmatch(text) is None:
        raise ConfigurationError(f"{name} must be a number optionally followed by k, m or g.")
    arguments.extend((flag, text))


def _file(value: str | Path, name: str) -> str:
    return str(_cli.existing_file(value, name))


def _files(values: Sequence[str | Path], name: str) -> list[str]:
    if isinstance(values, (str, Path)) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of files.")
    return [_file(item, name) for item in values]


# ---------------------------------------------------------------------------
# Running programs
# ---------------------------------------------------------------------------


def _stream(
    executable: str,
    arguments: Sequence[str],
    directory: Path,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: float,
) -> int:
    """Run one program without a shell, streaming standard output and standard error into files."""
    timeout = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        try:
            process = subprocess.Popen(
                [executable, *arguments],
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                cwd=directory,
                shell=False,
                env={**os.environ, "LANG": "C", "LC_ALL": "C"},
                start_new_session=os.name == "posix",
            )
        except OSError as exc:
            raise BackendExecutionError(
                "External command could not be started.",
                code="MINIMAP2_START_FAILED",
                context={"error": str(exc)},
            ) from exc
        deadline = time.monotonic() + timeout
        try:
            while process.poll() is None:
                if time.monotonic() > deadline:
                    raise BackendTimeoutError(f"{Path(executable).name} exceeded its time limit.")
                if (
                    stdout_path.stat().st_size > _MAX_ARTIFACT_BYTES
                    or stderr_path.stat().st_size > _MAX_LOG_BYTES
                ):
                    raise BackendExecutionError(
                        f"{Path(executable).name} exceeded its output size limit."
                    )
                time.sleep(0.02)
        finally:
            if process.poll() is None:
                if os.name == "posix":
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(process.pid, signal.SIGKILL)
                else:
                    process.kill()
            process.wait()
    return int(process.returncode)


def _tail(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[-_cli.LOG_TAIL_CHARS :]


def _version(name: str, explicit: str | None, hint: str, scratch: Path) -> str | None:
    probe = _cli.run(
        name,
        ["--version"],
        cwd=scratch,
        expected=[],
        explicit=explicit,
        hint=hint,
        timeout_seconds=60.0,
    )
    return probe["version"] if isinstance(probe["version"], str) else None


def _resolve_paftools(paftools_path: str | None, k8_path: str | None) -> tuple[str, str]:
    """The k8 shell and the paftools.js script (explicit path, then PATH)."""
    k8 = _cli.program("k8", explicit=k8_path, hint=_PAFTOOLS_HINT)
    script = _cli.program("paftools.js", explicit=paftools_path, hint=_PAFTOOLS_HINT)
    return k8, script


def _paftools_version(k8: str, script: str, scratch: Path) -> str | None:
    probe = _cli.run(
        "k8",
        [script, "version"],
        cwd=scratch,
        expected=[],
        explicit=k8,
        hint=_PAFTOOLS_HINT,
        timeout_seconds=60.0,
    )
    text = str(probe["log_tail"]).strip()
    return text.splitlines()[0].strip() if text else None


def _result(
    operation: str,
    backend: str,
    version: str | None,
    command: Sequence[str],
    directory: Path,
    outputs: dict[str, Path],
    log: Path,
    **extra: object,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": backend,
        "operation": operation,
        "version": version,
        "command": list(command),
        "output_directory": str(directory),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "log_tail": _tail(log),
        "provenance": {
            "tool": "minimap2",
            "tested_version": _TESTED_MINIMAP2,
            "citation": _CITATION,
        },
    }
    result.update(extra)
    return result


def _paftools(
    subcommand: str,
    arguments: Sequence[str],
    output_directory: str | Path,
    output_name: str,
    *,
    operation: str,
    paftools_path: str | None,
    k8_path: str | None,
    timeout_seconds: float,
    script: str | None = None,
    **extra: object,
) -> dict[str, Any]:
    """Run ``k8 paftools.js <subcommand>`` (or another k8 script) with standard output streamed to a file."""
    k8, paftools_script = (
        _resolve_paftools(paftools_path, k8_path)
        if script is None
        else (
            _cli.program("k8", explicit=k8_path, hint=_PAFTOOLS_HINT),
            script,
        )
    )
    directory = _cli.new_directory(output_directory, "output_directory")
    output = directory / output_name
    log = directory / "tool.log"
    command = [paftools_script, *([subcommand] if subcommand else []), *arguments]
    return_code = _stream(k8, command, directory, output, log, timeout_seconds)
    if return_code != 0:
        raise BackendExecutionError(
            f"{operation} failed (exit code {return_code}).",
            code="PAFTOOLS_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": _tail(log)},
        )
    version = _paftools_version(k8, paftools_script, directory) if script is None else None
    return _result(
        operation,
        "paftools.js" if script is None else Path(script).name,
        version,
        ["k8", Path(paftools_script).name, *([subcommand] if subcommand else []), *arguments],
        directory,
        {"output": output, "log": log},
        log,
        **extra,
    )


def _table(path: Path, limit: int = _ROW_LIMIT) -> tuple[list[list[str]], int]:
    """Tab-separated lines of a result file (header and comment lines excluded) and their total count."""
    rows: list[list[str]] = []
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            total += 1
            if len(rows) < limit:
                rows.append(line.rstrip("\n").split("\t"))
    return rows, total


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")[:_REPORT_LIMIT]


def _number(value: str) -> Any:
    for kind in (int, float):
        try:
            return kind(value)
        except ValueError:
            continue
    return value


def _slug(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def _colon_metrics(text: str) -> dict[str, Any]:
    """``label: value`` lines (optionally starting with ``#``) as a dictionary with snake_case keys.

    ``6 (100.00%)`` becomes ``label: 6`` and ``label_percent: 100.0``; comma-separated values become lists.
    """
    metrics: dict[str, Any] = {}
    for line in text.splitlines():
        match = re.match(r"^#?\s*(.+?):\s*(\S.*)$", line)
        if not match:
            continue
        key, raw = _slug(match.group(1)), match.group(2).strip()
        counted = re.fullmatch(r"(\d+) \(([\d.]+)%\)", raw)
        if counted:
            metrics[key], metrics[f"{key}_percent"] = int(counted.group(1)), float(counted.group(2))
        elif "," in raw:
            metrics[key] = [_number(item) for item in raw.split(",")]
        else:
            metrics[key] = _number(raw)
    return metrics


# ---------------------------------------------------------------------------
# minimap2 aligner
# ---------------------------------------------------------------------------


def minimap2_map_reads(
    reference: str | Path,
    query: str | Path,
    output_directory: str | Path,
    *,
    preset: Literal[
        "map-ont",
        "map-hifi",
        "map-pb",
        "lr:hq",
        "splice",
        "splice:hq",
        "splice:sr",
        "sr",
        "asm5",
        "asm10",
        "asm20",
        "ava-pb",
        "ava-ont",
    ] = "map-ont",
    mate: str | Path | None = None,
    output_format: Literal["paf", "sam"] = "paf",
    base_level_alignment: bool = False,
    cs_tag: Literal["short", "long"] | None = None,
    md_tag: bool = False,
    eqx: bool = False,
    kmer_size: int | None = None,
    window_size: int | None = None,
    max_intron_length: int | str | None = None,
    splice_strand: Literal["f", "r", "b", "n"] | None = None,
    junction_bed: str | Path | None = None,
    secondary: bool | None = None,
    max_secondary: int | None = None,
    soft_clip_supplementary: bool = False,
    threads: int | None = None,
    minimap2_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align reads or sequences to a reference with minimap2 (``minimap2 -x <preset>``).

    ``preset`` is the ``-x`` preset: ``map-ont``, ``map-hifi``, ``map-pb``, ``lr:hq``, ``sr``, ``asm5``/``asm10``/``asm20``,
    the spliced RNA or cDNA presets ``splice``, ``splice:hq`` and ``splice:sr`` (``splice_strand`` is ``-u``,
    ``max_intron_length`` is ``-G`` and ``junction_bed`` is ``--junc-bed``), and the all-vs-all read overlap presets
    ``ava-pb`` and ``ava-ont`` (pass the same file as ``reference`` and ``query``). ``mate`` is the second file of a
    paired short-read library (``sr``). ``output_format='sam'`` adds ``-a``; ``base_level_alignment`` adds ``-c`` (CIGAR in
    PAF). ``cs_tag``/``md_tag``/``eqx`` add ``--cs``, ``--MD`` and ``--eqx``. Writes ``alignments.paf`` or
    ``alignments.sam`` plus ``minimap2.log``; ``alignment_count`` counts the alignment lines and, for PAF, ``records`` holds
    the first rows (12 standard columns plus a ``tags`` dictionary).
    """
    _cli.choice(preset, _PRESETS, "preset")
    _cli.choice(output_format, ("paf", "sam"), "output_format")
    target = _file(reference, "reference")
    reads = _file(query, "query")
    second = None if mate is None else _file(mate, "mate")
    if splice_strand is not None:
        _cli.choice(splice_strand, ("f", "r", "b", "n"), "splice_strand")
    executable = _cli.program("minimap2", explicit=minimap2_path, hint=_MINIMAP2_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    output = directory / f"alignments.{output_format}"
    log = directory / "minimap2.log"
    arguments = ["-x", preset]
    _switch(arguments, "-a", output_format == "sam", "output_format")
    _switch(arguments, "-c", base_level_alignment, "base_level_alignment")
    if cs_tag is not None:
        arguments.append(
            "--cs=long" if _cli.choice(cs_tag, ("short", "long"), "cs_tag") == "long" else "--cs"
        )
    _switch(arguments, "--MD", md_tag, "md_tag")
    _switch(arguments, "--eqx", eqx, "eqx")
    _opt_int(arguments, "-k", kmer_size, "kmer_size", minimum=1, maximum=28)
    _opt_int(arguments, "-w", window_size, "window_size", minimum=1, maximum=255)
    _opt_size(arguments, "-G", max_intron_length, "max_intron_length")
    _opt_text(arguments, "-u", splice_strand, "splice_strand")
    _opt_file(arguments, "--junc-bed", junction_bed, "junction_bed")
    if secondary is not None:
        arguments.append("--secondary=yes" if _flag(secondary, "secondary") else "--secondary=no")
    _opt_int(arguments, "-N", max_secondary, "max_secondary", minimum=0)
    _switch(arguments, "-Y", soft_clip_supplementary, "soft_clip_supplementary")
    _opt_int(arguments, "-t", threads, "threads", minimum=1, maximum=1024)
    arguments.extend((target, reads))
    if second is not None:
        arguments.append(second)
    return_code = _stream(executable, arguments, directory, output, log, timeout_seconds)
    if return_code != 0:
        raise BackendExecutionError(
            f"minimap2 failed (exit code {return_code}).",
            code="MINIMAP2_EXECUTION_FAILED",
            context={"log_tail": _tail(log)},
        )
    version = _version("minimap2", minimap2_path, _MINIMAP2_HINT, directory)
    count = 0
    records: list[dict[str, Any]] = []
    with output.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            count += 1
            if output_format == "paf" and len(records) < _ROW_LIMIT:
                fields = line.rstrip("\n").split("\t")
                row: dict[str, Any] = {
                    name: _number(value) for name, value in zip(_PAF_COLUMNS, fields, strict=False)
                }
                row["tags"] = {
                    tag[:2]: tag[5:] for tag in fields[12:] if len(tag) > 5 and tag[2] == ":"
                }
                records.append(row)
    return _result(
        "map_reads",
        "minimap2",
        version,
        ["minimap2", *arguments],
        directory,
        {"alignments": output, "log": log},
        log,
        preset=preset,
        alignment_count=count,
        **({"records": records} if output_format == "paf" else {}),
    )


# ---------------------------------------------------------------------------
# paftools.js: format conversion
# ---------------------------------------------------------------------------


def minimap2_view_paf(
    paf: str | Path,
    output_directory: str | Path,
    *,
    output_format: Literal["aln", "maf", "lastz-cigar"] = "aln",
    line_length: int | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Render PAF alignments as BLAST-like text, MAF or lastz-cigar (``paftools.js view``).

    The PAF must carry the ``cs`` tag (``minimap2 -c --cs``). ``output_format`` is ``aln`` (BLAST-like, ``-f aln``), ``maf``
    or ``lastz-cigar``; ``line_length`` is ``-l`` for the BLAST-like layout (default 80). Writes ``view.txt``.
    """
    _cli.choice(output_format, ("aln", "maf", "lastz-cigar"), "output_format")
    arguments = ["-f", output_format]
    _opt_int(arguments, "-l", line_length, "line_length", minimum=1)
    arguments.append(_file(paf, "paf"))
    result = _paftools(
        "view",
        arguments,
        output_directory,
        "view.txt",
        operation="view_paf",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["output_format"] = output_format
    return result


def minimap2_sam_to_paf(
    sam: str | Path,
    output_directory: str | Path,
    *,
    alignments: Literal["all", "primary_and_supplementary", "primary"] = "all",
    include_unmapped: bool = False,
    long_cs: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert SAM alignments to PAF (``paftools.js sam2paf``).

    ``alignments`` selects ``all`` records, ``primary_and_supplementary`` (``-p``) or ``primary`` only (``-P``);
    ``include_unmapped`` is ``-U`` and ``long_cs`` writes the long form of the cs tag (``-L``). Writes ``converted.paf``.
    """
    _cli.choice(alignments, ("all", "primary_and_supplementary", "primary"), "alignments")
    arguments: list[str] = []
    if alignments == "primary_and_supplementary":
        arguments.append("-p")
    elif alignments == "primary":
        arguments.append("-P")
    _switch(arguments, "-U", include_unmapped, "include_unmapped")
    _switch(arguments, "-L", long_cs, "long_cs")
    arguments.append(_file(sam, "sam"))
    result = _paftools(
        "sam2paf",
        arguments,
        output_directory,
        "converted.paf",
        operation="sam_to_paf",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["row_count"] = _table(Path(result["outputs"]["output"]), 5)
    return result


def minimap2_delta_to_paf(
    delta: str | Path,
    output_directory: str | Path,
    *,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert MUMmer ``.delta`` alignments to PAF (``paftools.js delta2paf``). Writes ``converted.paf``."""
    result = _paftools(
        "delta2paf",
        [_file(delta, "delta")],
        output_directory,
        "converted.paf",
        operation="delta_to_paf",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["row_count"] = _table(Path(result["outputs"]["output"]), 5)
    return result


def minimap2_spliced_alignments_to_bed12(
    alignments: str | Path,
    output_directory: str | Path,
    *,
    keep_multiple: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert spliced alignments in PAF (with the ``cg`` tag) or SAM to BED12 gene structures (``paftools.js splice2bed``).

    ``keep_multiple`` keeps secondary mappings (``-m``). Writes ``structures.bed`` and parses it into ``records`` (BED12
    columns).
    """
    arguments: list[str] = []
    _switch(arguments, "-m", keep_multiple, "keep_multiple")
    arguments.append(_file(alignments, "alignments"))
    result = _paftools(
        "splice2bed",
        arguments,
        output_directory,
        "structures.bed",
        operation="spliced_alignments_to_bed12",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _bed_rows(
        Path(result["outputs"]["output"]), BED12_COLUMNS
    )
    return result


BED12_COLUMNS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "block_count",
    "block_sizes",
    "block_starts",
)


def _bed_rows(path: Path, names: Sequence[str]) -> tuple[list[dict[str, Any]], int]:
    rows, total = _table(path)
    return [
        {name: _number(value) for name, value in zip(names, row, strict=False)} for row in rows
    ], total


def minimap2_annotation_to_bed12(
    annotation: str | Path,
    output_directory: str | Path,
    *,
    junctions_only: bool = False,
    short_names: bool = False,
    fai: str | Path | None = None,
    ensembl_canonical_only: bool = False,
    output_gff: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GTF or GENCODE-style GFF3 gene annotation to BED12 (``paftools.js gff2bed``).

    GFF3 records need ``transcript_id=`` attributes as written by GENCODE/Ensembl; plain ``Parent=`` GFF3 is rejected by
    the tool. ``junctions_only`` writes a junction BED (``-j``), ``short_names`` shortens the names (``-s``), ``fai`` is a
    FASTA index for chromosome name conversion (``-u``), ``ensembl_canonical_only`` keeps transcripts tagged
    ``Ensembl_canonical`` (``-e``) and ``output_gff`` writes GFF instead of BED (``-g``, with ``fai``). Writes
    ``annotation.bed`` (``annotation.gff`` for ``output_gff``).
    """
    arguments: list[str] = []
    _switch(arguments, "-j", junctions_only, "junctions_only")
    _switch(arguments, "-s", short_names, "short_names")
    _opt_file(arguments, "-u", fai, "fai")
    _switch(arguments, "-e", ensembl_canonical_only, "ensembl_canonical_only")
    _switch(arguments, "-g", output_gff, "output_gff")
    if output_gff and fai is None:
        raise ConfigurationError("output_gff needs fai (paftools.js gff2bed -g is used with -u).")
    arguments.append(_file(annotation, "annotation"))
    name = "annotation.gff" if output_gff else "annotation.bed"
    result = _paftools(
        "gff2bed",
        arguments,
        output_directory,
        name,
        operation="annotation_to_bed12",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["record_count"] = _table(Path(result["outputs"]["output"]), 50)
    return result


def minimap2_gff3_to_junctions(
    gff3: str | Path,
    output_directory: str | Path,
    *,
    feature: str | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract intron junctions from a GFF3 file as BED (``paftools.js gff2junc``).

    Introns are the gaps between consecutive ``feature`` records (tool default ``CDS``; use ``exon`` for exon gaps) that
    share a ``Parent=``. Writes ``junctions.bed`` (chrom, start, end, transcript, score, strand; parsed into ``records``).
    """
    arguments: list[str] = []
    _opt_text(arguments, "-f", feature, "feature")
    arguments.append(_file(gff3, "gff3"))
    result = _paftools(
        "gff2junc",
        arguments,
        output_directory,
        "junctions.bed",
        operation="gff3_to_junctions",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _bed_rows(
        Path(result["outputs"]["output"]), ("chrom", "start", "end", "name", "score", "strand")
    )
    return result


def minimap2_longcs_to_sequences(
    paf: str | Path,
    output_directory: str | Path,
    *,
    query: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reconstruct the aligned target (default) or query (``query=True``, ``-q``) sequences from a long-cs PAF.

    Runs ``paftools.js longcs2seq``; the PAF needs ``--cs=long``. Writes ``sequences.fasta``.
    """
    arguments: list[str] = []
    _switch(arguments, "-q", query, "query")
    arguments.append(_file(paf, "paf"))
    result = _paftools(
        "longcs2seq",
        arguments,
        output_directory,
        "sequences.fasta",
        operation="longcs_to_sequences",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    with Path(result["outputs"]["output"]).open(encoding="utf-8") as handle:
        result["record_count"] = sum(1 for line in handle if line.startswith(">"))
    return result


def minimap2_paf_to_gff3(
    paf: str | Path,
    output_directory: str | Path,
    *,
    amino_acid: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert miniprot protein-to-genome PAF to GFF3 gene models (``paftools.js paf2gff``).

    ``amino_acid`` is ``-a`` (use the amino-acid cs). Writes ``models.gff3``.
    """
    arguments: list[str] = []
    _switch(arguments, "-a", amino_acid, "amino_acid")
    arguments.append(_file(paf, "paf"))
    result = _paftools(
        "paf2gff",
        arguments,
        output_directory,
        "models.gff3",
        operation="paf_to_gff3",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["record_count"] = _table(Path(result["outputs"]["output"]), 50)
    return result


def minimap2_convert_simulated_reads(
    simulator: Literal["mason2", "pbsim", "badread"],
    input_paths: Sequence[str | Path],
    output_directory: str | Path,
    *,
    reference_fai: str | Path | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert simulated reads to FASTQ/FASTA whose names carry the true coordinates.

    ``simulator='mason2'`` converts a mason2 SAM to FASTQ (``paftools.js mason2fq``); ``'pbsim'`` converts one or more PBSIM
    MAF files to FASTA (``pbsim2fq``, needs ``reference_fai``); ``'badread'`` converts a Badread FASTQ to FASTA
    (``badread2fa``, needs ``reference_fai``). Writes ``reads.fastq`` (mason2) or ``reads.fasta``; the names have the form
    ``name!chrom!start!end!strand`` that ``minimap2_evaluate_mapping`` reads.
    """
    _cli.choice(simulator, ("mason2", "pbsim", "badread"), "simulator")
    paths = _files(input_paths, "input_paths")
    if simulator == "mason2":
        if len(paths) != 1 or reference_fai is not None:
            raise ConfigurationError(
                "simulator='mason2' takes exactly one SAM file and no reference_fai."
            )
        subcommand, arguments, name = "mason2fq", paths, "reads.fastq"
    else:
        if reference_fai is None:
            raise ConfigurationError(f"simulator={simulator!r} needs reference_fai.")
        if simulator == "badread" and len(paths) != 1:
            raise ConfigurationError("simulator='badread' takes exactly one FASTQ file.")
        subcommand = "pbsim2fq" if simulator == "pbsim" else "badread2fa"
        arguments, name = [_file(reference_fai, "reference_fai"), *paths], "reads.fasta"
    result = _paftools(
        subcommand,
        arguments,
        output_directory,
        name,
        operation="convert_simulated_reads",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        simulator=simulator,
    )
    with Path(result["outputs"]["output"]).open(encoding="utf-8") as handle:
        lines = sum(1 for _ in handle)
    result["record_count"] = lines // 4 if simulator == "mason2" else lines // 2
    return result


def minimap2_simulated_names_to_bed(
    names_file: str | Path,
    output_directory: str | Path,
    *,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert simulated read names (``name!chrom!start!end!strand``, one per line) to BED truth (``paftools.js sim2bed``).

    Writes ``truth.bed`` (chrom, start, end, name, score, strand; parsed into ``records``).
    """
    result = _paftools(
        "sim2bed",
        [_file(names_file, "names_file")],
        output_directory,
        "truth.bed",
        operation="simulated_names_to_bed",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"], result["record_count"] = _bed_rows(
        Path(result["outputs"]["output"]), ("chrom", "start", "end", "name", "score", "strand")
    )
    return result


# ---------------------------------------------------------------------------
# paftools.js: alignment and assembly statistics
# ---------------------------------------------------------------------------


def minimap2_alignment_stats(
    alignments: str | Path,
    output_directory: str | Path,
    *,
    per_alignment_cs: bool = False,
    gap_output_length: int | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary statistics of PAF or SAM alignments (``paftools.js stat``).

    Reports mapped sequences, primary and secondary alignments, bases, substitutions and indel counts by size.
    ``per_alignment_cs`` is ``-c`` (use the cs tag of PAF) and ``gap_output_length`` is ``-l`` (print gaps of at least
    that length). Writes ``stats.txt``; ``metrics`` holds the numbers.
    """
    arguments: list[str] = []
    _switch(arguments, "-c", per_alignment_cs, "per_alignment_cs")
    _opt_int(arguments, "-l", gap_output_length, "gap_output_length", minimum=0)
    arguments.append(_file(alignments, "alignments"))
    result = _paftools(
        "stat",
        arguments,
        output_directory,
        "stats.txt",
        operation="alignment_stats",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    text = _text(Path(result["outputs"]["output"]))
    result["metrics"] = {
        _slug(label): int(value)
        for label, value in re.findall(r"^Number of (.+?): (\d+)$", text, re.MULTILINE)
    }
    return result


def minimap2_assembly_stats(
    fai: str | Path,
    alignments: Sequence[str | Path],
    output_directory: str | Path,
    *,
    min_query_length: int | None = None,
    min_block_length: int | None = None,
    max_divergence: float | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assembly contiguity statistics from assembly-to-reference alignments (``paftools.js asmstat``).

    ``fai`` is the reference FASTA index and ``alignments`` one PAF per assembly. ``min_query_length`` is ``-q`` (default
    0), ``min_block_length`` ``-l`` (10000) and ``max_divergence`` ``-d`` (0.01). Writes ``assembly_stats.tsv``;
    ``metrics`` maps each metric (Length, NG50, NGA50, ...) to one value per assembly.
    """
    arguments: list[str] = []
    _opt_int(arguments, "-q", min_query_length, "min_query_length", minimum=0)
    _opt_int(arguments, "-l", min_block_length, "min_block_length", minimum=0)
    _opt_num(arguments, "-d", max_divergence, "max_divergence", minimum=0.0)
    arguments.extend([_file(fai, "fai"), *_files(alignments, "alignments")])
    result = _paftools(
        "asmstat",
        arguments,
        output_directory,
        "assembly_stats.tsv",
        operation="assembly_stats",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    rows, _ = _table(Path(result["outputs"]["output"]), 10_000)
    result["metrics"] = {
        row[0]: [_number(value) for value in row[1:]]
        for row in rows
        if row[0] != "Metric" and len(row) > 1
    }
    return result


def minimap2_assembly_gene_completeness(
    reference_paf: str | Path,
    assembly_pafs: Sequence[str | Path],
    output_directory: str | Path,
    *,
    min_identity: float | None = None,
    min_coverage: float | None = None,
    autosomes_only: bool = False,
    show_incomplete: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gene completeness of assemblies from spliced alignments of the same transcripts (``paftools.js asmgene``).

    ``reference_paf`` aligns the transcripts (cDNA) to the reference, each of ``assembly_pafs`` to an assembly
    (``minimap2 -cx splice``). ``min_identity`` is ``-i`` (0.99), ``min_coverage`` ``-c`` (0.99), ``autosomes_only`` ``-a``
    and ``show_incomplete`` prints fragmented and missing genes (``-e``). Writes ``completeness.tsv``; ``metrics`` maps
    full_sgl, full_dup, frag, part50+, ... to one count per input.
    """
    arguments: list[str] = []
    _opt_num(arguments, "-i", min_identity, "min_identity", minimum=0.0, maximum=1.0)
    _opt_num(arguments, "-c", min_coverage, "min_coverage", minimum=0.0, maximum=1.0)
    _switch(arguments, "-a", autosomes_only, "autosomes_only")
    _switch(arguments, "-e", show_incomplete, "show_incomplete")
    arguments.extend(
        [_file(reference_paf, "reference_paf"), *_files(assembly_pafs, "assembly_pafs")]
    )
    result = _paftools(
        "asmgene",
        arguments,
        output_directory,
        "completeness.tsv",
        operation="assembly_gene_completeness",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    metrics: dict[str, list[Any]] = {}
    for line in _text(Path(result["outputs"]["output"])).splitlines():
        fields = line.split("\t")
        if fields[0] == "X" and len(fields) > 2:
            metrics[fields[1]] = [_number(value) for value in fields[2:]]
    result["metrics"] = metrics
    return result


def minimap2_detect_misjoins(
    paf: str | Path,
    output_directory: str | Path,
    *,
    centromere_bed: str | Path | None = None,
    centromere_overlap_ratio: float | None = None,
    min_block_length: int | str | None = None,
    max_gap: int | str | None = None,
    show_non_centromeric: bool = False,
    show_long_blocks: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Detect large-scale misjoins in assembly-to-reference alignments (``paftools.js misjoin``).

    ``centromere_bed`` is ``-c``, ``centromere_overlap_ratio`` ``-r`` (0.5), ``min_block_length`` ``-l`` (1m) and ``max_gap``
    ``-g`` (1m) accept ``1000000`` or ``1m``; ``show_non_centromeric`` is ``-e`` and ``show_long_blocks`` is ``-p``. Writes
    ``misjoins.txt``; ``metrics`` holds the counts (inter-chromosomal misjoins, gaps, candidate inversions).
    """
    arguments: list[str] = []
    _opt_file(arguments, "-c", centromere_bed, "centromere_bed")
    _opt_num(
        arguments,
        "-r",
        centromere_overlap_ratio,
        "centromere_overlap_ratio",
        minimum=0.0,
        maximum=1.0,
    )
    _opt_size(arguments, "-l", min_block_length, "min_block_length")
    _opt_size(arguments, "-g", max_gap, "max_gap")
    _switch(arguments, "-e", show_non_centromeric, "show_non_centromeric")
    _switch(arguments, "-p", show_long_blocks, "show_long_blocks")
    arguments.append(_file(paf, "paf"))
    result = _paftools(
        "misjoin",
        arguments,
        output_directory,
        "misjoins.txt",
        operation="detect_misjoins",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["metrics"] = _colon_metrics(_text(Path(result["outputs"]["output"])))
    return result


def minimap2_liftover_coordinates(
    paf: str | Path,
    query_bed: str | Path,
    output_directory: str | Path,
    *,
    min_mapping_quality: int | None = None,
    min_alignment_length: int | None = None,
    max_divergence: float | None = None,
    merge_intervals: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift BED intervals between sequences through PAF alignments (``paftools.js liftover``, a simplistic liftOver).

    The PAF needs the ``cg`` tag (``minimap2 -c``); ``query_bed`` names intervals on the PAF's query sequences. ``-q``
    ``min_mapping_quality`` (5), ``-l`` ``min_alignment_length`` (50000), ``-d`` ``max_divergence`` (>=1 disables) and
    ``-m`` ``merge_intervals``. Writes ``lifted.bed``.
    """
    arguments: list[str] = []
    _opt_int(arguments, "-q", min_mapping_quality, "min_mapping_quality", minimum=0)
    _opt_int(arguments, "-l", min_alignment_length, "min_alignment_length", minimum=0)
    _opt_num(arguments, "-d", max_divergence, "max_divergence", minimum=0.0)
    _switch(arguments, "-m", merge_intervals, "merge_intervals")
    arguments.extend((_file(paf, "paf"), _file(query_bed, "query_bed")))
    result = _paftools(
        "liftover",
        arguments,
        output_directory,
        "lifted.bed",
        operation="liftover_coordinates",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["record_count"] = _table(Path(result["outputs"]["output"]))
    return result


# ---------------------------------------------------------------------------
# paftools.js: variants and VCF
# ---------------------------------------------------------------------------


def minimap2_call_assembly_variants(
    paf: str | Path,
    output_directory: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    sample_name: str | None = None,
    min_coverage_length: int | None = None,
    min_call_length: int | None = None,
    min_mapping_quality: int | None = None,
    gap_threshold: int | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call variants from assembly-to-reference alignments with the cs tag (``paftools.js call``).

    The PAF (``minimap2 -cx asm5 --cs``) is sorted by target name and start first, as the tool requires
    (``sort -k6,6 -k8,8n``); the sorted copy is kept as ``sorted.paf``. With ``reference_fasta`` (``-f``) the calls are
    written as VCF (``variants.vcf``, sample name ``sample_name`` = ``-s``); without it ``variants.txt`` holds the tool's
    text calls. ``min_coverage_length`` is ``-l`` (10000), ``min_call_length`` ``-L`` (50000), ``min_mapping_quality`` ``-q``
    (5) and ``gap_threshold`` ``-g`` (50). The coverage and variant-count summary goes to the log and is parsed into
    ``summary``.
    """
    source = Path(_file(paf, "paf"))
    reference = None if reference_fasta is None else _file(reference_fasta, "reference_fasta")
    if sample_name is not None and reference is None:
        raise ConfigurationError("sample_name needs reference_fasta (VCF output).")
    k8, script = _resolve_paftools(paftools_path, k8_path)
    directory = _cli.new_directory(output_directory, "output_directory")
    lines = [
        line
        for line in source.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]
    parsed: list[tuple[bytes, int, str]] = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 12:
            raise ConfigurationError("paf must have at least 12 tab-separated columns per line.")
        parsed.append((fields[5].encode(), int(fields[7]), line))
    parsed.sort()
    sorted_path = directory / "sorted.paf"
    sorted_path.write_text("".join(f"{item[2]}\n" for item in parsed), encoding="utf-8")
    arguments: list[str] = []
    _opt_int(arguments, "-l", min_coverage_length, "min_coverage_length", minimum=0)
    _opt_int(arguments, "-L", min_call_length, "min_call_length", minimum=0)
    _opt_int(arguments, "-q", min_mapping_quality, "min_mapping_quality", minimum=0)
    _opt_int(arguments, "-g", gap_threshold, "gap_threshold", minimum=0)
    if reference is not None:
        arguments.extend(("-f", reference))
    _opt_text(arguments, "-s", sample_name, "sample_name")
    arguments.append(str(sorted_path))
    name = "variants.vcf" if reference is not None else "variants.txt"
    output, log = directory / name, directory / "tool.log"
    return_code = _stream(k8, [script, "call", *arguments], directory, output, log, timeout_seconds)
    if return_code != 0:
        raise BackendExecutionError(
            f"call_assembly_variants failed (exit code {return_code}).",
            code="PAFTOOLS_EXECUTION_FAILED",
            context={"log_tail": _tail(log)},
        )
    summary: dict[str, Any] = {}
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^(\d+) (.+?)(?:; ts/tv = (\S+))?$", line.strip())
        if match:
            summary[_slug(match.group(2))] = int(match.group(1))
            if match.group(3) is not None:
                summary["ts_tv"] = _number(match.group(3))
    return _result(
        "call_assembly_variants",
        "paftools.js",
        _paftools_version(k8, script, directory),
        ["k8", Path(script).name, "call", *arguments],
        directory,
        {"output": output, "sorted_paf": sorted_path, "log": log},
        log,
        summary=summary,
    )


def minimap2_pair_haplotype_vcf(
    pair_vcf: str | Path,
    output_directory: str | Path,
    *,
    male: bool = False,
    genome_version: Literal["37"] | None = None,
    sample_name: str | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the VCF records of two haplotypes into one phased diploid VCF (``paftools.js vcfpair``).

    The input has one sample column per haplotype. ``male`` is ``-m`` and requires ``genome_version='37'`` (``-g``) to treat
    the pseudo-autosomal regions of chrX; ``sample_name`` is ``-s`` (tool default ``syndip``). Writes ``diploid.vcf``.
    """
    arguments: list[str] = []
    _switch(arguments, "-m", male, "male")
    if genome_version is not None:
        arguments.extend(("-g", _cli.choice(genome_version, ("37",), "genome_version")))
    if male and genome_version is None:
        raise ConfigurationError("male=True needs genome_version='37'.")
    _opt_text(arguments, "-s", sample_name, "sample_name")
    arguments.append(_file(pair_vcf, "pair_vcf"))
    result = _paftools(
        "vcfpair",
        arguments,
        output_directory,
        "diploid.vcf",
        operation="pair_haplotype_vcf",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    with Path(result["outputs"]["output"]).open(encoding="utf-8") as handle:
        result["record_count"] = sum(
            1 for line in handle if line.strip() and not line.startswith("#")
        )
    return result


def minimap2_vcf_stats(
    vcf: str | Path,
    output_directory: str | Path,
    *,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary statistics of a VCF: substitutions, transitions/transversions and indels by size (``paftools.js vcfstat``).

    Writes ``vcf_stats.txt``; ``metrics`` holds the numbers (for example ``substitutions``, ``ts_tv``, ``insertions``).
    """
    result = _paftools(
        "vcfstat",
        [_file(vcf, "vcf")],
        output_directory,
        "vcf_stats.txt",
        operation="vcf_stats",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["metrics"] = _colon_metrics(_text(Path(result["outputs"]["output"])))
    return result


def minimap2_select_vcf_by_length(
    vcf: str | Path,
    output_directory: str | Path,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the VCF records whose variant length lies in a range (``paftools.js vcfsel``).

    ``min_length`` is ``-l`` and ``max_length`` is ``-L``; the header is kept. Writes ``selected.vcf``.
    """
    arguments: list[str] = []
    _opt_int(arguments, "-l", min_length, "min_length", minimum=0)
    _opt_int(arguments, "-L", max_length, "max_length", minimum=0)
    arguments.append(_file(vcf, "vcf"))
    result = _paftools(
        "vcfsel",
        arguments,
        output_directory,
        "selected.vcf",
        operation="select_vcf_by_length",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    with Path(result["outputs"]["output"]).open(encoding="utf-8") as handle:
        result["record_count"] = sum(
            1 for line in handle if line.strip() and not line.startswith("#")
        )
    return result


def minimap2_compare_sv_callsets(
    base_vcf: str | Path,
    call_vcf: str | Path,
    output_directory: str | Path,
    *,
    confident_regions_bed: str | Path | None = None,
    min_discard_length: int | None = None,
    min_sv_length: int | None = None,
    max_sv_length: int | None = None,
    window_size: int | None = None,
    max_allele_difference: float | None = None,
    show_errors: bool = False,
    show_matches: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two structural-variant callsets in VCF (``paftools.js sveval``).

    ``confident_regions_bed`` is ``-r``, ``min_discard_length`` ``-f`` (30), ``min_sv_length`` ``-i`` (50), ``max_sv_length``
    ``-x`` (100000), ``window_size`` ``-w`` (500), ``max_allele_difference`` ``-d`` (0.5), ``show_errors`` ``-e`` and ``show_matches`` ``-p``. Writes
    ``sveval.txt`` (the tool's report lines, in ``rows``).
    """
    arguments: list[str] = []
    _opt_file(arguments, "-r", confident_regions_bed, "confident_regions_bed")
    _opt_int(arguments, "-f", min_discard_length, "min_discard_length", minimum=0)
    _opt_int(arguments, "-i", min_sv_length, "min_sv_length", minimum=0)
    _opt_int(arguments, "-x", max_sv_length, "max_sv_length", minimum=0)
    _opt_int(arguments, "-w", window_size, "window_size", minimum=0)
    _opt_num(arguments, "-d", max_allele_difference, "max_allele_difference", minimum=0.0)
    _switch(arguments, "-e", show_errors, "show_errors")
    _switch(arguments, "-p", show_matches, "show_matches")
    arguments.extend((_file(base_vcf, "base_vcf"), _file(call_vcf, "call_vcf")))
    result = _paftools(
        "sveval",
        arguments,
        output_directory,
        "sveval.txt",
        operation="compare_sv_callsets",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["row_count"] = _table(Path(result["outputs"]["output"]))
    return result


# ---------------------------------------------------------------------------
# paftools.js: mapping and annotation evaluation
# ---------------------------------------------------------------------------


def minimap2_evaluate_mapping(
    alignments: str | Path,
    output_directory: str | Path,
    *,
    min_overlap_ratio: float | None = None,
    print_wrong_mapq: int | None = None,
    mode: Literal[0, 1, 2] | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate mapping accuracy against the truth encoded in simulated read names (``paftools.js mapeval``).

    Read names must be ``name!chrom!start!end!strand`` (PBSIM/Badread via ``minimap2_convert_simulated_reads``) or the
    mason2 paired-end form. ``min_overlap_ratio`` is ``-r`` (0.1), ``print_wrong_mapq`` ``-Q`` prints wrong mappings with at
    least that mapping quality, ``mode`` is ``-m`` (0 longest alignment only, 1 first only, 2 all primary). Writes
    ``mapeval.tsv``; ``records`` has one row per mapping-quality threshold (``mapq``, ``reads`` and ``wrong`` at that
    threshold, ``cumulative_error_rate`` and ``cumulative_reads`` for mapQ at least that value) and ``unmapped`` is the
    number of unmapped reads.
    """
    arguments: list[str] = []
    _opt_num(arguments, "-r", min_overlap_ratio, "min_overlap_ratio", minimum=0.0, maximum=1.0)
    _opt_int(arguments, "-Q", print_wrong_mapq, "print_wrong_mapq", minimum=0)
    if mode is not None:
        arguments.extend(("-m", str(_cli.integer(mode, "mode", minimum=0, maximum=2))))
    arguments.append(_file(alignments, "alignments"))
    result = _paftools(
        "mapeval",
        arguments,
        output_directory,
        "mapeval.tsv",
        operation="evaluate_mapping",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    rows, result["row_count"] = _table(Path(result["outputs"]["output"]))
    result["rows"] = rows
    result["records"] = [
        {
            "mapq": int(row[1]),
            "reads": int(row[2]),
            "wrong": int(row[3]),
            "cumulative_error_rate": float(row[4]),
            "cumulative_reads": int(row[5]),
        }
        for row in rows
        if row[0] == "Q" and len(row) >= 6
    ]
    unmapped = [int(row[1]) for row in rows if row[0] == "U" and len(row) >= 2]
    result["unmapped"] = unmapped[0] if unmapped else None
    return result


def minimap2_compare_paf(
    base_paf: str | Path,
    test_paf: str | Path,
    output_directory: str | Path,
    *,
    min_mapping_quality: int | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two PAF alignment files base by base (``paftools.js pafcmp``).

    ``min_mapping_quality`` is ``-q`` (10). Writes ``pafcmp.txt``; ``metrics`` counts the base alignments that are correct,
    wrong, missing or additional in ``test_paf`` relative to ``base_paf``.
    """
    arguments: list[str] = []
    _opt_int(arguments, "-q", min_mapping_quality, "min_mapping_quality", minimum=0)
    arguments.extend((_file(base_paf, "base_paf"), _file(test_paf, "test_paf")))
    result = _paftools(
        "pafcmp",
        arguments,
        output_directory,
        "pafcmp.txt",
        operation="compare_paf",
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    metrics: dict[str, int] = {}
    for line in _text(Path(result["outputs"]["output"])).splitlines():
        match = re.match(r"^X\t(\d+) (.+)$", line)
        if match:
            metrics[_slug(match.group(2))] = int(match.group(1))
    result["metrics"] = metrics
    return result


def _evaluate_annotation(
    subcommand: str,
    operation: str,
    annotation: str | Path,
    alignments: str | Path,
    output_directory: str | Path,
    flags: list[str],
    paftools_path: str | None,
    k8_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    arguments = [*flags, _file(annotation, "annotation"), _file(alignments, "alignments")]
    result = _paftools(
        subcommand,
        arguments,
        output_directory,
        f"{subcommand}.txt",
        operation=operation,
        paftools_path=paftools_path,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"], result["row_count"] = _table(Path(result["outputs"]["output"]))
    result["report_text"] = _text(Path(result["outputs"]["output"]))
    result["metrics"] = _colon_metrics(result["report_text"])
    return result


def minimap2_evaluate_junctions(
    annotation_gtf: str | Path,
    alignments: str | Path,
    output_directory: str | Path,
    *,
    input_format: Literal["sam", "miniprot_paf", "bed"] = "sam",
    tolerance: int | None = None,
    standard_chromosomes_only: bool = False,
    first_alignment_only: bool = False,
    print_overlapping: bool = False,
    print_errors: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate splice junctions of spliced alignments against a gene annotation (``paftools.js junceval``).

    ``alignments`` is SAM by default; ``input_format='miniprot_paf'`` (``-a``) or ``'bed'`` (``-b``) change the reader.
    ``tolerance`` is ``-l`` (junction position tolerance, 0 = exact), ``standard_chromosomes_only`` is ``-c``,
    ``first_alignment_only`` is ``-1``, ``print_overlapping`` is ``-p`` and ``print_errors`` is ``-e``. Writes
    ``junceval.txt`` (report in ``report_text`` and ``rows``).
    """
    _cli.choice(input_format, ("sam", "miniprot_paf", "bed"), "input_format")
    flags: list[str] = []
    _opt_int(flags, "-l", tolerance, "tolerance", minimum=0)
    _switch(flags, "-p", print_overlapping, "print_overlapping")
    _switch(flags, "-e", print_errors, "print_errors")
    _switch(flags, "-c", standard_chromosomes_only, "standard_chromosomes_only")
    _switch(flags, "-1", first_alignment_only, "first_alignment_only")
    if input_format == "miniprot_paf":
        flags.append("-a")
    elif input_format == "bed":
        flags.append("-b")
    return _evaluate_annotation(
        "junceval",
        "evaluate_junctions",
        annotation_gtf,
        alignments,
        output_directory,
        flags,
        paftools_path,
        k8_path,
        timeout_seconds,
    )


def minimap2_evaluate_exons(
    annotation_gtf: str | Path,
    alignments: str | Path,
    output_directory: str | Path,
    *,
    input_format: Literal["sam", "miniprot_paf", "bed"] = "sam",
    tolerance: int | None = None,
    coding_only: bool = False,
    standard_chromosomes_only: bool = False,
    first_alignment_only: bool = False,
    skip_first_exon: bool = False,
    skip_terminal_exons: bool = False,
    base_level: bool = False,
    print_overlapping: bool = False,
    print_errors: bool = False,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate exons of spliced alignments against a gene annotation (``paftools.js exoneval``).

    Options as for ``minimap2_evaluate_junctions`` plus ``coding_only`` (``-d``, coding regions instead of exons),
    ``skip_first_exon`` (``-f``, miniprot mode), ``skip_terminal_exons`` (``-t``) and ``base_level`` (``-s``, base-level
    sensitivity and specificity, more memory). Writes ``exoneval.txt``.
    """
    _cli.choice(input_format, ("sam", "miniprot_paf", "bed"), "input_format")
    flags: list[str] = []
    _opt_int(flags, "-l", tolerance, "tolerance", minimum=0)
    _switch(flags, "-d", coding_only, "coding_only")
    _switch(flags, "-p", print_overlapping, "print_overlapping")
    _switch(flags, "-e", print_errors, "print_errors")
    _switch(flags, "-c", standard_chromosomes_only, "standard_chromosomes_only")
    _switch(flags, "-1", first_alignment_only, "first_alignment_only")
    _switch(flags, "-f", skip_first_exon, "skip_first_exon")
    _switch(flags, "-t", skip_terminal_exons, "skip_terminal_exons")
    _switch(flags, "-s", base_level, "base_level")
    if input_format == "miniprot_paf":
        flags.append("-a")
    elif input_format == "bed":
        flags.append("-b")
    return _evaluate_annotation(
        "exoneval",
        "evaluate_exons",
        annotation_gtf,
        alignments,
        output_directory,
        flags,
        paftools_path,
        k8_path,
        timeout_seconds,
    )


def minimap2_evaluate_overlaps(
    read_to_reference_paf: str | Path,
    overlap_paf: str | Path,
    output_directory: str | Path,
    *,
    min_overlap_length: int | None = None,
    min_mapping_quality: int | None = None,
    min_mapped_fraction: float | None = None,
    paftools_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate read-overlap sensitivity against the mapping of the reads to a reference (``paftools.js ov-eval``).

    ``read_to_reference_paf`` is sorted by target name and start first (``sort -k6,6 -k8,8n``; the sorted copy is kept as
    ``sorted_to_ref.paf``); ``overlap_paf`` is the all-vs-all overlap PAF (``minimap2_map_reads`` with ``ava-pb`` or
    ``ava-ont``). ``min_overlap_length`` is ``-l`` (2000), ``min_mapping_quality`` ``-q`` (10) and ``min_mapped_fraction``
    ``-f`` (0.95). Writes ``ov_eval.txt`` (``rows`` and ``report_text``).
    """
    source = Path(_file(read_to_reference_paf, "read_to_reference_paf"))
    overlaps = _file(overlap_paf, "overlap_paf")
    k8, script = _resolve_paftools(paftools_path, k8_path)
    directory = _cli.new_directory(output_directory, "output_directory")
    parsed: list[tuple[bytes, int, str]] = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            raise ConfigurationError(
                "read_to_reference_paf must have at least 12 tab-separated columns per line."
            )
        parsed.append((fields[5].encode(), int(fields[7]), line))
    parsed.sort()
    sorted_path = directory / "sorted_to_ref.paf"
    sorted_path.write_text("".join(f"{item[2]}\n" for item in parsed), encoding="utf-8")
    arguments: list[str] = []
    _opt_int(arguments, "-l", min_overlap_length, "min_overlap_length", minimum=0)
    _opt_int(arguments, "-q", min_mapping_quality, "min_mapping_quality", minimum=0)
    _opt_num(arguments, "-f", min_mapped_fraction, "min_mapped_fraction", minimum=0.0, maximum=1.0)
    arguments.extend((str(sorted_path), overlaps))
    output, log = directory / "ov_eval.txt", directory / "tool.log"
    return_code = _stream(
        k8, [script, "ov-eval", *arguments], directory, output, log, timeout_seconds
    )
    if return_code != 0:
        raise BackendExecutionError(
            f"evaluate_overlaps failed (exit code {return_code}).",
            code="PAFTOOLS_EXECUTION_FAILED",
            context={"log_tail": _tail(log)},
        )
    rows, total = _table(output)
    report = _text(output)
    metrics: dict[str, Any] = {}
    for pattern, key in (
        (r"(\d+) overlaps inferred", "inferred_overlaps"),
        (r"(\d+) missed", "missed_overlaps"),
        (r"([\d.]+)% sensitivity", "sensitivity_percent"),
    ):
        match = re.search(pattern, report)
        if match:
            metrics[key] = _number(match.group(1))
    return _result(
        "evaluate_overlaps",
        "paftools.js",
        _paftools_version(k8, script, directory),
        ["k8", Path(script).name, "ov-eval", *arguments],
        directory,
        {"output": output, "sorted_to_ref_paf": sorted_path, "log": log},
        log,
        rows=rows,
        row_count=total,
        report_text=report,
        metrics=metrics,
    )


def minimap2_cluster_overlap_reads(
    overlap_paf: str | Path,
    output_directory: str | Path,
    *,
    min_coverage: float | None = None,
    max_divergence: float | None = None,
    max_difference: float | None = None,
    pafcluster_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster reads from an all-vs-all overlap PAF into read groups (``pafcluster.js`` from minimap2's ``misc/`` folder).

    ``overlap_paf`` comes from ``minimap2_map_reads`` with ``ava-pb`` or ``ava-ont``. ``min_coverage`` is ``-c`` (0.9),
    ``max_divergence`` ``-d`` (0.015) and ``max_difference`` ``-e`` (20000). Writes ``clusters.txt`` in the tool's own record format
    (``SD`` seed, ``CL`` cluster length and ``//`` terminated blocks; also in ``rows``).
    """
    script = _cli.program("pafcluster.js", explicit=pafcluster_path, hint=_PAFCLUSTER_HINT)
    arguments: list[str] = []
    _opt_num(arguments, "-c", min_coverage, "min_coverage", minimum=0.0, maximum=1.0)
    _opt_num(arguments, "-d", max_divergence, "max_divergence", minimum=0.0)
    _opt_num(arguments, "-e", max_difference, "max_difference", minimum=0.0)
    arguments.append(_file(overlap_paf, "overlap_paf"))
    result = _paftools(
        "",
        arguments,
        output_directory,
        "clusters.txt",
        operation="cluster_overlap_reads",
        paftools_path=None,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        script=script,
    )
    result["rows"], result["row_count"] = _table(Path(result["outputs"]["output"]))
    return result


__all__ = [
    "minimap2_alignment_stats",
    "minimap2_annotation_to_bed12",
    "minimap2_assembly_gene_completeness",
    "minimap2_assembly_stats",
    "minimap2_call_assembly_variants",
    "minimap2_cluster_overlap_reads",
    "minimap2_compare_paf",
    "minimap2_compare_sv_callsets",
    "minimap2_convert_simulated_reads",
    "minimap2_delta_to_paf",
    "minimap2_detect_misjoins",
    "minimap2_evaluate_exons",
    "minimap2_evaluate_junctions",
    "minimap2_evaluate_mapping",
    "minimap2_evaluate_overlaps",
    "minimap2_gff3_to_junctions",
    "minimap2_liftover_coordinates",
    "minimap2_longcs_to_sequences",
    "minimap2_map_reads",
    "minimap2_paf_to_gff3",
    "minimap2_pair_haplotype_vcf",
    "minimap2_sam_to_paf",
    "minimap2_select_vcf_by_length",
    "minimap2_simulated_names_to_bed",
    "minimap2_spliced_alignments_to_bed12",
    "minimap2_vcf_stats",
    "minimap2_view_paf",
]

_register_public_verbose(globals(), tuple(__all__))
