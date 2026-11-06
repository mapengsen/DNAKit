"""BWA command-line facades: index build stages, backtrack alignment (aln, samse, sampe), SMEM seeds, paired-end
merging, MEM-length histograms and the bwa.kit scripts (ALT post-processing, HLA typing, reference generation).

The user installs BWA (``conda install -c bioconda bwa``) and, for the ``bwakit_*`` functions, bwa.kit
(``conda install -c bioconda bwakit``, or the unpacked release tarball); BWA 0.7.19 and bwa.kit 0.7.15 are the tested
versions. Each function runs one fixed command with an explicit typed signature mapped to fixed flags; no shell is
used and the program is resolved from an explicit path, then PATH. The bwa.kit directory is resolved from an explicit
path, then from ``run-bwamem`` on PATH (the conda package keeps the scripts in ``share/bwakit-*``). Outputs are **new**
files (an existing file is never overwritten). Results are JSON-serialisable dictionaries: ``status``, ``backend``,
``version``, ``command``, ``return_code``, ``outputs``, ``log_tail`` plus small parsed summaries.

``bwa mem`` is not repeated here; these functions cover the other BWA subcommands and the bwa.kit scripts.
"""

from __future__ import annotations

import re
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p13_cli as _p13
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_PROGRAM = "bwa"
_HINT = (
    "Install BWA 0.7.x (for example `conda install -c bioconda bwa`, or build it from "
    "https://github.com/lh3/bwa), then put it on PATH or pass bwa_path."
)
_KIT_HINT = (
    "Install bwa.kit (for example `conda install -c bioconda bwakit`, or unpack the bwakit release tarball from "
    "https://sourceforge.net/projects/bio-bwa/files/bwakit/) and pass bwakit_directory, or put run-bwamem on PATH."
)
_VERSION_FLAGS: tuple[
    str, ...
] = ()  # bwa prints its usage, with the version, when called without arguments
_STAGES = ("fa2pac", "pac2bwt", "pac2bwtgen", "bwtupdate", "bwt2sa")
_REFERENCE_SETS = ("hs38", "hs38a", "hs38DH", "hs37", "hs37d5")
_GENE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check(bwa_path: str | None) -> None:
    _cli.program(_PROGRAM, explicit=bwa_path, hint=_HINT)


def _run(
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    bwa_path: str | None,
    timeout_seconds: float,
    operation: str,
) -> dict[str, Any]:
    # bwa has no --version flag: its usage text (printed without arguments) carries the version.
    result = _p13.run(
        _PROGRAM,
        arguments,
        cwd=cwd,
        expected=expected,
        explicit=bwa_path,
        hint=_HINT,
        version_flags=_VERSION_FLAGS,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = operation
    return result


def _stream(
    arguments: Sequence[str],
    output: Path,
    bwa_path: str | None,
    timeout_seconds: float,
    operation: str,
) -> dict[str, Any]:
    result = _p13.run_to_file(
        _PROGRAM,
        arguments,
        output=output,
        explicit=bwa_path,
        hint=_HINT,
        version_flags=_VERSION_FLAGS,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = operation
    return result


def _index_prefix(value: str, name: str = "index_prefix") -> Path:
    """A BWA index given by its prefix: ``<prefix>.bwt`` must exist."""
    prefix = Path(_cli.text(value, name)).expanduser().resolve()
    if not Path(f"{prefix}.bwt").is_file():
        raise ConfigurationError(
            f"{name} must be a BWA index prefix with an existing <prefix>.bwt: {prefix}"
        )
    return prefix


def _optional_int(flag: str, value: int | None, name: str, *, minimum: int = 0) -> list[str]:
    return [] if value is None else [flag, str(_cli.integer(value, name, minimum=minimum))]


def _switch(on: bool, flag: str, name: str) -> list[str]:
    if not isinstance(on, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return [flag] if on else []


def _read_group(value: str | None) -> list[str]:
    if value is None:
        return []
    text = _cli.text(value, "read_group")
    if not text.startswith("@RG"):
        raise ConfigurationError(
            "read_group must be a read group header line such as '@RG\\tID:foo\\tSM:bar'."
        )
    return ["-r", text]


# ---------------------------------------------------------------------------
# bwa index build stages
# ---------------------------------------------------------------------------


def bwa_index_stage(
    stage: str,
    input_path: str,
    output_path: str,
    *,
    forward_only: bool = False,
    sa_interval: int = 32,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one stage of BWA's index construction by hand (``bwa fa2pac|pac2bwt|pac2bwtgen|bwtupdate|bwt2sa``).

    ``bwa index`` chains these stages; running them separately gives access to the intermediate files.
    ``stage="fa2pac"`` packs a FASTA file into 2-bit ``<output_path>.pac`` plus ``.ann`` and ``.amb`` (``output_path``
    is the output prefix; ``forward_only`` is ``-f``). ``stage="pac2bwt"`` builds the BWT from the packed reference
    with BWA's default algorithm and ``stage="pac2bwtgen"`` with the alternative algorithm (``bwtgen``); the result is in
    the old format and needs ``bwtupdate`` before BWA can use it. ``stage="bwtupdate"`` converts a ``.bwt`` file to the
    current format (``input_path`` is left untouched; the converted copy is ``output_path``). ``stage="bwt2sa"`` builds
    the suffix array ``output_path`` from a current-format BWT (``sa_interval`` is ``-i``, default 32).
    """
    kind = _cli.choice(stage, _STAGES, "stage")
    source = _cli.existing_file(input_path, "input_path")
    interval = _cli.integer(sa_interval, "sa_interval", minimum=1)
    if kind != "fa2pac" and forward_only:
        raise ConfigurationError("forward_only only applies to stage='fa2pac'.")
    _check(bwa_path)
    if kind == "fa2pac":
        prefix = Path(output_path).expanduser().resolve()
        if not prefix.parent.is_dir():
            raise ConfigurationError(f"The folder for output_path does not exist: {prefix.parent}")
        targets = [Path(f"{prefix}{suffix}") for suffix in (".pac", ".ann", ".amb")]
        clash = [str(path) for path in targets if path.exists()]
        if clash:
            raise ConfigurationError(
                f"Output files already exist and will not be overwritten: {clash}"
            )
        arguments = [
            "fa2pac",
            *_switch(forward_only, "-f", "forward_only"),
            str(source),
            str(prefix),
        ]
        return _run(
            arguments,
            cwd=prefix.parent,
            expected=targets,
            bwa_path=bwa_path,
            timeout_seconds=timeout_seconds,
            operation="bwa fa2pac",
        )
    target = _cli.new_file(output_path, "output_path")
    if kind == "bwtupdate":
        shutil.copyfile(source, target)
        result = _run(
            ["bwtupdate", str(target)],
            cwd=target.parent,
            expected=[target],
            bwa_path=bwa_path,
            timeout_seconds=timeout_seconds,
            operation="bwa bwtupdate",
        )
        return result
    arguments = {
        "pac2bwt": ["pac2bwt", str(source), str(target)],
        "pac2bwtgen": ["pac2bwtgen", str(source), str(target)],
        "bwt2sa": ["bwt2sa", "-i", str(interval), str(source), str(target)],
    }[kind]
    return _run(
        arguments,
        cwd=target.parent,
        expected=[target],
        bwa_path=bwa_path,
        timeout_seconds=timeout_seconds,
        operation=f"bwa {kind}",
    )


# ---------------------------------------------------------------------------
# aln, samse, sampe
# ---------------------------------------------------------------------------


def bwa_aln(
    index_prefix: str,
    reads_path: str,
    output_path: str,
    *,
    max_diff: float = 0.04,
    max_gap_opens: float = 1,
    max_gap_extensions: int = -1,
    indel_end_skip: int = 5,
    max_long_deletion_occurrences: int = 10,
    seed_length: int = 32,
    max_seed_diff: int = 2,
    max_queue_entries: int = 2_000_000,
    threads: int = 1,
    mismatch_penalty: int = 3,
    gap_open_penalty: int = 11,
    gap_extension_penalty: int = 4,
    max_best_hits: int = 30,
    trim_quality: int = 0,
    barcode_length: int | None = None,
    log_scaled_gap_penalty: bool = False,
    non_iterative: bool = False,
    illumina_1_3: bool = False,
    bam_input: bool = False,
    bam_read: str | None = None,
    filter_casava: bool = False,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """BWA-backtrack alignment of short reads, writing suffix-array intervals (``bwa aln``).

    The defaults are those of ``bwa aln``: ``max_diff`` (``-n``, an integer number of differences or a fraction for the
    missing probability), ``max_gap_opens`` (``-o``), ``max_gap_extensions`` (``-e``, -1 disables long gaps),
    ``indel_end_skip`` (``-i``), ``max_long_deletion_occurrences`` (``-d``), ``seed_length`` (``-l``),
    ``max_seed_diff`` (``-k``), ``max_queue_entries`` (``-m``), ``threads`` (``-t``), the penalties (``-M``, ``-O``,
    ``-E``), ``max_best_hits`` (``-R``), ``trim_quality`` (``-q``), ``barcode_length`` (``-B``),
    ``log_scaled_gap_penalty`` (``-L``), ``non_iterative`` (``-N``), ``illumina_1_3`` (``-I``), ``bam_input`` (``-b``;
    ``bam_read`` then selects ``"single"`` (``-0``), ``"first"`` (``-1``) or ``"second"`` (``-2``) reads) and
    ``filter_casava`` (``-Y``). The ``.sai`` intervals in ``output_path`` are turned into SAM by
    :func:`bwa_sai_to_sam`. ``index_prefix`` is a ``bwa index`` prefix.
    """
    prefix = _index_prefix(index_prefix)
    reads = _cli.existing_file(reads_path, "reads_path")
    args = ["aln", "-n", f"{_cli.number(max_diff, 'max_diff', minimum=0):g}"]
    args += ["-o", f"{_cli.number(max_gap_opens, 'max_gap_opens', minimum=0):g}"]
    args += ["-e", str(_cli.integer(max_gap_extensions, "max_gap_extensions", minimum=-1))]
    args += ["-i", str(_cli.integer(indel_end_skip, "indel_end_skip", minimum=0))]
    args += [
        "-d",
        str(
            _cli.integer(max_long_deletion_occurrences, "max_long_deletion_occurrences", minimum=0)
        ),
    ]
    args += ["-l", str(_cli.integer(seed_length, "seed_length", minimum=1))]
    args += ["-k", str(_cli.integer(max_seed_diff, "max_seed_diff", minimum=0))]
    args += ["-m", str(_cli.integer(max_queue_entries, "max_queue_entries", minimum=1))]
    args += ["-t", str(_cli.integer(threads, "threads", minimum=1, maximum=256))]
    args += ["-M", str(_cli.integer(mismatch_penalty, "mismatch_penalty", minimum=0))]
    args += ["-O", str(_cli.integer(gap_open_penalty, "gap_open_penalty", minimum=0))]
    args += ["-E", str(_cli.integer(gap_extension_penalty, "gap_extension_penalty", minimum=0))]
    args += ["-R", str(_cli.integer(max_best_hits, "max_best_hits", minimum=0))]
    args += ["-q", str(_cli.integer(trim_quality, "trim_quality", minimum=0))]
    args += _optional_int("-B", barcode_length, "barcode_length", minimum=1)
    args += _switch(log_scaled_gap_penalty, "-L", "log_scaled_gap_penalty")
    args += _switch(non_iterative, "-N", "non_iterative")
    args += _switch(illumina_1_3, "-I", "illumina_1_3")
    args += _switch(filter_casava, "-Y", "filter_casava")
    if bam_read is not None and not bam_input:
        raise ConfigurationError("bam_read only applies together with bam_input=True.")
    if bam_input:
        args.append("-b")
        if bam_read is not None:
            args.append(
                {"single": "-0", "first": "-1", "second": "-2"}[
                    _cli.choice(bam_read, ("single", "first", "second"), "bam_read")
                ]
            )
    _check(bwa_path)
    out = _cli.new_file(output_path, "output_path")
    args += ["-f", str(out), str(prefix), str(reads)]
    result = _run(
        args,
        cwd=out.parent,
        expected=[out],
        bwa_path=bwa_path,
        timeout_seconds=timeout_seconds,
        operation="bwa aln",
    )
    return result


def bwa_sai_to_sam(
    index_prefix: str,
    sai_paths: Sequence[str],
    reads_paths: Sequence[str],
    output_path: str,
    *,
    read_group: str | None = None,
    max_hits: int | None = None,
    max_insert_size: int | None = None,
    max_occurrences_per_end: int | None = None,
    max_discordant_hits: int | None = None,
    chimeric_rate: float | None = None,
    preload_index: bool = False,
    disable_smith_waterman: bool = False,
    disable_insert_size_estimate: bool = False,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate SAM from the intervals of :func:`bwa_aln`: single-end (``bwa samse``) or paired-end (``bwa sampe``).

    One ``.sai`` and one read file give single-end output (``samse``); two of each, in the order of read 1 and read 2,
    give paired-end output (``sampe``). ``read_group`` is a header line such as ``@RG\\tID:foo\\tSM:bar`` (``-r``).
    ``max_hits`` is ``-n``: the maximum number of alternative hits written in the XA tag for single-end reads (BWA's
    default 3) or the maximum hits output for paired reads (default 3). The paired-end only options are
    ``max_insert_size`` (``-a``, default 500), ``max_occurrences_per_end`` (``-o``, default 100000),
    ``max_discordant_hits`` (``-N``, default 10), ``chimeric_rate`` (``-c``, default 1e-5), ``preload_index`` (``-P``),
    ``disable_smith_waterman`` (``-s``) and ``disable_insert_size_estimate`` (``-A``).
    """
    prefix = _index_prefix(index_prefix)
    if isinstance(sai_paths, str) or isinstance(reads_paths, str):
        raise ConfigurationError("sai_paths and reads_paths must be lists of files.")
    sai = [_cli.existing_file(path, "sai_paths") for path in sai_paths]
    reads = [_cli.existing_file(path, "reads_paths") for path in reads_paths]
    if len(sai) != len(reads) or len(sai) not in (1, 2):
        raise ConfigurationError(
            "Give one .sai and one read file (single-end) or two of each (paired-end)."
        )
    paired = len(sai) == 2
    paired_only = {
        "max_insert_size": max_insert_size,
        "max_occurrences_per_end": max_occurrences_per_end,
        "max_discordant_hits": max_discordant_hits,
        "chimeric_rate": chimeric_rate,
    }
    if not paired and (
        any(value is not None for value in paired_only.values())
        or preload_index
        or disable_smith_waterman
        or disable_insert_size_estimate
    ):
        raise ConfigurationError(
            "max_insert_size, max_occurrences_per_end, max_discordant_hits, chimeric_rate, preload_index, "
            "disable_smith_waterman and disable_insert_size_estimate only apply to paired-end input."
        )
    args = ["sampe" if paired else "samse"]
    args += _optional_int("-n", max_hits, "max_hits", minimum=1)
    args += _read_group(read_group)
    if paired:
        args += _optional_int("-a", max_insert_size, "max_insert_size", minimum=1)
        args += _optional_int("-o", max_occurrences_per_end, "max_occurrences_per_end", minimum=1)
        args += _optional_int("-N", max_discordant_hits, "max_discordant_hits", minimum=1)
        if chimeric_rate is not None:
            args += ["-c", f"{_cli.number(chimeric_rate, 'chimeric_rate', minimum=0.0):g}"]
        args += _switch(preload_index, "-P", "preload_index")
        args += _switch(disable_smith_waterman, "-s", "disable_smith_waterman")
        args += _switch(disable_insert_size_estimate, "-A", "disable_insert_size_estimate")
    _check(bwa_path)
    out = _cli.new_file(output_path, "output_path")
    args += [
        "-f",
        str(out),
        str(prefix),
        *[str(path) for path in sai],
        *[str(path) for path in reads],
    ]
    result = _run(
        args,
        cwd=out.parent,
        expected=[out],
        bwa_path=bwa_path,
        timeout_seconds=timeout_seconds,
        operation="bwa sampe" if paired else "bwa samse",
    )
    result["paired"] = paired
    result["alignment_count"] = _p13.count_lines(out, skip_prefixes=("@",)) if out.is_file() else 0
    return result


# ---------------------------------------------------------------------------
# fastmap, pemerge, maxk
# ---------------------------------------------------------------------------


def _parse_fastmap(path: Path, limit: int) -> tuple[list[dict[str, Any]], int]:
    queries: list[dict[str, Any]] = []
    total = 0
    current: dict[str, Any] | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if fields[0] == "SQ" and len(fields) >= 3:
            total += 1
            current = (
                {"query": fields[1], "length": int(fields[2]), "matches": []}
                if total <= limit
                else None
            )
            if current is not None:
                queries.append(current)
        elif fields[0] == "EM" and current is not None and len(fields) >= 4:
            current["matches"].append(
                {
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "hits": int(fields[3]),
                    "positions": [position for position in fields[4:] if position],
                }
            )
    return queries, total


def bwa_fastmap(
    index_prefix: str,
    reads_path: str,
    output_path: str,
    *,
    min_smem_length: int = 17,
    max_interval_size: int = 20,
    min_interval_size: int = 1,
    max_mem_length: int = 2_147_483_647,
    stop_interval_size: int = 0,
    max_queries_returned: int = 1000,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the super-maximal exact matches (SMEMs) of sequences against a BWA index (``bwa fastmap``).

    ``min_smem_length`` is ``-l``, ``max_interval_size`` is ``-w`` (coordinates are listed only for matches with at most
    this many hits), ``min_interval_size`` is ``-i``, ``max_mem_length`` is ``-L`` and ``stop_interval_size`` is ``-I``.
    ``output_path`` receives BWA's text (``SQ`` line per query, then ``EM`` lines ``start end hits positions`` and
    ``//``); the first ``max_queries_returned`` queries are also returned parsed in ``queries``.
    """
    prefix = _index_prefix(index_prefix)
    reads = _cli.existing_file(reads_path, "reads_path")
    limit = _cli.integer(max_queries_returned, "max_queries_returned", minimum=1)
    args = ["fastmap", "-l", str(_cli.integer(min_smem_length, "min_smem_length", minimum=1))]
    args += ["-w", str(_cli.integer(max_interval_size, "max_interval_size", minimum=1))]
    args += ["-i", str(_cli.integer(min_interval_size, "min_interval_size", minimum=1))]
    args += ["-L", str(_cli.integer(max_mem_length, "max_mem_length", minimum=1))]
    args += ["-I", str(_cli.integer(stop_interval_size, "stop_interval_size", minimum=0))]
    args += [str(prefix), str(reads)]
    _check(bwa_path)
    out = _cli.new_file(output_path, "output_path")
    result = _stream(args, out, bwa_path, timeout_seconds, "bwa fastmap")
    queries, total = _parse_fastmap(out, limit) if out.is_file() else ([], 0)
    result["query_count"] = total
    result["queries"] = queries
    return result


def bwa_pemerge(
    read1_path: str,
    output_path: str,
    *,
    read2_path: str | None = None,
    output: str = "all",
    threads: int = 1,
    min_overlap: int = 10,
    max_error_sum: int = 70,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge overlapping paired-end reads into single reads (``bwa pemerge``, experimental in BWA).

    ``read1_path`` and ``read2_path`` are the two FASTQ files, or ``read1_path`` alone is an interleaved FASTQ.
    ``output`` is ``"all"`` (merged reads plus the pairs that could not be merged, default), ``"merged"`` (``-m``) or
    ``"unmerged"`` (``-u``). ``min_overlap`` is ``-T`` and ``max_error_sum`` is ``-Q``. The FASTQ is written to
    ``output_path``; merged reads carry ``merged`` in their name line.
    """
    read1 = _cli.existing_file(read1_path, "read1_path")
    mode = _cli.choice(output, ("all", "merged", "unmerged"), "output")
    args = ["pemerge", "-t", str(_cli.integer(threads, "threads", minimum=1, maximum=256))]
    args += ["-T", str(_cli.integer(min_overlap, "min_overlap", minimum=1))]
    args += ["-Q", str(_cli.integer(max_error_sum, "max_error_sum", minimum=0))]
    args += {"all": [], "merged": ["-m"], "unmerged": ["-u"]}[mode]
    args.append(str(read1))
    if read2_path is not None:
        args.append(str(_cli.existing_file(read2_path, "read2_path")))
    _check(bwa_path)
    out = _cli.new_file(output_path, "output_path")
    result = _stream(args, out, bwa_path, timeout_seconds, "bwa pemerge")
    result["merged_pairs"] = (
        len(re.findall(r"\bmerged\b", " ".join(_header_lines(out)))) if out.is_file() else 0
    )
    result["record_count"] = _p13.count_lines(out) // 4 if out.is_file() else 0
    return result


def _header_lines(path: Path) -> list[str]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for index, line in enumerate(handle) if index % 4 == 0]


def bwa_maxk(
    bwt_path: str,
    sequences_path: str,
    output_path: str,
    *,
    self_match: bool = False,
    bwa_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of maximal-exact-match lengths of sequences against an index (``bwa maxk``).

    For every query position the longest super-maximal exact match covering it is taken; the output has 256 lines
    ``length<TAB>number of positions`` (matches of 255 or more bases are capped at 255). ``bwt_path`` is the ``.bwt``
    file of a ``bwa index`` (BWA opens that file, not the prefix). ``self_match`` is ``-s``: matches are required to
    occur at least twice in the index, which is how a sequence that is itself indexed is compared with the rest of
    the index. ``histogram`` lists the lengths with a non-zero count.
    """
    bwt = _cli.existing_file(bwt_path, "bwt_path")
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    args = ["maxk", *_switch(self_match, "-s", "self_match"), str(bwt), str(sequences)]
    _check(bwa_path)
    out = _cli.new_file(output_path, "output_path")
    result = _stream(args, out, bwa_path, timeout_seconds, "bwa maxk")
    histogram: list[dict[str, int]] = []
    if out.is_file():
        for line in out.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split("\t")
            if (
                len(fields) == 2
                and fields[0].isdigit()
                and fields[1].isdigit()
                and int(fields[1]) > 0
            ):
                histogram.append({"length": int(fields[0]), "positions": int(fields[1])})
    result["histogram"] = histogram
    return result


# ---------------------------------------------------------------------------
# bwa.kit scripts
# ---------------------------------------------------------------------------


def _bwakit_directory(explicit: str | None, required: Sequence[str]) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit).expanduser())
    else:
        found = shutil.which("run-bwamem")
        if found is not None:
            real = Path(found).resolve()
            candidates.append(real.parent)
            candidates += sorted((real.parent.parent / "share").glob("bwakit-*"))
    for directory in candidates:
        if directory.is_dir() and all((directory / name).exists() for name in required):
            return directory.resolve()
    raise BackendUnavailableError(
        f"The bwa.kit directory with {', '.join(required)} was not found.", hint=_KIT_HINT
    )


def _kit_version(kit: Path) -> str | None:
    """bwa.kit version from its directory name (``bwakit-0.7.15-1``) or the top of its README."""
    match = re.search(r"\d+(?:\.\d+)+", kit.name)
    return match.group(0) if match is not None else None


def bwakit_postalt(
    alt_path: str,
    alignments_path: str,
    output_path: str,
    *,
    hla_prefix: str | None = None,
    min_alt_fraction_for_mapq: float = 1.0,
    bwakit_directory: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Post-process ``bwa mem`` alignments for ALT contigs, decoys and HLA genes (bwa.kit ``bwa-postalt.js``).

    ``alt_path`` is the ``<ref>.alt`` file (ALT contig to primary assembly alignments in SAM) written next to a
    reference prepared with :func:`bwakit_generate_reference`; ``alignments_path`` is the SAM from ``bwa mem`` on
    that reference. The script lifts the ALT hits to the primary assembly, groups them and re-estimates the mapping
    quality across groups; the new SAM goes to ``output_path``. ``hla_prefix`` additionally writes the reads that
    match the HLA genes to ``<hla_prefix>.HLA-A.fq`` and so on (``-p``); ``min_alt_fraction_for_mapq`` is ``-r``.
    """
    alt = _cli.existing_file(alt_path, "alt_path")
    alignments = _cli.existing_file(alignments_path, "alignments_path")
    fraction = _cli.number(
        min_alt_fraction_for_mapq, "min_alt_fraction_for_mapq", minimum=0.0, maximum=1.0
    )
    kit = _bwakit_directory(bwakit_directory, ("k8", "bwa-postalt.js"))
    args = [str(kit / "bwa-postalt.js"), "-r", f"{fraction:g}"]
    hla_files: list[Path] = []
    if hla_prefix is not None:
        stem = Path(_cli.text(hla_prefix, "hla_prefix")).expanduser().resolve()
        if not stem.parent.is_dir():
            raise ConfigurationError(f"The folder for hla_prefix does not exist: {stem.parent}")
        args += ["-p", str(stem)]
        hla_files = [
            Path(f"{stem}.{gene}.fq")
            for gene in ("HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")
        ]
        clash = [str(path) for path in hla_files if path.exists()]
        if clash:
            raise ConfigurationError(
                f"Output files already exist and will not be overwritten: {clash}"
            )
    args += [str(alt), str(alignments)]
    out = _cli.new_file(output_path, "output_path")
    result = _p13.run_to_file(
        str(kit / "k8"),
        args,
        output=out,
        explicit=str(kit / "k8"),
        hint=_KIT_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["version"] = _kit_version(kit)
    result["operation"] = "bwa-postalt.js"
    result["backend"] = "bwakit"
    result["bwakit_directory"] = str(kit)
    result["hla_files"] = [str(path) for path in hla_files if path.is_file()]
    result["alignment_count"] = _p13.count_lines(out, skip_prefixes=("@",)) if out.is_file() else 0
    return result


def bwakit_type_hla(
    prefix: str,
    gene: str,
    *,
    from_contigs: bool = False,
    max_genotypes: int = 50,
    bwakit_directory: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """HLA typing of one gene from the reads (or contigs) that match it (bwa.kit ``typeHLA.sh``).

    ``prefix`` and ``gene`` name the input ``<prefix>.<gene>.fq`` (for example ``sample.hla.HLA-A.fq``, the file that
    ``run-bwamem -H`` or :func:`bwakit_postalt` with ``hla_prefix`` writes). By default the reads are assembled
    de novo (fermi2) and the exons are mapped to the contigs; with ``from_contigs=True`` (``-A``) the input is already
    a set of assembled contigs and the assembly step is skipped. The genotypes are written to ``<prefix>.<gene>.gt``
    (which must not exist yet) and the best ``max_genotypes`` are returned in ``genotypes`` as allele 1, allele 2,
    mismatches on the primary exons, mismatches on the other exons and number of exons used; the first one is the
    typing call. Temporary files are removed by the script, except ``<prefix>.<gene>.mag.gz`` of the de novo assembly.
    """
    stem = Path(_cli.text(prefix, "prefix")).expanduser().resolve()
    name = _cli.text(gene, "gene")
    if not _GENE.match(name):
        raise ConfigurationError("gene must be a gene name such as HLA-A.")
    limit = _cli.integer(max_genotypes, "max_genotypes", minimum=1)
    reads = Path(f"{stem}.{name}.fq")
    if not reads.is_file():
        raise ConfigurationError(f"Input file does not exist: {reads}")
    genotypes_file = Path(f"{stem}.{name}.gt")
    if genotypes_file.exists():
        raise ConfigurationError(
            f"Output already exists and will not be overwritten: {genotypes_file}"
        )
    kit = _bwakit_directory(bwakit_directory, ("typeHLA.sh", "k8", "resource-human-HLA"))
    args = ["-A"] if _switch(from_contigs, "-A", "from_contigs") else []
    args += [str(stem), name]
    result = _p13.run(
        str(kit / "typeHLA.sh"),
        args,
        cwd=stem.parent,
        expected=[genotypes_file],
        explicit=str(kit / "typeHLA.sh"),
        hint=_KIT_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["version"] = _kit_version(kit)
    result["operation"] = "typeHLA.sh"
    result["backend"] = "bwakit"
    result["bwakit_directory"] = str(kit)
    rows: list[dict[str, Any]] = []
    total = 0
    if genotypes_file.is_file():
        for line in genotypes_file.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split("\t")
            if fields[0] == "GT" and len(fields) >= 6:
                total += 1
                if total <= limit:
                    rows.append(
                        {
                            "allele1": fields[1],
                            "allele2": fields[2],
                            "mismatches_primary_exons": int(fields[3]),
                            "mismatches_other_exons": int(fields[4]),
                            "exons": int(fields[5]),
                        }
                    )
    result["genotype_count"] = total
    result["genotypes"] = rows
    result["typed"] = total > 0
    return result


def bwakit_generate_reference(
    analysis_set: str,
    output_directory: str,
    *,
    bwakit_directory: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and assemble a human reference analysis set with bwa.kit (``run-gen-ref``).

    ``analysis_set`` is ``"hs38"`` (GRCh38 primary assembly and EBV), ``"hs38a"`` (plus ALT contigs), ``"hs38DH"`` (plus
    decoys and HLA genes; recommended for GRCh38), ``"hs37"`` or ``"hs37d5"`` (GRCh37 with decoys). The script
    downloads about 1-3 GB from NCBI (network access and ``wget`` are required) and writes ``<analysis_set>.fa`` (and,
    for the ALT-aware sets, ``<analysis_set>.fa.alt``) into the new ``output_directory``; the reference still has to be
    indexed with ``bwa index``.
    """
    name = _cli.choice(analysis_set, _REFERENCE_SETS, "analysis_set")
    kit = _bwakit_directory(bwakit_directory, ("run-gen-ref",))
    if shutil.which("wget") is None:
        raise BackendUnavailableError(
            "wget was not found on PATH.", hint="Install wget; run-gen-ref downloads with it."
        )
    directory = _cli.new_directory(output_directory, "output_directory")
    expected = [directory / f"{name}.fa"]
    result = _p13.run(
        str(kit / "run-gen-ref"),
        [name],
        cwd=directory,
        expected=expected,
        explicit=str(kit / "run-gen-ref"),
        hint=_KIT_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["version"] = _kit_version(kit)
    result["operation"] = "run-gen-ref"
    result["backend"] = "bwakit"
    result["bwakit_directory"] = str(kit)
    alt = directory / f"{name}.fa.alt"
    result["alt_file"] = str(alt) if alt.is_file() else None
    return result


# ---------------------------------------------------------------------------
# qualfa2fq.pl
# ---------------------------------------------------------------------------


def bwa_qual_fasta_to_fastq(
    fasta_path: str,
    quality_path: str,
    output_path: str,
    *,
    script_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine a FASTA file and a quality-value file into FASTQ (BWA's ``qualfa2fq.pl``).

    ``quality_path`` is a FASTA-like file whose records hold whitespace-separated Phred scores (as written by the old
    454 and SOLiD pipelines); each score is converted to the character with code score + 33. Plain or ``.gz`` files
    are accepted. The script ships with BWA (``qualfa2fq.pl``; the conda package installs it on PATH) and needs
    perl; pass ``script_path`` for another location.
    """
    fasta = _cli.existing_file(fasta_path, "fasta_path")
    quality = _cli.existing_file(quality_path, "quality_path")
    _cli.program("qualfa2fq.pl", explicit=script_path, hint=_HINT)
    out = _cli.new_file(output_path, "output_path")
    result = _p13.run_to_file(
        "qualfa2fq.pl",
        [str(fasta), str(quality)],
        output=out,
        explicit=script_path,
        hint=_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = "qualfa2fq.pl"
    result["record_count"] = _fastq_records(out)
    return result


def _fastq_records(path: Path) -> int:
    """Records of a FASTQ whose quality line may be wrapped (qualfa2fq.pl wraps at 60): count the ``+`` lines."""
    if not path.is_file():
        return 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.rstrip("\n") == "+")


__all__ = [
    "bwa_aln",
    "bwa_fastmap",
    "bwa_index_stage",
    "bwa_maxk",
    "bwa_pemerge",
    "bwa_qual_fasta_to_fastq",
    "bwa_sai_to_sam",
    "bwakit_generate_reference",
    "bwakit_postalt",
    "bwakit_type_hla",
]

_register_public_verbose(globals(), tuple(__all__))
