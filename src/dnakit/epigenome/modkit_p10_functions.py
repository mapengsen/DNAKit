"""Oxford Nanopore modkit command-line facades (modified-base tables, bedMethyl, modBAM, motifs and DMRs).

The user installs modkit (a release binary from https://github.com/nanoporetech/modkit/releases, or ``conda install -c
bioconda modkit``) and either puts it on PATH or passes ``modkit_path``; 0.6.4 is the tested version. Each public function
runs one modkit subcommand with an explicit typed signature mapped to a fixed set of long flags: no shell, no generic
"run any command", positional arguments always before options (modkit's multi-value options such as ``--modified-bases``
would otherwise swallow them). Outputs are new explicit paths or fresh directories, so an existing result is never
overwritten. Results are JSON-serialisable dictionaries with status, backend, version, command, return code, existing
outputs, a bounded log tail and, for tables, a row count and a short preview. An option that is left at ``None`` or
``False`` is not passed, so modkit's own defaults apply (the documented defaults of modkit 0.6.4 are named in each
docstring where they matter). Not every option of a subcommand is exposed; the ones that change the inputs, the
algorithm and the outputs are.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p10_cli as _p10
from dnakit.exceptions import ConfigurationError

_PROGRAM = "modkit"
_TESTED_VERSION = "0.6.4"
_CITATION = "https://github.com/nanoporetech/modkit"
_HINT = (
    "Install modkit 0.6.4 or later (a release binary from https://github.com/nanoporetech/modkit/releases, "
    "or `conda install -c bioconda modkit`), then put it on PATH or pass modkit_path."
)
_MOTIF = re.compile(r"[A-Za-z]+")
_MOD_CODE = re.compile(r"[A-Za-z0-9:]+")
_THRESHOLD = re.compile(r"[A-Za-z0-9:]+:(?:\d+\.?\d*|\.\d+)")
_EDGE = re.compile(r"\d+(?:,\d+)?")
_NAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*")


# ---------------------------------------------------------------------------
# Argument helpers
# ---------------------------------------------------------------------------


class _Args:
    """An argument vector: positionals first, then flags in the order they are added."""

    def __init__(self, *positional: str) -> None:
        self.items: list[str] = list(positional)

    def switch(self, flag: str, on: bool, name: str) -> None:
        if not isinstance(on, bool):
            raise ConfigurationError(f"{name} must be True or False.")
        if on:
            self.items.append(flag)

    def value(self, flag: str, value: object | None) -> None:
        if value is not None:
            self.items += [flag, str(value)]

    def many(self, flag: str, values: Sequence[str] | None) -> None:
        if values:
            self.items += [flag, *values]


def _file(path: str | Path, name: str) -> str:
    return str(_cli.existing_file(path, name))


def _optional_file(path: str | Path | None, name: str) -> str | None:
    return None if path is None else _file(path, name)


def _optional_int(value: int | None, name: str, *, minimum: int = 0) -> int | None:
    return None if value is None else _cli.integer(value, name, minimum=minimum)


def _fraction(value: float | None, name: str) -> float | None:
    return None if value is None else _cli.number(value, name, minimum=0.0, maximum=1.0)


def _pattern(value: str | None, name: str, expression: re.Pattern[str]) -> str | None:
    if value is None:
        return None
    text = _cli.text(value, name)
    if expression.fullmatch(text) is None:
        raise ConfigurationError(f"{name} has an invalid format: {text!r}")
    return text


def _codes(values: Sequence[str] | None, name: str) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    out = []
    for value in values:
        text = _cli.text(value, name)
        if _MOD_CODE.fullmatch(text) is None:
            raise ConfigurationError(
                f"{name} entries may contain only letters, digits and ':': {text!r}"
            )
        out.append(text)
    return out


def _thresholds(values: Sequence[str] | None, name: str) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list such as ['m:0.8'].")
    out = []
    for value in values:
        text = _cli.text(value, name)
        if _THRESHOLD.fullmatch(text) is None:
            raise ConfigurationError(f"{name} entries must look like 'm:0.8': {text!r}")
        out.append(text)
    return out


def _region(value: str | None) -> str | None:
    return None if value is None else _cli.text(value, "region")


def _edge_filter(value: str | None) -> str | None:
    return _pattern(value, "edge_filter", _EDGE)


def _motif_pair(motif: str | None, offset: int, flag: str, args: _Args) -> None:
    if motif is None:
        return
    seq = _cli.text(motif, "motif")
    if _MOTIF.fullmatch(seq) is None:
        raise ConfigurationError(f"motif must be letters only (IUPAC codes allowed): {seq!r}")
    position = _cli.integer(offset, "motif_offset", minimum=0)
    if position >= len(seq):
        raise ConfigurationError("motif_offset must be smaller than the motif length.")
    args.items += [flag, seq, str(position)]


def _triples(values: Sequence[Sequence[str]] | None, name: str) -> list[tuple[str, str, str]]:
    """Known motifs as (sequence, offset, modification code) triples."""
    out: list[tuple[str, str, str]] = []
    for item in values or ():
        if isinstance(item, str) or len(item) != 3:
            raise ConfigurationError(f"each entry of {name} must be [sequence, offset, mod_code].")
        seq = _cli.text(item[0], name)
        if _MOTIF.fullmatch(seq) is None:
            raise ConfigurationError(f"{name} sequence must be letters only: {seq!r}")
        offset = str(
            _cli.integer(int(item[1]) if not isinstance(item[1], int) else item[1], name, minimum=0)
        )
        code = _cli.text(item[2], name)
        if _MOD_CODE.fullmatch(code) is None:
            raise ConfigurationError(f"{name} mod code is invalid: {code!r}")
        out.append((seq, offset, code))
    return out


def _context(value: Sequence[int] | None, name: str) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, str) or len(value) != 2:
        raise ConfigurationError(f"{name} must be two integers: (before, after).")
    return [
        str(_cli.integer(value[0], name, minimum=0)),
        str(_cli.integer(value[1], name, minimum=0)),
    ]


def _name(value: str, name: str) -> str:
    text = _cli.text(value, name)
    if _NAME.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} may contain only letters, digits, '_', '.' and '-': {text!r}"
        )
    return text


def _filtering(
    args: _Args,
    percentile: float | None,
    threshold: float | None,
    mod_thresholds: Sequence[str] | None,
    no_filtering: bool,
    *,
    percentile_flag: str = "--filter-percentile",
) -> None:
    """The pass/fail threshold options shared by the BAM based subcommands."""
    args.switch("--no-filtering", no_filtering, "no_filtering")
    if percentile is not None:
        args.value(
            percentile_flag,
            _cli.number(percentile, "filter_percentile", minimum=0.0, maximum=100.0),
        )
    if threshold is not None:
        args.value(
            "--filter-threshold",
            _cli.number(threshold, "filter_threshold", minimum=0.0, maximum=1.0),
        )
    for entry in _thresholds(mod_thresholds, "mod_thresholds") or ():
        args.value("--mod-threshold", entry)


def _sampling(args: _Args, num_reads: int | None, fraction: float | None, seed: int | None) -> None:
    args.value("--num-reads", _optional_int(num_reads, "num_reads", minimum=1))
    args.value("--sampling-frac", _fraction(fraction, "sampling_frac"))
    args.value("--seed", _optional_int(seed, "seed"))


def _threads(args: _Args, threads: int | None, flag: str = "--threads") -> None:
    args.value(flag, _optional_int(threads, "threads", minimum=1))


def _modkit(
    subcommand: Sequence[str],
    args: _Args,
    *,
    cwd: Path,
    expected: Sequence[Path],
    modkit_path: str | None,
    progress: bool = True,
) -> dict[str, Any]:
    arguments = [*subcommand, *args.items]
    if progress:
        arguments.append("--suppress-progress")
    return _cli.run(
        _PROGRAM, arguments, cwd=cwd, expected=list(expected), explicit=modkit_path, hint=_HINT
    )


def _finish(
    operation: str,
    rendered: dict[str, Any],
    parse: Callable[[], dict[str, Any]] | None = None,
    **extra: Any,
) -> dict[str, Any]:
    """Add the operation name, provenance and (only after a completed run) the parsed output and extra keys."""
    result = dict(rendered)
    result["operation"] = operation
    result["provenance"] = _p10.provenance(
        "modkit", _TESTED_VERSION, rendered.get("version"), _CITATION
    )
    if result["status"] == "completed":
        result.update(extra)
        if parse is not None:
            result.update(parse())
    return result


def _check_executable(modkit_path: str | None) -> None:
    """Resolve modkit before anything is created on disk, so a missing program leaves no output behind."""
    _cli.program(_PROGRAM, explicit=modkit_path, hint=_HINT)


def _new_dir(path: str | Path, name: str, modkit_path: str | None) -> Path:
    """A fresh output directory, created only after modkit has been found."""
    _check_executable(modkit_path)
    return _cli.new_directory(path, name)


def _new_output(path: str | Path, name: str = "output_path") -> Path:
    return _cli.new_file(path, name)


def _bedmethyl_preview(path: Path) -> dict[str, Any]:
    return {"record_count": _p10.count_records(path, comment="#")}


# ---------------------------------------------------------------------------
# pileup and pileup-hemi
# ---------------------------------------------------------------------------


def modkit_pileup(
    bam_path: str,
    output_path: str,
    *,
    modified_bases: Sequence[str] | None = None,
    reference: str | None = None,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    combine_mods: bool = False,
    combine_strands: bool = False,
    mask: bool = False,
    region: str | None = None,
    include_bed: str | None = None,
    edge_filter: str | None = None,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    no_filtering: bool = False,
    max_depth: int | None = None,
    num_reads: int | None = None,
    sampling_frac: float | None = None,
    seed: int | None = None,
    duplex: bool = False,
    header: bool = False,
    bgzf: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tabulate base-modification calls per genomic position into a bedMethyl file (``modkit pileup``).

    ``bam_path`` is a sorted, indexed modBAM. With ``cpg=True`` or a ``motif`` (sequence plus 0-based
    ``motif_offset``) only those sites are reported and ``reference`` (a FASTA with .fai) is required.
    ``modified_bases`` restricts the output to modifications such as ``['5mC']`` or ``['C:m']``;
    ``combine_mods`` sums all modification codes of a base; ``combine_strands`` sums both strands of a
    motif onto the positive strand. Unless ``no_filtering`` is set, modkit estimates a pass threshold at
    ``filter_percentile`` (default 10) from a sample of the reads. ``bgzf=True`` writes bgzip output, which you
    index with ``tabix -p bed`` before using it as input of ``modkit_stats``, ``modkit_localize``,
    ``modkit_bedmethyl_merge`` or ``modkit_dmr``. Result key: ``record_count``.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_path)
    ref = _optional_file(reference, "reference")
    bed = _optional_file(include_bed, "include_bed")
    if (cpg or motif is not None) and ref is None:
        raise ConfigurationError("reference is required with cpg=True or a motif.")
    args = _Args(bam, str(out))
    args.switch("--duplex", duplex, "duplex")
    args.switch("--bgzf", bgzf, "bgzf")
    args.switch("--header", header, "header")
    args.switch("--cpg", cpg, "cpg")
    _motif_pair(motif, motif_offset, "--motif", args)
    args.value("--reference", ref)
    args.switch("--mask", mask, "mask")
    args.switch("--combine-mods", combine_mods, "combine_mods")
    args.switch("--combine-strands", combine_strands, "combine_strands")
    args.value("--region", _region(region))
    args.value("--include-bed", bed)
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.value("--max-depth", _optional_int(max_depth, "max_depth", minimum=1))
    _filtering(args, filter_percentile, filter_threshold, mod_thresholds, no_filtering)
    _sampling(args, num_reads, sampling_frac, seed)
    _threads(args, threads)
    codes = _codes(modified_bases, "modified_bases")
    args.many("--modified-bases", codes)
    _check_executable(modkit_path)
    rendered = _modkit(["pileup"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path)
    return _finish("modkit_pileup", rendered, parse=lambda: _bedmethyl_preview(out))


def modkit_pileup_hemi(
    bam_path: str,
    reference: str,
    output_path: str,
    *,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    combine_mods: bool = False,
    mask: bool = False,
    ignore: str | None = None,
    region: str | None = None,
    include_bed: str | None = None,
    edge_filter: str | None = None,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    no_filtering: bool = False,
    max_depth: int | None = None,
    num_reads: int | None = None,
    sampling_frac: float | None = None,
    seed: int | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tabulate double-stranded (hemi-)methylation patterns across motif positions (``modkit pileup-hemi``).

    One of ``cpg=True`` (shorthand for motif ``CG`` offset 0) or ``motif`` with ``motif_offset`` is required;
    ``reference`` is a FASTA with .fai. ``ignore`` drops one modification code. The output is a bedMethyl
    table whose rows describe the pattern on both strands of each motif site. Result key: ``record_count``.
    """
    bam = _file(bam_path, "bam_path")
    ref = _file(reference, "reference")
    out = _new_output(output_path)
    if cpg == (motif is not None):
        raise ConfigurationError("Give exactly one of cpg=True or motif.")
    args = _Args(bam)
    args.value("--ref", ref)
    args.value("--out-bed", out)
    args.switch("--cpg", cpg, "cpg")
    _motif_pair(motif, motif_offset, "--motif", args)
    args.value("--ignore", _pattern(ignore, "ignore", _MOD_CODE))
    args.switch("--mask", mask, "mask")
    args.switch("--combine-mods", combine_mods, "combine_mods")
    args.value("--region", _region(region))
    args.value("--include-bed", _optional_file(include_bed, "include_bed"))
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.value("--max-depth", _optional_int(max_depth, "max_depth", minimum=1))
    _filtering(args, filter_percentile, filter_threshold, mod_thresholds, no_filtering)
    _sampling(args, num_reads, sampling_frac, seed)
    _threads(args, threads)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["pileup-hemi"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish("modkit_pileup_hemi", rendered, parse=lambda: _bedmethyl_preview(out))


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------


def modkit_validate(
    bam_and_bed: Sequence[Sequence[str]],
    output_path: str,
    *,
    canonical_base: str | None = None,
    ignore: str | None = None,
    edge_filter: str | None = None,
    invert_edge_filter: bool = False,
    min_alignment_identity: float | None = None,
    min_alignment_length: int | None = None,
    filter_quantile: float | None = None,
    filter_threshold: float | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate modified-base calls against ground-truth positions (``modkit validate``).

    ``bam_and_bed`` is a list of ``[modBAM, truth BED]`` pairs, one per sample; the name column of each BED
    holds the modification code (or ``-`` for a canonical-base position). ``canonical_base`` (A, C, G or T)
    is needed when the truth has only canonical sites or ChEBI codes. ``output_path`` receives modkit's
    machine-readable summary table. Result keys: ``row_count``, ``rows`` (preview) and ``columns``.
    """
    if isinstance(bam_and_bed, str) or not bam_and_bed:
        raise ConfigurationError("bam_and_bed must be a non-empty list of [bam, bed] pairs.")
    out = _new_output(output_path)
    args = _Args()
    for pair in bam_and_bed:
        if isinstance(pair, str) or len(pair) != 2:
            raise ConfigurationError("each entry of bam_and_bed must be [bam_path, bed_path].")
        args.items += [
            "--bam-and-bed",
            _file(pair[0], "bam_and_bed"),
            _file(pair[1], "bam_and_bed"),
        ]
    if canonical_base is not None:
        args.value(
            "--canonical-base", _cli.choice(canonical_base, ("A", "C", "G", "T"), "canonical_base")
        )
    args.value("--ignore", _pattern(ignore, "ignore", _MOD_CODE))
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.switch("--invert-edge-filter", invert_edge_filter, "invert_edge_filter")
    if min_alignment_identity is not None:
        args.value(
            "--min-identity",
            _cli.number(min_alignment_identity, "min_alignment_identity", minimum=0.0, maximum=1.0),
        )
    args.value(
        "--min-length", _optional_int(min_alignment_length, "min_alignment_length", minimum=1)
    )
    args.value("--filter-quantile", _fraction(filter_quantile, "filter_quantile"))
    args.value("--filter-threshold", _fraction(filter_threshold, "filter_threshold"))
    _threads(args, threads)
    args.value("--out-filepath", out)
    _check_executable(modkit_path)
    rendered = _modkit(["validate"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path)
    return _finish(
        "modkit_validate",
        rendered,
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


# ---------------------------------------------------------------------------
# motif subcommands
# ---------------------------------------------------------------------------


def modkit_motif_bed(
    fasta_path: str,
    motif: str,
    offset: int,
    output_path: str,
    *,
    mask: bool = False,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BED file with every location of a sequence motif in a FASTA (``modkit motif bed``).

    ``motif`` may use IUPAC codes; ``offset`` is the 0-based position inside the motif that the BED interval
    marks (``CG`` with offset 0 gives the C of every CpG on the forward strand). ``mask=True`` respects
    soft-masking. Result key: ``record_count``.
    """
    fasta = _file(fasta_path, "fasta_path")
    out = _new_output(output_path)
    seq = _cli.text(motif, "motif")
    if _MOTIF.fullmatch(seq) is None:
        raise ConfigurationError(f"motif must be letters only (IUPAC codes allowed): {seq!r}")
    position = _cli.integer(offset, "offset", minimum=0)
    if position >= len(seq):
        raise ConfigurationError("offset must be smaller than the motif length.")
    args = _Args("bed", fasta, seq, str(position))
    args.switch("--mask", mask, "mask")
    executable = _cli.program(_PROGRAM, explicit=modkit_path, hint=_HINT)
    rendered = _p10.stdout_to_file(
        _PROGRAM,
        ["motif", *args.items],
        output=out,
        cwd=out.parent,
        explicit=modkit_path,
        hint=_HINT,
        program_path=executable,
    )
    return _finish(
        "modkit_motif_bed",
        rendered,
        parse=lambda: {"record_count": _p10.count_records(out, comment=None)},
    )


def _motif_common(
    args: _Args,
    bedmethyl: str,
    reference: str,
    contig: str | None,
    known_motifs: Sequence[Sequence[str]] | None,
    known_motifs_table: str | None,
    min_coverage: int | None,
    context_size: Sequence[int] | None,
    low_threshold: float | None,
    high_threshold: float | None,
    threads: int | None,
) -> None:
    args.value("--in-bedmethyl", _file(bedmethyl, "bedmethyl_path"))
    args.value("--ref", _file(reference, "reference"))
    args.value("--contig", None if contig is None else _cli.text(contig, "contig"))
    for seq, offset, code in _triples(known_motifs, "known_motifs"):
        args.items += ["--known-motif", seq, offset, code]
    args.value("--known-motifs-table", _optional_file(known_motifs_table, "known_motifs_table"))
    args.value("--min-coverage", _optional_int(min_coverage, "min_coverage", minimum=1))
    context = _context(context_size, "context_size")
    if context:
        args.items += ["--context-size", *context]
    args.value("--low-thresh", _fraction(low_threshold, "low_threshold"))
    args.value("--high-thresh", _fraction(high_threshold, "high_threshold"))
    _threads(args, threads)


def modkit_motif_search(
    bedmethyl_path: str,
    reference: str,
    output_table: str,
    *,
    contig: str | None = None,
    mod_codes: Sequence[str] | None = None,
    known_motifs: Sequence[Sequence[str]] | None = None,
    known_motifs_table: str | None = None,
    context_size: Sequence[int] | None = None,
    init_context_size: Sequence[int] | None = None,
    min_coverage: int | None = None,
    min_sites: int | None = None,
    min_log_odds: float | None = None,
    min_frac_mod: float | None = None,
    low_threshold: float | None = None,
    high_threshold: float | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a genome for modification-enriched sequence motifs from a bedMethyl table (``modkit motif search``).

    ``bedmethyl_path`` is a bgzip bedMethyl with a tabix index and ``reference`` a FASTA with .fai. The
    machine-readable table of enriched motifs is written to ``output_table``. ``known_motifs`` are
    ``[sequence, offset, mod_code]`` triples whose statistics are reported next to the discovered ones;
    ``mod_codes`` limits the search to those modification codes. Result keys: ``row_count``, ``columns``, ``rows``.
    """
    out = _new_output(output_table, "output_table")
    args = _Args()
    _motif_common(
        args,
        bedmethyl_path,
        reference,
        contig,
        known_motifs,
        known_motifs_table,
        min_coverage,
        context_size,
        low_threshold,
        high_threshold,
        threads,
    )
    init = _context(init_context_size, "init_context_size")
    if init:
        args.items += ["--init-context-size", *init]
    args.value("--min-sites", _optional_int(min_sites, "min_sites", minimum=1))
    if min_log_odds is not None:
        args.value("--min-log-odds", _cli.number(min_log_odds, "min_log_odds"))
    args.value("--min-frac-mod", _fraction(min_frac_mod, "min_frac_mod"))
    args.value("--out-table", out)
    for code in _codes(mod_codes, "mod_codes") or ():
        args.value("--mod-code", code)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["motif", "search"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish(
        "modkit_motif_search",
        rendered,
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


def modkit_motif_evaluate(
    bedmethyl_path: str,
    reference: str,
    output_table: str,
    *,
    known_motifs: Sequence[Sequence[str]] | None = None,
    known_motifs_table: str | None = None,
    contig: str | None = None,
    context_size: Sequence[int] | None = None,
    min_coverage: int | None = None,
    low_threshold: float | None = None,
    high_threshold: float | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute enrichment statistics of given motifs from a bedMethyl table (``modkit motif evaluate``).

    The motifs are ``known_motifs`` (``[sequence, offset, mod_code]`` triples) and/or a ``known_motifs_table``
    TSV (``mod_code``, ``motif_seq``, ``offset``). At least one of them is required. Result keys: ``row_count``,
    ``columns``, ``rows``.
    """
    if not known_motifs and known_motifs_table is None:
        raise ConfigurationError("Give known_motifs or known_motifs_table.")
    out = _new_output(output_table, "output_table")
    args = _Args()
    _motif_common(
        args,
        bedmethyl_path,
        reference,
        contig,
        known_motifs,
        known_motifs_table,
        min_coverage,
        context_size,
        low_threshold,
        high_threshold,
        threads,
    )
    args.value("--out", out)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["motif", "evaluate"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish(
        "modkit_motif_evaluate",
        rendered,
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


def modkit_motif_refine(
    bedmethyl_path: str,
    reference: str,
    output_table: str,
    *,
    known_motifs: Sequence[Sequence[str]] | None = None,
    known_motifs_table: str | None = None,
    contig: str | None = None,
    context_size: Sequence[int] | None = None,
    min_coverage: int | None = None,
    min_sites: int | None = None,
    min_log_odds: float | None = None,
    min_frac_mod: float | None = None,
    min_refine_frac_mod: float | None = None,
    min_refine_sites: int | None = None,
    low_threshold: float | None = None,
    high_threshold: float | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Refine a list of known motifs with a bedMethyl table (``modkit motif refine``).

    ``known_motifs`` (``[sequence, offset, mod_code]``) and/or ``known_motifs_table`` give the motifs to refine
    (at least one is required). ``min_refine_frac_mod`` and ``min_refine_sites`` bound the refinement.
    Result keys: ``row_count``, ``columns``, ``rows``.
    """
    if not known_motifs and known_motifs_table is None:
        raise ConfigurationError("Give known_motifs or known_motifs_table.")
    out = _new_output(output_table, "output_table")
    args = _Args()
    _motif_common(
        args,
        bedmethyl_path,
        reference,
        contig,
        known_motifs,
        known_motifs_table,
        min_coverage,
        context_size,
        low_threshold,
        high_threshold,
        threads,
    )
    args.value("--min-sites", _optional_int(min_sites, "min_sites", minimum=1))
    if min_log_odds is not None:
        args.value("--min-log-odds", _cli.number(min_log_odds, "min_log_odds"))
    args.value("--min-frac-mod", _fraction(min_frac_mod, "min_frac_mod"))
    args.value("--min-refine-frac-mod", _fraction(min_refine_frac_mod, "min_refine_frac_mod"))
    args.value("--min-refine-sites", _optional_int(min_refine_sites, "min_refine_sites", minimum=1))
    args.value("--out", out)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["motif", "refine"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish(
        "modkit_motif_refine",
        rendered,
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


# ---------------------------------------------------------------------------
# entropy, localize, stats
# ---------------------------------------------------------------------------


def modkit_entropy(
    bam_paths: Sequence[str],
    reference: str,
    output_path: str,
    *,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    base: str | None = None,
    regions_bed: str | None = None,
    num_positions: int | None = None,
    window_size: int | None = None,
    min_coverage: int | None = None,
    max_filtered_positions: int | None = None,
    combine_strands: bool = False,
    mask: bool = False,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    no_filtering: bool = False,
    num_reads: int | None = None,
    header: bool = False,
    drop_zeros: bool = False,
    prefix: str | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Methylation entropy over genomic windows of ``num_positions`` modified sites (``modkit entropy``).

    One or more modBAMs (``bam_paths``, pooled when several) and a FASTA ``reference`` are read. Choose the
    sites with ``cpg=True``, a ``motif`` (+ ``motif_offset``) or a ``base``. Without ``regions_bed`` the
    windows go to the new BED file ``output_path``; with ``regions_bed`` ``output_path`` is a new directory
    that receives one BED per region (files start with ``prefix``) plus modkit's regions summary. modkit's
    defaults are 4 positions in a window of at most 50 bp. Result keys: ``record_count`` (BED mode) or
    ``files`` (directory mode).
    """
    if isinstance(bam_paths, str) or not bam_paths:
        raise ConfigurationError("bam_paths must be a non-empty list of modBAM files.")
    bams = [_file(path, "bam_paths") for path in bam_paths]
    ref = _file(reference, "reference")
    regions = _optional_file(regions_bed, "regions_bed")
    if regions is None:
        out = _new_output(output_path)
        directory_mode = False
    else:
        out = _new_dir(output_path, "output_path", modkit_path)
        directory_mode = True
    args = _Args()
    for bam in bams:
        args.value("--in-bam", bam)
    args.value("--ref", ref)
    args.switch("--cpg", cpg, "cpg")
    _motif_pair(motif, motif_offset, "--motif", args)
    args.value("--base", None if base is None else _cli.choice(base, ("A", "C", "G", "T"), "base"))
    args.value("--regions", regions)
    args.value("--num-positions", _optional_int(num_positions, "num_positions", minimum=1))
    args.value("--window-size", _optional_int(window_size, "window_size", minimum=1))
    args.value("--min-coverage", _optional_int(min_coverage, "min_coverage", minimum=1))
    args.value(
        "--max-filtered-positions",
        _optional_int(max_filtered_positions, "max_filtered_positions", minimum=0),
    )
    args.switch("--combine-strands", combine_strands, "combine_strands")
    args.switch("--mask", mask, "mask")
    _filtering(args, filter_percentile, filter_threshold, mod_thresholds, no_filtering)
    args.value("--num-reads", _optional_int(num_reads, "num_reads", minimum=1))
    args.switch("--header", header, "header")
    args.switch("--drop-zeros", drop_zeros, "drop_zeros")
    args.value("--prefix", None if prefix is None else _name(prefix, "prefix"))
    _threads(args, threads)
    args.value("--out-bed", out)
    _check_executable(modkit_path)
    if directory_mode:
        rendered = _modkit(["entropy"], args, cwd=out, expected=[], modkit_path=modkit_path)
        files = sorted(str(path) for path in out.iterdir())
        rendered["outputs"] = files
        if rendered["return_code"] == 0 and not files:
            rendered["status"] = "failed"
        return _finish("modkit_entropy", rendered, files=files)
    rendered = _modkit(["entropy"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path)
    return _finish(
        "modkit_entropy",
        rendered,
        parse=lambda: {"record_count": _p10.count_records(out, comment="#")},
    )


def modkit_localize(
    bedmethyl_path: str,
    regions_bed: str,
    genome_sizes: str,
    output_path: str,
    *,
    window: int | None = None,
    stranded: str | None = None,
    stranded_features: bool = False,
    min_coverage: int | None = None,
    chart_path: str | None = None,
    chart_name: str | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Aggregate modification levels in windows centred on genomic features (``modkit localize``; alias ``localise``).

    ``bedmethyl_path`` must be bgzip-compressed with a tabix index (``<file>.tbi``); ``regions_bed`` holds the
    features, ``genome_sizes`` a two-column chromosome sizes file. ``window`` is the number of bases on each
    side (modkit default 2000). ``chart_path`` optionally writes an HTML chart. Result keys: ``row_count``,
    ``columns``, ``rows`` (the position-by-position table), ``chart`` (path or None).
    """
    bed = _file(bedmethyl_path, "bedmethyl_path")
    out = _new_output(output_path)
    chart = None if chart_path is None else _new_output(chart_path, "chart_path")
    args = _Args(bed)
    args.value("--regions", _file(regions_bed, "regions_bed"))
    args.value("--genome-sizes", _file(genome_sizes, "genome_sizes"))
    args.value("--window", _optional_int(window, "window", minimum=1))
    if stranded is not None:
        args.value("--stranded", _cli.choice(stranded, ("true", "false", "both"), "stranded"))
    args.switch("--stranded-features", stranded_features, "stranded_features")
    args.value("--min-coverage", _optional_int(min_coverage, "min_coverage", minimum=1))
    args.value("--chart", chart)
    args.value("--name", None if chart_name is None else _cli.text(chart_name, "chart_name"))
    _threads(args, threads)
    args.value("--out-file", out)
    _check_executable(modkit_path)
    expected = [out] if chart is None else [out, chart]
    rendered = _modkit(
        ["localize"],
        args,
        cwd=out.parent,
        expected=expected,
        modkit_path=modkit_path,
        progress=False,
    )
    return _finish(
        "modkit_localize",
        rendered,
        chart=None if chart is None else str(chart),
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


def modkit_stats(
    bedmethyl_path: str,
    regions_bed: str,
    output_path: str,
    *,
    mod_codes: str | None = None,
    min_coverage: int | None = None,
    no_header: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Modification level (fraction modified) over each region of a BED file (``modkit stats``).

    ``bedmethyl_path`` must be bgzip-compressed with a tabix index; ``mod_codes`` is a comma-separated list
    of modification codes (default: all codes found); ``min_coverage`` keeps records with at least that valid
    coverage (modkit default 1). Result keys: ``row_count``, ``columns`` (None with ``no_header``), ``rows``.
    """
    bed = _file(bedmethyl_path, "bedmethyl_path")
    out = _new_output(output_path)
    args = _Args(bed)
    args.value("--regions", _file(regions_bed, "regions_bed"))
    args.value("--mod-codes", _pattern(mod_codes, "mod_codes", re.compile(r"[A-Za-z0-9:,]+")))
    args.value("--min-coverage", _optional_int(min_coverage, "min_coverage", minimum=1))
    args.switch("--no-header", no_header, "no_header")
    _threads(args, threads)
    args.value("--out-table", out)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["stats"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path, progress=False
    )
    return _finish(
        "modkit_stats",
        rendered,
        parse=lambda: _p10.table_preview(out, header=not no_header, comment=None),
    )


# ---------------------------------------------------------------------------
# bedmethyl utilities
# ---------------------------------------------------------------------------


def modkit_bedmethyl_merge(
    bedmethyl_paths: Sequence[str],
    genome_sizes: str,
    output_path: str,
    *,
    header: bool = False,
    min_samples: int | None = None,
    min_sample_coverage: int | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Outer-join two or more bedMethyl files, summing the counts of records that overlap (``modkit bedmethyl merge``).

    Inputs are bgzip bedMethyl files with tabix indexes; ``genome_sizes`` is a two-column chromosome sizes
    file. ``min_samples`` and ``min_sample_coverage`` drop records seen in too few samples or with too little
    coverage. Result key: ``record_count``.
    """
    if isinstance(bedmethyl_paths, str) or len(bedmethyl_paths) < 2:
        raise ConfigurationError("bedmethyl_paths must list at least two bedMethyl files.")
    beds = [_file(path, "bedmethyl_paths") for path in bedmethyl_paths]
    out = _new_output(output_path)
    args = _Args(*beds)
    args.value("--genome-sizes", _file(genome_sizes, "genome_sizes"))
    args.switch("--header", header, "header")
    args.value("--min-samples", _optional_int(min_samples, "min_samples", minimum=1))
    args.value(
        "--min-sample-coverage",
        _optional_int(min_sample_coverage, "min_sample_coverage", minimum=1),
    )
    _threads(args, threads)
    args.value("--out-bed", out)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["bedmethyl", "merge"],
        args,
        cwd=out.parent,
        expected=[out],
        modkit_path=modkit_path,
        progress=False,
    )
    return _finish("modkit_bedmethyl_merge", rendered, parse=lambda: _bedmethyl_preview(out))


def modkit_bedmethyl_to_bigwig(
    bedmethyl_path: str,
    output_path: str,
    mod_codes: Sequence[str],
    *,
    chrom_sizes: str | None = None,
    header_bam: str | None = None,
    negative_strand_values: bool = False,
    nzooms: int | None = None,
    uncompressed: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a bedMethyl table into a BigWig track of fraction modified (``modkit bedmethyl tobigwig``).

    ``mod_codes`` lists the modification codes to write (one track per code is written with the code appended to
    the file name when several are given). The chromosome sizes come from ``chrom_sizes`` (two-column file) or
    from the header of ``header_bam``; exactly one of them is required. ``negative_strand_values`` writes
    negative-strand records as negative values.
    """
    bed = _file(bedmethyl_path, "bedmethyl_path")
    out = _new_output(output_path)
    if (chrom_sizes is None) == (header_bam is None):
        raise ConfigurationError("Give exactly one of chrom_sizes or header_bam.")
    codes = _codes(mod_codes, "mod_codes")
    args = _Args(bed, str(out))
    args.value("--sizes", _optional_file(chrom_sizes, "chrom_sizes"))
    args.value("--header", _optional_file(header_bam, "header_bam"))
    for code in codes or ():
        args.value("--mod-codes", code)
    args.switch("--negative-strand-values", negative_strand_values, "negative_strand_values")
    args.value("--nzooms", _optional_int(nzooms, "nzooms", minimum=0))
    args.switch("--uncompressed", uncompressed, "uncompressed")
    args.value("--nthreads", _optional_int(threads, "threads", minimum=1))
    _check_executable(modkit_path)
    rendered = _modkit(
        ["bedmethyl", "tobigwig"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish("modkit_bedmethyl_to_bigwig", rendered)


def modkit_bedmethyl_map_to_genome(
    bedmethyl_path: str,
    output_path: str,
    gtf_path: str,
    transcript_id: str,
    *,
    transcript_version: int | None = None,
    ignore_version: bool = False,
    header: bool = False,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map a transcriptome-aligned bedMethyl file onto genome coordinates with a GTF (``modkit bedmethyl map-to-genome``).

    The records of one transcript (``transcript_id``, optionally with ``transcript_version``; ``ignore_version``
    matches without the version suffix) are lifted to genome coordinates using the exons of ``gtf_path``.
    Result key: ``record_count``.
    """
    bed = _file(bedmethyl_path, "bedmethyl_path")
    out = _new_output(output_path)
    args = _Args(bed, str(out))
    args.value("--transcript-id", _cli.text(transcript_id, "transcript_id"))
    args.value(
        "--transcript-version", _optional_int(transcript_version, "transcript_version", minimum=0)
    )
    args.switch("--ignore-version", ignore_version, "ignore_version")
    args.value("--gtf", _file(gtf_path, "gtf_path"))
    args.switch("--header", header, "header")
    _check_executable(modkit_path)
    rendered = _modkit(
        ["bedmethyl", "map-to-genome"],
        args,
        cwd=out.parent,
        expected=[out],
        modkit_path=modkit_path,
        progress=False,
    )
    return _finish(
        "modkit_bedmethyl_map_to_genome", rendered, parse=lambda: _bedmethyl_preview(out)
    )


# ---------------------------------------------------------------------------
# modbam utilities
# ---------------------------------------------------------------------------


def modkit_modbam_check_tags(
    bam_path: str,
    output_dir: str,
    *,
    permissive: bool = False,
    region: str | None = None,
    num_reads: int | None = None,
    head: int | None = None,
    mapped_only: bool = False,
    allow_non_primary: bool = False,
    ignore_index: bool = False,
    prefix: str | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate the MM and ML base-modification tags of a BAM (``modkit modbam check-tags``).

    modkit exits non-zero when invalid records are found unless ``permissive`` is set; the status is then
    ``failed`` and ``log_tail`` holds the reason. ``output_dir`` is a new directory that receives the failure
    tables modkit writes (files start with ``prefix``). Result key: ``files``.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_dir(output_dir, "output_dir", modkit_path)
    args = _Args(bam)
    args.switch("--permissive", permissive, "permissive")
    args.value("--out-dir", out)
    args.value("--prefix", None if prefix is None else _name(prefix, "prefix"))
    _threads(args, threads)
    args.switch("--ignore-index", ignore_index, "ignore_index")
    args.value("--num-reads", _optional_int(num_reads, "num_reads", minimum=1))
    args.switch("--allow-non-primary", allow_non_primary, "allow_non_primary")
    args.switch("--mapped-only", mapped_only, "mapped_only")
    args.value("--region", _region(region))
    args.value("--head", _optional_int(head, "head", minimum=1))
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "check-tags"], args, cwd=out, expected=[], modkit_path=modkit_path
    )
    files = sorted(str(path) for path in out.iterdir())
    rendered["outputs"] = files
    return _finish("modkit_modbam_check_tags", rendered, files=files)


def modkit_modbam_adjust_mods(
    bam_path: str,
    output_bam: str,
    *,
    ignore: str | None = None,
    convert_from: str | None = None,
    convert_to: str | None = None,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    discard_motifs: bool = False,
    edge_filter: str | None = None,
    invert_edge_filter: bool = False,
    filter_probs: bool = False,
    mapped_only: bool = False,
    num_reads: int | None = None,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    output_sam: bool = False,
    fast_forward: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Adjust the base-modification calls of a modBAM and write a new BAM (``modkit modbam adjust-mods``).

    ``ignore`` removes calls of one modification code; ``convert_from`` and ``convert_to`` (given together)
    rename one code to another (collapsing, for example ``h`` into ``m``); ``cpg`` or ``motif`` with
    ``discard_motifs`` keeps or discards calls by sequence context; ``edge_filter``/``filter_probs`` and the
    thresholds drop calls near read ends or below the pass threshold. ``output_sam=True`` writes SAM
    instead of BAM.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_bam, "output_bam")
    if (convert_from is None) != (convert_to is None):
        raise ConfigurationError("convert_from and convert_to must be given together.")
    args = _Args(bam, str(out))
    args.switch("--ff", fast_forward, "fast_forward")
    args.switch("--output-sam", output_sam, "output_sam")
    args.value("--ignore", _pattern(ignore, "ignore", _MOD_CODE))
    if convert_from is not None:
        args.items += [
            "--convert",
            _pattern(convert_from, "convert_from", _MOD_CODE) or "",
            _pattern(convert_to, "convert_to", _MOD_CODE) or "",
        ]
    _motif_pair(motif, motif_offset, "--motif", args)
    args.switch("--cpg", cpg, "cpg")
    args.switch("--discard-motifs", discard_motifs, "discard_motifs")
    _threads(args, threads)
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.switch("--invert-edge-filter", invert_edge_filter, "invert_edge_filter")
    args.switch("--filter-probs", filter_probs, "filter_probs")
    args.switch("--mapped-only", mapped_only, "mapped_only")
    args.value("--num-reads", _optional_int(num_reads, "num_reads", minimum=1))
    _filtering(args, filter_percentile, filter_threshold, mod_thresholds, False)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "adjust-mods"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish("modkit_modbam_adjust_mods", rendered)


def modkit_modbam_update_tags(
    bam_path: str,
    output_bam: str,
    *,
    mode: str | None = None,
    no_implicit_probs: bool = False,
    output_sam: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rename lower-case ``Mm``/``Ml`` tags to ``MM``/``ML`` and set the mode flag (``modkit modbam update-tags``).

    ``mode`` is ``explicit`` (the ``?`` flag) or ``implicit`` (the ``.`` flag); ``no_implicit_probs`` stops modkit
    from adding implicit canonical calls. ``output_sam=True`` writes SAM.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_bam, "output_bam")
    args = _Args(bam, str(out))
    args.value(
        "--mode", None if mode is None else _cli.choice(mode, ("explicit", "implicit"), "mode")
    )
    args.switch("--no-implicit-probs", no_implicit_probs, "no_implicit_probs")
    args.switch("--output-sam", output_sam, "output_sam")
    _threads(args, threads)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "update-tags"],
        args,
        cwd=out.parent,
        expected=[out],
        modkit_path=modkit_path,
        progress=False,
    )
    return _finish("modkit_modbam_update_tags", rendered)


def modkit_modbam_sample_probs(
    bam_path: str,
    output_dir: str,
    *,
    percentiles: str | None = None,
    quantiles: str | None = None,
    histogram: bool = False,
    reference: str | None = None,
    motif: str | None = None,
    motif_offset: int = 0,
    region: str | None = None,
    include_bed: str | None = None,
    edge_filter: str | None = None,
    use_non_primary: bool = False,
    num_reads: int | None = None,
    sampling_frac: float | None = None,
    no_sampling: bool = False,
    seed: int | None = None,
    ignore_index: bool = False,
    prefix: str | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the distribution of base-modification probabilities in a modBAM (``modkit modbam sample-probs``).

    Writes modkit's probability and threshold tables (and, with ``histogram=True``, probability histograms) into the
    new directory ``output_dir``. ``percentiles`` is a comma-separated list such as ``"10,50,90"`` (modkit
    default) and ``quantiles`` one such as ``"0.1,0.5,0.9"``. Reads are sampled (``num_reads`` default 10042)
    unless ``no_sampling``. Result key: ``files``.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_dir(output_dir, "output_dir", modkit_path)
    number_list = re.compile(r"\d+(?:\.\d+)?(?:,\d+(?:\.\d+)?)*")
    args = _Args(bam)
    args.value("--percentiles", _pattern(percentiles, "percentiles", number_list))
    args.value("--quantiles", _pattern(quantiles, "quantiles", number_list))
    args.value("--out-dir", out)
    args.value("--prefix", None if prefix is None else _name(prefix, "prefix"))
    args.switch("--hist", histogram, "histogram")
    args.value("--reference", _optional_file(reference, "reference"))
    _motif_pair(motif, motif_offset, "--motif", args)
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.value("--region", _region(region))
    args.value("--include-bed", _optional_file(include_bed, "include_bed"))
    args.switch("--use-non-primary", use_non_primary, "use_non_primary")
    _sampling(args, num_reads, sampling_frac, seed)
    args.switch("--no-sampling", no_sampling, "no_sampling")
    args.switch("--ignore-index", ignore_index, "ignore_index")
    _threads(args, threads)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "sample-probs"], args, cwd=out, expected=[], modkit_path=modkit_path
    )
    files = sorted(str(path) for path in out.iterdir())
    rendered["outputs"] = files
    if rendered["return_code"] == 0 and not files:
        rendered["status"] = "failed"
    return _finish("modkit_modbam_sample_probs", rendered, files=files)


def modkit_modbam_summary(
    bam_path: str,
    output_path: str,
    *,
    tsv: bool = False,
    combine_mods: bool = False,
    cpg: bool = False,
    reference: str | None = None,
    motif: str | None = None,
    motif_offset: int = 0,
    region: str | None = None,
    include_bed: str | None = None,
    edge_filter: str | None = None,
    mapped_only: bool = False,
    matched_only: bool = False,
    num_reads: int | None = None,
    sampling_frac: float | None = None,
    seed: int | None = None,
    filter_quantile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    ignore_index: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the modification tags present in a BAM and basic call statistics (``modkit modbam summary``).

    Writes modkit's default totals table (``#`` lines) and calls table to ``output_path``, or one tab-separated table with
    ``tsv=True``. Reads are sampled (``num_reads`` default 10042) unless ``num_reads`` is raised;
    ``cpg``/``motif`` need ``reference``. Result keys: ``row_count``, ``rows`` (preview of the call rows).
    """
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_path)
    ref = _optional_file(reference, "reference")
    if (cpg or motif is not None) and ref is None:
        raise ConfigurationError("reference is required with cpg=True or a motif.")
    args = _Args(bam)
    args.switch("--ignore-index", ignore_index, "ignore_index")
    args.value("--reference", ref)
    args.switch("--cpg", cpg, "cpg")
    args.switch("--tsv", tsv, "tsv")
    args.switch("--combine-mods", combine_mods, "combine_mods")
    _sampling(args, num_reads, sampling_frac, seed)
    if filter_quantile is not None:
        args.value(
            "--filter-quantile",
            _cli.number(filter_quantile, "filter_quantile", minimum=0.0, maximum=1.0),
        )
    if filter_threshold is not None:
        args.value(
            "--filter-threshold",
            _cli.number(filter_threshold, "filter_threshold", minimum=0.0, maximum=1.0),
        )
    for entry in _thresholds(mod_thresholds, "mod_thresholds") or ():
        args.value("--mod-threshold", entry)
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.value("--include-bed", _optional_file(include_bed, "include_bed"))
    args.switch("--mapped-only", mapped_only, "mapped_only")
    args.switch("--matched-only", matched_only, "matched_only")
    _motif_pair(motif, motif_offset, "--motif", args)
    args.value("--region", _region(region))
    _threads(args, threads)
    args.items.append("--suppress-progress")
    executable = _cli.program(_PROGRAM, explicit=modkit_path, hint=_HINT)
    rendered = _p10.stdout_to_file(
        _PROGRAM,
        ["modbam", "summary", *args.items],
        output=out,
        cwd=out.parent,
        explicit=modkit_path,
        hint=_HINT,
        program_path=executable,
    )
    return _finish(
        "modkit_modbam_summary",
        rendered,
        parse=lambda: _p10.table_preview(out, header=tsv, comment="#"),
    )


def modkit_modbam_call_mods(
    bam_path: str,
    output_bam: str,
    *,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    discard_motifs: bool = False,
    edge_filter: str | None = None,
    invert_edge_filter: bool = False,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    no_filtering: bool = False,
    num_reads: int | None = None,
    sampling_frac: float | None = None,
    seed: int | None = None,
    output_sam: bool = False,
    fast_forward: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Threshold the modification probabilities and write a modBAM of hard calls (``modkit modbam call-mods``).

    Each base is set to probability 100 percent when a modification is called and to 0 percent when it is called
    canonical, using the same pass-threshold estimation as ``modkit_pileup`` (``filter_percentile`` default 10).
    ``cpg``/``motif`` restrict (or, with ``discard_motifs``, exclude) sequence contexts. ``output_sam=True`` writes SAM.
    """
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_bam, "output_bam")
    args = _Args(bam, str(out))
    args.switch("--ff", fast_forward, "fast_forward")
    _threads(args, threads)
    _sampling(args, num_reads, sampling_frac, seed)
    _filtering(args, filter_percentile, filter_threshold, mod_thresholds, no_filtering)
    args.value("--edge-filter", _edge_filter(edge_filter))
    args.switch("--invert-edge-filter", invert_edge_filter, "invert_edge_filter")
    _motif_pair(motif, motif_offset, "--motif", args)
    args.switch("--cpg", cpg, "cpg")
    args.switch("--discard-motifs", discard_motifs, "discard_motifs")
    args.switch("--output-sam", output_sam, "output_sam")
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "call-mods"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish("modkit_modbam_call_mods", rendered)


def modkit_modbam_repair(
    donor_bam: str,
    acceptor_bam: str,
    output_bam: str,
    *,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy correct MM/ML tags from a donor modBAM onto the reads of an acceptor BAM (``modkit modbam repair``).

    Both BAMs must be sorted by read name and the donor must contain every read of the acceptor (extra donor
    reads are allowed). Typical use: restore tags lost by a trimming or alignment step.
    """
    donor = _file(donor_bam, "donor_bam")
    acceptor = _file(acceptor_bam, "acceptor_bam")
    out = _new_output(output_bam, "output_bam")
    args = _Args()
    args.value("--donor-bam", donor)
    args.value("--acceptor-bam", acceptor)
    args.value("--output-bam", out)
    _threads(args, threads)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["modbam", "repair"],
        args,
        cwd=out.parent,
        expected=[out],
        modkit_path=modkit_path,
        progress=False,
    )
    return _finish("modkit_modbam_repair", rendered)


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

_EXTRACT_TABLES = ("full", "calls", "read-stats")
_EXTRACT_SHARED = frozenset(
    {
        "reference",
        "region",
        "include_bed",
        "exclude_bed",
        "motif",
        "cpg",
        "annotate_motifs",
        "mask",
        "ignore",
        "edge_filter",
        "invert_edge_filter",
        "mapped_only",
        "allow_non_primary",
        "num_reads",
        "kmer_size",
        "bgzf",
        "no_headers",
        "ignore_implicit",
    }
)
_EXTRACT_CALLS_ONLY = frozenset({"filter_percentile", "pass_only", "no_filtering", "sampling_frac"})


def modkit_extract(
    bam_path: str,
    output_path: str,
    *,
    table: str = "full",
    reference: str | None = None,
    region: str | None = None,
    include_bed: str | None = None,
    exclude_bed: str | None = None,
    cpg: bool = False,
    motif: str | None = None,
    motif_offset: int = 0,
    annotate_motifs: bool = False,
    mask: bool = False,
    ignore: str | None = None,
    edge_filter: str | None = None,
    invert_edge_filter: bool = False,
    mapped_only: bool = False,
    allow_non_primary: bool = False,
    num_reads: int | None = None,
    kmer_size: int | None = None,
    filter_percentile: float | None = None,
    filter_threshold: float | None = None,
    mod_thresholds: Sequence[str] | None = None,
    no_filtering: bool = False,
    pass_only: bool = False,
    sampling_frac: float | None = None,
    mod_codes: Sequence[str] | None = None,
    ignore_implicit: bool = False,
    ignore_index: bool = False,
    bgzf: bool = False,
    no_headers: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract read-level base-modification information from a modBAM into a table (``modkit extract``).

    ``table`` selects the subcommand: ``"full"`` (``extract full``) gives the probability of every modification
    call per read position; ``"calls"`` (``extract calls``) gives hard calls per read position using the pass threshold
    of ``modkit_pileup`` (``pass_only`` keeps only passing calls); ``"read-stats"`` (``extract read-stats``) gives one row per
    read with modified/unmodified counts per code (``mod_codes``, such as ``['C:m']``, is required; only
    ``filter_threshold`` and ``mod_thresholds`` apply). Options of another table raise ConfigurationError. Result keys: ``row_count``,
    ``columns``, ``rows`` (a preview of the first rows).
    """
    kind = _cli.choice(table, _EXTRACT_TABLES, "table")
    bam = _file(bam_path, "bam_path")
    out = _new_output(output_path)
    given = {
        "reference": reference,
        "region": region,
        "include_bed": include_bed,
        "exclude_bed": exclude_bed,
        "motif": motif,
        "cpg": cpg,
        "annotate_motifs": annotate_motifs,
        "mask": mask,
        "ignore": ignore,
        "edge_filter": edge_filter,
        "invert_edge_filter": invert_edge_filter,
        "mapped_only": mapped_only,
        "allow_non_primary": allow_non_primary,
        "num_reads": num_reads,
        "kmer_size": kmer_size,
        "bgzf": bgzf,
        "no_headers": no_headers,
        "ignore_implicit": ignore_implicit,
        "ignore_index": ignore_index,
        "filter_percentile": filter_percentile,
        "pass_only": pass_only,
        "no_filtering": no_filtering,
        "sampling_frac": sampling_frac,
    }
    used = {name for name, value in given.items() if value not in (None, False)}
    if kind == "read-stats":
        bad = sorted(used)
        if not mod_codes:
            raise ConfigurationError("table='read-stats' needs mod_codes such as ['C:m'].")
        if bad:
            raise ConfigurationError(f"table='read-stats' does not take: {', '.join(bad)}.")
    elif kind == "full":
        bad = sorted(used & _EXTRACT_CALLS_ONLY)
        if bad or filter_threshold is not None or mod_thresholds is not None:
            raise ConfigurationError(
                "table='full' does not take the threshold options "
                f"({', '.join(bad + (['filter_threshold'] if filter_threshold is not None else []) + (['mod_thresholds'] if mod_thresholds is not None else []))}); use table='calls'."
            )
    if kind != "read-stats" and mod_codes is not None:
        raise ConfigurationError("mod_codes only applies to table='read-stats'.")
    ref = _optional_file(reference, "reference")
    if (cpg or motif is not None or annotate_motifs) and ref is None:
        raise ConfigurationError("reference is required with cpg=True, a motif or annotate_motifs.")
    args = _Args(bam, str(out))
    if kind == "read-stats":
        args.value("--filter-threshold", _fraction(filter_threshold, "filter_threshold"))
        for entry in _thresholds(mod_thresholds, "mod_thresholds") or ():
            args.value("--mod-thresholds", entry)
        args.many("--mod-codes", _codes(mod_codes, "mod_codes"))
    else:
        args.value("--reference", ref)
        args.switch("--ignore-index", ignore_index, "ignore_index")
        args.switch("--bgzf", bgzf, "bgzf")
        args.switch("--no-headers", no_headers, "no_headers")
        args.value("--kmer-size", _optional_int(kmer_size, "kmer_size", minimum=1))
        args.switch("--mapped-only", mapped_only, "mapped_only")
        args.switch("--allow-non-primary", allow_non_primary, "allow_non_primary")
        if kind == "full":
            args.value("--num-reads", _optional_int(num_reads, "num_reads", minimum=1))
        args.value("--region", _region(region))
        args.value("--include-bed", _optional_file(include_bed, "include_bed"))
        args.value("--exclude-bed", _optional_file(exclude_bed, "exclude_bed"))
        args.value("--edge-filter", _edge_filter(edge_filter))
        args.switch("--invert-edge-filter", invert_edge_filter, "invert_edge_filter")
        args.switch("--ignore-implicit", ignore_implicit, "ignore_implicit")
        _motif_pair(motif, motif_offset, "--motif", args)
        args.switch("--annotate-motifs", annotate_motifs, "annotate_motifs")
        args.switch("--cpg", cpg, "cpg")
        args.switch("--mask", mask, "mask")
        args.value("--ignore", _pattern(ignore, "ignore", _MOD_CODE))
        if kind == "calls":
            args.switch("--pass-only", pass_only, "pass_only")
            _filtering(args, filter_percentile, filter_threshold, mod_thresholds, no_filtering)
            args.value("--sample-num-reads", _optional_int(num_reads, "num_reads", minimum=1))
            args.value("--sampling-frac", _fraction(sampling_frac, "sampling_frac"))
    _threads(args, threads)
    _check_executable(modkit_path)
    rendered = _modkit(
        ["extract", kind], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish(
        "modkit_extract",
        rendered,
        table=kind,
        parse=lambda: _p10.table_preview(
            out, delimiter="," if kind == "read-stats" else "\t", header=not no_headers, comment="#"
        ),
    )


# ---------------------------------------------------------------------------
# dmr
# ---------------------------------------------------------------------------


def modkit_dmr(
    samples: dict[str, str],
    regions_bed: str,
    reference: str,
    output_path: str,
    *,
    mode: str = "pair",
    base: Sequence[str] = ("C",),
    assign_code: Sequence[str] | None = None,
    single_code: str | None = None,
    mask: bool = False,
    min_valid_coverage: int | None = None,
    segment: str | None = None,
    reference_sample: str | None = None,
    prefix: str | None = None,
    header: bool = False,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Differential methylation test over regions between samples (``modkit dmr pair`` and ``modkit dmr multi``).

    ``samples`` maps a sample name to a bgzip bedMethyl with a tabix index. ``mode="pair"`` needs exactly two
    samples (the first is the control ``-a``, the second the experiment ``-b``) and writes the BED file
    ``output_path``; ``mode="multi"`` compares all pairs of two or more samples (``reference_sample`` restricts
    the comparisons to those against that sample) and writes one BED per pair into the new directory
    ``output_path`` (file names start with ``prefix``). ``regions_bed`` holds the regions to test (not needed
    with ``segment``, a BED output path for modkit's segmentation, only in pair mode). ``base`` is the
    primary base(s) to analyse (modkit default C); ``assign_code`` entries such as ``h:m`` assign codes to
    the compared modification; ``single_code`` compares a single code. Result keys: ``record_count`` (pair) or
    ``files`` (multi).
    """
    kind = _cli.choice(mode, ("pair", "multi"), "mode")
    if not isinstance(samples, dict) or not samples:
        raise ConfigurationError("samples must be a non-empty {name: bedmethyl_path} mapping.")
    paths = {_name(name, "sample name"): _file(path, "samples") for name, path in samples.items()}
    if kind == "pair" and len(paths) != 2:
        raise ConfigurationError("mode='pair' needs exactly two samples.")
    if kind == "multi" and len(paths) < 2:
        raise ConfigurationError("mode='multi' needs at least two samples.")
    regions = _file(regions_bed, "regions_bed")
    ref = _file(reference, "reference")
    if kind == "pair":
        out = _new_output(output_path)
        directory = False
    else:
        if segment is not None:
            raise ConfigurationError("segment only applies to mode='pair'.")
        out = _new_dir(output_path, "output_path", modkit_path)
        directory = True
    args = _Args()
    if kind == "pair":
        (control, experiment) = list(paths.values())
        args.value("-a", control)
        args.value("-b", experiment)
        args.value("--out-path", out)
        if segment is not None:
            args.value("--segment", _new_output(segment, "segment"))
        args.switch("--header", header, "header")
    else:
        for name, path in paths.items():
            args.items += ["--sample", path, name]
        args.value("--out-dir", out)
        args.value("--prefix", None if prefix is None else _name(prefix, "prefix"))
        if reference_sample is not None and reference_sample not in paths:
            raise ConfigurationError("reference_sample must be one of the sample names.")
        args.value("--ref-sample", reference_sample)
        args.switch("--header", header, "header")
    args.value("--regions-bed", regions)
    args.value("--ref", ref)
    for code in _codes(list(base), "base") or ():
        args.value("--base", code)
    for entry in _codes(assign_code, "assign_code") or ():
        args.value("--assign-code", entry)
    args.value("--single-code", _pattern(single_code, "single_code", _MOD_CODE))
    args.switch("--mask", mask, "mask")
    args.value(
        "--min-valid-coverage", _optional_int(min_valid_coverage, "min_valid_coverage", minimum=1)
    )
    _threads(args, threads)
    _check_executable(modkit_path)
    if directory:
        rendered = _modkit(["dmr", "multi"], args, cwd=out, expected=[], modkit_path=modkit_path)
        files = sorted(str(path) for path in out.iterdir())
        rendered["outputs"] = files
        if rendered["return_code"] == 0 and not files:
            rendered["status"] = "failed"
        return _finish("modkit_dmr", rendered, mode=kind, files=files)
    expected = [out] + (
        [Path(args.items[args.items.index("--segment") + 1])] if segment is not None else []
    )
    rendered = _modkit(
        ["dmr", "pair"], args, cwd=out.parent, expected=expected, modkit_path=modkit_path
    )
    return _finish(
        "modkit_dmr",
        rendered,
        mode=kind,
        parse=lambda: {"record_count": _p10.count_records(out, comment="#")},
    )


def modkit_dmr_isoform(
    bedmethyl_path: str,
    gtf_path: str,
    output_path: str,
    *,
    gene_id: str | None = None,
    gene_name: str | None = None,
    gene_strand: str | None = None,
    gene_chrom: str | None = None,
    min_valid_coverage: int | None = None,
    single_code: str | None = None,
    full: bool = False,
    ignore_version: bool = False,
    plot_path: str | None = None,
    min_score: float | None = None,
    max_pval: float | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare methylation across the transcript isoforms of genes from a transcriptome-aligned bedMethyl (``modkit dmr isoform``).

    ``bedmethyl_path`` is bgzip with a tabix index, ``gtf_path`` the annotation. ``gene_id`` / ``gene_name`` (with
    ``gene_strand`` and ``gene_chrom``) select one gene; otherwise all genes are tested. ``full=True`` writes the
    full per-site table. ``plot_path`` optionally writes modkit's isoform plot (SVG) for the selected gene.
    Result keys: ``row_count``, ``columns``, ``rows``, ``plot``.
    """
    bed = _file(bedmethyl_path, "bedmethyl_path")
    gtf = _file(gtf_path, "gtf_path")
    out = _new_output(output_path)
    plot = None if plot_path is None else _new_output(plot_path, "plot_path")
    args = _Args(bed, str(out))
    args.value("--gtf", gtf)
    args.value(
        "--min-valid-coverage", _optional_int(min_valid_coverage, "min_valid_coverage", minimum=1)
    )
    args.value("--single-code", _pattern(single_code, "single_code", _MOD_CODE))
    args.switch("--full", full, "full")
    args.switch("--ignore-version", ignore_version, "ignore_version")
    args.value("--gene-id", None if gene_id is None else _cli.text(gene_id, "gene_id"))
    args.value("--gene-name", None if gene_name is None else _cli.text(gene_name, "gene_name"))
    args.value(
        "--gene-strand",
        None if gene_strand is None else _cli.choice(gene_strand, ("+", "-"), "gene_strand"),
    )
    args.value("--gene-chrom", None if gene_chrom is None else _cli.text(gene_chrom, "gene_chrom"))
    _threads(args, threads)
    args.value("--plot", plot)
    if min_score is not None:
        args.value("--min-score", _cli.number(min_score, "min_score"))
    args.value(
        "--max-pval",
        None if max_pval is None else _cli.number(max_pval, "max_pval", minimum=0.0, maximum=1.0),
    )
    _check_executable(modkit_path)
    expected = [out] if plot is None else [out, plot]
    rendered = _modkit(
        ["dmr", "isoform"], args, cwd=out.parent, expected=expected, modkit_path=modkit_path
    )
    return _finish(
        "modkit_dmr_isoform",
        rendered,
        plot=None if plot is None else str(plot),
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


def modkit_dmr_compare_tx_sites(
    condition_a: str,
    condition_b: str,
    gtf_path: str,
    output_path: str,
    *,
    min_valid_coverage: int | None = None,
    single_code: str | None = None,
    full: bool = False,
    ignore_version: bool = False,
    plot_path: str | None = None,
    top_k: int | None = None,
    min_effect_size: float | None = None,
    plot_title: str | None = None,
    threads: int | None = None,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare per-site methylation of transcripts between two conditions (``modkit dmr compare-tx-sites``; alias ``compare-transcriptome-sites``).

    ``condition_a`` and ``condition_b`` are transcriptome-aligned bgzip bedMethyl files with tabix indexes and
    ``gtf_path`` the annotation. The per-gene/transcript comparison table goes to ``output_path`` (``full=True`` adds all
    columns). ``plot_path`` optionally writes a volcano-style plot of the ``top_k`` genes with at least
    ``min_effect_size``. Result keys: ``row_count``, ``columns``, ``rows``, ``plot``.
    """
    a = _file(condition_a, "condition_a")
    b = _file(condition_b, "condition_b")
    gtf = _file(gtf_path, "gtf_path")
    out = _new_output(output_path)
    plot = None if plot_path is None else _new_output(plot_path, "plot_path")
    args = _Args()
    args.value("-a", a)
    args.value("-b", b)
    args.value("--out", out)
    args.switch("--full", full, "full")
    args.value("--gtf", gtf)
    args.switch("--ignore-version", ignore_version, "ignore_version")
    args.value(
        "--min-valid-coverage", _optional_int(min_valid_coverage, "min_valid_coverage", minimum=1)
    )
    args.value("--single-code", _pattern(single_code, "single_code", _MOD_CODE))
    _threads(args, threads)
    args.value("--plot", plot)
    args.value("--top-k", _optional_int(top_k, "top_k", minimum=1))
    if min_effect_size is not None:
        args.value(
            "--min-effect-size", _cli.number(min_effect_size, "min_effect_size", minimum=0.0)
        )
    args.value("--plot-title", None if plot_title is None else _cli.text(plot_title, "plot_title"))
    _check_executable(modkit_path)
    expected = [out] if plot is None else [out, plot]
    rendered = _modkit(
        ["dmr", "compare-tx-sites"],
        args,
        cwd=out.parent,
        expected=expected,
        modkit_path=modkit_path,
    )
    return _finish(
        "modkit_dmr_compare_tx_sites",
        rendered,
        plot=None if plot is None else str(plot),
        parse=lambda: _p10.table_preview(out, header=True, comment=None),
    )


# ---------------------------------------------------------------------------
# open-chromatin
# ---------------------------------------------------------------------------


def modkit_open_chromatin_predict(
    bam_path: str,
    model_path: str,
    output_path: str,
    *,
    region: str | None = None,
    include_bed: str | None = None,
    min_coverage: int | None = None,
    threshold: float | None = None,
    batch_size: int | None = None,
    step_size: int | None = None,
    no_header: bool = False,
    dry_run: bool = False,
    modkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict open-chromatin regions from exogenous 6mA (Fiber-seq / DAF-style) signal in a modBAM (``modkit open-chromatin predict``).

    ``model_path`` is a modkit open-chromatin model directory (it holds ``model_config.json``; the release archive ships
    one under ``models/``). ``threshold`` is the score cut-off and
    ``min_coverage`` the minimum read coverage. With ``dry_run=True`` modkit only checks the inputs and
    the model and writes the (possibly empty) output. Result keys: ``row_count``.
    """
    bam = _file(bam_path, "bam_path")
    model_directory = Path(model_path).expanduser().resolve()
    if not model_directory.is_dir():
        raise ConfigurationError(
            f"model_path must be an existing model directory: {model_directory}"
        )
    model = str(model_directory)
    out = _new_output(output_path)
    args = _Args(bam)
    args.value("--model", model)
    args.value("--batch-size", _optional_int(batch_size, "batch_size", minimum=1))
    args.value("--step-size", _optional_int(step_size, "step_size", minimum=1))
    args.value("--min-coverage", _optional_int(min_coverage, "min_coverage", minimum=1))
    args.value(
        "--threshold",
        None
        if threshold is None
        else _cli.number(threshold, "threshold", minimum=0.0, maximum=1.0),
    )
    args.value("--output", out)
    args.value("--include-bed", _optional_file(include_bed, "include_bed"))
    args.value("--region", _region(region))
    args.switch("--no-header", no_header, "no_header")
    args.switch("--dryrun", dry_run, "dry_run")
    _check_executable(modkit_path)
    rendered = _modkit(
        ["open-chromatin", "predict"], args, cwd=out.parent, expected=[out], modkit_path=modkit_path
    )
    return _finish(
        "modkit_open_chromatin_predict",
        rendered,
        parse=lambda: {"record_count": _p10.count_records(out, comment="#")},
    )


__all__ = [
    "modkit_bedmethyl_map_to_genome",
    "modkit_bedmethyl_merge",
    "modkit_bedmethyl_to_bigwig",
    "modkit_dmr",
    "modkit_dmr_compare_tx_sites",
    "modkit_dmr_isoform",
    "modkit_entropy",
    "modkit_extract",
    "modkit_localize",
    "modkit_modbam_adjust_mods",
    "modkit_modbam_call_mods",
    "modkit_modbam_check_tags",
    "modkit_modbam_repair",
    "modkit_modbam_sample_probs",
    "modkit_modbam_summary",
    "modkit_modbam_update_tags",
    "modkit_motif_bed",
    "modkit_motif_evaluate",
    "modkit_motif_refine",
    "modkit_motif_search",
    "modkit_open_chromatin_predict",
    "modkit_pileup",
    "modkit_pileup_hemi",
    "modkit_stats",
    "modkit_validate",
]

_register_public_verbose(globals(), tuple(__all__))
