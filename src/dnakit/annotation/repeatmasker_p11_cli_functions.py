"""RepeatMasker program and utility-script facades (pack P11): fixed argument lists, no shell, bounded execution.

Each public function runs one RepeatMasker program (``RepeatMasker``, ``ProcessRepeats``, ``DateRepeats``, ``DupMasker``,
``RepeatProteinMask``) or one script of its ``util/`` folder (4.2.6 is the tested version). The user installs and
configures RepeatMasker (``conda install -c bioconda repeatmasker``, or the tarball from
https://www.repeatmasker.org/RepeatMasker/ followed by ``perl configure``) and either puts the programs on PATH or passes
``repeatmasker_dir``, the installation folder (the programs are looked up there, the ``util`` scripts in its ``util``
subfolder). Perl scripts whose ``#!`` interpreter does not exist (several ``util`` scripts name ``/usr/local/bin/perl``) are
run with ``perl`` from PATH. The programs write their results next to their input, so inputs are copied into a fresh
``output_directory`` (never overwritten) and the results are written there; standard output is streamed to a file or to
``tool.log``. Results are JSON-serialisable dictionaries with status, version, command, ``outputs`` (every file produced),
small parsed results and a bounded log tail. Optional parameters default to ``None`` / ``False`` and are then not
passed, so RepeatMasker's own defaults apply. Searches with ``species`` need the FamDB library configured in the
RepeatMasker installation; searches with ``library`` need only a search engine (RMBlast) and TRF.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p11_cli as _p
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_TESTED_VERSION = "4.2.6"
_CITATION = "https://www.repeatmasker.org/"
_HINT = (
    "Install and configure RepeatMasker 4.2 or later (for example `conda install -c bioconda repeatmasker`, or the "
    "tarball from https://www.repeatmasker.org/RepeatMasker/ followed by `perl configure`), then put its programs on PATH "
    "or pass repeatmasker_dir (the installation folder)."
)
_ENGINES = ("crossmatch", "wublast", "abblast", "ncbi", "rmblast", "hmmer")
_SPEED = {"quick": "-q", "rush": "-qq", "slow": "-s"}
_IS_MODES = {"report": None, "clip": "-is_clip", "skip": "-no_is"}
_MASK_STYLES = {"n": None, "x": "-x", "lowercase_all": "-small", "soft": "-xsmall"}
_VERSION = re.compile(r"\d+(?:\.\d+)+")


# ---------------------------------------------------------------------------
# Locating and running programs
# ---------------------------------------------------------------------------


def _locate(name: str, repeatmasker_dir: str | Path | None, *, util: bool = False) -> str:
    """The script ``name`` in the installation folder (``util/`` subfolder for utilities) or on PATH."""
    if repeatmasker_dir is not None:
        base = Path(repeatmasker_dir).expanduser().resolve()
        candidate = base / "util" / name if util else base / name
        if not candidate.is_file():
            raise BackendUnavailableError(
                f"{name} was not found in the RepeatMasker folder {base}.", hint=_HINT
            )
        return str(candidate)
    return _cli.program(name, explicit=None, hint=_HINT)


def _launcher(script: str) -> list[str]:
    """The command prefix that runs ``script``: itself when its ``#!`` interpreter exists, otherwise ``perl script``."""
    path = Path(script)
    try:
        first = path.open("rb").readline().decode("utf-8", errors="replace")
    except OSError:
        first = ""
    if first.startswith("#!"):
        parts = first[2:].split()
        if parts:
            interpreter = parts[0]
            if interpreter.endswith("/env") and len(parts) > 1:
                if shutil.which(parts[1]) is not None:
                    return [script]
            elif Path(interpreter).is_file():
                return [script]
        if "perl" in first:
            return [
                _cli.program("perl", explicit=None, hint="Install Perl and put it on PATH."),
                script,
            ]
    return [script]


def _repeatmasker_version(repeatmasker_dir: str | Path | None) -> str | None:
    """The RepeatMasker version (``RepeatMasker -v``), or None when the main program cannot be run."""
    try:
        script = _locate("RepeatMasker", repeatmasker_dir)
    except BackendUnavailableError:
        return None
    return _probe_version(_launcher(script))


def _probe_version(command: Sequence[str]) -> str | None:
    """Run ``<command> -v`` and read the version number from its output."""
    try:
        completed = subprocess.run(
            [*command, "-v"],
            capture_output=True,
            text=True,
            timeout=120,
            stdin=subprocess.DEVNULL,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = _VERSION.search((completed.stdout or "") + (completed.stderr or ""))
    return match.group(0) if match else None


def _run(
    operation: str,
    script: str,
    arguments: Sequence[str],
    directory: Path,
    *,
    repeatmasker_dir: str | Path | None,
    timeout_seconds: float,
    stdout_name: str | None = None,
    before: set[str] | None = None,
    **extra: object,
) -> dict[str, Any]:
    """Run one script inside ``directory`` and return the result with every file that appeared in it."""
    launch = _launcher(script)
    command = [*launch[1:], *arguments] if len(launch) > 1 else list(arguments)
    executable = launch[0]
    _p.run_checked(
        operation,
        executable,
        command,
        directory,
        stdout_name=stdout_name,
        timeout_seconds=timeout_seconds,
    )
    skip = set(before or ()) | {_p.LOG_NAME}
    outputs = {
        path.name: path
        for path in sorted(directory.iterdir())
        if path.is_file() and path.name not in skip
    }
    if stdout_name is not None and stdout_name not in skip:
        outputs[stdout_name] = directory / stdout_name
    return _p.result(
        operation,
        Path(script).name,
        _repeatmasker_version(repeatmasker_dir),
        [Path(script).name, *arguments],
        directory,
        outputs,
        tool="RepeatMasker",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        **extra,
    )


def _stage(source: str | Path, directory: Path, name: str, *, new_name: str | None = None) -> Path:
    """Copy an input file into the output directory (the programs write next to their input)."""
    path = Path(_p.file_path(source, name))
    target = directory / (new_name or path.name)
    if target.exists():
        raise ConfigurationError(
            f"Two inputs would both be copied to {target.name}; rename one of them."
        )
    shutil.copy2(path, target)
    return target


def _listing(directory: Path) -> set[str]:
    return {path.name for path in directory.iterdir()}


def _file_name(path: Path) -> str:
    if path.name.startswith("-"):
        raise ConfigurationError(f"File names must not start with '-': {path.name}")
    return path.name


def _out_table(path: Path) -> dict[str, Any]:
    """Annotation count and per-class counts of a RepeatMasker ``.out`` file."""
    classes: dict[str, int] = {}
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if index < 3 or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 11:
                continue
            total += 1
            classes[fields[10]] = classes.get(fields[10], 0) + 1
    return {"annotations": total, "class_counts": classes}


def _table_summary(path: Path) -> dict[str, Any]:
    """Sequence count, length, GC and masked fraction from the header of a RepeatMasker ``.tbl`` file."""
    text = path.read_text(encoding="utf-8", errors="replace")
    summary: dict[str, Any] = {}
    for key, pattern in (
        ("sequences", r"^sequences:\s+(\d+)"),
        ("total_length", r"^total length:\s+(\d+) bp"),
        ("gc_percent", r"^GC level:\s+([\d.]+) %"),
        ("masked_bases", r"^bases masked:\s+(\d+) bp"),
        ("masked_percent", r"^bases masked:\s+\d+ bp \(\s*([\d.]+) %\)"),
    ):
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            summary[key] = float(match.group(1)) if "." in match.group(1) else int(match.group(1))
    return summary


def _engine_arguments(arguments: list[str], engine: str | None) -> None:
    if engine is not None:
        arguments.extend(("-engine", _cli.choice(engine, _ENGINES, "engine")))


# ---------------------------------------------------------------------------
# RepeatMasker itself
# ---------------------------------------------------------------------------


def mask_repeats(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    species: str | None = None,
    library: str | Path | None = None,
    search_engine: str | None = None,
    speed: str | None = None,
    parallel: int | None = None,
    cutoff: int | None = None,
    mask_low_complexity: bool = True,
    mask_interspersed: bool = True,
    screen_small_rnas: bool = True,
    alu_only: bool = False,
    max_divergence: float | None = None,
    gc_level: int | None = None,
    calculate_gc: bool = False,
    max_fragment_length: int | None = None,
    skip_excision: bool = False,
    insertion_sequences: str = "report",
    mask_style: str = "n",
    lowercase_ambiguous: bool = False,
    alignments: bool = False,
    alignments_in_repeat_orientation: bool = False,
    polymorphic_simple_repeats: bool = False,
    ace_output: bool = False,
    gff_output: bool = False,
    html_output: bool = False,
    untouched_annotation: bool = False,
    crossmatch_output: bool = False,
    no_id_column: bool = False,
    exclude_n_runs: bool = False,
    uncurated_families: bool = False,
    cache_dir: str | Path | None = None,
    famdb_data_dir: str | Path | None = None,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Screen DNA for interspersed repeats and low-complexity sequence: annotation and masked FASTA (RepeatMasker).

    ``species`` (``-species``, a clade or species name resolved with FamDB; RepeatMasker's default is human) or
    ``library`` (``-lib``, a FASTA of repeat consensus sequences, ``>name#class/subclass``) selects the repeats searched;
    ``search_engine`` is ``-engine`` (``crossmatch``, ``wublast``, ``abblast``, ``ncbi``, ``rmblast`` or ``hmmer``).
    ``speed`` is ``quick`` (``-q``), ``rush`` (``-qq``) or ``slow`` (``-s``); ``parallel`` is ``-pa``; ``cutoff`` is the
    score cut-off for custom libraries (``-cutoff``, default 225). What is masked: ``mask_low_complexity=False``
    is ``-nolow``, ``mask_interspersed=False`` is ``-noint`` (only simple repeats and low complexity), ``screen_small_rnas=False``
    is ``-norna``, ``alu_only`` is ``-alu`` and ``max_divergence`` is ``-div``. Search details: ``gc_level`` (``-gc``),
    ``calculate_gc`` (``-gccalc``), ``max_fragment_length`` (``-frag``), ``skip_excision`` (``-nocut``).
    ``insertion_sequences`` controls the bacterial insertion-element check: ``report`` (default), ``clip`` (``-is_clip``,
    clip them out before the analysis) or ``skip`` (``-no_is``). Output: ``mask_style`` ``n`` (repeats replaced by N,
    default), ``x`` (``-x``), ``lowercase_all`` (``-small``) or ``soft`` (``-xsmall``, repeats lower case);
    ``lowercase_ambiguous`` (``-lcambig``); extra files: ``alignments`` (``-a``, ``.align``),
    ``alignments_in_repeat_orientation`` (``-inv``), ``polymorphic_simple_repeats`` (``-poly``), ``ace_output``,
    ``gff_output`` (``-gff``), ``html_output``, ``untouched_annotation`` (``-u``), ``crossmatch_output`` (``-xm``),
    ``no_id_column``, ``exclude_n_runs`` (``-excln``), ``uncurated_families``, ``cache_dir`` and ``famdb_data_dir``.
    Result keys: ``summary`` (sequences, total length, GC, masked bases and percent from ``.tbl``), ``annotations`` and
    ``class_counts`` (from ``.out``) and ``outputs`` (``<name>.out``, ``<name>.masked``, ``<name>.tbl``, ``<name>.cat``
    and the optional files).
    """
    source = Path(_p.file_path(fasta, "fasta"))
    arguments: list[str] = []
    if species is not None and library is not None:
        raise ConfigurationError("Give species or library, not both.")
    _p.opt_text(arguments, "-species", species, "species")
    _p.opt_file(arguments, "-lib", library, "library")
    _engine_arguments(arguments, search_engine)
    if speed is not None:
        arguments.append(_SPEED[_cli.choice(speed, tuple(_SPEED), "speed")])
    _p.opt_int(arguments, "-pa", parallel, "parallel", minimum=1, maximum=1024)
    _p.opt_int(arguments, "-cutoff", cutoff, "cutoff", minimum=1)
    if not _p.flag(mask_low_complexity, "mask_low_complexity"):
        arguments.append("-nolow")
    if not _p.flag(mask_interspersed, "mask_interspersed"):
        arguments.append("-noint")
    if not _p.flag(screen_small_rnas, "screen_small_rnas"):
        arguments.append("-norna")
    _p.switch(arguments, "-alu", alu_only, "alu_only")
    _p.opt_num(arguments, "-div", max_divergence, "max_divergence", minimum=0.0, maximum=100.0)
    _p.opt_int(arguments, "-gc", gc_level, "gc_level", minimum=0, maximum=100)
    _p.switch(arguments, "-gccalc", calculate_gc, "calculate_gc")
    _p.opt_int(arguments, "-frag", max_fragment_length, "max_fragment_length", minimum=1)
    _p.switch(arguments, "-nocut", skip_excision, "skip_excision")
    is_flag = _IS_MODES[_cli.choice(insertion_sequences, tuple(_IS_MODES), "insertion_sequences")]
    if is_flag is not None:
        arguments.append(is_flag)
    style_flag = _MASK_STYLES[_cli.choice(mask_style, tuple(_MASK_STYLES), "mask_style")]
    if style_flag is not None:
        arguments.append(style_flag)
    _p.switch(arguments, "-lcambig", lowercase_ambiguous, "lowercase_ambiguous")
    _p.switch(arguments, "-a", alignments, "alignments")
    _p.switch(
        arguments, "-inv", alignments_in_repeat_orientation, "alignments_in_repeat_orientation"
    )
    _p.switch(arguments, "-poly", polymorphic_simple_repeats, "polymorphic_simple_repeats")
    _p.switch(arguments, "-ace", ace_output, "ace_output")
    _p.switch(arguments, "-gff", gff_output, "gff_output")
    _p.switch(arguments, "-html", html_output, "html_output")
    _p.switch(arguments, "-u", untouched_annotation, "untouched_annotation")
    _p.switch(arguments, "-xm", crossmatch_output, "crossmatch_output")
    _p.switch(arguments, "-no_id", no_id_column, "no_id_column")
    _p.switch(arguments, "-excln", exclude_n_runs, "exclude_n_runs")
    _p.switch(arguments, "-uncurated", uncurated_families, "uncurated_families")
    _p.opt_dir(arguments, "-cache_dir", cache_dir, "cache_dir")
    _p.opt_dir(arguments, "-famdb_data_dir", famdb_data_dir, "famdb_data_dir")
    script = _locate("RepeatMasker", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "mask_repeats",
        script,
        [*arguments, "-dir", str(directory), str(source)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
    )
    outputs = result["outputs"]
    if f"{source.name}.out" in outputs:
        result.update(_out_table(Path(outputs[f"{source.name}.out"])))
    if f"{source.name}.tbl" in outputs:
        result["summary"] = _table_summary(Path(outputs[f"{source.name}.tbl"]))
    return result


def clip_bacterial_insertion_sequences(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip bacterial insertion-sequence (IS) elements out of DNA, for example E. coli contamination (``RepeatMasker -is_only``).

    RepeatMasker searches the sequences for E. coli IS elements (the IS library is part of FamDB), removes the elements
    that are flanked by a short duplicated target site and writes the remaining sequence to ``<name>.withoutIS``; the
    matches that are only reported (not clipped) go to ``<name>.alert``. Nothing else is masked or annotated. Result
    keys: ``clipped_fasta`` (path or None when no IS element was found) and ``outputs``.
    """
    source = Path(_p.file_path(fasta, "fasta"))
    script = _locate("RepeatMasker", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "clip_bacterial_insertion_sequences",
        script,
        ["-is_only", "-dir", str(directory), str(source)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
    )
    result["clipped_fasta"] = result["outputs"].get(f"{source.name}.withoutIS")
    return result


def process_repeat_hits(
    cat_file: str | Path,
    output_directory: str | Path,
    *,
    species: str | None = None,
    library: str | Path | None = None,
    mask_low_complexity: bool = True,
    mask_interspersed: bool = True,
    lowercase_ambiguous: bool = False,
    untouched_annotation: bool = False,
    crossmatch_output: bool = False,
    ace_output: bool = False,
    gff_output: bool = False,
    polymorphic_simple_repeats: bool = False,
    no_id_column: bool = False,
    exclude_n_runs: bool = False,
    align_l1_orf2: bool = False,
    alignments: bool = False,
    mask_source: str | Path | None = None,
    mask_style: str = "n",
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Post-process raw RepeatMasker search hits (``.cat``) into the repeat annotation (``.out``) and summary (``.tbl``).

    ``cat_file`` is the ``.cat`` (or ``.cat.gz``) file of a RepeatMasker run. ProcessRepeats resolves overlapping hits,
    joins fragments of one element and writes the annotation again; use it to recreate or reformat results without
    repeating the search. ``species`` (``-species``, RepeatMasker's default is human) or ``library`` (``-lib``) must match
    the original run. ``mask_low_complexity=False`` is ``-nolow``, ``mask_interspersed=False`` is ``-noint``;
    ``lowercase_ambiguous`` is ``-lcambig``; extra files: ``untouched_annotation`` (``-u``), ``crossmatch_output`` (``-xm``),
    ``ace_output``, ``gff_output`` (``-gff``), ``polymorphic_simple_repeats`` (``-poly``), ``alignments`` (``-a``, needs the
    alignment-bearing ``.cat``); ``no_id_column`` (``-no_id``), ``exclude_n_runs`` (``-excln``), ``align_l1_orf2``
    (``-orf2``). With ``mask_source`` (a FASTA file, ``-maskSource``) the sequence is masked in ``mask_style`` ``n``,
    ``x`` or ``soft``. Result keys: ``annotations`` and ``class_counts`` (from ``.out``) and ``summary`` (from ``.tbl``).
    """
    source = Path(_p.file_path(cat_file, "cat_file"))
    if species is not None and library is not None:
        raise ConfigurationError("Give species or library, not both.")
    arguments: list[str] = []
    _p.opt_text(arguments, "-species", species, "species")
    _p.opt_file(arguments, "-lib", library, "library")
    if not _p.flag(mask_low_complexity, "mask_low_complexity"):
        arguments.append("-nolow")
    if not _p.flag(mask_interspersed, "mask_interspersed"):
        arguments.append("-noint")
    _p.switch(arguments, "-lcambig", lowercase_ambiguous, "lowercase_ambiguous")
    _p.switch(arguments, "-u", untouched_annotation, "untouched_annotation")
    _p.switch(arguments, "-xm", crossmatch_output, "crossmatch_output")
    _p.switch(arguments, "-ace", ace_output, "ace_output")
    _p.switch(arguments, "-gff", gff_output, "gff_output")
    _p.switch(arguments, "-poly", polymorphic_simple_repeats, "polymorphic_simple_repeats")
    _p.switch(arguments, "-no_id", no_id_column, "no_id_column")
    _p.switch(arguments, "-excln", exclude_n_runs, "exclude_n_runs")
    _p.switch(arguments, "-orf2", align_l1_orf2, "align_l1_orf2")
    _p.switch(arguments, "-a", alignments, "alignments")
    _p.opt_file(arguments, "-maskSource", mask_source, "mask_source")
    style = _cli.choice(mask_style, ("n", "x", "soft"), "mask_style")
    if style != "n":
        arguments.append({"x": "-x", "soft": "-xsmall"}[style])
    script = _locate("ProcessRepeats", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "cat_file")
    before = _listing(directory)
    result = _run(
        "process_repeat_hits",
        script,
        [*arguments, _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    stem = re.sub(r"\.cat(\.gz)?$", "", staged.name)
    outputs = result["outputs"]
    if f"{stem}.out" in outputs:
        result.update(_out_table(Path(outputs[f"{stem}.out"])))
    if f"{stem}.tbl" in outputs:
        result["summary"] = _table_summary(Path(outputs[f"{stem}.tbl"]))
    return result


def date_repeats(
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    query_species: str,
    comparison_species: Sequence[str],
    sequence_file: str | Path | None = None,
    mask_species: str | None = None,
    mask_uncertain: bool = False,
    mask_low_complexity: bool = True,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Label repeats as lineage-specific or ancestral relative to other species, and optionally mask the lineage-specific ones.

    ``repeat_annotation`` is the RepeatMasker ``.out`` file of the analysed sequence from ``query_species``;
    ``comparison_species`` are the other (mammalian) species or clades, resolved with FamDB. DateRepeats writes a copy
    of the annotation with one added column per comparison species telling whether the repeat is expected to be present
    there (``<name>.out_<species>...``). With ``mask_species`` (one of the comparison species) it also writes the sequence
    file with all repeats that are in the query but absent from that species masked, for genome alignments; this needs
    ``sequence_file``, the FASTA that the annotation was made from. ``mask_uncertain`` (``-aggressive``) also masks the repeats
    whose age is unclear; ``mask_low_complexity=False`` (``-nolow``) leaves satellites, simple repeats and low-complexity
    DNA unmasked. FamDB must be configured in the RepeatMasker installation.
    """
    annotation = Path(_p.file_path(repeat_annotation, "repeat_annotation"))
    others = (
        [_cli.text(item, "comparison_species") for item in comparison_species]
        if not isinstance(comparison_species, str)
        else []
    )
    if not others:
        raise ConfigurationError("comparison_species must be a non-empty list of species names.")
    arguments = ["-query", _cli.text(query_species, "query_species")]
    for item in others:
        arguments.extend(("-comp", item))
    if mask_species is not None:
        if sequence_file is None:
            raise ConfigurationError(
                "mask_species needs sequence_file (the FASTA of the annotated sequence)."
            )
        arguments.extend(("-mask", _cli.text(mask_species, "mask_species")))
    elif sequence_file is not None:
        raise ConfigurationError("sequence_file is only used together with mask_species.")
    _p.switch(arguments, "-aggressive", mask_uncertain, "mask_uncertain")
    if not _p.flag(mask_low_complexity, "mask_low_complexity"):
        arguments.append("-nolow")
    script = _locate("DateRepeats", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(annotation, directory, "repeat_annotation")
    if sequence_file is not None:
        # DateRepeats looks for the sequence file next to the annotation, named like it without ".out".
        _stage(
            sequence_file, directory, "sequence_file", new_name=re.sub(r"\.out$", "", staged.name)
        )
    before = _listing(directory)
    return _run(
        "date_repeats",
        script,
        [_file_name(staged), *arguments],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )


def annotate_segmental_duplications(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    repeatmasker_out: str | Path | None = None,
    max_divergence: float | None = None,
    max_gap_width: int | None = None,
    force_search: bool = False,
    alignments: bool = False,
    gff_output: bool = False,
    search_engine: str | None = None,
    parallel: int | None = None,
    seeds_only: bool = False,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate segmental duplications (duplicons) by searching DNA against the duplicon library (DupMasker).

    DupMasker repeat-masks the sequence (human libraries), searches the masked sequence against the duplicon library,
    builds realignment regions for each duplicon seed, realigns them on the unmasked DNA and joins the duplicons; the
    annotation is ``<name>.duplicons``. An existing RepeatMasker annotation can be given as ``repeatmasker_out`` so that
    the masking step is not repeated. ``max_divergence`` (``-maxDiv``) filters duplicon seeds, ``max_gap_width``
    (``-maxWidth``, default 300) is the largest non-seed gap that is realigned, ``force_search`` (``-forceSearch``) repeats
    searches whose results exist, ``alignments`` (``-align``) stores alignments, ``gff_output`` (``-gff``) adds a GFF file,
    ``search_engine`` is ``abblast``, ``wublast`` or ``ncbi``, ``parallel`` is ``-pa`` and ``seeds_only`` (``-dupout``) stops
    after the seed search. Result key: ``duplicons`` (path of the annotation or None).
    """
    source = Path(_p.file_path(fasta, "fasta"))
    arguments: list[str] = []
    _p.opt_num(arguments, "-maxDiv", max_divergence, "max_divergence", minimum=0.0)
    _p.opt_int(arguments, "-maxWidth", max_gap_width, "max_gap_width", minimum=1)
    _p.switch(arguments, "-forceSearch", force_search, "force_search")
    _p.switch(arguments, "-align", alignments, "alignments")
    _p.switch(arguments, "-gff", gff_output, "gff_output")
    if search_engine is not None:
        arguments.extend(
            ("-engine", _cli.choice(search_engine, ("abblast", "wublast", "ncbi"), "search_engine"))
        )
    _p.opt_int(arguments, "-pa", parallel, "parallel", minimum=1, maximum=1024)
    _p.switch(arguments, "-dupout", seeds_only, "seeds_only")
    script = _locate("DupMasker", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "fasta")
    if repeatmasker_out is not None:
        _stage(repeatmasker_out, directory, "repeatmasker_out", new_name=f"{staged.name}.out")
    before = _listing(directory)
    result = _run(
        "annotate_segmental_duplications",
        script,
        [*arguments, _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    result["duplicons"] = result["outputs"].get(f"{staged.name}.duplicons")
    return result


def mask_repeat_proteins(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    search_engine: str | None = None,
    p_value: float | None = None,
    min_score: int | None = None,
    word_size: int | None = None,
    query_statistics_length: int | None = None,
    max_amino_acid_distance: int | None = None,
    mask_low_complexity: bool = True,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask repeat-derived protein-coding fragments in DNA by a translated search against repeat proteins (RepeatProteinMask).

    The DNA is compared (BLASTX) with the RepeatPeps database of repeat-encoded proteins; matching regions are masked
    and annotated. Low-complexity and simple repeats are always searched first. ``search_engine`` is ``ncbi``,
    ``abblast`` or ``wublast``; ``p_value`` (``-pvalue``) and ``min_score`` (``-minscore``) are thresholds (none by default);
    ``word_size`` is ``-wordsize`` (default 3), ``query_statistics_length`` is ``-queryStatLen``, ``max_amino_acid_distance``
    is ``-maxAADist`` (default 333) and ``mask_low_complexity=False`` is ``-noLowSimple``. The annotation is ``<name>.annot``
    and the masked sequence ``<name>.masked``; RepeatProteinMask reads ``Libraries/RepeatPeps.lib`` of the installation.
    """
    source = Path(_p.file_path(fasta, "fasta"))
    arguments: list[str] = []
    if search_engine is not None:
        arguments.extend(
            ("-engine", _cli.choice(search_engine, ("ncbi", "abblast", "wublast"), "search_engine"))
        )
    _p.opt_num(arguments, "-pvalue", p_value, "p_value", minimum=0.0)
    _p.opt_int(arguments, "-minscore", min_score, "min_score", minimum=0)
    _p.opt_int(arguments, "-wordsize", word_size, "word_size", minimum=1)
    _p.opt_int(
        arguments, "-queryStatLen", query_statistics_length, "query_statistics_length", minimum=1
    )
    _p.opt_int(
        arguments, "-maxAADist", max_amino_acid_distance, "max_amino_acid_distance", minimum=1
    )
    if not _p.flag(mask_low_complexity, "mask_low_complexity"):
        arguments.append("-noLowSimple")
    script = _locate("RepeatProteinMask", repeatmasker_dir)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "fasta")
    before = _listing(directory)
    return _run(
        "mask_repeat_proteins",
        script,
        [*arguments, _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )


def mask_tandem_repeats(
    fasta: str | Path,
    output_directory: str | Path,
    *,
    soft_mask: bool = False,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask simple tandem repeats in DNA with Tandem Repeats Finder (RepeatMasker's ``util/trfMask``).

    TRF (configured in the RepeatMasker installation) finds the tandem repeats; the masked sequence is written to
    ``<name>.masked`` (repeats replaced by N, or lower case with ``soft_mask=True`` which is ``-xsmall``), with the
    annotation in ``<name>.out`` and the alignments in ``<name>.align``. Result key: ``summary`` (the line trfMask
    prints: sequence count, total length, GC level and masked bases).
    """
    source = Path(_p.file_path(fasta, "fasta"))
    arguments = ["-xsmall"] if _p.flag(soft_mask, "soft_mask") else []
    script = _locate("trfMask", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "fasta")
    before = _listing(directory)
    result = _run(
        "mask_tandem_repeats",
        script,
        [*arguments, _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    match = re.search(
        r"seqCount = (\d+), totaSeqLen = (\d+), GC Level = ([\d.]+), bp masked = (\d+)",
        result["log_tail"],
    )
    if match:
        result["summary"] = {
            "sequences": int(match.group(1)),
            "total_length": int(match.group(2)),
            "gc_percent": float(match.group(3)),
            "masked_bases": int(match.group(4)),
        }
    return result


# ---------------------------------------------------------------------------
# Utilities that read or convert RepeatMasker results
# ---------------------------------------------------------------------------


def convert_repeats_to_gff3(
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "repeats.gff3",
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a RepeatMasker annotation (``.out``) to GFF3 (``util/rmOutToGFF3.pl``).

    Every annotation becomes a ``dispersed_repeat`` feature (source ``RepeatMasker``, score = Smith-Waterman score,
    ``ID`` = the element ID, ``Target`` = the repeat family and its coordinates). Result key: ``features``.
    """
    annotation = _p.file_path(repeat_annotation, "repeat_annotation")
    name = _p.file_stem(output_name, "output_name")
    script = _locate("rmOutToGFF3.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "convert_repeats_to_gff3",
        script,
        [annotation],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        stdout_name=name,
    )
    result["features"] = sum(
        1
        for line in Path(result["outputs"][name]).read_text().splitlines()
        if line and line[0] != "#"
    )
    return result


def convert_repeats_to_bed(
    repeat_file: str | Path,
    output_directory: str | Path,
    *,
    output_prefix: str | None = None,
    split_by: str | None = None,
    sort_by: str | None = None,
    min_length: int | None = None,
    min_hits: int | None = None,
    max_divergence: float | None = None,
    min_divergence: float | None = None,
    overlap_resolution: str | None = None,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a repeat annotation (``.out`` or ``.align``) to BED with class, length and divergence filters (``util/RM2Bed.py``).

    The BED file(s) are named ``<prefix>_rm.bed`` (or ``<prefix>_<group>.bed`` with ``split_by``). Options: ``output_prefix``
    (``-p``, default the first field of the input name), ``split_by`` (``-sp``: ``name``, ``family``, ``class`` or
    ``subclass``), ``sort_by`` (``-s``: ``size``, ``name``, ``family``, ``class``, ``subclass`` or ``diverge``),
    ``min_length`` (``-m``), ``min_hits`` (``-n``, with ``split_by``), ``max_divergence`` (``-dmax``) and ``min_divergence``
    (``-dmin``), and ``overlap_resolution`` (``-o``: ``higher_score``, ``longer_element`` or ``lower_divergence``, which
    assigns every base to one element). The script needs the pandas package. Result key: ``bed_files`` (paths) and
    ``intervals`` (total rows).
    """
    source = Path(_p.file_path(repeat_file, "repeat_file"))
    arguments: list[str] = []
    _p.opt_text(arguments, "-p", output_prefix, "output_prefix")
    if split_by is not None:
        arguments.extend(
            ("-sp", _cli.choice(split_by, ("name", "family", "class", "subclass"), "split_by"))
        )
    if sort_by is not None:
        arguments.extend(
            (
                "-s",
                _cli.choice(
                    sort_by, ("size", "name", "family", "class", "subclass", "diverge"), "sort_by"
                ),
            )
        )
    _p.opt_int(arguments, "-m", min_length, "min_length", minimum=0)
    _p.opt_int(arguments, "-n", min_hits, "min_hits", minimum=1)
    _p.opt_num(arguments, "-dmax", max_divergence, "max_divergence", minimum=0.0)
    _p.opt_num(arguments, "-dmin", min_divergence, "min_divergence", minimum=0.0)
    if overlap_resolution is not None:
        arguments.extend(
            (
                "-o",
                _cli.choice(
                    overlap_resolution,
                    ("higher_score", "longer_element", "lower_divergence"),
                    "overlap_resolution",
                ),
            )
        )
    script = _locate("RM2Bed.py", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "repeat_file")
    before = _listing(directory)
    result = _run(
        "convert_repeats_to_bed",
        script,
        [*arguments, "-d", str(directory), _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    beds = [path for key, path in result["outputs"].items() if key.endswith(".bed")]
    result["bed_files"] = list(beds)
    result["intervals"] = sum(_p.count_lines(Path(path)) for path in beds)
    return result


def export_repeats_ucsc_tables(
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    alignments: str | Path | None = None,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export a repeat annotation as UCSC ``rmsk`` tables (``util/rmToUCSCTables.pl``).

    Writes ``<name>.out.tsv`` (the annotation as the UCSC ``rmsk`` table) and ``<name>.join.tsv`` (the table joined with
    the element IDs) and, when ``alignments`` (the ``.align`` file, optionally gzipped) is given, ``<name>.align.tsv``; the
    alignment table needs ``intersectBed`` (bedtools) on PATH. Result keys: ``tables`` (file name to row count).
    """
    annotation = Path(_p.file_path(repeat_annotation, "repeat_annotation"))
    script = _locate("rmToUCSCTables.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(annotation, directory, "repeat_annotation")
    arguments = ["-out", _file_name(staged)]
    if alignments is not None:
        arguments.extend(("-align", _file_name(_stage(alignments, directory, "alignments"))))
    before = _listing(directory)
    result = _run(
        "export_repeats_ucsc_tables",
        script,
        arguments,
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    result["tables"] = {
        key: _p.count_lines(Path(path))
        for key, path in result["outputs"].items()
        if key.endswith(".tsv")
    }
    return result


def build_repeat_track_hub(
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    genome: str,
    hub_name: str,
    alignments: str | Path | None = None,
    chrom_sizes: str | Path | None = None,
    only_bigbed: bool = False,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a UCSC bigRmsk track hub (bigBed files and hub text files) from a repeat annotation (``util/rmToTrackHub.pl``).

    ``genome`` is the UCSC genome (assembly) name, ``hub_name`` the hub folder name (no spaces), created inside
    ``output_directory``. The script needs the UCSC program ``bedToBigBed`` on PATH. The chromosome sizes are downloaded
    from UCSC (network access) unless ``chrom_sizes`` (a two-column file) is given. ``alignments`` (``.align``) adds the
    alignment bigBed. ``only_bigbed`` (``-only_beds``) writes only the bigBed files, not ``hub.txt``, ``genomes.txt``,
    ``trackDB.txt`` and the HTML stubs. Result keys: ``hub_directory`` and ``hub_files`` (file names in the hub folder).
    """
    annotation = Path(_p.file_path(repeat_annotation, "repeat_annotation"))
    name = _p.file_stem(hub_name, "hub_name")
    assembly = _p.file_stem(genome, "genome")
    arguments = ["-genome", assembly, "-hubname", name]
    sizes = None if chrom_sizes is None else _p.file_path(chrom_sizes, "chrom_sizes")
    if sizes is not None:
        arguments.extend(("-chromsizes", sizes))
    _p.switch(arguments, "-only_beds", only_bigbed, "only_bigbed")
    script = _locate("rmToTrackHub.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(annotation, directory, "repeat_annotation")
    arguments.extend(("-out", _file_name(staged)))
    if alignments is not None:
        arguments.extend(("-align", _file_name(_stage(alignments, directory, "alignments"))))
    before = _listing(directory)
    result = _run(
        "build_repeat_track_hub",
        script,
        arguments,
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    hub = directory / name
    result["hub_directory"] = str(hub)
    result["hub_files"] = sorted(path.name for path in hub.iterdir()) if hub.is_dir() else []
    return result


def extract_repeat_sequences(
    fasta: str | Path,
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "repeats.fa",
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the sequence of every annotated repeat element as FASTA (``util/rmOut2Fasta.pl``).

    ``fasta`` is the file with all sequences named in the RepeatMasker ``.out`` annotation. Each record is named
    ``<sequence>:<start>-<end> ( <repeat>:<repeat start>-<repeat end> orient=<+|-> )`` and is the genomic sequence of that
    annotation, in the orientation of the genome. Result key: ``sequences``.
    """
    sequences = _p.file_path(fasta, "fasta")
    annotation = _p.file_path(repeat_annotation, "repeat_annotation")
    name = _p.file_stem(output_name, "output_name")
    script = _locate("rmOut2Fasta.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "extract_repeat_sequences",
        script,
        ["-fasta", sequences, "-out", annotation],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        stdout_name=name,
    )
    result["sequences"] = _p.count_lines(Path(result["outputs"][name]), ">")
    return result


def mask_sequences_from_annotation(
    fasta: str | Path,
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    soft_mask: bool = False,
    min_divergence: float | None = None,
    max_divergence: float | None = None,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hard- or soft-mask a FASTA file from a RepeatMasker annotation (``util/maskFile.pl``).

    Repeats listed in the ``.out`` file are replaced by N (default) or written in lower case (``soft_mask``, the script's
    ``-softmask``). ``min_divergence`` and ``max_divergence`` (``-minDivergence`` / ``-maxDivergence``, percent) mask only
    the repeats within that divergence range. The result is ``<name>.masked`` in the output folder. Result keys:
    ``masked_fasta``, ``masked_bases`` (N or lower-case bases in the output).
    """
    source = Path(_p.file_path(fasta, "fasta"))
    annotation = _p.file_path(repeat_annotation, "repeat_annotation")
    arguments = ["-annotations", annotation]
    _p.switch(arguments, "-softmask", soft_mask, "soft_mask")
    _p.opt_num(arguments, "-minDivergence", min_divergence, "min_divergence", minimum=0.0)
    _p.opt_num(arguments, "-maxDivergence", max_divergence, "max_divergence", minimum=0.0)
    script = _locate("maskFile.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "fasta")
    before = _listing(directory)
    result = _run(
        "mask_sequences_from_annotation",
        script,
        ["-fasta", _file_name(staged), *arguments],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    masked = Path(result["outputs"][f"{staged.name}.masked"])
    result["masked_fasta"] = str(masked)
    result["masked_bases"] = _count_masked(masked, lowercase=soft_mask)
    return result


def _count_masked(path: Path, *, lowercase: bool) -> int:
    total = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            body = line.strip()
            total += (
                sum(1 for character in body if character.islower())
                if lowercase
                else body.upper().count("N")
            )
    return total


def summarize_repeat_annotation(
    repeat_annotation: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "summary.txt",
    min_divergence: float | None = None,
    max_divergence: float | None = None,
    species: str | None = None,
    genome: str | Path | None = None,
    use_absolute_genome_size: bool = False,
    library_dir: str | Path | None = None,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a repeat annotation (``.out``) per family and class: counts, bases and genome fraction (``util/buildSummary.pl``).

    ``min_divergence`` and ``max_divergence`` (``-minDiv`` / ``-maxDiv``, a fraction 0-1) filter the annotations.
    ``species`` (``-species``) splits families into lineage-specific and ancestral (for runs made with ``-species``; needs
    FamDB). ``genome`` (``-genome``) is a ``.2bit`` file or a two-column ``sequence<TAB>non-N size`` table that sets the
    genome size and the sequences counted; ``use_absolute_genome_size`` is ``-useAbsoluteGenomeSize``. ``library_dir`` is
    ``-libdir``. buildSummary.pl requires the UCSC program ``twoBitInfo`` on PATH. Result key: ``report`` (the first
    4000 characters of the summary).
    """
    annotation = _p.file_path(repeat_annotation, "repeat_annotation")
    name = _p.file_stem(output_name, "output_name")
    arguments: list[str] = []
    _p.opt_num(arguments, "-minDiv", min_divergence, "min_divergence", minimum=0.0, maximum=1.0)
    _p.opt_num(arguments, "-maxDiv", max_divergence, "max_divergence", minimum=0.0, maximum=1.0)
    _p.opt_text(arguments, "-species", species, "species")
    _p.opt_file(arguments, "-genome", genome, "genome")
    _p.switch(
        arguments, "-useAbsoluteGenomeSize", use_absolute_genome_size, "use_absolute_genome_size"
    )
    _p.opt_dir(arguments, "-libdir", library_dir, "library_dir")
    script = _locate("buildSummary.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "summarize_repeat_annotation",
        script,
        [*arguments, annotation],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        stdout_name=name,
    )
    result["report"] = Path(result["outputs"][name]).read_text(encoding="utf-8", errors="replace")[
        :4000
    ]
    return result


def recalculate_repeat_divergence(
    alignment_file: str | Path,
    output_directory: str | Path,
    *,
    summary_name: str = "repeats.divsum",
    write_alignments: bool = False,
    adjust_cpg: bool = True,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recalculate the divergence of every repeat alignment and summarise it per family (``util/calcDivergenceFromAlign.pl``).

    ``alignment_file`` is the ``.align`` (optionally gzipped) file of a RepeatMasker run made with ``-a``. The divergence is
    the Kimura 2-parameter distance, with CpG sites of the consensus counted specially (two transitions as one, one
    transition as 1/10) unless ``adjust_cpg=False`` (``-noCpGMod``). The summary table (``-s``, the input of
    ``plot_repeat_landscape``) is written to ``summary_name``; ``write_alignments`` (``-a``) also writes
    ``<name>.with_divergence.align`` with the divergence on each alignment. Result key: ``families``, one row per
    repeat family with ``class``, ``repeat``, ``absolute_length``, ``well_characterised_length`` and ``kimura_percent``.
    """
    source = Path(_p.file_path(alignment_file, "alignment_file"))
    summary = _p.file_stem(summary_name, "summary_name")
    script = _locate("calcDivergenceFromAlign.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "alignment_file")
    arguments = ["-s", summary]
    if _p.flag(write_alignments, "write_alignments"):
        arguments.extend(("-a", "with_divergence.align"))
    if not _p.flag(adjust_cpg, "adjust_cpg"):
        arguments.append("-noCpGMod")
    before = _listing(directory)
    result = _run(
        "recalculate_repeat_divergence",
        script,
        [*arguments, _file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    families: list[dict[str, Any]] = []
    in_table = False
    for line in (
        Path(result["outputs"][summary]).read_text(encoding="utf-8", errors="replace").splitlines()
    ):
        if line.startswith("Class\tRepeat"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.strip() or line.startswith("-----"):
            if families:
                break
            continue
        fields = line.split("\t")
        if len(fields) >= 5:
            families.append(
                {
                    "class": fields[0],
                    "repeat": fields[1],
                    "absolute_length": _p.number(fields[2]),
                    "well_characterised_length": _p.number(fields[3]),
                    "kimura_percent": _p.number(fields[4]),
                }
            )
    result["families"] = families
    return result


def plot_repeat_landscape(
    divergence_summary: str | Path,
    output_directory: str | Path,
    *,
    genome_size: int | None = None,
    two_bit_file: str | Path | None = None,
    title: str | None = None,
    javascript_only: bool = False,
    output_name: str = "landscape.html",
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Plot the repeat landscape: genome fraction occupied by each repeat class at each divergence level (``util/createRepeatLandscape.pl``).

    ``divergence_summary`` is the ``.divsum`` file of ``recalculate_repeat_divergence``. The genome size for the
    percentages is ``genome_size`` (``-g``, bp) or taken from ``two_bit_file`` (``-twoBit``, needs the UCSC program
    ``twoBitInfo``). The plot is a stacked-area Google Charts page written to ``output_name`` (HTML; the page loads
    the charting library from google.com when opened) or only its JavaScript with ``javascript_only`` (``-j``); ``title``
    is the graph title. Result keys: ``classes`` (repeat classes) and ``landscape`` (one row per divergence level with
    ``divergence`` and the percentage of the genome per class), which is the data behind the plot.
    """
    summary = _p.file_path(divergence_summary, "divergence_summary")
    if (genome_size is None) == (two_bit_file is None):
        raise ConfigurationError("Give exactly one of genome_size and two_bit_file.")
    name = _p.file_stem(output_name, "output_name")
    arguments = ["-div", summary]
    if genome_size is not None:
        arguments.extend(("-g", str(_cli.integer(genome_size, "genome_size", minimum=1))))
    else:
        arguments.extend(("-twoBit", _p.file_path(two_bit_file, "two_bit_file")))  # type: ignore[arg-type]
    _p.opt_text(arguments, "-t", title, "title")
    _p.switch(arguments, "-j", javascript_only, "javascript_only")
    script = _locate("createRepeatLandscape.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    result = _run(
        "plot_repeat_landscape",
        script,
        arguments,
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        stdout_name=name,
    )
    page = Path(result["outputs"][name]).read_text(encoding="utf-8", errors="replace")
    classes = re.findall(r"data\.addColumn\('number', '([^']*)'\)", page)
    rows = []
    for match in re.finditer(r"^\s*\[\s*'([^']*)',([^\]]*)\],?\s*$", page, re.MULTILINE):
        values = [float(item) for item in match.group(2).split(",") if item.strip()]
        if len(values) == len(classes):
            rows.append(
                {"divergence": _p.number(match.group(1)), **dict(zip(classes, values, strict=True))}
            )
    result["classes"] = classes
    result["landscape"] = rows
    return result


def plot_duplicon_map(
    duplicons_file: str | Path,
    output_directory: str | Path,
    *,
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw the DupMasker duplicon map of every sequence as SVG (``util/dupliconToSVG.pl``).

    ``duplicons_file`` is the ``.duplicons`` annotation written by ``annotate_segmental_duplications`` (the script works
    with the human duplicon library). One SVG file per sequence is written, named ``<file>.<n>.svg``. Result key:
    ``svg_files`` (paths).
    """
    source = Path(_p.file_path(duplicons_file, "duplicons_file"))
    script = _locate("dupliconToSVG.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    staged = _stage(source, directory, "duplicons_file")
    before = _listing(directory)
    result = _run(
        "plot_duplicon_map",
        script,
        [_file_name(staged)],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    result["svg_files"] = [path for key, path in result["outputs"].items() if key.endswith(".svg")]
    return result


def combine_repeat_annotations(
    first_annotation: str | Path,
    first_alignments: str | Path,
    second_annotation: str | Path,
    second_alignments: str | Path,
    output_directory: str | Path,
    *,
    output_prefix: str = "combined",
    repeatmasker_dir: str | Path | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge two RepeatMasker annotation sets (``.out`` and ``.align``) into one sorted set (``util/combineRMFiles.pl``).

    Each set is a RepeatMasker ``.out`` annotation with its ``.align`` alignment file (a run with ``-a``). The second set
    must come from searching sequence pre-masked with the first, so the annotations do not overlap; element IDs are
    renumbered and the entries are placed in sequence and position order. The result is ``<output_prefix>.out`` and
    ``<output_prefix>.align``. Result keys: ``annotations`` and ``class_counts`` of the combined ``.out``.
    """
    prefix = _p.file_stem(output_prefix, "output_prefix")
    script = _locate("combineRMFiles.pl", repeatmasker_dir, util=True)
    directory = _cli.new_directory(output_directory, "output_directory")
    _stage(first_annotation, directory, "first_annotation", new_name="first.out")
    _stage(first_alignments, directory, "first_alignments", new_name="first.align")
    _stage(second_annotation, directory, "second_annotation", new_name="second.out")
    _stage(second_alignments, directory, "second_alignments", new_name="second.align")
    before = _listing(directory)
    result = _run(
        "combine_repeat_annotations",
        script,
        ["first", "second", prefix],
        directory,
        repeatmasker_dir=repeatmasker_dir,
        timeout_seconds=timeout_seconds,
        before=before,
    )
    if f"{prefix}.out" in result["outputs"]:
        result.update(_out_table(Path(result["outputs"][f"{prefix}.out"])))
    return result


__all__ = [
    "annotate_segmental_duplications",
    "build_repeat_track_hub",
    "clip_bacterial_insertion_sequences",
    "combine_repeat_annotations",
    "convert_repeats_to_bed",
    "convert_repeats_to_gff3",
    "date_repeats",
    "export_repeats_ucsc_tables",
    "extract_repeat_sequences",
    "mask_repeat_proteins",
    "mask_repeats",
    "mask_sequences_from_annotation",
    "mask_tandem_repeats",
    "plot_duplicon_map",
    "plot_repeat_landscape",
    "process_repeat_hits",
    "recalculate_repeat_divergence",
    "summarize_repeat_annotation",
]
_register_public_verbose(globals())
