"""pairtools command-line operations as typed DNAKit facades: Hi-C pair (.pairs/.pairsam) processing.

The user installs pairtools (``pip install pairtools``, or ``conda install -c bioconda pairtools``). Each function runs
one pairtools subcommand with an argument list built from explicit typed parameters mapped to fixed flags; no shell is
used; flags are passed only when the caller sets them, so pairtools' own defaults apply; outputs are new explicit
paths; results are JSON-serialisable dictionaries with status, backend, version, command, outputs and a bounded log tail.
"""

from __future__ import annotations

import functools
import importlib
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = (
    "Install pairtools (tested with 1.1.2: python -m pip install pairtools, or conda install -c bioconda pairtools), "
    "or pass pairtools_path."
)
_MAX_TAIL = 4000
_TIMEOUT_SECONDS = 3600.0
_SUFFIX = re.compile(r"^[A-Za-z0-9_.\-]+$")


def _program(pairtools_path: str | None) -> str:
    """Resolve the pairtools executable: explicit path, then PATH, then the interpreter's own bin folder."""
    if pairtools_path is not None:
        candidate = Path(pairtools_path).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"pairtools_path does not exist: {candidate}", hint=_HINT)
        return str(candidate)
    found = shutil.which("pairtools")
    if found is None:
        local = Path(sys.executable).parent / "pairtools"
        if local.is_file():
            found = str(local)
    if found is None:
        raise BackendUnavailableError("pairtools was not found on PATH.", hint=_HINT)
    return found


@functools.lru_cache(maxsize=8)
def _version(program: str) -> str:
    completed = subprocess.run(
        [program, "--version"],
        capture_output=True,
        text=True,
        timeout=120,
        shell=False,
        check=False,
    )
    return (completed.stdout or completed.stderr).strip()


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _flag(name: str, value: Any) -> list[str]:
    """A flag with its value, or nothing when the caller left the pairtools default in place."""
    return [] if value is None else [name, str(value)]


def _switch(name: str, enabled: bool) -> list[str]:
    return [name] if enabled else []


def _run(
    subcommand: Sequence[str],
    arguments: Sequence[str],
    outputs: Sequence[Path],
    pairtools_path: str | None,
) -> dict[str, Any]:
    program = _program(pairtools_path)
    args = [*subcommand, *arguments]
    cwd = outputs[0].parent if outputs else Path.cwd()
    completed = subprocess.run(
        [program, *args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=_TIMEOUT_SECONDS,
        shell=False,
        check=False,
        cwd=cwd,
    )
    produced = [p for p in outputs if p.is_file()]
    status = (
        "completed" if completed.returncode == 0 and len(produced) == len(outputs) else "failed"
    )
    return {
        "status": status,
        "backend": "pairtools",
        "version": _version(program),
        "command": ["pairtools", *args],
        "return_code": completed.returncode,
        "outputs": [str(p) for p in produced],
        "log_tail": (completed.stderr or "")[-_MAX_TAIL:],
    }


def _optional_outputs(*pairs: tuple[str, str | None]) -> tuple[list[str], list[Path]]:
    """Flags for optional new output files, with their resolved paths."""
    flags: list[str] = []
    paths: list[Path] = []
    for name, value in pairs:
        if value is not None:
            target = _new(value, name.lstrip("-").replace("-", "_"))
            flags += [name, str(target)]
            paths.append(target)
    return flags, paths


def _upstream(module: str) -> Any:
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(f"{module} could not be imported: {exc}", hint=_HINT) from exc


def pairtools_dedup(
    pairs_path: str,
    output_path: str,
    *,
    output_dups: str | None = None,
    output_unmapped: str | None = None,
    output_stats: str | None = None,
    max_mismatch: int | None = None,
    method: str | None = None,
    mark_dups: bool = False,
    n_proc: int | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find and remove (or mark) PCR and optical duplicate pairs in a .pairs file (pairtools dedup)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    flags, extra = _optional_outputs(
        ("--output-dups", output_dups),
        ("--output-unmapped", output_unmapped),
        ("--output-stats", output_stats),
    )
    arguments = [
        "-o",
        str(target),
        *flags,
        *_flag("--max-mismatch", max_mismatch),
        *_flag("--method", method),
        *_switch("--mark-dups", mark_dups),
        *_flag("--n-proc", n_proc),
        str(source),
    ]
    return _run(["dedup"], arguments, [target, *extra], pairtools_path)


def pairtools_filterbycov(
    pairs_path: str,
    output_path: str,
    *,
    output_highcov: str | None = None,
    output_unmapped: str | None = None,
    output_stats: str | None = None,
    max_cov: int | None = None,
    max_dist: int | None = None,
    method: str | None = None,
    mark_multi: bool = False,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove pairs from genomic regions of high coverage (pairtools filterbycov)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    flags, extra = _optional_outputs(
        ("--output-highcov", output_highcov),
        ("--output-unmapped", output_unmapped),
        ("--output-stats", output_stats),
    )
    arguments = [
        "-o",
        str(target),
        *flags,
        *_flag("--max-cov", max_cov),
        *_flag("--max-dist", max_dist),
        *_flag("--method", method),
        *_switch("--mark-multi", mark_multi),
        str(source),
    ]
    return _run(["filterbycov"], arguments, [target, *extra], pairtools_path)


def pairtools_flip(
    pairs_path: str,
    chroms_path: str,
    output_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flip pairs so that the matrix is upper-triangular, using chromosome order from a chromsizes file (pairtools flip)."""
    source = _existing(pairs_path, "pairs_path")
    chroms = _existing(chroms_path, "chroms_path")
    target = _new(output_path, "output_path")
    return _run(
        ["flip"], ["-c", str(chroms), "-o", str(target), str(source)], [target], pairtools_path
    )


def pairtools_header_generate(
    pairs_path: str,
    output_path: str,
    *,
    chroms_path: str | None = None,
    sam_path: str | None = None,
    columns: str | None = None,
    assembly: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate or extend the header of a .pairs/.pairsam file from a chromsizes file or a SAM header (pairtools header generate)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    chroms = _existing(chroms_path, "chroms_path") if chroms_path is not None else None
    sam = _existing(sam_path, "sam_path") if sam_path is not None else None
    arguments = [
        "-o",
        str(target),
        *_flag("--chroms-path", chroms),
        *_flag("--sam-path", sam),
        *_flag("--columns", columns),
        *_flag("--assembly", assembly),
        str(source),
    ]
    return _run(["header", "generate"], arguments, [target], pairtools_path)


def pairtools_header_transfer(
    pairs_path: str,
    reference_path: str,
    output_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy the header of a reference .pairs file onto another .pairs file (pairtools header transfer)."""
    source = _existing(pairs_path, "pairs_path")
    reference = _existing(reference_path, "reference_path")
    target = _new(output_path, "output_path")
    return _run(
        ["header", "transfer"],
        ["-o", str(target), "-r", str(reference), str(source)],
        [target],
        pairtools_path,
    )


def pairtools_header_set_columns(
    pairs_path: str,
    columns: str,
    output_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set the #columns header of a .pairs file to the given space-separated column names (pairtools header set-columns)."""
    if not columns.strip():
        raise ConfigurationError("columns must list at least one column name.")
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    return _run(
        ["header", "set-columns"],
        ["-o", str(target), "-c", columns, str(source)],
        [target],
        pairtools_path,
    )


def pairtools_header_validate_columns(
    pairs_path: str,
    reference_path: str,
    output_path: str,
    *,
    reference_columns: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the column layout of a .pairs file matches a reference file (pairtools header validate-columns)."""
    source = _existing(pairs_path, "pairs_path")
    reference = _existing(reference_path, "reference_path")
    target = _new(output_path, "output_path")
    arguments = [
        "-o",
        str(target),
        "-r",
        str(reference),
        *_flag("-c", reference_columns),
        str(source),
    ]
    return _run(["header", "validate-columns"], arguments, [target], pairtools_path)


def pairtools_markasdup(
    pairsam_path: str,
    output_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag every pair of a .pairsam file as a duplicate (pairtools markasdup)."""
    source = _existing(pairsam_path, "pairsam_path")
    target = _new(output_path, "output_path")
    return _run(["markasdup"], ["-o", str(target), str(source)], [target], pairtools_path)


def pairtools_merge(
    pairs_paths: Sequence[str],
    output_path: str,
    *,
    max_nmerge: int | None = None,
    tmpdir: str | None = None,
    memory: str | None = None,
    nproc: int | None = None,
    keep_first_header: bool = False,
    concatenate: bool = False,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several sorted .pairs files into one sorted file (pairtools merge)."""
    if not pairs_paths:
        raise ConfigurationError("pairs_paths must list at least one .pairs file.")
    sources = [str(_existing(p, "pairs_paths")) for p in pairs_paths]
    target = _new(output_path, "output_path")
    arguments = [
        "-o",
        str(target),
        *_flag("--max-nmerge", max_nmerge),
        *_flag("--tmpdir", tmpdir),
        *_flag("--memory", memory),
        *_flag("--nproc", nproc),
        *_switch("--keep-first-header", keep_first_header),
        *_switch("--concatenate", concatenate),
        *sources,
    ]
    return _run(["merge"], arguments, [target], pairtools_path)


def pairtools_parse(
    sam_path: str,
    chroms_path: str,
    output_path: str,
    *,
    assembly: str | None = None,
    min_mapq: int | None = None,
    max_molecule_size: int | None = None,
    add_columns: str | None = None,
    output_stats: str | None = None,
    output_parsed_alignments: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find ligation pairs in name-sorted SAM alignments and write a .pairs file (pairtools parse)."""
    sam = _existing(sam_path, "sam_path")
    chroms = _existing(chroms_path, "chroms_path")
    target = _new(output_path, "output_path")
    flags, extra = _optional_outputs(
        ("--output-stats", output_stats), ("--output-parsed-alignments", output_parsed_alignments)
    )
    arguments = [
        "-c",
        str(chroms),
        "-o",
        str(target),
        *_flag("--assembly", assembly),
        *_flag("--min-mapq", min_mapq),
        *_flag("--max-molecule-size", max_molecule_size),
        *_flag("--add-columns", add_columns),
        *flags,
        str(sam),
    ]
    return _run(["parse"], arguments, [target, *extra], pairtools_path)


def pairtools_parse2(
    sam_path: str,
    chroms_path: str,
    output_path: str,
    *,
    assembly: str | None = None,
    min_mapq: int | None = None,
    single_end: bool = False,
    output_stats: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract ligation pairs from SAM/BAM alignments with the parse2 walk algorithm (pairtools parse2)."""
    sam = _existing(sam_path, "sam_path")
    chroms = _existing(chroms_path, "chroms_path")
    target = _new(output_path, "output_path")
    flags, extra = _optional_outputs(("--output-stats", output_stats))
    arguments = [
        "-c",
        str(chroms),
        "-o",
        str(target),
        *_flag("--assembly", assembly),
        *_flag("--min-mapq", min_mapq),
        *_switch("--single-end", single_end),
        *flags,
        str(sam),
    ]
    return _run(["parse2"], arguments, [target, *extra], pairtools_path)


def pairtools_phase(
    pairs_path: str,
    output_path: str,
    *,
    phase_suffixes: Sequence[str] | None = None,
    min_mapq: int | None = None,
    report_scores: bool = False,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase pairs of a diploid genome using chromosome suffixes such as -1 and -2 (pairtools phase)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    if phase_suffixes is not None and (
        not phase_suffixes or any(not _SUFFIX.match(s) for s in phase_suffixes)
    ):
        raise ConfigurationError(
            "phase_suffixes must be one or more chromosome suffixes such as '-1' and '-2'."
        )
    # The positional input comes first: --phase-suffixes takes several values and would swallow it otherwise.
    arguments = [
        str(source),
        "-o",
        str(target),
        *_flag("--min-mapq", min_mapq),
        *_switch("--report-scores", report_scores),
    ]
    if phase_suffixes is not None:
        arguments += ["--phase-suffixes", *phase_suffixes]
    return _run(["phase"], arguments, [target], pairtools_path)


def pairtools_restrict(
    pairs_path: str,
    frags_path: str,
    output_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign restriction fragments to both ends of each pair from a fragment BED file (pairtools restrict)."""
    source = _existing(pairs_path, "pairs_path")
    frags = _existing(frags_path, "frags_path")
    target = _new(output_path, "output_path")
    return _run(
        ["restrict"], ["-f", str(frags), "-o", str(target), str(source)], [target], pairtools_path
    )


def pairtools_sample(
    pairs_path: str,
    fraction: float,
    output_path: str,
    *,
    seed: int | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep a random fraction of the pairs of a .pairs file (pairtools sample)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    if (
        isinstance(fraction, bool)
        or not isinstance(fraction, (int, float))
        or not 0 < fraction <= 1
    ):
        raise ConfigurationError("fraction must be a number in (0, 1].")
    arguments = [str(fraction), str(source), "-o", str(target), *_flag("--seed", seed)]
    return _run(["sample"], arguments, [target], pairtools_path)


def pairtools_scaling(
    input_path: str,
    output_path: str,
    *,
    view: str | None = None,
    dist_range: str | None = None,
    n_dist_bins_decade: int | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate contact-probability scaling (P(s)) of pairs, optionally restricted to a view (pairtools scaling)."""
    source = _existing(input_path, "input_path")
    target = _new(output_path, "output_path")
    arguments = [
        "-o",
        str(target),
        *_flag("--view", view),
        *_flag("--dist-range", dist_range),
        *_flag("--n-dist-bins-decade", n_dist_bins_decade),
        str(source),
    ]
    return _run(["scaling"], arguments, [target], pairtools_path)


def pairtools_select(
    condition: str,
    pairs_path: str,
    output_path: str,
    *,
    output_rest: str | None = None,
    chrom_subset: str | None = None,
    remove_columns: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep pairs that satisfy a Python condition on their columns, such as chrom1 == chrom2 (pairtools select)."""
    if not condition.strip():
        raise ConfigurationError("condition must be a non-empty expression.")
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    subset = _existing(chrom_subset, "chrom_subset") if chrom_subset is not None else None
    flags, extra = _optional_outputs(("--output-rest", output_rest))
    arguments = [
        condition,
        str(source),
        "-o",
        str(target),
        *flags,
        *_flag("--chrom-subset", subset),
        *_flag("--remove-columns", remove_columns),
    ]
    return _run(["select"], arguments, [target, *extra], pairtools_path)


def pairtools_sort(
    pairs_path: str,
    output_path: str,
    *,
    c1: str | None = None,
    c2: str | None = None,
    p1: str | None = None,
    p2: str | None = None,
    pair_type: str | None = None,
    nproc: int | None = None,
    tmpdir: str | None = None,
    memory: str | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort a .pairs or .pairsam file by chromosome and position columns (pairtools sort)."""
    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    arguments = [
        str(source),
        "-o",
        str(target),
        *_flag("--c1", c1),
        *_flag("--c2", c2),
        *_flag("--p1", p1),
        *_flag("--p2", p2),
        *_flag("--pt", pair_type),
        *_flag("--nproc", nproc),
        *_flag("--tmpdir", tmpdir),
        *_flag("--memory", memory),
    ]
    return _run(["sort"], arguments, [target], pairtools_path)


def pairtools_split(
    pairsam_path: str,
    output_pairs_path: str,
    output_sam_path: str,
    *,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a .pairsam file into a .pairs file and a SAM file (pairtools split)."""
    source = _existing(pairsam_path, "pairsam_path")
    pairs_target = _new(output_pairs_path, "output_pairs_path")
    sam_target = _new(output_sam_path, "output_sam_path")
    arguments = [str(source), "--output-pairs", str(pairs_target), "--output-sam", str(sam_target)]
    return _run(["split"], arguments, [pairs_target, sam_target], pairtools_path)


def pairtools_stats(
    input_paths: Sequence[str],
    output_path: str,
    *,
    merge: bool = False,
    n_dist_bins_decade: int | None = None,
    pairtools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute pair statistics of one .pairs file, or with merge=True merge several stats files (pairtools stats)."""
    if not input_paths:
        raise ConfigurationError("input_paths must list at least one file.")
    if not merge and len(input_paths) != 1:
        raise ConfigurationError("Without merge=True, stats takes exactly one .pairs file.")
    sources = [str(_existing(p, "input_paths")) for p in input_paths]
    target = _new(output_path, "output_path")
    arguments = [
        *sources,
        "-o",
        str(target),
        *_switch("--merge", merge),
        *_flag("--n-dist-bins-decade", n_dist_bins_decade),
    ]
    return _run(["stats"], arguments, [target], pairtools_path)


def pairtools_mark_split_pair_as_dup(columns: Sequence[str], *, verbose: bool = True) -> list[str]:
    """Mark one parsed .pairsam record (list of its columns) as a duplicate: pair type DD and its SAM records tagged."""
    if not columns or any(not isinstance(c, str) for c in columns):
        raise ConfigurationError(
            "columns must be a non-empty list of strings (one .pairsam record)."
        )
    module = _upstream("pairtools.lib.dedup")
    return list(module.mark_split_pair_as_dup(list(columns)))


__all__ = [
    "pairtools_dedup",
    "pairtools_filterbycov",
    "pairtools_flip",
    "pairtools_header_generate",
    "pairtools_header_set_columns",
    "pairtools_header_transfer",
    "pairtools_header_validate_columns",
    "pairtools_mark_split_pair_as_dup",
    "pairtools_markasdup",
    "pairtools_merge",
    "pairtools_parse",
    "pairtools_parse2",
    "pairtools_phase",
    "pairtools_restrict",
    "pairtools_sample",
    "pairtools_scaling",
    "pairtools_select",
    "pairtools_sort",
    "pairtools_split",
    "pairtools_stats",
]
_register_public_verbose(globals())
