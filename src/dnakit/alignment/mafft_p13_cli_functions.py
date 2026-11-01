"""MAFFT command-line facades: strategy-specific multiple alignment, adding sequences to an alignment, seed
alignments, profile-profile alignment, partition-tree alignment of large sets, sparse-core alignment and pairwise
distances.

The user installs MAFFT (``conda install -c bioconda mafft``, ``apt install mafft`` or the release package); 7.526 is
the tested version. Each function runs one MAFFT program with an explicit typed signature mapped to fixed flags; no
shell is used and the program is resolved from an explicit path, then PATH. MAFFT writes the alignment to standard
output, which is streamed to a **new** output file (an existing file is never overwritten); its progress messages
(standard error) are kept in ``log_tail``. Results are JSON-serialisable dictionaries: ``status``, ``backend``,
``version``, ``command``, ``return_code``, ``outputs``, ``log_tail`` and the number of aligned sequences.

The default ``mafft`` strategy is already available as ``dnakit.alignment.align_sequences(method="mafft")``; the
strategies here are the named MAFFT programs (``fftnsi``, ``nwns``, ``nwnsi``, ``linsi``, ``ginsi``, ``einsi``,
``qinsi``, ``xinsi``) that MAFFT installs next to ``mafft``.
"""

from __future__ import annotations

import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p13_cli as _p13
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = (
    "Install MAFFT 7.5 or later (for example `conda install -c bioconda mafft` or `apt install mafft`), then put it "
    "on PATH or pass mafft_path."
)
_STRATEGIES = ("fftnsi", "nwns", "nwnsi", "linsi", "ginsi", "einsi", "qinsi", "xinsi")
_ADD_MODES = ("add", "addfull", "addlong", "addfragments")
_PARTITION_TREES = ("parttree", "dpparttree", "fastaparttree")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _program(name: str, mafft_path: str | None) -> str:
    """Resolve ``name``; with ``mafft_path`` (the ``mafft`` executable) the sibling program in its folder is used."""
    if mafft_path is not None:
        base = Path(mafft_path).expanduser()
        if not base.is_file():
            raise BackendUnavailableError(f"mafft executable does not exist: {base}", hint=_HINT)
        sibling = base.parent / name
        if sibling.is_file():
            return str(sibling)
    return _cli.program(name, explicit=None, hint=_HINT)


def _common_options(
    threads: int,
    sequence_type: str | None,
    gap_open_penalty: float | None,
    gap_extension_offset: float | None,
    max_iterations: int | None,
    retree: int | None,
    reorder: bool,
    adjust_direction: str | None,
    output_format: str,
) -> list[str]:
    args = ["--quiet", "--thread", str(_cli.integer(threads, "threads", minimum=-1, maximum=256))]
    if sequence_type is not None:
        args.append({"nucleotide": "--nuc", "amino": "--amino"}[_cli.choice(sequence_type, ("nucleotide", "amino"), "sequence_type")])
    if gap_open_penalty is not None:
        args += ["--op", f"{_cli.number(gap_open_penalty, 'gap_open_penalty', minimum=0.0):g}"]
    if gap_extension_offset is not None:
        args += ["--ep", f"{_cli.number(gap_extension_offset, 'gap_extension_offset', minimum=0.0):g}"]
    if max_iterations is not None:
        args += ["--maxiterate", str(_cli.integer(max_iterations, "max_iterations", minimum=0))]
    if retree is not None:
        args += ["--retree", str(_cli.integer(retree, "retree", minimum=1))]
    if not isinstance(reorder, bool):
        raise ConfigurationError("reorder must be True or False.")
    if reorder:
        args.append("--reorder")
    if adjust_direction is not None:
        args.append(
            {"fast": "--adjustdirection", "accurate": "--adjustdirectionaccurately"}[
                _cli.choice(adjust_direction, ("fast", "accurate"), "adjust_direction")
            ]
        )
    if _cli.choice(output_format, ("fasta", "clustal"), "output_format") == "clustal":
        args.append("--clustalout")
    return args


def _stream(
    name: str,
    arguments: Sequence[str],
    output: Path,
    mafft_path: str | None,
    timeout_seconds: float,
    operation: str,
    *,
    environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    executable = _program(name, mafft_path)
    result = _p13.run_to_file(
        name,
        arguments,
        output=output,
        explicit=executable,
        hint=_HINT,
        version_flags=("--version",),
        timeout_seconds=timeout_seconds,
        environment=environment,
    )
    result["operation"] = operation
    result["sequence_count"] = _sequence_count(output)
    return result


def _sequence_count(path: Path) -> int:
    """Number of aligned sequences: FASTA ``>`` records, or the first block of a CLUSTAL alignment."""
    if not path.is_file():
        return 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if first.startswith("CLUSTAL"):
            names: set[str] = set()
            for line in handle:
                if line.strip() and not line.startswith(" ") and not line.startswith("CLUSTAL"):
                    names.add(line.split()[0])
            return len(names)
        count = 1 if first.startswith(">") else 0
        return count + sum(1 for line in handle if line.startswith(">"))


def _require_existing(mafft_path: str | None, name: str) -> None:
    _program(name, mafft_path)


# ---------------------------------------------------------------------------
# Strategy alignment
# ---------------------------------------------------------------------------


def mafft_align(
    sequences_path: str,
    output_path: str,
    *,
    strategy: str = "linsi",
    threads: int = 1,
    sequence_type: str | None = None,
    gap_open_penalty: float | None = None,
    gap_extension_offset: float | None = None,
    max_iterations: int | None = None,
    retree: int | None = None,
    reorder: bool = False,
    adjust_direction: str | None = None,
    output_format: str = "fasta",
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Multiple sequence alignment with one of MAFFT's named strategies (``mafft-<strategy>``).

    ``strategy`` selects the program and with it the options MAFFT sets: ``"fftnsi"`` (FFT-NS-i: FFT approximation,
    two guide-tree cycles, 2 refinement iterations), ``"nwns"`` (NW-NS-2: progressive alignment with
    Needleman-Wunsch distances, no FFT; use ``retree=1`` for NW-NS-1), ``"nwnsi"`` (NW-NS-i, iterative refinement),
    ``"linsi"`` (L-INS-i: local pairwise alignments, 1000 refinement iterations; the default here, for up to a few
    hundred sequences), ``"ginsi"`` (G-INS-i: global pairwise alignments), ``"einsi"`` (E-INS-i: generalised affine gap
    pairwise alignments), ``"qinsi"`` (Q-INS-i: McCaskill base-pair probabilities with global distances, for RNA) and
    ``"xinsi"`` (X-INS-i: McCaskill probabilities with SCARNA distances, for RNA; needs MAFFT's extensions).
    The optional arguments are the matching MAFFT options and are only passed when given: ``threads``
    (``--thread``), ``sequence_type`` (``"nucleotide"`` for ``--nuc`` or ``"amino"`` for ``--amino``),
    ``gap_open_penalty`` (``--op``), ``gap_extension_offset`` (``--ep``), ``max_iterations`` (``--maxiterate``),
    ``retree`` (``--retree``), ``reorder`` (``--reorder``), ``adjust_direction`` (``"fast"`` or ``"accurate"`` reverse
    complement detection) and ``output_format`` (``"fasta"`` or ``"clustal"``). The sequences stay in input order
    unless ``reorder`` is set.
    """
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    name = f"mafft-{_cli.choice(strategy, _STRATEGIES, 'strategy')}"
    args = _common_options(
        threads, sequence_type, gap_open_penalty, gap_extension_offset, max_iterations, retree, reorder, adjust_direction, output_format
    )
    args.append(str(sequences))
    _require_existing(mafft_path, name)
    out = _cli.new_file(output_path, "output_path")
    result = _stream(name, args, out, mafft_path, timeout_seconds, f"mafft {strategy}")
    result["strategy"] = strategy
    return result


# ---------------------------------------------------------------------------
# Adding sequences, seed alignments, profiles
# ---------------------------------------------------------------------------


def mafft_add_sequences(
    alignment_path: str,
    new_sequences_path: str,
    output_path: str,
    *,
    mode: str = "add",
    keep_length: bool = False,
    six_mer_pair: bool = False,
    threads: int = 1,
    sequence_type: str | None = None,
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert new sequences into an existing alignment without realigning it (``mafft --add`` family).

    ``mode`` is ``"add"`` (``--add``: full-length sequences), ``"addfull"`` (``--addfull``: MAFFT's separate
    full-length mode), ``"addlong"`` (``--addlong``: long sequences that are not full-length, such as contigs, which
    MAFFT aligns in a memory-saving way) or ``"addfragments"`` (``--addfragments``: short fragments such as reads or
    amplicons). In every mode the columns of the existing alignment are kept.
    ``keep_length`` (``--keeplength``) leaves the alignment length unchanged (insertions in the new sequences are
    dropped); ``six_mer_pair`` (``--6merpair``) uses the fast 6-mer based pairwise step recommended for many fragments.
    The output holds the existing sequences first and the new ones after them.
    """
    existing = _cli.existing_file(alignment_path, "alignment_path")
    new = _cli.existing_file(new_sequences_path, "new_sequences_path")
    flag = _cli.choice(mode, _ADD_MODES, "mode")
    args = ["--quiet", "--thread", str(_cli.integer(threads, "threads", minimum=-1, maximum=256))]
    if sequence_type is not None:
        args.append({"nucleotide": "--nuc", "amino": "--amino"}[_cli.choice(sequence_type, ("nucleotide", "amino"), "sequence_type")])
    for on, switch, label in ((keep_length, "--keeplength", "keep_length"), (six_mer_pair, "--6merpair", "six_mer_pair")):
        if not isinstance(on, bool):
            raise ConfigurationError(f"{label} must be True or False.")
        if on:
            args.append(switch)
    args += [f"--{flag}", str(new), str(existing)]
    _require_existing(mafft_path, "mafft")
    out = _cli.new_file(output_path, "output_path")
    result = _stream("mafft", args, out, mafft_path, timeout_seconds, f"mafft --{flag}")
    result["mode"] = flag
    return result


def mafft_seed_alignment(
    sequences_path: str,
    seed_alignment_paths: Sequence[str],
    output_path: str,
    *,
    threads: int = 1,
    sequence_type: str | None = None,
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align sequences while keeping the columns of one or more fixed seed alignments (``mafft --seed``).

    Each seed alignment file (a FASTA alignment) is given to MAFFT as ``--seed <file>``; its sequences are aligned as
    in the seed and the remaining sequences of ``sequences_path`` are aligned to them. The output contains the seed
    sequences and the sequences of ``sequences_path``.
    """
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    if isinstance(seed_alignment_paths, str) or not isinstance(seed_alignment_paths, Sequence) or not seed_alignment_paths:
        raise ConfigurationError("seed_alignment_paths must be a non-empty list of files.")
    seeds = [_cli.existing_file(path, "seed_alignment_paths") for path in seed_alignment_paths]
    args = ["--quiet", "--thread", str(_cli.integer(threads, "threads", minimum=-1, maximum=256))]
    if sequence_type is not None:
        args.append({"nucleotide": "--nuc", "amino": "--amino"}[_cli.choice(sequence_type, ("nucleotide", "amino"), "sequence_type")])
    for seed in seeds:
        args += ["--seed", str(seed)]
    args.append(str(sequences))
    _require_existing(mafft_path, "mafft")
    out = _cli.new_file(output_path, "output_path")
    return _stream("mafft", args, out, mafft_path, timeout_seconds, "mafft --seed")


def mafft_profile_alignment(
    alignment1_path: str,
    alignment2_path: str,
    output_path: str,
    *,
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two existing alignments to each other as profiles (``mafft-profile``).

    The columns of each input alignment are kept intact and gaps are inserted between columns to match the two
    profiles; the output has the sequences of the first alignment followed by those of the second.
    """
    first = _cli.existing_file(alignment1_path, "alignment1_path")
    second = _cli.existing_file(alignment2_path, "alignment2_path")
    _require_existing(mafft_path, "mafft-profile")
    out = _cli.new_file(output_path, "output_path")
    return _stream(
        "mafft-profile", [str(first), str(second)], out, mafft_path, timeout_seconds, "mafft-profile"
    )


# ---------------------------------------------------------------------------
# Large data sets
# ---------------------------------------------------------------------------


def mafft_partition_tree_align(
    sequences_path: str,
    output_path: str,
    *,
    tree_method: str = "parttree",
    partition_size: int | None = None,
    retree: int = 2,
    threads: int = 1,
    sequence_type: str | None = None,
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Progressive alignment guided by a partitioned tree, for tens of thousands of sequences (``mafft --parttree`` family).

    ``tree_method`` is ``"parttree"`` (``--parttree``: fast 6-mer based partition tree), ``"dpparttree"``
    (``--dpparttree``: partition tree whose group sub-alignments use dynamic programming) or ``"fastaparttree"``
    (``--fastaparttree``: partition tree from FASTA-style sequence similarity). ``partition_size`` is
    ``--partsize`` (the number of sequences per partition, MAFFT default 50) and ``retree`` is ``--retree``
    (MAFFT's own recommendation for these methods is 2).
    """
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    method = _cli.choice(tree_method, _PARTITION_TREES, "tree_method")
    args = ["--quiet", "--thread", str(_cli.integer(threads, "threads", minimum=-1, maximum=256))]
    if sequence_type is not None:
        args.append({"nucleotide": "--nuc", "amino": "--amino"}[_cli.choice(sequence_type, ("nucleotide", "amino"), "sequence_type")])
    args += [f"--{method}", "--retree", str(_cli.integer(retree, "retree", minimum=1))]
    if partition_size is not None:
        args += ["--partsize", str(_cli.integer(partition_size, "partition_size", minimum=2))]
    args.append(str(sequences))
    _require_existing(mafft_path, "mafft")
    out = _cli.new_file(output_path, "output_path")
    result = _stream("mafft", args, out, mafft_path, timeout_seconds, f"mafft --{method}")
    result["tree_method"] = method
    return result


def mafft_sparsecore_align(
    sequences_path: str,
    output_path: str,
    *,
    core_size: int = 500,
    core_candidates: str | int = "50%",
    seed: int = 0,
    core_options: str | None = None,
    add_options: str | None = None,
    keep_input_order: bool = False,
    mark_core: bool = False,
    mafft_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Large-scale alignment by aligning a sparse core first and adding the rest (``mafft-sparsecore.rb``).

    A random ``core_size`` subset of the ``core_candidates`` (the upper 50 % by length by default, or an integer
    count) is aligned with ``--globalpair --maxiterate 100`` (extra MAFFT options for this stage can be added with
    ``core_options``) and the remaining sequences are added to that core alignment (``add_options``). ``seed`` is the
    random seed of the core choice, ``keep_input_order`` keeps the input order instead of reordering by similarity and
    ``mark_core`` marks the core sequences with ``*`` in their names. Needs ruby.
    """
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    size = _cli.integer(core_size, "core_size", minimum=2)
    random_seed = _cli.integer(seed, "seed", minimum=0)
    candidates = str(core_candidates)
    if not candidates or candidates.startswith("-") or any(char.isspace() for char in candidates):
        raise ConfigurationError("core_candidates must be an integer count or a percentage such as '50%'.")
    mafft = _program("mafft", mafft_path)
    script = _program("mafft-sparsecore.rb", mafft_path)
    if shutil.which("ruby") is None:
        raise BackendUnavailableError(
            "ruby was not found on PATH; mafft-sparsecore.rb is a Ruby script.", hint="Install ruby (apt install ruby)."
        )
    args = ["-i", str(sequences), "-m", mafft, "-p", str(size), "-n", candidates, "-s", str(random_seed)]
    for flag, value, label in (("-C", core_options, "core_options"), ("-A", add_options, "add_options")):
        if value is not None:
            args += [flag, _cli.text(value, label)]
    if not isinstance(keep_input_order, bool) or not isinstance(mark_core, bool):
        raise ConfigurationError("keep_input_order and mark_core must be True or False.")
    if keep_input_order:
        args += ["-o", "i"]
    if mark_core:
        args.append("-M")
    out = _cli.new_file(output_path, "output_path")
    result = _p13.run_to_file(
        "mafft-sparsecore.rb",
        args,
        output=out,
        explicit=script,
        hint=_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["version"] = _p13.version_of(mafft, ("--version",))
    result["operation"] = "mafft-sparsecore.rb"
    result["sequence_count"] = _sequence_count(out)
    return result


# ---------------------------------------------------------------------------
# Pairwise distances
# ---------------------------------------------------------------------------


def mafft_pairwise_distances(
    sequences_path: str,
    output_path: str,
    *,
    sequence_type: str | None = None,
    max_rows_returned: int = 10_000,
    mafft_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise sequence distances computed from pairwise alignments (``mafft-distance``).

    ``output_path`` receives one line ``i-j d=<distance> l=<length i>,<length j>`` per sequence pair (1-based input
    order); the first ``max_rows_returned`` are also returned in ``distances`` as ``i``, ``j``, ``distance``,
    ``length_i`` and ``length_j``. ``sequence_type`` is ``None`` (detected), ``"nucleotide"`` (``-D``) or ``"amino"``
    (``-P``).
    """
    sequences = _cli.existing_file(sequences_path, "sequences_path")
    limit = _cli.integer(max_rows_returned, "max_rows_returned", minimum=1)
    args: list[str] = []
    if sequence_type is not None:
        args.append({"nucleotide": "-D", "amino": "-P"}[_cli.choice(sequence_type, ("nucleotide", "amino"), "sequence_type")])
    args.append(str(sequences))
    _require_existing(mafft_path, "mafft-distance")
    out = _cli.new_file(output_path, "output_path")
    executable = _program("mafft-distance", mafft_path)
    result = _p13.run_to_file(
        "mafft-distance",
        args,
        output=out,
        explicit=executable,
        hint=_HINT,
        version_flags=None,
        timeout_seconds=timeout_seconds,
    )
    result["version"] = _p13.version_of(_program("mafft", mafft_path), ("--version",))
    result["operation"] = "mafft-distance"
    rows: list[dict[str, Any]] = []
    total = 0
    if out.is_file():
        for line in out.read_text(encoding="utf-8", errors="replace").splitlines():
            parts = line.split()
            if len(parts) == 3 and "-" in parts[0] and parts[1].startswith("d=") and parts[2].startswith("l="):
                total += 1
                if total <= limit:
                    first, second = parts[0].split("-", 1)
                    length_i, length_j = parts[2][2:].split(",")
                    rows.append(
                        {
                            "i": int(first),
                            "j": int(second),
                            "distance": float(parts[1][2:]),
                            "length_i": int(length_i),
                            "length_j": int(length_j),
                        }
                    )
    result["pair_count"] = total
    result["distances"] = rows
    return result


__all__ = [
    "mafft_add_sequences",
    "mafft_align",
    "mafft_pairwise_distances",
    "mafft_partition_tree_align",
    "mafft_profile_alignment",
    "mafft_seed_alignment",
    "mafft_sparsecore_align",
]

_register_public_verbose(globals(), tuple(__all__))
