"""Primer3 operations not covered elsewhere (tested with primer3-py 2.3.1 / libprimer3 2.6.1 and Primer3 2.6.1):
primer design for every ``PRIMER_TASK`` (check, list, sequencing, cloning, discriminative), local alignment of
oligos with ``ntdpal`` and k-mer based repeat masking with ``primer3_masker``.

Design uses the ``primer3-py`` bindings (extra ``bio-primer3``). ``ntdpal`` and ``primer3_masker`` are Primer3
command-line programs (``conda install -c bioconda primer3``); no shell is used and the argument lists are fixed.
Oligo melting temperature, hairpin and dimer thermodynamics are already provided by ``calculate_primer3_tm``,
``calculate_primer3_hairpin`` and ``calculate_primer3_dimer``.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli, _p14_cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_HINT = "Install Primer3 2.6 or later (for example `conda install -c bioconda primer3`) and put it on PATH or pass primer3_path."
_TASKS = (
    "generic",
    "check_primers",
    "pick_primer_list",
    "pick_sequencing_primers",
    "pick_cloning_primers",
    "pick_discriminative_primers",
)
_SETTING = re.compile(r"^PRIMER_[A-Z0-9_]+$")
_KEY = re.compile(r"^PRIMER_(LEFT|RIGHT|INTERNAL|PAIR)_(\d+)_([A-Z0-9_]+)$")
_DNA = re.compile(r"^[ACGTNacgtn]+$")
_ALIGN_LINE = re.compile(
    r"score=(-?[\d.]+) len=(\d+)(?: align_end_1=(\d+) align_end_2=(\d+))?(?: \|(.*)\|)?"
)


def _region(value: Sequence[int] | None, name: str, length: int) -> list[int] | None:
    if value is None:
        return None
    if isinstance(value, str) or len(value) != 2:
        raise ConfigurationError(f"{name} must be [start, length] (0-based).")
    start = _gap_cli.integer(value[0], f"{name} start", minimum=0)
    size = _gap_cli.integer(value[1], f"{name} length", minimum=1)
    if start + size > length:
        raise ConfigurationError(f"{name} lies outside the template.")
    return [start, size]


def _oligo(value: str | None, name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not _DNA.fullmatch(value):
        raise ConfigurationError(f"{name} must be a string of A, C, G, T, N.")
    return value


def primer3_design_by_task(
    template: str,
    *,
    task: str = "generic",
    left_primer: str | None = None,
    right_primer: str | None = None,
    internal_oligo: str | None = None,
    target: Sequence[int] | None = None,
    included_region: Sequence[int] | None = None,
    excluded_regions: Sequence[Sequence[int]] | None = None,
    product_size_range: Sequence[Sequence[int]] | None = None,
    num_return: int | None = None,
    settings: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run Primer3 for any task with primer3-py (``primer3_core PRIMER_TASK=...``).

    ``task`` is 'generic' (primer pairs), 'check_primers' (evaluate the given ``left_primer`` / ``right_primer`` /
    ``internal_oligo``), 'pick_primer_list' (ranked single primers without pairing), 'pick_sequencing_primers'
    (primers around ``target``), 'pick_cloning_primers' (primers at the ends of ``included_region``) or
    'pick_discriminative_primers' (primers whose 3' end sits on the 1-base ``target``, to tell alleles apart).
    Coordinates are 0-based ``[start, length]``. ``settings`` holds further Primer3 tags (``PRIMER_OPT_TM`` and so
    on; unset tags keep Primer3's defaults and ``PRIMER_TASK`` is set by ``task``). Returns the oligos
    (``left``/``right``/``internal``, each with sequence, 0-based start, length, tm, gc_percent, penalty and the
    thermodynamic fields), the ``pairs`` and Primer3's ``explain`` counters.
    """
    primer3 = require("primer3", "bio-primer3")
    if not isinstance(template, str) or not _DNA.fullmatch(template):
        raise ConfigurationError("template must be a string of A, C, G, T, N.")
    chosen = _gap_cli.choice(task, _TASKS, "task")
    seq_args: dict[str, Any] = {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": template}
    for key, value in (
        ("SEQUENCE_PRIMER", _oligo(left_primer, "left_primer")),
        ("SEQUENCE_PRIMER_REVCOMP", _oligo(right_primer, "right_primer")),
        ("SEQUENCE_INTERNAL_OLIGO", _oligo(internal_oligo, "internal_oligo")),
    ):
        if value is not None:
            seq_args[key] = value
    for key, value, name in (
        ("SEQUENCE_TARGET", target, "target"),
        ("SEQUENCE_INCLUDED_REGION", included_region, "included_region"),
    ):
        region = _region(value, name, len(template))
        if region is not None:
            seq_args[key] = region
    if excluded_regions is not None:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [
            _region(item, "excluded_regions", len(template)) for item in excluded_regions
        ]
    global_args: dict[str, Any] = {}
    for key, value in (settings or {}).items():
        if not _SETTING.fullmatch(key) or key == "PRIMER_TASK":
            raise ConfigurationError(
                f"settings key not allowed: {key!r} (use PRIMER_* tags; the task has its own parameter)."
            )
        global_args[key] = value
    global_args["PRIMER_TASK"] = chosen
    if num_return is not None:
        global_args["PRIMER_NUM_RETURN"] = _gap_cli.integer(
            num_return, "num_return", minimum=1, maximum=1000
        )
    if product_size_range is not None:
        global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [
            list(map(int, item)) for item in product_size_range
        ]
    try:
        raw = primer3.bindings.design_primers(seq_args, global_args)
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"primer3 rejected the input: {exc}") from exc
    result: dict[str, Any] = {
        "status": "completed",
        "backend": "primer3-py",
        "version": getattr(primer3, "__version__", None),
        "task": chosen,
    }
    groups: dict[str, dict[int, dict[str, Any]]] = {
        "LEFT": {},
        "RIGHT": {},
        "INTERNAL": {},
        "PAIR": {},
    }
    explain: dict[str, str] = {}
    for key, value in raw.items():
        match = _KEY.match(key)
        if match is not None:
            groups[match.group(1)].setdefault(int(match.group(2)), {})[match.group(3).lower()] = (
                value
            )
        elif key.endswith("_EXPLAIN"):
            explain[key] = value
    for name, rows in groups.items():
        items = []
        for index in sorted(rows):
            row = rows[index]
            row.pop("", None)
            items.append(row)
        result["pairs" if name == "PAIR" else name.lower()] = items
    for kind in ("left", "right", "internal"):
        for index, row in enumerate(result[kind]):
            raw_pos = raw.get(f"PRIMER_{kind.upper()}_{index}")
            if raw_pos is not None:
                row["start"], row["length"] = int(raw_pos[0]), int(raw_pos[1])
    result["explain"] = explain
    result["error"] = raw.get("PRIMER_ERROR") or None
    result["warning"] = raw.get("PRIMER_WARNING") or None
    if result["error"]:
        result["status"] = "failed"
    return result


def ntdpal_align_oligos(
    sequence_1: str,
    sequence_2: str,
    *,
    mode: str = "local",
    gap_open_penalty: float | None = None,
    gap_extend_penalty: float | None = None,
    max_gap: int | None = None,
    gc_weighted_scores: bool = False,
    ambiguity_codes: bool = False,
    primer3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align two short DNA sequences with Primer3's dynamic-programming aligner (``ntdpal``).

    ``mode`` is 'global' (mode g), 'global_end_anchored' (G, also written e; the alignment must reach the end of
    both sequences), 'local' (l) or 'local_end_anchored' (L). ``gap_open_penalty`` and ``gap_extend_penalty`` are
    subtracted from the score, ``max_gap`` is the longest allowed gap (default 3), ``gc_weighted_scores`` scores G/C
    matches 3, A/T 2 and mismatches -0.5 instead of +1/-1, and ``ambiguity_codes`` enables IUPAC codes. Returns the
    ``score``, ``length`` and the aligned position pairs ``pairs`` (index in sequence 1, index in sequence 2, 0-based);
    for local modes also the end positions ``end_1`` / ``end_2``.
    """
    codes = {"global": "g", "global_end_anchored": "G", "local": "l", "local_end_anchored": "L"}
    letter = codes[_gap_cli.choice(mode, tuple(codes), "mode")]
    first = _p14_cli.nucleotides(sequence_1, "sequence_1")
    second = _p14_cli.nucleotides(sequence_2, "sequence_2")
    arguments: list[str] = []
    if gap_open_penalty is not None:
        arguments += ["-g", str(_gap_cli.number(gap_open_penalty, "gap_open_penalty", minimum=0.0))]
    if gap_extend_penalty is not None:
        arguments += [
            "-l",
            str(_gap_cli.number(gap_extend_penalty, "gap_extend_penalty", minimum=0.0)),
        ]
    if max_gap is not None:
        arguments += ["-m", str(_gap_cli.integer(max_gap, "max_gap", minimum=0))]
    if _p14_cli.flag(gc_weighted_scores, "gc_weighted_scores"):
        arguments.append("-h")
    if _p14_cli.flag(ambiguity_codes, "ambiguity_codes"):
        arguments.append("-a")
    if letter in ("l", "L"):
        arguments.append("-e")
    arguments += [first, second, letter]
    with tempfile.TemporaryDirectory(prefix="dnakit_ntdpal_") as folder:
        work = Path(folder)
        result = _p14_cli.run(
            "ntdpal",
            arguments,
            cwd=work,
            expected=[],
            explicit=_exe("ntdpal", primer3_path),
            hint=_HINT,
            version_args=("-h",),
            stdout_file=work / "out.txt",
        )
        text = (work / "out.txt").read_text(encoding="utf-8", errors="replace").strip()
    result["outputs"] = []
    match = _ALIGN_LINE.search(text)
    if result["status"] == "completed" and match is not None:
        result["score"] = float(match.group(1))
        result["length"] = int(match.group(2))
        if match.group(3) is not None:
            result["end_1"], result["end_2"] = int(match.group(3)), int(match.group(4))
        pairs = (match.group(5) or "").split("|")
        result["pairs"] = [[int(a) for a in item.split(",")] for item in pairs if item]
        result["mode"] = mode
    elif result["status"] == "completed":
        result["status"] = "failed"
    return result


def _exe(name: str, primer3_path: str | None) -> str | None:
    if primer3_path is None:
        return None
    folder = Path(primer3_path).expanduser()
    if not folder.is_dir():
        raise ConfigurationError(
            f"primer3_path must be a folder that contains the Primer3 programs: {folder}"
        )
    return str(folder / name)


def primer3_mask_repeats(
    input_fasta: str,
    output: str,
    *,
    kmer_lists_path: str,
    list_prefix: str = "homo_sapiens",
    probability_cutoff: float | None = None,
    absolute_cutoff: int | None = None,
    kmer_list: str | None = None,
    mask_5_prime: int | None = None,
    mask_3_prime: int | None = None,
    masking_character: str | None = None,
    soft_mask: bool = False,
    masking_direction: str | None = None,
    primer3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask repeat k-mers in FASTA sequences so that primers avoid them (``primer3_masker``).

    By default each position is masked when the probability that its k-mer comes from a repeat reaches
    ``probability_cutoff`` (Primer3 default 0.1), using the list files ``<list_prefix>_11.list`` and
    ``<list_prefix>_16.list`` in ``kmer_lists_path`` (``homo_sapiens`` or ``mus_musculus``; the lists ship with the
    Primer3 source, folder ``kmer_lists``). With ``absolute_cutoff`` the count-based mode is used instead: a k-mer is
    masked when it occurs at least that many times in ``kmer_list`` (path of one list file). ``mask_5_prime`` and
    ``mask_3_prime`` give how many bases are masked beyond a flagged k-mer (defaults 1 and 0), ``masking_character``
    replaces them (default N), ``soft_mask`` writes lower case instead, and ``masking_direction`` is 'fwd', 'rev' or
    'both'. Writes the masked FASTA to a new file.
    """
    source = _gap_cli.existing_file(input_fasta, "input_fasta")
    output_path = _gap_cli.new_file(output, "output")
    lists = _p14_cli.existing_directory(kmer_lists_path, "kmer_lists_path")
    arguments = ["-lh", str(lists) + "/", "-lp", _gap_cli.text(list_prefix, "list_prefix")]
    if (absolute_cutoff is None) != (kmer_list is None):
        raise ConfigurationError("absolute_cutoff and kmer_list must be given together.")
    if absolute_cutoff is not None:
        if probability_cutoff is not None:
            raise ConfigurationError("probability_cutoff cannot be combined with absolute_cutoff.")
        arguments += [
            "-a",
            str(_gap_cli.integer(absolute_cutoff, "absolute_cutoff", minimum=1)),
            "-l",
            str(_gap_cli.existing_file(str(kmer_list), "kmer_list")),
        ]
    elif probability_cutoff is not None:
        arguments += [
            "-p",
            str(
                _gap_cli.number(probability_cutoff, "probability_cutoff", minimum=0.0, maximum=1.0)
            ),
        ]
    if mask_5_prime is not None:
        arguments += ["-m5", str(_gap_cli.integer(mask_5_prime, "mask_5_prime", minimum=0))]
    if mask_3_prime is not None:
        arguments += ["-m3", str(_gap_cli.integer(mask_3_prime, "mask_3_prime", minimum=0))]
    if masking_character is not None:
        if (
            not isinstance(masking_character, str)
            or len(masking_character) != 1
            or not masking_character.isalpha()
        ):
            raise ConfigurationError("masking_character must be one letter.")
        arguments += ["-c", masking_character]
    if _p14_cli.flag(soft_mask, "soft_mask"):
        arguments.append("-s")
    if masking_direction is not None:
        arguments += [
            "-d",
            _gap_cli.choice(masking_direction, ("fwd", "rev", "both"), "masking_direction"),
        ]
    arguments.append(str(source))
    result = _p14_cli.run(
        "primer3_masker",
        arguments,
        cwd=output_path.parent,
        expected=[],
        explicit=_exe("primer3_masker", primer3_path),
        hint=_HINT,
        version_args=("-h",),
        stdout_file=output_path,
    )
    result["mode"] = "absolute_count" if absolute_cutoff is not None else "probability"
    return result


__all__ = ["ntdpal_align_oligos", "primer3_design_by_task", "primer3_mask_repeats"]
_register_public_verbose(globals())
