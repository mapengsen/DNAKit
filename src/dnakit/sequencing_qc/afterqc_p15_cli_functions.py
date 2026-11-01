"""AfterQC command-line facades: paired-read overlap analysis, debubble and barcode handling.

AfterQC (https://github.com/OpenGene/AfterQC, tested against 0.9.7) is a Python 2 program: the user installs it (for example
``conda install -c bioconda afterqc``, which brings Python 2.7) and puts ``after.py`` on PATH or passes ``afterqc_path``. Each
function runs ``after.py`` with a fixed argument whitelist and no shell, in a new output folder holding the ``good``, ``bad``
and ``QC`` sub-folders that AfterQC writes. Results are JSON-serialisable dictionaries (status, backend, version, command,
output folders and files, log tail).

The plain filtering and trimming options (fixed or automatic end trimming, poly-X, unqualified-base and N-base limits) are
``dnakit.sequencing_qc.short_reads.run_afterqc``; this module adds the other AfterQC operations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p15_cli as _cli
from dnakit.exceptions import ConfigurationError

_TOOL = _cli.Tool(
    name="afterqc",
    executable="after.py",
    hint=(
        "Install AfterQC 0.9.7 (a Python 2 program; for example `conda install -c bioconda afterqc`) and put after.py on PATH or "
        "pass afterqc_path."
    ),
    tested_version="0.9.7",
    citation="https://doi.org/10.1186/s12859-017-1469-3",
    version_arguments=("--version",),
    version_pattern=r"(\d+\.\d+\.\d+)",
)


def _listing(folder: Path) -> list[str]:
    return sorted(str(path) for path in folder.rglob("*") if path.is_file()) if folder.is_dir() else []


def _quality_options(
    trim_front: int, trim_tail: int, qualified_quality_phred: int, unqualified_base_limit: int, poly_size_limit: int,
    n_base_limit: int, seq_len_req: int, gzip_output: bool,
) -> list[str]:
    arguments = [
        "-f", str(_cli.integer(trim_front, "trim_front", minimum=-1)),
        "-t", str(_cli.integer(trim_tail, "trim_tail", minimum=-1)),
        "-q", str(_cli.integer(qualified_quality_phred, "qualified_quality_phred", minimum=0, maximum=93)),
        "-u", str(_cli.integer(unqualified_base_limit, "unqualified_base_limit", minimum=0)),
        "-p", str(_cli.integer(poly_size_limit, "poly_size_limit", minimum=0)),
        "-n", str(_cli.integer(n_base_limit, "n_base_limit", minimum=0)),
        "-s", str(_cli.integer(seq_len_req, "seq_len_req", minimum=0)),
    ]
    if _cli.flag(gzip_output, "gzip_output"):
        arguments.append("--gzip")
    return arguments


def afterqc_overlap_analysis(
    read1: str,
    read2: str,
    output_directory: str,
    *,
    correct_mismatches: bool = True,
    mask_mismatch: bool = False,
    analyse_overlap: bool = True,
    store_overlap: bool = False,
    trim_front: int = -1,
    trim_tail: int = -1,
    qualified_quality_phred: int = 15,
    unqualified_base_limit: int = 60,
    poly_size_limit: int = 35,
    n_base_limit: int = 5,
    seq_len_req: int = 35,
    gzip_output: bool = False,
    afterqc_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim adapters, correct mismatched bases and store the overlap of read pairs by overlap analysis (AfterQC).

    AfterQC aligns the two mates of every pair on their overlap. Where the mates overlap it trims the adapter that sticks out
    (``analyse_overlap``, AfterQC's default; ``False`` is ``--no_overlap``), corrects bases that differ between the mates by
    keeping the one with the higher quality (``correct_mismatches``, default; ``False`` is ``--no_correction``), sets the
    quality of unresolved mismatching bases to 0 (``mask_mismatch``, ``--mask_mismatch``) and writes only the overlapping part
    of good pairs to the ``overlap`` folder (``store_overlap``, ``--store_overlap on``). The remaining options are AfterQC's
    filters and keep its defaults (``trim_front`` and ``trim_tail`` of -1 detect the trimming automatically). Output goes to the
    new folder ``output_directory`` with the sub-folders ``good``, ``bad``, ``QC`` and, when stored, ``overlap``; the result
    lists the files in ``good_files``, ``bad_files``, ``qc_files`` and ``overlap_files``.
    """
    first = _cli.existing_file(read1, "read1")
    second = _cli.existing_file(read2, "read2")
    if first == second:
        raise ConfigurationError("read1 and read2 must be different files.")
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = ["-1", str(first), "-2", str(second), "-g", str(out / "good"), "-b", str(out / "bad"), "-r", str(out / "QC")]
    arguments += _quality_options(
        trim_front, trim_tail, qualified_quality_phred, unqualified_base_limit, poly_size_limit, n_base_limit, seq_len_req,
        gzip_output,
    )
    if not _cli.flag(analyse_overlap, "analyse_overlap"):
        arguments.append("--no_overlap")
    if not _cli.flag(correct_mismatches, "correct_mismatches"):
        arguments.append("--no_correction")
    if _cli.flag(mask_mismatch, "mask_mismatch"):
        arguments.append("--mask_mismatch")
    if _cli.flag(store_overlap, "store_overlap"):
        arguments += ["--store_overlap", "on", "--overlap_output_folder", str(out / "overlap")]
    record = _cli.invoke(_TOOL, arguments, cwd=out, expected=[out / "good"], explicit=afterqc_path, timeout_seconds=timeout_seconds)
    return _cli.finish(
        record, "overlap_analysis_of_read_pairs", output_directory=str(out), good_files=_listing(out / "good"),
        bad_files=_listing(out / "bad"), qc_files=_listing(out / "QC"), overlap_files=_listing(out / "overlap"),
        analyse_overlap=analyse_overlap, correct_mismatches=correct_mismatches, store_overlap=store_overlap,
    )


def afterqc_debubble(
    input_directory: str,
    output_directory: str,
    *,
    draw: bool = True,
    trim_front: int = -1,
    trim_tail: int = -1,
    qualified_quality_phred: int = 15,
    unqualified_base_limit: int = 60,
    poly_size_limit: int = 35,
    n_base_limit: int = 5,
    seq_len_req: int = 35,
    gzip_output: bool = False,
    afterqc_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove reads from sequencer bubbles (air bubbles in the flow cell) before filtering (``after.py --debubble``).

    AfterQC maps the tile and position of every read of the Illumina FASTQ files in ``input_directory`` (files whose names
    contain ``R1``/``R2``), finds circular regions with an abnormal density of bad reads and moves the reads inside them to the
    ``bad`` folder. ``draw`` writes the bubble pictures (needs the Pillow imaging package in AfterQC's Python). The circles
    found are in ``debubble/circles.csv`` (returned as ``circles``); the cleaned reads are in ``good`` and the removed ones in
    ``bad``. The other options are AfterQC's filters with their defaults. All results go to the new folder ``output_directory``.
    """
    source = _cli.existing_directory(input_directory, "input_directory")
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = [
        "-d", str(source), "-g", str(out / "good"), "-b", str(out / "bad"), "-r", str(out / "QC"), "--debubble",
        "--debubble_dir", str(out / "debubble"), "--draw", "on" if _cli.flag(draw, "draw") else "off",
    ]
    arguments += _quality_options(
        trim_front, trim_tail, qualified_quality_phred, unqualified_base_limit, poly_size_limit, n_base_limit, seq_len_req,
        gzip_output,
    )
    record = _cli.invoke(_TOOL, arguments, cwd=out, expected=[out / "good"], explicit=afterqc_path, timeout_seconds=timeout_seconds)
    circles = out / "debubble" / "circles.csv"
    return _cli.finish(
        record, "remove_sequencer_bubble_reads", output_directory=str(out), good_files=_listing(out / "good"),
        bad_files=_listing(out / "bad"), debubble_files=_listing(out / "debubble"),
        circles=_cli.read_table(circles, delimiter=",", header=False) if circles.is_file() else [],
    )


def afterqc_move_barcodes_to_read_names(
    read1: str,
    output_directory: str,
    *,
    barcode_length: int = 12,
    barcode_verify: str = "CAGTA",
    barcode_flag: str = "barcode",
    read2: str | None = None,
    qualified_quality_phred: int = 15,
    unqualified_base_limit: int = 60,
    poly_size_limit: int = 35,
    n_base_limit: int = 5,
    seq_len_req: int = 35,
    gzip_output: bool = False,
    afterqc_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Move single-molecule barcodes from the start of reads into the read names (``after.py --barcode on``).

    For barcoded sequencing AfterQC reads the first ``barcode_length`` bases as the barcode, checks that the next bases equal
    ``barcode_verify`` and writes the barcode into the read name (and removes it from the sequence); reads whose verify sequence
    does not match go to ``bad``. AfterQC treats a file as barcoded when its name contains ``barcode_flag``, so ``read1`` (and
    ``read2``) must contain that word. Front trimming is switched off for barcoded files, as in AfterQC. The other options are
    AfterQC's filters. Output goes to the new folder ``output_directory`` (``good``, ``bad``, ``QC``).
    """
    first = _cli.existing_file(read1, "read1")
    flag = _cli.text(barcode_flag, "barcode_flag")
    if flag not in first.name:
        raise ConfigurationError(f"AfterQC treats a file as barcoded only when its name contains {flag!r}: {first.name}")
    out = _cli.new_directory(output_directory, "output_directory")
    arguments = ["-1", str(first)]
    if read2 is not None:
        arguments += ["-2", str(_cli.existing_file(read2, "read2"))]
    arguments += ["-g", str(out / "good"), "-b", str(out / "bad"), "-r", str(out / "QC"), "--barcode", "on"]
    arguments += ["--barcode_length", str(_cli.integer(barcode_length, "barcode_length", minimum=1))]
    arguments += ["--barcode_verify", _cli.text(barcode_verify, "barcode_verify"), "--barcode_flag", flag]
    arguments += _quality_options(
        -1, -1, qualified_quality_phred, unqualified_base_limit, poly_size_limit, n_base_limit, seq_len_req, gzip_output
    )
    record = _cli.invoke(_TOOL, arguments, cwd=out, expected=[out / "good"], explicit=afterqc_path, timeout_seconds=timeout_seconds)
    return _cli.finish(
        record, "move_barcodes_to_read_names", output_directory=str(out), good_files=_listing(out / "good"),
        bad_files=_listing(out / "bad"), qc_files=_listing(out / "QC"),
    )


__all__ = ["afterqc_debubble", "afterqc_move_barcodes_to_read_names", "afterqc_overlap_analysis"]
_register_public_verbose(globals())
