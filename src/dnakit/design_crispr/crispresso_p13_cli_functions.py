"""CRISPResso2 genome-editing analysis: command-line facades and ports of the repository helper scripts.

The user installs CRISPResso2 (``conda install -c bioconda crispresso2``, which also brings bowtie2, samtools and
fastp); 2.3.4 is the tested version. The ``crispresso_*`` facades run one CRISPResso2 console program
(``CRISPResso``, ``CRISPRessoBatch``, ``CRISPRessoPooled``, ``CRISPRessoWGS``, ``CRISPRessoCompare``,
``CRISPRessoPooledWGSCompare``, ``CRISPRessoAggregate``) with typed arguments mapped to fixed flags; no shell is used and
the program is resolved from an explicit path, then PATH. Options are passed only when given, so CRISPResso's own
defaults apply. Output folders are fresh directories. Results are JSON-serialisable dictionaries with ``status``,
``backend``, ``version``, ``command``, ``return_code``, ``outputs``, ``log_tail`` and parsed summaries.

The six ``crispresso_*`` helpers at the end port CRISPResso2's ``scripts/`` folder (which the package does not install).
Reading FASTQ files and counting bases need no extra package; the ones that read a finished CRISPResso run or draw
allele tables import the ``CRISPResso2`` Python package installed through Bioconda.

Version note: installation hints now point to the upstream Bioconda distribution.
"""

from __future__ import annotations

import gzip
import re
import zipfile
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _common
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p13_cli as _p13
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_HINT = (
    "Install CRISPResso2 with `conda install -c conda-forge -c bioconda crispresso2`, "
    "then put its programs on PATH or pass crispresso_path."
)
_SEQ = re.compile(r"^[ACGTNacgtn]+$")


def _seqs(value: str | Sequence[str] | None, name: str) -> str | None:
    if value is None:
        return None
    items = [value] if isinstance(value, str) else list(value)
    if not items:
        raise ConfigurationError(f"{name} must not be empty.")
    for item in items:
        if not isinstance(item, str) or not _SEQ.match(item):
            raise ConfigurationError(f"{name} must contain only A, C, G, T and N.")
    return ",".join(items)


def _names(value: str | Sequence[str] | None, name: str) -> str | None:
    if value is None:
        return None
    items = [value] if isinstance(value, str) else list(value)
    for item in items:
        _cli.text(item, name)
        if "," in item:
            raise ConfigurationError(f"{name} must not contain commas.")
    return ",".join(items)


def _opt(flag: str, value: object | None) -> list[str]:
    return [] if value is None else [flag, str(value)]


def _on(flag: str, value: bool, name: str) -> list[str]:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return [flag] if value else []


def _table(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        return []
    names = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        row: dict[str, Any] = {}
        for key, value in zip(names, line.split("\t"), strict=False):
            try:
                row[key] = int(value) if re.fullmatch(r"-?\d+", value) else float(value)
            except ValueError:
                row[key] = value
        rows.append(row)
    return rows


def _execute(
    program: str,
    args: Sequence[str],
    out: Path,
    crispresso_path: str | None,
    timeout_seconds: float,
    operation: str,
    *,
    cwd: Path | None = None,
) -> dict[str, Any]:
    result = _p13.run(
        program,
        args,
        cwd=cwd or out.parent,
        expected=[],
        explicit=crispresso_path,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = operation
    files = sorted(str(p) for p in out.rglob("*") if p.is_file()) if out.is_dir() else []
    result["outputs"] = files
    runs = sorted(p for p in out.glob("CRISPResso*_on_*") if p.is_dir()) if out.is_dir() else []
    result["run_folders"] = [str(p) for p in runs]
    if result["return_code"] == 0 and not runs:
        result["status"] = "failed"
    return result


def _start(program: str, crispresso_path: str | None, out_dir: str) -> Path:
    _cli.program(program, explicit=crispresso_path, hint=_HINT)
    out = Path(out_dir).expanduser().resolve()
    if out.exists():
        raise ConfigurationError(f"output_folder must not exist yet: {out}")
    if not out.parent.is_dir():
        raise ConfigurationError(f"The folder for output_folder does not exist: {out.parent}")
    return out


# ---------------------------------------------------------------------------
# CRISPResso: one amplicon set, every editing mode
# ---------------------------------------------------------------------------


def crispresso_quantify_editing(
    fastq_r1_path: str,
    output_folder: str,
    *,
    amplicon_sequences: str | Sequence[str],
    fastq_r2_path: str | None = None,
    amplicon_names: str | Sequence[str] | None = None,
    guide_sequences: str | Sequence[str] | None = None,
    guide_names: str | Sequence[str] | None = None,
    expected_hdr_amplicon_sequence: str | None = None,
    coding_sequences: str | Sequence[str] | None = None,
    base_editor_output: bool = False,
    conversion_nuc_from: str | None = None,
    conversion_nuc_to: str | None = None,
    prime_editing_pegrna_spacer_seq: str | None = None,
    prime_editing_pegrna_extension_seq: str | None = None,
    prime_editing_pegrna_scaffold_seq: str | None = None,
    prime_editing_nicking_guide_seq: str | None = None,
    prime_editing_extension_quantification_window_size: int | None = None,
    prime_editing_scaffold_min_match_length: int | None = None,
    name: str | None = None,
    quantification_window_size: int | None = None,
    quantification_window_center: int | None = None,
    quantification_window_coordinates: str | None = None,
    exclude_bp_from_left: int | None = None,
    exclude_bp_from_right: int | None = None,
    min_average_read_quality: int | None = None,
    min_single_bp_quality: int | None = None,
    min_bp_quality_or_n: int | None = None,
    min_identity_score: float | None = None,
    ignore_substitutions: bool = False,
    ignore_insertions: bool = False,
    ignore_deletions: bool = False,
    discard_indel_reads: bool = False,
    trim_sequences: bool = False,
    write_detailed_allele_table: bool = False,
    fastq_output: bool = False,
    bam_output: bool = False,
    vcf_output: bool = False,
    suppress_plots: bool = False,
    suppress_report: bool = False,
    keep_intermediate: bool = False,
    n_processes: int | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Quantify genome-editing outcomes at one or more amplicons from sequencing reads (``CRISPResso``).

    One function covers every CRISPResso editing mode: plain indel quantification (``amplicon_sequences`` with optional
    ``guide_sequences``); HDR outcomes against a donor amplicon (``expected_hdr_amplicon_sequence``, ``-e``);
    frameshift analysis of edits in coding sequences (``coding_sequences``, ``-c``); base-editing conversion analysis
    (``base_editor_output=True`` with ``conversion_nuc_from``/``conversion_nuc_to``, for example ``"C"``/``"T"``); and
    prime-editing analysis (``prime_editing_pegrna_spacer_seq`` plus extension, scaffold and nicking-guide sequences and
    the two window/match-length parameters). Sequences may be one string or a list (joined with commas as CRISPResso
    expects). ``output_folder`` is created and holds ``CRISPResso_on_<name>``; ``quantification`` returns its
    ``CRISPResso_quantification_of_editing_frequency.txt`` rows and ``mapping_statistics`` the mapping counts.
    The remaining arguments map one to one to CRISPResso options (``-w``, ``-wc``, ``-qwc``, ``-q``, ``-s``,
    ``--min_bp_quality_or_N``, ``--default_min_aln_score``, ``--ignore_*``, ``--discard_indel_reads``,
    ``--trim_sequences``, ``--write_detailed_allele_table``, ``--fastq_output``, ``--bam_output``, ``--vcf_output``,
    ``--suppress_plots``, ``--suppress_report``, ``--keep_intermediate``, ``-p``).
    """
    r1 = _cli.existing_file(fastq_r1_path, "fastq_r1_path")
    args = ["-r1", str(r1)]
    if fastq_r2_path is not None:
        args += ["-r2", str(_cli.existing_file(fastq_r2_path, "fastq_r2_path"))]
    args += _opt("-a", _seqs(amplicon_sequences, "amplicon_sequences"))
    args += _opt("-an", _names(amplicon_names, "amplicon_names"))
    args += _opt("-g", _seqs(guide_sequences, "guide_sequences"))
    args += _opt("-gn", _names(guide_names, "guide_names"))
    args += _opt("-e", _seqs(expected_hdr_amplicon_sequence, "expected_hdr_amplicon_sequence"))
    args += _opt("-c", _seqs(coding_sequences, "coding_sequences"))
    if base_editor_output:
        if conversion_nuc_from is None or conversion_nuc_to is None:
            raise ConfigurationError("base_editor_output needs conversion_nuc_from and conversion_nuc_to.")
    elif conversion_nuc_from is not None or conversion_nuc_to is not None:
        raise ConfigurationError("conversion_nuc_from/to only apply with base_editor_output=True.")
    args += _on("--base_editor_output", base_editor_output, "base_editor_output")
    for flag, value in (("--conversion_nuc_from", conversion_nuc_from), ("--conversion_nuc_to", conversion_nuc_to)):
        args += _opt(flag, None if value is None else _cli.choice(value, ("A", "C", "G", "T"), flag))
    args += _opt("--prime_editing_pegRNA_spacer_seq", _seqs(prime_editing_pegrna_spacer_seq, "pegrna spacer"))
    args += _opt("--prime_editing_pegRNA_extension_seq", _seqs(prime_editing_pegrna_extension_seq, "pegrna extension"))
    args += _opt("--prime_editing_pegRNA_scaffold_seq", _seqs(prime_editing_pegrna_scaffold_seq, "pegrna scaffold"))
    args += _opt("--prime_editing_nicking_guide_seq", _seqs(prime_editing_nicking_guide_seq, "nicking guide"))
    for flag, value, label in (
        ("--prime_editing_pegRNA_extension_quantification_window_size", prime_editing_extension_quantification_window_size, "window"),
        ("--prime_editing_pegRNA_scaffold_min_match_length", prime_editing_scaffold_min_match_length, "match length"),
        ("-w", quantification_window_size, "quantification_window_size"),
        ("-wc", quantification_window_center, "quantification_window_center"),
        ("--exclude_bp_from_left", exclude_bp_from_left, "exclude_bp_from_left"),
        ("--exclude_bp_from_right", exclude_bp_from_right, "exclude_bp_from_right"),
        ("-q", min_average_read_quality, "min_average_read_quality"),
        ("-s", min_single_bp_quality, "min_single_bp_quality"),
        ("--min_bp_quality_or_N", min_bp_quality_or_n, "min_bp_quality_or_n"),
    ):
        args += _opt(flag, None if value is None else _cli.integer(value, label, minimum=-1 if flag == "-wc" else 0))
    if quantification_window_coordinates is not None:
        args += ["-qwc", _cli.text(quantification_window_coordinates, "quantification_window_coordinates")]
    if min_identity_score is not None:
        args += ["--default_min_aln_score", f"{_cli.number(min_identity_score, 'min_identity_score', minimum=0.0, maximum=100.0):g}"]
    if name is not None:
        args += ["-n", _cli.text(name, "name")]
    for flag, value, label in (
        ("--ignore_substitutions", ignore_substitutions, "ignore_substitutions"),
        ("--ignore_insertions", ignore_insertions, "ignore_insertions"),
        ("--ignore_deletions", ignore_deletions, "ignore_deletions"),
        ("--discard_indel_reads", discard_indel_reads, "discard_indel_reads"),
        ("--trim_sequences", trim_sequences, "trim_sequences"),
        ("--write_detailed_allele_table", write_detailed_allele_table, "write_detailed_allele_table"),
        ("--fastq_output", fastq_output, "fastq_output"),
        ("--bam_output", bam_output, "bam_output"),
        ("--vcf_output", vcf_output, "vcf_output"),
        ("--suppress_plots", suppress_plots, "suppress_plots"),
        ("--suppress_report", suppress_report, "suppress_report"),
        ("--keep_intermediate", keep_intermediate, "keep_intermediate"),
    ):
        args += _on(flag, value, label)
    args += _opt("-p", None if n_processes is None else _cli.integer(n_processes, "n_processes", minimum=1, maximum=256))
    out = _start("CRISPResso", crispresso_path, output_folder)
    out.mkdir()
    result = _execute("CRISPResso", [*args, "-o", str(out)], out, crispresso_path, timeout_seconds, "CRISPResso")
    run = Path(result["run_folders"][0]) if result["run_folders"] else None
    result["quantification"] = _table(run / "CRISPResso_quantification_of_editing_frequency.txt") if run else []
    stats = _table(run / "CRISPResso_mapping_statistics.txt") if run else []
    result["mapping_statistics"] = stats[0] if stats else {}
    return result


def crispresso_batch(
    batch_settings_path: str,
    output_folder: str,
    *,
    batch_name: str | None = None,
    amplicon_sequences: str | Sequence[str] | None = None,
    guide_sequences: str | Sequence[str] | None = None,
    quantification_window_size: int | None = None,
    suppress_plots: bool = False,
    suppress_batch_summary_plots: bool = False,
    n_processes: int | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run CRISPResso on many samples listed in a tab-separated settings file (``CRISPRessoBatch``).

    The header row of ``batch_settings_path`` holds CRISPResso parameter names (``name``, ``fastq_r1``, ``fastq_r2``,
    ``amplicon_seq``, ``guide_seq`` ...), one sample per row; ``amplicon_sequences``/``guide_sequences`` set defaults
    for samples that leave them empty. The per-sample folders and the batch summary are written to ``output_folder``.
    """
    settings = _cli.existing_file(batch_settings_path, "batch_settings_path")
    args = ["-bs", str(settings)]
    args += _opt("-a", _seqs(amplicon_sequences, "amplicon_sequences"))
    args += _opt("-g", _seqs(guide_sequences, "guide_sequences"))
    args += _opt("-w", None if quantification_window_size is None else _cli.integer(quantification_window_size, "quantification_window_size", minimum=0))
    args += _opt("-n", None if batch_name is None else _cli.text(batch_name, "batch_name"))
    args += _on("--suppress_plots", suppress_plots, "suppress_plots") + _on("--suppress_batch_summary_plots", suppress_batch_summary_plots, "suppress_batch_summary_plots")
    args += _opt("-p", None if n_processes is None else _cli.integer(n_processes, "n_processes", minimum=1, maximum=256))
    out = _start("CRISPRessoBatch", crispresso_path, output_folder)
    out.mkdir()
    return _execute("CRISPRessoBatch", [*args, "-o", str(out)], out, crispresso_path, timeout_seconds, "CRISPRessoBatch")


def crispresso_pooled(
    fastq_r1_path: str,
    output_folder: str,
    *,
    amplicons_file_path: str | None = None,
    genome_index: str | None = None,
    gene_annotations_path: str | None = None,
    fastq_r2_path: str | None = None,
    name: str | None = None,
    min_reads_to_use_region: int | None = None,
    demultiplex_only_at_amplicons: bool = False,
    guide_sequences: str | Sequence[str] | None = None,
    quantification_window_size: int | None = None,
    suppress_plots: bool = False,
    n_processes: int | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Demultiplex and quantify pooled amplicon sequencing (``CRISPRessoPooled``).

    Give ``amplicons_file_path`` (tab-delimited ``name``, ``amplicon_seq``, optional guide and window columns), a
    bowtie2 ``genome_index`` (``-x``; reads are mapped and regions with enough reads are analysed, optionally
    annotated with ``gene_annotations_path``), or both. ``min_reads_to_use_region`` is CRISPResso's threshold
    (default 1000) and ``demultiplex_only_at_amplicons`` restricts demultiplexing to the listed amplicons.
    """
    if amplicons_file_path is None and genome_index is None:
        raise ConfigurationError("Give amplicons_file_path, genome_index or both.")
    args = ["-r1", str(_cli.existing_file(fastq_r1_path, "fastq_r1_path"))]
    if fastq_r2_path is not None:
        args += ["-r2", str(_cli.existing_file(fastq_r2_path, "fastq_r2_path"))]
    if amplicons_file_path is not None:
        args += ["-f", str(_cli.existing_file(amplicons_file_path, "amplicons_file_path"))]
    args += _opt("-x", None if genome_index is None else _cli.text(genome_index, "genome_index"))
    if gene_annotations_path is not None:
        args += ["--gene_annotations", str(_cli.existing_file(gene_annotations_path, "gene_annotations_path"))]
    args += _opt("-n", None if name is None else _cli.text(name, "name"))
    args += _opt("--min_reads_to_use_region", None if min_reads_to_use_region is None else _cli.integer(min_reads_to_use_region, "min_reads_to_use_region", minimum=1))
    args += _on("--demultiplex_only_at_amplicons", demultiplex_only_at_amplicons, "demultiplex_only_at_amplicons")
    args += _opt("-g", _seqs(guide_sequences, "guide_sequences"))
    args += _opt("-w", None if quantification_window_size is None else _cli.integer(quantification_window_size, "quantification_window_size", minimum=0))
    args += _on("--suppress_plots", suppress_plots, "suppress_plots")
    args += _opt("-p", None if n_processes is None else _cli.integer(n_processes, "n_processes", minimum=1, maximum=256))
    out = _start("CRISPRessoPooled", crispresso_path, output_folder)
    out.mkdir()
    return _execute("CRISPRessoPooled", [*args, "-o", str(out)], out, crispresso_path, timeout_seconds, "CRISPRessoPooled")


def crispresso_wgs(
    bam_path: str,
    region_file_path: str,
    reference_fasta_path: str,
    output_folder: str,
    *,
    name: str | None = None,
    gene_annotations_path: str | None = None,
    min_reads_to_use_region: int | None = None,
    guide_sequences: str | Sequence[str] | None = None,
    quantification_window_size: int | None = None,
    suppress_plots: bool = False,
    n_processes: int | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Quantify editing at BED regions of a whole-genome BAM (``CRISPRessoWGS``).

    ``bam_path`` is the sorted, indexed WGS BAM, ``region_file_path`` a BED file of the regions to analyse and
    ``reference_fasta_path`` the reference the reads were aligned to.
    """
    args = ["-b", str(_cli.existing_file(bam_path, "bam_path")), "-f", str(_cli.existing_file(region_file_path, "region_file_path"))]
    args += ["-r", str(_cli.existing_file(reference_fasta_path, "reference_fasta_path"))]
    args += _opt("-n", None if name is None else _cli.text(name, "name"))
    if gene_annotations_path is not None:
        args += ["--gene_annotations", str(_cli.existing_file(gene_annotations_path, "gene_annotations_path"))]
    args += _opt("--min_reads_to_use_region", None if min_reads_to_use_region is None else _cli.integer(min_reads_to_use_region, "min_reads_to_use_region", minimum=1))
    args += _opt("-g", _seqs(guide_sequences, "guide_sequences"))
    args += _opt("-w", None if quantification_window_size is None else _cli.integer(quantification_window_size, "quantification_window_size", minimum=0))
    args += _on("--suppress_plots", suppress_plots, "suppress_plots")
    args += _opt("-p", None if n_processes is None else _cli.integer(n_processes, "n_processes", minimum=1, maximum=256))
    out = _start("CRISPRessoWGS", crispresso_path, output_folder)
    out.mkdir()
    return _execute("CRISPRessoWGS", [*args, "-o", str(out)], out, crispresso_path, timeout_seconds, "CRISPRessoWGS")


def _compare(
    program: str,
    folder1: str,
    folder2: str,
    output_folder: str,
    sample_names: Sequence[str] | None,
    name: str | None,
    qvalue: float | None,
    crispresso_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    first = Path(folder1).expanduser().resolve()
    second = Path(folder2).expanduser().resolve()
    for folder in (first, second):
        if not folder.is_dir():
            raise ConfigurationError(f"Not a CRISPResso output folder: {folder}")
    args = [str(first), str(second)]
    if sample_names is not None:
        if isinstance(sample_names, str) or len(sample_names) != 2:
            raise ConfigurationError("sample_names must hold exactly two names.")
        args += ["-n1", _cli.text(sample_names[0], "sample_names"), "-n2", _cli.text(sample_names[1], "sample_names")]
    args += _opt("-n", None if name is None else _cli.text(name, "name"))
    args += _opt("--reported_qvalue_cutoff", None if qvalue is None else f"{_cli.number(qvalue, 'reported_qvalue_cutoff', minimum=0.0, maximum=1.0):g}")
    out = _start(program, crispresso_path, output_folder)
    out.mkdir()
    result = _p13.run(
        program, [*args, "-o", str(out)], cwd=out, expected=[], explicit=crispresso_path, hint=_HINT, timeout_seconds=timeout_seconds
    )
    result["operation"] = program
    result["outputs"] = sorted(str(p) for p in out.rglob("*") if p.is_file())
    result["run_folders"] = [str(p) for p in sorted(out.glob("*_on_*")) if p.is_dir()]
    if result["return_code"] == 0 and not result["run_folders"]:
        result["status"] = "failed"
    return result


def crispresso_compare(
    folder1: str,
    folder2: str,
    output_folder: str,
    *,
    sample_names: Sequence[str] | None = None,
    name: str | None = None,
    reported_qvalue_cutoff: float | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the editing outcomes of two finished CRISPResso analyses (``CRISPRessoCompare``).

    ``folder1`` and ``folder2`` are ``CRISPResso_on_<name>`` folders of runs on the same amplicon; the comparison
    plots and tables (with q-values below ``reported_qvalue_cutoff``) go to ``output_folder``. The two analyses must
    have been run with plots and the report enabled.
    """
    return _compare("CRISPRessoCompare", folder1, folder2, output_folder, sample_names, name, reported_qvalue_cutoff, crispresso_path, timeout_seconds)


def crispresso_pooled_wgs_compare(
    folder1: str,
    folder2: str,
    output_folder: str,
    *,
    sample_names: Sequence[str] | None = None,
    name: str | None = None,
    reported_qvalue_cutoff: float | None = None,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two CRISPRessoPooled or CRISPRessoWGS analyses region by region (``CRISPRessoPooledWGSCompare``)."""
    return _compare("CRISPRessoPooledWGSCompare", folder1, folder2, output_folder, sample_names, name, reported_qvalue_cutoff, crispresso_path, timeout_seconds)


def crispresso_aggregate(
    search_directory: str,
    name: str,
    *,
    folder_prefixes: Sequence[str] | None = None,
    folder_suffix: str | None = None,
    min_reads_for_inclusion: int | None = None,
    suppress_plots: bool = False,
    crispresso_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Aggregate the editing summaries of many CRISPResso runs into one report (``CRISPRessoAggregate``).

    CRISPResso folders in ``search_directory`` whose names start with one of ``folder_prefixes`` and/or end with
    ``folder_suffix`` are combined; the summary ``CRISPRessoAggregate_on_<name>`` is written into the same directory.
    """
    directory = Path(search_directory).expanduser().resolve()
    if not directory.is_dir():
        raise ConfigurationError(f"search_directory does not exist: {directory}")
    if not folder_prefixes and folder_suffix is None:
        raise ConfigurationError("Give folder_prefixes, folder_suffix or both.")
    args = ["-n", _cli.text(name, "name")]
    for prefix in folder_prefixes or ():
        args += ["-p", _cli.text(prefix, "folder_prefixes")]
    args += _opt("-s", None if folder_suffix is None else _cli.text(folder_suffix, "folder_suffix"))
    args += _opt("--min_reads_for_inclusion", None if min_reads_for_inclusion is None else _cli.integer(min_reads_for_inclusion, "min_reads_for_inclusion", minimum=0))
    args += _on("--suppress_plots", suppress_plots, "suppress_plots")
    _cli.program("CRISPRessoAggregate", explicit=crispresso_path, hint=_HINT)
    target = directory / f"CRISPRessoAggregate_on_{name}"
    if target.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {target}")
    result = _p13.run(
        "CRISPRessoAggregate", args, cwd=directory, expected=[], explicit=crispresso_path, hint=_HINT, timeout_seconds=timeout_seconds
    )
    result["operation"] = "CRISPRessoAggregate"
    result["outputs"] = sorted(str(p) for p in target.rglob("*") if p.is_file()) if target.is_dir() else []
    if result["return_code"] == 0 and not target.is_dir():
        result["status"] = "failed"
    return result


# ---------------------------------------------------------------------------
# Ports of CRISPResso2 scripts/
# ---------------------------------------------------------------------------


def _open_text(path: Path, mode: str) -> Any:
    return gzip.open(path, mode + "t") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def crispresso_filter_reads_by_sequence(
    fastq_r1_path: str,
    output_r1_path: str,
    *,
    fastq_r2_path: str | None = None,
    output_r2_path: str | None = None,
    include_sequences: Sequence[str] = (),
    exclude_sequences: Sequence[str] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep FASTQ reads that contain (or do not contain) given sequences (port of ``filterReadsOnSequencePresence.py``).

    A read is kept when it contains at least one of ``include_sequences`` (all reads pass when the list is empty) and
    none of ``exclude_sequences``; with ``fastq_r2_path`` a pair is judged on both mates and written to both outputs.
    ``.gz`` inputs and outputs are supported. The original single-end writer omits the newline after each quality line,
    which glues records together; this port writes a proper FASTQ. Returns the read counts and the per-sequence counts.
    """
    inc = [_cli.text(s, "include_sequences") for s in include_sequences]
    exc = [_cli.text(s, "exclude_sequences") for s in exclude_sequences]
    if not inc and not exc:
        raise ConfigurationError("Give include_sequences or exclude_sequences.")
    r1 = _cli.existing_file(fastq_r1_path, "fastq_r1_path")
    r2 = _cli.existing_file(fastq_r2_path, "fastq_r2_path") if fastq_r2_path else None
    if (r2 is None) != (output_r2_path is None):
        raise ConfigurationError("fastq_r2_path and output_r2_path must be given together.")
    o1 = _cli.new_file(output_r1_path, "output_r1_path")
    o2 = _cli.new_file(output_r2_path, "output_r2_path") if output_r2_path else None
    total = written = 0
    inc_count: dict[str, int] = defaultdict(int)
    exc_count: dict[str, int] = defaultdict(int)
    handles = [_open_text(r1, "r"), _open_text(o1, "w")]
    if r2 is not None and o2 is not None:
        handles += [_open_text(r2, "r"), _open_text(o2, "w")]
    try:
        while True:
            rec1 = [handles[0].readline().rstrip() for _ in range(4)]
            rec2 = [handles[2].readline().rstrip() for _ in range(4)] if r2 is not None else None
            if not rec1[0]:
                break
            total += 1
            seqs = [rec1[1]] + ([rec2[1]] if rec2 else [])
            keep = not inc
            for s in inc:
                if any(s in x for x in seqs):
                    keep = True
                    inc_count[s] += 1
            drop = False
            for s in exc:
                if any(s in x for x in seqs):
                    drop = True
                    exc_count[s] += 1
            if keep and not drop:
                written += 1
                handles[1].write("\n".join(rec1) + "\n")
                if rec2:
                    handles[3].write("\n".join(rec2) + "\n")
    finally:
        for handle in handles:
            handle.close()
    return {
        "status": "completed",
        "operation": "filterReadsOnSequencePresence",
        "outputs": [str(p) for p in (o1, o2) if p is not None],
        "reads_read": total,
        "reads_written": written,
        "reads_containing_included": {s: inc_count[s] for s in inc},
        "reads_containing_excluded": {s: exc_count[s] for s in exc},
    }


def _plus_dict(line: str) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for element in line.rstrip("\n").split(" "):
        parts = element.split("=")
        out[parts[0]] = parts[1] if len(parts) > 1 else None
    return out


def crispresso_count_high_quality_bases(
    annotated_fastq_path: str,
    position: int,
    original_base: str,
    target_base: str,
    *,
    reference_name: str = "Reference",
    quality_cutoff: int = 30,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count high- and low-quality bases at one reference position in a CRISPResso-annotated FASTQ (port of ``countHighQualityBases.py``).

    ``annotated_fastq_path`` is the ``CRISPResso_output.fastq.gz`` written with ``--fastq_output``; ``position`` is the
    0-based reference position. Reads whose ``ALN`` is another reference are counted as ``unaligned``. A base counts as
    high quality when its Phred score (ASCII - 33) is at least ``quality_cutoff`` (default 30).
    """
    path = _cli.existing_file(annotated_fastq_path, "annotated_fastq_path")
    pos = _cli.integer(position, "position", minimum=0)
    cutoff = _cli.integer(quality_cutoff, "quality_cutoff", minimum=0)
    counts = dict.fromkeys(
        ("total_reads", "unaligned", "original_high", "original_low", "target_high", "target_low", "deletion", "other_high", "other_low"), 0
    )
    with _open_text(path, "r") as handle:
        while True:
            header, _seq, plus, qual = (handle.readline() for _ in range(4))
            if not header:
                break
            counts["total_reads"] += 1
            info = _plus_dict(plus)
            if "ALN" not in info:
                raise ConfigurationError(f"ALN not in plus line: {plus.strip()}")
            if info["ALN"] != reference_name:
                counts["unaligned"] += 1
                continue
            details = {}
            for element in (info.get("ALN_DETAILS") or "").split("&"):
                fields = element.split(",")
                if len(fields) != 4:
                    raise ConfigurationError(f"unexpected ALN_DETAILS: {fields}")
                details[fields[0]] = fields
            read_aln, ref_aln = details[reference_name][1], details[reference_name][2]
            read_pos = -1
            in_aln: list[int] = []
            in_read: list[int] = []
            for i, (r, f) in enumerate(zip(read_aln, ref_aln, strict=False)):
                if r != "-":
                    read_pos += 1
                if f != "-":
                    in_aln.append(i)
                    in_read.append(read_pos)
            base = read_aln[in_aln[pos]]
            high = ord(qual[in_read[pos]]) - 33 >= cutoff
            if base == original_base:
                counts["original_high" if high else "original_low"] += 1
            elif base == target_base:
                counts["target_high" if high else "target_low"] += 1
            elif base == "-":
                counts["deletion"] += 1
            else:
                counts["other_high" if high else "other_low"] += 1
    return {"status": "completed", "operation": "countHighQualityBases", "counts": counts}


def _require_crispresso(module: str) -> Any:
    try:
        return _common.require(module, "crispresso2")
    except BackendUnavailableError as exc:
        raise BackendUnavailableError(exc.message, hint=_HINT) from exc


def _load_run(folder: str) -> tuple[Any, Path, Any, Any]:
    shared = _require_crispresso("CRISPResso2.CRISPRessoShared")
    pandas = _common.require("pandas", "bio-analysis")
    directory = Path(folder).expanduser().resolve()
    if not directory.is_dir():
        raise ConfigurationError(f"Not a CRISPResso output folder: {directory}")
    info = shared.load_crispresso_info(str(directory))
    return info, directory, pandas, shared


def _alleles(info: dict[str, Any], directory: Path, pandas: Any) -> Any:
    running = info["running_info"]
    with zipfile.ZipFile(directory / running["allele_frequency_table_zip_filename"]) as archive, archive.open(
        running["allele_frequency_table_filename"]
    ) as handle:
        return pandas.read_csv(handle, sep="\t")


def crispresso_count_sgrna_specific_edits(crispresso_folder: str, *, verbose: bool = True) -> dict[str, Any]:
    """Count the edits that fall in each sgRNA's window in a finished CRISPResso run (port of ``count_sgRNA_specific_edits.py``).

    The run must have used ``--write_detailed_allele_table`` and CRISPResso2 2.2.15 or later. Returns, per reference,
    the number of modified and total reads and the read counts per category (``UNMODIFIED`` or ``MODIFIED <guide> + ...``).
    """
    info, directory, pandas, _ = _load_run(crispresso_folder)
    parts = info["running_info"]["version"].split(".")
    if (int(parts[0]), int(parts[1]), int(parts[2])) < (2, 2, 15):
        raise ConfigurationError("The CRISPResso run needs CRISPResso2 2.2.15 or later.")
    if not info["running_info"]["args"].write_detailed_allele_table:
        raise ConfigurationError("The CRISPResso run must use --write_detailed_allele_table.")
    refs = info["results"]["refs"]
    include: dict[str, dict[str, Any]] = {}
    order: dict[str, list[str]] = {}
    for ref in info["results"]["ref_names"]:
        include[ref], order[ref] = {}, []
        for idx, guide in enumerate(refs[ref]["sgRNA_names"]):
            label = guide or refs[ref]["sgRNA_orig_sequences"][idx]
            if label in include[ref]:
                label = f"{label}_{refs[ref]['sgRNA_intervals'][idx][0]}"
            include[ref][label] = refs[ref]["sgRNA_include_idxs"][idx]
            order[ref].append(label)
    df = _alleles(info, directory, pandas)
    total: dict[str, int] = defaultdict(int)
    modified: dict[str, int] = defaultdict(int)
    categories: dict[str, dict[str, int]] = {ref: defaultdict(int) for ref in include}

    def positions(text: str) -> list[int]:
        return [int(x) for x in text[1:-1].split(",")] if text != "[]" else []

    for _, row in df.iterrows():
        ref = row["Reference_Name"]
        if ref not in include:
            continue
        reads = int(row["#Reads"])
        total[ref] += reads
        hit = set(positions(row["insertion_positions"])) | set(positions(row["deletion_positions"])) | set(positions(row["substitution_positions"]))
        guides = [g for g in order[ref] if any(i in hit for i in include[ref][g])]
        if guides:
            modified[ref] += reads
            categories[ref]["MODIFIED " + " + ".join(guides)] += reads
        else:
            categories[ref]["UNMODIFIED"] += reads
    return {
        "status": "completed",
        "operation": "count_sgRNA_specific_edits",
        "references": {
            ref: {"modified_reads": modified[ref], "total_reads": total[ref], "categories": dict(categories[ref])} for ref in include
        },
    }


def crispresso_reconstitute_reads(crispresso_folder: str, output_fastq_path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Regenerate a FASTQ from the allele frequency table of a finished CRISPResso run (port of ``reconstituteReads.py``).

    Each allele is written once per read (gaps removed, names ``SEQ_ID_<allele>_<copy>``, quality ``A``) so the run can be
    re-analysed. ``reads_written`` is compared with the number of aligned reads of the original run.
    """
    info, directory, pandas, _ = _load_run(crispresso_folder)
    out = _cli.new_file(output_fastq_path, "output_fastq_path")
    df = _alleles(info, directory, pandas)
    if df.shape[0] == 0:
        raise ConfigurationError("No alleles found in allele frequency table.")
    for column in ("#Reads", "Aligned_Sequence", "Reference_Sequence"):
        if column not in df.columns:
            raise ConfigurationError(f"Cannot parse allele frequency table, missing {column} column.")
    reads = 0
    with out.open("w", encoding="utf-8") as handle:
        for allele, (_, row) in enumerate(df.iterrows()):
            seq = row["Aligned_Sequence"].replace("-", "")
            for copy in range(int(row["#Reads"])):
                handle.write(f"@SEQ_ID_{allele}_{copy}\n{seq}\n+\n{'A' * len(seq)}\n")
                reads += 1
    stats = info["running_info"]["alignment_stats"]
    expected = stats["N_CACHED_ALN"] + stats["N_COMPUTED_ALN"]
    return {
        "status": "completed",
        "operation": "reconstituteReads",
        "outputs": [str(out)],
        "alleles_written": int(df.shape[0]),
        "reads_written": reads,
        "reads_aligned_in_original_run": int(expected),
        "matches_original": reads == expected,
    }


def _plot_alleles(
    crispresso_folder: str,
    output_root: str,
    mode: str,
    min_frequency: float | None,
    max_rows: int | None,
    plot_cut_point: bool,
    save_png: bool,
    plot_left: int,
    plot_right: int,
    plot_center: int | None,
) -> dict[str, Any]:
    info, directory, pandas, shared = _load_run(crispresso_folder)
    plot = _require_crispresso("CRISPResso2.plots.CRISPRessoPlot")
    args = info["running_info"]["args"]
    if not args.write_detailed_allele_table:
        raise ConfigurationError("The CRISPResso run must use --write_detailed_allele_table.")
    min_freq = args.min_frequency_alleles_around_cut_to_plot if min_frequency is None else min_frequency
    rows = args.max_rows_alleles_around_cut_to_plot if max_rows is None else max_rows
    root = Path(output_root).expanduser().resolve()
    if not root.parent.is_dir():
        raise ConfigurationError(f"The folder for output_root does not exist: {root.parent}")
    df = _alleles(info, directory, pandas)
    df["ref_positions"] = df["ref_positions"].apply(lambda v: [int(x) for x in v[1:-1].split(",")])
    if mode == "ambiguous":
        df = df[df["Reference_Name"].str.contains("AMBIGUOUS")]
    produced: list[str] = []
    refs = info["results"]["refs"]

    def window(frame: Any, cut: int, left: int, right: int) -> Any:
        if frame.shape[0] == 0:
            return frame
        length = len(frame["Reference_Sequence"].iloc[0].replace("-", ""))
        if cut + right + 1 > length:
            raise ConfigurationError(f"The plotting window cannot extend past the end of the amplicon ({length} bp).")
        if mode == "ambiguous":
            return shared.get_dataframe_around_cut_asymmetrical(frame, cut, left, right)
        rowsout = []
        for _, r in frame.iterrows():
            idx = r["ref_positions"].index(cut)
            sl = slice(idx - left + 1, idx + right + 1)
            rowsout.append((r["Aligned_Sequence"][sl], r["Reference_Sequence"][sl], r["Read_Status"] == "UNMODIFIED", r["n_deleted"], r["n_inserted"], r["n_mutated"], r["#Reads"], r["%Reads"]))
        out = pandas.DataFrame(rowsout, columns=["Aligned_Sequence", "Reference_Sequence", "Unedited", "n_deleted", "n_inserted", "n_mutated", "#Reads", "%Reads"])
        out = out.groupby(["Aligned_Sequence", "Reference_Sequence", "Unedited", "n_deleted", "n_inserted", "n_mutated"]).sum().reset_index().set_index("Aligned_Sequence")
        out = out.sort_values(by=["#Reads", "Aligned_Sequence", "Reference_Sequence"], ascending=[False, True, True])
        out["Unedited"] = out["Unedited"] > 0
        return out

    for ref in info["results"]["ref_names"]:
        ref_info = refs[ref]
        targets = [("custom", plot_center, plot_center, None)] if plot_center is not None else [
            (g or seq, cut, pc, None) for g, seq, cut, pc in zip(ref_info["sgRNA_names"], ref_info["sgRNA_sequences"], ref_info["sgRNA_cut_points"], ref_info["sgRNA_plot_cut_points"], strict=True)
        ]
        for label, cut, pcut, _unused in targets:
            if mode == "ambiguous":
                left = right = max(1, args.plot_window_size)
                frame = df.loc[df["Reference_Name"] == "AMBIGUOUS_" + ref]
            else:
                left, right, frame = plot_left, plot_right, df
            around = window(frame, cut, left, right)
            if around.shape[0] < 1:
                continue
            start = cut - left
            intervals = [(a - start - 1, b - start - 1) for a, b in ref_info["sgRNA_intervals"]]
            seq = ref_info["sequence"][cut - left + 1 : cut + right + 1]
            filename = f"{root}_{ref}_{label}"
            kwargs: dict[str, Any] = {} if mode == "ambiguous" else {"cut_point_ind": cut - start}
            plot.plot_alleles_table(
                seq, df_alleles=around, fig_filename_root=filename, MIN_FREQUENCY=min_freq, MAX_N_ROWS=rows, SAVE_ALSO_PNG=save_png,
                plot_cut_point=plot_cut_point, sgRNA_intervals=intervals, sgRNA_names=ref_info["sgRNA_names"],
                sgRNA_mismatches=ref_info["sgRNA_mismatches"], annotate_wildtype_allele=args.annotate_wildtype_allele, **kwargs,
            )
            produced += [str(p) for p in sorted(root.parent.glob(Path(filename).name + ".*"))]
    return {
        "status": "completed" if produced else "failed",
        "outputs": produced,
        "plot_count": len({Path(p).stem for p in produced}),
        "reads_in_run": int(df["#Reads"].sum()) if df.shape[0] else 0,
    }


def crispresso_plot_ambiguous_alleles(
    crispresso_folder: str,
    output_root: str,
    *,
    min_frequency: float | None = None,
    max_rows: int | None = None,
    plot_cut_point: bool = True,
    save_png: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Allele-table plot of the ambiguous reads of a finished CRISPResso run (port of ``plotAmbiguous.py``).

    Ambiguous reads align equally well to several amplicons; they are plotted against the first reference, one plot per
    reference and guide, to ``<output_root>_<reference>_<guide>.pdf`` (and ``.png`` with ``save_png``). The run needs
    ``--write_detailed_allele_table``. ``min_frequency`` and ``max_rows`` default to the run's own settings. The original
    script appends each plot name to the previous one; this port names every plot independently.
    """
    result = _plot_alleles(crispresso_folder, output_root, "ambiguous", min_frequency, max_rows, plot_cut_point, save_png, 0, 0, None)
    result["operation"] = "plotAmbiguous"
    return result


def crispresso_plot_allele_window(
    crispresso_folder: str,
    output_root: str,
    *,
    plot_left: int = 20,
    plot_right: int = 20,
    plot_center: int | None = None,
    min_frequency: float | None = None,
    max_rows: int | None = None,
    plot_cut_point: bool = True,
    save_png: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Allele-table plot over a custom window around each cut site (port of ``plotCustomAllelePlot.py``).

    ``plot_left`` and ``plot_right`` bases around every sgRNA cut point (default 20 each) or around ``plot_center`` (a
    reference position; then no per-guide plots are made). Plots go to ``<output_root>_<reference>_<guide>.pdf`` (and
    ``.png`` with ``save_png``). The run needs ``--write_detailed_allele_table``.
    """
    left = _cli.integer(plot_left, "plot_left", minimum=1)
    right = _cli.integer(plot_right, "plot_right", minimum=1)
    result = _plot_alleles(crispresso_folder, output_root, "custom", min_frequency, max_rows, plot_cut_point, save_png, left, right, plot_center)
    result["operation"] = "plotCustomAllelePlot"
    return result


__all__ = [
    "crispresso_aggregate",
    "crispresso_batch",
    "crispresso_compare",
    "crispresso_count_high_quality_bases",
    "crispresso_count_sgrna_specific_edits",
    "crispresso_filter_reads_by_sequence",
    "crispresso_plot_allele_window",
    "crispresso_plot_ambiguous_alleles",
    "crispresso_pooled",
    "crispresso_pooled_wgs_compare",
    "crispresso_quantify_editing",
    "crispresso_reconstitute_reads",
    "crispresso_wgs",
]

_register_public_verbose(globals(), tuple(__all__))
