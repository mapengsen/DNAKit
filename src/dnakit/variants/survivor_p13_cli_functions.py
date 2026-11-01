"""SURVIVOR command-line facades: SV simulation and evaluation, SV VCF filtering, statistics and format conversion.

The user installs SURVIVOR (``conda install -c bioconda survivor``, or build it from
https://github.com/fritzsedlazeck/SURVIVOR); 1.0.7 is the tested version. Each function runs one SURVIVOR command with
an explicit typed signature; SURVIVOR's positional arguments are mapped to named parameters, no shell is used and the
program is resolved from an explicit path, then PATH. Outputs are **new** files (an existing file is never overwritten).
Results are JSON-serialisable dictionaries: ``status``, ``backend``, ``version``, ``command``, ``return_code``,
``outputs``, ``log_tail`` plus the parsed summary of the command where it prints one.

Merging SV VCF files is already provided by ``dnakit.variants.merge_structural_variants_survivor`` and is not repeated
here. ``genComp`` and ``parent_phasing`` take the merged multi-sample VCF written by that function. ``select_sample`` is not
wrapped: it segfaults after the first ranking step in SURVIVOR 1.0.7.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p13_cli as _p13
from dnakit.exceptions import ConfigurationError

_PROGRAM = "SURVIVOR"
_HINT = (
    "Install SURVIVOR 1.0.7 or later (for example `conda install -c bioconda survivor`, or build it from "
    "https://github.com/fritzsedlazeck/SURVIVOR), then put it on PATH or pass survivor_path."
)
_SIM_TYPES = (
    ("DUPLICATION", "duplication", (100, 10000, 3)),
    ("INDEL", "indel", (20, 500, 1)),
    ("TRANSLOCATION", "translocation", (1000, 3000, 2)),
    ("INVERSION", "inversion", (600, 800, 4)),
    ("INV_del", "inversion_deletion", (600, 800, 2)),
    ("INV_dup", "inversion_duplication", (600, 800, 2)),
)
_OVERALL = re.compile(
    r"Overall:\s+(\d+)\s+([\d/]+)\s+([\d/]+)\s+([\d/]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)"
)
_SV_TYPES = ("DEL", "DUP", "INV", "TRA", "INS")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check(survivor_path: str | None) -> None:
    """Resolve the program before anything is created, so a missing program leaves no files behind."""
    _cli.program(_PROGRAM, explicit=survivor_path, hint=_HINT)


def _run(
    command: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    survivor_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    result = _p13.run(
        _PROGRAM,
        [command, *arguments],
        cwd=cwd,
        expected=expected,
        explicit=survivor_path,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = f"SURVIVOR {command}"
    return result


def _size(value: int | None, name: str) -> int:
    """SURVIVOR uses -1 for a disabled size/support limit; ``None`` means disabled."""
    return -1 if value is None else _cli.integer(value, name, minimum=0)


def _vcf_records(path: Path) -> int:
    return _p13.count_lines(path, skip_prefixes=("#",)) if path.is_file() else 0


def _sim_parameter_text(values: dict[str, Sequence[int]]) -> str:
    lines = ["PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!"]
    for key, _, _ in _SIM_TYPES:
        low, high, number = values[key]
        lines += [
            f"{key}_minimum_length: {low}",
            f"{key}_maximum_length: {high}",
            f"{key}_number: {number}",
        ]
    return "\n".join(lines) + "\n"


def _triplet(value: Sequence[int], name: str) -> tuple[int, int, int]:
    if isinstance(value, str) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigurationError(f"{name} must be [minimum_length, maximum_length, number].")
    low = _cli.integer(value[0], f"{name} minimum_length", minimum=1)
    high = _cli.integer(value[1], f"{name} maximum_length", minimum=1)
    number = _cli.integer(value[2], f"{name} number", minimum=0)
    if low > high:
        raise ConfigurationError(f"{name} minimum_length must not exceed maximum_length.")
    return low, high, number


# ---------------------------------------------------------------------------
# Simulation and evaluation
# ---------------------------------------------------------------------------


def survivor_simulate_sv_genome(
    reference_path: str,
    output_prefix: str,
    *,
    duplication: Sequence[int] = (100, 10000, 3),
    indel: Sequence[int] = (20, 500, 1),
    translocation: Sequence[int] = (1000, 3000, 2),
    inversion: Sequence[int] = (600, 800, 4),
    inversion_deletion: Sequence[int] = (600, 800, 2),
    inversion_duplication: Sequence[int] = (600, 800, 2),
    snp_frequency: float = 0.0,
    coordinates_for: str = "simulated_reads",
    parameter_file: str | None = None,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate SVs and SNPs into a reference genome (``SURVIVOR simSV``).

    Every SV class is ``[minimum_length, maximum_length, number]``; the defaults are the values of the parameter file
    that SURVIVOR itself generates. ``snp_frequency`` (0-1) is the SNP mutation rate. ``coordinates_for`` is
    ``"simulated_reads"`` (0: reads will be simulated from the altered genome and mapped to the original reference) or
    ``"real_reads"`` (1: real reads will be mapped to the altered genome). Alternatively ``parameter_file`` is an
    existing SURVIVOR parameter file that replaces the six class arguments. SURVIVOR has no seed, so the placement of
    the events differs between runs. Writes ``<prefix>.fasta`` (altered genome), ``<prefix>.bed`` (events),
    ``<prefix>.insertions.fa``, ``<prefix>.vcf`` and, unless ``parameter_file`` is given, ``<prefix>.parameters.txt``.
    """
    reference = _cli.existing_file(reference_path, "reference_path")
    prefix = Path(output_prefix).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for output_prefix does not exist: {prefix.parent}")
    rate = _cli.number(snp_frequency, "snp_frequency", minimum=0.0, maximum=1.0)
    mode = _cli.choice(coordinates_for, ("simulated_reads", "real_reads"), "coordinates_for")
    outputs = [Path(f"{prefix}{suffix}") for suffix in (".fasta", ".bed", ".insertions.fa", ".vcf")]
    parameters = Path(f"{prefix}.parameters.txt")
    targets = outputs + ([parameters] if parameter_file is None else [])
    clash = [str(path) for path in targets if path.exists()]
    if clash:
        raise ConfigurationError(f"Output files already exist and will not be overwritten: {clash}")
    supplied = {
        "DUPLICATION": duplication,
        "INDEL": indel,
        "TRANSLOCATION": translocation,
        "INVERSION": inversion,
        "INV_del": inversion_deletion,
        "INV_dup": inversion_duplication,
    }
    checked = {key: _triplet(value, name) for (key, name, _), value in zip(_SIM_TYPES, supplied.values(), strict=True)}
    chosen = _cli.existing_file(parameter_file, "parameter_file") if parameter_file is not None else None
    _check(survivor_path)
    if chosen is None:
        parameters.write_text(_sim_parameter_text(checked), encoding="utf-8")
        chosen = parameters
    args = [str(reference), str(chosen), repr(rate), "0" if mode == "simulated_reads" else "1", str(prefix)]
    result = _run(
        "simSV",
        args,
        cwd=prefix.parent,
        expected=outputs,
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    if parameter_file is None and parameters.is_file():
        result["outputs"].append(str(parameters))
    result["parameter_file"] = str(chosen)
    result["event_count"] = _p13.count_lines(outputs[1]) if outputs[1].is_file() else 0
    return result


def _samtools_text(
    alignments: Path, reference: Path | None, samtools_path: str | None, work: Path, timeout_seconds: float
) -> Path:
    """Write the alignment records (no header) as SAM text; BAM/CRAM and MD-tag creation need samtools."""
    hint = "Install samtools (for example `conda install -c bioconda samtools`) or pass samtools_path."
    target = work / "alignments.sam"
    suffix = alignments.suffix.lower()
    if suffix in (".bam", ".cram") or reference is not None:
        if reference is not None:
            args = ["calmd", str(alignments), str(reference)]
        else:
            args = ["view", str(alignments)]
        with tempfile.TemporaryDirectory(dir=work) as scratch:
            raw = Path(scratch) / "records.sam"
            result = _p13.run_to_file(
                "samtools",
                args,
                output=raw,
                explicit=samtools_path,
                hint=hint,
                timeout_seconds=timeout_seconds,
            )
            if result["status"] != "completed":
                raise ConfigurationError(f"samtools could not read the alignments: {result['log_tail'][-500:]}")
            _strip_header(raw, target)
        return target
    _strip_header(alignments, target)
    return target


def _strip_header(source: Path, target: Path) -> None:
    with source.open(encoding="utf-8", errors="replace") as reader, target.open("w", encoding="utf-8") as writer:
        for line in reader:
            if not line.startswith("@"):
                writer.write(line)


def survivor_scan_read_error_profile(
    alignments_path: str,
    output_path: str,
    *,
    min_read_length: int = 1000,
    reference_path: str | None = None,
    survivor_path: str | None = None,
    samtools_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract a per-position error profile from aligned long reads (``SURVIVOR scanreads``).

    ``alignments_path`` is a SAM, BAM or CRAM file whose records carry the ``MD`` tag (SURVIVOR stops without it);
    with ``reference_path`` the MD tags are first computed with ``samtools calmd``. Alignments shorter than
    ``min_read_length`` bases are ignored. The profile (read-length distribution and, per read position, the
    probabilities of stop, match, mismatch, insertion and deletion) is the input of :func:`survivor_simulate_long_reads`
    and is written to ``output_path``. BAM/CRAM input and ``reference_path`` need samtools.
    """
    alignments = _cli.existing_file(alignments_path, "alignments_path")
    reference = _cli.existing_file(reference_path, "reference_path") if reference_path is not None else None
    minimum = _cli.integer(min_read_length, "min_read_length", minimum=1)
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    with tempfile.TemporaryDirectory(dir=out.parent) as scratch:
        records = _samtools_text(alignments, reference, samtools_path, Path(scratch), timeout_seconds)
        result = _p13.run_to_file(
            _PROGRAM,
            ["scanreads", str(minimum), str(out)],
            output=None,
            cwd=out.parent,
            stdin_path=records,
            extra_outputs=[out],
            explicit=survivor_path,
            hint=_HINT,
            timeout_seconds=timeout_seconds,
        )
    result["operation"] = "SURVIVOR scanreads"
    valid = re.search(r"Number of valid reads:\s*(\d+)", result["log_tail"])
    result["valid_reads"] = int(valid.group(1)) if valid else None
    return result


def survivor_simulate_long_reads(
    reference_path: str,
    error_profile_path: str,
    output_path: str,
    *,
    coverage: float = 2.0,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate PacBio- or ONT-like long reads from a reference (``SURVIVOR simreads``).

    ``error_profile_path`` is a profile written by :func:`survivor_scan_read_error_profile` (or one of SURVIVOR's
    shipped Oxford Nanopore or PacBio profiles); read lengths and error rates come from it. ``coverage`` is the
    simulated depth. The reads are written as FASTA to ``output_path``; each read name holds the chromosome, the
    left-most start and the strand. SURVIVOR has no seed, so the reads differ between runs.
    """
    reference = _cli.existing_file(reference_path, "reference_path")
    profile = _cli.existing_file(error_profile_path, "error_profile_path")
    depth = _cli.number(coverage, "coverage", minimum=0.01)
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _run(
        "simreads",
        [str(reference), str(profile), f"{depth:g}", str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    match = re.search(r"Num of reads:\s*(\d+)", result["log_tail"])
    result["read_count"] = int(match.group(1)) if match else None
    return result


def survivor_evaluate_sv_calls(
    vcf_path: str,
    simulated_bed_path: str,
    output_prefix: str,
    *,
    max_distance: int = 10,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate an SV call set against the events written by :func:`survivor_simulate_sv_genome` (``SURVIVOR eval``).

    Calls must match the simulated type and lie within ``max_distance`` bases of the simulated breakpoints. The
    correctly found calls go to ``<prefix>_right.vcf`` and the calls that were not simulated to
    ``<prefix>addition.vcf``. ``evaluation`` holds the simulated total, and per type (DEL, DUP, INV, TRA, INS) the
    found, missing and additional counts, the sensitivity and the false discovery rate that SURVIVOR prints.
    """
    vcf = _cli.existing_file(vcf_path, "vcf_path")
    bed = _cli.existing_file(simulated_bed_path, "simulated_bed_path")
    distance = _cli.integer(max_distance, "max_distance", minimum=0)
    prefix = Path(output_prefix).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for output_prefix does not exist: {prefix.parent}")
    outputs = [Path(f"{prefix}_right.vcf"), Path(f"{prefix}addition.vcf")]
    clash = [str(path) for path in outputs if path.exists()]
    if clash:
        raise ConfigurationError(f"Output files already exist and will not be overwritten: {clash}")
    _check(survivor_path)
    result = _run(
        "eval",
        [str(vcf), str(bed), str(distance), str(prefix)],
        cwd=prefix.parent,
        expected=[],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["outputs"] = [str(path) for path in outputs if path.is_file()]
    match = _OVERALL.search(result["log_tail"])
    if match is None:
        result["status"] = "failed"
        result["evaluation"] = None
        return result

    def per_type(text: str) -> dict[str, int]:
        return dict(zip(_SV_TYPES, (int(part) for part in text.split("/")), strict=False))

    result["evaluation"] = {
        "simulated": int(match.group(1)),
        "found": per_type(match.group(2)),
        "missing": per_type(match.group(3)),
        "additional": per_type(match.group(4)),
        "sensitivity": float(match.group(5)),
        "false_discovery_rate": float(match.group(6)),
    }
    return result


# ---------------------------------------------------------------------------
# Filtering, statistics, annotation
# ---------------------------------------------------------------------------


def survivor_filter_sv_vcf(
    vcf_path: str,
    output_path: str,
    *,
    ignore_regions_bed: str | None = None,
    min_size: int | None = None,
    max_size: int | None = None,
    min_allele_frequency: float = 0.0,
    min_read_support: int | None = None,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter an SV VCF by size, regions to ignore, allele frequency and read support (``SURVIVOR filter``).

    ``ignore_regions_bed`` removes SVs overlapping the regions (``None`` disables it); ``min_size``/``max_size`` bound
    the SV length and ``min_read_support`` is the minimum ``RE`` read count (``None`` disables each, SURVIVOR's ``-1``);
    ``min_allele_frequency`` (0-1) is the minimum allele frequency. The records before and after filtering are counted
    in ``records_in`` and ``records_out``.
    """
    vcf = _cli.existing_file(vcf_path, "vcf_path")
    regions = str(_cli.existing_file(ignore_regions_bed, "ignore_regions_bed")) if ignore_regions_bed else "NA"
    low = _size(min_size, "min_size")
    high = _size(max_size, "max_size")
    if low >= 0 and high >= 0 and low > high:
        raise ConfigurationError("min_size must not exceed max_size.")
    frequency = _cli.number(min_allele_frequency, "min_allele_frequency", minimum=0.0, maximum=1.0)
    support = _size(min_read_support, "min_read_support")
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _run(
        "filter",
        [str(vcf), regions, str(low), str(high), repr(frequency), str(support), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["records_in"] = _vcf_records(vcf)
    result["records_out"] = _vcf_records(out)
    return result


def survivor_sv_vcf_stats(
    vcf_path: str,
    summary_path: str,
    *,
    min_size: int | None = None,
    max_size: int | None = None,
    min_read_support: int | None = None,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary statistics of an SV VCF: counts by type and by length class (``SURVIVOR stats``).

    ``min_size``, ``max_size`` and ``min_read_support`` restrict the counted SVs (``None`` disables each). The length
    table (rows 0-50 bp, 50-100 bp, 100-1000 bp, 1000-10000 bp ...; columns Del, Dup, Inv, INS, TRA, UNK) is written
    to ``summary_path`` and returned as ``length_table``; ``totals`` has the number of SVs per type as printed by
    SURVIVOR.
    """
    vcf = _cli.existing_file(vcf_path, "vcf_path")
    low = _size(min_size, "min_size")
    high = _size(max_size, "max_size")
    support = _size(min_read_support, "min_read_support")
    _check(survivor_path)
    out = _cli.new_file(summary_path, "summary_path")
    result = _run(
        "stats",
        [str(vcf), str(low), str(high), str(support), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["totals"] = _stats_totals(result["log_tail"])
    result["length_table"] = _length_table(out) if out.is_file() else []
    return result


def _stats_totals(text: str) -> dict[str, int]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.startswith("Tot\t") and index + 1 < len(lines):
            names = line.split("\t")
            values = lines[index + 1].split("\t")
            return {name: int(value) for name, value in zip(names, values, strict=False) if value.strip().isdigit()}
    return {}


def _length_table(path: Path) -> list[dict[str, Any]]:
    lines = [line for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if not lines:
        return []
    names = lines[0].split("\t")
    rows: list[dict[str, Any]] = []
    for line in lines[1:]:
        fields = line.split("\t")
        row: dict[str, Any] = {names[0]: fields[0]}
        for name, value in zip(names[1:], fields[1:], strict=False):
            row[name] = int(value) if value.strip().lstrip("-").isdigit() else value
        rows.append(row)
    return rows


def survivor_annotate_sv_with_mummer(
    vcf_path: str,
    show_diff_paths: Sequence[str],
    output_path: str,
    *,
    max_distance: int = 100,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate an SV VCF with the breakpoints found by MUMmer (``SURVIVOR compMUMMer``).

    ``show_diff_paths`` are the output files of MUMmer's ``show-diff`` (default output with its header line, not the
    ``-H`` form); SURVIVOR is given them through a list file written next to the output as
    ``<output>.show_diff_files.txt``. SVs whose breakpoints lie within ``max_distance`` bases of a MUMmer breakpoint
    are annotated. ``changed_records`` counts the records that differ from the input VCF.
    """
    vcf = _cli.existing_file(vcf_path, "vcf_path")
    if isinstance(show_diff_paths, str) or not isinstance(show_diff_paths, Sequence) or not show_diff_paths:
        raise ConfigurationError("show_diff_paths must be a non-empty list of files.")
    diffs = [_cli.existing_file(path, "show_diff_paths") for path in show_diff_paths]
    distance = _cli.integer(max_distance, "max_distance", minimum=0)
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    listing = _cli.new_file(f"{out}.show_diff_files.txt", "show_diff list file")
    listing.write_text("".join(f"{path}\n" for path in diffs), encoding="utf-8")
    result = _run(
        "compMUMMer",
        [str(vcf), str(listing), str(distance), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["list_file"] = str(listing)
    if out.is_file():
        before = set(vcf.read_text(encoding="utf-8", errors="replace").splitlines())
        after = [line for line in out.read_text(encoding="utf-8", errors="replace").splitlines() if not line.startswith("#")]
        result["changed_records"] = sum(1 for line in after if line not in before)
    return result


# ---------------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------------


def survivor_sv_vcf_to_bed(
    vcf_path: str,
    output_path: str,
    *,
    min_size: int | None = None,
    max_size: int | None = None,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an SV VCF to a BEDPE-style BED file with an optional size range (``SURVIVOR vcftobed``).

    Each line holds both breakpoints (chromosome, start, end for each), the SV id, the type string and the strands.
    ``min_size``/``max_size`` select the SV length (``None`` disables each, SURVIVOR's ``-1``).
    """
    vcf = _cli.existing_file(vcf_path, "vcf_path")
    low = _size(min_size, "min_size")
    high = _size(max_size, "max_size")
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _run(
        "vcftobed",
        [str(vcf), str(low), str(high), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _p13.count_lines(out) if out.is_file() else 0
    return result


def survivor_convert_to_sv_vcf(
    input_path: str,
    output_path: str,
    *,
    input_format: str = "bed",
    sv_type: str | None = None,
    min_size: int = 0,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert SV calls of another program or format to an SV VCF (``SURVIVOR bedtovcf/smaptovcf/bedpetovcf/convertAssemblytics``).

    ``input_format`` selects the SURVIVOR command: ``"bed"`` (``bedtovcf``: ``chr<TAB>start<TAB>stop`` lines that are all
    of the SV type ``sv_type``, for example ``DEL``, ``INS``, ``DUP``, ``INV``), ``"smap"`` (``smaptovcf``: Bionano
    SMAP, beta), ``"bedpe"`` (``bedpetovcf``: Lumpy BEDPE, beta) or ``"assemblytics"`` (``convertAssemblytics``: the
    Assemblytics BED, keeping SVs of at least ``min_size`` bases). ``sv_type`` is required for ``"bed"`` only and
    ``min_size`` applies to ``"assemblytics"`` only.
    """
    source = _cli.existing_file(input_path, "input_path")
    kind = _cli.choice(input_format, ("bed", "smap", "bedpe", "assemblytics"), "input_format")
    if sv_type is not None and kind != "bed":
        raise ConfigurationError("sv_type only applies to input_format='bed'.")
    if kind == "bed":
        if sv_type is None:
            raise ConfigurationError("sv_type is required for input_format='bed'.")
        command, arguments = "bedtovcf", [str(source), _cli.text(sv_type, "sv_type"), ""]
    elif kind == "smap":
        command, arguments = "smaptovcf", [str(source), ""]
    elif kind == "bedpe":
        command, arguments = "bedpetovcf", [str(source), ""]
    else:
        command = "convertAssemblytics"
        arguments = [str(source), str(_cli.integer(min_size, "min_size", minimum=0)), ""]
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    arguments[-1] = str(out)
    result = _run(
        command,
        arguments,
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["input_format"] = kind
    result["record_count"] = _vcf_records(out)
    return result


def survivor_hapcut_to_phased_vcf(
    snp_vcf_path: str,
    hapcut_path: str,
    output_path: str,
    *,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a phased VCF from a HapCUT2 final haplotype file and the SNP VCF given to HapCUT2 (``SURVIVOR hapcuttovcf``).

    The phased SNPs get ``0|1`` style genotypes and a ``PS`` phase-set identifier from the HapCUT2 blocks; SNPs that
    HapCUT2 did not phase are copied unchanged.
    """
    snps = _cli.existing_file(snp_vcf_path, "snp_vcf_path")
    haplotypes = _cli.existing_file(hapcut_path, "hapcut_path")
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _run(
        "hapcuttovcf",
        [str(snps), str(haplotypes), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["record_count"] = _vcf_records(out)
    return result


def survivor_low_mapping_quality_regions(
    coverage_path: str,
    output_path: str,
    *,
    border_size: int = 10,
    min_coverage: int = 2,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster a MQ0 coverage table into BED regions to filter SVs in poorly mappable sequence (``SURVIVOR bincov``).

    ``coverage_path`` is ``samtools depth`` output of the low-mapping-quality reads; positions with coverage above
    ``min_coverage`` that are at most ``border_size`` bases apart are merged into one region. The BED goes to
    ``output_path`` and can be given to :func:`survivor_filter_sv_vcf` as ``ignore_regions_bed``.
    """
    coverage = _cli.existing_file(coverage_path, "coverage_path")
    border = _cli.integer(border_size, "border_size", minimum=0)
    minimum = _cli.integer(min_coverage, "min_coverage", minimum=0)
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _p13.run_to_file(
        _PROGRAM,
        ["bincov", str(coverage), str(border), str(minimum)],
        output=out,
        explicit=survivor_path,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    result["operation"] = "SURVIVOR bincov"
    result["region_count"] = _p13.count_lines(out) if out.is_file() else 0
    return result


# ---------------------------------------------------------------------------
# Multi-sample analysis of a merged VCF
# ---------------------------------------------------------------------------


def survivor_pairwise_overlap_matrix(
    merged_vcf_path: str,
    output_path: str,
    *,
    normalize: bool = False,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise SV overlap matrix across the samples of a merged multi-sample VCF (``SURVIVOR genComp``).

    Entry (i, j) counts the SVs supported by both samples; with ``normalize=True`` the counts are normalised by SURVIVOR.
    The matrix is written to ``output_path`` and returned as ``matrix``.
    """
    vcf = _cli.existing_file(merged_vcf_path, "merged_vcf_path")
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    flag = "1" if _flag(normalize, "normalize") else "0"
    result = _run(
        "genComp",
        [str(vcf), flag, str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["matrix"] = _matrix(out) if out.is_file() else []
    return result


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _matrix(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.strip():
            rows.append([float(value) for value in line.split()])
    return rows


def survivor_parental_phasing(
    merged_parental_vcf_path: str,
    hapcut_path: str,
    gatk_path: str,
    proband_snp_vcf_path: str,
    output_path: str,
    *,
    survivor_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign the SNPs of a proband to the parental haplotypes for SV phasing (``SURVIVOR parent_phasing``).

    Takes the merged parental VCF, the HapCUT2 final haplotype file, the GATK phasing file and the SNP VCF of the
    proband, and writes SURVIVOR's phasing matrix (one tab-separated line per SNP) to ``output_path``; the lines are
    returned in ``rows``.
    """
    parental = _cli.existing_file(merged_parental_vcf_path, "merged_parental_vcf_path")
    haplotypes = _cli.existing_file(hapcut_path, "hapcut_path")
    gatk = _cli.existing_file(gatk_path, "gatk_path")
    proband = _cli.existing_file(proband_snp_vcf_path, "proband_snp_vcf_path")
    _check(survivor_path)
    out = _cli.new_file(output_path, "output_path")
    result = _run(
        "parent_phasing",
        [str(parental), str(haplotypes), str(gatk), str(proband), str(out)],
        cwd=out.parent,
        expected=[out],
        survivor_path=survivor_path,
        timeout_seconds=timeout_seconds,
    )
    result["rows"] = (
        [line.split("\t") for line in out.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
        if out.is_file()
        else []
    )
    return result


__all__ = [
    "survivor_annotate_sv_with_mummer",
    "survivor_convert_to_sv_vcf",
    "survivor_evaluate_sv_calls",
    "survivor_filter_sv_vcf",
    "survivor_hapcut_to_phased_vcf",
    "survivor_low_mapping_quality_regions",
    "survivor_pairwise_overlap_matrix",
    "survivor_parental_phasing",
    "survivor_scan_read_error_profile",
    "survivor_simulate_long_reads",
    "survivor_simulate_sv_genome",
    "survivor_sv_vcf_stats",
    "survivor_sv_vcf_to_bed",
]

_register_public_verbose(globals(), tuple(__all__))
