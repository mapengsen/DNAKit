"""hifiasm command-line facade: haplotype-resolved de novo assembly with every assembly mode and report option.

The user installs hifiasm (``conda install -c bioconda hifiasm``, or build it from https://github.com/chhylp123/hifiasm);
0.25.0 is the tested version. ``hifiasm_assemble`` runs one ``hifiasm`` command whose mode is chosen by its arguments:
HiFi reads (default), ONT simplex reads (``ont``), HiFi plus ultra-long ONT reads (``ultra_long_reads``), Hi-C
partitioning (``hic_reads1``/``hic_reads2``), trio binning (``paternal_kmers``/``maternal_kmers``), diploid dual
scaffolding, purge-dups style duplicate haplotig purging, and the read-correction, overlap, breakpoint, misjoin,
telomere and chimera options. Every option is passed only when it is given, so hifiasm's own defaults apply otherwise.
No shell is used and the program is resolved from an explicit path, then PATH. The result is a JSON-serialisable
dictionary: ``status``, ``backend``, ``version``, ``command``, ``return_code``, ``outputs`` (every file hifiasm wrote
with the prefix), the assembly graphs, the BED reports and the contig summaries.

The simpler ``dnakit.assembly.assemble_hifiasm`` (HiFi reads, default options, primary contig graph) is kept as it is.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError

_PROGRAM = "hifiasm"
_HINT = (
    "Install hifiasm 0.25 or later (for example `conda install -c bioconda hifiasm`, or build it from "
    "https://github.com/chhylp123/hifiasm), then put it on PATH or pass hifiasm_path."
)


def _files(values: Sequence[str] | str | None, name: str, *, minimum: int = 1) -> list[Path]:
    if values is None:
        return []
    items = [values] if isinstance(values, str) else list(values)
    if len(items) < minimum:
        raise ConfigurationError(f"{name} must hold at least {minimum} file(s).")
    return [_cli.existing_file(item, name) for item in items]


def _optional_int(flag: str, value: int | None, name: str, *, minimum: int | None = 0) -> list[str]:
    return [] if value is None else [flag, str(_cli.integer(value, name, minimum=minimum))]


def _optional_float(flag: str, value: float | None, name: str, *, minimum: float = 0.0, maximum: float | None = None) -> list[str]:
    return [] if value is None else [flag, f"{_cli.number(value, name, minimum=minimum, maximum=maximum):g}"]


def _switch(on: bool, flag: str, name: str) -> list[str]:
    if not isinstance(on, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return [flag] if on else []


def _graph_summary(path: Path) -> dict[str, Any]:
    contigs = 0
    bases = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("S\t"):
                fields = line.rstrip("\n").split("\t")
                contigs += 1
                if len(fields) > 2 and fields[2] != "*":
                    bases += len(fields[2])
                else:
                    length = next((f[5:] for f in fields[3:] if f.startswith("LN:i:")), "0")
                    bases += int(length)
    return {"file": path.name, "segments": contigs, "total_bp": bases}


def hifiasm_assemble(
    reads_paths: Sequence[str],
    output_prefix: str,
    *,
    ont: bool = False,
    ultra_long_paths: Sequence[str] | None = None,
    hic_reads1_paths: Sequence[str] | None = None,
    hic_reads2_paths: Sequence[str] | None = None,
    paternal_kmers_path: str | None = None,
    maternal_kmers_path: str | None = None,
    paternal_reads_list_path: str | None = None,
    maternal_reads_list_path: str | None = None,
    dual_scaffold: bool = False,
    scaffold_max_gap: int | None = None,
    purge_level: int | None = None,
    purge_similarity: float | None = None,
    purge_min_overlaps: int | None = None,
    haplotype_count: int | None = None,
    write_error_corrected_reads: bool = False,
    write_overlaps_paf: bool = False,
    low_quality_percent: int | None = None,
    break_low_coverage: int | None = None,
    break_high_coverage: int | None = None,
    break_exact_overlap_rate: float | None = None,
    misjoin_min_size: int | None = None,
    telomere_motif: str | None = None,
    chimeric_read_support: int | None = None,
    chimeric_flank_length: int | None = None,
    min_ont_read_length: int | None = None,
    min_ont_mean_quality: int | None = None,
    min_ultra_long_length: int | None = None,
    ultra_long_error_rate: float | None = None,
    primary_assembly: bool = False,
    hom_coverage: int | None = None,
    genome_size: str | None = None,
    bloom_filter_bits: int | None = None,
    threads: int = 1,
    hifiasm_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """De novo assembly with hifiasm in any of its modes (``hifiasm``).

    ``reads_paths`` are HiFi reads in FASTA/FASTQ (optionally gzipped), or ONT simplex R10 FASTQ reads with
    ``ont=True`` (``--ont``). ``output_prefix`` is the output prefix (``-o``); hifiasm writes
    ``<prefix>.*`` files next to it and the call refuses to run when such files already exist.

    Modes (combinable as hifiasm allows): ``ultra_long_paths`` integrates ultra-long ONT reads for a hybrid assembly
    (``--ul``; ``min_ultra_long_length`` is ``--ul-cut``, ``ultra_long_error_rate`` is ``--ul-rate``);
    ``hic_reads1_paths`` and ``hic_reads2_paths`` give Hi-C partitioning of haplotypes (``--h1``/``--h2``;
    ``misjoin_min_size`` is ``--l-msjoin``, the minimum size of unitigs checked for misjoins, 0 disables);
    ``paternal_kmers_path`` and ``maternal_kmers_path`` are ``yak count`` k-mer dumps for trio binning (``-1``/``-2``;
    the optional read-name lists are ``-3``/``-4``); ``dual_scaffold`` outputs the self-scaffolded diploid assemblies
    (``--dual-scaf``; ``scaffold_max_gap`` is ``--scaf-gap``).

    Duplicate haplotig purging: ``purge_level`` is ``-l`` (0 none, 1 light, 2/3 aggressive), ``purge_similarity`` is
    ``-s``, ``purge_min_overlaps`` is ``-O`` and ``haplotype_count`` is ``--n-hap``. Outputs and reports:
    ``write_error_corrected_reads`` writes the corrected reads as FASTA (``--write-ec``, ``<prefix>.ec.fa``),
    ``write_overlaps_paf`` writes the read overlaps (``--write-paf``, ``<prefix>.ovlp.paf``) and
    ``low_quality_percent`` sets the BED report of contig regions with at least that percentage of inconsistent reads
    (``--lowQ``, 0 disables). Contig breaking: ``break_low_coverage`` (``--b-cov``) and ``break_high_coverage``
    (``--h-cov``) break contigs at positions with fewer / more times the coverage, using
    ``break_exact_overlap_rate`` (``--m-rate``). ``telomere_motif`` preserves telomeric contig ends (``--telo-m``,
    ``CCCTAA`` for human). ONT mode filters: ``chimeric_read_support`` is ``--chem-c`` (reads with at most this many
    supporting reads are masked as chimeric), ``chimeric_flank_length`` is ``--chem-f``, ``min_ont_read_length`` is
    ``--rl-cut`` and ``min_ont_mean_quality`` is ``--sc-cut``. General: ``primary_assembly`` (``--primary``),
    ``hom_coverage`` (``--hom-cov``), ``genome_size`` (``--hg-size``, for example ``"5m"``), ``bloom_filter_bits`` (``-f``;
    hifiasm's default of 37 reserves 16 GB of memory, ``0`` disables the filter for small genomes) and ``threads``
    (``-t``).
    """
    reads = _files(list(reads_paths) if not isinstance(reads_paths, str) else [reads_paths], "reads_paths")
    prefix = Path(output_prefix).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder for output_prefix does not exist: {prefix.parent}")
    existing = sorted(path.name for path in prefix.parent.glob(f"{prefix.name}.*"))
    if existing:
        raise ConfigurationError(f"Files with this output_prefix already exist and will not be overwritten: {existing[:5]}")
    ultra_long = _files(ultra_long_paths, "ultra_long_paths")
    hic1 = _files(hic_reads1_paths, "hic_reads1_paths")
    hic2 = _files(hic_reads2_paths, "hic_reads2_paths")
    if bool(hic1) != bool(hic2):
        raise ConfigurationError("hic_reads1_paths and hic_reads2_paths must be given together.")
    if (paternal_kmers_path is None) != (maternal_kmers_path is None):
        raise ConfigurationError("paternal_kmers_path and maternal_kmers_path must be given together.")
    if (paternal_reads_list_path is None) != (maternal_reads_list_path is None):
        raise ConfigurationError("paternal_reads_list_path and maternal_reads_list_path must be given together.")
    if hic1 and paternal_kmers_path is not None:
        raise ConfigurationError("Hi-C partitioning and trio binning cannot be combined.")
    if purge_level is not None and not 0 <= purge_level <= 3:
        raise ConfigurationError("purge_level must be 0, 1, 2 or 3.")

    args = ["-t", str(_cli.integer(threads, "threads", minimum=1, maximum=256)), "-o", str(prefix)]
    args += _switch(ont, "--ont", "ont")
    args += _optional_int("-f", bloom_filter_bits, "bloom_filter_bits", minimum=0)
    args += _optional_int("--hom-cov", hom_coverage, "hom_coverage", minimum=1)
    if genome_size is not None:
        args += ["--hg-size", _cli.text(genome_size, "genome_size")]
    args += _switch(primary_assembly, "--primary", "primary_assembly")
    if ultra_long:
        args += ["--ul", ",".join(str(path) for path in ultra_long)]
    args += _optional_float("--ul-rate", ultra_long_error_rate, "ultra_long_error_rate", maximum=1.0)
    args += _optional_int("--ul-cut", min_ultra_long_length, "min_ultra_long_length")
    if hic1:
        args += ["--h1", ",".join(str(path) for path in hic1), "--h2", ",".join(str(path) for path in hic2)]
    args += _optional_int("--l-msjoin", misjoin_min_size, "misjoin_min_size")
    if paternal_kmers_path is not None and maternal_kmers_path is not None:
        args += ["-1", str(_cli.existing_file(paternal_kmers_path, "paternal_kmers_path"))]
        args += ["-2", str(_cli.existing_file(maternal_kmers_path, "maternal_kmers_path"))]
    if paternal_reads_list_path is not None and maternal_reads_list_path is not None:
        args += ["-3", str(_cli.existing_file(paternal_reads_list_path, "paternal_reads_list_path"))]
        args += ["-4", str(_cli.existing_file(maternal_reads_list_path, "maternal_reads_list_path"))]
    args += _switch(dual_scaffold, "--dual-scaf", "dual_scaffold")
    args += _optional_int("--scaf-gap", scaffold_max_gap, "scaffold_max_gap")
    args += _optional_int("-l", purge_level, "purge_level")
    args += _optional_float("-s", purge_similarity, "purge_similarity", maximum=1.0)
    args += _optional_int("-O", purge_min_overlaps, "purge_min_overlaps", minimum=1)
    args += _optional_int("--n-hap", haplotype_count, "haplotype_count", minimum=1)
    args += _switch(write_error_corrected_reads, "--write-ec", "write_error_corrected_reads")
    args += _switch(write_overlaps_paf, "--write-paf", "write_overlaps_paf")
    args += _optional_int("--lowQ", low_quality_percent, "low_quality_percent")
    args += _optional_int("--b-cov", break_low_coverage, "break_low_coverage")
    args += _optional_int("--h-cov", break_high_coverage, "break_high_coverage", minimum=-1)
    args += _optional_float("--m-rate", break_exact_overlap_rate, "break_exact_overlap_rate", maximum=1.0)
    if telomere_motif is not None:
        args += ["--telo-m", _cli.text(telomere_motif, "telomere_motif")]
    args += _optional_int("--chem-c", chimeric_read_support, "chimeric_read_support")
    args += _optional_int("--chem-f", chimeric_flank_length, "chimeric_flank_length", minimum=1)
    args += _optional_int("--rl-cut", min_ont_read_length, "min_ont_read_length")
    args += _optional_int("--sc-cut", min_ont_mean_quality, "min_ont_mean_quality")
    args += [str(path) for path in reads]

    _cli.program(_PROGRAM, explicit=hifiasm_path, hint=_HINT)
    result = _cli.run(
        _PROGRAM,
        args,
        cwd=prefix.parent,
        expected=[],
        explicit=hifiasm_path,
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    produced = sorted(prefix.parent.glob(f"{prefix.name}.*"))
    result["operation"] = "hifiasm"
    result["outputs"] = [str(path) for path in produced]
    graphs = [path for path in produced if path.name.endswith(".gfa") and ".noseq." not in path.name]
    result["assembly_graphs"] = [str(path) for path in graphs]
    result["bed_reports"] = [str(path) for path in produced if path.suffix == ".bed"]
    result["error_corrected_reads"] = next((str(path) for path in produced if path.name.endswith(".ec.fa")), None)
    result["overlaps_paf"] = next((str(path) for path in produced if path.name.endswith(".ovlp.paf")), None)
    result["contig_summary"] = [_graph_summary(path) for path in graphs if path.name.endswith("p_ctg.gfa")]
    if result["return_code"] == 0 and not graphs:
        result["status"] = "failed"
    elif result["return_code"] == 0:
        result["status"] = "completed"
    return result


__all__ = ["hifiasm_assemble"]

_register_public_verbose(globals(), tuple(__all__))
