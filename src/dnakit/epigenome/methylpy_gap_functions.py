"""methylpy operations as typed DNAKit facades: bisulfite-read alignment, methylation calling and DMR finding.

The user installs methylpy and the aligner/tool executables (samtools, bowtie2 or bowtie, cutadapt, bgzip and tabix,
wigToBigWig, picard). Each function calls methylpy's own function for the operation, with the upstream defaults kept:
only the arguments the caller sets are forwarded. DNAKit adds three things: the executables are checked first, so a
missing program raises BackendUnavailableError with an install hint (upstream would call sys.exit); upstream failures
are raised as BackendExecutionError; and the result is a JSON-serialisable dictionary with status, backend, version,
function name and output paths.

One upstream step is not called: methylpy's single-end quality_trim runs cutadapt through shell=True with interpolated
strings. methylpy_quality_trim_se runs the same cutadapt command as an argument list with no shell.
"""

from __future__ import annotations

import functools
import importlib
import re
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_MAX_TAIL = 4000
_ADAPTER = re.compile(r"^[ACGTNRYSWKMBDHVU]+$", re.IGNORECASE)
_INSTALL = {
    "samtools": "conda install -c bioconda samtools, or apt install samtools; put it on PATH or pass path_to_samtools",
    "bowtie2": "conda install -c bioconda bowtie2, or apt install bowtie2; put it on PATH or pass path_to_aligner",
    "bowtie2-build": "conda install -c bioconda bowtie2 (provides bowtie2-build); pass path_to_aligner if it is not on PATH",
    "bowtie": "conda install -c bioconda bowtie; put it on PATH or pass path_to_aligner",
    "bowtie-build": "conda install -c bioconda bowtie (provides bowtie-build); pass path_to_aligner if it is not on PATH",
    "cutadapt": "python -m pip install cutadapt, or conda install -c bioconda cutadapt; put it on PATH or pass path_to_cutadapt",
    "bgzip": "conda install -c bioconda htslib (provides bgzip and tabix); put them on PATH or pass path_to_bgzip",
    "tabix": "conda install -c bioconda htslib (provides bgzip and tabix); put them on PATH or pass path_to_tabix",
    "wigToBigWig": "conda install -c bioconda ucsc-wigtobigwig; put it on PATH or pass path_to_wigToBigWig",
    "picard": "conda install -c bioconda picard, and Java; pass path_to_picard if it is not on PATH",
}
_GSL_HINT = (
    "The bundled run_rms_tests.out needs the GSL shared libraries (libgsl.so.0 and libgslcblas.so.0): "
    "install GSL (apt install libgsl-dev, or conda install -c conda-forge gsl)."
)
_METHYLPY_HINT = "Install methylpy (python -m pip install methylpy) in the same environment as DNAKit."


def _upstream(module: str) -> Any:
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(f"{module} could not be imported: {exc}", hint=_METHYLPY_HINT) from exc


@functools.lru_cache(maxsize=1)
def _version() -> str:
    try:
        return metadata.version("methylpy")
    except metadata.PackageNotFoundError:
        return "unknown"


def _require(program: str, prefix: str | None = None) -> str:
    """Check an executable the way methylpy finds it: prefix + name when a directory prefix is given, else PATH."""
    if prefix:
        directory = prefix if prefix.endswith("/") else prefix + "/"
        candidate = Path(directory + program)
        if not candidate.is_file():
            raise BackendUnavailableError(f"{program} was not found at {candidate}.", hint=_INSTALL.get(program, ""))
        return str(candidate)
    found = shutil.which(program)
    if found is None and (Path(sys.executable).parent / program).is_file():
        found = str(Path(sys.executable).parent / program)
    if found is None:
        raise BackendUnavailableError(
            f"{program} was not found on PATH.", hint=_INSTALL.get(program, f"Install {program} and put it on PATH.")
        )
    return found


def _existing(path: str | Path, name: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return str(resolved)


def _existing_all(paths: Sequence[str | Path], name: str) -> list[str]:
    if isinstance(paths, (str, Path)) or not paths:
        raise ConfigurationError(f"{name} must be a non-empty list of files.")
    return [_existing(p, name) for p in paths]


def _new(path: str | Path, name: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return str(resolved)


def _forward(arguments: Mapping[str, Any]) -> dict[str, Any]:
    """Keep only the arguments the caller set, so upstream defaults apply to the rest."""
    return {key: value for key, value in arguments.items() if key != "verbose" and value is not None}


def _call(
    module: str,
    function: str,
    arguments: Mapping[str, Any],
    *,
    programs: Sequence[tuple[str, str | None]] = (),
    outputs: Sequence[str] = (),
    hint: str | None = None,
) -> dict[str, Any]:
    for program, prefix in programs:
        _require(program, prefix)
    target = getattr(_upstream(module), function)
    forwarded = _forward(arguments)
    try:
        target(**forwarded)
    except SystemExit as exc:
        raise BackendExecutionError(f"methylpy {function} stopped: {exc.code}") from exc
    except Exception as exc:
        message = f"methylpy {function} failed: {exc}"
        raise BackendExecutionError(message + (f" {hint}" if hint else "")) from exc
    return {
        "status": "completed",
        "backend": "methylpy",
        "version": _version(),
        "function": f"{module}.{function}",
        "outputs": [str(p) for p in outputs if Path(p).exists()],
        "arguments": sorted(forwarded),
    }


def _bundled_rms_binary() -> None:
    module = _upstream("methylpy.DMRfind")
    binary = Path(module.__file__).resolve().parent / "run_rms_tests.out"
    if not binary.is_file():
        raise BackendUnavailableError("methylpy does not ship run_rms_tests.out in this installation.", hint=_METHYLPY_HINT)


def methylpy_dmrfind(
    allc_files: Sequence[str],
    output_prefix: str,
    *,
    samples: Sequence[str] | None = None,
    mc_type: str | None = None,
    chroms: Sequence[str] | None = None,
    num_procs: int | None = None,
    min_cov: int | None = None,
    dmr_max_dist: int | None = None,
    sig_cutoff: float | None = None,
    num_sims: int | None = None,
    num_sig_tests: int | None = None,
    seed: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find differentially methylated regions between samples from allc files (methylpy DMRfind, bundled permutation binary)."""
    arguments = dict(locals())
    _bundled_rms_binary()
    arguments["allc_files"] = _existing_all(allc_files, "allc_files")
    return _call("methylpy.DMRfind", "DMRfind", arguments, hint=_GSL_HINT)


def methylpy_run_rms_tests(
    files: Sequence[str],
    output: str,
    samples: Sequence[str],
    *,
    min_cov: int | None = None,
    num_sims: int | str | None = None,
    num_sig_tests: int | str | None = None,
    seed: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the root-mean-square permutation tests of DMRfind through methylpy's bundled C++ binary."""
    arguments = dict(locals())
    _bundled_rms_binary()
    arguments["files"] = _existing_all(files, "files")
    return _call("methylpy.DMRfind", "run_rms_tests", arguments, hint=_GSL_HINT)


def methylpy_bgzip_allc_file(
    allc_file: str,
    *,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    buffer_line_number: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compress an allc methylation-call file with bgzip and index it with tabix (methylpy bgzip_allc_file)."""
    arguments = dict(locals())
    source = _existing(allc_file, "allc_file")
    arguments["allc_file"] = source
    return _call(
        "methylpy.utilities",
        "bgzip_allc_file",
        arguments,
        programs=[("bgzip", path_to_bgzip), ("tabix", path_to_tabix)],
    )


def methylpy_convert_allc_to_bigwig(
    input_allc_file: str,
    output_file: str,
    reference_fasta: str,
    *,
    mc_type: str | None = None,
    bin_size: int | None = None,
    path_to_wigToBigWig: str | None = None,
    path_to_samtools: str | None = None,
    min_bin_sites: int | None = None,
    min_bin_cov: int | None = None,
    max_site_cov: int | None = None,
    min_site_cov: int | None = None,
    remove_chr_prefix: bool | None = None,
    add_chr_prefix: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an allc file to a binned methylation bigWig track (methylpy convert_allc_to_bigwig: samtools + wigToBigWig)."""
    arguments = dict(locals())
    arguments["input_allc_file"] = _existing(input_allc_file, "input_allc_file")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.utilities",
        "convert_allc_to_bigwig",
        arguments,
        programs=[("samtools", path_to_samtools), ("wigToBigWig", path_to_wigToBigWig)],
        outputs=[output_file],
    )


def methylpy_quality_trim_se(
    inputf: Sequence[str],
    output: Sequence[str],
    *,
    quality_base: int | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    adapter_seq: str = "AGATCGGAAGAGCACACGTCTG",
    input_format: str | None = None,
    error_rate: float | None = None,
    max_adapter_removal: int | None = None,
    overlap_length: int | None = None,
    zero_cap: bool = False,
    path_to_cutadapt: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim adapters and low-quality ends of single-end bisulfite reads with cutadapt, one command per file (no shell)."""
    sources = _existing_all(inputf, "inputf")
    targets = [_new(o, "output") for o in output]
    if len(targets) != len(sources):
        raise ConfigurationError("Provide one output file for each input file.")
    if input_format is not None and input_format not in ("fasta", "fastq", "sra-fastq"):
        raise ConfigurationError("input_format must be 'fasta', 'fastq' or 'sra-fastq'.")
    if not _ADAPTER.match(adapter_seq):
        raise ConfigurationError("adapter_seq must be a DNA sequence.")
    program = _require("cutadapt", path_to_cutadapt)
    base = [program, "--quiet"]
    if zero_cap:
        base.append("-z")
    if input_format:
        base += ["-f", input_format]
    if error_rate:
        base += ["-e", str(error_rate)]
    if max_adapter_removal:
        base += ["-n", str(max_adapter_removal)]
    if overlap_length:
        base += ["-O", str(overlap_length)]
    if min_read_len:
        base += ["-m", str(min_read_len)]
    if min_qual_score:
        base += ["-q", str(min_qual_score)]
    if quality_base:
        base += [f"--quality-base={quality_base}"]
    base += ["-a", adapter_seq.upper()]
    logs: list[str] = []
    for source, target in zip(sources, targets):
        completed = subprocess.run(
            [*base, "-o", target, source], capture_output=True, text=True, shell=False, check=False, timeout=3600
        )
        logs.append((completed.stderr or "")[-_MAX_TAIL:])
        if completed.returncode != 0:
            raise BackendExecutionError(f"cutadapt failed on {source} with exit code {completed.returncode}.")
    return {
        "status": "completed",
        "backend": "methylpy",
        "version": _version(),
        "function": "methylpy.call_mc_se.quality_trim (argument-list cutadapt)",
        "outputs": [t for t in targets if Path(t).is_file()],
        "log_tail": logs,
    }


def methylpy_quality_trim_pe(
    inputf_read1: str,
    outputf_read1: str,
    inputf_read2: str,
    outputf_read2: str,
    *,
    quality_base: int | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    adapter_seq_read1: str | None = None,
    adapter_seq_read2: str | None = None,
    num_procs: int | None = None,
    input_format: str | None = None,
    error_rate: float | None = None,
    max_adapter_removal: int | None = None,
    overlap_length: int | None = None,
    zero_cap: bool | None = None,
    path_to_cutadapt: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim adapters and low-quality ends of paired-end bisulfite reads (methylpy quality_trim_pe, cutadapt).

    Upstream quality_trim_pe mis-handles plain string paths (it assigns the list conversion to the wrong names and then
    iterates the characters of each path), so the facade always passes one-element lists.
    """
    arguments = dict(locals())
    arguments["inputf_read1"] = [_existing(inputf_read1, "inputf_read1")]
    arguments["inputf_read2"] = [_existing(inputf_read2, "inputf_read2")]
    arguments["outputf_read1"] = [_new(outputf_read1, "outputf_read1")]
    arguments["outputf_read2"] = [_new(outputf_read2, "outputf_read2")]
    return _call(
        "methylpy.call_mc_pe",
        "quality_trim_pe",
        arguments,
        programs=[("cutadapt", path_to_cutadapt)],
        outputs=[outputf_read1, outputf_read2],
    )


def methylpy_se_run_methylation_pipeline(
    read_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    libraries: Sequence[str] | None = None,
    unmethylated_control: str | None = None,
    path_to_output: str | None = None,
    sig_cutoff: float | None = None,
    num_procs: int | None = None,
    sort_mem: str | None = None,
    num_upstr_bases: int | None = None,
    num_downstr_bases: int | None = None,
    generate_allc_file: bool | None = None,
    compress_output: bool | None = None,
    bgzip: bool | None = None,
    binom_test: bool | None = None,
    min_cov: int | None = None,
    trim_reads: bool | None = None,
    aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    remove_clonal: bool | None = None,
    remove_chr_prefix: bool | None = None,
    pbat: bool | None = None,
    path_to_samtools: str | None = None,
    path_to_aligner: str | None = None,
    path_to_cutadapt: str | None = None,
    path_to_picard: str | None = None,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    adapter_seq: str | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the single-end bisulfite pipeline: align to converted references, call methylated sites (methylpy run_methylation_pipeline)."""
    arguments = dict(locals())
    arguments["read_files"] = _existing_all(read_files, "read_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    programs: list[tuple[str, str | None]] = [
        ("samtools", path_to_samtools),
        (aligner or "bowtie2", path_to_aligner),
    ]
    if trim_reads is not False:
        programs.append(("cutadapt", path_to_cutadapt))
    if bgzip:
        programs += [("bgzip", path_to_bgzip), ("tabix", path_to_tabix)]
    return _call("methylpy.call_mc_se", "run_methylation_pipeline", arguments, programs=programs)


def methylpy_se_run_mapping(
    current_library: str,
    library_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    path_to_output: str | None = None,
    path_to_samtools: str | None = None,
    path_to_aligner: str | None = None,
    aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    pbat: bool | None = None,
    num_procs: int | None = None,
    trim_reads: bool | None = None,
    path_to_cutadapt: str | None = None,
    adapter_seq: str | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    keep_temp_files: bool | None = None,
    sort_mem: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align one library of single-end bisulfite reads to the converted references with bowtie2 or bowtie (methylpy run_mapping)."""
    arguments = dict(locals())
    arguments["library_files"] = _existing_all(library_files, "library_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    programs: list[tuple[str, str | None]] = [("samtools", path_to_samtools), (aligner or "bowtie2", path_to_aligner)]
    if trim_reads is not False:
        programs.append(("cutadapt", path_to_cutadapt))
    return _call("methylpy.call_mc_se", "run_mapping", arguments, programs=programs)


def methylpy_se_run_alignment(
    current_library: str,
    library_read_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    path_to_output: str | None = None,
    path_to_samtools: str | None = None,
    aligner: str | None = None,
    path_to_aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    num_procs: int | None = None,
    keep_temp_files: bool | None = None,
    sort_mem: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align single-end reads of one library on both converted references and merge the strands (methylpy run_alignment)."""
    arguments = dict(locals())
    arguments["library_read_files"] = _existing_all(library_read_files, "library_read_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    return _call(
        "methylpy.call_mc_se",
        "run_alignment",
        arguments,
        programs=[("samtools", path_to_samtools), (aligner or "bowtie2", path_to_aligner)],
    )


def methylpy_se_build_ref(
    input_files: Sequence[str],
    output: str,
    *,
    buffsize: int | None = None,
    aligner: str | None = None,
    path_to_aligner: str | None = None,
    num_procs: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the C-to-T and G-to-A converted reference FASTA files and their aligner indexes (methylpy build_ref)."""
    arguments = dict(locals())
    arguments["input_files"] = _existing_all(input_files, "input_files")
    builder = "bowtie-build" if (aligner or "bowtie2").lower() == "bowtie" else "bowtie2-build"
    return _call("methylpy.call_mc_se", "build_ref", arguments, programs=[(builder, path_to_aligner)])


def _pipeline_programs(
    aligner: str | None,
    path_to_aligner: str | None,
    path_to_samtools: str | None,
    trim_reads: bool | None,
    path_to_cutadapt: str | None,
    bgzip: bool | None,
    path_to_bgzip: str | None,
    path_to_tabix: str | None,
) -> list[tuple[str, str | None]]:
    programs: list[tuple[str, str | None]] = [("samtools", path_to_samtools), (aligner or "bowtie2", path_to_aligner)]
    if trim_reads is not False:
        programs.append(("cutadapt", path_to_cutadapt))
    if bgzip:
        programs += [("bgzip", path_to_bgzip), ("tabix", path_to_tabix)]
    return programs


def _calling_programs(
    path_to_samtools: str | None, bgzip: bool | None, path_to_bgzip: str | None, path_to_tabix: str | None
) -> list[tuple[str, str | None]]:
    programs: list[tuple[str, str | None]] = [("samtools", path_to_samtools)]
    if bgzip:
        programs += [("bgzip", path_to_bgzip), ("tabix", path_to_tabix)]
    return programs


def methylpy_remove_clonal_bam(
    input_bam: str,
    output_bam: str,
    metric: str,
    *,
    is_pe: bool | None = None,
    path_to_picard: str | None = None,
    java_options: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove PCR-duplicate (clonal) reads from a BAM file with Picard MarkDuplicates, run through java (methylpy remove_clonal_bam)."""
    arguments = dict(locals())
    source = _existing(input_bam, "input_bam")
    target = _new(output_bam, "output_bam")
    if not path_to_picard:
        raise ConfigurationError(
            "path_to_picard must name the folder that contains picard.jar (methylpy reads path_to_picard + '/picard.jar')."
        )
    jar = Path(path_to_picard) / "picard.jar"
    if not jar.is_file():
        raise BackendUnavailableError(f"picard.jar was not found at {jar}.", hint=_INSTALL["picard"])
    arguments.update(input_bam=source, output_bam=target)
    return _call("methylpy.call_mc_se", "remove_clonal_bam", arguments, programs=[("java", None)], outputs=[target])


def methylpy_se_find_multi_mappers(
    inputf: str,
    output: str,
    *,
    num_procs: int | None = None,
    min_mapq: int | None = None,
    keep_temp_files: bool | None = None,
    append: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split single-end alignments into uniquely and multiply mapped reads (methylpy find_multi_mappers, reads SAM output of the aligner)."""
    arguments = dict(locals())
    arguments["inputf"] = _existing(inputf, "inputf")
    return _call("methylpy.call_mc_se", "find_multi_mappers", arguments, outputs=[output])


def methylpy_se_merge_sorted_multimap(
    current_library: str,
    files: Sequence[str],
    prefix: str,
    reference_fasta: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sorted single-end alignment files of one library and keep the uniquely mapping reads (methylpy merge_sorted_multimap)."""
    arguments = dict(locals())
    arguments["files"] = _existing_all(files, "files")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call("methylpy.call_mc_se", "merge_sorted_multimap", arguments, programs=[("samtools", path_to_samtools)])


def methylpy_se_merge_sorted_multimap_max_mapq(
    current_library: str,
    files: Sequence[str],
    prefix: str,
    reference_fasta: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sorted single-end alignment files and keep the best-MAPQ alignment of multiply mapped reads (methylpy merge_sorted_multimap_max_mapq)."""
    arguments = dict(locals())
    arguments["files"] = _existing_all(files, "files")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.call_mc_se", "merge_sorted_multimap_max_mapq", arguments, programs=[("samtools", path_to_samtools)]
    )


def methylpy_se_merge_bam_files(
    input_files: Sequence[str],
    output: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several BAM files into one with a corrected header (methylpy merge_bam_files, samtools merge)."""
    arguments = dict(locals())
    arguments["input_files"] = _existing_all(input_files, "input_files")
    return _call("methylpy.call_mc_se", "merge_bam_files", arguments, programs=[("samtools", path_to_samtools)], outputs=[output])


def methylpy_se_call_methylated_sites(
    inputf: str,
    sample: str,
    reference_fasta: str,
    *,
    unmethylated_control: str | None = None,
    sig_cutoff: float | None = None,
    num_procs: int | None = None,
    num_upstr_bases: int | None = None,
    num_downstr_bases: int | None = None,
    generate_mpileup_file: bool | None = None,
    compress_output: bool | None = None,
    bgzip: bool | None = None,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    min_mapq: int | None = None,
    min_cov: int | None = None,
    binom_test: bool | None = None,
    path_to_samtools: str | None = None,
    remove_chr_prefix: bool | None = None,
    sort_mem: str | None = None,
    add_snp_info: bool | None = None,
    min_base_quality: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call methylated cytosines from a sorted single-end BAM with a binomial test against the reference (methylpy call_methylated_sites)."""
    arguments = dict(locals())
    arguments["inputf"] = _existing(inputf, "inputf")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.call_mc_se",
        "call_methylated_sites",
        arguments,
        programs=_calling_programs(path_to_samtools, bgzip, path_to_bgzip, path_to_tabix),
    )


def methylpy_se_call_methylated_sites_with_snp_info(
    inputf: str,
    sample: str,
    reference_fasta: str,
    *,
    unmethylated_control: float | None = None,
    sig_cutoff: float | None = None,
    num_procs: int | None = None,
    num_upstr_bases: int | None = None,
    num_downstr_bases: int | None = None,
    generate_mpileup_file: bool | None = None,
    compress_output: bool | None = None,
    bgzip: bool | None = None,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    min_mapq: int | None = None,
    min_cov: int | None = None,
    binom_test: bool | None = None,
    path_to_samtools: str | None = None,
    remove_chr_prefix: bool | None = None,
    sort_mem: str | None = None,
    add_snp_info: bool | None = None,
    min_base_quality: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call methylated cytosines of a single-end BAM and write the SNP information at each covered site (methylpy call_methylated_sites_with_SNP_info)."""
    arguments = dict(locals())
    arguments["inputf"] = _existing(inputf, "inputf")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.call_mc_se",
        "call_methylated_sites_with_SNP_info",
        arguments,
        programs=_calling_programs(path_to_samtools, bgzip, path_to_bgzip, path_to_tabix),
    )


def methylpy_pe_run_methylation_pipeline(
    read1_files: Sequence[str],
    read2_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    libraries: Sequence[str] | None = None,
    unmethylated_control: str | None = None,
    path_to_output: str | None = None,
    sig_cutoff: float | None = None,
    num_procs: int | None = None,
    sort_mem: str | None = None,
    num_upstr_bases: int | None = None,
    num_downstr_bases: int | None = None,
    generate_allc_file: bool | None = None,
    compress_output: bool | None = None,
    bgzip: bool | None = None,
    binom_test: bool | None = None,
    min_cov: int | None = None,
    trim_reads: bool | None = None,
    aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    remove_clonal: bool | None = None,
    remove_chr_prefix: bool | None = None,
    pbat: bool | None = None,
    path_to_samtools: str | None = None,
    path_to_aligner: str | None = None,
    path_to_cutadapt: str | None = None,
    path_to_picard: str | None = None,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    adapter_seq_read1: str | None = None,
    adapter_seq_read2: str | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the paired-end bisulfite pipeline: align to converted references, call methylated sites (methylpy run_methylation_pipeline_pe)."""
    arguments = dict(locals())
    arguments["read1_files"] = _existing_all(read1_files, "read1_files")
    arguments["read2_files"] = _existing_all(read2_files, "read2_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    programs = _pipeline_programs(aligner, path_to_aligner, path_to_samtools, trim_reads, path_to_cutadapt,
                                  bgzip, path_to_bgzip, path_to_tabix)
    return _call("methylpy.call_mc_pe", "run_methylation_pipeline_pe", arguments, programs=programs)


def methylpy_pe_run_mapping(
    current_library: str,
    library_read1_files: Sequence[str],
    library_read2_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    path_to_output: str | None = None,
    path_to_samtools: str | None = None,
    path_to_aligner: str | None = None,
    aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    num_procs: int | None = None,
    trim_reads: bool | None = None,
    path_to_cutadapt: str | None = None,
    adapter_seq_read1: str | None = None,
    adapter_seq_read2: str | None = None,
    min_qual_score: int | None = None,
    min_read_len: int | None = None,
    keep_temp_files: bool | None = None,
    sort_mem: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align one library of paired-end bisulfite reads to the converted references with bowtie2 or bowtie (methylpy run_mapping_pe)."""
    arguments = dict(locals())
    arguments["library_read1_files"] = _existing_all(library_read1_files, "library_read1_files")
    arguments["library_read2_files"] = _existing_all(library_read2_files, "library_read2_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    programs: list[tuple[str, str | None]] = [("samtools", path_to_samtools), (aligner or "bowtie2", path_to_aligner)]
    if trim_reads is not False:
        programs.append(("cutadapt", path_to_cutadapt))
    return _call("methylpy.call_mc_pe", "run_mapping_pe", arguments, programs=programs)


def methylpy_pe_run_alignment(
    current_library: str,
    library_read1_files: Sequence[str],
    library_read2_files: Sequence[str],
    sample: str,
    forward_reference: str,
    reverse_reference: str,
    reference_fasta: str,
    *,
    path_to_output: str | None = None,
    path_to_samtools: str | None = None,
    aligner: str | None = None,
    aligner_options: Sequence[str] | None = None,
    path_to_aligner: str | None = None,
    merge_by_max_mapq: bool | None = None,
    min_mapq: int | None = None,
    num_procs: int | None = None,
    keep_temp_files: bool | None = None,
    sort_mem: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align paired-end reads of one library on both converted references and merge the strands (methylpy run_alignment_pe)."""
    arguments = dict(locals())
    arguments["library_read1_files"] = _existing_all(library_read1_files, "library_read1_files")
    arguments["library_read2_files"] = _existing_all(library_read2_files, "library_read2_files")
    for name in ("forward_reference", "reverse_reference", "reference_fasta"):
        arguments[name] = _existing(arguments[name], name)
    return _call(
        "methylpy.call_mc_pe",
        "run_alignment_pe",
        arguments,
        programs=[("samtools", path_to_samtools), (aligner or "bowtie2", path_to_aligner)],
    )


def methylpy_pe_find_multi_mappers(
    inputf: str,
    output: str,
    *,
    num_procs: int | None = None,
    min_mapq: int | None = None,
    keep_temp_files: bool | None = None,
    append: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split paired-end alignments into uniquely and multiply mapped reads (methylpy find_multi_mappers_pe)."""
    arguments = dict(locals())
    arguments["inputf"] = _existing(inputf, "inputf")
    return _call("methylpy.call_mc_pe", "find_multi_mappers_pe", arguments, outputs=[output])


def methylpy_pe_merge_sorted_multimap(
    current_library: str,
    files: Sequence[str],
    prefix: str,
    reference_fasta: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sorted paired-end alignment files of one library and keep the uniquely mapping read pairs (methylpy merge_sorted_multimap_pe)."""
    arguments = dict(locals())
    arguments["files"] = _existing_all(files, "files")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call("methylpy.call_mc_pe", "merge_sorted_multimap_pe", arguments, programs=[("samtools", path_to_samtools)])


def methylpy_pe_merge_sorted_multimap_max_mapq(
    current_library: str,
    files: Sequence[str],
    prefix: str,
    reference_fasta: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sorted paired-end alignment files and keep the best-MAPQ pair of multiply mapped reads (methylpy merge_sorted_multimap_pe_max_mapq)."""
    arguments = dict(locals())
    arguments["files"] = _existing_all(files, "files")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.call_mc_pe", "merge_sorted_multimap_pe_max_mapq", arguments, programs=[("samtools", path_to_samtools)]
    )


def methylpy_pe_flip_read2_strand(
    input_file: str,
    output_file: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flip the strand of read 2 in every mapped read pair so that both reads report the original strand (methylpy flip_read2_strand)."""
    arguments = dict(locals())
    arguments["input_file"] = _existing(input_file, "input_file")
    return _call("methylpy.call_mc_pe", "flip_read2_strand", arguments, programs=[("samtools", path_to_samtools)],
                 outputs=[output_file])


def methylpy_pe_merge_bam_files(
    input_files: Sequence[str],
    output: str,
    *,
    path_to_samtools: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge several paired-end BAM files into one with a corrected header (methylpy merge_bam_files of call_mc_pe)."""
    arguments = dict(locals())
    arguments["input_files"] = _existing_all(input_files, "input_files")
    return _call("methylpy.call_mc_pe", "merge_bam_files", arguments, programs=[("samtools", path_to_samtools)],
                 outputs=[output])


def methylpy_pe_call_methylated_sites(
    inputf: str,
    sample: str,
    reference_fasta: str,
    *,
    unmethylated_control: str | None = None,
    sig_cutoff: float | None = None,
    num_procs: int | None = None,
    num_upstr_bases: int | None = None,
    num_downstr_bases: int | None = None,
    generate_mpileup_file: bool | None = None,
    compress_output: bool | None = None,
    bgzip: bool | None = None,
    path_to_bgzip: str | None = None,
    path_to_tabix: str | None = None,
    buffer_line_number: int | None = None,
    min_mapq: int | None = None,
    min_cov: int | None = None,
    binom_test: bool | None = None,
    path_to_samtools: str | None = None,
    remove_chr_prefix: bool | None = None,
    add_snp_info: bool | None = None,
    sort_mem: str | None = None,
    path_to_files: str | None = None,
    min_base_quality: int | None = None,
    keep_temp_files: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call methylated cytosines from a sorted paired-end BAM with a binomial test against the reference (methylpy call_methylated_sites_pe)."""
    arguments = dict(locals())
    arguments["inputf"] = _existing(inputf, "inputf")
    arguments["reference_fasta"] = _existing(reference_fasta, "reference_fasta")
    return _call(
        "methylpy.call_mc_pe",
        "call_methylated_sites_pe",
        arguments,
        programs=_calling_programs(path_to_samtools, bgzip, path_to_bgzip, path_to_tabix),
    )


__all__ = [
    "methylpy_bgzip_allc_file",
    "methylpy_convert_allc_to_bigwig",
    "methylpy_dmrfind",
    "methylpy_pe_call_methylated_sites",
    "methylpy_pe_find_multi_mappers",
    "methylpy_pe_flip_read2_strand",
    "methylpy_pe_merge_bam_files",
    "methylpy_pe_merge_sorted_multimap",
    "methylpy_pe_merge_sorted_multimap_max_mapq",
    "methylpy_pe_run_alignment",
    "methylpy_pe_run_mapping",
    "methylpy_pe_run_methylation_pipeline",
    "methylpy_quality_trim_pe",
    "methylpy_quality_trim_se",
    "methylpy_remove_clonal_bam",
    "methylpy_run_rms_tests",
    "methylpy_se_build_ref",
    "methylpy_se_call_methylated_sites",
    "methylpy_se_call_methylated_sites_with_snp_info",
    "methylpy_se_find_multi_mappers",
    "methylpy_se_merge_bam_files",
    "methylpy_se_merge_sorted_multimap",
    "methylpy_se_merge_sorted_multimap_max_mapq",
    "methylpy_se_run_alignment",
    "methylpy_se_run_mapping",
    "methylpy_se_run_methylation_pipeline",
]
_register_public_verbose(globals())
