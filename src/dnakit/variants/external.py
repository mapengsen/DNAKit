"""JSON-safe facades for caller-installed variant command-line tools."""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import (
    ExecutionConfig,
    artifacts,
    input_path,
    positive,
    run,
    run_to_file,
    workspace,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools.serialization import serialize_for_agent


def call_freebayes(
    reference: str | Path,
    alignments: str | Path | Sequence[str | Path],
    output_directory: str | Path,
    *,
    ploidy: int = 2,
    min_alternate_fraction: float | None = None,
    min_alternate_count: int | None = None,
    program: str = "freebayes",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call small variants from one or more indexed BAM/CRAM files with FreeBayes."""

    reference_path = input_path(reference)
    alignment_values = (
        (alignments,) if isinstance(alignments, (str, Path)) else tuple(alignments)
    )
    if not alignment_values:
        raise ConfigurationError("alignments must contain at least one BAM or CRAM path.")
    alignment_paths = tuple(input_path(path) for path in alignment_values)
    positive(ploidy, "ploidy")
    _validate_fraction(min_alternate_fraction, "min_alternate_fraction")
    if min_alternate_count is not None:
        positive(min_alternate_count, "min_alternate_count")
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / "variants.vcf"
    arguments = ["-f", str(reference_path), "-p", str(ploidy)]
    if min_alternate_fraction is not None:
        arguments.extend(("--min-alternate-fraction", str(min_alternate_fraction)))
    if min_alternate_count is not None:
        arguments.extend(("--min-alternate-count", str(min_alternate_count)))
    arguments.extend(str(path) for path in alignment_paths)
    log = run_to_file(program, arguments, target, config)
    return _envelope(
        "call_freebayes",
        "freebayes",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "reference": reference_path,
            "alignments": alignment_paths,
            "ploidy": ploidy,
            "log": log,
        },
        verbose=verbose,
    )


def load_gemini_database(
    source_vcf: str | Path,
    output_directory: str | Path,
    *,
    annotation_type: str | None = None,
    cores: int = 1,
    program: str = "gemini",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Load a VCF into a caller-installed GEMINI variant database.

    GEMINI annotation assets and any external annotation downloads remain the
    caller's responsibility.  The adapter only exposes the documented
    ``gemini load`` boundary and never downloads those assets implicitly.
    """

    source_path = input_path(source_vcf)
    _validate_program(program, "program")
    _validate_integer(cores, "cores", minimum=1, maximum=256)
    if annotation_type is not None:
        _validate_simple_text(annotation_type, "annotation_type", maximum=64)
    config = ExecutionConfig(
        threads=cores,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    database = directory / "variants.db"
    arguments = ["load", "-v", str(source_path)]
    if annotation_type is not None:
        arguments.extend(("-t", annotation_type))
    if cores > 1:
        arguments.extend(("--cores", str(cores)))
    arguments.append(str(database))
    log = run(program, arguments, directory, config, outputs=(database,))
    return _envelope(
        "load_gemini_database",
        "GEMINI",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source_vcf": source_path,
            "database": database,
            "annotation_type": annotation_type,
            "cores": cores,
            "log": log,
        },
        verbose=verbose,
    )


def query_gemini_database(
    database: str | Path,
    query: str,
    output_directory: str | Path,
    *,
    header: bool = True,
    result_format: str | None = None,
    program: str = "gemini",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Run a bounded GEMINI SQL query and capture its tabular stdout.

    ``query`` is passed as one argv value to GEMINI; no shell is used.  The
    optional format is restricted to the documented default and
    ``sampledetail`` views.
    """

    database_path = input_path(database)
    if not isinstance(query, str) or not query.strip() or len(query) > 32_768 or "\x00" in query:
        raise ConfigurationError("query must be non-empty SQL text of at most 32768 characters.")
    _validate_program(program, "program")
    _validate_bool(header, "header")
    if result_format is not None:
        result_format = _validate_choice(result_format, "result_format", {"default", "sampledetail"})
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    output = directory / "query.tsv"
    arguments = ["query"]
    if header:
        arguments.append("--header")
    if result_format is not None:
        arguments.extend(("--format", result_format))
    arguments.extend(("-q", query, str(database_path)))
    log = run_to_file(program, arguments, output, config)
    return _envelope(
        "query_gemini_database",
        "GEMINI",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "database": database_path,
            "query": query,
            "header": header,
            "result_format": result_format,
            "output": output,
            "log": log,
        },
        verbose=verbose,
    )


def annotate_vcf_ensembl_vep(
    source_vcf: str | Path,
    output_directory: str | Path,
    *,
    species: str = "homo_sapiens",
    assembly: str | None = None,
    cache_directory: str | Path | None = None,
    fasta: str | Path | None = None,
    offline: bool = True,
    output_format: str = "vcf",
    fork: int = 1,
    program: str = "vep",
    timeout_seconds: float = 7_200.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Annotate a VCF with caller-installed Ensembl Variant Effect Predictor.

    VEP caches, transcript/reference resources and any remote database access
    remain the caller's responsibility.  The adapter exposes a small local
    boundary and never installs or downloads those resources implicitly.
    """

    source_path = input_path(source_vcf)
    _validate_program(program, "program")
    _validate_simple_text(species, "species", maximum=128)
    if assembly is not None:
        _validate_simple_text(assembly, "assembly", maximum=64)
    _validate_bool(offline, "offline")
    _validate_choice(output_format, "output_format", {"vcf", "tab", "json"})
    _validate_integer(fork, "fork", minimum=1, maximum=32)
    cache_path: Path | None = None
    if cache_directory is not None:
        cache_path = Path(cache_directory).expanduser().resolve()
        if not cache_path.is_dir():
            raise ConfigurationError(f"cache_directory must be an existing directory: {cache_path}")
    fasta_path = input_path(fasta) if fasta is not None else None
    config = ExecutionConfig(
        threads=max(1, fork),
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    output = directory / {
        "vcf": "annotated.vcf",
        "tab": "annotated.tsv",
        "json": "annotated.json",
    }[output_format]
    arguments = [
        "--input_file",
        str(source_path),
        "--output_file",
        str(output),
        "--species",
        species,
    ]
    if output_format == "vcf":
        arguments.append("--vcf")
    elif output_format == "tab":
        arguments.append("--tab")
    else:
        arguments.append("--json")
    if assembly is not None:
        arguments.extend(("--assembly", assembly))
    if offline:
        arguments.append("--offline")
    if cache_path is not None:
        arguments.extend(("--cache", "--dir_cache", str(cache_path)))
    if fasta_path is not None:
        arguments.extend(("--fasta", str(fasta_path)))
    if fork > 1:
        arguments.extend(("--fork", str(fork)))
    log = run(program, arguments, directory, config, outputs=(output,))
    return _envelope(
        "annotate_vcf_ensembl_vep",
        "Ensembl VEP",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source_vcf": source_path,
            "output": output,
            "species": species,
            "assembly": assembly,
            "cache_directory": cache_path,
            "fasta": fasta_path,
            "offline": offline,
            "output_format": output_format,
            "fork": fork,
            "log": log,
        },
        verbose=verbose,
    )


def call_gatk_haplotypecaller(
    reference: str | Path,
    alignment: str | Path,
    output_directory: str | Path,
    *,
    emit_gvcf: bool = False,
    output_mode: str = "EMIT_VARIANTS_ONLY",
    min_confidence: float | None = None,
    ploidy: int | None = None,
    pair_hmm_threads: int | None = None,
    interval: str | Path | None = None,
    java_max_heap_gb: int | None = None,
    program: str = "gatk",
    timeout_seconds: float = 7_200.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call germline variants with caller-installed GATK4 HaplotypeCaller.

    This is a bounded GATK4 boundary for one reference, one coordinate-sorted
    BAM/CRAM and one VCF/GVCF output. It does not bundle GATK, Java, reference
    indexes or the broader GATK tool collection.
    """

    reference_path = input_path(reference)
    alignment_path = input_path(alignment)
    interval_path = input_path(interval) if interval is not None else None
    _validate_program(program, "program")
    _validate_bool(emit_gvcf, "emit_gvcf")
    output_mode = _validate_choice(
        output_mode,
        "output_mode",
        {"EMIT_VARIANTS_ONLY", "EMIT_ALL_CONFIDENT_SITES", "EMIT_ALL_ACTIVE_SITES"},
    )
    _validate_number(min_confidence, "min_confidence", minimum=0.0, maximum=10_000.0)
    _validate_integer(ploidy, "ploidy", minimum=1, maximum=64)
    _validate_integer(pair_hmm_threads, "pair_hmm_threads", minimum=1, maximum=256)
    _validate_integer(java_max_heap_gb, "java_max_heap_gb", minimum=1, maximum=512)
    config = ExecutionConfig(
        threads=pair_hmm_threads or 1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / ("variants.g.vcf.gz" if emit_gvcf else "variants.vcf.gz")
    arguments: list[str] = []
    if java_max_heap_gb is not None:
        arguments.extend(("--java-options", f"-Xmx{java_max_heap_gb}g"))
    arguments.extend(("HaplotypeCaller", "-R", str(reference_path), "-I", str(alignment_path)))
    arguments.extend(("-O", str(target), "--output-mode", output_mode))
    if emit_gvcf:
        arguments.extend(("-ERC", "GVCF"))
    if min_confidence is not None:
        arguments.extend(("--standard-min-confidence-threshold-for-calling", str(min_confidence)))
    if ploidy is not None:
        arguments.extend(("--sample-ploidy", str(ploidy)))
    if pair_hmm_threads is not None:
        arguments.extend(("--native-pair-hmm-threads", str(pair_hmm_threads)))
    if interval_path is not None:
        arguments.extend(("-L", str(interval_path)))
    log = run(program, arguments, directory, config, outputs=(target,))
    return _envelope(
        "call_gatk_haplotypecaller",
        "gatk4",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "reference": reference_path,
            "alignment": alignment_path,
            "interval": interval_path,
            "output": target,
            "emit_gvcf": emit_gvcf,
            "log": log,
        },
        verbose=verbose,
    )


def call_deepvariant(
    reference: str | Path,
    alignment: str | Path,
    output_directory: str | Path,
    *,
    model_type: str = "WGS",
    num_shards: int | None = None,
    emit_gvcf: bool = False,
    program: str = "run_deepvariant",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call variants with caller-installed DeepVariant.

    This bounded facade exposes the documented reference, indexed BAM/CRAM,
    model and output VCF/GVCF boundary. DeepVariant models, indexes and its
    container or binary remain caller-installed dependencies.
    """

    reference_path = input_path(reference)
    alignment_path = input_path(alignment)
    _validate_program(program, "program")
    model_type = _validate_choice(model_type, "model_type", {"WGS", "WES", "PACBIO", "ONT_R104"})
    _validate_bool(emit_gvcf, "emit_gvcf")
    _validate_integer(num_shards, "num_shards", minimum=1, maximum=512)
    config = ExecutionConfig(
        threads=num_shards or 1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    vcf = directory / "output.vcf.gz"
    gvcf = directory / "output.g.vcf.gz"
    arguments = [
        f"--model_type={model_type}",
        f"--ref={reference_path}",
        f"--reads={alignment_path}",
        f"--output_vcf={vcf}",
    ]
    outputs = [vcf]
    if emit_gvcf:
        arguments.append(f"--output_gvcf={gvcf}")
        outputs.append(gvcf)
    if num_shards is not None:
        arguments.append(f"--num_shards={num_shards}")
    log = run(program, arguments, directory, config, outputs=tuple(outputs))
    return _envelope(
        "call_deepvariant",
        "deepvariant",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "reference": reference_path,
            "alignment": alignment_path,
            "model_type": model_type,
            "output_vcf": vcf,
            "output_gvcf": gvcf if emit_gvcf else None,
            "num_shards": num_shards,
            "log": log,
        },
        verbose=verbose,
    )


def filter_vcf_vcftools(
    source: str | Path,
    output_directory: str | Path,
    *,
    gzipped: bool = False,
    remove_indels: bool = False,
    remove_snps: bool = False,
    max_missing: float | None = None,
    min_allele_frequency: float | None = None,
    program: str = "vcftools",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Apply a bounded VCFtools filter set and write a recoded VCF."""

    source_path = input_path(source)
    _validate_fraction(max_missing, "max_missing")
    _validate_fraction(min_allele_frequency, "min_allele_frequency")
    if remove_indels and remove_snps:
        raise ConfigurationError("remove_indels and remove_snps cannot both be true.")
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    prefix = directory / "filtered"
    output_vcf = Path(str(prefix) + ".recode.vcf")
    arguments = ["--gzvcf" if gzipped else "--vcf", str(source_path)]
    if remove_indels:
        arguments.append("--remove-indels")
    if remove_snps:
        arguments.append("--remove-snps")
    if max_missing is not None:
        arguments.extend(("--max-missing", str(max_missing)))
    if min_allele_frequency is not None:
        arguments.extend(("--maf", str(min_allele_frequency)))
    arguments.extend(("--recode", "--recode-INFO-all", "--out", str(prefix)))
    log = run(program, arguments, directory, config, outputs=(output_vcf,))
    return _envelope(
        "filter_vcf_vcftools",
        "vcftools",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "output": output_vcf,
            "log": log,
        },
        verbose=verbose,
    )


def annotate_vcf_vcfanno(
    source: str | Path,
    configuration: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    lua_script: str | Path | None = None,
    ends: bool = False,
    permissive_overlap: bool = False,
    program: str = "vcfanno",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Annotate a VCF with caller-provided vcfanno TOML configuration."""

    source_path = input_path(source)
    configuration_path = input_path(configuration)
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be an integer in [1, 256].")
    if not isinstance(ends, bool) or not isinstance(permissive_overlap, bool):
        raise ConfigurationError("ends and permissive_overlap must be booleans.")
    lua_path = input_path(lua_script) if lua_script is not None else None
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    target = directory / "annotated.vcf"
    arguments = ["-p", str(threads)]
    if ends:
        arguments.append("-ends")
    if permissive_overlap:
        arguments.append("-permissive-overlap")
    if lua_path is not None:
        arguments.extend(("-lua", str(lua_path)))
    arguments.extend((str(configuration_path), str(source_path)))
    log = run_to_file(program, arguments, target, config)
    return _envelope(
        "annotate_vcf_vcfanno",
        "vcfanno",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "source": source_path,
            "configuration": configuration_path,
            "output": target,
            "threads": threads,
            "lua_script": lua_path,
            "log": log,
        },
        verbose=verbose,
    )


def merge_structural_variants_survivor(
    sources: Sequence[str | Path],
    output_directory: str | Path,
    *,
    max_distance: int = 100,
    min_support: int = 1,
    type_support: bool = True,
    strand_support: bool = True,
    min_sv_size: int = 30,
    program: str = "SURVIVOR",
    timeout_seconds: float = 3_600.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Merge structural-variant VCFs with the bounded SURVIVOR merge mode."""

    source_values = tuple(sources)
    if len(source_values) < 2:
        raise ConfigurationError("sources must contain at least two VCF paths.")
    source_paths = tuple(input_path(path) for path in source_values)
    _validate_bounded_integer(max_distance, "max_distance", minimum=0, maximum=10_000_000)
    _validate_bounded_integer(min_support, "min_support", minimum=1, maximum=10_000)
    if not isinstance(type_support, bool) or not isinstance(strand_support, bool):
        raise ConfigurationError("type_support and strand_support must be booleans.")
    _validate_bounded_integer(min_sv_size, "min_sv_size", minimum=0, maximum=10_000_000)
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    source_list = directory / "survivor_inputs.txt"
    source_list.write_text("\n".join(str(path) for path in source_paths) + "\n", encoding="utf-8")
    target = directory / "merged.vcf"
    arguments = [
        "merge",
        str(source_list),
        str(max_distance),
        str(min_support),
        "1" if type_support else "0",
        "1" if strand_support else "0",
        "0",
        str(min_sv_size),
        str(target),
    ]
    log = run(program, arguments, directory, config, outputs=(target,))
    return _envelope(
        "merge_structural_variants_survivor",
        "SURVIVOR",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "sources": source_paths,
            "source_list": source_list,
            "output": target,
            "max_distance": max_distance,
            "min_support": min_support,
            "type_support": type_support,
            "strand_support": strand_support,
            "min_sv_size": min_sv_size,
            "log": log,
        },
        verbose=verbose,
    )


def call_cutesv(
    alignment: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    threads: int = 1,
    sample: str | None = None,
    program: str = "cuteSV",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call long-read structural variants with caller-installed cuteSV."""

    alignment_path = input_path(alignment)
    reference_path = input_path(reference)
    _validate_program(program, "program")
    _validate_integer(threads, "threads", minimum=1, maximum=256)
    if sample is not None:
        _validate_simple_text(sample, "sample", maximum=256)
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    output_vcf = directory / "variants.vcf"
    work_directory = directory / "work"
    work_directory.mkdir(parents=True, exist_ok=True)
    arguments = [
        str(alignment_path),
        str(reference_path),
        str(output_vcf),
        str(work_directory),
        "--threads",
        str(threads),
    ]
    if sample is not None:
        arguments.extend(("--sample", sample))
    log = run(program, arguments, directory, config, outputs=(output_vcf,))
    return _envelope(
        "call_cutesv",
        "cuteSV",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "alignment": alignment_path,
            "reference": reference_path,
            "output": output_vcf,
            "work_directory": work_directory,
            "threads": threads,
            "sample": sample,
            "log": log,
        },
        verbose=verbose,
    )


def call_sniffles2(
    alignment: str | Path,
    output_directory: str | Path,
    *,
    reference: str | Path | None = None,
    threads: int = 1,
    sample: str | None = None,
    min_sv_length: int | None = None,
    program: str = "sniffles",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Call long-read structural variants with caller-installed Sniffles2."""

    alignment_path = input_path(alignment)
    reference_path = input_path(reference) if reference is not None else None
    _validate_program(program, "program")
    _validate_integer(threads, "threads", minimum=1, maximum=256)
    _validate_integer(min_sv_length, "min_sv_length", minimum=1, maximum=10_000_000)
    if sample is not None:
        _validate_simple_text(sample, "sample", maximum=256)
    config = ExecutionConfig(
        threads=threads,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    output_vcf = directory / "variants.vcf"
    arguments = [
        "--input",
        str(alignment_path),
        "--vcf",
        str(output_vcf),
        "--threads",
        str(threads),
    ]
    if reference_path is not None:
        arguments.extend(("--reference", str(reference_path)))
    if sample is not None:
        arguments.extend(("--sample-id", sample))
    if min_sv_length is not None:
        arguments.extend(("--minsvlen", str(min_sv_length)))
    log = run(program, arguments, directory, config, outputs=(output_vcf,))
    return _envelope(
        "call_sniffles2",
        "sniffles2",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "alignment": alignment_path,
            "reference": reference_path,
            "output": output_vcf,
            "threads": threads,
            "sample": sample,
            "min_sv_length": min_sv_length,
            "log": log,
        },
        verbose=verbose,
    )


def call_delly(
    alignment: str | Path,
    reference: str | Path,
    output_directory: str | Path,
    *,
    mode: str = "sr",
    program: str = "delly",
    timeout_seconds: float = 86_400.0,
    max_output_bytes: int = 1_000_000,
    max_artifact_bytes: int = 10_000_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Discover short-read structural variants with caller-installed DELLY."""

    alignment_path = input_path(alignment)
    reference_path = input_path(reference)
    _validate_program(program, "program")
    mode = _validate_choice(mode, "mode", {"sr", "lr"})
    config = ExecutionConfig(
        threads=1,
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        max_artifact_bytes=max_artifact_bytes,
    )
    directory = workspace(output_directory)
    output_vcf = directory / "variants.vcf"
    arguments = [mode, "-g", str(reference_path), str(alignment_path)]
    log = run_to_file(program, arguments, output_vcf, config)
    return _envelope(
        "call_delly",
        "DELLY",
        {
            "directory": directory,
            "artifacts": artifacts(directory),
            "alignment": alignment_path,
            "reference": reference_path,
            "output": output_vcf,
            "mode": mode,
            "log": log,
        },
        verbose=verbose,
    )


def _envelope(
    operation: str,
    backend: str,
    result: dict[str, object],
    *,
    verbose: bool,
) -> dict[str, object]:
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": operation,
            "backend": backend,
            "result": result,
            "provenance": {
                "implementation": "dnakit.variants.external",
                "execution": "caller-installed executable; no shell or implicit download",
            },
        },
        verbose=verbose,
    )


def _validate_fraction(value: float | None, name: str) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be between 0 and 1.")
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ConfigurationError(f"{name} must be between 0 and 1.")


def _validate_number(
    value: float | None, name: str, *, minimum: float, maximum: float
) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    if not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")


def _validate_choice(value: str, name: str, choices: set[str]) -> str:
    if not isinstance(value, str) or value not in choices:
        raise ConfigurationError(f"{name} must be one of {sorted(choices)}.")
    return value


def _validate_bool(value: bool, name: str) -> None:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")


def _validate_integer(value: int | None, name: str, *, minimum: int, maximum: int) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")


def _validate_program(value: str, name: str) -> None:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > 4096
        or any(char in value for char in "\x00\r\n")
    ):
        raise ConfigurationError(f"{name} must be a bounded executable name or path.")


def _validate_bounded_integer(value: int, name: str, *, minimum: int, maximum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")


def _validate_simple_text(value: str, name: str, *, maximum: int) -> None:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > maximum
        or any(char.isspace() for char in value)
        or any(char in value for char in "/\\\x00\r\n")
    ):
        raise ConfigurationError(f"{name} must be a simple non-empty text value.")


__all__ = [
    "annotate_vcf_ensembl_vep",
    "annotate_vcf_vcfanno",
    "call_cutesv",
    "call_deepvariant",
    "call_delly",
    "call_freebayes",
    "call_gatk_haplotypecaller",
    "call_sniffles2",
    "filter_vcf_vcftools",
    "load_gemini_database",
    "merge_structural_variants_survivor",
    "query_gemini_database",
]

_register_public_verbose(globals(), tuple(__all__))
