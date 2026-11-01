"""Picard quality-control, metrics and fingerprinting tools as typed DNAKit facades (integration
pack P08, part 1).

The user installs Picard (``conda install -c bioconda picard``, or ``picard.jar`` plus a Java
runtime). Each function runs one Picard tool with an explicit, typed ``KEY=value`` argument list (no
shell). Picard is found from ``picard_path`` (a ``picard`` launcher or a ``picard.jar``), then
``picard`` on PATH, then the ``PICARD_JAR`` environment variable; a missing program raises
``BackendUnavailableError`` with an install hint. Outputs are new explicit paths: an existing output
is never overwritten. Results are JSON-serialisable dictionaries with ``status`` (``completed`` only
when the exit code is accepted and every expected output exists), ``backend``, ``tool``,
``version``, ``command``, ``return_code``, ``outputs``, ``log_tail`` and, where the tool writes a
Picard metrics file, the parsed ``metrics``.

Optional arguments left as ``None`` are not passed, so the defaults are exactly Picard's own
(documented per function where they matter). Tested with Picard 3.5.0. Charts (``chart_output``)
need ``Rscript`` on PATH, as in Picard itself.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _picard_cli as _p
from dnakit.exceptions import ConfigurationError

_DATA_TYPES = ("READGROUP", "LIBRARY", "SAMPLE", "FILE")
_ACCUMULATION_LEVELS = ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP")
_STRAND_SPECIFICITY = (
    "NONE",
    "FIRST_READ_TRANSCRIPTION_STRAND",
    "SECOND_READ_TRANSCRIPTION_STRAND",
)
_ILLUMINA_DATA_TYPES = (
    "Position",
    "BaseCalls",
    "QualityScores",
    "PF",
    "Barcodes",
    "Locs",
    "Clocs",
    "Pos",
    "Filter",
    "CBCL",
)


def _prefix(value: str | Path, name: str, suffixes: Sequence[str]) -> Path:
    """An output prefix: its folder must exist and none of the files it will create may exist."""
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a path prefix.")
    base = Path(value).expanduser().resolve()
    if not base.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {base.parent}")
    for suffix in suffixes:
        if Path(f"{base}{suffix}").exists():
            raise ConfigurationError(f"{name} would overwrite an existing file: {base}{suffix}")
    return base


def _with(prefix: Path, suffixes: Sequence[str]) -> dict[str, Path]:
    return {suffix.lstrip("."): Path(f"{prefix}{suffix}") for suffix in suffixes}


def _sam_common(
    validation_stringency: str | None, reference_sequence: str | Path | None
) -> list[str]:
    reference = _p.optional_file(reference_sequence, "reference_sequence")
    return _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency)) + _p.kv(
        "REFERENCE_SEQUENCE", reference
    )


def _metrics_of(result: dict[str, Any], path: Path, key: str = "metrics") -> dict[str, Any]:
    if path.is_file():
        parsed = _p.read_metrics(path)
        result[key] = parsed["metrics"]
        result[f"{key}_class"] = parsed["metrics_class"]
        if parsed["histogram"]:
            result["histogram"] = parsed["histogram"]
    return result


# ---------------------------------------------------------------------------
# Merging and formatting metrics files
# ---------------------------------------------------------------------------


def picard_accumulate_metrics(
    input_paths: Sequence[str],
    output_path: str,
    *,
    metric_type: str = "quality_yield",
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge sharded Picard metric files into one (``AccumulateQualityYieldMetrics`` or
    ``AccumulateVariantCallingMetrics``).

    ``metric_type="quality_yield"`` merges ``QualityYieldMetrics`` files (``input_paths`` are the
    files, ``output_path`` is the merged file). ``metric_type="variant_calling"`` merges Variant
    Calling metrics shards: ``input_paths`` and ``output_path`` are file prefixes without the
    ``.variant_calling_detail_metrics`` / ``.variant_calling_summary_metrics`` extensions, and both
    merged files are written.
    """
    kind = _p.choice(metric_type, ("quality_yield", "variant_calling"), "metric_type")
    if kind == "quality_yield":
        sources = _p.existing_files(input_paths, "input_paths")
        target = _p.new_file(output_path, "output_path")
        args = [f"INPUT={path}" for path in sources] + [f"OUTPUT={target}"]
        result = _p.run(
            "AccumulateQualityYieldMetrics",
            args,
            outputs={"output": target},
            picard_path=picard_path,
            timeout_seconds=timeout_seconds,
        )
        return _metrics_of(result, target)
    if isinstance(input_paths, str) or not input_paths:
        raise ConfigurationError("input_paths must be a non-empty list of metric file prefixes.")
    prefixes = []
    for value in input_paths:
        base = Path(value).expanduser().resolve()
        _p.existing_file(f"{base}.variant_calling_detail_metrics", "input_paths")
        _p.existing_file(f"{base}.variant_calling_summary_metrics", "input_paths")
        prefixes.append(base)
    target_prefix = _prefix(
        output_path,
        "output_path",
        (".variant_calling_detail_metrics", ".variant_calling_summary_metrics"),
    )
    outputs = _with(
        target_prefix, (".variant_calling_detail_metrics", ".variant_calling_summary_metrics")
    )
    args = [f"INPUT={path}" for path in prefixes] + [f"OUTPUT={target_prefix}"]
    result = _p.run(
        "AccumulateVariantCallingMetrics",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    _metrics_of(result, outputs["variant_calling_summary_metrics"], "summary_metrics")
    _metrics_of(result, outputs["variant_calling_detail_metrics"], "detail_metrics")
    return result


def picard_format_contamination_metrics(
    input_path: str,
    output_prefix: str,
    *,
    source: str = "bafregress",
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format a contamination estimator's output as a Picard metrics file.

    ``source="bafregress"`` runs ``CreateBafRegressMetricsFile`` (writes
    ``<prefix>.bafregress_metrics``); ``source="verifyidintensity"`` runs
    ``CreateVerifyIDIntensityContaminationMetricsFile`` (writes
    ``<prefix>.verifyidintensity_metrics``). ``input_path`` is the estimator's captured standard
    output.
    """
    kind = _p.choice(source, ("bafregress", "verifyidintensity"), "source")
    tool = {
        "bafregress": "CreateBafRegressMetricsFile",
        "verifyidintensity": "CreateVerifyIDIntensityContaminationMetricsFile",
    }[kind]
    suffix = f".{kind}_metrics"
    source_file = _p.existing_file(input_path, "input_path")
    base = _prefix(output_prefix, "output_prefix", (suffix,))
    outputs = _with(base, (suffix,))
    result = _p.run(
        tool,
        [f"INPUT={source_file}", f"OUTPUT={base}"],
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, outputs[f"{kind}_metrics"])


def picard_convert_artifact_to_oxog(
    input_base: str,
    output_base: str,
    *,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert ``CollectSequencingArtifactMetrics`` output into ``CollectOxoGMetrics`` format
    (``ConvertSequencingArtifactToOxoG``).

    ``input_base`` is the prefix used for the artifact metrics
    (``<base>.pre_adapter_detail_metrics`` and ``<base>.bait_bias_detail_metrics`` must exist);
    ``output_base`` is the prefix of the new ``<base>.oxog_metrics`` file.
    """
    base = Path(input_base).expanduser().resolve()
    _p.existing_file(f"{base}.pre_adapter_detail_metrics", "input_base")
    _p.existing_file(f"{base}.bait_bias_detail_metrics", "input_base")
    target = _prefix(output_base, "output_base", (".oxog_metrics",))
    outputs = _with(target, (".oxog_metrics",))
    result = _p.run(
        "ConvertSequencingArtifactToOxoG",
        [f"INPUT_BASE={base}", f"OUTPUT_BASE={target}"],
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, outputs["oxog_metrics"])


# ---------------------------------------------------------------------------
# Alignment-based metrics
# ---------------------------------------------------------------------------


def picard_base_distribution_by_cycle(
    input_path: str,
    output_path: str,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    chart_output: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    reference_sequence: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide distribution per read cycle of a SAM/BAM/CRAM file
    (``CollectBaseDistributionByCycle``).

    Picard's defaults apply (all reads, ``assume_sorted`` true). ``chart_output`` is an optional new
    ``.pdf`` and needs R.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    chart = _p.optional_new_file(chart_output, "chart_output")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.kv(
        "ALIGNED_READS_ONLY",
        None
        if aligned_reads_only is None
        else _p.boolean(aligned_reads_only, "aligned_reads_only"),
    )
    args += _p.kv(
        "PF_READS_ONLY",
        None if pf_reads_only is None else _p.boolean(pf_reads_only, "pf_reads_only"),
    )
    args += _p.kv("CHART_OUTPUT", chart)
    args += _p.kv(
        "ASSUME_SORTED",
        None if assume_sorted is None else _p.boolean(assume_sorted, "assume_sorted"),
    )
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"output": target, **({"chart": chart} if chart is not None else {})}
    result = _p.run(
        "CollectBaseDistributionByCycle",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        parsed = _p.read_metrics(target)
        result["cycles"] = parsed["metrics"]
        result["cycle_count"] = len(parsed["metrics"])
    return result


def picard_duplicate_metrics(
    input_path: str,
    metrics_file: str,
    *,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    reference_sequence: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Duplication metrics read from an already duplicate-marked SAM/BAM/CRAM file
    (``CollectDuplicateMetrics``).

    Picard notes that the optical-duplicate count is always 0 here, so ``ESTIMATED_LIBRARY_SIZE`` is
    slightly off.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(metrics_file, "metrics_file")
    args = [f"INPUT={source}", f"METRICS_FILE={target}"]
    args += _p.kv(
        "ASSUME_SORTED",
        None if assume_sorted is None else _p.boolean(assume_sorted, "assume_sorted"),
    )
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "CollectDuplicateMetrics",
        args,
        outputs={"metrics_file": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


def picard_oxog_metrics(
    input_path: str,
    output_path: str,
    reference_sequence: str,
    *,
    context_size: int | None = None,
    contexts: Sequence[str] | None = None,
    dbsnp: str | None = None,
    intervals: str | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_quality_score: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_non_pf_reads: bool | None = None,
    use_oq: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """8-oxoguanine (CpCG to CpCA) oxidation error metrics per sequence context
    (``CollectOxoGMetrics``).

    Picard's defaults apply (context size 1, mapping quality 30, base quality 20, insert size 60 to
    600).
    """
    source = _p.existing_file(input_path, "input_path")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}", f"REFERENCE_SEQUENCE={reference}"]
    args += _p.kv(
        "CONTEXT_SIZE",
        None if context_size is None else _p.integer(context_size, "context_size", minimum=0),
    )
    args += _p.kv_many("CONTEXTS", _p.texts(contexts, "contexts"))
    args += _p.kv("DB_SNP", _p.optional_file(dbsnp, "dbsnp"))
    args += _p.kv("INTERVALS", _p.optional_file(intervals, "intervals"))
    args += _p.kv(
        "MINIMUM_MAPPING_QUALITY",
        None
        if minimum_mapping_quality is None
        else _p.integer(minimum_mapping_quality, "minimum_mapping_quality", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_QUALITY_SCORE",
        None
        if minimum_quality_score is None
        else _p.integer(minimum_quality_score, "minimum_quality_score", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_INSERT_SIZE",
        None
        if minimum_insert_size is None
        else _p.integer(minimum_insert_size, "minimum_insert_size", minimum=0),
    )
    args += _p.kv(
        "MAXIMUM_INSERT_SIZE",
        None
        if maximum_insert_size is None
        else _p.integer(maximum_insert_size, "maximum_insert_size", minimum=0),
    )
    args += _p.kv(
        "INCLUDE_NON_PF_READS",
        None
        if include_non_pf_reads is None
        else _p.boolean(include_non_pf_reads, "include_non_pf_reads"),
    )
    args += _p.kv("USE_OQ", None if use_oq is None else _p.boolean(use_oq, "use_oq"))
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    result = _p.run(
        "CollectOxoGMetrics",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


def picard_snvq_yield_metrics(
    input_path: str,
    output_path: str,
    *,
    alternate_quality_attribute: str | None = None,
    include_bq_histogram: bool | None = None,
    include_secondary_alignments: bool | None = None,
    include_supplemental_alignments: bool | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    reference_sequence: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Yield of bases passing quality thresholds from SNV-quality reads
    (``CollectQualityYieldMetricsSNVQ``).

    The input reads must carry the per-base SNV quality tags ``qa``, ``qc``, ``qg`` and ``qt`` (as
    written by Ultima Genomics pipelines); without them Picard stops with an error and the result
    has ``status`` ``failed``.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.kv(
        "ALTERNATE_QUALITY_ATTRIBUTE",
        None
        if alternate_quality_attribute is None
        else _p.text(alternate_quality_attribute, "alternate_quality_attribute"),
    )
    args += _p.kv(
        "INCLUDE_BQ_HISTOGRAM",
        None
        if include_bq_histogram is None
        else _p.boolean(include_bq_histogram, "include_bq_histogram"),
    )
    args += _p.kv(
        "INCLUDE_SECONDARY_ALIGNMENTS",
        None
        if include_secondary_alignments is None
        else _p.boolean(include_secondary_alignments, "include_secondary_alignments"),
    )
    args += _p.kv(
        "INCLUDE_SUPPLEMENTAL_ALIGNMENTS",
        None
        if include_supplemental_alignments is None
        else _p.boolean(include_supplemental_alignments, "include_supplemental_alignments"),
    )
    args += _p.kv(
        "ASSUME_SORTED",
        None if assume_sorted is None else _p.boolean(assume_sorted, "assume_sorted"),
    )
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "CollectQualityYieldMetricsSNVQ",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


def picard_rnaseq_metrics(
    input_path: str,
    output_path: str,
    ref_flat: str,
    strand_specificity: str,
    *,
    ribosomal_intervals: str | None = None,
    chart_output: str | None = None,
    end_bias_bases: int | None = None,
    minimum_length: int | None = None,
    rrna_fragment_percentage: float | None = None,
    ignore_sequence: Sequence[str] | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    reference_sequence: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """RNA-seq coverage and transcript-region metrics from aligned reads (``CollectRnaSeqMetrics``).

    ``ref_flat`` is a UCSC refFlat gene annotation. ``strand_specificity`` is ``NONE``,
    ``FIRST_READ_TRANSCRIPTION_STRAND`` or ``SECOND_READ_TRANSCRIPTION_STRAND``.
    ``metric_accumulation_level`` takes ``ALL_READS`` (Picard's default), ``SAMPLE``, ``LIBRARY``
    and ``READ_GROUP``. ``chart_output`` needs R.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    annotation = _p.existing_file(ref_flat, "ref_flat")
    strand = _p.choice(strand_specificity, _STRAND_SPECIFICITY, "strand_specificity")
    chart = _p.optional_new_file(chart_output, "chart_output")
    args = [
        f"INPUT={source}",
        f"OUTPUT={target}",
        f"REF_FLAT={annotation}",
        f"STRAND_SPECIFICITY={strand}",
    ]
    args += _p.kv(
        "RIBOSOMAL_INTERVALS", _p.optional_file(ribosomal_intervals, "ribosomal_intervals")
    )
    args += _p.kv("CHART_OUTPUT", chart)
    args += _p.kv(
        "END_BIAS_BASES",
        None if end_bias_bases is None else _p.integer(end_bias_bases, "end_bias_bases", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_LENGTH",
        None if minimum_length is None else _p.integer(minimum_length, "minimum_length", minimum=0),
    )
    args += _p.kv(
        "RRNA_FRAGMENT_PERCENTAGE",
        None
        if rrna_fragment_percentage is None
        else _p.number(rrna_fragment_percentage, "rrna_fragment_percentage", minimum=0, maximum=1),
    )
    args += _p.kv_many("IGNORE_SEQUENCE", _p.texts(ignore_sequence, "ignore_sequence"))
    args += _p.kv_many(
        "METRIC_ACCUMULATION_LEVEL",
        _p.choices(metric_accumulation_level, _ACCUMULATION_LEVELS, "metric_accumulation_level"),
    )
    args += _p.kv(
        "ASSUME_SORTED",
        None if assume_sorted is None else _p.boolean(assume_sorted, "assume_sorted"),
    )
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"output": target, **({"chart": chart} if chart is not None else {})}
    result = _p.run(
        "CollectRnaSeqMetrics",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


def picard_sam_error_metrics(
    input_path: str,
    output_prefix: str,
    reference_sequence: str,
    *,
    error_metrics: Sequence[str] | None = None,
    vcf: str | None = None,
    intervals: Sequence[str] | None = None,
    interval_iterator: bool | None = None,
    min_base_q: int | None = None,
    min_mapping_q: int | None = None,
    prior_q: int | None = None,
    long_homopolymer: int | None = None,
    location_bin_size: int | None = None,
    max_loci: int | None = None,
    probability: float | None = None,
    file_extension: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Base-error metrics stratified by context, cycle, quality and other covariates
    (``CollectSamErrorMetrics``).

    One file ``<output_prefix>.<error>_by_<stratifier>`` is written per ``ERROR[:STRATIFIER...]``
    request; they are listed in ``files``. ``error_metrics`` entries are added to Picard's default
    list (``ERROR``, ``ERROR:BASE_QUALITY``, ``ERROR:INSERT_LENGTH``, ... ``INDEL_ERROR``), so
    repeating a default entry is rejected by Picard.
    """
    source = _p.existing_file(input_path, "input_path")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    base = _prefix(output_prefix, "output_prefix", (".error_by_all",))
    args = [f"INPUT={source}", f"OUTPUT={base}", f"REFERENCE_SEQUENCE={reference}"]
    args += _p.kv_many("ERROR_METRICS", _p.texts(error_metrics, "error_metrics"))
    args += _p.kv("VCF", _p.optional_file(vcf, "vcf"))
    args += _p.kv_many(
        "INTERVALS", [_p.existing_file(path, "intervals") for path in (intervals or [])]
    )
    args += _p.kv(
        "INTERVAL_ITERATOR",
        None if interval_iterator is None else _p.boolean(interval_iterator, "interval_iterator"),
    )
    args += _p.kv(
        "MIN_BASE_Q",
        None if min_base_q is None else _p.integer(min_base_q, "min_base_q", minimum=0),
    )
    args += _p.kv(
        "MIN_MAPPING_Q",
        None if min_mapping_q is None else _p.integer(min_mapping_q, "min_mapping_q", minimum=0),
    )
    args += _p.kv("PRIOR_Q", None if prior_q is None else _p.integer(prior_q, "prior_q", minimum=0))
    args += _p.kv(
        "LONG_HOMOPOLYMER",
        None
        if long_homopolymer is None
        else _p.integer(long_homopolymer, "long_homopolymer", minimum=1),
    )
    args += _p.kv(
        "LOCATION_BIN_SIZE",
        None
        if location_bin_size is None
        else _p.integer(location_bin_size, "location_bin_size", minimum=1),
    )
    args += _p.kv(
        "MAX_LOCI", None if max_loci is None else _p.integer(max_loci, "max_loci", minimum=0)
    )
    args += _p.kv(
        "PROBABILITY",
        None
        if probability is None
        else _p.number(probability, "probability", minimum=0, maximum=1),
    )
    args += _p.kv(
        "FILE_EXTENSION",
        None if file_extension is None else _p.text(file_extension, "file_extension"),
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    required = {"error_by_all": Path(f"{base}.error_by_all")}
    result = _p.run(
        "CollectSamErrorMetrics",
        args,
        outputs=required,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    files = sorted(path for path in base.parent.glob(f"{base.name}.*") if path.is_file())
    result["files"] = [str(path) for path in files]
    return result


def picard_sequencing_artifact_metrics(
    input_path: str,
    output_prefix: str,
    reference_sequence: str,
    *,
    context_size: int | None = None,
    contexts_to_print: Sequence[str] | None = None,
    dbsnp: str | None = None,
    intervals: str | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_quality_score: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_duplicates: bool | None = None,
    include_non_pf_reads: bool | None = None,
    include_unpaired: bool | None = None,
    tandem_reads: bool | None = None,
    use_oq: bool | None = None,
    assume_sorted: bool | None = None,
    file_extension: str | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Library-preparation (pre-adapter) and bait-bias substitution artifacts
    (``CollectSequencingArtifactMetrics``).

    Writes ``<prefix>.pre_adapter_detail_metrics``, ``.pre_adapter_summary_metrics``,
    ``.bait_bias_detail_metrics``, ``.bait_bias_summary_metrics`` and ``.error_summary_metrics``.
    Picard's defaults apply (context size 1, mapping quality 30, base quality 20, insert size 60 to
    600, duplicates and non-PF reads excluded).
    """
    source = _p.existing_file(input_path, "input_path")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    suffixes: tuple[str, ...] = (
        ".pre_adapter_detail_metrics",
        ".pre_adapter_summary_metrics",
        ".bait_bias_detail_metrics",
        ".bait_bias_summary_metrics",
        ".error_summary_metrics",
    )
    extension = None if file_extension is None else _p.text(file_extension, "file_extension")
    if extension is not None:
        suffixes = tuple(f"{suffix}{extension}" for suffix in suffixes)
    base = _prefix(output_prefix, "output_prefix", suffixes)
    args = [f"INPUT={source}", f"OUTPUT={base}", f"REFERENCE_SEQUENCE={reference}"]
    args += _p.kv(
        "CONTEXT_SIZE",
        None if context_size is None else _p.integer(context_size, "context_size", minimum=0),
    )
    args += _p.kv_many("CONTEXTS_TO_PRINT", _p.texts(contexts_to_print, "contexts_to_print"))
    args += _p.kv("DB_SNP", _p.optional_file(dbsnp, "dbsnp"))
    args += _p.kv("INTERVALS", _p.optional_file(intervals, "intervals"))
    args += _p.kv(
        "MINIMUM_MAPPING_QUALITY",
        None
        if minimum_mapping_quality is None
        else _p.integer(minimum_mapping_quality, "minimum_mapping_quality", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_QUALITY_SCORE",
        None
        if minimum_quality_score is None
        else _p.integer(minimum_quality_score, "minimum_quality_score", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_INSERT_SIZE",
        None
        if minimum_insert_size is None
        else _p.integer(minimum_insert_size, "minimum_insert_size", minimum=0),
    )
    args += _p.kv(
        "MAXIMUM_INSERT_SIZE",
        None
        if maximum_insert_size is None
        else _p.integer(maximum_insert_size, "maximum_insert_size", minimum=0),
    )
    args += _p.kv(
        "INCLUDE_DUPLICATES",
        None
        if include_duplicates is None
        else _p.boolean(include_duplicates, "include_duplicates"),
    )
    args += _p.kv(
        "INCLUDE_NON_PF_READS",
        None
        if include_non_pf_reads is None
        else _p.boolean(include_non_pf_reads, "include_non_pf_reads"),
    )
    args += _p.kv(
        "INCLUDE_UNPAIRED",
        None if include_unpaired is None else _p.boolean(include_unpaired, "include_unpaired"),
    )
    args += _p.kv(
        "TANDEM_READS", None if tandem_reads is None else _p.boolean(tandem_reads, "tandem_reads")
    )
    args += _p.kv("USE_OQ", None if use_oq is None else _p.boolean(use_oq, "use_oq"))
    args += _p.kv(
        "ASSUME_SORTED",
        None if assume_sorted is None else _p.boolean(assume_sorted, "assume_sorted"),
    )
    args += _p.kv("FILE_EXTENSION", extension)
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=0),
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    outputs = _with(base, suffixes)
    result = _p.run(
        "CollectSequencingArtifactMetrics",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    summary_key = next(key for key in outputs if key.startswith("pre_adapter_summary_metrics"))
    _metrics_of(result, outputs[summary_key], "pre_adapter_summary")
    return result


def picard_wgs_metrics_nonzero_coverage(
    input_path: str,
    output_path: str,
    reference_sequence: str,
    *,
    chart_output: str | None = None,
    theoretical_sensitivity_output: str | None = None,
    allele_fraction: Sequence[float] | None = None,
    coverage_cap: int | None = None,
    count_unpaired: bool | None = None,
    include_bq_histogram: bool | None = None,
    intervals: str | None = None,
    locus_accumulation_cap: int | None = None,
    minimum_base_quality: int | None = None,
    minimum_mapping_quality: int | None = None,
    read_length: int | None = None,
    sample_size: int | None = None,
    use_fast_algorithm: bool | None = None,
    stop_after: int | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Whole-genome coverage metrics, also restricted to bases with non-zero coverage
    (``CollectWgsMetricsWithNonZeroCoverage``).

    Picard marks this tool EXPERIMENTAL. Defaults: coverage cap 250, mapping quality 20, base
    quality 20, read length 150, sample size 10000. ``chart_output`` needs R. The ``metrics`` rows
    are ``WHOLE_GENOME`` and ``NON_ZERO_REGIONS``.
    """
    source = _p.existing_file(input_path, "input_path")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    target = _p.new_file(output_path, "output_path")
    chart = _p.optional_new_file(chart_output, "chart_output")
    sensitivity = _p.optional_new_file(
        theoretical_sensitivity_output, "theoretical_sensitivity_output"
    )
    args = [f"INPUT={source}", f"OUTPUT={target}", f"REFERENCE_SEQUENCE={reference}"]
    args += _p.kv("CHART_OUTPUT", chart)
    args += _p.kv("THEORETICAL_SENSITIVITY_OUTPUT", sensitivity)
    args += _p.kv_many(
        "ALLELE_FRACTION",
        [
            _p.number(value, "allele_fraction", minimum=0, maximum=1)
            for value in (allele_fraction or [])
        ],
    )
    args += _p.kv(
        "COVERAGE_CAP",
        None if coverage_cap is None else _p.integer(coverage_cap, "coverage_cap", minimum=1),
    )
    args += _p.kv(
        "COUNT_UNPAIRED",
        None if count_unpaired is None else _p.boolean(count_unpaired, "count_unpaired"),
    )
    args += _p.kv(
        "INCLUDE_BQ_HISTOGRAM",
        None
        if include_bq_histogram is None
        else _p.boolean(include_bq_histogram, "include_bq_histogram"),
    )
    args += _p.kv("INTERVALS", _p.optional_file(intervals, "intervals"))
    args += _p.kv(
        "LOCUS_ACCUMULATION_CAP",
        None
        if locus_accumulation_cap is None
        else _p.integer(locus_accumulation_cap, "locus_accumulation_cap", minimum=1),
    )
    args += _p.kv(
        "MINIMUM_BASE_QUALITY",
        None
        if minimum_base_quality is None
        else _p.integer(minimum_base_quality, "minimum_base_quality", minimum=0),
    )
    args += _p.kv(
        "MINIMUM_MAPPING_QUALITY",
        None
        if minimum_mapping_quality is None
        else _p.integer(minimum_mapping_quality, "minimum_mapping_quality", minimum=0),
    )
    args += _p.kv(
        "READ_LENGTH",
        None if read_length is None else _p.integer(read_length, "read_length", minimum=1),
    )
    args += _p.kv(
        "SAMPLE_SIZE",
        None if sample_size is None else _p.integer(sample_size, "sample_size", minimum=1),
    )
    args += _p.kv(
        "USE_FAST_ALGORITHM",
        None
        if use_fast_algorithm is None
        else _p.boolean(use_fast_algorithm, "use_fast_algorithm"),
    )
    args += _p.kv(
        "STOP_AFTER",
        None if stop_after is None else _p.integer(stop_after, "stop_after", minimum=-1),
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    outputs = {"output": target}
    if chart is not None:
        outputs["chart"] = chart
    if sensitivity is not None:
        outputs["theoretical_sensitivity_output"] = sensitivity
    result = _p.run(
        "CollectWgsMetricsWithNonZeroCoverage",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


# ---------------------------------------------------------------------------
# Variant-calling metrics
# ---------------------------------------------------------------------------


def picard_variant_calling_metrics(
    input_path: str,
    dbsnp: str,
    output_prefix: str,
    *,
    gvcf_input: bool | None = None,
    sequence_dictionary: str | None = None,
    target_intervals: str | None = None,
    thread_count: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-sample and aggregate variant-calling metrics of a VCF (``CollectVariantCallingMetrics``).

    ``dbsnp`` is a reference dbSNP file (dbSNP or VCF format, indexed). Writes
    ``<prefix>.variant_calling_detail_metrics`` and ``<prefix>.variant_calling_summary_metrics``.
    """
    source = _p.existing_file(input_path, "input_path")
    known = _p.existing_file(dbsnp, "dbsnp")
    suffixes = (".variant_calling_detail_metrics", ".variant_calling_summary_metrics")
    base = _prefix(output_prefix, "output_prefix", suffixes)
    args = [f"INPUT={source}", f"DBSNP={known}", f"OUTPUT={base}"]
    args += _p.kv(
        "GVCF_INPUT", None if gvcf_input is None else _p.boolean(gvcf_input, "gvcf_input")
    )
    args += _p.kv(
        "SEQUENCE_DICTIONARY", _p.optional_file(sequence_dictionary, "sequence_dictionary")
    )
    args += _p.kv("TARGET_INTERVALS", _p.optional_file(target_intervals, "target_intervals"))
    args += _p.kv(
        "THREAD_COUNT",
        None if thread_count is None else _p.integer(thread_count, "thread_count", minimum=1),
    )
    outputs = _with(base, suffixes)
    result = _p.run(
        "CollectVariantCallingMetrics",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    _metrics_of(result, outputs["variant_calling_summary_metrics"], "summary_metrics")
    _metrics_of(result, outputs["variant_calling_detail_metrics"], "detail_metrics")
    return result


def picard_arrays_variant_calling_metrics(
    input_path: str,
    dbsnp: str,
    output_prefix: str,
    *,
    call_rate_pf_threshold: float | None = None,
    sequence_dictionary: str | None = None,
    num_processors: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summary and per-sample metrics of a genotyping-array VCF written by GtcToVcf
    (``CollectArraysVariantCallingMetrics``).

    Writes ``<prefix>.arrays_variant_calling_detail_metrics``,
    ``.arrays_variant_calling_summary_metrics`` and ``.arrays_control_code_summary_metrics``. The
    call-rate threshold for an autocall pass defaults to 0.98.
    """
    source = _p.existing_file(input_path, "input_path")
    known = _p.existing_file(dbsnp, "dbsnp")
    suffixes = (
        ".arrays_variant_calling_detail_metrics",
        ".arrays_variant_calling_summary_metrics",
        ".arrays_control_code_summary_metrics",
    )
    base = _prefix(output_prefix, "output_prefix", suffixes)
    args = [f"INPUT={source}", f"DBSNP={known}", f"OUTPUT={base}"]
    args += _p.kv(
        "CALL_RATE_PF_THRESHOLD",
        None
        if call_rate_pf_threshold is None
        else _p.number(call_rate_pf_threshold, "call_rate_pf_threshold", minimum=0, maximum=1),
    )
    args += _p.kv(
        "SEQUENCE_DICTIONARY", _p.optional_file(sequence_dictionary, "sequence_dictionary")
    )
    args += _p.kv(
        "NUM_PROCESSORS",
        None if num_processors is None else _p.integer(num_processors, "num_processors"),
    )
    outputs = _with(base, suffixes)
    result = _p.run(
        "CollectArraysVariantCallingMetrics",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        required=[
            "arrays_variant_calling_detail_metrics",
            "arrays_variant_calling_summary_metrics",
        ],
    )
    _metrics_of(result, outputs["arrays_variant_calling_summary_metrics"], "summary_metrics")
    _metrics_of(result, outputs["arrays_variant_calling_detail_metrics"], "detail_metrics")
    return result


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def picard_read_group_checksum(
    input_path: str,
    output_path: str,
    *,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hash code of the read groups in a SAM/BAM header (``CalculateReadGroupChecksum``).

    Adding or removing read groups changes the code. The hash is written to ``output_path`` and
    returned as ``checksum``.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    result = _p.run(
        "CalculateReadGroupChecksum",
        [f"INPUT={source}", f"OUTPUT={target}"],
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["checksum"] = target.read_text(encoding="utf-8").strip()
    return result


def picard_check_terminator_block(
    input_path: str,
    *,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the last BGZF block of a block-compressed file such as a BAM is the well-formed
    terminator (``CheckTerminatorBlock``).

    Picard exits with code 100 when the terminator is missing or malformed; that is a verdict, not a
    failure, so ``status`` is ``completed`` and ``valid`` is False.
    """
    source = _p.existing_file(input_path, "input_path")
    result = _p.run(
        "CheckTerminatorBlock",
        [f"INPUT={source}"],
        outputs={},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        ok_codes=(0, 100),
        cwd=source.parent,
    )
    result["valid"] = result["return_code"] == 0
    return result


def picard_check_illumina_directory(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    *,
    data_types: Sequence[str] | None = None,
    tile_numbers: Sequence[int] | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check an Illumina BaseCalls directory for the files every tile and cycle needs
    (``CheckIlluminaDirectory``).

    ``read_structure`` is the integer/letter description of a run (``B`` sample barcode, ``M``
    molecular barcode, ``T`` template, ``S`` skip), for example ``76T8B76T``. ``valid`` is true when
    Picard finds no missing or empty file; otherwise ``status`` is ``failed`` and ``log_tail`` lists
    what is missing. Files are never created (Picard's ``FAKE_FILES`` and ``LINK_LOCS`` are not
    exposed).
    """
    directory = _p.existing_directory(basecalls_dir, "basecalls_dir")
    if isinstance(lanes, (str, bytes)) or not lanes:
        raise ConfigurationError("lanes must be a non-empty list of lane numbers.")
    args = [f"BASECALLS_DIR={directory}"]
    args += _p.kv_many("LANES", [_p.integer(lane, "lanes", minimum=1) for lane in lanes])
    args += _p.kv("READ_STRUCTURE", _p.text(read_structure, "read_structure"))
    args += _p.kv_many("DATA_TYPES", _p.choices(data_types, _ILLUMINA_DATA_TYPES, "data_types"))
    args += _p.kv_many(
        "TILE_NUMBERS",
        [_p.integer(tile, "tile_numbers", minimum=1) for tile in (tile_numbers or [])],
    )
    result = _p.run(
        "CheckIlluminaDirectory",
        args,
        outputs={},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        cwd=directory,
    )
    result["valid"] = result["status"] == "completed"
    return result


# ---------------------------------------------------------------------------
# Fingerprinting
# ---------------------------------------------------------------------------


def picard_fingerprint_metrics(
    input_paths: Sequence[str],
    haplotype_map: str,
    output_path: str,
    *,
    calculate_by: str | None = None,
    genotype_lod_threshold: float | None = None,
    number_of_sampling: int | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Viability statistics of single fingerprints (``CalculateFingerprintMetrics``).

    Reports p-values under the null that each fingerprint is a non-contaminated diploid human, plus
    the discriminatory power. ``calculate_by`` is ``READGROUP`` (Picard's default), ``LIBRARY``,
    ``SAMPLE`` or ``FILE``.
    """
    sources = _p.existing_files(input_paths, "input_paths")
    haplotypes = _p.existing_file(haplotype_map, "haplotype_map")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={path}" for path in sources] + [
        f"HAPLOTYPE_MAP={haplotypes}",
        f"OUTPUT={target}",
    ]
    args += _p.kv(
        "CALCULATE_BY",
        None if calculate_by is None else _p.choice(calculate_by, _DATA_TYPES, "calculate_by"),
    )
    args += _p.kv(
        "GENOTYPE_LOD_THRESHOLD",
        None
        if genotype_lod_threshold is None
        else _p.number(genotype_lod_threshold, "genotype_lod_threshold"),
    )
    args += _p.kv(
        "NUMBER_OF_SAMPLING",
        None
        if number_of_sampling is None
        else _p.integer(number_of_sampling, "number_of_sampling", minimum=1),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "CalculateFingerprintMetrics",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


def picard_extract_fingerprint(
    input_path: str,
    haplotype_map: str,
    output_path: str,
    reference_sequence: str,
    *,
    sample_alias: str | None = None,
    contamination: float | None = None,
    extract_contamination: bool | None = None,
    extract_non_representatives_too: bool | None = None,
    locus_max_reads: int | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the fingerprint genotype likelihoods of a file into a VCF (``ExtractFingerprint``).

    The likelihoods are the PL values at the fingerprint sites.

    ``reference_sequence`` is required by Picard for this tool."""
    source = _p.existing_file(input_path, "input_path")
    haplotypes = _p.existing_file(haplotype_map, "haplotype_map")
    target = _p.new_file(output_path, "output_path")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    args = [
        f"INPUT={source}",
        f"HAPLOTYPE_MAP={haplotypes}",
        f"OUTPUT={target}",
        f"REFERENCE_SEQUENCE={reference}",
    ]
    args += _p.kv(
        "SAMPLE_ALIAS", None if sample_alias is None else _p.text(sample_alias, "sample_alias")
    )
    args += _p.kv(
        "CONTAMINATION",
        None
        if contamination is None
        else _p.number(contamination, "contamination", minimum=0, maximum=1),
    )
    args += _p.kv(
        "EXTRACT_CONTAMINATION",
        None
        if extract_contamination is None
        else _p.boolean(extract_contamination, "extract_contamination"),
    )
    args += _p.kv(
        "EXTRACT_NON_REPRESENTATIVES_TOO",
        None
        if extract_non_representatives_too is None
        else _p.boolean(extract_non_representatives_too, "extract_non_representatives_too"),
    )
    args += _p.kv(
        "LOCUS_MAX_READS",
        None
        if locus_max_reads is None
        else _p.integer(locus_max_reads, "locus_max_reads", minimum=1),
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    result = _p.run(
        "ExtractFingerprint",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    if target.is_file():
        result["sites"] = _p.count_lines(target, skip_prefix="#")
    return result


def picard_check_fingerprint(
    input_path: str,
    genotypes: str,
    haplotype_map: str,
    output_prefix: str,
    *,
    expected_sample_alias: str | None = None,
    observed_sample_alias: str | None = None,
    genotype_lod_threshold: float | None = None,
    ignore_read_groups: bool | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that a SAM/BAM/CRAM or VCF matches the expected sample's known genotypes
    (``CheckFingerprint``).

    ``genotypes`` is a VCF of known genotypes; ``haplotype_map`` lists the fingerprint SNPs. Writes
    ``<prefix>.fingerprinting_summary_metrics`` and ``<prefix>.fingerprinting_detail_metrics``. Exit
    code 1 (expected sample not found in the genotypes) and 2 (no valid checks) are reported as
    ``failed``.
    """
    source = _p.existing_file(input_path, "input_path")
    known = _p.existing_file(genotypes, "genotypes")
    haplotypes = _p.existing_file(haplotype_map, "haplotype_map")
    suffixes = (".fingerprinting_summary_metrics", ".fingerprinting_detail_metrics")
    base = _prefix(output_prefix, "output_prefix", suffixes)
    args = [
        f"INPUT={source}",
        f"GENOTYPES={known}",
        f"HAPLOTYPE_MAP={haplotypes}",
        f"OUTPUT={base}",
    ]
    args += _p.kv(
        "EXPECTED_SAMPLE_ALIAS",
        None
        if expected_sample_alias is None
        else _p.text(expected_sample_alias, "expected_sample_alias"),
    )
    args += _p.kv(
        "OBSERVED_SAMPLE_ALIAS",
        None
        if observed_sample_alias is None
        else _p.text(observed_sample_alias, "observed_sample_alias"),
    )
    args += _p.kv(
        "GENOTYPE_LOD_THRESHOLD",
        None
        if genotype_lod_threshold is None
        else _p.number(genotype_lod_threshold, "genotype_lod_threshold"),
    )
    args += _p.kv(
        "IGNORE_READ_GROUPS",
        None
        if ignore_read_groups is None
        else _p.boolean(ignore_read_groups, "ignore_read_groups"),
    )
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = _with(base, suffixes)
    result = _p.run(
        "CheckFingerprint",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    _metrics_of(result, outputs["fingerprinting_summary_metrics"], "summary_metrics")
    _metrics_of(result, outputs["fingerprinting_detail_metrics"], "detail_metrics")
    return result


def picard_crosscheck_fingerprints(
    input_paths: Sequence[str],
    haplotype_map: str,
    output_path: str,
    *,
    crosscheck_by: str | None = None,
    second_input_paths: Sequence[str] | None = None,
    crosscheck_mode: str | None = None,
    lod_threshold: float | None = None,
    expect_all_groups_to_match: bool | None = None,
    output_errors_only: bool | None = None,
    allow_duplicate_reads: bool | None = None,
    calculate_tumor_aware_results: bool | None = None,
    matrix_output: str | None = None,
    num_threads: int | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that groups of reads or variants come from the same individual
    (``CrosscheckFingerprints``).

    ``crosscheck_by`` rolls fingerprints up to ``READGROUP`` (Picard's default), ``LIBRARY``,
    ``SAMPLE`` or ``FILE``. ``second_input_paths`` compares each sample of the first set with the
    same sample of the second set (``crosscheck_mode`` is ``CHECK_SAME_SAMPLE`` by default, or
    ``CHECK_ALL_OTHERS``). Picard exits with code 1 when a group relates to the others unexpectedly;
    the metrics file is still written, so ``status`` stays ``completed`` and ``mismatch_detected``
    is true. Exit code 2 (no valid checks) is a failure.
    """
    sources = _p.existing_files(input_paths, "input_paths")
    haplotypes = _p.existing_file(haplotype_map, "haplotype_map")
    target = _p.new_file(output_path, "output_path")
    matrix = _p.optional_new_file(matrix_output, "matrix_output")
    second = (
        _p.existing_files(second_input_paths, "second_input_paths") if second_input_paths else []
    )
    args = [f"INPUT={path}" for path in sources] + [
        f"HAPLOTYPE_MAP={haplotypes}",
        f"OUTPUT={target}",
    ]
    args += [f"SECOND_INPUT={path}" for path in second]
    args += _p.kv(
        "CROSSCHECK_BY",
        None if crosscheck_by is None else _p.choice(crosscheck_by, _DATA_TYPES, "crosscheck_by"),
    )
    args += _p.kv(
        "CROSSCHECK_MODE",
        None
        if crosscheck_mode is None
        else _p.choice(
            crosscheck_mode, ("CHECK_SAME_SAMPLE", "CHECK_ALL_OTHERS"), "crosscheck_mode"
        ),
    )
    args += _p.kv(
        "LOD_THRESHOLD",
        None if lod_threshold is None else _p.number(lod_threshold, "lod_threshold"),
    )
    args += _p.kv(
        "EXPECT_ALL_GROUPS_TO_MATCH",
        None
        if expect_all_groups_to_match is None
        else _p.boolean(expect_all_groups_to_match, "expect_all_groups_to_match"),
    )
    args += _p.kv(
        "OUTPUT_ERRORS_ONLY",
        None
        if output_errors_only is None
        else _p.boolean(output_errors_only, "output_errors_only"),
    )
    args += _p.kv(
        "ALLOW_DUPLICATE_READS",
        None
        if allow_duplicate_reads is None
        else _p.boolean(allow_duplicate_reads, "allow_duplicate_reads"),
    )
    args += _p.kv(
        "CALCULATE_TUMOR_AWARE_RESULTS",
        None
        if calculate_tumor_aware_results is None
        else _p.boolean(calculate_tumor_aware_results, "calculate_tumor_aware_results"),
    )
    args += _p.kv("MATRIX_OUTPUT", matrix)
    args += _p.kv(
        "NUM_THREADS",
        None if num_threads is None else _p.integer(num_threads, "num_threads", minimum=1),
    )
    args += ["EXIT_CODE_WHEN_MISMATCH=1", "EXIT_CODE_WHEN_NO_VALID_CHECKS=2"]
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"output": target, **({"matrix": matrix} if matrix is not None else {})}
    result = _p.run(
        "CrosscheckFingerprints",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        ok_codes=(0, 1),
    )
    result["mismatch_detected"] = result["return_code"] == 1
    return _metrics_of(result, target)


def picard_cluster_crosscheck_metrics(
    input_path: str,
    output_path: str,
    *,
    lod_threshold: float | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster the groups of a ``CrosscheckFingerprints`` metrics file by LOD score
    (``ClusterCrosscheckMetrics``).

    Two groups are connected when their LOD score exceeds ``lod_threshold`` (Picard's default is
    0.0); the resulting clusters are written to ``output_path`` and parsed in ``metrics``.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.kv(
        "LOD_THRESHOLD",
        None if lod_threshold is None else _p.number(lod_threshold, "lod_threshold"),
    )
    result = _p.run(
        "ClusterCrosscheckMetrics",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics_of(result, target)


__all__ = [
    "picard_accumulate_metrics",
    "picard_arrays_variant_calling_metrics",
    "picard_base_distribution_by_cycle",
    "picard_check_fingerprint",
    "picard_check_illumina_directory",
    "picard_check_terminator_block",
    "picard_cluster_crosscheck_metrics",
    "picard_convert_artifact_to_oxog",
    "picard_crosscheck_fingerprints",
    "picard_duplicate_metrics",
    "picard_extract_fingerprint",
    "picard_fingerprint_metrics",
    "picard_format_contamination_metrics",
    "picard_oxog_metrics",
    "picard_read_group_checksum",
    "picard_rnaseq_metrics",
    "picard_sam_error_metrics",
    "picard_sequencing_artifact_metrics",
    "picard_snvq_yield_metrics",
    "picard_variant_calling_metrics",
    "picard_wgs_metrics_nonzero_coverage",
]
_register_public_verbose(globals())
