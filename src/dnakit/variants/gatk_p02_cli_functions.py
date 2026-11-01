"""GATK command-line facades (pack P02): VCF filtering, recalibration, scoring and VCF utility
tools.

Each public function runs one GATK tool (4.7.0.0 is the tested version) through the shared facade
layer ``dnakit.bioinformatics._gap_cli`` (``dnakit.backends.external.execute_bounded_command``):
fixed argument whitelist per tool, no shell, the ``gatk`` launcher resolved from an explicit path,
then PATH. The user installs GATK (conda ``bioconda::gatk4`` or the release zip from
https://github.com/broadinstitute/gatk/releases) and Java 17; ``gatk_path`` points at the ``gatk``
launcher when it is not on PATH. Outputs are written to new explicit paths (an existing output is
refused, nothing is overwritten) and GATK's own index files (``.tbi``, ``.idx``) are reported next
to them.

Results are JSON-serialisable dictionaries: ``status`` (``completed`` only when the exit code is 0
and every expected output exists), ``backend``, ``version``, ``command``, ``return_code``,
``outputs``, a bounded ``log_tail``, the GATK ``tool`` name and tool-specific parsed values (record
counts, evaluation tables, validation verdict). Only arguments the caller sets are passed on, so
every unset argument keeps GATK's own default; booleans whose GATK default is true are Python
``True`` by default.

Tools that are Picard programs inside GATK (FixVcfHeader, LiftoverVcf, MakeSitesOnlyVcf, SplitVcfs,
VcfFormatConverter, VcfToIntervalList) and walkers (the other tools) behave as documented by the
GATK tool index. NVScoreVariants and the PYTHON_IFOREST backend of the variant-annotation-model
tools (the GATK default) start a Python process: activate the GATK conda environment
(``gatkcondaenv.yml`` in the release) or otherwise put a ``python`` with its packages first on PATH.
Feature files that GATK queries by interval (resources, dbSNP, comparison and mask files) must be
indexed (``gatk IndexFeatureFile``).
"""

from __future__ import annotations

import gzip
import os
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError

_PROGRAM = "gatk"
_TESTED_VERSION = "4.7.0.0"
_HINT = (
    "Install GATK 4.6 or later with Java 17 (conda install -c bioconda gatk4, or the release zip "
    "from https://github.com/broadinstitute/gatk/releases), then put the `gatk` launcher on PATH "
    "or pass gatk_path."
)
_VCF_SUFFIXES = (".vcf", ".vcf.gz", ".bcf")
_TAG = re.compile(r"[A-Za-z0-9_.+-]+")
_MODES = ("SNP", "INDEL", "BOTH")
_VARIANT_TYPES = ("SNP", "INDEL")
_MODEL_BACKENDS = ("JAVA_BGMM", "PYTHON_IFOREST", "PYTHON_SCRIPT")
_MATCHING_STRATEGIES = (
    "START_POSITION",
    "START_POSITION_AND_GIVEN_REPRESENTATION",
    "START_POSITION_AND_MINIMAL_REPRESENTATION",
)
_VALIDATION_TYPES = ("ALL", "REF", "IDS", "ALLELES", "CHR_COUNTS")
_ID_METHODS = ("CONCAT_ALL", "FIRST", "LAST", "NONE")
_THRESHOLD_STRATEGIES = ("CONSTANT", "FALSE_DISCOVERY_RATE", "OPTIMAL_F_SCORE")
_TENSOR_TYPES = ("reference", "read_tensor")
_REPORT_TABLE = re.compile(r"^#:GATKTable:(?P<name>[^:\s]+):(?P<description>.*)$")
_REPORT_FORMAT = re.compile(r"^#:GATKTable:\d+:\d+:")


# ---------------------------------------------------------------------------
# Argument builder and shared helpers
# ---------------------------------------------------------------------------


class _Args:
    """Collects the GATK argument vector; every method validates its value before adding it."""

    def __init__(self) -> None:
        self.items: list[str] = []

    def on(self, flag: str, value: bool, name: str) -> None:
        """A boolean whose GATK default is false: passed (as ``true``) only when set."""
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be True or False.")
        if value:
            self.items += [flag, "true"]

    def off(self, flag: str, value: bool, name: str) -> None:
        """A boolean whose GATK default is true: passed (as ``false``) only when switched off."""
        if not isinstance(value, bool):
            raise ConfigurationError(f"{name} must be True or False.")
        if not value:
            self.items += [flag, "false"]

    def text(self, flag: str, value: str | None, name: str) -> None:
        if value is not None:
            self.items += [flag, _cli.text(value, name)]

    def integer(
        self,
        flag: str,
        value: int | None,
        name: str,
        *,
        minimum: int | None = None,
        maximum: int | None = None,
    ) -> None:
        if value is not None:
            self.items += [flag, str(_cli.integer(value, name, minimum=minimum, maximum=maximum))]

    def number(
        self,
        flag: str,
        value: float | None,
        name: str,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
    ) -> None:
        if value is not None:
            self.items += [flag, repr(_cli.number(value, name, minimum=minimum, maximum=maximum))]

    def choice(self, flag: str, value: str | None, options: Sequence[str], name: str) -> None:
        if value is not None:
            self.items += [flag, _cli.choice(value, options, name)]

    def file(self, flag: str, value: str | Path | None, name: str) -> None:
        if value is not None:
            self.items += [flag, str(_cli.existing_file(value, name))]

    def many_text(self, flag: str, values: Sequence[str] | None, name: str) -> None:
        for value in _sequence(values, name):
            self.items += [flag, _cli.text(value, name)]

    def many_choice(
        self, flag: str, values: Sequence[str] | None, options: Sequence[str], name: str
    ) -> None:
        for value in _sequence(values, name):
            self.items += [flag, _cli.choice(value, options, name)]

    def many_file(self, flag: str, values: Sequence[str | Path] | None, name: str) -> None:
        for value in _sequence(values, name):
            self.items += [flag, str(_cli.existing_file(value, name))]

    def features(
        self, flag: str, values: Sequence[str | Path | Mapping[str, Any]] | None, name: str
    ) -> None:
        """Feature inputs (resources, comparison sets): a path, or ``{path, name, tags}`` for a
        tagged input."""
        for value in _sequence(values, name):
            if isinstance(value, Mapping):
                path = _cli.existing_file(value.get("path", ""), f"{name} path")
                label = _tag(value.get("name"), f"{name} name") if value.get("name") else None
                tags = value.get("tags") or {}
                if not isinstance(tags, Mapping):
                    raise ConfigurationError(f"{name} tags must be a mapping of tag name to value.")
                parts = [] if label is None else [label]
                parts += [
                    f"{_tag(key, f'{name} tag')}={_tag_value(tag)}" for key, tag in tags.items()
                ]
                self.items += [f"{flag}:{','.join(parts)}" if parts else flag, str(path)]
            else:
                self.items += [flag, str(_cli.existing_file(value, name))]

    def feature(self, flag: str, value: str | Path | Mapping[str, Any] | None, name: str) -> None:
        if value is not None:
            self.features(flag, [value], name)

    def intervals(
        self, intervals: Sequence[str] | None, exclude_intervals: Sequence[str] | None
    ) -> None:
        for flag, values, name in (
            ("-L", intervals, "intervals"),
            ("-XL", exclude_intervals, "exclude_intervals"),
        ):
            for value in _sequence(values, name):
                self.items += [flag, _interval(value, name)]


def _sequence(values: Sequence[Any] | None, name: str) -> list[Any]:
    if values is None:
        return []
    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    return list(values)


def _tag(value: object, name: str) -> str:
    if not isinstance(value, str) or _TAG.fullmatch(value) is None:
        raise ConfigurationError(f"{name} may contain only letters, digits and '_', '.', '+', '-'.")
    return value


def _tag_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    return _tag(value, "tag value")


def _interval(value: object, name: str) -> str:
    """An interval expression (``chr1:100-200``) or the path of an interval file (BED,
    interval_list, VCF)."""
    text = _cli.text(value, name)
    candidate = Path(text).expanduser()
    return str(candidate.resolve()) if candidate.is_file() else text


def _output_vcf(value: str | Path, name: str = "output_path") -> Path:
    path = _cli.new_file(value, name)
    if not str(path).endswith(_VCF_SUFFIXES):
        raise ConfigurationError(
            f"{name} must end with .vcf, .vcf.gz or .bcf (GATK picks the format from it)."
        )
    return path


def _input_vcf(value: str | Path, name: str = "vcf") -> Path:
    path = _cli.existing_file(value, name)
    if not str(path).endswith(_VCF_SUFFIXES):
        raise ConfigurationError(f"{name} must end with .vcf, .vcf.gz or .bcf.")
    return path


def _common_folder(paths: Sequence[Path]) -> Path:
    """The folder that contains every output (the executor requires monitored outputs to lie inside
    cwd)."""
    return Path(os.path.commonpath([str(path.parent) for path in paths]))


def _index_files(paths: Sequence[Path]) -> list[str]:
    found: list[str] = []
    for path in paths:
        for suffix in (".tbi", ".idx", ".csi"):
            candidate = Path(f"{path}{suffix}")
            if candidate.is_file():
                found.append(str(candidate))
    return found


def _gatk(
    tool: str,
    args: _Args,
    *,
    cwd: Path,
    expected: Sequence[Path],
    gatk_path: str | None,
    java_max_heap_gb: int | None,
    timeout_seconds: float,
    operation: str,
) -> dict[str, Any]:
    """Run ``gatk [--java-options -Xmx..] <tool> <args>`` and add the tool name, operation and index
    files."""
    heap: list[str] = []
    if java_max_heap_gb is not None:
        gigabytes = _cli.integer(java_max_heap_gb, "java_max_heap_gb", minimum=1, maximum=1024)
        heap = ["--java-options", f"-Xmx{gigabytes}g"]
    seconds = _cli.number(timeout_seconds, "timeout_seconds", minimum=1.0, maximum=86_400.0)
    result = _cli.run(
        _PROGRAM,
        [*heap, tool, *args.items],
        cwd=cwd,
        expected=list(expected),
        explicit=gatk_path,
        hint=_HINT,
        timeout_seconds=seconds,
    )
    result["tool"] = tool
    result["operation"] = operation
    result["tested_version"] = _TESTED_VERSION
    result["index_files"] = _index_files(expected)
    return result


def _vcf_summary(path: Path) -> dict[str, Any] | None:
    """Record count, FILTER value counts and samples of a VCF (plain or bgzipped); ``None`` for BCF
    or unreadable."""
    name = str(path)
    if name.endswith(".bcf") or not path.is_file():
        return None
    opener = gzip.open if name.endswith(".gz") else open
    records = 0
    filters: Counter[str] = Counter()
    samples: list[str] = []
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if line.startswith("#CHROM"):
                    samples = fields[9:]
                elif len(fields) >= 7:
                    records += 1
                    filters[fields[6]] += 1
    except (OSError, EOFError):
        return None
    return {"records": records, "filters": dict(sorted(filters.items())), "samples": samples}


def _with_summary(result: dict[str, Any], output: Path, key: str = "summary") -> dict[str, Any]:
    result["output"] = str(output) if output.is_file() else None
    result[key] = _vcf_summary(output) if result["status"] == "completed" else None
    return result


def _problem_lines(log: str) -> list[str]:
    return [
        line.strip()
        for line in log.splitlines()
        if re.search(r"\b(ERROR|WARN|Exception)\b|USER ERROR", line)
    ][:20]


# ---------------------------------------------------------------------------
# Validation and evaluation
# ---------------------------------------------------------------------------


def validate_vcf_gatk(
    vcf: str,
    *,
    reference: str | None = None,
    dbsnp: str | None = None,
    validate_gvcf: bool = False,
    fail_gvcf_on_overlap: bool = False,
    do_not_validate_filtered_records: bool = False,
    warn_on_errors: bool = False,
    exclude_validation_types: Sequence[str] | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a VCF with GATK's strict criteria (gatk ValidateVariants).

    Checks the reference allele against ``reference`` (REF), IDs against ``dbsnp`` (IDS), allele and
    chromosome counts (ALLELES, CHR_COUNTS); ``exclude_validation_types`` lists checks to skip (ALL,
    REF, IDS, ALLELES, CHR_COUNTS). ``validate_gvcf`` validates a GVCF and ``fail_gvcf_on_overlap``
    rejects overlapping GVCF positions. ``warn_on_errors`` only warns instead of stopping at the
    first error. Result keys: ``valid`` (exit code 0 and no error) and ``problems`` (warning and
    error lines of the log). The tool writes no output file; for an unindexed plain VCF GATK may
    create a ``.idx`` index next to the input.
    """
    path = _input_vcf(vcf)
    args = _Args()
    args.items += ["-V", str(path)]
    args.file("-R", reference, "reference")
    args.file("--dbsnp", dbsnp, "dbsnp")
    args.on("--validate-GVCF", validate_gvcf, "validate_gvcf")
    args.on("--fail-gvcf-on-overlap", fail_gvcf_on_overlap, "fail_gvcf_on_overlap")
    args.on(
        "--do-not-validate-filtered-records",
        do_not_validate_filtered_records,
        "do_not_validate_filtered_records",
    )
    args.on("--warn-on-errors", warn_on_errors, "warn_on_errors")
    args.many_choice(
        "--validation-type-to-exclude",
        exclude_validation_types,
        _VALIDATION_TYPES,
        "exclude_validation_types",
    )
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "ValidateVariants",
        args,
        cwd=path.parent,
        expected=[],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="validate_vcf_gatk",
    )
    result["valid"] = result["return_code"] == 0
    result["problems"] = _problem_lines(result["log_tail"])
    return result


def _parse_gatk_report(text: str) -> dict[str, dict[str, Any]]:
    """Parse a GATKReport (the output format of VariantEval) into ``{table: {description, columns,
    rows}}``."""
    tables: dict[str, dict[str, Any]] = {}
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        match = _REPORT_TABLE.match(line)
        if match is not None and not _REPORT_FORMAT.match(line):
            index += 1
            columns = lines[index].split() if index < len(lines) else []
            rows: list[list[object]] = []
            index += 1
            while index < len(lines) and lines[index].strip():
                rows.append([_cell(value) for value in lines[index].split()])
                index += 1
            tables[match.group("name")] = {
                "description": match.group("description"),
                "columns": columns,
                "rows": rows,
            }
        index += 1
    return tables


def _cell(value: str) -> object:
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def evaluate_variants_gatk(
    eval_vcfs: Sequence[str],
    reference: str,
    output_path: str,
    *,
    comparison_vcfs: Sequence[str] | None = None,
    dbsnp: str | None = None,
    gold_standard: str | None = None,
    eval_modules: Sequence[str] | None = None,
    stratification_modules: Sequence[str] | None = None,
    use_all_standard_modules: bool = True,
    use_all_standard_stratifications: bool = True,
    known_names: Sequence[str] | None = None,
    merge_evals: bool = False,
    require_strict_allele_match: bool = False,
    keep_ac0: bool = False,
    samples: Sequence[str] | None = None,
    sample_ploidy: int | None = None,
    pedigree: str | None = None,
    ancestral_alignments: str | None = None,
    known_cnvs: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """General-purpose variant evaluation: % in dbSNP, genotype concordance, Ti/Tv and more (gatk
    VariantEval).

    ``eval_vcfs`` are the call sets to evaluate, ``comparison_vcfs`` the sets to compare with and
    ``dbsnp`` / ``gold_standard`` the known-variation resources. ``eval_modules`` (-EV) and
    ``stratification_modules`` (-ST) add modules; switch ``use_all_standard_modules`` /
    ``use_all_standard_stratifications`` off to run only the modules listed (-no-ev / -no-st). The
    GATKReport written to ``output_path`` is parsed into ``tables`` (``{name: {description, columns,
    rows}}``). VariantEval is a beta GATK tool.
    """
    evals = _sequence(eval_vcfs, "eval_vcfs")
    if not evals:
        raise ConfigurationError("eval_vcfs must list at least one VCF.")
    output = _cli.new_file(output_path, "output_path")
    args = _Args()
    args.many_file("--eval", evals, "eval_vcfs")
    args.items += ["-R", str(_cli.existing_file(reference, "reference"))]
    args.items += ["-O", str(output)]
    args.many_file("--comparison", comparison_vcfs, "comparison_vcfs")
    args.file("--dbsnp", dbsnp, "dbsnp")
    args.file("--gold-standard", gold_standard, "gold_standard")
    args.many_text("-EV", eval_modules, "eval_modules")
    args.many_text("-ST", stratification_modules, "stratification_modules")
    args.on(
        "--do-not-use-all-standard-modules",
        not use_all_standard_modules,
        "use_all_standard_modules",
    )
    args.on(
        "--do-not-use-all-standard-stratifications",
        not use_all_standard_stratifications,
        "use_all_standard_stratifications",
    )
    args.many_text("--knownNames", known_names, "known_names")
    args.on("--merge-evals", merge_evals, "merge_evals")
    args.on(
        "--require-strict-allele-match", require_strict_allele_match, "require_strict_allele_match"
    )
    args.on("--keep-ac0", keep_ac0, "keep_ac0")
    args.many_text("--sample", samples, "samples")
    args.integer("--sample-ploidy", sample_ploidy, "sample_ploidy", minimum=1, maximum=64)
    args.file("--pedigree", pedigree, "pedigree")
    args.file("--ancestral-alignments", ancestral_alignments, "ancestral_alignments")
    args.file("--known-cnvs", known_cnvs, "known_cnvs")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "VariantEval",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="evaluate_variants_gatk",
    )
    result["output"] = str(output) if output.is_file() else None
    result["tables"] = (
        _parse_gatk_report(output.read_text(encoding="utf-8", errors="replace"))
        if result["status"] == "completed"
        else {}
    )
    return result


# ---------------------------------------------------------------------------
# Recalibration (VQSR) and tranche filtering
# ---------------------------------------------------------------------------


def recalibrate_variants_vqsr(
    vcfs: Sequence[str],
    resources: Sequence[Mapping[str, Any]],
    annotations: Sequence[str],
    recal_path: str,
    tranches_path: str,
    *,
    mode: str = "SNP",
    truth_sensitivity_tranches: Sequence[float] | None = None,
    max_gaussians: int | None = None,
    max_negative_gaussians: int | None = None,
    max_iterations: int | None = None,
    k_means_iterations: int | None = None,
    minimum_bad_variants: int | None = None,
    maximum_training_variants: int | None = None,
    bad_lod_score_cutoff: float | None = None,
    target_titv: float | None = None,
    standard_deviation_threshold: float | None = None,
    output_model: str | None = None,
    use_allele_specific_annotations: bool = False,
    trust_all_polymorphic: bool = False,
    ignore_all_filters: bool = False,
    ignore_filters: Sequence[str] | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a recalibration model that scores variant quality for filtering (gatk
    VariantRecalibrator).

    ``resources`` are mappings ``{"path", "name", "tags"}`` where ``tags`` carry the VQSR roles, for
    example ``{"path": "hapmap.vcf", "name": "hapmap", "tags": {"known": False, "training": True,
    "truth": True, "prior": 15.0}}``. ``annotations`` are the INFO annotations used by the model
    (-an), ``mode`` is SNP, INDEL or BOTH. Writes the recal table to ``recal_path`` and the tranches
    file to ``tranches_path``; both feed ``apply_vqsr``. The R plot script is not requested.
    """
    inputs = _sequence(vcfs, "vcfs")
    names = _sequence(annotations, "annotations")
    if not inputs or not names or not _sequence(resources, "resources"):
        raise ConfigurationError(
            "vcfs, resources and annotations must each list at least one item."
        )
    recal = _cli.new_file(recal_path, "recal_path")
    tranches = _cli.new_file(tranches_path, "tranches_path")
    model = _cli.new_file(output_model, "output_model") if output_model is not None else None
    args = _Args()
    args.many_file("-V", inputs, "vcfs")
    args.features("--resource", resources, "resources")
    args.many_text("-an", names, "annotations")
    args.items += ["-O", str(recal), "--tranches-file", str(tranches)]
    args.choice("--mode", mode, _MODES, "mode")
    for level in _sequence(truth_sensitivity_tranches, "truth_sensitivity_tranches"):
        args.items += [
            "-tranche",
            repr(_cli.number(level, "truth_sensitivity_tranches", minimum=0.0, maximum=100.0)),
        ]
    args.integer("--max-gaussians", max_gaussians, "max_gaussians", minimum=1)
    args.integer(
        "--max-negative-gaussians", max_negative_gaussians, "max_negative_gaussians", minimum=1
    )
    args.integer("--max-iterations", max_iterations, "max_iterations", minimum=1)
    args.integer("--k-means-iterations", k_means_iterations, "k_means_iterations", minimum=1)
    args.integer("--minimum-bad-variants", minimum_bad_variants, "minimum_bad_variants", minimum=0)
    args.integer(
        "--maximum-training-variants",
        maximum_training_variants,
        "maximum_training_variants",
        minimum=1,
    )
    args.number("--bad-lod-score-cutoff", bad_lod_score_cutoff, "bad_lod_score_cutoff")
    args.number("--target-titv", target_titv, "target_titv", minimum=0.0)
    args.number(
        "--standard-deviation-threshold",
        standard_deviation_threshold,
        "standard_deviation_threshold",
        minimum=0.0,
    )
    if model is not None:
        args.items += ["--output-model", str(model)]
    args.on("-AS", use_allele_specific_annotations, "use_allele_specific_annotations")
    args.on("--trust-all-polymorphic", trust_all_polymorphic, "trust_all_polymorphic")
    args.on("--ignore-all-filters", ignore_all_filters, "ignore_all_filters")
    args.many_text("--ignore-filter", ignore_filters, "ignore_filters")
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    expected = [recal, tranches] + ([model] if model is not None else [])
    result = _gatk(
        "VariantRecalibrator",
        args,
        cwd=_common_folder(expected),
        expected=expected,
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="recalibrate_variants_vqsr",
    )
    result["recal_file"] = str(recal) if recal.is_file() else None
    result["tranches_file"] = str(tranches) if tranches.is_file() else None
    result["model_file"] = str(model) if model is not None and model.is_file() else None
    return result


def apply_vqsr(
    vcf: str,
    recal_file: str,
    output_path: str,
    *,
    tranches_file: str | None = None,
    truth_sensitivity_filter_level: float | None = None,
    lod_score_cutoff: float | None = None,
    mode: str = "SNP",
    use_allele_specific_annotations: bool = False,
    exclude_filtered: bool = False,
    ignore_all_filters: bool = False,
    ignore_filters: Sequence[str] | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply a score cutoff to filter variants based on a recalibration table (gatk ApplyVQSR).

    ``recal_file`` and ``tranches_file`` come from ``recalibrate_variants_vqsr``. Filtering starts
    at ``truth_sensitivity_filter_level`` (percent, needs the tranches file) or at
    ``lod_score_cutoff``; ``mode`` is SNP, INDEL or BOTH. Each variant is annotated with its VQSLOD;
    ``exclude_filtered`` drops filtered sites. Result keys: ``output`` and ``summary`` (record
    count, FILTER value counts, samples).
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.items += ["--recal-file", str(_cli.existing_file(recal_file, "recal_file"))]
    args.file("--tranches-file", tranches_file, "tranches_file")
    args.number(
        "--truth-sensitivity-filter-level",
        truth_sensitivity_filter_level,
        "truth_sensitivity_filter_level",
        minimum=0.0,
        maximum=100.0,
    )
    args.number("--lod-score-cutoff", lod_score_cutoff, "lod_score_cutoff")
    args.choice("--mode", mode, _MODES, "mode")
    args.on("-AS", use_allele_specific_annotations, "use_allele_specific_annotations")
    args.on("--exclude-filtered", exclude_filtered, "exclude_filtered")
    args.on("--ignore-all-filters", ignore_all_filters, "ignore_all_filters")
    args.many_text("--ignore-filter", ignore_filters, "ignore_filters")
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "ApplyVQSR",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="apply_vqsr",
    )
    return _with_summary(result, output)


def filter_variant_tranches(
    vcf: str,
    resources: Sequence[str | Mapping[str, Any]],
    output_path: str,
    *,
    snp_tranches: Sequence[float] | None = None,
    indel_tranches: Sequence[float] | None = None,
    info_key: str | None = None,
    invalidate_previous_filters: bool = False,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply tranche filtering from a score in the INFO field and truth resources (gatk
    FilterVariantTranches).

    ``resources`` are validated VCFs of known common variation; ``snp_tranches`` and
    ``indel_tranches`` are sensitivity levels in percent (for example 99.9). ``info_key`` names the
    score annotation (GATK default CNN_2D). ``invalidate_previous_filters`` clears existing FILTER
    values. Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    if not _sequence(resources, "resources"):
        raise ConfigurationError("resources must list at least one VCF.")
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.features("--resource", resources, "resources")
    for flag, values, name in (
        ("--snp-tranche", snp_tranches, "snp_tranches"),
        ("--indel-tranche", indel_tranches, "indel_tranches"),
    ):
        for level in _sequence(values, name):
            args.items += [flag, repr(_cli.number(level, name, minimum=0.0, maximum=100.0))]
    args.text("--info-key", info_key, "info_key")
    args.on(
        "--invalidate-previous-filters", invalidate_previous_filters, "invalidate_previous_filters"
    )
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "FilterVariantTranches",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_variant_tranches",
    )
    return _with_summary(result, output)


# ---------------------------------------------------------------------------
# Variant-annotation models (extract, train, score) and the CNN scorer
# ---------------------------------------------------------------------------


def _prefix(value: str, name: str) -> Path:
    """An output file prefix: the folder must exist and no file of that prefix may exist yet."""
    path = Path(_cli.text(value, name)).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    clashes = sorted(p.name for p in path.parent.glob(f"{path.name}.*"))
    if clashes:
        raise ConfigurationError(f"{name} would overwrite existing files: {clashes[:3]}")
    return path


def _prefixed_outputs(prefix: Path) -> list[str]:
    return sorted(str(p) for p in prefix.parent.glob(f"{prefix.name}.*") if p.is_file())


def extract_variant_annotations(
    vcf: str,
    annotations: Sequence[str],
    output_prefix: str,
    *,
    resources: Sequence[str | Mapping[str, Any]] | None = None,
    modes: Sequence[str] | None = None,
    maximum_number_of_unlabeled_variants: int | None = None,
    reservoir_sampling_random_seed: int | None = None,
    resource_matching_strategy: str | None = None,
    omit_alleles_in_hdf5: bool = False,
    do_not_gzip_vcf_output: bool = False,
    do_not_trust_all_polymorphic: bool = False,
    ignore_all_filters: bool = False,
    ignore_filters: Sequence[str] | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract site-level annotations, labels and metadata from a VCF to HDF5 files (gatk
    ExtractVariantAnnotations).

    ``annotations`` are the INFO annotations to extract. ``resources`` are labelling VCFs given as
    paths or as ``{"path", "name", "tags"}`` mappings (tags such as ``{"training": True}``).
    ``modes`` selects SNP and/or INDEL (GATK extracts both by default). Output files start with
    ``output_prefix`` (``<prefix>.annot.hdf5``, ``<prefix>.vcf.gz``, and the unlabeled / per-mode
    HDF5 files). Result key ``prefix_outputs`` lists every file written with that prefix.
    ExtractVariantAnnotations is a beta GATK tool.
    """
    source = _input_vcf(vcf)
    prefix = _prefix(output_prefix, "output_prefix")
    names = _sequence(annotations, "annotations")
    if not names:
        raise ConfigurationError("annotations must list at least one annotation.")
    args = _Args()
    args.items += ["-V", str(source), "-O", str(prefix)]
    args.many_text("-A", names, "annotations")
    args.features("--resource", resources, "resources")
    args.many_choice("--mode", modes, _VARIANT_TYPES, "modes")
    args.integer(
        "--maximum-number-of-unlabeled-variants",
        maximum_number_of_unlabeled_variants,
        "maximum_number_of_unlabeled_variants",
        minimum=0,
    )
    args.integer(
        "--reservoir-sampling-random-seed",
        reservoir_sampling_random_seed,
        "reservoir_sampling_random_seed",
    )
    args.choice(
        "--resource-matching-strategy",
        resource_matching_strategy,
        _MATCHING_STRATEGIES,
        "resource_matching_strategy",
    )
    args.on("--omit-alleles-in-hdf5", omit_alleles_in_hdf5, "omit_alleles_in_hdf5")
    args.on("--do-not-gzip-vcf-output", do_not_gzip_vcf_output, "do_not_gzip_vcf_output")
    args.on(
        "--do-not-trust-all-polymorphic",
        do_not_trust_all_polymorphic,
        "do_not_trust_all_polymorphic",
    )
    args.on("--ignore-all-filters", ignore_all_filters, "ignore_all_filters")
    args.many_text("--ignore-filter", ignore_filters, "ignore_filters")
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    annotation_hdf5 = Path(f"{prefix}.annot.hdf5")
    result = _gatk(
        "ExtractVariantAnnotations",
        args,
        cwd=prefix.parent,
        expected=[annotation_hdf5],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="extract_variant_annotations",
    )
    result["annotations_hdf5"] = str(annotation_hdf5) if annotation_hdf5.is_file() else None
    result["prefix_outputs"] = _prefixed_outputs(prefix)
    return result


def train_variant_annotations_model(
    annotations_hdf5: str,
    output_prefix: str,
    *,
    unlabeled_annotations_hdf5: str | None = None,
    model_backend: str = "PYTHON_IFOREST",
    hyperparameters_json: str | None = None,
    python_script: str | None = None,
    modes: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train a model that scores variant calls from site-level annotations (gatk
    TrainVariantAnnotationsModel).

    ``annotations_hdf5`` is the labeled file written by ``extract_variant_annotations``; with
    ``unlabeled_annotations_hdf5`` a positive-unlabeled approach is used, otherwise positive-only.
    The backend is ``PYTHON_IFOREST`` (GATK default; scikit-learn IsolationForest, positive-only,
    needs the GATK Python environment), ``PYTHON_SCRIPT`` (with ``python_script``) or ``JAVA_BGMM``
    (needs ``hyperparameters_json``; GATK 4.7.0.0 has not implemented it yet and stops with a
    NotImplementedException). Output files start with ``output_prefix`` (for example
    ``<prefix>.snp.scorer.pkl`` and ``<prefix>.snp.trainingScores.hdf5``); result key
    ``prefix_outputs`` lists them.
    """
    labeled = _cli.existing_file(annotations_hdf5, "annotations_hdf5")
    prefix = _prefix(output_prefix, "output_prefix")
    backend = _cli.choice(model_backend, _MODEL_BACKENDS, "model_backend")
    if backend == "PYTHON_SCRIPT" and python_script is None:
        raise ConfigurationError("model_backend PYTHON_SCRIPT needs python_script.")
    if backend == "JAVA_BGMM" and hyperparameters_json is None:
        raise ConfigurationError("model_backend JAVA_BGMM needs hyperparameters_json.")
    args = _Args()
    args.items += [
        "--annotations-hdf5",
        str(labeled),
        "-O",
        str(prefix),
        "--model-backend",
        backend,
    ]
    args.file(
        "--unlabeled-annotations-hdf5", unlabeled_annotations_hdf5, "unlabeled_annotations_hdf5"
    )
    args.file("--hyperparameters-json", hyperparameters_json, "hyperparameters_json")
    args.file("--python-script", python_script, "python_script")
    args.many_choice("--mode", modes, _VARIANT_TYPES, "modes")
    result = _gatk(
        "TrainVariantAnnotationsModel",
        args,
        cwd=prefix.parent,
        expected=[],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="train_variant_annotations_model",
    )
    result["prefix_outputs"] = _prefixed_outputs(prefix)
    if result["status"] == "completed" and not result["prefix_outputs"]:
        result["status"] = "failed"
    return result


def score_variant_annotations(
    vcf: str,
    annotations: Sequence[str],
    model_prefix: str,
    output_prefix: str,
    *,
    resources: Sequence[str | Mapping[str, Any]] | None = None,
    modes: Sequence[str] | None = None,
    model_backend: str = "PYTHON_IFOREST",
    python_script: str | None = None,
    snp_calibration_sensitivity_threshold: float | None = None,
    indel_calibration_sensitivity_threshold: float | None = None,
    score_key: str | None = None,
    calibration_sensitivity_key: str | None = None,
    low_score_filter_name: str | None = None,
    snp_key: str | None = None,
    double_format: str | None = None,
    resource_matching_strategy: str | None = None,
    omit_alleles_in_hdf5: bool = False,
    do_not_gzip_vcf_output: bool = False,
    do_not_trust_all_polymorphic: bool = False,
    ignore_all_filters: bool = False,
    ignore_filters: Sequence[str] | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score variant calls from site-level annotations with a previously trained model (gatk
    ScoreVariantAnnotations).

    ``model_prefix`` is the ``output_prefix`` given to ``train_variant_annotations_model`` and
    ``model_backend`` must match the backend used there (JAVA_BGMM, PYTHON_IFOREST or
    PYTHON_SCRIPT). ``resources`` carrying a calibration role (tags such as ``{"calibration":
    True}``) allow SCORE to be turned into a calibration sensitivity;
    ``snp_/indel_calibration_sensitivity_threshold`` hard-filter sites at or above that sensitivity.
    Output files start with ``output_prefix`` (``<prefix>.vcf.gz`` with SCORE and, when calibrated,
    CALIBRATION_SENSITIVITY, and the annotation HDF5 files). Result keys: ``scored_vcf`` and
    ``summary``.
    """
    source = _input_vcf(vcf)
    prefix = _prefix(output_prefix, "output_prefix")
    names = _sequence(annotations, "annotations")
    if not names:
        raise ConfigurationError("annotations must list at least one annotation.")
    backend = _cli.choice(model_backend, _MODEL_BACKENDS, "model_backend")
    if backend == "PYTHON_SCRIPT" and python_script is None:
        raise ConfigurationError("model_backend PYTHON_SCRIPT needs python_script.")
    args = _Args()
    args.items += [
        "-V",
        str(source),
        "-O",
        str(prefix),
        "--model-prefix",
        _cli.text(model_prefix, "model_prefix"),
    ]
    args.items += ["--model-backend", backend]
    args.many_text("-A", names, "annotations")
    args.features("--resource", resources, "resources")
    args.many_choice("--mode", modes, _VARIANT_TYPES, "modes")
    args.file("--python-script", python_script, "python_script")
    args.number(
        "--snp-calibration-sensitivity-threshold",
        snp_calibration_sensitivity_threshold,
        "snp_calibration_sensitivity_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    args.number(
        "--indel-calibration-sensitivity-threshold",
        indel_calibration_sensitivity_threshold,
        "indel_calibration_sensitivity_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    args.text("--score-key", score_key, "score_key")
    args.text(
        "--calibration-sensitivity-key", calibration_sensitivity_key, "calibration_sensitivity_key"
    )
    args.text("--low-score-filter-name", low_score_filter_name, "low_score_filter_name")
    args.text("--snp-key", snp_key, "snp_key")
    args.text("--double-format", double_format, "double_format")
    args.choice(
        "--resource-matching-strategy",
        resource_matching_strategy,
        _MATCHING_STRATEGIES,
        "resource_matching_strategy",
    )
    args.on("--omit-alleles-in-hdf5", omit_alleles_in_hdf5, "omit_alleles_in_hdf5")
    args.on("--do-not-gzip-vcf-output", do_not_gzip_vcf_output, "do_not_gzip_vcf_output")
    args.on(
        "--do-not-trust-all-polymorphic",
        do_not_trust_all_polymorphic,
        "do_not_trust_all_polymorphic",
    )
    args.on("--ignore-all-filters", ignore_all_filters, "ignore_all_filters")
    args.many_text("--ignore-filter", ignore_filters, "ignore_filters")
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    scored = Path(f"{prefix}.vcf" if do_not_gzip_vcf_output else f"{prefix}.vcf.gz")
    result = _gatk(
        "ScoreVariantAnnotations",
        args,
        cwd=prefix.parent,
        expected=[scored],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="score_variant_annotations",
    )
    result["scored_vcf"] = str(scored) if scored.is_file() else None
    result["summary"] = _vcf_summary(scored) if result["status"] == "completed" else None
    result["prefix_outputs"] = _prefixed_outputs(prefix)
    return result


def score_variants_with_cnn(
    vcf: str,
    reference: str,
    output_path: str,
    *,
    reads: str | None = None,
    tensor_type: str = "reference",
    batch_size: int | None = None,
    accelerator: str | None = None,
    random_seed: int | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate a VCF with scores of a PyTorch convolutional neural network (gatk NVScoreVariants).

    ``tensor_type`` is ``reference`` (1D reference tensors) or ``read_tensor`` (2D tensors; needs
    ``reads``, a BAM file). ``accelerator`` is auto, cpu, cuda, mps, tpu and so on. GATK launches a
    Python process, so the GATK conda environment (PyTorch, PyTorch Lightning, Biopython and the
    ``scorevariants`` package of ``gatkPythonPackageArchive.zip``; ``gatkcondaenv.yml`` in the GATK
    release) must be active, or the archive folder must be on PYTHONPATH; with PyTorch 2.6 or
    later also set ``TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1`` because the shipped weights are full
    pickles. The 1D model reads the INFO annotations MQ, DP, SOR, FS, QD, MQRankSum and
    ReadPosRankSum, which must be defined in the VCF header. The score is written to INFO as
    CNN_1D (``reference``) or CNN_2D (``read_tensor``). NVScoreVariants is a beta GATK tool and the
    CNN weights ship inside the GATK jar. Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    kind = _cli.choice(tensor_type, _TENSOR_TYPES, "tensor_type")
    if kind == "read_tensor" and reads is None:
        raise ConfigurationError("tensor_type read_tensor needs reads (a BAM file).")
    args = _Args()
    args.items += [
        "-V",
        str(source),
        "-R",
        str(_cli.existing_file(reference, "reference")),
        "-O",
        str(output),
    ]
    args.items += ["--tensor-type", kind]
    args.file("-I", reads, "reads")
    args.integer("--batch-size", batch_size, "batch_size", minimum=1)
    args.text("--accelerator", accelerator, "accelerator")
    args.integer("--random-seed", random_seed, "random_seed")
    result = _gatk(
        "NVScoreVariants",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="score_variants_with_cnn",
    )
    return _with_summary(result, output)


# ---------------------------------------------------------------------------
# Somatic and mitochondrial filtering
# ---------------------------------------------------------------------------


def create_somatic_panel_of_normals(
    genomicsdb_workspace: str,
    output_path: str,
    *,
    reference: str | None = None,
    germline_resource: str | None = None,
    min_sample_count: int | None = None,
    max_germline_probability: float | None = None,
    call_genotypes: bool = False,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Make a panel of normals for Mutect2 from a GenomicsDB workspace of normal-sample calls (gatk
    CreateSomaticPanelOfNormals).

    ``genomicsdb_workspace`` is the folder written by ``gatk GenomicsDBImport`` over the Mutect2
    ``--max-mnp-distance 0`` VCFs of the normals (passed to GATK as ``gendb://<folder>``). A site
    enters the panel when at least ``min_sample_count`` samples carry it (GATK default 2);
    ``germline_resource`` (a gnomAD-like VCF with allele fractions) and ``max_germline_probability``
    drop probable germline variants. CreateSomaticPanelOfNormals is a beta GATK tool. Result keys:
    ``output`` and ``summary``.
    """
    workspace = Path(_cli.text(genomicsdb_workspace, "genomicsdb_workspace")).expanduser().resolve()
    if not workspace.is_dir():
        raise ConfigurationError(f"genomicsdb_workspace must be an existing folder: {workspace}")
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", f"gendb://{workspace}", "-O", str(output)]
    args.file("-R", reference, "reference")
    args.file("--germline-resource", germline_resource, "germline_resource")
    args.integer("--min-sample-count", min_sample_count, "min_sample_count", minimum=1)
    args.number(
        "--max-germline-probability",
        max_germline_probability,
        "max_germline_probability",
        minimum=0.0,
        maximum=1.0,
    )
    args.on("--call-genotypes", call_genotypes, "call_genotypes")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "CreateSomaticPanelOfNormals",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="create_somatic_panel_of_normals",
    )
    return _with_summary(result, output)


def filter_mutect_calls(
    vcf: str,
    reference: str,
    output_path: str,
    *,
    stats: str | None = None,
    filtering_stats_path: str | None = None,
    contamination_tables: Sequence[str] | None = None,
    contamination_estimate: float | None = None,
    tumor_segmentation: Sequence[str] | None = None,
    orientation_bias_priors: Sequence[str] | None = None,
    threshold_strategy: str | None = None,
    false_discovery_rate: float | None = None,
    f_score_beta: float | None = None,
    initial_threshold: float | None = None,
    min_allele_fraction: float | None = None,
    unique_alt_read_count: int | None = None,
    min_median_base_quality: int | None = None,
    min_median_mapping_quality: int | None = None,
    min_median_read_position: int | None = None,
    max_events_in_region: int | None = None,
    max_alt_allele_count: int | None = None,
    normal_p_value_threshold: float | None = None,
    microbial_mode: bool = False,
    mitochondria_mode: bool = False,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter somatic SNVs and indels called by Mutect2 (gatk FilterMutectCalls).

    ``vcf`` is the Mutect2 VCF and ``stats`` its ``.stats`` file (GATK default: ``<vcf>.stats``).
    Contamination comes from ``contamination_tables`` / ``tumor_segmentation``
    (CalculateContamination) or ``contamination_estimate``; ``orientation_bias_priors`` are the
    ``.tar.gz`` priors of LearnReadOrientationModel. ``threshold_strategy`` is OPTIMAL_F_SCORE (GATK
    default), FALSE_DISCOVERY_RATE or CONSTANT. ``microbial_mode`` and ``mitochondria_mode`` switch
    to GATK's microbial and mitochondrial filter defaults. Other filter parameters keep GATK's
    defaults unless set. Result keys: ``output``, ``summary`` (FILTER value counts) and
    ``filtering_stats``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    stats_out = (
        _cli.new_file(filtering_stats_path, "filtering_stats_path")
        if filtering_stats_path is not None
        else None
    )
    args = _Args()
    args.items += [
        "-V",
        str(source),
        "-R",
        str(_cli.existing_file(reference, "reference")),
        "-O",
        str(output),
    ]
    args.file("--stats", stats, "stats")
    if stats_out is not None:
        args.items += ["--filtering-stats", str(stats_out)]
    args.many_file("--contamination-table", contamination_tables, "contamination_tables")
    args.number(
        "--contamination-estimate",
        contamination_estimate,
        "contamination_estimate",
        minimum=0.0,
        maximum=1.0,
    )
    args.many_file("--tumor-segmentation", tumor_segmentation, "tumor_segmentation")
    args.many_file(
        "--orientation-bias-artifact-priors", orientation_bias_priors, "orientation_bias_priors"
    )
    args.choice(
        "--threshold-strategy", threshold_strategy, _THRESHOLD_STRATEGIES, "threshold_strategy"
    )
    args.number(
        "--false-discovery-rate",
        false_discovery_rate,
        "false_discovery_rate",
        minimum=0.0,
        maximum=1.0,
    )
    args.number("--f-score-beta", f_score_beta, "f_score_beta", minimum=0.0)
    args.number(
        "--initial-threshold", initial_threshold, "initial_threshold", minimum=0.0, maximum=1.0
    )
    args.number(
        "--min-allele-fraction",
        min_allele_fraction,
        "min_allele_fraction",
        minimum=0.0,
        maximum=1.0,
    )
    args.integer(
        "--unique-alt-read-count", unique_alt_read_count, "unique_alt_read_count", minimum=0
    )
    args.integer("--min-median-base-quality", min_median_base_quality, "min_median_base_quality")
    args.integer(
        "--min-median-mapping-quality", min_median_mapping_quality, "min_median_mapping_quality"
    )
    args.integer("--min-median-read-position", min_median_read_position, "min_median_read_position")
    args.integer("--max-events-in-region", max_events_in_region, "max_events_in_region", minimum=0)
    args.integer("--max-alt-allele-count", max_alt_allele_count, "max_alt_allele_count", minimum=1)
    args.number(
        "--normal-p-value-threshold",
        normal_p_value_threshold,
        "normal_p_value_threshold",
        minimum=0.0,
        maximum=1.0,
    )
    args.on("--microbial-mode", microbial_mode, "microbial_mode")
    args.on("--mitochondria-mode", mitochondria_mode, "mitochondria_mode")
    args.intervals(intervals, exclude_intervals)
    expected = [output] + ([stats_out] if stats_out is not None else [])
    result = _gatk(
        "FilterMutectCalls",
        args,
        cwd=_common_folder(expected),
        expected=expected,
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_mutect_calls",
    )
    result["filtering_stats"] = (
        str(stats_out) if stats_out is not None and stats_out.is_file() else None
    )
    return _with_summary(result, output)


def filter_alignment_artifacts(
    vcf: str,
    reads: str,
    reference: str,
    bwa_mem_index_image: str,
    output_path: str,
    *,
    bam_output_path: str | None = None,
    dont_skip_filtered_variants: bool = False,
    kmer_size: int | None = None,
    drop_ratio: float | None = None,
    min_seed_length: int | None = None,
    seed_split_factor: float | None = None,
    indel_start_tolerance: int | None = None,
    max_reasonable_fragment_length: int | None = None,
    min_aligner_score_difference_per_base: float | None = None,
    min_mismatch_difference_per_base: float | None = None,
    num_regular_contigs: int | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter alignment artifacts from a VCF callset by realigning the reads around each variant
    (gatk FilterAlignmentArtifacts).

    ``vcf`` is a (Mutect2) callset, ``reads`` the BAM it was called from, ``bwa_mem_index_image``
    the image written by ``gatk BwaMemIndexImageCreator`` for ``reference``. Variants whose reads
    realign better elsewhere get the ``alignment-artifact`` FILTER. ``dont_skip_filtered_variants``
    also realigns variants that are already filtered; ``bam_output_path`` writes the assembled
    haplotypes. FilterAlignmentArtifacts is an experimental GATK tool. Result keys: ``output`` and
    ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    bam_out = (
        _cli.new_file(bam_output_path, "bam_output_path") if bam_output_path is not None else None
    )
    args = _Args()
    args.items += ["-V", str(source), "-I", str(_cli.existing_file(reads, "reads"))]
    args.items += ["-R", str(_cli.existing_file(reference, "reference")), "-O", str(output)]
    args.items += [
        "--bwa-mem-index-image",
        str(_cli.existing_file(bwa_mem_index_image, "bwa_mem_index_image")),
    ]
    if bam_out is not None:
        args.items += ["--bam-output", str(bam_out)]
    args.on(
        "--dont-skip-filtered-variants", dont_skip_filtered_variants, "dont_skip_filtered_variants"
    )
    args.integer("--kmer-size", kmer_size, "kmer_size", minimum=1)
    args.number("--drop-ratio", drop_ratio, "drop_ratio", minimum=0.0, maximum=1.0)
    args.integer("--minimum-seed-length", min_seed_length, "min_seed_length", minimum=1)
    args.number("--seed-split-factor", seed_split_factor, "seed_split_factor", minimum=0.0)
    args.integer(
        "--indel-start-tolerance", indel_start_tolerance, "indel_start_tolerance", minimum=0
    )
    args.integer(
        "--max-reasonable-fragment-length",
        max_reasonable_fragment_length,
        "max_reasonable_fragment_length",
        minimum=1,
    )
    args.number(
        "--min-aligner-score-difference-per-base",
        min_aligner_score_difference_per_base,
        "min_aligner_score_difference_per_base",
        minimum=0.0,
    )
    args.number(
        "--min-mismatch-difference-per-base",
        min_mismatch_difference_per_base,
        "min_mismatch_difference_per_base",
        minimum=0.0,
    )
    args.integer("--num-regular-contigs", num_regular_contigs, "num_regular_contigs", minimum=1)
    args.intervals(intervals, exclude_intervals)
    expected = [output] + ([bam_out] if bam_out is not None else [])
    result = _gatk(
        "FilterAlignmentArtifacts",
        args,
        cwd=_common_folder(expected),
        expected=expected,
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_alignment_artifacts",
    )
    result["bam_output"] = str(bam_out) if bam_out is not None and bam_out.is_file() else None
    return _with_summary(result, output)


def filter_mt_low_heteroplasmy(
    vcf: str,
    output_path: str,
    *,
    low_het_threshold: float | None = None,
    max_allowed_low_hets: int | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter all low-heteroplasmy sites when too many of them pass the other filters (gatk
    MTLowHeteroplasmyFilterTool).

    A site is low heteroplasmy when its allele fraction is below ``low_het_threshold`` (GATK default
    0.1); when more than ``max_allowed_low_hets`` (default 3) such sites pass every other filter,
    all of them are filtered. Meant for mitochondrial Mutect2 output. Result keys: ``output`` and
    ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.number(
        "--low-het-threshold", low_het_threshold, "low_het_threshold", minimum=0.0, maximum=1.0
    )
    args.integer("--max-allowed-low-hets", max_allowed_low_hets, "max_allowed_low_hets", minimum=0)
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "MTLowHeteroplasmyFilterTool",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_mt_low_heteroplasmy",
    )
    return _with_summary(result, output)


def filter_numt_variants(
    vcf: str,
    output_path: str,
    *,
    autosomal_coverage: float | None = None,
    max_numt_autosomal_copies: float | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag alleles that may be nuclear mitochondrial segments (NuMTs) from coverage and allele
    depth (gatk NuMTFilterTool).

    Uses the median autosomal coverage (``autosomal_coverage``, GATK default 0.0) and the allele
    depth: alleles whose depth is below ``max_numt_autosomal_copies`` (default 4.0) times half the
    autosomal coverage get the ``possible_numt`` FILTER. Meant for mitochondrial Mutect2 output.
    Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.number("--autosomal-coverage", autosomal_coverage, "autosomal_coverage", minimum=0.0)
    args.number(
        "--max-numt-autosomal-copies",
        max_numt_autosomal_copies,
        "max_numt_autosomal_copies",
        minimum=0.0,
    )
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "NuMTFilterTool",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_numt_variants",
    )
    return _with_summary(result, output)


# ---------------------------------------------------------------------------
# Hard filtering and annotation
# ---------------------------------------------------------------------------


def filter_vcf_variants_gatk(
    vcf: str,
    output_path: str,
    *,
    filter_expressions: Sequence[str] | None = None,
    filter_names: Sequence[str] | None = None,
    genotype_filter_expressions: Sequence[str] | None = None,
    genotype_filter_names: Sequence[str] | None = None,
    invert_filter_expression: bool = False,
    invert_genotype_filter_expression: bool = False,
    missing_values_evaluate_as_failing: bool = False,
    set_filtered_genotype_to_no_call: bool = False,
    invalidate_previous_filters: bool = False,
    mask: str | None = None,
    mask_name: str | None = None,
    mask_extension: int | None = None,
    mask_description: str | None = None,
    filter_not_in_mask: bool = False,
    cluster_size: int | None = None,
    cluster_window_size: int | None = None,
    apply_allele_specific_filters: bool = False,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter variant calls on INFO and/or FORMAT annotations by JEXL expressions (gatk
    VariantFiltration).

    ``filter_expressions`` (for example ``"QD < 2.0"``) test INFO fields and
    ``genotype_filter_expressions`` test FORMAT fields; ``filter_names`` / ``genotype_filter_names``
    give the label written to FILTER for each expression (the lists must have the same length). A
    ``mask`` file filters records overlapping it (label ``mask_name``, default Mask);
    ``cluster_size`` SNPs within ``cluster_window_size`` bases are filtered as a cluster. Failing
    records are marked, not removed. Result keys: ``output`` and ``summary`` (FILTER value counts).
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    expressions = _sequence(filter_expressions, "filter_expressions")
    names = _sequence(filter_names, "filter_names")
    genotype_expressions = _sequence(genotype_filter_expressions, "genotype_filter_expressions")
    genotype_names = _sequence(genotype_filter_names, "genotype_filter_names")
    if len(expressions) != len(names) or len(genotype_expressions) != len(genotype_names):
        raise ConfigurationError(
            "each filter expression needs exactly one filter name (likewise for genotype filters)."
        )
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    for expression, name in zip(expressions, names, strict=True):
        args.items += [
            "--filter-expression",
            _cli.text(expression, "filter_expressions"),
            "--filter-name",
            _cli.text(name, "filter_names"),
        ]
    for expression, name in zip(genotype_expressions, genotype_names, strict=True):
        args.items += [
            "--genotype-filter-expression",
            _cli.text(expression, "genotype_filter_expressions"),
            "--genotype-filter-name",
            _cli.text(name, "genotype_filter_names"),
        ]
    args.on("--invert-filter-expression", invert_filter_expression, "invert_filter_expression")
    args.on(
        "--invert-genotype-filter-expression",
        invert_genotype_filter_expression,
        "invert_genotype_filter_expression",
    )
    args.on(
        "--missing-values-evaluate-as-failing",
        missing_values_evaluate_as_failing,
        "missing_values_evaluate_as_failing",
    )
    args.on(
        "--set-filtered-genotype-to-no-call",
        set_filtered_genotype_to_no_call,
        "set_filtered_genotype_to_no_call",
    )
    args.on(
        "--invalidate-previous-filters", invalidate_previous_filters, "invalidate_previous_filters"
    )
    args.file("--mask", mask, "mask")
    args.text("--mask-name", mask_name, "mask_name")
    args.integer("--mask-extension", mask_extension, "mask_extension", minimum=0)
    args.text("--mask-description", mask_description, "mask_description")
    args.on("--filter-not-in-mask", filter_not_in_mask, "filter_not_in_mask")
    args.integer("--cluster-size", cluster_size, "cluster_size", minimum=2)
    args.integer("--cluster-window-size", cluster_window_size, "cluster_window_size", minimum=0)
    args.on(
        "--apply-allele-specific-filters",
        apply_allele_specific_filters,
        "apply_allele_specific_filters",
    )
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "VariantFiltration",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="filter_vcf_variants_gatk",
    )
    return _with_summary(result, output)


def annotate_vcf_gatk(
    vcf: str,
    output_path: str,
    *,
    annotations: Sequence[str] | None = None,
    annotation_groups: Sequence[str] | None = None,
    annotations_to_exclude: Sequence[str] | None = None,
    expressions: Sequence[str] | None = None,
    resources: Sequence[str | Mapping[str, Any]] | None = None,
    resource_allele_concordance: bool = False,
    dbsnp: str | None = None,
    reads: Sequence[str] | None = None,
    pedigree: str | None = None,
    founder_ids: Sequence[str] | None = None,
    enable_all_annotations: bool = False,
    disable_tool_default_annotations: bool = False,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    exclude_intervals: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add annotations to the records of a VCF (gatk VariantAnnotator).

    ``annotations`` (-A, for example ``Coverage``, ``QualByDepth``) and ``annotation_groups`` (-G,
    for example ``StandardAnnotation``) choose what is added, ``annotations_to_exclude`` (-AX)
    removes some. Annotations that need reads (depth, mapping quality) require ``reads`` (BAM/CRAM
    files) and usually ``reference``. ``resources`` add INFO fields from other VCFs (``{"path",
    "name"}``), ``dbsnp`` fills the ID column, and ``expressions`` (-E, ``<resource>.<INFO field>``)
    choose which resource fields are copied. Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.many_text("-A", annotations, "annotations")
    args.many_text("-G", annotation_groups, "annotation_groups")
    args.many_text("-AX", annotations_to_exclude, "annotations_to_exclude")
    args.many_text("-E", expressions, "expressions")
    args.features("--resource", resources, "resources")
    args.on(
        "--resource-allele-concordance", resource_allele_concordance, "resource_allele_concordance"
    )
    args.file("--dbsnp", dbsnp, "dbsnp")
    args.many_file("-I", reads, "reads")
    args.file("--pedigree", pedigree, "pedigree")
    args.many_text("--founder-id", founder_ids, "founder_ids")
    args.on("--enable-all-annotations", enable_all_annotations, "enable_all_annotations")
    args.on(
        "--disable-tool-default-annotations",
        disable_tool_default_annotations,
        "disable_tool_default_annotations",
    )
    args.file("-R", reference, "reference")
    args.intervals(intervals, exclude_intervals)
    result = _gatk(
        "VariantAnnotator",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="annotate_vcf_gatk",
    )
    return _with_summary(result, output)


def print_variants_spark(
    vcf: str,
    output_path: str,
    *,
    reference: str | None = None,
    variant_shard_size: int | None = None,
    use_nio: bool = False,
    spark_master: str | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print (rewrite) the variants of a VCF as one file through the GATK Spark engine (gatk
    PrintVariantsSpark).

    Runs in local Spark mode unless ``spark_master`` names a cluster (for example ``local[4]``); the
    file is read in shards of ``variant_shard_size`` bases (GATK default 10000) and written as a
    single ``output_path``. GATK 4.7.0.0 cannot serialise genotype columns between Spark tasks (Kryo
    error), so the tool works on sites-only VCFs; use ``make_sites_only_vcf_gatk`` first. Result
    keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.file("-R", reference, "reference")
    args.integer("--variantShardSize", variant_shard_size, "variant_shard_size", minimum=1)
    args.on("--use-nio", use_nio, "use_nio")
    if spark_master is not None:
        args.items += [
            "--",
            "--spark-runner",
            "LOCAL" if spark_master.startswith("local") else "SPARK",
        ]
        args.items += ["--spark-master", _cli.text(spark_master, "spark_master")]
    result = _gatk(
        "PrintVariantsSpark",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="print_variants_spark",
    )
    return _with_summary(result, output)


# ---------------------------------------------------------------------------
# VCF header, dictionary, conversion and site-level utilities
# ---------------------------------------------------------------------------


def fix_vcf_header_gatk(
    vcf: str,
    output_path: str,
    *,
    header: str | None = None,
    enforce_same_samples: bool = True,
    check_first_n_records: int | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace or fix a VCF header (gatk FixVcfHeader).

    With ``header`` (a VCF whose header is used) the input header is replaced;
    ``enforce_same_samples`` (GATK default true) requires the same samples in the same order.
    Without ``header`` the tool reads the records and adds every FILTER, INFO and FORMAT definition
    that the records use but the header lacks, with a dummy description; ``check_first_n_records``
    limits that scan (GATK default: all). Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["--INPUT", str(source), "--OUTPUT", str(output)]
    args.file("--HEADER", header, "header")
    args.off("--ENFORCE_SAME_SAMPLES", enforce_same_samples, "enforce_same_samples")
    args.integer(
        "--CHECK_FIRST_N_RECORDS", check_first_n_records, "check_first_n_records", minimum=-1
    )
    result = _gatk(
        "FixVcfHeader",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="fix_vcf_header_gatk",
    )
    return _with_summary(result, output)


def update_vcf_sequence_dictionary_gatk(
    vcf: str,
    output_path: str,
    *,
    source_dictionary: str | None = None,
    reference: str | None = None,
    replace: bool = False,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Update the sequence dictionary (contig lines) of a VCF from a dictionary, variant, alignment
    or reference file (gatk UpdateVCFSequenceDictionary).

    ``source_dictionary`` is a ``.dict``, VCF, BAM or FASTA file; alternatively ``reference``
    supplies the dictionary. ``replace`` forces replacement when the VCF already has contig lines.
    The dictionary must be valid for every variant of the VCF. Result keys: ``output`` and
    ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    if source_dictionary is None and reference is None:
        raise ConfigurationError("Give source_dictionary or reference as the dictionary source.")
    args = _Args()
    args.items += ["-V", str(source), "-O", str(output)]
    args.file("--source-dictionary", source_dictionary, "source_dictionary")
    args.file("-R", reference, "reference")
    args.on("--replace", replace, "replace")
    result = _gatk(
        "UpdateVCFSequenceDictionary",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="update_vcf_sequence_dictionary_gatk",
    )
    return _with_summary(result, output)


def liftover_vcf_gatk(
    vcf: str,
    chain_file: str,
    reference: str,
    output_path: str,
    reject_path: str,
    *,
    min_match: float | None = None,
    recover_swapped_ref_alt: bool = False,
    write_original_position: bool = False,
    write_original_alleles: bool = False,
    tags_to_drop: Sequence[str] | None = None,
    tags_to_reverse: Sequence[str] | None = None,
    allow_missing_fields_in_header: bool = False,
    disable_sort: bool = False,
    warn_on_missing_contig: bool = False,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift a VCF over from one reference build to another with a UCSC chain file (gatk
    LiftoverVcf).

    ``chain_file`` is the chain of the old-to-new build, ``reference`` the FASTA (with ``.dict``) of
    the TARGET build. Lifted records go to ``output_path`` (sorted and indexed unless
    ``disable_sort``) and records that cannot be lifted to ``reject_path``. ``min_match`` is the
    minimum fraction of the variant that must map (GATK default 1.0); ``recover_swapped_ref_alt``
    rescues bi-allelic SNPs whose alleles swapped between builds (``tags_to_reverse`` such as AF are
    transformed x -> 1-x, ``tags_to_drop`` are removed). Result keys: ``output``, ``reject``,
    ``summary`` (lifted) and ``reject_summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    reject = _output_vcf(reject_path, "reject_path")
    if output == reject:
        raise ConfigurationError("output_path and reject_path must differ.")
    args = _Args()
    args.items += ["--INPUT", str(source), "--OUTPUT", str(output), "--REJECT", str(reject)]
    args.items += ["--CHAIN", str(_cli.existing_file(chain_file, "chain_file"))]
    args.items += ["--REFERENCE_SEQUENCE", str(_cli.existing_file(reference, "reference"))]
    args.number("--LIFTOVER_MIN_MATCH", min_match, "min_match", minimum=0.0, maximum=1.0)
    args.on("--RECOVER_SWAPPED_REF_ALT", recover_swapped_ref_alt, "recover_swapped_ref_alt")
    args.on("--WRITE_ORIGINAL_POSITION", write_original_position, "write_original_position")
    args.on("--WRITE_ORIGINAL_ALLELES", write_original_alleles, "write_original_alleles")
    args.many_text("--TAGS_TO_DROP", tags_to_drop, "tags_to_drop")
    args.many_text("--TAGS_TO_REVERSE", tags_to_reverse, "tags_to_reverse")
    args.on(
        "--ALLOW_MISSING_FIELDS_IN_HEADER",
        allow_missing_fields_in_header,
        "allow_missing_fields_in_header",
    )
    args.on("--DISABLE_SORT", disable_sort, "disable_sort")
    args.on("--WARN_ON_MISSING_CONTIG", warn_on_missing_contig, "warn_on_missing_contig")
    result = _gatk(
        "LiftoverVcf",
        args,
        cwd=_common_folder([output, reject]),
        expected=[output, reject],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="liftover_vcf_gatk",
    )
    result["reject"] = str(reject) if reject.is_file() else None
    result["reject_summary"] = _vcf_summary(reject) if result["status"] == "completed" else None
    return _with_summary(result, output)


def make_sites_only_vcf_gatk(
    vcf: str,
    output_path: str,
    *,
    samples: Sequence[str] | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Drop all genotype columns of a VCF/BCF but keep the site-level information, including AN and
    AF (gatk MakeSitesOnlyVcf).

    ``samples`` restricts the site-level information to the listed samples (-S) before genotypes are
    removed. The output format follows the extension of ``output_path`` (.vcf, .vcf.gz, .bcf).
    Result keys: ``output`` and ``summary``.
    """
    source = _input_vcf(vcf)
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["--INPUT", str(source), "--OUTPUT", str(output)]
    args.many_text("--SAMPLE", samples, "samples")
    result = _gatk(
        "MakeSitesOnlyVcf",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="make_sites_only_vcf_gatk",
    )
    return _with_summary(result, output)


def split_vcf_snps_indels_gatk(
    vcf: str,
    snp_output_path: str,
    indel_output_path: str,
    *,
    strict: bool = True,
    sequence_dictionary: str | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a VCF/BCF into a SNP file and an INDEL file with identical headers (gatk SplitVcfs).

    With ``strict`` (GATK default true) a record that is neither SNP nor indel (for example a
    structural variant or mixed type) stops the run; switch it off to drop such records.
    ``sequence_dictionary`` replaces the dictionary of the input. Result keys: ``snp_output``,
    ``indel_output``, ``snp_summary`` and ``indel_summary``.
    """
    source = _input_vcf(vcf)
    snps = _output_vcf(snp_output_path, "snp_output_path")
    indels = _output_vcf(indel_output_path, "indel_output_path")
    if snps == indels:
        raise ConfigurationError("snp_output_path and indel_output_path must differ.")
    args = _Args()
    args.items += ["--INPUT", str(source), "--SNP_OUTPUT", str(snps), "--INDEL_OUTPUT", str(indels)]
    args.off("--STRICT", strict, "strict")
    args.file("--SEQUENCE_DICTIONARY", sequence_dictionary, "sequence_dictionary")
    result = _gatk(
        "SplitVcfs",
        args,
        cwd=_common_folder([snps, indels]),
        expected=[snps, indels],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="split_vcf_snps_indels_gatk",
    )
    done = result["status"] == "completed"
    result["snp_output"] = str(snps) if snps.is_file() else None
    result["indel_output"] = str(indels) if indels.is_file() else None
    result["snp_summary"] = _vcf_summary(snps) if done else None
    result["indel_summary"] = _vcf_summary(indels) if done else None
    return result


def convert_vcf_format_gatk(
    input_path: str,
    output_path: str,
    *,
    require_index: bool = True,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a VCF to BCF or a BCF to VCF (gatk VcfFormatConverter).

    The input and output formats are taken from the file extensions (.vcf, .vcf.gz, .bcf).
    ``require_index`` (GATK default true) fails when the input has no index. Result keys: ``output``
    and ``summary`` (record counts are only available for VCF outputs).
    """
    source = _input_vcf(input_path, "input_path")
    output = _output_vcf(output_path)
    args = _Args()
    args.items += ["--INPUT", str(source), "--OUTPUT", str(output)]
    args.off("--REQUIRE_INDEX", require_index, "require_index")
    result = _gatk(
        "VcfFormatConverter",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="convert_vcf_format_gatk",
    )
    return _with_summary(result, output)


def vcf_to_interval_list_gatk(
    vcf: str,
    output_path: str,
    *,
    include_filtered: bool = False,
    variant_id_method: str | None = None,
    gatk_path: str | None = None,
    java_max_heap_gb: int | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a Picard interval list from the records of a VCF/BCF (gatk VcfToIntervalList).

    Filtered records are skipped unless ``include_filtered``. ``variant_id_method`` controls the
    interval names: CONCAT_ALL (GATK default: all variant IDs joined, or ``interval-<n>``), FIRST,
    LAST or NONE. The interval list file is read back into ``intervals`` as ``[contig, start, end,
    strand, name]`` rows.
    """
    source = _input_vcf(vcf)
    output = _cli.new_file(output_path, "output_path")
    args = _Args()
    args.items += ["--INPUT", str(source), "--OUTPUT", str(output)]
    args.on("--INCLUDE_FILTERED", include_filtered, "include_filtered")
    args.choice("--VARIANT_ID_METHOD", variant_id_method, _ID_METHODS, "variant_id_method")
    result = _gatk(
        "VcfToIntervalList",
        args,
        cwd=output.parent,
        expected=[output],
        gatk_path=gatk_path,
        java_max_heap_gb=java_max_heap_gb,
        timeout_seconds=timeout_seconds,
        operation="vcf_to_interval_list_gatk",
    )
    result["output"] = str(output) if output.is_file() else None
    rows: list[list[object]] = []
    if result["status"] == "completed":
        for line in output.read_text(encoding="utf-8", errors="replace").splitlines():
            if line and not line.startswith("@"):
                fields = line.split("\t")
                rows.append([fields[0], int(fields[1]), int(fields[2]), *fields[3:]])
    result["intervals"] = rows
    return result


__all__ = [
    "annotate_vcf_gatk",
    "apply_vqsr",
    "convert_vcf_format_gatk",
    "create_somatic_panel_of_normals",
    "evaluate_variants_gatk",
    "extract_variant_annotations",
    "filter_alignment_artifacts",
    "filter_mt_low_heteroplasmy",
    "filter_mutect_calls",
    "filter_numt_variants",
    "filter_variant_tranches",
    "filter_vcf_variants_gatk",
    "fix_vcf_header_gatk",
    "liftover_vcf_gatk",
    "make_sites_only_vcf_gatk",
    "print_variants_spark",
    "recalibrate_variants_vqsr",
    "score_variant_annotations",
    "score_variants_with_cnn",
    "split_vcf_snps_indels_gatk",
    "train_variant_annotations_model",
    "update_vcf_sequence_dictionary_gatk",
    "validate_vcf_gatk",
    "vcf_to_interval_list_gatk",
]
_register_public_verbose(globals())
