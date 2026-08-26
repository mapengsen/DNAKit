"""Bounded local execution of strictly whitelisted DNAKit workflows."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Any, Literal, cast

from dnakit.core import ArtifactRef, DNAAlphabet, DNARecord, DNASet, Issue, IssueSeverity
from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.provenance import Provenance
from dnakit.datasets import (
    DeduplicationConfig,
    SplitConfig,
    SplitResult,
    deduplicate,
    split,
)
from dnakit.descriptors import (
    base_composition,
    gc_at_content,
    length_features,
    shannon_entropy,
)
from dnakit.exceptions import ConfigurationError, InputFormatError, SequenceError
from dnakit.fingerprints import kmer_fingerprint
from dnakit.io import ReadConfig, WriteConfig, read, write
from dnakit.standardize import (
    DatasetValidationConfig,
    NormalizationConfig,
    ValidationConfig,
    normalize,
    validate,
)
from dnakit.visualization import SaveConfig, build_html_report, save_html_report
from dnakit.workflows.manifest import (
    RunManifestBuilder,
    artifact_from_path,
    load_manifest,
    save_manifest,
)
from dnakit.workflows.schema import LoadedWorkflow, WorkflowSpec, WorkflowStep, load_workflow

WorkflowProgressStatus = Literal["planned", "started", "succeeded", "failed", "skipped", "blocked"]
WorkflowStepStatus = Literal["planned", "succeeded", "failed", "skipped", "blocked"]
WorkflowRunStatus = Literal["dry-run", "succeeded", "failed"]
ProgressCallback = Callable[["WorkflowProgress"], None]

_OUTPUT_MARKER = ".dnakit-workflow-output-v1"
_MANIFEST_NAME = "run-manifest.json"
_EXECUTION_SCHEMA = "dnakit-workflow-execution-v1"
_MEDIA_TYPES = {
    "fasta": "text/x-fasta",
    "fastq": "text/x-fastq",
    "csv": "text/csv",
    "tsv": "text/tab-separated-values",
    "json": "application/json",
    "jsonl": "application/x-ndjson",
    "genbank": "text/x-genbank",
    "html": "text/html",
}
_FORMAT_ALIASES = {
    "fa": "fasta",
    "fna": "fasta",
    "fq": "fastq",
    "gb": "genbank",
    "gbk": "genbank",
    "ndjson": "jsonl",
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _plain(value: object) -> object:
    return to_json_compatible(value)


def _plain_dict(value: Mapping[str, object]) -> dict[str, object]:
    return cast(dict[str, object], _plain(value))


def _canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            _plain(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Workflow state must remain finite and JSON-compatible.",
            code="INVALID_WORKFLOW_STATE",
        ) from exc


def _sha256(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _error_payload(error: Exception) -> dict[str, object]:
    raw_code = getattr(error, "code", "WORKFLOW_STEP_ERROR")
    code = raw_code if isinstance(raw_code, str) and raw_code else "WORKFLOW_STEP_ERROR"
    payload: dict[str, object] = {
        "type": type(error).__name__,
        "code": code,
        "message": str(error),
    }
    context = getattr(error, "context", None)
    if isinstance(context, Mapping) and context:
        payload["context"] = _plain(context)
    hint = getattr(error, "hint", None)
    if isinstance(hint, str) and hint:
        payload["hint"] = hint
    return payload


@dataclass(frozen=True, slots=True)
class WorkflowProgress:
    """Library-level progress event; DNAKit never prints it itself."""

    run_id: str
    step_id: str
    operation: str
    index: int
    total: int
    status: WorkflowProgressStatus
    message: str


@dataclass(frozen=True, slots=True, init=False)
class WorkflowStepResult:
    """Auditable outcome for one planned workflow step."""

    id: str
    operation: str
    input_ref: str
    status: WorkflowStepStatus
    params: FrozenDict
    started_at: str | None
    finished_at: str | None
    duration_seconds: float
    input_sha256: str | None
    record_count: int | None
    summary: FrozenDict
    artifacts: tuple[FrozenDict, ...]
    error: FrozenDict | None

    def __init__(
        self,
        id: str,
        operation: str,
        input_ref: str,
        status: WorkflowStepStatus,
        params: Mapping[str, object],
        *,
        started_at: str | None = None,
        finished_at: str | None = None,
        duration_seconds: float = 0.0,
        input_sha256: str | None = None,
        record_count: int | None = None,
        summary: Mapping[str, object] | None = None,
        artifacts: tuple[Mapping[str, object], ...] = (),
        error: Mapping[str, object] | None = None,
    ) -> None:
        if not id or not operation or not input_ref:
            raise ConfigurationError("Workflow step result identifiers must be non-empty.")
        if status not in {"planned", "succeeded", "failed", "skipped", "blocked"}:
            raise ConfigurationError("Workflow step result status is invalid.")
        if duration_seconds < 0:
            raise ConfigurationError("Workflow step duration must be non-negative.")
        if record_count is not None and (
            isinstance(record_count, bool) or not isinstance(record_count, int) or record_count < 0
        ):
            raise ConfigurationError("Workflow step record_count must be non-negative or None.")
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "operation", operation)
        object.__setattr__(self, "input_ref", input_ref)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "params", freeze_mapping(params))
        object.__setattr__(self, "started_at", started_at)
        object.__setattr__(self, "finished_at", finished_at)
        object.__setattr__(self, "duration_seconds", float(duration_seconds))
        object.__setattr__(self, "input_sha256", input_sha256)
        object.__setattr__(self, "record_count", record_count)
        object.__setattr__(self, "summary", freeze_mapping(summary))
        object.__setattr__(
            self,
            "artifacts",
            tuple(freeze_mapping(artifact) for artifact in artifacts),
        )
        object.__setattr__(self, "error", None if error is None else freeze_mapping(error))

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], _plain(self))


@dataclass(frozen=True, slots=True)
class WorkflowRunResult:
    """Final local workflow outcome and persisted artifact references."""

    run_id: str
    status: WorkflowRunStatus
    dry_run: bool
    resumed: bool
    spec_sha256: str
    input_sha256: str
    output_dir: str
    manifest_path: str | None
    manifest_artifact: ArtifactRef | None
    steps: tuple[WorkflowStepResult, ...]
    artifacts: tuple[ArtifactRef, ...]

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], _plain(self))


@dataclass(frozen=True, slots=True)
class _State:
    records: DNASet
    results: Mapping[str, object]
    split_result: SplitResult | None = None


@dataclass(frozen=True, slots=True)
class _OperationOutcome:
    state: _State
    summary: Mapping[str, object]
    artifacts: tuple[ArtifactRef, ...] = ()


@dataclass(frozen=True, slots=True)
class _ResumeContext:
    entries: Mapping[str, Mapping[str, object]]
    prior_steps: tuple[Mapping[str, object], ...]


def _emit(
    callback: ProgressCallback | None,
    *,
    run_id: str,
    step_id: str,
    operation: str,
    index: int,
    total: int,
    status: WorkflowProgressStatus,
    message: str,
) -> None:
    if callback is None:
        return
    try:
        callback(WorkflowProgress(run_id, step_id, operation, index, total, status, message))
    except Exception:
        # Progress is advisory and must never change scientific output or leave
        # an already-persisted step in an ambiguous failed state.
        return


def _best_effort_emit(
    callback: ProgressCallback | None,
    *,
    run_id: str,
    step_id: str,
    operation: str,
    index: int,
    total: int,
    status: WorkflowProgressStatus,
    message: str,
) -> None:
    _emit(
        callback,
        run_id=run_id,
        step_id=step_id,
        operation=operation,
        index=index,
        total=total,
        status=status,
        message=message,
    )


def _prepare_output_dir(path: Path, *, create: bool) -> None:
    marker = path / _OUTPUT_MARKER
    if path.exists():
        if path.is_symlink() or not path.is_dir():
            raise ConfigurationError(
                "Workflow output_dir must be a non-symlink directory.",
                code="UNSAFE_WORKFLOW_OUTPUT_DIR",
                context={"output_dir": str(path)},
            )
        nonempty = next(path.iterdir(), None) is not None
        if nonempty and (not marker.is_file() or marker.is_symlink()):
            raise ConfigurationError(
                "Existing non-empty output_dir lacks the DNAKit workflow marker.",
                code="UNSAFE_WORKFLOW_OUTPUT_DIR",
                context={"output_dir": str(path)},
                hint="Choose an empty dedicated directory.",
            )
    elif create:
        path.mkdir(parents=True, exist_ok=False)
    if create and not marker.exists():
        try:
            with marker.open("x", encoding="utf-8") as handle:
                handle.write(_EXECUTION_SCHEMA + "\n")
                handle.flush()
                os.fsync(handle.fileno())
        except FileExistsError as exc:
            if not marker.is_file() or marker.is_symlink():
                raise ConfigurationError(
                    "Workflow output marker is unsafe.",
                    code="UNSAFE_WORKFLOW_OUTPUT_DIR",
                ) from exc


def _ensure_target(root: Path, relative: str, *, create_parent: bool) -> Path:
    raw = Path(relative)
    candidate = root.joinpath(*raw.parts)
    resolved = candidate.resolve(strict=False)
    if not resolved.is_relative_to(root) or resolved == root:
        raise ConfigurationError(
            "Workflow output target escapes output_dir.",
            code="UNSAFE_WORKFLOW_OUTPUT_PATH",
            context={"target": relative},
        )
    current = root
    for part in raw.parts[:-1]:
        current /= part
        if current.exists() and current.is_symlink():
            raise ConfigurationError(
                "Workflow output target traverses a symbolic link.",
                code="UNSAFE_WORKFLOW_OUTPUT_PATH",
                context={"path": str(current)},
            )
    if create_parent:
        resolved.parent.mkdir(parents=True, exist_ok=True)
    if resolved.parent.resolve() != resolved.parent or not resolved.parent.is_relative_to(root):
        raise ConfigurationError(
            "Workflow output parent is unsafe.", code="UNSAFE_WORKFLOW_OUTPUT_PATH"
        )
    if resolved.exists() and resolved.is_symlink():
        raise ConfigurationError(
            "Workflow output target must not be a symbolic link.",
            code="UNSAFE_WORKFLOW_OUTPUT_PATH",
            context={"target": relative},
        )
    return resolved


def _input_media_type(plan: LoadedWorkflow) -> str:
    format_name = plan.spec.input.format
    if format_name is None:
        name = plan.input_path.name.lower()
        if name.endswith(".gz"):
            name = name[:-3]
        suffix = Path(name).suffix.lower().lstrip(".")
        format_name = _FORMAT_ALIASES.get(suffix, suffix)
    return _MEDIA_TYPES.get(format_name, "application/octet-stream")


def _read_input(plan: LoadedWorkflow) -> DNASet:
    config = ReadConfig(
        alphabet=DNAAlphabet(plan.spec.input.alphabet),
        compression=cast(Any, plan.spec.input.compression),
        max_records=plan.spec.limits.max_records,
        max_sequence_symbols=plan.spec.limits.max_sequence_length,
        max_input_bytes=plan.spec.limits.max_input_bytes,
    )
    records: list[DNARecord] = []
    total = 0
    with read(
        plan.input_path,
        format=plan.spec.input.format,
        config=config,
    ) as source:
        for index, record in enumerate(source):
            symbol_length = record.sequence.symbol_length
            coordinate_span = record.sequence.coordinate_span
            if symbol_length > plan.spec.limits.max_sequence_length or (
                coordinate_span is not None
                and coordinate_span > plan.spec.limits.max_sequence_length
            ):
                raise ConfigurationError(
                    "Input record exceeds limits.max_sequence_length.",
                    code="WORKFLOW_SEQUENCE_LIMIT",
                    context={
                        "index": index,
                        "record_id": record.id,
                        "symbol_length": symbol_length,
                        "coordinate_span": coordinate_span,
                        "maximum": plan.spec.limits.max_sequence_length,
                    },
                )
            total += symbol_length if coordinate_span is None else coordinate_span
            if total > plan.spec.limits.max_total_bases:
                raise ConfigurationError(
                    "Input exceeds limits.max_total_bases.",
                    code="WORKFLOW_TOTAL_BASE_LIMIT",
                    context={"observed": total, "maximum": plan.spec.limits.max_total_bases},
                )
            records.append(record)
    return DNASet.from_records(records, source=str(plan.input_path))


def _params(step: WorkflowStep) -> dict[str, object]:
    return _plain_dict(step.params)


def _with_result(
    state: _State,
    step_id: str,
    summary: Mapping[str, object],
    *,
    records: DNASet | None = None,
    split_result: SplitResult | Literal["preserve"] | None = "preserve",
) -> _State:
    results = dict(state.results)
    results[step_id] = _plain_dict(summary)
    resolved_split = state.split_result if split_result == "preserve" else split_result
    return _State(state.records if records is None else records, results, resolved_split)


def _check_result_size(value: object, spec: WorkflowSpec, *, step_id: str) -> None:
    size = len(_canonical_bytes(value))
    if size > spec.limits.max_result_bytes:
        raise ConfigurationError(
            "Workflow step result exceeds limits.max_result_bytes.",
            code="WORKFLOW_RESULT_LIMIT",
            context={
                "step_id": step_id,
                "byte_size": size,
                "maximum": spec.limits.max_result_bytes,
            },
        )


def _normalize_operation(
    step: WorkflowStep, state: _State, spec: WorkflowSpec
) -> _OperationOutcome:
    config = NormalizationConfig(**cast(dict[str, Any], _params(step)))
    records: list[DNARecord] = []
    audit: list[dict[str, object]] = []
    for record in state.records:
        result = normalize(record.sequence, config=config)
        if result.sequence is None:
            raise SequenceError(
                "Workflow normalization did not produce a valid DNA sequence.",
                code="WORKFLOW_NORMALIZATION_FAILED",
                context={
                    "record_id": record.id,
                    "issue_codes": [issue.code for issue in result.issues],
                },
            )
        if result.sequence.symbol_length != record.sequence.symbol_length and (
            record.features or record.letter_annotations
        ):
            raise SequenceError(
                "Length-changing normalization cannot preserve record annotations.",
                code="WORKFLOW_ANNOTATION_LENGTH_CHANGED",
                context={"record_id": record.id},
            )
        records.append(
            DNARecord(
                result.sequence,
                record.id,
                description=record.description,
                features=record.features,
                metadata=record.metadata,
                letter_annotations=record.letter_annotations,
            )
        )
        audit.append(
            {
                "record_id": record.id,
                "modified": result.was_modified,
                "change_count": len(result.changes),
                "issue_codes": [issue.code for issue in result.issues],
            }
        )
    normalized = DNASet(
        records,
        name=state.records.name,
        source=state.records.source,
        version=state.records.version,
        metadata=state.records.metadata,
    )
    summary = {
        "record_count": len(normalized),
        "modified_record_count": sum(bool(row["modified"]) for row in audit),
        "records": audit,
    }
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(
        _with_result(state, step.id, summary, records=normalized, split_result=None), summary
    )


def _validation_operation(
    step: WorkflowStep, state: _State, spec: WorkflowSpec
) -> _OperationOutcome:
    params = _params(step)
    fail_on_invalid = cast(bool, params.pop("fail_on_invalid", True))
    require_unique_ids = cast(bool, params.pop("require_unique_ids", True))
    collect_record_reports = cast(bool, params.pop("collect_record_reports", True))
    record_config = ValidationConfig(**cast(dict[str, Any], params))
    report = validate(
        state.records,
        config=DatasetValidationConfig(
            record_config,
            require_unique_ids=require_unique_ids,
            collect_record_reports=collect_record_reports,
        ),
    )
    invalid_ids = (
        [
            item.record_id
            for item in report.record_reports or ()
            if not item.is_valid and item.record_id is not None
        ]
        if collect_record_reports
        else []
    )
    summary = {
        "record_count": report.record_count,
        "is_valid": report.is_valid,
        "ids_unique": report.ids_unique,
        "duplicate_ids": [item.id for item in report.duplicate_ids],
        "invalid_record_ids": invalid_ids,
        "issue_codes": [item.code for item in report.issues],
    }
    if fail_on_invalid and not report.is_valid:
        raise ConfigurationError(
            "Workflow validation failed.",
            code="WORKFLOW_VALIDATION_FAILED",
            context=summary,
            hint="Set fail_on_invalid=false only when downstream handling is intentional.",
        )
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(_with_result(state, step.id, summary), summary)


def _descriptor_operation(
    step: WorkflowStep, state: _State, spec: WorkflowSpec
) -> _OperationOutcome:
    params = _params(step)
    metrics = cast(list[str], params.get("metrics", ["length", "gc"]))
    ambiguity_policy = cast(str, params.get("ambiguity_policy", "error"))
    if len(state.records) * len(metrics) * 256 > spec.limits.max_result_bytes:
        raise ConfigurationError(
            "Workflow descriptor table would exceed limits.max_result_bytes.",
            code="WORKFLOW_RESULT_LIMIT",
            context={"step_id": step.id, "record_count": len(state.records)},
        )
    rows: list[dict[str, object]] = []
    for record in state.records:
        row: dict[str, object] = {"record_id": record.id}
        for metric in metrics:
            if metric == "length":
                payload = length_features(record).to_dict()
            elif metric == "gc":
                payload = gc_at_content(record, ambiguity_policy=ambiguity_policy).to_dict()
            elif metric == "composition":
                payload = base_composition(record, ambiguity_policy=ambiguity_policy).to_dict()
            else:
                payload = shannon_entropy(record, ambiguity_policy=ambiguity_policy).to_dict()
            row[metric] = payload
        rows.append(row)
    summary = {"record_count": len(rows), "metrics": metrics, "rows": rows}
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(_with_result(state, step.id, summary), summary)


def _fingerprint_operation(
    step: WorkflowStep, state: _State, spec: WorkflowSpec
) -> _OperationOutcome:
    params = _params(step)
    params.setdefault("k", 3)
    params.setdefault("max_dimension", spec.limits.max_fingerprint_dimension)
    k = cast(int, params["k"])
    raw_dimension = 4**k
    canonical = cast(bool, params.get("canonical", False))
    dimension = (
        (raw_dimension + (0 if k % 2 else 4 ** (k // 2))) // 2 if canonical else raw_dimension
    )
    representation = cast(str, params.get("representation", "dense"))
    schema_bytes = dimension * (k + 4)
    value_bytes = (
        dimension * max(1, len(state.records)) * 2
        if representation == "dense"
        else sum(record.sequence.symbol_length for record in state.records) * (k + 16)
    )
    if schema_bytes + value_bytes > spec.limits.max_result_bytes:
        raise ConfigurationError(
            "Workflow fingerprint result would exceed limits.max_result_bytes.",
            code="WORKFLOW_RESULT_LIMIT",
            context={
                "step_id": step.id,
                "estimated_minimum_bytes": schema_bytes + value_bytes,
                "maximum": spec.limits.max_result_bytes,
            },
        )
    rows: list[dict[str, object]] = []
    feature_names: tuple[str, ...] | None = None
    for record in state.records:
        result = kmer_fingerprint(record, **cast(dict[str, Any], params))
        if feature_names is None:
            feature_names = result.feature_names
        rows.append(
            {
                "record_id": record.id,
                "schema_version": result.schema_version,
                "dimension": len(result.feature_names),
                "observation_count": result.observation_count,
                "values": result.values,
            }
        )
    summary = {
        "record_count": len(rows),
        "k": params["k"],
        "canonical": params.get("canonical", False),
        "mode": params.get("mode", "count"),
        "representation": params.get("representation", "dense"),
        "feature_names": () if feature_names is None else feature_names,
        "rows": rows,
    }
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(_with_result(state, step.id, summary), summary)


def _deduplicate_operation(
    step: WorkflowStep, state: _State, spec: WorkflowSpec
) -> _OperationOutcome:
    params = _params(step)
    equivalence = cast(str, params.pop("equivalence", "exact"))
    result = deduplicate(
        state.records,
        equivalence=cast(Any, equivalence),
        config=DeduplicationConfig(**cast(dict[str, Any], params)),
    )
    summary = result.to_dict()
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(
        _with_result(state, step.id, summary, records=result.records, split_result=None), summary
    )


def _split_operation(step: WorkflowStep, state: _State, spec: WorkflowSpec) -> _OperationOutcome:
    params = _params(step)
    result = split(
        state.records,
        config=SplitConfig(seed=spec.seed, **cast(dict[str, Any], params)),
    )
    summary = result.to_dict()
    _check_result_size(summary, spec, step_id=step.id)
    return _OperationOutcome(_with_result(state, step.id, summary, split_result=result), summary)


def _selected_records(state: _State, subset: object, *, step_id: str) -> DNASet:
    if subset is None:
        return state.records
    if not isinstance(subset, str):
        raise ConfigurationError("Workflow subset must be a string.")
    if state.split_result is None:
        raise ConfigurationError(
            "Workflow subset selection requires an upstream split step.",
            code="WORKFLOW_SPLIT_REQUIRED",
            context={"step_id": step_id, "subset": subset},
        )
    try:
        return state.split_result.get(subset)
    except KeyError as exc:
        raise ConfigurationError(
            "Workflow subset name is not present in the upstream split.",
            code="UNKNOWN_WORKFLOW_SUBSET",
            context={"step_id": step_id, "subset": subset},
        ) from exc


def _output_format(target: str, explicit: object) -> str:
    if isinstance(explicit, str):
        return explicit
    name = target[:-3] if target.lower().endswith(".gz") else target
    suffix = Path(name).suffix.lower().lstrip(".")
    return _FORMAT_ALIASES.get(suffix, suffix)


def _publish_stage(stage: Path, target: Path, *, overwrite: bool) -> bool:
    existed = target.exists()
    if existed and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing workflow output: {target}")
    if overwrite:
        os.replace(stage, target)
    else:
        os.link(stage, target)
        stage.unlink()
    return existed


def _write_operation(
    step: WorkflowStep,
    state: _State,
    spec: WorkflowSpec,
    output_root: Path,
    available_output_bytes: int,
) -> _OperationOutcome:
    params = _params(step)
    relative_target = cast(str, params["target"])
    records = _selected_records(state, params.get("subset"), step_id=step.id)
    format_name = _output_format(relative_target, params.get("format"))
    compression = cast(str, params.get("compression", "auto"))
    if compression == "auto":
        compression = "gzip" if relative_target.lower().endswith(".gz") else "none"
    preliminary_summary = {
        "target": relative_target,
        "format": format_name,
        "record_count": len(records),
        "byte_count": spec.limits.max_output_bytes,
        "sha256": "0" * 64,
        "overwritten": False,
        "subset": params.get("subset"),
    }
    _check_result_size({**state.results, step.id: preliminary_summary}, spec, step_id=step.id)
    target = _ensure_target(output_root, relative_target, create_parent=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.{step.id}.", suffix=".stage", dir=target.parent
    )
    os.close(descriptor)
    stage = Path(temporary_name)
    stage.unlink()
    try:
        try:
            result = write(
                records,
                stage,
                format=format_name,
                config=WriteConfig(
                    overwrite=False,
                    create_parents=False,
                    compression=cast(Any, compression),
                    compression_level=cast(int, params.get("compression_level", 6)),
                    line_width=cast(int, params.get("line_width", 80)),
                    max_output_bytes=available_output_bytes,
                ),
            )
        except InputFormatError as exc:
            if exc.code != "OUTPUT_BYTE_LIMIT_EXCEEDED":
                raise
            raise ConfigurationError(
                "Workflow outputs exceed the total limits.max_output_bytes.",
                code="WORKFLOW_OUTPUT_LIMIT",
                context={"step_id": step.id, "remaining": available_output_bytes},
            ) from exc
        byte_count = result.byte_count or 0
        if byte_count > available_output_bytes:
            raise ConfigurationError(
                "Workflow outputs exceed the total limits.max_output_bytes.",
                code="WORKFLOW_OUTPUT_LIMIT",
                context={
                    "step_id": step.id,
                    "byte_size": byte_count,
                    "remaining": available_output_bytes,
                },
            )
        overwritten = _publish_stage(stage, target, overwrite=spec.overwrite)
    finally:
        stage.unlink(missing_ok=True)
    artifact = artifact_from_path(
        target,
        media_type=_MEDIA_TYPES[format_name],
        schema_version="dnakit-io-v1",
    )
    summary = {
        "target": relative_target,
        "format": format_name,
        "record_count": len(records),
        "byte_count": artifact.byte_size,
        "sha256": artifact.sha256,
        "overwritten": overwritten,
        "subset": params.get("subset"),
    }
    return _OperationOutcome(_with_result(state, step.id, summary), summary, (artifact,))


def _report_operation(
    step: WorkflowStep,
    state: _State,
    spec: WorkflowSpec,
    output_root: Path,
    available_output_bytes: int,
) -> _OperationOutcome:
    params = _params(step)
    relative_target = cast(str, params["target"])
    records = _selected_records(state, params.get("subset"), step_id=step.id)
    max_result_bytes = cast(int, params.get("max_result_bytes", spec.limits.max_result_bytes))
    preliminary_summary = {
        "target": relative_target,
        "format": "html",
        "record_count": len(records),
        "result_names": tuple(sorted(state.results)),
        "byte_count": spec.limits.max_output_bytes,
        "sha256": "0" * 64,
        "overwritten": False,
        "subset": params.get("subset"),
    }
    _check_result_size({**state.results, step.id: preliminary_summary}, spec, step_id=step.id)
    target = _ensure_target(output_root, relative_target, create_parent=True)
    artifact = build_html_report(
        records,
        results=state.results,
        title=cast(str, params.get("title", "DNAKit report")),
        max_records=cast(int, params.get("max_records", spec.limits.max_records)),
        max_result_bytes=max_result_bytes,
        max_total_sequence_symbols=spec.limits.max_total_bases,
        max_output_bytes=available_output_bytes,
    )
    byte_size = len(artifact.html.encode("utf-8"))
    if byte_size > available_output_bytes:
        raise ConfigurationError(
            "Workflow outputs exceed the total limits.max_output_bytes.",
            code="WORKFLOW_OUTPUT_LIMIT",
            context={
                "step_id": step.id,
                "byte_size": byte_size,
                "remaining": available_output_bytes,
            },
        )
    saved = save_html_report(
        artifact,
        target,
        config=SaveConfig(overwrite=spec.overwrite, create_parents=True),
    )
    summary = {
        "target": relative_target,
        "format": "html",
        "record_count": len(records),
        "result_names": artifact.result_names,
        "byte_count": saved.target_artifact.byte_size,
        "sha256": saved.target_artifact.sha256,
        "overwritten": saved.overwritten,
        "subset": params.get("subset"),
    }
    return _OperationOutcome(
        _with_result(state, step.id, summary), summary, (saved.target_artifact,)
    )


def _execute_operation(
    step: WorkflowStep,
    state: _State,
    spec: WorkflowSpec,
    output_root: Path,
    available_output_bytes: int,
) -> _OperationOutcome:
    if step.operation == "normalize":
        return _normalize_operation(step, state, spec)
    if step.operation == "validate":
        return _validation_operation(step, state, spec)
    if step.operation == "descriptors":
        return _descriptor_operation(step, state, spec)
    if step.operation == "fingerprint":
        return _fingerprint_operation(step, state, spec)
    if step.operation == "deduplicate":
        return _deduplicate_operation(step, state, spec)
    if step.operation == "split":
        return _split_operation(step, state, spec)
    if step.operation == "write":
        return _write_operation(step, state, spec, output_root, available_output_bytes)
    return _report_operation(step, state, spec, output_root, available_output_bytes)


def _state_digest(state: _State, step: WorkflowStep) -> str:
    record_digest = hashlib.sha256()
    for record in state.records:
        encoded = _canonical_bytes(record)
        record_digest.update(len(encoded).to_bytes(8, "big"))
        record_digest.update(encoded)
    split_payload: object = None if state.split_result is None else state.split_result.to_dict()
    return _sha256(
        {
            "records_sha256": record_digest.hexdigest(),
            "record_count": len(state.records),
            "results_sha256": _sha256(state.results),
            "split_sha256": _sha256(split_payload),
            "step": step.to_dict(),
        }
    )


def _artifact_summary(artifact: ArtifactRef, output_root: Path) -> dict[str, object]:
    absolute = (Path.cwd() / artifact.relative_path).resolve()
    if not absolute.is_relative_to(output_root):
        raise ConfigurationError(
            "Workflow output artifact escaped output_dir.",
            code="UNSAFE_WORKFLOW_OUTPUT_PATH",
            context={"path": str(absolute)},
        )
    return {
        "path": absolute.relative_to(output_root).as_posix(),
        "media_type": artifact.media_type,
        "schema_version": artifact.schema_version,
        "sha256": artifact.sha256,
        "byte_size": artifact.byte_size,
    }


def _resume_entry_map(steps: object) -> dict[str, Mapping[str, object]]:
    result: dict[str, Mapping[str, object]] = {}
    if not isinstance(steps, list):
        return result
    for raw in steps:
        if not isinstance(raw, Mapping):
            continue
        step_id = raw.get("id")
        status = raw.get("status")
        if isinstance(step_id, str) and status in {"succeeded", "skipped"}:
            result[step_id] = cast(Mapping[str, object], raw)
    return result


def _load_resume_context(
    manifest_path: Path,
    *,
    spec_sha256: str,
    input_sha256: str,
) -> _ResumeContext:
    payload = load_manifest(manifest_path)
    resolved = payload.get("resolved_config")
    if not isinstance(resolved, dict) or resolved.get("schema_version") != _EXECUTION_SCHEMA:
        raise ConfigurationError(
            "Existing manifest is not a DNAKit workflow execution manifest.",
            code="WORKFLOW_RESUME_MISMATCH",
        )
    if resolved.get("spec_sha256") != spec_sha256:
        raise ConfigurationError(
            "Workflow configuration hash differs from the resumable manifest.",
            code="WORKFLOW_RESUME_MISMATCH",
        )
    execution = resolved.get("execution")
    if not isinstance(execution, dict) or execution.get("input_sha256") != input_sha256:
        raise ConfigurationError(
            "Workflow input hash differs from the resumable manifest.",
            code="WORKFLOW_RESUME_MISMATCH",
        )
    current = execution.get("steps", [])
    prior = execution.get("resume_candidates", [])
    entries = _resume_entry_map(prior)
    entries.update(_resume_entry_map(current))
    prior_steps = tuple(
        cast(Mapping[str, object], item)
        for item in [
            *(prior if isinstance(prior, list) else []),
            *(current if isinstance(current, list) else []),
        ]
        if isinstance(item, Mapping)
    )
    return _ResumeContext(entries, prior_steps)


def _verified_resume_artifacts(
    entry: Mapping[str, object],
    *,
    step: WorkflowStep,
    input_sha256: str,
    output_root: Path,
) -> tuple[tuple[ArtifactRef, ...], Mapping[str, object]]:
    if entry.get("operation") != step.operation or entry.get("input_sha256") != input_sha256:
        raise ConfigurationError(
            "A resumable step no longer has the same operation or input state hash.",
            code="WORKFLOW_RESUME_MISMATCH",
            context={"step_id": step.id},
        )
    raw_artifacts = entry.get("artifacts")
    if not isinstance(raw_artifacts, list) or len(raw_artifacts) != 1:
        raise ConfigurationError(
            "A resumable output step must have exactly one artifact audit.",
            code="WORKFLOW_RESUME_MISMATCH",
            context={"step_id": step.id},
        )
    expected_target = step.params.get("target")
    verified: list[ArtifactRef] = []
    for raw in raw_artifacts:
        if not isinstance(raw, Mapping):
            raise ConfigurationError(
                "A resumable artifact audit is invalid.", code="WORKFLOW_RESUME_MISMATCH"
            )
        path_value = raw.get("path")
        media_type = raw.get("media_type")
        schema_version = raw.get("schema_version")
        expected_hash = raw.get("sha256")
        expected_size = raw.get("byte_size")
        if not all(
            isinstance(value, str) and value
            for value in (path_value, media_type, schema_version, expected_hash)
        ):
            raise ConfigurationError(
                "A resumable artifact audit is incomplete.", code="WORKFLOW_RESUME_MISMATCH"
            )
        if (
            isinstance(expected_size, bool)
            or not isinstance(expected_size, int)
            or expected_size < 0
        ):
            raise ConfigurationError(
                "A resumable artifact byte size is invalid.",
                code="WORKFLOW_RESUME_MISMATCH",
                context={"step_id": step.id},
            )
        if path_value != expected_target:
            raise ConfigurationError(
                "A resumable artifact path differs from the configured output target.",
                code="WORKFLOW_RESUME_MISMATCH",
                context={
                    "step_id": step.id,
                    "configured_target": expected_target,
                    "manifest_target": path_value,
                },
            )
        target = _ensure_target(output_root, cast(str, path_value), create_parent=False)
        if target.stat().st_size != expected_size:
            raise ConfigurationError(
                "A resumable output artifact failed size validation.",
                code="WORKFLOW_RESUME_ARTIFACT_MISMATCH",
                context={
                    "step_id": step.id,
                    "path": cast(str, path_value),
                    "expected_size": expected_size,
                    "actual_size": target.stat().st_size,
                },
            )
        artifact = artifact_from_path(
            target,
            media_type=cast(str, media_type),
            schema_version=cast(str, schema_version),
            max_bytes=max(1, expected_size),
        )
        if artifact.sha256 != expected_hash or artifact.byte_size != expected_size:
            raise ConfigurationError(
                "A resumable output artifact failed checksum validation.",
                code="WORKFLOW_RESUME_ARTIFACT_MISMATCH",
                context={
                    "step_id": step.id,
                    "path": cast(str, path_value),
                    "expected_sha256": expected_hash,
                    "actual_sha256": artifact.sha256,
                },
            )
        verified.append(artifact)
    summary = entry.get("summary")
    if not isinstance(summary, Mapping):
        summary = {"resume": "verified"}
    return tuple(verified), cast(Mapping[str, object], summary)


def _workflow_provenance() -> Provenance:
    try:
        yaml_version = version("PyYAML")
    except PackageNotFoundError:
        yaml_version = "unknown"
    return Provenance(dependency_versions={"PyYAML": yaml_version})


def _resolved_manifest_config(
    plan: LoadedWorkflow,
    *,
    input_sha256: str,
    resume_requested: bool,
    steps: tuple[WorkflowStepResult, ...],
    resume_candidates: tuple[Mapping[str, object], ...],
) -> dict[str, object]:
    return {
        "schema_version": _EXECUTION_SCHEMA,
        "spec_sha256": plan.spec.sha256,
        "workflow": plan.spec.to_dict(),
        "runtime": {
            "resume_requested": resume_requested,
            "network_access": False,
            "arbitrary_command_execution": False,
        },
        "execution": {
            "input_sha256": input_sha256,
            "steps": [step.to_dict() for step in steps],
            "resume_candidates": [_plain_dict(item) for item in resume_candidates],
        },
    }


def _save_execution_manifest(
    plan: LoadedWorkflow,
    *,
    input_artifact: ArtifactRef,
    config_artifact: ArtifactRef,
    resume_requested: bool,
    steps: tuple[WorkflowStepResult, ...],
    resume_candidates: tuple[Mapping[str, object], ...],
    outputs: tuple[ArtifactRef, ...],
    issues: tuple[Issue, ...],
    status: Literal["running", "succeeded", "failed"],
    started_at: str,
) -> ArtifactRef:
    builder = RunManifestBuilder(
        plan.spec.run_id,
        ("dnakit", "workflow", str(plan.config_path)),
        _resolved_manifest_config(
            plan,
            input_sha256=input_artifact.sha256,
            resume_requested=resume_requested,
            steps=steps,
            resume_candidates=resume_candidates,
        ),
        seed=plan.spec.seed,
        seed_derivation="workflow master seed applied to deterministic dataset steps",
        provenance=_workflow_provenance(),
        started_at=started_at,
    )
    builder.add_input(config_artifact)
    builder.add_input(input_artifact)
    for artifact in outputs:
        builder.add_output(artifact)
    for issue in issues:
        builder.add_issue(issue)
    manifest = builder.build(status=status)
    target = plan.output_dir / _MANIFEST_NAME
    return save_manifest(manifest, target, overwrite=target.exists())


def _failed_issue(error: Exception, step_id: str | None) -> Issue:
    payload = _error_payload(error)
    details: dict[str, object] = {"error": payload}
    if step_id is not None:
        details["step_id"] = step_id
    return Issue(
        cast(str, payload["code"]),
        IssueSeverity.ERROR,
        f"Workflow execution failed: {error}",
        details=details,
    )


def _planned_result(step: WorkflowStep) -> WorkflowStepResult:
    return WorkflowStepResult(
        step.id,
        step.operation,
        step.input_ref,
        "planned",
        _params(step),
        summary={"execution": "not started"},
    )


def run_workflow(
    path: str | os.PathLike[str],
    *,
    dry_run: bool = False,
    resume: bool = False,
    progress: ProgressCallback | None = None,
) -> WorkflowRunResult:
    """Execute one bounded local JSON/YAML workflow without commands or network access.

    ``dry_run`` validates the complete schema, paths, input artifact, dependencies,
    and operation parameters without creating ``output_dir``.  ``resume`` reruns
    deterministic in-memory computation but skips a prior ``write`` or ``report``
    only after validating the workflow hash, input-state hash, output path, byte
    size, and SHA-256 digest.
    """

    if not isinstance(dry_run, bool) or not isinstance(resume, bool):
        raise ConfigurationError("dry_run and resume must be booleans.")
    if dry_run and resume:
        raise ConfigurationError(
            "dry_run and resume are mutually exclusive.",
            code="INVALID_WORKFLOW_MODE",
        )
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.")
    plan = load_workflow(path)
    _prepare_output_dir(plan.output_dir, create=False)
    config_artifact = artifact_from_path(
        plan.config_path,
        media_type=(
            "application/json" if plan.config_path.suffix.lower() == ".json" else "application/yaml"
        ),
        schema_version=_EXECUTION_SCHEMA,
        max_bytes=1_000_000,
    )
    if (
        config_artifact.sha256 != plan.config_file_sha256
        or config_artifact.byte_size != plan.config_file_size
    ):
        raise InputFormatError(
            "Workflow configuration changed after it was parsed.",
            code="WORKFLOW_CONFIG_CHANGED",
            context={"path": str(plan.config_path)},
        )
    input_artifact = artifact_from_path(
        plan.input_path,
        media_type=_input_media_type(plan),
        schema_version="dnakit-workflow-input-v1",
        max_bytes=plan.spec.limits.max_input_bytes,
    )
    if dry_run:
        for step in plan.spec.steps:
            if step.operation not in {"write", "report"}:
                continue
            target = _ensure_target(
                plan.output_dir,
                cast(str, step.params["target"]),
                create_parent=False,
            )
            if target.exists() and not plan.spec.overwrite:
                raise FileExistsError(f"Workflow output exists and overwrite is disabled: {target}")
        planned = tuple(_planned_result(step) for step in plan.spec.steps)
        for index, step in enumerate(plan.spec.steps, start=1):
            _emit(
                progress,
                run_id=plan.spec.run_id,
                step_id=step.id,
                operation=step.operation,
                index=index,
                total=len(plan.spec.steps),
                status="planned",
                message="validated and planned",
            )
        return WorkflowRunResult(
            plan.spec.run_id,
            "dry-run",
            True,
            False,
            plan.spec.sha256,
            input_artifact.sha256,
            str(plan.output_dir),
            None,
            None,
            planned,
            (),
        )

    _prepare_output_dir(plan.output_dir, create=True)
    manifest_path = plan.output_dir / _MANIFEST_NAME
    resume_context = _ResumeContext({}, ())
    if manifest_path.exists():
        if resume:
            resume_context = _load_resume_context(
                manifest_path,
                spec_sha256=plan.spec.sha256,
                input_sha256=input_artifact.sha256,
            )
        elif not plan.spec.overwrite:
            raise FileExistsError(
                "Workflow manifest already exists; enable resume or explicit overwrite: "
                f"{manifest_path}"
            )

    started_at = _utc_now()
    step_results: list[WorkflowStepResult] = []
    outputs: list[ArtifactRef] = []
    issues: list[Issue] = []
    manifest_artifact = _save_execution_manifest(
        plan,
        input_artifact=input_artifact,
        config_artifact=config_artifact,
        resume_requested=resume,
        steps=(),
        resume_candidates=resume_context.prior_steps,
        outputs=(),
        issues=(),
        status="running",
        started_at=started_at,
    )

    try:
        _emit(
            progress,
            run_id=plan.spec.run_id,
            step_id="input",
            operation="read",
            index=0,
            total=len(plan.spec.steps),
            status="started",
            message="reading input artifact",
        )
        records = _read_input(plan)
        verified_input = artifact_from_path(
            plan.input_path,
            media_type=_input_media_type(plan),
            schema_version="dnakit-workflow-input-v1",
            max_bytes=plan.spec.limits.max_input_bytes,
        )
        if (
            verified_input.sha256 != input_artifact.sha256
            or verified_input.byte_size != input_artifact.byte_size
        ):
            raise InputFormatError(
                "Workflow input changed between hashing and parsing.",
                code="WORKFLOW_INPUT_CHANGED",
                context={"path": str(plan.input_path)},
            )
        _emit(
            progress,
            run_id=plan.spec.run_id,
            step_id="input",
            operation="read",
            index=0,
            total=len(plan.spec.steps),
            status="succeeded",
            message=f"read {len(records)} records",
        )
    except Exception as error:
        issue = _failed_issue(error, None)
        issues.append(issue)
        manifest_artifact = _save_execution_manifest(
            plan,
            input_artifact=input_artifact,
            config_artifact=config_artifact,
            resume_requested=resume,
            steps=(),
            resume_candidates=resume_context.prior_steps,
            outputs=(),
            issues=tuple(issues),
            status="failed",
            started_at=started_at,
        )
        if plan.spec.error_policy == "raise":
            raise
        return WorkflowRunResult(
            plan.spec.run_id,
            "failed",
            False,
            resume,
            plan.spec.sha256,
            input_artifact.sha256,
            str(plan.output_dir),
            str(manifest_path),
            manifest_artifact,
            (),
            (),
        )

    states: dict[str, _State] = {"input": _State(records, {})}
    total_steps = len(plan.spec.steps)
    resumed_any = False
    for index, step in enumerate(plan.spec.steps, start=1):
        state = states.get(step.input_ref)
        if state is None:
            result = WorkflowStepResult(
                step.id,
                step.operation,
                step.input_ref,
                "blocked",
                _params(step),
                summary={"reason": "upstream step did not produce a state"},
                error={
                    "code": "WORKFLOW_UPSTREAM_FAILED",
                    "message": "The referenced upstream step failed or was blocked.",
                },
            )
            step_results.append(result)
            _best_effort_emit(
                progress,
                run_id=plan.spec.run_id,
                step_id=step.id,
                operation=step.operation,
                index=index,
                total=total_steps,
                status="blocked",
                message="upstream state unavailable",
            )
            continue

        step_started_at = _utc_now()
        timer = perf_counter()
        input_digest = _state_digest(state, step)
        try:
            _emit(
                progress,
                run_id=plan.spec.run_id,
                step_id=step.id,
                operation=step.operation,
                index=index,
                total=total_steps,
                status="started",
                message="step started",
            )
            resume_entry = resume_context.entries.get(step.id) if resume else None
            if resume_entry is not None and step.operation in {"write", "report"}:
                artifacts, summary = _verified_resume_artifacts(
                    resume_entry,
                    step=step,
                    input_sha256=input_digest,
                    output_root=plan.output_dir,
                )
                outcome = _OperationOutcome(
                    _with_result(state, step.id, summary), summary, artifacts
                )
                status: WorkflowStepStatus = "skipped"
                resumed_any = True
            else:
                used_output_bytes = sum(artifact.byte_size for artifact in outputs)
                outcome = _execute_operation(
                    step,
                    state,
                    plan.spec,
                    plan.output_dir,
                    plan.spec.limits.max_output_bytes - used_output_bytes,
                )
                status = "succeeded"
            if (
                sum(artifact.byte_size for artifact in outputs)
                + sum(artifact.byte_size for artifact in outcome.artifacts)
                > plan.spec.limits.max_output_bytes
            ):
                raise ConfigurationError(
                    "Workflow outputs exceed the total limits.max_output_bytes.",
                    code="WORKFLOW_OUTPUT_LIMIT",
                    context={"step_id": step.id},
                )
            _check_result_size(outcome.state.results, plan.spec, step_id=step.id)
            states[step.id] = outcome.state
            outputs.extend(outcome.artifacts)
            elapsed = perf_counter() - timer
            artifact_summaries = tuple(
                _artifact_summary(artifact, plan.output_dir) for artifact in outcome.artifacts
            )
            result = WorkflowStepResult(
                step.id,
                step.operation,
                step.input_ref,
                status,
                _params(step),
                started_at=step_started_at,
                finished_at=_utc_now(),
                duration_seconds=elapsed,
                input_sha256=input_digest,
                record_count=len(outcome.state.records),
                summary=outcome.summary,
                artifacts=artifact_summaries,
            )
            step_results.append(result)
            _emit(
                progress,
                run_id=plan.spec.run_id,
                step_id=step.id,
                operation=step.operation,
                index=index,
                total=total_steps,
                status=status,
                message="verified and skipped" if status == "skipped" else "step succeeded",
            )
            manifest_artifact = _save_execution_manifest(
                plan,
                input_artifact=input_artifact,
                config_artifact=config_artifact,
                resume_requested=resume,
                steps=tuple(step_results),
                resume_candidates=resume_context.prior_steps,
                outputs=tuple(outputs),
                issues=tuple(issues),
                status="running",
                started_at=started_at,
            )
        except Exception as error:
            elapsed = perf_counter() - timer
            payload = _error_payload(error)
            result = WorkflowStepResult(
                step.id,
                step.operation,
                step.input_ref,
                "failed",
                _params(step),
                started_at=step_started_at,
                finished_at=_utc_now(),
                duration_seconds=elapsed,
                input_sha256=input_digest,
                record_count=len(state.records),
                error=payload,
            )
            step_results.append(result)
            issues.append(_failed_issue(error, step.id))
            _best_effort_emit(
                progress,
                run_id=plan.spec.run_id,
                step_id=step.id,
                operation=step.operation,
                index=index,
                total=total_steps,
                status="failed",
                message=str(error),
            )
            if plan.spec.error_policy == "raise":
                manifest_artifact = _save_execution_manifest(
                    plan,
                    input_artifact=input_artifact,
                    config_artifact=config_artifact,
                    resume_requested=resume,
                    steps=tuple(step_results),
                    resume_candidates=resume_context.prior_steps,
                    outputs=tuple(outputs),
                    issues=tuple(issues),
                    status="failed",
                    started_at=started_at,
                )
                raise

    failed = any(result.status in {"failed", "blocked"} for result in step_results)
    final_status: WorkflowRunStatus = "failed" if failed else "succeeded"
    manifest_artifact = _save_execution_manifest(
        plan,
        input_artifact=input_artifact,
        config_artifact=config_artifact,
        resume_requested=resume,
        steps=tuple(step_results),
        resume_candidates=resume_context.prior_steps if failed else (),
        outputs=tuple(outputs),
        issues=tuple(issues),
        status="failed" if failed else "succeeded",
        started_at=started_at,
    )
    return WorkflowRunResult(
        plan.spec.run_id,
        final_status,
        False,
        resumed_any,
        plan.spec.sha256,
        input_artifact.sha256,
        str(plan.output_dir),
        str(manifest_path),
        manifest_artifact,
        tuple(step_results),
        tuple(outputs),
    )


__all__ = [
    "ProgressCallback",
    "WorkflowProgress",
    "WorkflowRunResult",
    "WorkflowStepResult",
    "run_workflow",
]
