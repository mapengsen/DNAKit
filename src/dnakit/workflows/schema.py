"""Strict, non-executable schema for local DNAKit workflows."""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast

import yaml  # type: ignore[import-untyped]

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.datasets import DeduplicationConfig, SplitConfig
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.standardize import (
    DatasetValidationConfig,
    NormalizationConfig,
    ValidationConfig,
)

WorkflowErrorPolicy = Literal["raise", "collect"]
WorkflowOperation = Literal[
    "normalize",
    "validate",
    "descriptors",
    "fingerprint",
    "deduplicate",
    "split",
    "write",
    "report",
]

_SCHEMA_VERSION = "dnakit-workflow-v1"
_MAX_CONFIG_BYTES = 1_000_000
_MAX_CONFIG_DEPTH = 32
_MAX_CONFIG_NODES = 100_000
_SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_OPERATIONS: frozenset[str] = frozenset(
    {
        "normalize",
        "validate",
        "descriptors",
        "fingerprint",
        "deduplicate",
        "split",
        "write",
        "report",
    }
)
_INPUT_FORMATS = frozenset({"fasta", "fastq", "csv", "tsv", "json", "jsonl", "genbank"})
_OUTPUT_FORMATS = _INPUT_FORMATS
_RESERVED_OUTPUTS = frozenset({".dnakit-workflow-output-v1", "run-manifest.json"})

_STEP_PARAMETERS: Mapping[str, frozenset[str]] = {
    "normalize": frozenset(
        {
            "alphabet",
            "uppercase",
            "remove_whitespace",
            "remove_invisible",
            "removable_separators",
            "keep_ambiguous",
            "keep_u",
            "keep_other",
            "u_policy",
            "ambiguity_policy",
            "ambiguity_mask",
            "allow_gaps",
            "raise_on_error",
        }
    ),
    "validate": frozenset(
        {
            "alphabet",
            "allow_empty",
            "allow_gaps",
            "allow_unknown_gap_length",
            "min_length",
            "max_length",
            "sequence_length",
            "max_ambiguity_fraction",
            "ambiguity_denominator_includes_gap",
            "required_metadata_fields",
            "required_letter_annotations",
            "check_phred_quality",
            "minimum_phred",
            "maximum_phred",
            "minimum_mean_phred",
            "require_unique_ids",
            "collect_record_reports",
            "fail_on_invalid",
        }
    ),
    "descriptors": frozenset({"metrics", "ambiguity_policy"}),
    "fingerprint": frozenset(
        {
            "k",
            "canonical",
            "mode",
            "representation",
            "ambiguity_policy",
            "overlapping",
            "cross_gaps",
            "max_dimension",
        }
    ),
    "deduplicate": frozenset(
        {
            "equivalence",
            "representative_policy",
            "conflict_field",
            "conflict_policy",
            "merge_metadata",
        }
    ),
    "split": frozenset(
        {
            "method",
            "ratios",
            "shuffle",
            "preserve_order",
            "metadata_key",
            "missing_metadata_policy",
            "similarity_k",
            "similarity_threshold",
            "similarity_ambiguity_policy",
            "similarity_gap_policy",
            "max_pairwise_records",
        }
    ),
    "write": frozenset(
        {"target", "format", "subset", "compression", "compression_level", "line_width"}
    ),
    "report": frozenset({"target", "subset", "title", "max_records", "max_result_bytes"}),
}


def _configuration_error(
    message: str,
    *,
    code: str = "INVALID_WORKFLOW_CONFIG",
    context: Mapping[str, object] | None = None,
) -> ConfigurationError:
    return ConfigurationError(message, code=code, context=context)


def _require_mapping(value: object, *, field: str) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise _configuration_error(f"Workflow {field} must be an object.", context={"field": field})
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str) or not key:
            raise _configuration_error(
                f"Workflow {field} keys must be non-empty strings.",
                context={"field": field},
            )
        result[key] = item
    return result


def _reject_unknown(values: Mapping[str, object], allowed: frozenset[str], *, field: str) -> None:
    unknown = sorted(set(values) - allowed)
    if unknown:
        raise _configuration_error(
            f"Workflow {field} contains unknown fields.",
            code="UNKNOWN_WORKFLOW_FIELD",
            context={"field": field, "unknown": unknown},
        )


def _require_string(value: object, *, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _configuration_error(
            f"Workflow {field} must be a non-empty string.", context={"field": field}
        )
    return value


def _require_bool(value: object, *, field: str) -> bool:
    if not isinstance(value, bool):
        raise _configuration_error(f"Workflow {field} must be a boolean.", context={"field": field})
    return value


def _positive_int(value: object, *, field: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0 or value > maximum:
        raise _configuration_error(
            f"Workflow {field} must be an integer in [1, {maximum}].",
            code="INVALID_WORKFLOW_LIMIT",
            context={"field": field, "value": value},
        )
    return value


def _safe_relative_output(value: object, *, field: str) -> str:
    text = _require_string(value, field=field)
    path = Path(text)
    if (
        "\x00" in text
        or "\\" in text
        or ":" in text
        or len(text) > 1_024
        or path.is_absolute()
        or text in {".", ".."}
        or any(part in {"", ".", ".."} for part in path.parts)
        or any(len(part) > 255 for part in path.parts)
        or path.name in _RESERVED_OUTPUTS
    ):
        raise _configuration_error(
            f"Workflow {field} must be a safe relative path below output_dir.",
            code="UNSAFE_WORKFLOW_OUTPUT_PATH",
            context={"field": field, "path": text},
        )
    return path.as_posix()


def _plain_mapping(values: Mapping[str, object]) -> dict[str, object]:
    return cast(dict[str, object], to_json_compatible(values))


@dataclass(frozen=True, slots=True)
class WorkflowLimits:
    """Hard resource ceilings applied before or during local execution."""

    max_steps: int = 50
    max_records: int = 100_000
    max_sequence_length: int = 10_000_000
    max_total_bases: int = 100_000_000
    max_input_bytes: int = 250_000_000
    max_fingerprint_dimension: int = 1_000_000
    max_result_bytes: int = 20_000_000
    max_output_bytes: int = 250_000_000

    def __post_init__(self) -> None:
        maxima = {
            "max_steps": 1_000,
            "max_records": 1_000_000,
            "max_sequence_length": 100_000_000,
            "max_total_bases": 1_000_000_000,
            "max_input_bytes": 2_000_000_000,
            "max_fingerprint_dimension": 1_000_000,
            "max_result_bytes": 100_000_000,
            "max_output_bytes": 2_000_000_000,
        }
        for field, maximum in maxima.items():
            _positive_int(getattr(self, field), field=f"limits.{field}", maximum=maximum)

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class WorkflowInput:
    """One explicit local sequence input artifact."""

    path: str
    format: str | None = None
    alphabet: str = "iupac"
    compression: str = "auto"

    def __post_init__(self) -> None:
        _require_string(self.path, field="input.path")
        if "\x00" in self.path:
            raise _configuration_error("Workflow input.path must not contain NUL.")
        if self.format is not None and self.format not in _INPUT_FORMATS:
            raise _configuration_error(
                "Workflow input.format is unsupported.",
                code="UNSUPPORTED_WORKFLOW_INPUT_FORMAT",
                context={"format": self.format, "allowed": sorted(_INPUT_FORMATS)},
            )
        if self.alphabet not in {"strict", "iupac"}:
            raise _configuration_error("Workflow input.alphabet must be 'strict' or 'iupac'.")
        if self.compression not in {"auto", "none", "gzip"}:
            raise _configuration_error("Workflow input.compression must be auto, none, or gzip.")


@dataclass(frozen=True, slots=True)
class WorkflowStep:
    """One whitelisted step referencing only an earlier workflow value."""

    id: str
    operation: WorkflowOperation
    input_ref: str
    params: FrozenDict

    def __init__(
        self,
        id: str,
        operation: WorkflowOperation | str,
        input_ref: str,
        params: Mapping[str, object] | None = None,
    ) -> None:
        if not isinstance(id, str) or _SAFE_IDENTIFIER.fullmatch(id) is None:
            raise _configuration_error(
                "Workflow step id must be a safe identifier of at most 64 characters.",
                context={"step_id": id},
            )
        if operation not in _OPERATIONS:
            raise _configuration_error(
                "Workflow step operation is not whitelisted.",
                code="UNKNOWN_WORKFLOW_OPERATION",
                context={"step_id": id, "operation": operation},
            )
        if not isinstance(input_ref, str) or (
            input_ref != "input" and _SAFE_IDENTIFIER.fullmatch(input_ref) is None
        ):
            raise _configuration_error(
                "Workflow step input must reference 'input' or a safe earlier step id.",
                context={"step_id": id, "input": input_ref},
            )
        resolved_params = (
            {} if params is None else _require_mapping(params, field=f"steps.{id}.params")
        )
        _reject_unknown(
            resolved_params,
            _STEP_PARAMETERS[operation],
            field=f"steps.{id}.params",
        )
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "operation", cast(WorkflowOperation, operation))
        object.__setattr__(self, "input_ref", input_ref)
        object.__setattr__(self, "params", freeze_mapping(resolved_params))

    def to_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "operation": self.operation,
            "input": self.input_ref,
            "params": _plain_mapping(self.params),
        }


@dataclass(frozen=True, slots=True)
class WorkflowSpec:
    """Fully resolved and validated workflow definition."""

    schema_version: str
    run_id: str
    input: WorkflowInput
    output_dir: str
    seed: int
    error_policy: WorkflowErrorPolicy
    overwrite: bool
    limits: WorkflowLimits
    steps: tuple[WorkflowStep, ...]

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA_VERSION:
            raise _configuration_error(
                "Unsupported workflow schema_version.",
                code="UNSUPPORTED_WORKFLOW_SCHEMA",
                context={"schema_version": self.schema_version, "expected": _SCHEMA_VERSION},
            )
        if not isinstance(self.run_id, str) or _SAFE_IDENTIFIER.fullmatch(self.run_id) is None:
            raise _configuration_error(
                "Workflow run_id must be a safe identifier of at most 64 characters."
            )
        if not isinstance(self.input, WorkflowInput):
            raise _configuration_error("Workflow input must be WorkflowInput.")
        _safe_relative_output(self.output_dir, field="output_dir")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise _configuration_error("Workflow seed must be an integer.")
        if self.error_policy not in {"raise", "collect"}:
            raise _configuration_error("Workflow error_policy must be 'raise' or 'collect'.")
        if not isinstance(self.overwrite, bool):
            raise _configuration_error("Workflow overwrite must be a boolean.")
        if not isinstance(self.limits, WorkflowLimits):
            raise _configuration_error("Workflow limits must be WorkflowLimits.")
        if (
            not isinstance(self.steps, tuple)
            or any(not isinstance(step, WorkflowStep) for step in self.steps)
            or not self.steps
            or len(self.steps) > self.limits.max_steps
        ):
            raise _configuration_error(
                "Workflow steps must contain WorkflowStep objects and stay within "
                "limits.max_steps.",
                code="WORKFLOW_STEP_LIMIT",
                context={"count": len(self.steps), "max_steps": self.limits.max_steps},
            )
        known = {"input"}
        targets: set[str] = set()
        for step in self.steps:
            if step.id in known:
                raise _configuration_error(
                    "Workflow step ids must be unique and cannot be 'input'.",
                    code="DUPLICATE_WORKFLOW_STEP",
                    context={"step_id": step.id},
                )
            if step.input_ref not in known:
                raise _configuration_error(
                    "Workflow steps may reference only the input or an earlier step.",
                    code="UNKNOWN_WORKFLOW_INPUT_REFERENCE",
                    context={"step_id": step.id, "input": step.input_ref},
                )
            known.add(step.id)
            if step.operation in {"write", "report"}:
                target = _safe_relative_output(
                    step.params.get("target"), field=f"steps.{step.id}.params.target"
                )
                if target in targets:
                    raise _configuration_error(
                        "Workflow output targets must be unique.",
                        code="DUPLICATE_WORKFLOW_OUTPUT",
                        context={"target": target},
                    )
                targets.add(target)
        _validate_step_parameters(self.steps, seed=self.seed, limits=self.limits)

    @property
    def sha256(self) -> str:
        encoded = json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "input": cast(dict[str, object], to_json_compatible(self.input)),
            "output_dir": self.output_dir,
            "seed": self.seed,
            "error_policy": self.error_policy,
            "overwrite": self.overwrite,
            "limits": self.limits.to_dict(),
            "steps": [step.to_dict() for step in self.steps],
        }


@dataclass(frozen=True, slots=True)
class LoadedWorkflow:
    """Validated workflow plus resolved local paths."""

    spec: WorkflowSpec
    config_path: Path
    input_path: Path
    output_dir: Path
    config_file_sha256: str
    config_file_size: int


def _validate_shape(value: object) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(value, 0)]
    while stack:
        item, depth = stack.pop()
        nodes += 1
        if nodes > _MAX_CONFIG_NODES or depth > _MAX_CONFIG_DEPTH:
            raise InputFormatError(
                "Workflow configuration exceeds structural limits.",
                code="WORKFLOW_CONFIG_STRUCTURE_LIMIT",
                context={
                    "max_nodes": _MAX_CONFIG_NODES,
                    "max_depth": _MAX_CONFIG_DEPTH,
                },
            )
        if isinstance(item, Mapping):
            stack.extend((key, depth + 1) for key in item)
            stack.extend((child, depth + 1) for child in item.values())
        elif isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray)):
            stack.extend((child, depth + 1) for child in item)


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise InputFormatError(
                "Workflow JSON contains a duplicate object key.",
                code="DUPLICATE_WORKFLOW_FIELD",
                context={"key": key},
            )
        result[key] = value
    return result


class _UniqueKeySafeLoader(yaml.SafeLoader):  # type: ignore[misc]
    """SafeLoader variant that rejects duplicate mapping keys."""


def _construct_unique_mapping(
    loader: yaml.SafeLoader, node: yaml.nodes.MappingNode, deep: bool = False
) -> dict[object, object]:
    loader.flatten_mapping(node)
    result: dict[object, object] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in result
        except TypeError as exc:
            raise InputFormatError(
                "Workflow YAML mapping keys must be hashable strings.",
                code="INVALID_WORKFLOW_CONFIG",
            ) from exc
        if duplicate:
            raise InputFormatError(
                "Workflow YAML contains a duplicate mapping key.",
                code="DUPLICATE_WORKFLOW_FIELD",
                context={"key": str(key)},
            )
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueKeySafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping
)


def _parse_payload(path: Path) -> tuple[dict[str, object], bytes]:
    if path.suffix.lower() not in {".json", ".yaml", ".yml"}:
        raise InputFormatError(
            "Workflow configuration must use JSON or YAML.",
            code="UNSUPPORTED_WORKFLOW_CONFIG_FORMAT",
            context={"path": str(path)},
        )
    if path.is_symlink() or not path.is_file():
        raise InputFormatError(
            "Workflow configuration must be a regular non-symlink file.",
            code="INVALID_WORKFLOW_CONFIG_PATH",
            context={"path": str(path)},
        )
    try:
        with path.open("rb") as handle:
            payload = handle.read(_MAX_CONFIG_BYTES + 1)
        if len(payload) > _MAX_CONFIG_BYTES:
            raise InputFormatError(
                "Workflow configuration exceeds the byte limit.",
                code="WORKFLOW_CONFIG_SIZE_LIMIT",
                context={"byte_size": len(payload), "max_bytes": _MAX_CONFIG_BYTES},
            )
        text = payload.decode("utf-8")
        if path.suffix.lower() == ".json":
            parsed: object = json.loads(text, object_pairs_hook=_unique_json_object)
        else:
            if any(
                isinstance(token, yaml.tokens.AliasToken)
                for token in yaml.scan(text, Loader=_UniqueKeySafeLoader)
            ):
                raise InputFormatError(
                    "Workflow YAML aliases are disabled.",
                    code="WORKFLOW_YAML_ALIAS_DISABLED",
                )
            parsed = yaml.load(text, Loader=_UniqueKeySafeLoader)
    except InputFormatError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, yaml.YAMLError, RecursionError) as exc:
        raise InputFormatError(
            "Workflow configuration could not be parsed.",
            code="INVALID_WORKFLOW_CONFIG",
            context={"path": str(path)},
        ) from exc
    _validate_shape(parsed)
    return _require_mapping(parsed, field="root"), payload


def _parse_limits(value: object) -> WorkflowLimits:
    if value is None:
        return WorkflowLimits()
    values = _require_mapping(value, field="limits")
    allowed = frozenset(
        {
            "max_steps",
            "max_records",
            "max_sequence_length",
            "max_total_bases",
            "max_input_bytes",
            "max_fingerprint_dimension",
            "max_result_bytes",
            "max_output_bytes",
        }
    )
    _reject_unknown(values, allowed, field="limits")
    defaults = WorkflowLimits()
    return WorkflowLimits(
        max_steps=cast(int, values.get("max_steps", defaults.max_steps)),
        max_records=cast(int, values.get("max_records", defaults.max_records)),
        max_sequence_length=cast(
            int, values.get("max_sequence_length", defaults.max_sequence_length)
        ),
        max_total_bases=cast(int, values.get("max_total_bases", defaults.max_total_bases)),
        max_input_bytes=cast(int, values.get("max_input_bytes", defaults.max_input_bytes)),
        max_fingerprint_dimension=cast(
            int,
            values.get("max_fingerprint_dimension", defaults.max_fingerprint_dimension),
        ),
        max_result_bytes=cast(int, values.get("max_result_bytes", defaults.max_result_bytes)),
        max_output_bytes=cast(int, values.get("max_output_bytes", defaults.max_output_bytes)),
    )


def _parse_input(value: object) -> WorkflowInput:
    values = _require_mapping(value, field="input")
    _reject_unknown(values, frozenset({"path", "format", "alphabet", "compression"}), field="input")
    if "path" not in values:
        raise _configuration_error("Workflow input.path is required.")
    format_value = values.get("format")
    if format_value is not None and not isinstance(format_value, str):
        raise _configuration_error("Workflow input.format must be a string or null.")
    return WorkflowInput(
        path=_require_string(values["path"], field="input.path"),
        format=format_value,
        alphabet=cast(str, values.get("alphabet", "iupac")),
        compression=cast(str, values.get("compression", "auto")),
    )


def _parse_steps(value: object) -> tuple[WorkflowStep, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise _configuration_error("Workflow steps must be an array.")
    steps: list[WorkflowStep] = []
    for index, raw_step in enumerate(value):
        values = _require_mapping(raw_step, field=f"steps[{index}]")
        _reject_unknown(
            values,
            frozenset({"id", "operation", "input", "params"}),
            field=f"steps[{index}]",
        )
        missing = sorted({"id", "operation", "input"} - values.keys())
        if missing:
            raise _configuration_error(
                "Workflow step is missing required fields.",
                context={"index": index, "missing": missing},
            )
        operation = _require_string(values["operation"], field=f"steps[{index}].operation")
        steps.append(
            WorkflowStep(
                _require_string(values["id"], field=f"steps[{index}].id"),
                operation,
                _require_string(values["input"], field=f"steps[{index}].input"),
                _require_mapping(values.get("params", {}), field=f"steps[{index}].params"),
            )
        )
    return tuple(steps)


def _validate_step_parameters(
    steps: tuple[WorkflowStep, ...], *, seed: int, limits: WorkflowLimits
) -> None:
    descriptor_metrics = {"length", "gc", "composition", "entropy"}
    for step in steps:
        params = _plain_mapping(step.params)
        if step.operation == "normalize":
            NormalizationConfig(**cast(dict[str, Any], params))
        elif step.operation == "validate":
            record_params = {
                key: value
                for key, value in params.items()
                if key not in {"require_unique_ids", "collect_record_reports", "fail_on_invalid"}
            }
            record_config = ValidationConfig(**cast(dict[str, Any], record_params))
            DatasetValidationConfig(
                record_config,
                require_unique_ids=cast(bool, params.get("require_unique_ids", True)),
                collect_record_reports=cast(bool, params.get("collect_record_reports", True)),
            )
            if "fail_on_invalid" in params:
                _require_bool(params["fail_on_invalid"], field=f"steps.{step.id}.fail_on_invalid")
        elif step.operation == "descriptors":
            metrics = params.get("metrics", ["length", "gc"])
            if not isinstance(metrics, Sequence) or isinstance(metrics, (str, bytes, bytearray)):
                raise _configuration_error("Descriptor metrics must be an array of names.")
            metric_names = tuple(metrics)
            if not metric_names or any(
                not isinstance(metric, str) or metric not in descriptor_metrics
                for metric in metric_names
            ):
                raise _configuration_error(
                    "Descriptor metrics contain an unsupported name.",
                    context={"step_id": step.id, "allowed": sorted(descriptor_metrics)},
                )
            if len(set(metric_names)) != len(metric_names):
                raise _configuration_error("Descriptor metrics must not contain duplicates.")
            ambiguity = params.get("ambiguity_policy", "error")
            if ambiguity not in {"error", "ignore"}:
                raise _configuration_error("Descriptor ambiguity_policy is invalid.")
        elif step.operation == "fingerprint":
            k = _positive_int(params.get("k", 3), field=f"steps.{step.id}.k", maximum=30)
            for flag in ("canonical", "overlapping", "cross_gaps"):
                if flag in params:
                    _require_bool(params[flag], field=f"steps.{step.id}.{flag}")
            if params.get("mode", "count") not in {"count", "binary", "presence", "frequency"}:
                raise _configuration_error("Fingerprint mode is invalid.")
            if params.get("representation", "dense") not in {"dense", "sparse"}:
                raise _configuration_error("Fingerprint representation is invalid.")
            if params.get("ambiguity_policy", "error") not in {"error", "ignore"}:
                raise _configuration_error("Fingerprint ambiguity_policy is invalid.")
            requested_dimension = _positive_int(
                params.get("max_dimension", limits.max_fingerprint_dimension),
                field=f"steps.{step.id}.max_dimension",
                maximum=limits.max_fingerprint_dimension,
            )
            canonical = cast(bool, params.get("canonical", False))
            raw_dimension = 4**k
            schema_dimension = (
                (raw_dimension + (0 if k % 2 else 4 ** (k // 2))) // 2
                if canonical
                else raw_dimension
            )
            if schema_dimension > requested_dimension:
                raise _configuration_error(
                    "Fingerprint schema exceeds the configured dimension limit.",
                    code="WORKFLOW_FINGERPRINT_LIMIT",
                    context={
                        "step_id": step.id,
                        "dimension": schema_dimension,
                        "max_dimension": requested_dimension,
                    },
                )
        elif step.operation == "deduplicate":
            equivalence = params.pop("equivalence", "exact")
            if equivalence not in {
                "exact",
                "reverse_complement",
                "circular",
                "circular_reverse_complement",
            }:
                raise _configuration_error("Deduplication equivalence is invalid.")
            DeduplicationConfig(**cast(dict[str, Any], params))
        elif step.operation == "split":
            SplitConfig(seed=seed, **cast(dict[str, Any], params))
        elif step.operation == "write":
            target = _safe_relative_output(params.get("target"), field=f"steps.{step.id}.target")
            format_value = params.get("format")
            if format_value is not None and format_value not in _OUTPUT_FORMATS:
                raise _configuration_error(
                    "Workflow write format is unsupported.",
                    code="UNSUPPORTED_WORKFLOW_OUTPUT_FORMAT",
                )
            if format_value is None:
                name = target[:-3] if target.lower().endswith(".gz") else target
                suffix = Path(name).suffix.lower().lstrip(".")
                aliases = {
                    "fa": "fasta",
                    "fna": "fasta",
                    "fq": "fastq",
                    "gb": "genbank",
                    "gbk": "genbank",
                    "ndjson": "jsonl",
                }
                if aliases.get(suffix, suffix) not in _OUTPUT_FORMATS:
                    raise _configuration_error(
                        "Workflow write target does not identify a supported format.",
                        code="UNSUPPORTED_WORKFLOW_OUTPUT_FORMAT",
                    )
            subset = params.get("subset")
            if subset is not None:
                _require_string(subset, field=f"steps.{step.id}.subset")
            if params.get("compression", "auto") not in {"auto", "none", "gzip"}:
                raise _configuration_error("Workflow write compression is invalid.")
            _positive_int(
                params.get("line_width", 80),
                field=f"steps.{step.id}.line_width",
                maximum=1_000_000,
            )
            compression_level = params.get("compression_level", 6)
            if (
                isinstance(compression_level, bool)
                or not isinstance(compression_level, int)
                or not 0 <= compression_level <= 9
            ):
                raise _configuration_error("Workflow write compression_level must be in [0, 9].")
        else:
            target = _safe_relative_output(params.get("target"), field=f"steps.{step.id}.target")
            if Path(target).suffix.lower() not in {".html", ".htm"}:
                raise _configuration_error(
                    "Workflow report target must end with .html or .htm.",
                    code="UNSUPPORTED_WORKFLOW_OUTPUT_FORMAT",
                )
            subset = params.get("subset")
            if subset is not None:
                _require_string(subset, field=f"steps.{step.id}.subset")
            title = _require_string(
                params.get("title", "DNAKit report"), field=f"steps.{step.id}.title"
            )
            if len(title) > 500:
                raise _configuration_error("Workflow report title is too long.")
            max_records = _positive_int(
                params.get("max_records", limits.max_records),
                field=f"steps.{step.id}.max_records",
                maximum=limits.max_records,
            )
            del max_records
            _positive_int(
                params.get("max_result_bytes", limits.max_result_bytes),
                field=f"steps.{step.id}.max_result_bytes",
                maximum=limits.max_result_bytes,
            )


def _resolve_output_dir(config_path: Path, raw: str) -> Path:
    normalized = _safe_relative_output(raw, field="output_dir")
    base = config_path.parent.resolve()
    candidate = base.joinpath(*Path(normalized).parts)
    resolved = candidate.resolve(strict=False)
    if not resolved.is_relative_to(base) or resolved == base:
        raise _configuration_error(
            "Workflow output_dir escapes the configuration directory.",
            code="UNSAFE_WORKFLOW_OUTPUT_DIR",
            context={"output_dir": raw},
        )
    current = base
    for part in Path(normalized).parts:
        current /= part
        if current.exists() and current.is_symlink():
            raise _configuration_error(
                "Workflow output_dir must not traverse symbolic links.",
                code="UNSAFE_WORKFLOW_OUTPUT_DIR",
                context={"path": str(current)},
            )
    return resolved


def load_workflow(path: str | os.PathLike[str]) -> LoadedWorkflow:
    """Load JSON/YAML without evaluating commands, Python, tags, or templates."""

    try:
        config_path = Path(path).expanduser().absolute()
    except TypeError as exc:
        raise _configuration_error("Workflow path must be a filesystem path.") from exc
    values, config_bytes = _parse_payload(config_path)
    _reject_unknown(
        values,
        frozenset(
            {
                "schema_version",
                "run_id",
                "input",
                "output_dir",
                "seed",
                "error_policy",
                "overwrite",
                "limits",
                "steps",
            }
        ),
        field="root",
    )
    missing = sorted({"schema_version", "run_id", "input", "output_dir", "steps"} - values.keys())
    if missing:
        raise _configuration_error(
            "Workflow configuration is missing required fields.", context={"missing": missing}
        )
    limits = _parse_limits(values.get("limits"))
    spec = WorkflowSpec(
        schema_version=_require_string(values["schema_version"], field="schema_version"),
        run_id=_require_string(values["run_id"], field="run_id"),
        input=_parse_input(values["input"]),
        output_dir=_safe_relative_output(values["output_dir"], field="output_dir"),
        seed=cast(int, values.get("seed", 0)),
        error_policy=cast(WorkflowErrorPolicy, values.get("error_policy", "raise")),
        overwrite=_require_bool(values.get("overwrite", False), field="overwrite"),
        limits=limits,
        steps=_parse_steps(values["steps"]),
    )
    input_path = Path(spec.input.path).expanduser()
    if not input_path.is_absolute():
        input_path = config_path.parent / input_path
    input_path = input_path.absolute()
    return LoadedWorkflow(
        spec=spec,
        config_path=config_path,
        input_path=input_path,
        output_dir=_resolve_output_dir(config_path, spec.output_dir),
        config_file_sha256=hashlib.sha256(config_bytes).hexdigest(),
        config_file_size=len(config_bytes),
    )


__all__ = [
    "LoadedWorkflow",
    "WorkflowErrorPolicy",
    "WorkflowInput",
    "WorkflowLimits",
    "WorkflowOperation",
    "WorkflowSpec",
    "WorkflowStep",
    "load_workflow",
]
