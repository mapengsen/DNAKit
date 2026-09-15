"""Configuration and immutable schemas for user-defined property models."""

from __future__ import annotations

import math
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal, TypeAlias
from urllib.parse import urlparse

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError
from dnakit.representations import get_embedding_model

from .providers import get_schema_llm_provider

TaskType: TypeAlias = Literal["classification", "regression"]
FineTuneMethod: TypeAlias = Literal["full", "linear_probe", "machine_learning"]
ClassificationMode: TypeAlias = Literal[
    "auto",
    "binary",
    "multiclass",
    "multilabel",
    "multioutput",
]
DataFormat: TypeAlias = Literal[
    "auto",
    "csv",
    "tsv",
    "delimited",
    "json",
    "jsonl",
    "yaml",
    "xml",
    "excel",
    "parquet",
    "sqlite",
    "fasta",
    "fastq",
    "genbank",
]

_DATA_FORMATS = {
    "auto",
    "csv",
    "tsv",
    "delimited",
    "json",
    "jsonl",
    "yaml",
    "xml",
    "excel",
    "parquet",
    "sqlite",
    "fasta",
    "fastq",
    "genbank",
}
_SEQUENCE_TRANSFORMS = {
    "strip",
    "uppercase",
    "remove_whitespace",
    "u_to_t",
    "non_iupac_to_n",
}
_ENVIRONMENT_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _single_line(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or any(character in value for character in "\r\n\x00")
    ):
        raise ConfigurationError(
            f"{name} must be non-empty single-line text.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
        )
    return value.strip()


def _optional_field(value: object, name: str) -> str | None:
    if value is None:
        return None
    return _single_line(value, name)


@dataclass(frozen=True, slots=True)
class DataMapping:
    """A declarative, non-executable mapping from source records to model rows.

    ``record_path`` and field names use dot-separated mapping keys.  Numeric list
    indexes, expressions, SQL, regular expressions, and executable code are not
    accepted.  A custom :class:`DataReader` can handle formats outside the built-in
    reader set while still returning ordinary record mappings.
    """

    format: DataFormat = "auto"
    sequence_field: str | None = "sequence"
    label_fields: tuple[str, ...] = ()
    id_field: str | None = None
    split_field: str | None = None
    record_path: str | None = None
    delimiter: str | None = None
    encoding: str = "utf-8"
    header_row: int = 0
    sheet_name: str | int | None = None
    table_name: str | None = None
    missing_values: tuple[str, ...] = ("", "NA", "N/A", "null", "None", ".")
    sequence_transforms: tuple[str, ...] = (
        "strip",
        "uppercase",
        "remove_whitespace",
        "u_to_t",
    )
    fasta_header_separator: str = "|"
    feature_fields: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        resolved_format = str(self.format).strip().lower()
        if resolved_format not in _DATA_FORMATS:
            raise ConfigurationError(
                "Unknown custom-property data format.",
                code="INVALID_CUSTOM_DATA_MAPPING",
                context={"format": self.format, "available": tuple(sorted(_DATA_FORMATS))},
            )
        object.__setattr__(self, "format", resolved_format)
        object.__setattr__(
            self,
            "sequence_field",
            _optional_field(self.sequence_field, "sequence_field"),
        )
        label_fields = tuple(_single_line(value, "label_fields") for value in self.label_fields)
        if len(set(label_fields)) != len(label_fields):
            raise ConfigurationError(
                "label_fields must be unique.", code="INVALID_CUSTOM_DATA_MAPPING"
            )
        object.__setattr__(self, "label_fields", label_fields)
        for name in ("id_field", "split_field", "record_path", "table_name"):
            object.__setattr__(self, name, _optional_field(getattr(self, name), name))
        features = tuple(_single_line(value, "feature_fields") for value in self.feature_fields)
        excluded = {*label_fields, self.id_field, self.split_field, self.sequence_field}
        if len(set(features)) != len(features) or set(features) & excluded:
            raise ConfigurationError(
                "feature_fields must be unique and separate from labels, ID, split and sequence.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )
        if self.sequence_field is None and not features:
            raise ConfigurationError(
                "A sequence_field or explicit numeric feature_fields is required.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )
        object.__setattr__(self, "feature_fields", features)
        if self.delimiter is not None and (
            not isinstance(self.delimiter, str)
            or len(self.delimiter) != 1
            or self.delimiter in "\r\n\x00"
        ):
            raise ConfigurationError(
                "delimiter must be one non-line-breaking character or None.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )
        object.__setattr__(self, "encoding", _single_line(self.encoding, "encoding"))
        if isinstance(self.header_row, bool) or not isinstance(self.header_row, int):
            raise ConfigurationError(
                "header_row must be a non-negative integer.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )
        if not 0 <= self.header_row <= 10_000:
            raise ConfigurationError(
                "header_row must be in [0, 10000].", code="INVALID_CUSTOM_DATA_MAPPING"
            )
        if self.sheet_name is not None and (
            isinstance(self.sheet_name, bool)
            or not isinstance(self.sheet_name, (str, int))
            or (isinstance(self.sheet_name, int) and self.sheet_name < 0)
            or (isinstance(self.sheet_name, str) and not self.sheet_name.strip())
        ):
            raise ConfigurationError(
                "sheet_name must be a non-empty name, a non-negative index, or None.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )
        missing_values = tuple(str(value) for value in self.missing_values)
        if len(set(missing_values)) != len(missing_values):
            raise ConfigurationError(
                "missing_values must be unique.", code="INVALID_CUSTOM_DATA_MAPPING"
            )
        object.__setattr__(self, "missing_values", missing_values)
        transforms = tuple(str(value).strip().lower() for value in self.sequence_transforms)
        unknown = sorted(set(transforms) - _SEQUENCE_TRANSFORMS)
        if unknown or len(set(transforms)) != len(transforms):
            raise ConfigurationError(
                "sequence_transforms contain unknown or duplicate operations.",
                code="INVALID_CUSTOM_DATA_MAPPING",
                context={
                    "unknown": tuple(unknown),
                    "available": tuple(sorted(_SEQUENCE_TRANSFORMS)),
                },
            )
        object.__setattr__(self, "sequence_transforms", transforms)
        if (
            not isinstance(self.fasta_header_separator, str)
            or len(self.fasta_header_separator) != 1
            or self.fasta_header_separator in "\r\n\x00"
        ):
            raise ConfigurationError(
                "fasta_header_separator must be one non-line-breaking character.",
                code="INVALID_CUSTOM_DATA_MAPPING",
            )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible mapping plan."""

        return {
            "format": self.format,
            "sequence_field": self.sequence_field,
            "feature_fields": list(self.feature_fields),
            "label_fields": list(self.label_fields),
            "id_field": self.id_field,
            "split_field": self.split_field,
            "record_path": self.record_path,
            "delimiter": self.delimiter,
            "encoding": self.encoding,
            "header_row": self.header_row,
            "sheet_name": self.sheet_name,
            "table_name": self.table_name,
            "missing_values": list(self.missing_values),
            "sequence_transforms": list(self.sequence_transforms),
            "fasta_header_separator": self.fasta_header_separator,
        }


@dataclass(frozen=True, slots=True)
class LLMSchemaConfig:
    """Configure schema analysis with provider presets or a custom compatible URL.

    New providers require an explicit model name. Omitting all provider/model
    options retains the original DeepSeek defaults. Keys are read only on use.
    """

    model: str | None = None
    base_url: str | None = None
    api_key: str | None = field(default=None, repr=False)
    api_key_env: str | None = None
    timeout_seconds: float = 60.0
    max_output_tokens: int = 2_000
    sample_rows: int = 8
    max_sample_characters: int = 12_000
    include_sample_values: bool = True
    provider: str = "deepseek"
    temperature: float | None = None
    token_limit_parameter: Literal["auto", "max_tokens", "max_completion_tokens"] = "auto"
    json_mode: bool = True
    extra_body: Mapping[str, object] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        preset = get_schema_llm_provider(self.provider)
        object.__setattr__(self, "provider", preset.name)
        model = self.model
        if model is None and preset.name == "deepseek":
            model = "deepseek-v4-flash"
        if model is None:
            raise ConfigurationError(
                "Specify a model name available in your provider account.",
                code="INVALID_LLM_SCHEMA_CONFIG",
                hint="Pass model='your-model-id'; provider presets do not select a paid model.",
            )
        object.__setattr__(self, "model", _single_line(model, "model"))
        requested_url = self.base_url if self.base_url is not None else preset.base_url
        if requested_url is None:
            raise ConfigurationError(
                "The custom provider requires an explicit base_url.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        base_url = _single_line(requested_url, "base_url").rstrip("/")
        parsed = urlparse(base_url)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ConfigurationError(
                "base_url must be an absolute HTTP(S) URL.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        if parsed.username is not None or parsed.password is not None:
            raise ConfigurationError(
                "base_url must not contain credentials.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        if parsed.query or parsed.fragment:
            raise ConfigurationError(
                "base_url must not contain a query string or fragment.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        object.__setattr__(self, "base_url", base_url)
        if self.api_key_env is None:
            object.__setattr__(self, "api_key_env", preset.api_key_env)
        if self.api_key is not None:
            try:
                api_key = _single_line(self.api_key, "api_key")
            except ConfigurationError as exc:
                raise ConfigurationError(
                    "api_key must be non-empty single-line text or None.",
                    code="INVALID_LLM_SCHEMA_CONFIG",
                ) from exc
            object.__setattr__(self, "api_key", api_key)
        if not isinstance(self.api_key_env, str) or not _ENVIRONMENT_NAME.fullmatch(
            self.api_key_env
        ):
            raise ConfigurationError(
                "api_key_env must be a valid environment-variable name.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        if (
            isinstance(self.timeout_seconds, bool)
            or not isinstance(self.timeout_seconds, (int, float))
            or not math.isfinite(self.timeout_seconds)
            or not 0 < self.timeout_seconds <= 600
        ):
            raise ConfigurationError(
                "timeout_seconds must be in (0, 600].", code="INVALID_LLM_SCHEMA_CONFIG"
            )
        for name, integer_value, lower, upper in (
            ("max_output_tokens", self.max_output_tokens, 128, 100_000),
            ("sample_rows", self.sample_rows, 1, 100),
            ("max_sample_characters", self.max_sample_characters, 256, 100_000),
        ):
            if (
                isinstance(integer_value, bool)
                or not isinstance(integer_value, int)
                or not lower <= integer_value <= upper
            ):
                raise ConfigurationError(
                    f"{name} must be in [{lower}, {upper}].",
                    code="INVALID_LLM_SCHEMA_CONFIG",
                )
        if not isinstance(self.include_sample_values, bool):
            raise ConfigurationError(
                "include_sample_values must be boolean.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        if self.temperature is not None and (
            isinstance(self.temperature, bool)
            or not isinstance(self.temperature, (int, float))
            or not math.isfinite(self.temperature)
            or not 0 <= self.temperature <= 2
        ):
            raise ConfigurationError(
                "temperature must be None or a finite number in [0, 2].",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        if not isinstance(self.token_limit_parameter, str) or self.token_limit_parameter not in {
            "auto",
            "max_tokens",
            "max_completion_tokens",
        }:
            raise ConfigurationError(
                "Unknown token_limit_parameter.", code="INVALID_LLM_SCHEMA_CONFIG"
            )
        if not isinstance(self.json_mode, bool):
            raise ConfigurationError("json_mode must be boolean.", code="INVALID_LLM_SCHEMA_CONFIG")
        protected = {
            "model",
            "messages",
            "stream",
            "response_format",
            "temperature",
            "max_tokens",
            "max_completion_tokens",
            "max_output_tokens",
            "api_key",
            "api_key_env",
            "authorization",
            "headers",
            "base_url",
        }
        if not isinstance(self.extra_body, Mapping) or any(
            not isinstance(key, str) or key.lower() in protected for key in self.extra_body
        ):
            raise ConfigurationError(
                "extra_body must be a JSON mapping without authentication or core request fields.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            )
        try:
            object.__setattr__(self, "extra_body", FrozenDict(self.extra_body))
        except ConfigurationError:
            raise ConfigurationError(
                "extra_body must contain only finite, non-recursive JSON values.",
                code="INVALID_LLM_SCHEMA_CONFIG",
            ) from None

    def resolved_api_key(self) -> str:
        """Return the explicit key or read it from the configured environment variable."""

        assert self.api_key_env is not None  # Resolved and validated during construction.
        value = self.api_key or os.environ.get(self.api_key_env)
        if (
            value is None
            or not value.strip()
            or any(character in value for character in "\r\n\x00")
        ):
            raise ConfigurationError(
                "No valid single-line API key is available for LLM schema analysis.",
                code="MISSING_LLM_API_KEY",
                hint=f"Set the {self.api_key_env} environment variable.",
            )
        return value.strip()

    @property
    def endpoint(self) -> str:
        """Return the normalized Chat Completions endpoint."""

        assert self.base_url is not None  # Resolved and validated during construction.
        if self.base_url.endswith("/chat/completions"):
            return self.base_url
        return f"{self.base_url}/chat/completions"

    @property
    def resolved_token_limit_parameter(self) -> str:
        """Use OpenAI's current token field while retaining compatible-provider support."""

        if self.token_limit_parameter != "auto":
            return self.token_limit_parameter
        if self.provider == "openai" or urlparse(self.endpoint).hostname == "api.openai.com":
            return "max_completion_tokens"
        return "max_tokens"


@dataclass(frozen=True, slots=True)
class CustomTrainingModel:
    """Fine-tuning capabilities for one registered DNA foundation model."""

    name: str
    display_name: str
    linear_probe: bool
    full_finetune: bool
    notes: str

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "linear_probe": self.linear_probe,
            "full_finetune": self.full_finetune,
            "notes": self.notes,
        }


def available_custom_training_models(*, verbose: bool = True) -> tuple[CustomTrainingModel, ...]:
    """Return training support for every registered DNA representation model."""

    from dnakit.representations import available_embedding_models

    models: list[CustomTrainingModel] = []
    for name in available_embedding_models():
        spec = get_embedding_model(name)
        full = spec.loader == "transformers"
        notes = (
            "Full fine-tuning and frozen linear probing are supported."
            if full
            else (
                "Frozen linear probing is supported; this backend does not expose a "
                "portable PyTorch full-fine-tuning path."
            )
        )
        models.append(CustomTrainingModel(name, spec.display_name, True, full, notes))
    return tuple(models)


def get_custom_training_model(name: str, *, verbose: bool = True) -> CustomTrainingModel:
    """Resolve fine-tuning capabilities for a registered DNA model."""

    canonical = get_embedding_model(name).name
    return next(item for item in available_custom_training_models() if item.name == canonical)


@dataclass(frozen=True, slots=True)
class CustomPropertyTrainingConfig:
    """Configure deep-learning or machine-learning property training.

    With ``finetune_method='machine_learning'``, ``model`` selects an estimator
    family, ``model_params`` configures it, and ``feature_source`` selects sequence
    features, frozen embeddings, or explicitly mapped numeric columns.

    Set ``cv_folds=10`` for complete cross-validation on a training-only source,
    followed by a final fit on all its rows. This ML-only mode ignores the three
    split ratios, rejects mapped split columns, and never reads an external test
    set. Fold membership and held-out predictions are saved locally.
    """

    task_type: TaskType
    model: str = "dnabert2"
    finetune_method: FineTuneMethod = "full"
    classification_mode: ClassificationMode = "auto"
    mapping: DataMapping | None = None
    llm: LLMSchemaConfig | None = None
    output_dir: str | os.PathLike[str] = "property-model"
    train_ratio: float = 0.8
    validation_ratio: float = 0.1
    test_ratio: float = 0.1
    seed: int = 42
    epochs: int = 10
    batch_size: int = 32
    workers: int = 6
    learning_rate: float | None = 1e-4
    beta1: float = 0.9
    beta2: float = 0.95
    weight_decay: float = 0.1
    early_stopping_patience: int = 3
    threshold: float = 0.5
    pooling: Literal["mean", "cls", "max", "last"] = "mean"
    max_length: int | None = None
    max_records: int = 1_000_000
    max_file_bytes: int = 1_000_000_000
    device: str = "auto"
    dtype: Literal["auto", "float32", "float16", "bfloat16"] = "auto"
    checkpoint_dir: str | os.PathLike[str] | None = None
    checkpoint_path: str | os.PathLike[str] | None = None
    model_source_path: str | os.PathLike[str] | None = None
    hf_token: str | None = field(default=None, repr=False)
    allow_remote_code: bool = False
    show_progress: bool = True
    overwrite: bool = False
    gradient_accumulation_steps: int = 1
    max_grad_norm: float = 1.0
    missing_label_policy: Literal["error", "drop"] = "error"
    stratify: bool = True
    feature_source: Literal["sequence", "embedding", "numeric"] = "sequence"
    embedding_model: str = "dnabert2"
    model_params: Mapping[str, object] = field(default_factory=dict)
    kmer_k: int = 3
    standardize_features: bool = True
    cv_folds: int | None = None

    def __post_init__(self) -> None:
        if self.task_type not in {"classification", "regression"}:
            raise ConfigurationError(
                "task_type must be 'classification' or 'regression'.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if self.finetune_method not in {"full", "linear_probe", "machine_learning"}:
            raise ConfigurationError(
                "finetune_method must be 'full', 'linear_probe' or 'machine_learning'.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if self.cv_folds is not None and (
            isinstance(self.cv_folds, bool)
            or not isinstance(self.cv_folds, int)
            or not 2 <= self.cv_folds <= 100
            or self.finetune_method != "machine_learning"
        ):
            raise ConfigurationError(
                "cv_folds must be an integer in [2, 100] for machine_learning training.",
                code="INVALID_CUSTOM_CV",
            )
        if self.finetune_method == "machine_learning":
            from .ml_models import validate_ml_config

            validate_ml_config(self)
        else:
            capability = get_custom_training_model(self.model)
            object.__setattr__(self, "model", capability.name)
        if self.finetune_method == "full" and not capability.full_finetune:
            raise ConfigurationError(
                "The selected model does not provide full fine-tuning in DNAKit.",
                code="CUSTOM_FULL_FINETUNE_UNSUPPORTED",
                context={"model": capability.name},
                hint="Choose finetune_method='linear_probe' or a model with full_finetune=True.",
            )
        if self.classification_mode not in {
            "auto",
            "binary",
            "multiclass",
            "multilabel",
            "multioutput",
        }:
            raise ConfigurationError(
                "Unknown classification_mode.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            )
        if self.task_type == "regression" and self.classification_mode != "auto":
            raise ConfigurationError(
                "classification_mode applies only to classification.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if self.mapping is not None and not isinstance(self.mapping, DataMapping):
            raise ConfigurationError(
                "mapping must be DataMapping or None.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if (
            self.mapping is not None
            and self.mapping.sequence_field is None
            and not (
                self.finetune_method == "machine_learning" and self.feature_source == "numeric"
            )
        ):
            raise ConfigurationError(
                "This feature source requires a mapped sequence field.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if self.llm is not None and not isinstance(self.llm, LLMSchemaConfig):
            raise ConfigurationError(
                "llm must be LLMSchemaConfig or None.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        try:
            output_dir = os.fspath(self.output_dir)
        except TypeError as exc:
            raise ConfigurationError(
                "output_dir must be path-like.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            ) from exc
        if not output_dir.strip():
            raise ConfigurationError(
                "output_dir must not be empty.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            )
        object.__setattr__(self, "output_dir", output_dir)
        ratios = (self.train_ratio, self.validation_ratio, self.test_ratio)
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or (value <= 0 if self.cv_folds is None else value < 0)
            for value in ratios
        ) or not math.isclose(sum(ratios), 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ConfigurationError(
                "Split ratios must sum to 1 and be positive (nonnegative in CV mode).",
                code="INVALID_CUSTOM_SPLIT",
            )
        for name, bounded_integer, lower, upper in (
            ("seed", self.seed, 0, 2**32 - 1),
            ("epochs", self.epochs, 1, 1_000_000),
            ("batch_size", self.batch_size, 1, 65_536),
            ("workers", self.workers, 0, 256),
            ("early_stopping_patience", self.early_stopping_patience, 0, 1_000_000),
            ("max_records", self.max_records, 3, 10_000_000),
            ("max_file_bytes", self.max_file_bytes, 1, 100_000_000_000),
            (
                "gradient_accumulation_steps",
                self.gradient_accumulation_steps,
                1,
                65_536,
            ),
        ):
            if (
                isinstance(bounded_integer, bool)
                or not isinstance(bounded_integer, int)
                or not lower <= bounded_integer <= upper
            ):
                raise ConfigurationError(
                    f"{name} must be in [{lower}, {upper}].",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                )
        numeric_options: tuple[tuple[str, int | float | None, bool], ...] = (
            ("learning_rate", self.learning_rate, False),
            ("weight_decay", self.weight_decay, True),
            ("max_grad_norm", self.max_grad_norm, False),
        )
        for name, numeric_value, allow_zero in numeric_options:
            if numeric_value is None:
                if name == "learning_rate":
                    continue
                raise ConfigurationError(
                    f"{name} must be a number.",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                )
            if (
                isinstance(numeric_value, bool)
                or not isinstance(numeric_value, (int, float))
                or not math.isfinite(numeric_value)
                or numeric_value < 0
                or (not allow_zero and numeric_value == 0)
            ):
                constraint = "non-negative" if allow_zero else "positive"
                raise ConfigurationError(
                    f"{name} must be a finite {constraint} number.",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                )
        for name, beta in (("beta1", self.beta1), ("beta2", self.beta2)):
            if (
                isinstance(beta, bool)
                or not isinstance(beta, (int, float))
                or not math.isfinite(beta)
                or not 0 <= beta < 1
            ):
                raise ConfigurationError(
                    f"{name} must be a finite number in [0, 1).",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                )
        if (
            isinstance(self.threshold, bool)
            or not isinstance(self.threshold, (int, float))
            or not math.isfinite(self.threshold)
            or not 0 < self.threshold < 1
        ):
            raise ConfigurationError(
                "threshold must be in (0, 1).", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            )
        if self.pooling not in {"mean", "cls", "max", "last"}:
            raise ConfigurationError(
                "Unknown pooling method.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            )
        if self.max_length is not None and (
            isinstance(self.max_length, bool)
            or not isinstance(self.max_length, int)
            or not 1 <= self.max_length <= 10_000_000
        ):
            raise ConfigurationError(
                "max_length must be None or an integer in [1, 10000000].",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if not isinstance(self.device, str) or not self.device.strip():
            raise ConfigurationError(
                "device must be non-empty text.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
            )
        if self.dtype not in {"auto", "float32", "float16", "bfloat16"}:
            raise ConfigurationError("Unknown dtype.", code="INVALID_CUSTOM_PROPERTY_CONFIG")
        for path_name in ("checkpoint_dir", "checkpoint_path", "model_source_path"):
            value = getattr(self, path_name)
            if value is None:
                continue
            try:
                path = os.fspath(value)
            except TypeError as exc:
                raise ConfigurationError(
                    f"{path_name} must be path-like or None.",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                ) from exc
            if not path.strip():
                raise ConfigurationError(
                    f"{path_name} must not be empty.",
                    code="INVALID_CUSTOM_PROPERTY_CONFIG",
                )
            object.__setattr__(self, path_name, path)
        if self.checkpoint_dir is not None and self.checkpoint_path is not None:
            raise ConfigurationError(
                "checkpoint_dir and checkpoint_path are mutually exclusive.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        for name in ("allow_remote_code", "show_progress", "overwrite", "stratify"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(
                    f"{name} must be boolean.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
                )
        if self.missing_label_policy not in {"error", "drop"}:
            raise ConfigurationError(
                "missing_label_policy must be 'error' or 'drop'.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )

    @property
    def resolved_learning_rate(self) -> float:
        """Return the configured rate, defaulting to 1e-4 for either training method."""

        if self.learning_rate is not None:
            return float(self.learning_rate)
        return 1e-4


@dataclass(frozen=True, slots=True)
class TargetSchema:
    """Serializable encoding contract learned from the training split."""

    task_type: TaskType
    mode: str
    label_names: tuple[str, ...]
    class_names: tuple[tuple[str, ...], ...] = ()
    regression_means: tuple[float, ...] = ()
    regression_scales: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if self.task_type not in {"classification", "regression"}:
            raise ConfigurationError(
                "The target task type is invalid.", code="INVALID_CUSTOM_TARGET_SCHEMA"
            )
        label_names = tuple(_single_line(value, "label_names") for value in self.label_names)
        if not label_names or len(set(label_names)) != len(label_names):
            raise ConfigurationError(
                "Target label names must be non-empty and unique.",
                code="INVALID_CUSTOM_TARGET_SCHEMA",
            )
        object.__setattr__(self, "label_names", label_names)
        class_names = tuple(
            tuple(_single_line(value, "class_names") for value in classes)
            for classes in self.class_names
        )
        if any(len(set(classes)) != len(classes) for classes in class_names):
            raise ConfigurationError(
                "Class names must be unique within every target.",
                code="INVALID_CUSTOM_TARGET_SCHEMA",
            )
        object.__setattr__(self, "class_names", class_names)
        means = tuple(float(value) for value in self.regression_means)
        scales = tuple(float(value) for value in self.regression_scales)
        if not all(math.isfinite(value) for value in (*means, *scales)):
            raise ConfigurationError(
                "Regression target normalization must be finite.",
                code="INVALID_CUSTOM_TARGET_SCHEMA",
            )
        object.__setattr__(self, "regression_means", means)
        object.__setattr__(self, "regression_scales", scales)

        if self.task_type == "regression":
            valid = (
                self.mode == "regression"
                and not class_names
                and len(means) == len(label_names)
                and len(scales) == len(label_names)
                and all(value > 0 for value in scales)
            )
        elif self.mode == "binary":
            valid = (
                len(label_names) == 1
                and len(class_names) == 1
                and len(class_names[0]) == 2
                and not means
                and not scales
            )
        elif self.mode == "multiclass":
            valid = (
                len(class_names) == 1
                and len(class_names[0]) >= 2
                and (len(label_names) == 1 or class_names[0] == label_names)
                and not means
                and not scales
            )
        elif self.mode == "multilabel":
            valid = (
                len(label_names) >= 2
                and len(class_names) == len(label_names)
                and all(classes == ("0", "1") for classes in class_names)
                and not means
                and not scales
            )
        elif self.mode == "multioutput":
            valid = (
                len(label_names) >= 2
                and len(class_names) == len(label_names)
                and all(len(classes) >= 2 for classes in class_names)
                and not means
                and not scales
            )
        else:
            valid = False
        if not valid:
            raise ConfigurationError(
                "The target schema fields are inconsistent.",
                code="INVALID_CUSTOM_TARGET_SCHEMA",
            )

    @property
    def output_dimension(self) -> int:
        if self.task_type == "regression" or self.mode == "multilabel":
            return len(self.label_names)
        return sum(len(values) for values in self.class_names)

    def to_dict(self) -> dict[str, object]:
        return {
            "task_type": self.task_type,
            "mode": self.mode,
            "label_names": list(self.label_names),
            "class_names": [list(values) for values in self.class_names],
            "regression_means": list(self.regression_means),
            "regression_scales": list(self.regression_scales),
            "output_dimension": self.output_dimension,
        }


@dataclass(frozen=True, slots=True)
class DataInspection:
    """Auditable source inspection and the selected mapping plan."""

    source_path: str
    detected_format: str
    source_size_bytes: int
    columns: tuple[str, ...]
    estimated_records: int | None
    mapping: DataMapping
    mapping_origin: Literal["explicit", "local", "llm"]
    recommendations: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    llm_model: str | None = None
    shared_sample_records: int = 0
    shared_sample_characters: int = 0

    def to_dict(self) -> dict[str, object]:
        return {
            "source_path": self.source_path,
            "detected_format": self.detected_format,
            "source_size_bytes": self.source_size_bytes,
            "columns": list(self.columns),
            "estimated_records": self.estimated_records,
            "mapping": self.mapping.to_dict(),
            "mapping_origin": self.mapping_origin,
            "recommendations": list(self.recommendations),
            "warnings": list(self.warnings),
            "llm_model": self.llm_model,
            "shared_sample_records": self.shared_sample_records,
            "shared_sample_characters": self.shared_sample_characters,
        }


@dataclass(frozen=True, slots=True)
class CustomPropertyTrainingResult:
    """Training metrics and paths for one user-defined property model."""

    model_name: str
    task_type: str
    finetune_method: str
    output_dir: str
    target_schema: TargetSchema
    split_counts: FrozenDict
    metrics: FrozenDict
    history: tuple[FrozenDict, ...]
    artifacts: FrozenDict
    data_inspection: DataInspection
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "model_name": self.model_name,
            "task_type": self.task_type,
            "finetune_method": self.finetune_method,
            "output_dir": self.output_dir,
            "target_schema": self.target_schema.to_dict(),
            "split_counts": to_json_compatible(self.split_counts),
            "metrics": to_json_compatible(self.metrics),
            "history": to_json_compatible(self.history),
            "artifacts": to_json_compatible(self.artifacts),
            "data_inspection": self.data_inspection.to_dict(),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class CustomPropertyPredictionResult:
    """Predictions produced by a saved user-defined property model."""

    model_dir: str
    record_ids: tuple[str, ...]
    label_names: tuple[str, ...]
    predictions: tuple[FrozenDict, ...]
    probabilities: tuple[FrozenDict, ...] = ()
    output_file: FrozenDict | None = None

    def to_dict(self) -> dict[str, object]:
        if self.output_file is not None:
            return {
                "model_dir": self.model_dir,
                **to_json_compatible(self.output_file),
            }
        return {
            "model_dir": self.model_dir,
            "record_ids": list(self.record_ids),
            "label_names": list(self.label_names),
            "predictions": to_json_compatible(self.predictions),
            "probabilities": to_json_compatible(self.probabilities),
        }


__all__ = [
    "ClassificationMode",
    "CustomPropertyPredictionResult",
    "CustomPropertyTrainingConfig",
    "CustomPropertyTrainingResult",
    "CustomTrainingModel",
    "DataFormat",
    "DataInspection",
    "DataMapping",
    "FineTuneMethod",
    "LLMSchemaConfig",
    "TargetSchema",
    "TaskType",
    "available_custom_training_models",
    "get_custom_training_model",
]

_register_public_verbose(globals())
