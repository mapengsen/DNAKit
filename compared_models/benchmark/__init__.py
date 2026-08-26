from .data import (
    SequenceBenchmarkCollator,
    SequenceBenchmarkDataset,
    build_sequence_dataloader,
    load_processed_task_dataframe,
    load_processed_task_meta,
    resolve_processed_task_names,
)
from .modeling import (
    MODEL_SPECS,
    FrozenSequenceLinearProbe,
    available_compared_models,
    build_linear_probe_model,
    get_model_spec,
    normalize_finetune_method,
    validate_finetune_method,
)
from .nt_preprocess import (
    available_nt_tasks,
    prepare_nt_task,
)

try:
    from .eqtl_preprocess import (
        available_eqtl_tasks,
        prepare_eqtl_task,
    )
except ModuleNotFoundError:  # pragma: no cover - optional dependency by env
    available_eqtl_tasks = None
    prepare_eqtl_task = None

try:
    from .genomic_benchmarks_preprocess import (
        available_genomic_benchmark_tasks,
        prepare_genomic_benchmark_task,
    )
except ModuleNotFoundError:  # pragma: no cover - optional dependency by env
    available_genomic_benchmark_tasks = None
    prepare_genomic_benchmark_task = None

__all__ = [
    "MODEL_SPECS",
    "FrozenSequenceLinearProbe",
    "SequenceBenchmarkCollator",
    "SequenceBenchmarkDataset",
    "available_compared_models",
    "available_eqtl_tasks",
    "available_genomic_benchmark_tasks",
    "available_nt_tasks",
    "build_linear_probe_model",
    "build_sequence_dataloader",
    "get_model_spec",
    "load_processed_task_dataframe",
    "load_processed_task_meta",
    "normalize_finetune_method",
    "prepare_eqtl_task",
    "prepare_genomic_benchmark_task",
    "prepare_nt_task",
    "resolve_processed_task_names",
    "validate_finetune_method",
]
