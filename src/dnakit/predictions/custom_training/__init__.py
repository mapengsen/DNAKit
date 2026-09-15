"""User-defined DNA sequence property training and inference."""

from dnakit._verbosity import register_public as _register_public_verbose

from .data import (
    DataReader,
    PropertyData,
    detect_data_format,
    inspect_property_data,
    load_property_data,
)
from .encoder import TrainableSequenceEncoder
from .llm import (
    SCHEMA_ANALYSIS_SYSTEM_PROMPT,
    OpenAICompatibleSchemaClient,
    SchemaLLMClient,
)
from .ml_models import (
    MachineLearningModel,
    available_machine_learning_models,
    get_machine_learning_model,
)
from .providers import (
    SchemaLLMProvider,
    available_schema_llm_providers,
    get_schema_llm_provider,
)
from .schema import (
    ClassificationMode,
    CustomPropertyPredictionResult,
    CustomPropertyTrainingConfig,
    CustomPropertyTrainingResult,
    CustomTrainingModel,
    DataFormat,
    DataInspection,
    DataMapping,
    FineTuneMethod,
    LLMSchemaConfig,
    TargetSchema,
    TaskType,
    available_custom_training_models,
    get_custom_training_model,
)
from .splitting import SplitAssignment, split_property_data
from .targets import compute_metrics
from .training import predict_custom_properties, train_custom_property_model

__all__ = [
    "SCHEMA_ANALYSIS_SYSTEM_PROMPT",
    "ClassificationMode",
    "CustomPropertyPredictionResult",
    "CustomPropertyTrainingConfig",
    "CustomPropertyTrainingResult",
    "CustomTrainingModel",
    "DataFormat",
    "DataInspection",
    "DataMapping",
    "DataReader",
    "FineTuneMethod",
    "LLMSchemaConfig",
    "MachineLearningModel",
    "OpenAICompatibleSchemaClient",
    "PropertyData",
    "SchemaLLMClient",
    "SchemaLLMProvider",
    "SplitAssignment",
    "TargetSchema",
    "TaskType",
    "TrainableSequenceEncoder",
    "available_custom_training_models",
    "available_machine_learning_models",
    "available_schema_llm_providers",
    "compute_metrics",
    "detect_data_format",
    "get_custom_training_model",
    "get_machine_learning_model",
    "get_schema_llm_provider",
    "inspect_property_data",
    "load_property_data",
    "predict_custom_properties",
    "split_property_data",
    "train_custom_property_model",
]

_register_public_verbose(globals())
