"""Lazy registry and parameter validation for supervised machine-learning models."""

from __future__ import annotations

import importlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

if TYPE_CHECKING:
    from .schema import CustomPropertyTrainingConfig


@dataclass(frozen=True, slots=True)
class MachineLearningModel:
    """An algorithm family and its supported supervised tasks."""

    name: str
    module: str
    classifier: str | None
    regressor: str | None
    nonnegative_features: bool = False

    @property
    def tasks(self) -> tuple[str, ...]:
        return tuple(
            task
            for task, implementation in (
                ("classification", self.classifier),
                ("regression", self.regressor),
            )
            if implementation is not None
        )

    @property
    def dependency(self) -> str:
        return "scikit-learn" if self.module.startswith("sklearn.") else self.module

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "tasks": list(self.tasks),
            "dependency": self.dependency,
            "nonnegative_features": self.nonnegative_features,
        }


_MODELS = (
    MachineLearningModel("logistic_regression", "sklearn.linear_model", "LogisticRegression", None),
    MachineLearningModel("linear_regression", "sklearn.linear_model", None, "LinearRegression"),
    MachineLearningModel("ridge", "sklearn.linear_model", "RidgeClassifier", "Ridge"),
    MachineLearningModel("lasso", "sklearn.linear_model", None, "Lasso"),
    MachineLearningModel("elastic_net", "sklearn.linear_model", None, "ElasticNet"),
    MachineLearningModel("bayesian_ridge", "sklearn.linear_model", None, "BayesianRidge"),
    MachineLearningModel("ard", "sklearn.linear_model", None, "ARDRegression"),
    MachineLearningModel("huber", "sklearn.linear_model", None, "HuberRegressor"),
    MachineLearningModel("ransac", "sklearn.linear_model", None, "RANSACRegressor"),
    MachineLearningModel("theil_sen", "sklearn.linear_model", None, "TheilSenRegressor"),
    MachineLearningModel("sgd", "sklearn.linear_model", "SGDClassifier", "SGDRegressor"),
    MachineLearningModel("perceptron", "sklearn.linear_model", "Perceptron", None),
    MachineLearningModel("svm", "sklearn.svm", "SVC", "SVR"),
    MachineLearningModel("linear_svm", "sklearn.svm", "LinearSVC", "LinearSVR"),
    MachineLearningModel("nu_svm", "sklearn.svm", "NuSVC", "NuSVR"),
    MachineLearningModel("kernel_ridge", "sklearn.kernel_ridge", None, "KernelRidge"),
    MachineLearningModel("knn", "sklearn.neighbors", "KNeighborsClassifier", "KNeighborsRegressor"),
    MachineLearningModel(
        "decision_tree", "sklearn.tree", "DecisionTreeClassifier", "DecisionTreeRegressor"
    ),
    MachineLearningModel("extra_tree", "sklearn.tree", "ExtraTreeClassifier", "ExtraTreeRegressor"),
    MachineLearningModel(
        "random_forest", "sklearn.ensemble", "RandomForestClassifier", "RandomForestRegressor"
    ),
    MachineLearningModel(
        "extra_trees", "sklearn.ensemble", "ExtraTreesClassifier", "ExtraTreesRegressor"
    ),
    MachineLearningModel(
        "gradient_boosting",
        "sklearn.ensemble",
        "GradientBoostingClassifier",
        "GradientBoostingRegressor",
    ),
    MachineLearningModel(
        "hist_gradient_boosting",
        "sklearn.ensemble",
        "HistGradientBoostingClassifier",
        "HistGradientBoostingRegressor",
    ),
    MachineLearningModel("adaboost", "sklearn.ensemble", "AdaBoostClassifier", "AdaBoostRegressor"),
    MachineLearningModel("bagging", "sklearn.ensemble", "BaggingClassifier", "BaggingRegressor"),
    MachineLearningModel("gaussian_nb", "sklearn.naive_bayes", "GaussianNB", None),
    MachineLearningModel("multinomial_nb", "sklearn.naive_bayes", "MultinomialNB", None, True),
    MachineLearningModel("bernoulli_nb", "sklearn.naive_bayes", "BernoulliNB", None),
    MachineLearningModel("complement_nb", "sklearn.naive_bayes", "ComplementNB", None, True),
    MachineLearningModel(
        "lda", "sklearn.discriminant_analysis", "LinearDiscriminantAnalysis", None
    ),
    MachineLearningModel(
        "qda", "sklearn.discriminant_analysis", "QuadraticDiscriminantAnalysis", None
    ),
    MachineLearningModel(
        "gaussian_process",
        "sklearn.gaussian_process",
        "GaussianProcessClassifier",
        "GaussianProcessRegressor",
    ),
    MachineLearningModel("mlp", "sklearn.neural_network", "MLPClassifier", "MLPRegressor"),
    MachineLearningModel("dummy", "sklearn.dummy", "DummyClassifier", "DummyRegressor"),
    MachineLearningModel("xgboost", "xgboost", "XGBClassifier", "XGBRegressor"),
    MachineLearningModel("lightgbm", "lightgbm", "LGBMClassifier", "LGBMRegressor"),
    MachineLearningModel("catboost", "catboost", "CatBoostClassifier", "CatBoostRegressor"),
)


def available_machine_learning_models() -> tuple[MachineLearningModel, ...]:
    """List algorithm capabilities without importing optional dependencies."""
    return _MODELS


def get_machine_learning_model(name: str) -> MachineLearningModel:
    """Resolve a supported algorithm family."""
    canonical = str(name).strip().lower().replace("-", "_")
    for model in _MODELS:
        if model.name == canonical:
            return model
    raise ConfigurationError(
        "Unknown machine-learning model.",
        code="UNKNOWN_MACHINE_LEARNING_MODEL",
        context={"model": name, "available": tuple(model.name for model in _MODELS)},
    )


def validate_ml_config(config: CustomPropertyTrainingConfig) -> None:
    """Validate the ML branch without importing estimator libraries."""
    from .schema import DataMapping

    spec = get_machine_learning_model(config.model)
    object.__setattr__(config, "model", spec.name)
    if config.task_type not in spec.tasks:
        raise ConfigurationError(
            "The selected machine-learning model does not support this task.",
            code="UNSUPPORTED_MACHINE_LEARNING_TASK",
            context={"model": spec.name, "tasks": spec.tasks},
        )
    if config.feature_source not in {"sequence", "embedding", "numeric"}:
        raise ConfigurationError("Unknown feature_source.", code="INVALID_CUSTOM_PROPERTY_CONFIG")
    if config.feature_source == "embedding":
        from dnakit.representations import get_embedding_model

        object.__setattr__(
            config, "embedding_model", get_embedding_model(config.embedding_model).name
        )
    if config.feature_source == "numeric" and (
        not isinstance(config.mapping, DataMapping) or not config.mapping.feature_fields
    ):
        raise ConfigurationError(
            "Numeric input requires explicit DataMapping.feature_fields and label_fields.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
        )
    if (
        isinstance(config.kmer_k, bool)
        or not isinstance(config.kmer_k, int)
        or not 1 <= config.kmer_k <= 6
    ):
        raise ConfigurationError("kmer_k must be in [1, 6].", code="INVALID_CUSTOM_PROPERTY_CONFIG")
    if not isinstance(config.standardize_features, bool):
        raise ConfigurationError(
            "standardize_features must be boolean.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
        )
    if not isinstance(config.model_params, Mapping) or not all(
        isinstance(key, str) for key in config.model_params
    ):
        raise ConfigurationError(
            "model_params must be a string-keyed mapping.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
        )
    try:
        # Deliberately accept only portable parameters, not executable estimator objects.
        raw = (
            to_json_compatible(config.model_params)
            if isinstance(config.model_params, FrozenDict)
            else dict(config.model_params)
        )
        payload = json.loads(json.dumps(raw, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "model_params must contain finite JSON values.", code="INVALID_CUSTOM_PROPERTY_CONFIG"
        ) from exc
    object.__setattr__(config, "model_params", FrozenDict(payload))


def require_ml_module(name: str) -> Any:
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        extra = name if name in {"xgboost", "lightgbm", "catboost"} else "ml"
        raise BackendUnavailableError(
            f"Machine-learning dependency {name!r} is not installed.",
            code="MISSING_MACHINE_LEARNING_DEPENDENCY",
            hint=f'Install with: python -m pip install "dnakit[{extra}]"',
        ) from exc


def create_estimator(config: CustomPropertyTrainingConfig) -> Any:
    spec = get_machine_learning_model(config.model)
    implementation = spec.classifier if config.task_type == "classification" else spec.regressor
    estimator_type = getattr(require_ml_module(spec.module), str(implementation))
    supported = estimator_type().get_params(deep=False)
    defaults: dict[str, object] = {}
    if "random_state" in supported:
        defaults["random_state"] = config.seed
    if "n_jobs" in supported:
        defaults["n_jobs"] = max(1, config.workers)
    if config.model in {"logistic_regression", "mlp"}:
        defaults["max_iter"] = 1000
    if config.model == "sgd" and config.task_type == "classification":
        defaults["loss"] = "log_loss"
    if config.model == "lightgbm":
        defaults["verbosity"] = -1
    if config.model == "catboost":
        defaults.update(random_seed=config.seed, thread_count=max(1, config.workers), verbose=False)
    params = dict(cast(Mapping[str, object], to_json_compatible(config.model_params)))
    # Targets use integer class indices; JSON object keys round-trip as text.
    for weight_parameter in ("class_weight", "class_weights"):
        weights = params.get(weight_parameter)
        if isinstance(weights, Mapping):
            try:
                if any(str(int(key)) != str(key) or int(key) < 0 for key in weights):
                    raise ValueError("invalid class index")
                params[weight_parameter] = {int(key): value for key, value in weights.items()}
            except (TypeError, ValueError) as exc:
                raise ConfigurationError(
                    "Class-weight dictionaries must use nonnegative encoded class indices.",
                    code="INVALID_MACHINE_LEARNING_PARAMS",
                ) from exc
    if spec.module.startswith("sklearn."):
        unknown = set(params) - set(supported)
        if unknown:
            raise ConfigurationError(
                "Unknown model_params for the selected estimator.",
                code="INVALID_MACHINE_LEARNING_PARAMS",
                context={"parameters": tuple(sorted(unknown)), "model": config.model},
            )
    defaults.update(params)
    if config.model == "catboost":
        # Keep CatBoost's auxiliary logs inside the explicitly managed artifact flow.
        defaults["allow_writing_files"] = False
        if "random_state" in params:
            defaults.pop("random_seed", None)
    try:
        return estimator_type(**defaults)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Invalid estimator parameters.", code="INVALID_MACHINE_LEARNING_PARAMS"
        ) from exc


__all__ = [
    "MachineLearningModel",
    "available_machine_learning_models",
    "get_machine_learning_model",
]
