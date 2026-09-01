"""Build JSON schemas and adapt JSON values to DNAKit public API inputs."""

from __future__ import annotations

import inspect
import os
import types
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence, Set
from dataclasses import dataclass
from enum import Enum
from pathlib import Path, PurePath
from typing import (
    Annotated,
    Any,
    Literal,
    TypeVar,
    Union,
    cast,
    get_args,
    get_origin,
    get_type_hints,
)

from dnakit.core import DNA, DNARecord, DNASequence, DNASet
from dnakit.exceptions import ConfigurationError

JSONSchema = dict[str, object]
_NONE_TYPE = type(None)
_MAX_SCHEMA_DEPTH = 8


@dataclass(frozen=True, slots=True)
class _AgentDictResult:
    values: Mapping[str, object]

    def to_dict(self) -> dict[str, object]:
        return dict(self.values)


@dataclass(frozen=True, slots=True)
class ParameterPlan:
    """One public callable parameter and its Agent-facing conversion contract."""

    name: str
    annotation: object
    kind: str
    required: bool
    supported: bool
    schema: JSONSchema
    default: object
    has_default: bool
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class CallPlan:
    """Resolved schema and argument-conversion plan for one public callable."""

    signature: str
    parameters: tuple[ParameterPlan, ...]
    input_schema: JSONSchema
    output_schema: JSONSchema
    accepts_keyword_arguments: bool
    compatible: bool
    incompatibilities: tuple[str, ...]


def _error(message: str, *, path: str, value: object | None = None) -> ConfigurationError:
    context: dict[str, object] = {"path": path}
    if value is not None:
        context["received_type"] = type(value).__name__
    return ConfigurationError(
        message,
        code="AGENT_ARGUMENT_TYPE",
        context=context,
        hint="Use describe_dnakit_tool to inspect the accepted JSON input schema.",
    )


def _resolve_hints(target: Callable[..., object] | type[object]) -> dict[str, object]:
    hint_target: object = target.__init__ if inspect.isclass(target) else target
    try:
        return dict(get_type_hints(hint_target, include_extras=True))
    except (NameError, TypeError):
        annotations = getattr(hint_target, "__annotations__", {})
        globals_mapping = getattr(hint_target, "__globals__", {})
        resolved: dict[str, object] = {}
        for name, annotation in annotations.items():
            if not isinstance(annotation, str):
                resolved[name] = annotation
                continue
            try:
                resolved[name] = eval(annotation, globals_mapping, {})
            except (NameError, TypeError, SyntaxError):
                resolved[name] = annotation
        return resolved


def _json_scalar_schema(values: tuple[object, ...]) -> JSONSchema:
    non_null = tuple(value for value in values if value is not None)
    types_seen = {
        "boolean"
        if isinstance(value, bool)
        else "integer"
        if isinstance(value, int)
        else "number"
        if isinstance(value, float)
        else "string"
        for value in non_null
    }
    schema: JSONSchema = {"enum": list(values)}
    if len(types_seen) == 1 and None not in values:
        schema["type"] = next(iter(types_seen))
    return schema


def _default_for_schema(value: object) -> object | None:
    if isinstance(value, Enum):
        return cast(object, value.value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        converted = [_default_for_schema(item) for item in value]
        if all(
            item is not None or original is None
            for item, original in zip(converted, value, strict=True)
        ):
            return converted
    if isinstance(value, Mapping) and all(isinstance(key, str) for key in value):
        converted_mapping = {str(key): _default_for_schema(item) for key, item in value.items()}
        if all(item is not None or value[key] is None for key, item in converted_mapping.items()):
            return converted_mapping
    return None


def _class_schema(
    annotation: type[object],
    *,
    seen: frozenset[type[object]],
    depth: int,
) -> tuple[bool, JSONSchema, str | None]:
    if annotation in seen or depth >= _MAX_SCHEMA_DEPTH:
        return True, {"type": "object", "x-python-type": annotation.__qualname__}, None
    try:
        signature = inspect.signature(annotation)
    except (TypeError, ValueError):
        return False, {}, f"{annotation.__qualname__} cannot be constructed from JSON"

    hints = _resolve_hints(annotation)
    properties: dict[str, object] = {}
    required: list[str] = []
    next_seen = seen | {annotation}
    for parameter in signature.parameters.values():
        if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
            continue
        if parameter.kind is inspect.Parameter.VAR_KEYWORD:
            continue
        child_annotation = hints.get(parameter.name, parameter.annotation)
        supported, child_schema, reason = _annotation_schema(
            child_annotation,
            seen=next_seen,
            depth=depth + 1,
        )
        is_required = parameter.default is inspect.Parameter.empty
        if not supported:
            if is_required:
                return (
                    False,
                    {},
                    f"{annotation.__qualname__}.{parameter.name}: {reason or 'unsupported type'}",
                )
            continue
        if parameter.default is not inspect.Parameter.empty:
            default = _default_for_schema(parameter.default)
            if default is not None or parameter.default is None:
                child_schema = dict(child_schema)
                child_schema["default"] = default
        properties[parameter.name] = child_schema
        if is_required:
            required.append(parameter.name)
    schema: JSONSchema = {
        "type": "object",
        "properties": properties,
        "additionalProperties": False,
        "x-python-type": f"{annotation.__module__}.{annotation.__qualname__}",
    }
    if required:
        schema["required"] = required
    return True, schema, None


def _annotation_schema(
    annotation: object,
    *,
    seen: frozenset[type[object]] = frozenset(),
    depth: int = 0,
) -> tuple[bool, JSONSchema, str | None]:
    if annotation is inspect.Parameter.empty or annotation is Any or annotation is object:
        return True, {}, None
    if annotation is None or annotation is _NONE_TYPE:
        return True, {"type": "null"}, None
    if isinstance(annotation, str):
        return False, {}, f"unresolved annotation {annotation!r}"
    if isinstance(annotation, TypeVar):
        if annotation.__bound__ is not None:
            return _annotation_schema(annotation.__bound__, seen=seen, depth=depth + 1)
        if annotation.__constraints__:
            schemas: list[JSONSchema] = []
            for constraint in annotation.__constraints__:
                supported, schema, reason = _annotation_schema(
                    constraint,
                    seen=seen,
                    depth=depth + 1,
                )
                if not supported:
                    return False, {}, reason
                schemas.append(schema)
            return True, {"anyOf": schemas}, None
        return True, {}, None

    origin = get_origin(annotation)
    arguments = get_args(annotation)
    if origin is Annotated:
        return _annotation_schema(arguments[0], seen=seen, depth=depth)
    if origin in (Union, types.UnionType):
        options: list[JSONSchema] = []
        reasons: list[str] = []
        supported_non_null = 0
        for option in arguments:
            supported, option_schema, reason = _annotation_schema(
                option,
                seen=seen,
                depth=depth + 1,
            )
            if supported:
                options.append(option_schema)
                if option not in (None, _NONE_TYPE):
                    supported_non_null += 1
            elif reason is not None:
                reasons.append(reason)
        if not options or (supported_non_null == 0 and reasons):
            return False, {}, "; ".join(reasons) or "union has no JSON-compatible branch"
        if len(options) == 1:
            return True, options[0], None
        return True, {"anyOf": options}, None
    if origin is Literal:
        return True, _json_scalar_schema(arguments), None
    if origin in (Callable,):
        return False, {}, "Python callbacks cannot be supplied through JSON"
    if origin in (Mapping, dict):
        key_annotation, value_annotation = arguments or (str, Any)
        key_supported = key_annotation in (str, Any, object)
        value_supported, value_schema, reason = _annotation_schema(
            value_annotation,
            seen=seen,
            depth=depth + 1,
        )
        if not key_supported or not value_supported:
            return False, {}, reason or "mapping keys must be strings"
        return True, {"type": "object", "additionalProperties": value_schema}, None
    if origin in (list, Sequence, Iterable, Iterator, set, frozenset, Set):
        item_annotation = arguments[0] if arguments else Any
        supported, item_schema, reason = _annotation_schema(
            item_annotation,
            seen=seen,
            depth=depth + 1,
        )
        if not supported:
            return False, {}, reason
        schema = {"type": "array", "items": item_schema}
        if origin in (set, frozenset, Set):
            schema["uniqueItems"] = True
        return True, schema, None
    if origin is tuple:
        if not arguments:
            return True, {"type": "array"}, None
        if len(arguments) == 2 and arguments[1] is Ellipsis:
            supported, item_schema, reason = _annotation_schema(
                arguments[0],
                seen=seen,
                depth=depth + 1,
            )
            if not supported:
                return False, {}, reason
            return True, {"type": "array", "items": item_schema}, None
        prefix: list[JSONSchema] = []
        for item in arguments:
            supported, item_schema, reason = _annotation_schema(
                item,
                seen=seen,
                depth=depth + 1,
            )
            if not supported:
                return False, {}, reason
            prefix.append(item_schema)
        return (
            True,
            {
                "type": "array",
                "prefixItems": prefix,
                "minItems": len(prefix),
                "maxItems": len(prefix),
            },
            None,
        )
    if origin is os.PathLike:
        return True, {"type": "string", "format": "path"}, None
    if origin is not None:
        origin_text = f"{origin!s}".lower()
        if "callable" in origin_text or "protocol" in origin_text or "io" in origin_text:
            return False, {}, f"{origin!s} requires a live Python object"
        return False, {}, f"unsupported generic type {origin!s}"

    if annotation is bool:
        return True, {"type": "boolean"}, None
    if annotation is int:
        return True, {"type": "integer"}, None
    if annotation is float:
        return True, {"type": "number"}, None
    if annotation is str:
        return True, {"type": "string"}, None
    if annotation is bytes:
        return True, {"type": "string", "contentEncoding": "utf-8"}, None
    if annotation is DNASequence:
        return (
            True,
            {
                "anyOf": [
                    {"type": "string", "description": "DNA sequence text"},
                    {
                        "type": "object",
                        "description": "A DNA record mapping accepted by dnakit.DNA",
                    },
                ]
            },
            None,
        )
    if annotation is DNARecord:
        return (
            True,
            {
                "type": "object",
                "required": ["sequence", "id"],
                "properties": {
                    "sequence": {"type": "string"},
                    "id": {"type": "string"},
                    "description": {"type": "string"},
                    "metadata": {"type": "object"},
                },
            },
            None,
        )
    if annotation is DNASet:
        return (
            True,
            {
                "anyOf": [
                    {"type": "array", "items": {"type": "object"}},
                    {"type": "object", "description": "DNA collection mapping"},
                ]
            },
            None,
        )
    if annotation is DNA:
        return True, {"description": "Any JSON input accepted by dnakit.DNA"}, None
    if inspect.isclass(annotation) and issubclass(annotation, Enum):
        return True, _json_scalar_schema(tuple(item.value for item in annotation)), None
    if inspect.isclass(annotation) and issubclass(annotation, PurePath):
        return True, {"type": "string", "format": "path"}, None
    if inspect.isclass(annotation) and issubclass(annotation, Mapping):
        return True, {"type": "object", "additionalProperties": {}}, None
    if inspect.isclass(annotation) and annotation.__module__.startswith("dnakit."):
        if annotation.__module__ == "dnakit.io.tables" and annotation.__qualname__ == "DictResult":
            return True, {"type": "object"}, None
        if inspect.isabstract(annotation) or getattr(annotation, "_is_protocol", False):
            return False, {}, f"{annotation.__qualname__} requires a live Python implementation"
        return _class_schema(annotation, seen=seen, depth=depth)
    return False, {}, f"{annotation!s} requires a live Python object"


def _format_annotation(annotation: object) -> str:
    if annotation is inspect.Parameter.empty:
        return "Any"
    if isinstance(annotation, str):
        return annotation
    return inspect.formatannotation(annotation)


def analyze_callable(function: Callable[..., object] | type[object]) -> CallPlan:
    """Resolve one callable into a bounded JSON input schema and conversion plan."""

    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError) as exc:
        return CallPlan(
            signature="(?)",
            parameters=(),
            input_schema={"type": "object"},
            output_schema={},
            accepts_keyword_arguments=False,
            compatible=False,
            incompatibilities=(f"signature unavailable: {exc}",),
        )
    hints = _resolve_hints(function)
    plans: list[ParameterPlan] = []
    properties: dict[str, object] = {}
    required_names: list[str] = []
    incompatibilities: list[str] = []
    accepts_keyword_arguments = False
    for parameter in signature.parameters.values():
        if parameter.kind is inspect.Parameter.VAR_KEYWORD:
            accepts_keyword_arguments = True
            continue
        if parameter.kind is inspect.Parameter.VAR_POSITIONAL:
            incompatibilities.append(f"{parameter.name}: variable positional arguments")
            continue
        annotation = hints.get(parameter.name, parameter.annotation)
        supported, schema, reason = _annotation_schema(annotation)
        required = parameter.default is inspect.Parameter.empty
        has_default = not required
        if supported and has_default:
            default = _default_for_schema(parameter.default)
            if default is not None or parameter.default is None:
                schema = dict(schema)
                schema["default"] = default
        plan = ParameterPlan(
            name=parameter.name,
            annotation=annotation,
            kind=parameter.kind.name,
            required=required,
            supported=supported,
            schema=schema,
            default=parameter.default,
            has_default=has_default,
            reason=reason,
        )
        plans.append(plan)
        if supported:
            properties[parameter.name] = schema
            if required:
                required_names.append(parameter.name)
        elif required:
            incompatibilities.append(f"{parameter.name}: {reason or 'unsupported input type'}")

    input_schema: JSONSchema = {
        "type": "object",
        "properties": properties,
        "additionalProperties": accepts_keyword_arguments,
    }
    if required_names:
        input_schema["required"] = required_names
    return_annotation = hints.get("return", signature.return_annotation)
    output_supported, output_schema, _ = _annotation_schema(return_annotation)
    if not output_supported:
        output_schema = {}
    return CallPlan(
        signature=str(signature),
        parameters=tuple(plans),
        input_schema=input_schema,
        output_schema=output_schema,
        accepts_keyword_arguments=accepts_keyword_arguments,
        compatible=not incompatibilities,
        incompatibilities=tuple(incompatibilities),
    )


def _construct_class(value: object, annotation: type[object], *, path: str) -> object:
    if isinstance(value, annotation):
        return value
    if not isinstance(value, Mapping):
        raise _error(f"{path} must be an object for {_format_annotation(annotation)}.", path=path)
    plan = analyze_callable(annotation)
    if not plan.compatible:
        raise _error(
            f"{_format_annotation(annotation)} cannot be constructed from JSON.",
            path=path,
        )
    positional, keywords = coerce_arguments(
        plan,
        {str(key): item for key, item in value.items()},
        path=path,
    )
    try:
        return annotation(*positional, **keywords)
    except ConfigurationError:
        raise
    except (TypeError, ValueError) as exc:
        raise _error(
            f"Could not construct {_format_annotation(annotation)}: {exc}",
            path=path,
        ) from exc


def coerce_value(value: object, annotation: object, *, path: str) -> object:
    """Convert one JSON-compatible value to a public DNAKit annotation."""

    if annotation is inspect.Parameter.empty or annotation is Any or annotation is object:
        return value
    if annotation is None or annotation is _NONE_TYPE:
        if value is None:
            return None
        raise _error(f"{path} must be null.", path=path, value=value)
    if isinstance(annotation, str):
        raise _error(f"{path} has an unsupported Python-only type.", path=path, value=value)
    if isinstance(annotation, TypeVar):
        if annotation.__bound__ is not None:
            return coerce_value(value, annotation.__bound__, path=path)
        if annotation.__constraints__:
            for constraint in annotation.__constraints__:
                try:
                    return coerce_value(value, constraint, path=path)
                except ConfigurationError:
                    continue
            raise _error(f"{path} does not match the constrained type.", path=path, value=value)
        return value

    origin = get_origin(annotation)
    arguments = get_args(annotation)
    if origin is Annotated:
        return coerce_value(value, arguments[0], path=path)
    if origin in (Union, types.UnionType):
        failures: list[str] = []
        for option in arguments:
            supported, _, _ = _annotation_schema(option)
            if not supported:
                continue
            try:
                return coerce_value(value, option, path=path)
            except ConfigurationError as exc:
                failures.append(exc.message)
        detail = failures[-1] if failures else "no JSON-compatible union branch"
        raise _error(f"{path} does not match an accepted type: {detail}", path=path, value=value)
    if origin is Literal:
        if value not in arguments:
            raise _error(f"{path} must be one of {arguments!r}.", path=path, value=value)
        return value
    if origin in (Mapping, dict):
        if not isinstance(value, Mapping):
            raise _error(f"{path} must be an object.", path=path, value=value)
        key_annotation, value_annotation = arguments or (str, Any)
        return {
            coerce_value(key, key_annotation, path=f"{path}.<key>"): coerce_value(
                item,
                value_annotation,
                path=f"{path}.{key}",
            )
            for key, item in value.items()
        }
    if origin in (list, Sequence, Iterable, Iterator, set, frozenset, Set):
        if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Iterable):
            raise _error(f"{path} must be an array.", path=path, value=value)
        item_annotation = arguments[0] if arguments else Any
        converted = [
            coerce_value(item, item_annotation, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
        if origin is list:
            return converted
        if origin in (set, Set):
            return set(converted)
        if origin is frozenset:
            return frozenset(converted)
        if origin is Iterator:
            return iter(converted)
        return tuple(converted)
    if origin is tuple:
        if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Iterable):
            raise _error(f"{path} must be an array.", path=path, value=value)
        values = list(value)
        if not arguments:
            return tuple(values)
        if len(arguments) == 2 and arguments[1] is Ellipsis:
            return tuple(
                coerce_value(item, arguments[0], path=f"{path}[{index}]")
                for index, item in enumerate(values)
            )
        if len(values) != len(arguments):
            raise _error(
                f"{path} must contain exactly {len(arguments)} items.",
                path=path,
                value=value,
            )
        return tuple(
            coerce_value(item, item_annotation, path=f"{path}[{index}]")
            for index, (item, item_annotation) in enumerate(zip(values, arguments, strict=True))
        )
    if origin is os.PathLike:
        if not isinstance(value, str):
            raise _error(f"{path} must be a filesystem path string.", path=path, value=value)
        return Path(value)

    if annotation is bool:
        if not isinstance(value, bool):
            raise _error(f"{path} must be a boolean.", path=path, value=value)
        return value
    if annotation is int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise _error(f"{path} must be an integer.", path=path, value=value)
        return value
    if annotation is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise _error(f"{path} must be a number.", path=path, value=value)
        return float(value)
    if annotation is str:
        if not isinstance(value, str):
            raise _error(f"{path} must be a string.", path=path, value=value)
        return value
    if annotation is bytes:
        if not isinstance(value, str):
            raise _error(f"{path} must be a UTF-8 string.", path=path, value=value)
        return value.encode("utf-8")
    if annotation is DNA:
        return DNA(value)
    if annotation is DNASequence:
        return DNA(value).sequence
    if annotation is DNARecord:
        return DNA(value).record
    if annotation is DNASet:
        return DNA(value).dataset
    if inspect.isclass(annotation) and issubclass(annotation, Enum):
        try:
            return annotation(value)
        except (TypeError, ValueError) as exc:
            raise _error(
                f"{path} is not a valid {_format_annotation(annotation)} value.",
                path=path,
                value=value,
            ) from exc
    if inspect.isclass(annotation) and issubclass(annotation, PurePath):
        if not isinstance(value, str):
            raise _error(f"{path} must be a filesystem path string.", path=path, value=value)
        return annotation(value)
    if inspect.isclass(annotation) and issubclass(annotation, Mapping):
        if not isinstance(value, Mapping):
            raise _error(f"{path} must be an object.", path=path, value=value)
        try:
            mapping_factory = cast(Callable[[object], object], annotation)
            return mapping_factory(value)
        except (TypeError, ValueError) as exc:
            raise _error(
                f"Could not construct {_format_annotation(annotation)}: {exc}",
                path=path,
                value=value,
            ) from exc
    if inspect.isclass(annotation) and annotation.__module__.startswith("dnakit."):
        if annotation.__module__ == "dnakit.io.tables" and annotation.__qualname__ == "DictResult":
            if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
                raise _error(f"{path} must be an object with string keys.", path=path, value=value)
            return _AgentDictResult({str(key): item for key, item in value.items()})
        return _construct_class(value, annotation, path=path)
    raise _error(f"{path} requires a live Python object.", path=path, value=value)


def coerce_arguments(
    plan: CallPlan,
    arguments: Mapping[str, object],
    *,
    path: str = "arguments",
) -> tuple[tuple[object, ...], dict[str, object]]:
    """Validate a JSON argument mapping and construct positional/keyword inputs."""

    if not isinstance(arguments, Mapping) or any(not isinstance(key, str) for key in arguments):
        raise ConfigurationError(
            "Agent tool arguments must be an object with string keys.",
            code="AGENT_ARGUMENTS_INVALID",
        )
    parameters = {parameter.name: parameter for parameter in plan.parameters}
    unsupported_supplied = sorted(
        name for name in arguments if name in parameters and not parameters[name].supported
    )
    if unsupported_supplied:
        raise ConfigurationError(
            "Python-only parameters cannot be supplied through the Agent interface.",
            code="AGENT_PARAMETER_UNSUPPORTED",
            context={"parameters": unsupported_supplied},
        )
    unknown = sorted(name for name in arguments if name not in parameters)
    if unknown and not plan.accepts_keyword_arguments:
        raise ConfigurationError(
            "Agent tool arguments contain unknown parameters.",
            code="AGENT_ARGUMENT_UNKNOWN",
            context={"parameters": unknown},
        )
    missing = sorted(
        parameter.name
        for parameter in plan.parameters
        if parameter.required and parameter.name not in arguments
    )
    if missing:
        raise ConfigurationError(
            "Agent tool arguments are missing required parameters.",
            code="AGENT_ARGUMENT_REQUIRED",
            context={"parameters": missing},
        )

    positional_plans = [
        parameter for parameter in plan.parameters if parameter.kind == "POSITIONAL_ONLY"
    ]
    last_positional = max(
        (index for index, parameter in enumerate(positional_plans) if parameter.name in arguments),
        default=-1,
    )
    positional: list[object] = []
    keywords: dict[str, object] = {}
    for index, parameter in enumerate(positional_plans):
        if index > last_positional:
            break
        raw = arguments.get(parameter.name, parameter.default)
        positional.append(coerce_value(raw, parameter.annotation, path=f"{path}.{parameter.name}"))
    for parameter in plan.parameters:
        if parameter.kind == "POSITIONAL_ONLY" or parameter.name not in arguments:
            continue
        keywords[parameter.name] = coerce_value(
            arguments[parameter.name],
            parameter.annotation,
            path=f"{path}.{parameter.name}",
        )
    if plan.accepts_keyword_arguments:
        for name in unknown:
            keywords[name] = arguments[name]
    return tuple(positional), keywords


__all__ = [
    "CallPlan",
    "JSONSchema",
    "ParameterPlan",
    "analyze_callable",
    "coerce_arguments",
    "coerce_value",
]
