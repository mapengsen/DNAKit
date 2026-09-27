"""Conservative, execution-free interface checks for implementation repairs.

Python checks compare declared top-level interfaces, not runtime monkey-patching
or semantic equivalence of annotation/default expressions. MCP checks preserve
JSON Schema constraints and ignore documentation only at schema locations.
"""

from __future__ import annotations

import ast
import json
from copy import deepcopy
from typing import Any, NoReturn

from dnakit.exceptions import ConfigurationError

_FUNCTIONS = (ast.FunctionDef, ast.AsyncFunctionDef)
_SCHEMA_MAPS = frozenset(
    {"properties", "patternProperties", "$defs", "definitions", "dependentSchemas"}
)
_SCHEMA_SINGLE = frozenset(
    {
        "items",
        "additionalItems",
        "additionalProperties",
        "contains",
        "not",
        "if",
        "then",
        "else",
        "propertyNames",
        "unevaluatedItems",
        "unevaluatedProperties",
        "contentSchema",
    }
)
_SCHEMA_LISTS = frozenset({"allOf", "anyOf", "oneOf", "prefixItems"})


def _breaking(reason: str) -> NoReturn:
    raise ConfigurationError(reason, code="TOOL_REPAIR_BREAKING_CHANGE")


def _expression(node: ast.AST | None) -> str | None:
    return ast.dump(node, include_attributes=False) if node is not None else None


def _parse(source: str, label: str) -> ast.Module:
    if not isinstance(source, str):
        _breaking(f"{label} source must be Python text.")
    try:
        return ast.parse(source, type_comments=True)
    except (SyntaxError, ValueError, TypeError, RecursionError) as exc:
        _breaking(f"{label} source cannot be parsed: {exc}.")


class _Bindings(ast.NodeVisitor):
    """Find protected names rebound at module scope, including control-flow blocks."""

    def __init__(self, protected: set[str], definitions: set[int]) -> None:
        self.protected = protected
        self.definitions = definitions
        self.rebound: set[str] = set()

    def visit_Name(self, node: ast.Name) -> None:
        if isinstance(node.ctx, (ast.Store, ast.Del)) and node.id in self.protected:
            self.rebound.add(node.id)

    def visit_FunctionDef(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        if id(node) not in self.definitions and node.name in self.protected:
            self.rebound.add(node.name)

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        if node.name in self.protected:
            self.rebound.add(node.name)

    def visit_Lambda(self, node: ast.Lambda) -> None:
        pass

    def visit_Import(self, node: ast.Import | ast.ImportFrom) -> None:
        for alias in node.names:
            name = alias.asname or alias.name.split(".")[0]
            if name in self.protected or alias.name == "*":
                self.rebound.add(name)

    visit_ImportFrom = visit_Import


def _functions(tree: ast.Module, function_name: str, label: str) -> dict[str, Any]:
    definitions = [node for node in tree.body if isinstance(node, _FUNCTIONS)]
    named = [node for node in definitions if node.name == function_name]
    if len(named) != 1:
        _breaking(
            f"{label} must define exactly one top-level function named {function_name!r}; "
            "class methods and assignment/import aliases are unsupported."
        )
    public: dict[str, Any] = {}
    for node in definitions:
        if not node.name.startswith("_") or node.name == function_name:
            if node.name in public:
                _breaking(f"{label} defines public function {node.name!r} more than once.")
            public[node.name] = node
    bindings = _Bindings(set(public), {id(node) for node in definitions})
    bindings.visit(tree)
    if bindings.rebound:
        _breaking(f"{label} rebinds public function names: {sorted(bindings.rebound)}.")
    return public


class _Generator(ast.NodeVisitor):
    def __init__(self) -> None:
        self.found = False

    def visit_Yield(self, node: ast.Yield | ast.YieldFrom) -> None:
        self.found = True

    visit_YieldFrom = visit_Yield

    def visit_FunctionDef(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        pass

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        pass

    def visit_Lambda(self, node: ast.Lambda) -> None:
        pass


def _is_generator(function: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    visitor = _Generator()
    for node in function.body:
        visitor.visit(node)
    return visitor.found


def _parameters(
    function: ast.FunctionDef | ast.AsyncFunctionDef,
) -> list[tuple[str, str, str | None, str | None]]:
    arguments = function.args
    positional = [*arguments.posonlyargs, *arguments.args]
    defaults: list[ast.expr | None] = [None] * (len(positional) - len(arguments.defaults))
    defaults.extend(arguments.defaults)
    result = [
        (
            argument.arg,
            "positional_only" if index < len(arguments.posonlyargs) else "positional_or_keyword",
            _expression(argument),
            _expression(default),
        )
        for index, (argument, default) in enumerate(zip(positional, defaults, strict=True))
    ]
    if arguments.vararg:
        argument = arguments.vararg
        result.append((argument.arg, "var_positional", _expression(argument), None))
    result.extend(
        (argument.arg, "keyword_only", _expression(argument), _expression(default))
        for argument, default in zip(arguments.kwonlyargs, arguments.kw_defaults, strict=True)
    )
    if arguments.kwarg:
        argument = arguments.kwarg
        result.append((argument.arg, "var_keyword", _expression(argument), None))
    names = [item[0] for item in result]
    if len(names) != len(set(names)):
        _breaking(f"Function {function.name!r} has duplicate parameter names.")
    return result


def _check_function(original: Any, candidate: Any) -> list[str]:
    name = original.name
    if type(original) is not type(candidate) or _is_generator(original) != _is_generator(candidate):
        _breaking(f"Function {name!r} cannot change async/generator kind.")
    for field in ("returns", "type_comment"):
        left, right = getattr(original, field), getattr(candidate, field)
        if (_expression(left) if isinstance(left, ast.AST) else left) != (
            _expression(right) if isinstance(right, ast.AST) else right
        ):
            _breaking(f"Function {name!r} cannot change its return annotation or type contract.")
    if [_expression(node) for node in original.decorator_list] != [
        _expression(node) for node in candidate.decorator_list
    ]:
        _breaking(f"Function {name!r} cannot change decorators during a compatible repair.")
    before, after = _parameters(original), _parameters(candidate)
    known = {item[0] for item in before}
    if [item for item in after if item[0] in known] != before:
        _breaking(
            f"Function {name!r} changed existing parameter names, order, kinds, "
            "annotations or defaults."
        )
    added = [item for item in after if item[0] not in known]
    ellipsis = _expression(ast.Constant(value=Ellipsis))
    for parameter, kind, _, default in added:
        if kind != "keyword_only" or default is None or default == ellipsis:
            _breaking(
                f"New parameter {parameter!r} must be keyword-only with a real default value."
            )
    return [item[0] for item in added]


def _exports(tree: ast.Module) -> list[str | None]:
    return [
        _expression(node)
        for node in tree.body
        if not isinstance(node, (*_FUNCTIONS, ast.ClassDef))
        and any(isinstance(child, ast.Name) and child.id == "__all__" for child in ast.walk(node))
    ]


def _public_aliases(tree: ast.Module) -> dict[str, str]:
    aliases = {node.name: node.name for node in tree.body if isinstance(node, _FUNCTIONS)}
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and isinstance(node.value, ast.Name):
            target = aliases.get(node.value.id)
            if target is not None:
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for item in targets:
                    if isinstance(item, ast.Name):
                        aliases[item.id] = target
    return {name: target for name, target in aliases.items() if not name.startswith("_")}


def check_python_compatibility(
    original_source: str,
    candidate_source: str,
    function_name: str,
) -> dict[str, object]:
    """Allow implementation repairs plus genuinely optional keyword-only additions."""
    if not isinstance(function_name, str) or not function_name.isidentifier():
        _breaking("function_name must name a top-level Python function, not a class method.")
    original = _parse(original_source, "Original")
    candidate = _parse(candidate_source, "Candidate")
    before = _functions(original, function_name, "Original")
    after = _functions(candidate, function_name, "Candidate")
    if (
        set(before) != set(after)
        or _exports(original) != _exports(candidate)
        or _public_aliases(original) != _public_aliases(candidate)
    ):
        _breaking(
            "A compatible repair cannot add, remove, rename or change exports of public functions."
        )
    additions: list[str] = []
    for name, function in before.items():
        added = _check_function(function, after[name])
        if name == function_name:
            additions = added
        elif added:
            _breaking(
                f"Optional parameters may only be added to the repaired function, not {name!r}."
            )
    return {"compatible": True, "added_optional_parameters": additions}


def _schema_contract(schema: object) -> object:
    if isinstance(schema, list):
        return [_schema_contract(item) for item in schema]
    if not isinstance(schema, dict):
        return deepcopy(schema)
    result: dict[str, Any] = {}
    for key, value in schema.items():
        if key in {"description", "examples"}:
            continue
        if key in _SCHEMA_MAPS and isinstance(value, dict):
            result[key] = {name: _schema_contract(child) for name, child in value.items()}
        elif key in _SCHEMA_SINGLE or key in _SCHEMA_LISTS:
            result[key] = _schema_contract(value)
        elif key == "dependencies" and isinstance(value, dict):
            result[key] = {
                name: _schema_contract(child) if isinstance(child, dict) else deepcopy(child)
                for name, child in value.items()
            }
        else:
            # enum, const and default contain business data, not nested schemas.
            result[key] = deepcopy(value)
    return result


def _json(value: object) -> str:
    try:
        return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)
    except (ValueError, TypeError, RecursionError):
        _breaking("MCP contracts must contain finite JSON values.")


def _difference(before: object, after: object, path: str = "tool") -> str:
    if isinstance(before, dict) and isinstance(after, dict):
        for key in dict.fromkeys([*before, *after]):
            selected = f"{path}.{key}"
            if key not in before or key not in after:
                return selected
            if _json(before[key]) != _json(after[key]):
                return _difference(before[key], after[key], selected)
    elif isinstance(before, list) and isinstance(after, list) and len(before) == len(after):
        for index, (left, right) in enumerate(zip(before, after, strict=True)):
            if _json(left) != _json(right):
                return _difference(left, right, f"{path}[{index}]")
    return path


def check_mcp_compatibility(
    original_tool: dict[str, Any],
    candidate_tool: dict[str, Any],
) -> dict[str, object]:
    """Preserve MCP contracts, allowing documentation and optional input additions."""
    if not isinstance(original_tool, dict) or not isinstance(candidate_tool, dict):
        _breaking("MCP tools must be JSON objects.")
    _json(original_tool)
    _json(candidate_tool)
    name = original_tool.get("name")
    if not isinstance(name, str) or not name or candidate_tool.get("name") != name:
        _breaking("The MCP tool name cannot change.")
    original = deepcopy(original_tool)
    candidate = deepcopy(candidate_tool)
    for tool in (original, candidate):
        tool.pop("description", None)
        tool.pop("examples", None)
        for field in ("inputSchema", "outputSchema"):
            if field in tool:
                tool[field] = _schema_contract(tool[field])
    before, after = original.get("inputSchema"), candidate.get("inputSchema")
    if not isinstance(before, dict) or not isinstance(after, dict):
        _breaking("MCP inputSchema must be an object schema.")
    old_properties = before.get("properties", {})
    new_properties = after.get("properties", {})
    if not isinstance(old_properties, dict) or not isinstance(new_properties, dict):
        _breaking("MCP input properties must be schema mappings.")
    removed = set(old_properties) - set(new_properties)
    if removed:
        _breaking(f"Existing MCP input parameters cannot be removed or renamed: {sorted(removed)}.")
    added = [parameter for parameter in new_properties if parameter not in old_properties]
    required = after.get("required", [])
    if not isinstance(required, list) or any(not isinstance(item, str) for item in required):
        _breaking("MCP required must be a list of parameter names.")
    for parameter in added:
        schema = new_properties[parameter]
        if parameter in required or not isinstance(schema, dict) or "default" not in schema:
            _breaking(f"New MCP parameter {parameter!r} must be optional and declare a default.")
        del new_properties[parameter]
    # An initially absent properties map may be introduced solely for additions.
    if "properties" not in before and not new_properties:
        after.pop("properties", None)
    if _json(original) != _json(candidate):
        _breaking(
            f"MCP contract changed at {_difference(original, candidate)}; "
            "only documentation and optional defaulted input additions may change."
        )
    return {"compatible": True, "added_optional_parameters": added}


__all__ = ["check_mcp_compatibility", "check_python_compatibility"]
