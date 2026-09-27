"""Compatibility repairs preserve public binding and schema contracts."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools.compatibility import check_mcp_compatibility, check_python_compatibility


def _python(
    parameters: str, *, returns: str = "int", body: str = "return 1", name: str = "measure"
) -> str:
    return f"def {name}({parameters}) -> {returns}:\n    {body}\n"


def test_python_implementation_fix_and_optional_keyword_parameter_are_compatible() -> None:
    original = _python(
        "sequence: str, /, scale: int = 1, *values: int, strict: bool = True, **options: str"
    )
    candidate = _python(
        "sequence: str, /, scale: int = 1, *values: int, strict: bool = True, "
        "limit: int | None = None, **options: str",
        body="return len(sequence) * scale",
    )
    assert check_python_compatibility(original, candidate, "measure") == {
        "compatible": True,
        "added_optional_parameters": ["limit"],
    }


@pytest.mark.parametrize(
    "candidate",
    [
        "",
        "renamed: int = 1",
        "value: str = 1",
        "value: int = 2",
        "value: int",
        "value: int = 1, extra: int = 0",
        "value: int = 1, *, extra: int",
        "value: int = 1, *, extra: int = ...",
        "value: int = 1, /",
        "*, value: int = 1",
        "value: int = 1, *extra",
        "value: int = 1, **extra",
    ],
)
def test_python_breaking_parameters_are_rejected(candidate: str) -> None:
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_python_compatibility(_python("value: int = 1"), _python(candidate), "measure")


@pytest.mark.parametrize(
    "before,after",
    [
        ("first, second", "second, first"),
        ("first, /, second", "first, second"),
        ("first, *args, **kwargs", "first, **kwargs"),
        ("first, *args, **kwargs", "first, *args"),
        ("first, *args, **kwargs", "first, *args, **renamed"),
        ("*, first=1, second=2", "*, second=2, first=1"),
        ("first=True", "first=1"),
        ("first={'description': 'old'}", "first={'description': 'new'}"),
    ],
)
def test_python_order_variadic_names_and_exact_defaults_are_preserved(
    before: str, after: str
) -> None:
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_python_compatibility(_python(before), _python(after), "measure")


@pytest.mark.parametrize(
    "candidate",
    [
        _python("value", returns="str"),
        _python("value", name="renamed"),
        "async " + _python("value"),
        _python("value", body="yield 1"),
        _python("value") + _python("value"),
        "class Wrapper:\n    def measure(self, value):\n        return 1\n",
        _python("value") + "measure = other\n",
        _python("value") + "if True:\n    measure = other\n",
        _python("value") + "public_alias = measure\n",
        _python("value") + _python("value", name="_helper") + "new_public_alias = _helper\n",
        _python("value") + "__all__.clear()\n",
        _python("value") + _python("value", name="new_public_function"),
        "@new_decorator\n" + _python("value"),
    ],
)
def test_python_return_kind_exports_and_target_identity_are_preserved(candidate: str) -> None:
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_python_compatibility(_python("value"), candidate, "measure")


def test_python_check_does_not_execute_imports_or_module_code(tmp_path: Path) -> None:
    marker = tmp_path / "must-not-exist"
    prefix = f"import missing_dependency\nopen({str(marker)!r}, 'w').write('executed')\n"
    assert check_python_compatibility(
        prefix + _python("value"), prefix + _python("value", body="return 2"), "measure"
    )["compatible"]
    assert not marker.exists()


def test_python_nested_generator_does_not_change_outer_function_kind() -> None:
    candidate = "def measure(value) -> int:\n    def helper():\n        yield value\n    return 2\n"
    assert check_python_compatibility(_python("value"), candidate, "measure")["compatible"]


def test_python_async_generator_can_be_repaired_without_changing_its_kind() -> None:
    original = "async def measure(value) -> int:\n    yield value\n"
    candidate = "async def measure(value, *, enabled=True) -> int:\n    yield value + 1\n"
    assert check_python_compatibility(original, candidate, "measure")[
        "added_optional_parameters"
    ] == ["enabled"]


def test_python_type_comments_are_part_of_the_contract() -> None:
    original = "def measure(value):  # type: (int) -> int\n    return 1\n"
    candidate = original.replace("(int) -> int", "(str) -> int")
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_python_compatibility(original, candidate, "measure")


def test_other_public_function_cannot_be_changed_during_target_repair() -> None:
    original = _python("value") + _python("other=1", name="existing")
    candidate = _python("value", body="return 2") + _python("other=2", name="existing")
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_python_compatibility(original, candidate, "measure")


def _mcp() -> dict[str, Any]:
    return {
        "name": "measure",
        "description": "Original documentation",
        "inputSchema": {
            "type": "object",
            "properties": {
                "sequence": {"type": "string", "minLength": 1, "description": "DNA text"},
                "scale": {"type": "integer", "default": 1, "minimum": 1},
            },
            "required": ["sequence"],
            "additionalProperties": False,
        },
        "outputSchema": {
            "type": "object",
            "properties": {"length": {"type": "integer"}},
            "required": ["length"],
        },
        "annotations": {"readOnlyHint": True, "destructiveHint": False},
    }


def test_mcp_allows_optional_defaulted_parameters_and_documentation_without_mutation() -> None:
    original = _mcp()
    candidate = deepcopy(original)
    candidate["description"] = "Clarified documentation"
    candidate["inputSchema"]["properties"]["sequence"]["description"] = "Canonical sequence text"
    candidate["inputSchema"]["properties"]["limit"] = {"type": ["integer", "null"], "default": None}
    candidate["outputSchema"]["properties"]["length"]["examples"] = [4]
    snapshots = deepcopy((original, candidate))
    assert check_mcp_compatibility(original, candidate) == {
        "compatible": True,
        "added_optional_parameters": ["limit"],
    }
    assert (original, candidate) == snapshots


@pytest.mark.parametrize(
    "change",
    [
        "name",
        "remove_parameter",
        "rename_parameter",
        "type",
        "default",
        "required",
        "remove_required",
        "new_required",
        "new_without_default",
        "input_constraint",
        "output_constraint",
        "output_type",
        "remove_output",
        "effect",
        "title",
    ],
)
def test_mcp_breaking_contract_changes_are_rejected(change: str) -> None:
    original = _mcp()
    candidate = deepcopy(original)
    schema = candidate["inputSchema"]
    properties = schema["properties"]
    if change == "name":
        candidate["name"] = "renamed"
    elif change == "remove_parameter":
        del properties["scale"]
    elif change == "rename_parameter":
        properties["factor"] = properties.pop("scale")
    elif change == "type":
        properties["scale"]["type"] = "number"
    elif change == "default":
        properties["scale"]["default"] = True
    elif change == "required":
        schema["required"].append("scale")
    elif change == "remove_required":
        schema["required"] = []
    elif change in {"new_required", "new_without_default"}:
        properties["limit"] = {"type": "integer"}
        if change == "new_required":
            properties["limit"]["default"] = 10
            schema["required"].append("limit")
    elif change == "input_constraint":
        schema["additionalProperties"] = True
    elif change == "output_constraint":
        candidate["outputSchema"]["required"] = []
    elif change == "output_type":
        candidate["outputSchema"]["properties"]["length"]["type"] = "string"
    elif change == "remove_output":
        del candidate["outputSchema"]
    elif change == "effect":
        candidate["annotations"]["destructiveHint"] = True
    else:
        candidate["title"] = "Changed functional metadata"
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_mcp_compatibility(original, candidate)


@pytest.mark.parametrize("keyword", ["default", "enum", "const"])
@pytest.mark.parametrize("business_key", ["description", "examples"])
def test_mcp_documentation_filter_preserves_business_values(
    keyword: str, business_key: str
) -> None:
    original = _mcp()
    value: Any = {business_key: "original business value"}
    original["inputSchema"]["properties"]["scale"][keyword] = (
        [value] if keyword == "enum" else value
    )
    candidate = deepcopy(original)
    changed = candidate["inputSchema"]["properties"]["scale"][keyword]
    (changed[0] if keyword == "enum" else changed)[business_key] = "modified business value"
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_mcp_compatibility(original, candidate)


def test_mcp_nested_schema_documentation_is_ignored_but_property_names_are_not() -> None:
    original = _mcp()
    original["outputSchema"] = {
        "type": "object",
        "properties": {"description": {"type": "string", "description": "Document text"}},
        "$defs": {"examples": {"type": "string", "description": "A type definition"}},
        "allOf": [{"description": "Branch text", "required": ["description"]}],
    }
    candidate = deepcopy(original)
    candidate["outputSchema"]["properties"]["description"]["description"] = "Clearer text"
    candidate["outputSchema"]["$defs"]["examples"]["description"] = "Clearer type definition"
    candidate["outputSchema"]["allOf"][0]["description"] = "Clearer branch"
    assert check_mcp_compatibility(original, candidate)["compatible"]
    del candidate["outputSchema"]["properties"]["description"]
    with pytest.raises(ConfigurationError, match="TOOL_REPAIR_BREAKING_CHANGE"):
        check_mcp_compatibility(original, candidate)


def test_mcp_optional_property_can_be_added_to_an_empty_object_contract() -> None:
    original: dict[str, Any] = {"name": "measure", "inputSchema": {"type": "object"}}
    candidate = deepcopy(original)
    candidate["inputSchema"]["properties"] = {"limit": {"type": "integer", "default": 0}}
    assert check_mcp_compatibility(original, candidate)["added_optional_parameters"] == ["limit"]


def test_mcp_error_identifies_the_changed_constraint() -> None:
    original = _mcp()
    candidate = deepcopy(original)
    candidate["outputSchema"]["properties"]["length"]["type"] = "string"
    with pytest.raises(ConfigurationError, match=r"outputSchema\.properties\.length\.type"):
        check_mcp_compatibility(original, candidate)
