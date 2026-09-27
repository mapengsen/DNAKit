"""Prepare reviewable contribution files without writing or publishing them."""

# Bilingual contribution documentation intentionally uses Chinese punctuation.
# ruff: noqa: RUF001

from __future__ import annotations

import hashlib
import json
import keyword
import re
from pathlib import Path, PureWindowsPath
from typing import Any
from urllib.parse import parse_qsl, urlsplit

from dnakit.exceptions import ConfigurationError

_EFFECTS = {"local", "network", "filesystem", "model", "external_process"}
_SENSITIVE = re.compile(
    r"(?:^|_)(?:api_?key|access_?token|refresh_?token|token|password|passwd|secret|"
    r"credentials?|authorization|private_?key)(?:$|_)",
    re.IGNORECASE,
)
_CHANGE_ID = re.compile(r"[a-z][a-z0-9_]{0,40}-[a-f0-9]{24,64}")


def _invalid(message: str) -> ConfigurationError:
    return ConfigurationError(message, code="TOOL_CONTRIBUTION_INVALID")


def _json(value: object) -> str:
    try:
        return json.dumps(value, ensure_ascii=True, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError, RecursionError) as exc:
        raise _invalid("Contribution cases must contain finite JSON values.") from exc


def _inspect_data(value: object) -> bool:
    """Reject recognizable credentials and flag machine-specific fixture paths."""
    if isinstance(value, dict):
        absolute = False
        for key, child in value.items():
            if not isinstance(key, str):
                raise _invalid("Contribution JSON objects require string keys.")
            if _SENSITIVE.search(key) and child not in (None, "", False):
                raise _invalid("Replace credentials with public test fixtures before contributing.")
            absolute = _inspect_data(child) or absolute
        return absolute
    if isinstance(value, list):
        return any([_inspect_data(child) for child in value])
    if isinstance(value, str):
        if "-----BEGIN " in value and "PRIVATE KEY-----" in value:
            raise _invalid("Private keys cannot appear in contribution material.")
        if value.startswith(("http://", "https://")):
            parts = urlsplit(value)
            if (
                parts.username
                or parts.password
                or any(_SENSITIVE.search(key) for key, _ in parse_qsl(parts.query))
            ):
                raise _invalid("Credential-bearing URLs cannot appear in contribution material.")
        return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()
    if value is None or isinstance(value, (bool, int, float)):
        return False
    raise _invalid("Contribution fixtures must use JSON objects, arrays and scalar values.")


def _string_list(value: object, label: str) -> list[str]:
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", item)
        for item in value
    ):
        raise _invalid(f"{label} must contain module names or effect names.")
    return sorted(set(value))


def _literal(value: str, *, indent: str = "    ") -> str:
    return "\n".join(
        indent + json.dumps(value[index : index + 36], ensure_ascii=True)
        for index in range(0, len(value), 36)
    )


def _test_source(tool_name: str, cases: list[dict[str, Any]], *, gated: bool) -> bytes:
    source = '''"""Fixed Python and Agent cases; review fixture data before submission."""

from __future__ import annotations

import importlib
import inspect
import json
import os
from copy import deepcopy
from typing import Any

import pytest

# Keep bounded generated JSON literals wrapped, including large sequence fixtures.
# fmt: off
__TOOL__
_CASES = json.loads(
__CASES__
)
# fmt: on
_GATED = __GATED__
pytestmark = pytest.mark.skipif(
    _GATED and os.getenv("DNAKIT_RUN_GENERATED_EFFECT_TESTS") != "1",
    reason="Review dependencies, fixture data and effects; explicitly opt in to these cases.",
)


def _invoke(driver: str, arguments: dict[str, Any]) -> object:
    from dnakit.exceptions import ConfigurationError
    from dnakit.tools.catalog import _build_spec
    from dnakit.tools.registry import AgentToolRegistry, default_tool_registry
    from dnakit.tools.serialization import serialize_for_agent

    module_name, function_name = _TOOL_NAME.rsplit(".", 1)
    function = getattr(importlib.import_module(module_name), function_name)
    values = deepcopy(arguments)
    if "verbose" in inspect.signature(function).parameters:
        values.setdefault("verbose", False)
    if driver == "agent":
        registry = default_tool_registry(verbose=False)
        try:
            registry.resolve(_TOOL_NAME)
        except ConfigurationError as error:
            if error.code != "AGENT_TOOL_NOT_FOUND":
                raise
            registry = AgentToolRegistry((_build_spec(module_name, function_name, function),))
        return registry.execute(_TOOL_NAME, values, allow_side_effects=_GATED)
    return serialize_for_agent(function(**values), verbose=False)


@pytest.mark.parametrize("driver", ["python", "agent"])
@pytest.mark.parametrize("case", _CASES)
def test_public_contract(driver: str, case: dict[str, Any]) -> None:
    if "expected_error" in case:
        with pytest.raises(Exception) as caught:
            _invoke(driver, case["arguments"])
        error = caught.value
        expected = case["expected_error"]
        if "type" in expected:
            assert expected["type"] in {
                type(error).__name__,
                f"{type(error).__module__}.{type(error).__name__}",
            }
        if "code" in expected:
            assert getattr(error, "code", None) == expected["code"]
        return
    observed = _invoke(driver, case["arguments"])
    assert json.dumps(observed, sort_keys=True, allow_nan=False) == json.dumps(
        case["expected"], sort_keys=True, allow_nan=False
    )
'''
    replacements = {
        "__TOOL__": (
            "_TOOL_NAME = " + json.dumps(tool_name)
            if len(tool_name) <= 80
            else "_TOOL_NAME = (\n" + _literal(tool_name) + "\n)"
        ),
        "__CASES__": _literal(_json(cases)),
        "__GATED__": repr(gated),
    }
    return re.sub(
        r"__TOOL__|__CASES__|__GATED__", lambda match: replacements[match.group()], source
    ).encode("utf-8")


def contribution_files(
    *,
    root: Path,
    tool_name: str,
    cases: list[dict[str, Any]],
    source_paths: list[Path],
    reason: str,
    change_id: str,
    kind: str = "generation",
    provenance: dict[str, Any] | None = None,
    effects: list[str] | None = None,
    dependencies: list[str] | None = None,
) -> dict[Path, bytes]:
    """Return proposed tests/docs/metadata; callers activate files and users submit PRs.

    Missing effect metadata is treated conservatively. Recognizable credentials
    are rejected; arbitrary private user data cannot be detected automatically,
    so fixture review is always required before sharing these local materials.
    No tool is imported or run, and no filesystem/network/git mutation occurs.
    """
    if (
        not isinstance(tool_name, str)
        or len(tool_name) > 150
        or re.fullmatch(r"dnakit(?:\.[A-Za-z][A-Za-z0-9_]*){2,}", tool_name) is None
        or any(keyword.iskeyword(part) for part in tool_name.split("."))
    ):
        raise _invalid("tool_name must be a full public DNAKit Python function name.")
    if not isinstance(change_id, str) or _CHANGE_ID.fullmatch(change_id) is None:
        raise _invalid("change_id must be a development kind followed by its hexadecimal digest.")
    if not isinstance(kind, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", kind):
        raise _invalid("kind must be an identifier.")
    if not isinstance(reason, str) or not reason.strip() or len(reason) > 16_000:
        raise _invalid("reason must describe the concrete problem in bounded text.")
    _inspect_data(reason)
    if not isinstance(cases, list) or not 1 <= len(cases) <= 50:
        raise _invalid("Provide 1 to 50 fixed public regression cases.")
    for case in cases:
        if (
            not isinstance(case, dict)
            or not isinstance(case.get("arguments"), dict)
            or set(case) - {"arguments", "expected", "expected_error"}
            or ("expected" in case) == ("expected_error" in case)
        ):
            raise _invalid("Each case needs arguments and exactly one expected result or error.")
        if "expected_error" in case:
            error = case["expected_error"]
            if (
                not isinstance(error, dict)
                or not error
                or set(error) - {"type", "code"}
                or any(not isinstance(value, str) or not value.strip() for value in error.values())
            ):
                raise _invalid("expected_error needs an exception type and/or stable error code.")
    serialized_cases = _json(cases)
    if len(serialized_cases.encode()) > 1_000_000:
        raise _invalid("Contribution cases exceed 1 MB; use small public fixtures.")
    absolute_fixtures = _inspect_data(cases)
    provenance = {} if provenance is None else provenance
    if not isinstance(provenance, dict):
        raise _invalid("provenance must be an object.")
    declared_effects = _string_list(["external_process"] if effects is None else effects, "effects")
    if (
        not declared_effects
        or set(declared_effects) - _EFFECTS
        or ("local" in declared_effects and len(declared_effects) > 1)
    ):
        raise _invalid("Declare known effects; local cannot be combined with other effects.")
    declared_dependencies = _string_list(
        provenance.get("dependencies", []) if dependencies is None else dependencies,
        "dependencies",
    )
    public_provenance: dict[str, str] = {}
    for field in ("upstream_url", "upstream_commit", "license"):
        if field not in provenance:
            continue
        value = provenance[field]
        if not isinstance(value, str) or not value.strip() or len(value) > 4000:
            raise _invalid(f"provenance.{field} must contain actual bounded text.")
        if _inspect_data(value):
            raise _invalid("Contribution provenance must not contain absolute local paths.")
        if field == "upstream_url":
            url = urlsplit(value)
            if url.scheme not in {"https", "http"} or not url.netloc or url.query or url.fragment:
                raise _invalid("upstream_url must be a public HTTP(S) source URL without secrets.")
        public_provenance[field] = value
    resolved_root = root.resolve()
    if not source_paths:
        raise _invalid("Provide at least one source artifact path.")
    sources: list[str] = []
    for path in source_paths:
        target = path if path.is_absolute() else resolved_root / path
        try:
            relative = target.resolve().relative_to(resolved_root).as_posix()
        except ValueError as exc:
            raise _invalid("Source artifacts must be inside the DNAKit checkout.") from exc
        if relative == "." or relative.startswith((".dnakit/", ".git/")):
            raise _invalid(
                "Contribution source paths cannot include runtime audits or git internals."
            )
        sources.append(relative)
    safe_name = tool_name.replace(".", "_")
    digest = hashlib.sha256(change_id.encode()).hexdigest()[:12]
    test_path = Path(f"others/tests/generated/test_{safe_name}_{digest}.py")
    zh_path = Path(f"docs/llm_tools/generated/{safe_name}.md")
    en_path = Path(f"docs/llm_tools/generated/{safe_name}.en.md")
    metadata_path = Path(f"registry/contributions/{change_id}.json")
    pr_path = metadata_path.with_suffix(".md")
    gated = declared_effects != ["local"] or absolute_fixtures
    prefix = "DNAKIT_RUN_GENERATED_EFFECT_TESTS=1 " if gated else ""
    commands = [
        f"{prefix}python -m pytest {test_path.as_posix()} -q",
        f"python -m ruff check {test_path.as_posix()}",
        f"python -m mypy {test_path.as_posix()}",
    ]
    first = cases[0]
    example = json.dumps(first, ensure_ascii=False, indent=2, allow_nan=False)
    arguments = json.dumps(first["arguments"], ensure_ascii=False, allow_nan=False)
    module, name = tool_name.rsplit(".", 1)
    python_example = (
        "import json\n"
        f"from {module} import {name}\n\n"
        f"arguments = json.loads({arguments!r})\n"
        f"result = {name}(**arguments)\n"
        "print(result)"
    )
    mcp_example = json.dumps(
        {"tool_name": tool_name, "arguments": first["arguments"], "allow_side_effects": gated},
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    )
    compatibility = (
        "Keep the tool/function name and all existing parameter names, positions, types and "
        "defaults. Do not remove parameters or add required parameters. New parameters must be "
        "keyword-only with defaults. Preserve successful existing calls and the documented "
        "result contract; add a fixed regression case for each repaired behavior."
    )
    review = (
        "Review fixtures for private data, credentials, absolute paths and environment assumptions "
        "before sharing. Automated checks cannot establish that arbitrary user data is public. "
        "Replace local inputs with small public/synthetic fixtures. Review source, effects, "
        "dependencies and licenses. Submit a PR yourself; maintainers review and merge. "
        "Preparing these files does not publish, push or create a PR."
    )
    gating = (
        "These tests are skipped by default because effects are non-local/unknown or fixtures "
        "contain absolute paths. After preparing dependencies, disposable outputs and public "
        "fixtures, explicitly set DNAKIT_RUN_GENERATED_EFFECT_TESTS=1."
        if gated
        else "Declared effects are local; these cases run by default without model/network calls."
    )
    dependency_text = ", ".join(declared_dependencies) or "No additional dependencies declared."
    discovery = (
        "For a custom-registered target, the MCP server must receive a registry containing this "
        "tool. Review public discovery registration before submitting; ordinary generated tools "
        "are registered automatically. The regression fallback selects only this custom target "
        "and does not expand the default catalog."
    )
    source_text = "\n".join(f"- `{path}`" for path in sorted(set(sources)))
    validation = "\n\n".join(f"```bash\n{command}\n```" for command in commands)
    english = (
        f"# {tool_name}\n\n{reason}\n\n"
        "## Inputs and expected outcome\n\n"
        f"```json\n{example}\n```\n\n"
        "Cases declaring expected_error verify the exception type/code instead of a result.\n\n"
        f"## Python API\n\n```python\n{python_example}\n```\n\n"
        "## MCP\n\nUse the existing DNAKit `call_dnakit_tool` endpoint; "
        "no new server is needed.\n\n"
        f"{discovery}\n\n"
        f"```json\n{mcp_example}\n```\n\n"
        f"## Validation and contribution\n\n{gating}\n\nDependencies: {dependency_text}\n\n"
        "Every fixed case checks both the Python API and AgentToolRegistry execution, comparing "
        "serialized results or expected errors. Each case can therefore execute twice. "
        "The following commands are instructions, not claims that checks have run.\n\n"
        f"{validation}\n\n{compatibility}\n\n{review}\n"
    )
    chinese = (
        f"# {tool_name}\n\n用途与待解决问题：{reason}\n\n"
        f"## 输入与预期结果\n\n```json\n{example}\n```\n\n"
        "expected_error 用例检查异常类型或稳定错误码，不要求正常返回。\n\n"
        f"## Python 调用\n\n```python\n{python_example}\n```\n\n"
        "## MCP 调用\n\n通过现有 DNAKit `call_dnakit_tool` 调用，无需创建新服务器。\n\n"
        "自定义注册目标的 MCP 服务也必须使用包含该工具的 registry。提交前需审查公开发现登记；"
        "常规 generated 工具已自动登记。回归测试仅为未登记目标建立单工具 registry，"
        "不会扩大默认工具目录。\n\n"
        f"```json\n{mcp_example}\n```\n\n"
        "## 验证与贡献\n\n"
        + (
            "存在非本地/未声明副作用或绝对路径用例，测试默认跳过。准备依赖、可丢弃输出和公开用例后，"
            "显式设置 `DNAKIT_RUN_GENERATED_EFFECT_TESTS=1`。\n\n"
            if gated
            else "声明为本地计算，用例默认运行，不调用模型或网络。\n\n"
        )
        + f"额外依赖：{dependency_text}\n\n"
        "每个固定用例分别调用 Python API 和 AgentToolRegistry，比较序列化结果或预期异常；"
        "同一用例可能实际执行两次。以下是待执行的验证命令，不代表已经通过。\n\n"
        f"{validation}\n\n"
        "保留函数/工具名与全部旧参数的名称、位置、类型和默认值；不得删除参数或新增必填参数。"
        "新参数必须是带默认值的 keyword-only 参数；保留原成功调用与返回契约，"
        "为每个修复行为添加固定回归用例。\n\n"
        "提交前审查源码、副作用、依赖、许可证与用例数据；替换私人数据、凭据、绝对路径和环境依赖，"
        "使用小型公开或合成用例。自动检查不能证明任意用户数据可公开。"
        "由用户手动提交 PR，维护者审查合并；生成这些文件不会发布、push 或创建 PR。\n"
    )
    provenance_lines = "\n".join(f"- {key}: {value}" for key, value in public_provenance.items())
    pr_body = (
        f"# {tool_name}: {kind}\n\n{reason}\n\n"
        "## Change\n\n"
        f"Source files:\n\n{source_text}\n\n"
        f"Adds {len(cases)} fixed cases for Python and Agent execution, bilingual usage, and "
        "contribution metadata. Source behavior is reviewed against the recorded cases.\n\n"
        f"## Compatibility\n\n{compatibility}\n\n{discovery}\n\n"
        f"## Validation commands (not yet claimed as run)\n\n{validation}\n\n{gating}\n\n"
        f"Dependencies: {dependency_text}\n\n"
        + (f"## Recorded provenance\n\n{provenance_lines}\n\n" if provenance_lines else "")
        + f"## Before submitting\n\n{review}\n"
    )
    files = {
        test_path: _test_source(tool_name, cases, gated=gated),
        zh_path: chinese.encode("utf-8"),
        en_path: english.encode("utf-8"),
        pr_path: pr_body.encode("utf-8"),
    }
    metadata = {
        "schema_version": 1,
        "change_id": change_id,
        "kind": kind,
        "tool_name": tool_name,
        "source_paths": sorted(set(sources)),
        "files": sorted([*(path.as_posix() for path in files), metadata_path.as_posix()]),
        "case_count": len(cases),
        "cases_sha256": hashlib.sha256(serialized_cases.encode()).hexdigest(),
        "effects": declared_effects,
        "dependencies": declared_dependencies,
        "provenance": public_provenance,
        "review_required": True,
        "review_items": ["fixture_privacy", "source_and_contract", "effects_and_dependencies"]
        + (["absolute_fixture_paths"] if absolute_fixtures else []),
        "tests_require_opt_in": gated,
        "validation_status": "commands_prepared_not_executed",
        "validation_commands": commands,
        "publication": "manual_pr_maintainer_merge",
    }
    files[metadata_path] = (json.dumps(metadata, ensure_ascii=False, indent=2) + "\n").encode()
    return {resolved_root / path: value for path, value in files.items()}
