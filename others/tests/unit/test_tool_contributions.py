"""Contribution bundles preserve executable cases without publishing private state."""

from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, cast

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools.contributions import contribution_files

_CHANGE = "generation-" + "a" * 24
_NAME = "dnakit.descriptors.canonical_kmer"


def _bundle(root: Path, **updates: Any) -> dict[Path, bytes]:
    arguments: dict[str, Any] = {
        "root": root,
        "tool_name": _NAME,
        "cases": [
            {"arguments": {"value": "TTTT"}, "expected": "AAAA"},
            {
                "arguments": {"value": "N"},
                "expected_error": {"type": "ConfigurationError", "code": "CONFIGURATION_ERROR"},
            },
        ],
        "source_paths": [root / "src/dnakit/descriptors/kmer.py"],
        "reason": "Preserve the canonical strand representation and reject ambiguous bases.",
        "change_id": _CHANGE,
        "effects": ["local"],
    }
    arguments.update(updates)
    return contribution_files(**arguments)


def _metadata(bundle: dict[Path, bytes]) -> dict[str, Any]:
    return cast(
        dict[str, Any],
        json.loads(next(value for path, value in bundle.items() if path.suffix == ".json")),
    )


def _run_generated(bundle: dict[Path, bytes], root: Path) -> subprocess.CompletedProcess[str]:
    test_path = next(path for path in bundle if path.suffix == ".py")
    test_path.parent.mkdir(parents=True, exist_ok=True)
    test_path.write_bytes(bundle[test_path])
    environment = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[3] / "src")}
    environment.pop("DNAKIT_RUN_GENERATED_EFFECT_TESTS", None)
    return subprocess.run(
        [sys.executable, "-m", "pytest", str(test_path), "-q", "-o", "addopts="],
        cwd=root,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=90,
    )


def test_bundle_is_pure_deterministic_and_contains_only_portable_metadata(tmp_path: Path) -> None:
    bundle = _bundle(
        tmp_path,
        provenance={
            "upstream_url": "https://example.org/author/project",
            "upstream_commit": "123abc",
            "license": "MIT",
            "dependencies": ["numpy"],
            "audit_dir": "/private/runtime/audit",
            "api_key": "never-copy-this",
        },
    )
    assert not any(tmp_path.iterdir())
    assert len(bundle) == 5
    metadata = _metadata(bundle)
    assert metadata["schema_version"] == 1 and metadata["change_id"] == _CHANGE
    assert metadata["source_paths"] == ["src/dnakit/descriptors/kmer.py"]
    assert metadata["dependencies"] == ["numpy"]
    assert metadata["review_required"] and not metadata["tests_require_opt_in"]
    assert metadata["case_count"] == 2
    assert metadata["validation_status"] == "commands_prepared_not_executed"
    assert set(metadata["files"]) == {path.relative_to(tmp_path).as_posix() for path in bundle}
    all_content = b"\n".join(bundle.values()).decode()
    assert str(tmp_path) not in all_content
    assert "never-copy-this" not in all_content and "/private/runtime/audit" not in all_content
    assert "call_dnakit_tool" in all_content and "maintainers review and merge" in all_content
    assert "123abc" in all_content and "MIT" in all_content
    assert _bundle(tmp_path) == _bundle(tmp_path)
    tree = ast.parse(next(value.decode() for path, value in bundle.items() if path.suffix == ".py"))
    assert not any(
        isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "exec"
        for node in ast.walk(tree)
    )


def test_generated_pytest_really_checks_python_agent_and_expected_errors(tmp_path: Path) -> None:
    result = _run_generated(_bundle(tmp_path), tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "4 passed" in result.stdout


def test_generated_pytest_executes_only_the_selected_custom_target(tmp_path: Path) -> None:
    module_directory = tmp_path / "src/dnakit"
    module_directory.mkdir(parents=True)
    source = module_directory / "contribution_custom.py"
    source.write_text(
        '"""Custom public function outside the default catalog."""\n\n'
        "def doubled(value: int) -> int:\n"
        '    """Double an integer while rejecting negative values."""\n'
        '    if value < 0:\n        raise ValueError("Expected a nonnegative integer")\n'
        "    return value * 2\n\n"
        '__all__ = ["doubled"]\n'
    )
    name = "dnakit.contribution_custom.doubled"
    bundle = _bundle(
        tmp_path,
        tool_name=name,
        source_paths=[source],
        cases=[
            {"arguments": {"value": 3}, "expected": 6},
            {"arguments": {"value": -1}, "expected_error": {"type": "ValueError"}},
        ],
    )
    test_directory = next(path.parent for path in bundle if path.suffix == ".py")
    test_directory.mkdir(parents=True)
    (test_directory / "conftest.py").write_text(
        "from pathlib import Path\n"
        "import dnakit\n"
        "from dnakit.tools.registry import default_tool_registry\n"
        'dnakit.__path__.append(str(Path(__file__).resolve().parents[3] / "src/dnakit"))\n'
        "def pytest_sessionfinish(session, exitstatus):\n"
        f"    assert all(t.name != {name!r} for t in default_tool_registry(verbose=False).tools)\n"
    )
    result = _run_generated(bundle, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "4 passed" in result.stdout
    english = next(value.decode() for path, value in bundle.items() if path.name.endswith(".en.md"))
    pr_body = bundle[tmp_path / f"registry/contributions/{_CHANGE}.md"].decode()
    assert "custom-registered target" in english and "custom-registered target" in pr_body
    assert "keyword-only with defaults" in pr_body


def test_generated_pytest_rejects_wrong_expected_result(tmp_path: Path) -> None:
    result = _run_generated(
        _bundle(tmp_path, cases=[{"arguments": {"value": "TTTT"}, "expected": "TTTT"}]),
        tmp_path,
    )
    assert result.returncode == 1
    assert "2 failed" in result.stdout


def test_template_like_data_and_unicode_are_preserved_as_data(tmp_path: Path) -> None:
    cases = [{"arguments": {"value": "__GATED__\n中文 'quoted'"}, "expected": "__CASES__"}]
    bundle = _bundle(tmp_path, cases=cases)
    tree = ast.parse(next(value for path, value in bundle.items() if path.suffix == ".py"))
    call = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and isinstance(node.targets[0], ast.Name)
        and node.targets[0].id == "_CASES"
    )
    assert isinstance(call, ast.Call)
    assert json.loads(ast.literal_eval(call.args[0])) == cases


@pytest.mark.parametrize("effects", [["network"], None])
def test_effectful_or_unknown_tests_skip_before_importing_the_tool(
    tmp_path: Path, effects: list[str] | None
) -> None:
    bundle = _bundle(tmp_path, tool_name="dnakit.uninstalled_backend.missing_tool", effects=effects)
    result = _run_generated(bundle, tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "4 skipped" in result.stdout
    assert _metadata(bundle)["tests_require_opt_in"]
    assert b"DNAKIT_RUN_GENERATED_EFFECT_TESTS=1" in next(
        content for path, content in bundle.items() if path.name.endswith(".en.md")
    )


def test_absolute_fixtures_are_flagged_and_do_not_leak_into_metadata(tmp_path: Path) -> None:
    fixture = str(tmp_path / "private-input.fasta")
    bundle = _bundle(tmp_path, cases=[{"arguments": {"input_path": fixture}, "expected": 1}])
    metadata = _metadata(bundle)
    assert metadata["tests_require_opt_in"]
    assert "absolute_fixture_paths" in metadata["review_items"]
    assert fixture not in json.dumps(metadata)
    assert "cases" not in metadata


@pytest.mark.parametrize(
    "arguments",
    [
        {"api_key": "credential"},
        {"nested": [{"access_token": "credential"}]},
        {"url": "https://user:password@example.org/data"},
        {"url": "https://example.org/data?token=credential"},
        {"key": "-----BEGIN RSA PRIVATE KEY-----\nprivate\n-----END RSA PRIVATE KEY-----"},
    ],
)
def test_recognizable_secrets_are_rejected_without_writes(
    tmp_path: Path, arguments: dict[str, Any]
) -> None:
    with pytest.raises(ConfigurationError, match="TOOL_CONTRIBUTION_INVALID"):
        _bundle(tmp_path, cases=[{"arguments": arguments, "expected": None}])
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize(
    "updates",
    [
        {"tool_name": "dnakit.foo.bar; import os"},
        {"change_id": "../../outside"},
        {"source_paths": [Path("../outside.py")]},
        {"source_paths": [Path(".dnakit/development/private.json")]},
        {"cases": [{"arguments": {}, "expected": float("nan")}]},
        {"cases": [{"arguments": {}, "expected": 1, "expected_error": {"type": "ValueError"}}]},
        {"cases": [{"arguments": {}, "expected_error": {}}]},
        {"effects": ["local", "network"]},
        {"provenance": {"upstream_url": "file:///private/source"}},
    ],
)
def test_invalid_or_nonportable_contracts_are_rejected(
    tmp_path: Path, updates: dict[str, Any]
) -> None:
    with pytest.raises(ConfigurationError, match="TOOL_CONTRIBUTION_INVALID"):
        _bundle(tmp_path, **updates)


def test_source_symlink_cannot_escape_checkout(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    (root / "external").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ConfigurationError, match="inside the DNAKit checkout"):
        _bundle(root, source_paths=[root / "external" / "source.py"])
