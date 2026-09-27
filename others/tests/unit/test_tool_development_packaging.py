"""Build and install a real wheel containing activated tools and optimized guidance.

Set DNAKIT_BUILD_PYTHON to an interpreter with setuptools>=77, wheel and build.
The test uses the current interpreter's runtime dependencies, copies real package
sources into a temporary checkout, and does not mutate or install into the repo.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import zipfile
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

import pytest

_NAME = "wheel_probe_length"
_FULL_NAME = f"dnakit.generated.{_NAME}"
_SKILL = "accessibility-aggregation"
_DESCRIPTION = "Count characters in a supplied DNA sequence."
_PARAMETER = "DNA sequence text whose character count is returned."
_SKILL_PREFIX = "Packaging verified workflow. "
_SOURCE = (
    '"""A deterministic generated-function packaging fixture."""\n\n'
    "def wheel_probe_length(sequence: str, *, verbose: bool = True) -> int:\n"
    '    """Return the number of characters in sequence."""\n'
    "    return len(sequence)\n"
)

_PREPARE = r"""
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
sys.path.insert(0, str(root / "src"))
from dnakit.tools import default_tool_registry
from dnakit.tools.contributions import contribution_files
from dnakit.tools.descriptions import append_description_version, description_contract
from dnakit.tools.development_support import DevelopmentJob
from dnakit.tools.optimization import _skill_target

tool_name, skill_id, description, parameter, skill_prefix = sys.argv[2:]
registry = default_tool_registry(verbose=False)
tool = registry.describe(tool_name)
store = append_description_version(
    {"schema_version": 1, "native": {}, "mcp": {}},
    kind="native", name=tool_name, server_id=None,
    contract=description_contract(tool),
    changes={"description": description, "parameters": {"sequence": parameter}},
    change_id="packaging-description-v1",
)
metadata = root / "src/dnakit/tools/_data/descriptions.json"
metadata.parent.mkdir(parents=True, exist_ok=True)
metadata.write_text(json.dumps(store, ensure_ascii=False), encoding="utf-8")
skill = root / "skills" / skill_id / "SKILL.md"
text = skill.read_text(encoding="utf-8")
candidate = text.replace("description: ", "description: " + skill_prefix, 1)
assert candidate != text
job = DevelopmentJob(root, "optimize_tool", skill_id,
                     {"reason": "Verify packaging of revised guidance"}, False)
_, _, writes = _skill_target(root, {"name": skill_id}, {"text": candidate}, job)
job.activate(writes)
assert skill.read_text(encoding="utf-8") == candidate
generated = root / "src/dnakit/generated"
entry = json.loads((generated / "manifest.json").read_text())["tools"][tool_name.rsplit(".", 1)[-1]]
materials = contribution_files(
    root=root, tool_name=tool_name,
    cases=[{"arguments": {"sequence": "ACGTAC"}, "expected": 6}],
    source_paths=[generated / entry["module"], generated / "manifest.json"],
    reason="Verify packaging of reviewable tool contributions.",
    change_id="generation-" + "a" * 24, effects=["local"], dependencies=[],
)
for path, content in materials.items():
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
"""

_VERIFY = r"""
import hashlib
import json
import sys
from pathlib import Path

target = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(target))
tool_name, skill_id, description, parameter, skill_prefix = sys.argv[2:]
import dnakit
from dnakit import generated
from dnakit.skills import get_skill, read_skill_file, search_skills
from dnakit.tools import default_tool_registry
import dnakit.tools.server as server_module
from dnakit.tools.development_support import _source_root
from dnakit.exceptions import ConfigurationError

assert Path(dnakit.__file__).resolve().is_relative_to(target)
assert Path(generated.__file__).resolve().is_relative_to(target)
assert generated.wheel_probe_length("ACGTAC") == 6
registry = default_tool_registry(verbose=False)
assert registry.execute(tool_name, {"sequence": "ACGTAC"}) == 6
assert any(item.name == tool_name for item in registry.search("wheel_probe_length"))
manifest = registry.describe(tool_name)
assert manifest["description"] == description
assert manifest["input_schema"]["properties"]["sequence"]["description"] == parameter

# Exercise DNAKit's actual MCP generic handlers without requiring FastMCP's
# optional transport dependency or opening any network connection.
class Registrar:
    def __init__(self, **settings):
        self.tools = {}
    def tool(self, *, name, description, annotations):
        def register(function):
            self.tools[name] = function
            return function
        return register

server_module._fastmcp_factory = lambda: Registrar
server = server_module.create_server(registry, verbose=False)
mcp_description = server.tools["describe_dnakit_tool"](tool_name)
mcp_call = server.tools["call_dnakit_tool"](tool_name, {"sequence": "ACGTAC"})
assert mcp_description == {"ok": True, "tool": manifest}
assert mcp_call["ok"] and mcp_call["result"] == 6
entry = get_skill(skill_id)
assert entry["license"] and entry["source_id"] and entry["upstream_url"]
assert entry["description"].startswith(skill_prefix)
assert any(item["id"] == skill_id for item in search_skills(skill_id)["skills"])
resource = read_skill_file(skill_id)
assert "description: " + skill_prefix in resource["content"]
mcp_resource = server.tools["read_dnakit_skill"](skill_id)
assert mcp_resource["ok"] and mcp_resource["content"] == resource["content"]
packaged_skill = target / "dnakit/skills/_data" / skill_id / "SKILL.md"
record = next(item for item in entry["files"] if item["path"] == "SKILL.md")
assert hashlib.sha256(packaged_skill.read_bytes()).hexdigest() == record["sha256"]
try:
    _source_root()
except ConfigurationError as exc:
    assert exc.code == "AGENT_DEVELOPMENT_SOURCE_REQUIRED"
else:
    raise AssertionError("Installed wheel must not be treated as a writable source checkout")
print(json.dumps({
    "imported_from": str(Path(dnakit.__file__).resolve()),
    "tool": tool_name, "python_result": 6, "registry_result": 6, "mcp_result": 6,
    "description": manifest["description"], "skill": skill_id,
    "skill_sha256": record["sha256"], "license": entry["license"],
    "source_id": entry["source_id"], "installed_development_rejected": True,
}, sort_keys=True))
"""


def _run(command: list[str], *, cwd: Path, log: Path, timeout: int = 300) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    log.write_text(completed.stdout, encoding="utf-8")
    assert completed.returncode == 0, f"See {log}\n{completed.stdout[-6000:]}"
    return completed.stdout


def _build_python() -> str:
    executable = os.environ.get("DNAKIT_BUILD_PYTHON", sys.executable)
    try:
        probe = subprocess.run(
            [
                executable,
                "-I",
                "-c",
                "import build, wheel; from importlib.metadata import version; "
                "from packaging.version import Version; "
                "assert Version(version('setuptools')) >= Version('77')",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
    except OSError:
        pytest.skip("DNAKIT_BUILD_PYTHON must point to an available build interpreter")
    if probe.returncode:
        pytest.skip("Wheel test requires setuptools>=77, wheel and build in DNAKIT_BUILD_PYTHON")
    return executable


def test_real_wheel_preserves_generated_functions_descriptions_and_skill_hashes(
    tmp_path: Path,
) -> None:
    builder = _build_python()
    root = Path(__file__).resolve().parents[3]
    project = tmp_path / "checkout"
    project.mkdir()
    print("[wheel validation] Copy real source and one cataloged skill", flush=True)
    shutil.copytree(
        root / "src",
        project / "src",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info"),
    )
    for name in (
        "pyproject.toml",
        "setup.py",
        "MANIFEST.in",
        "README.md",
        "README_CN.md",
        "CHANGELOG.md",
        "LICENSE",
        "DISCLAIMER.md",
        "THIRD_PARTY_NOTICES.md",
        "CONTRIBUTING.md",
        ".github/PULL_REQUEST_TEMPLATE.md",
    ):
        if (root / name).is_file():
            (project / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, project / name)
    shutil.copytree(root / "licenses", project / "licenses")
    catalog = cast(dict[str, Any], json.loads((root / "registry/catalog.json").read_text()))
    entry = deepcopy(next(item for item in catalog["skills"] if item["id"] == _SKILL))
    catalog["skills"] = [entry]
    catalog["sources"] = [item for item in catalog["sources"] if item["id"] == entry["source_id"]]
    catalog["categories"] = {entry["category"]: catalog["categories"][entry["category"]]}
    for resource in entry["files"]:
        source = root / "skills" / _SKILL / resource["path"]
        destination = project / "skills" / _SKILL / resource["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    (project / "registry").mkdir()
    (project / "registry/catalog.json").write_text(json.dumps(catalog), encoding="utf-8")

    digest = hashlib.sha256(_SOURCE.encode()).hexdigest()
    module = f"{_NAME}_{digest[:16]}.py"
    package = project / "src/dnakit/generated"
    (package / module).write_text(_SOURCE, encoding="utf-8")
    (package / "manifest.json").write_text(
        json.dumps(
            {
                "version": 1,
                "tools": {
                    _NAME: {
                        "module": module,
                        "sha256": digest,
                        "effect": "local",
                        "description": "Count DNA sequence characters.",
                        "contract": {"effects": ["local"]},
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    arguments = [_FULL_NAME, _SKILL, _DESCRIPTION, _PARAMETER, _SKILL_PREFIX]
    _run(
        [sys.executable, "-I", "-c", _PREPARE, str(project), *arguments],
        cwd=tmp_path,
        log=tmp_path / "prepare.log",
    )
    prepared_catalog = json.loads((project / "registry/catalog.json").read_text())
    prepared_entry = prepared_catalog["skills"][0]
    skill_record = next(item for item in prepared_entry["files"] if item["path"] == "SKILL.md")
    old_record = next(item for item in entry["files"] if item["path"] == "SKILL.md")
    assert skill_record["sha256"] != old_record["sha256"]
    assert prepared_entry["upstream_url"] == entry["upstream_url"]
    assert prepared_entry["license"] == entry["license"]

    print("[wheel validation] Build with the real setup.py and pyproject.toml", flush=True)
    dist = tmp_path / "dist"
    _run(
        [
            builder,
            "-I",
            "-m",
            "build",
            "--wheel",
            "--sdist",
            "--no-isolation",
            "--outdir",
            str(dist),
            str(project),
        ],
        cwd=tmp_path,
        log=tmp_path / "build.log",
    )
    wheels = list(dist.glob("dnakit-*.whl"))
    sdists = list(dist.glob("dnakit-*.tar.gz"))
    assert len(sdists) == 1
    with tarfile.open(sdists[0]) as source_archive:
        members = {name.partition("/")[2] for name in source_archive.getnames()}
        assert "CONTRIBUTING.md" in members
        assert ".github/PULL_REQUEST_TEMPLATE.md" in members
        materials = [
            path.relative_to(project).as_posix()
            for folder in (
                "others/tests/generated",
                "docs/llm_tools/generated",
                "registry/contributions",
            )
            for path in (project / folder).iterdir()
        ]
        assert len(materials) == 5 and all(path in members for path in materials)
        assert not any(".dnakit/" in path for path in members)
    assert len(wheels) == 1
    wheel = wheels[0]
    with zipfile.ZipFile(wheel) as archive:
        assert archive.read(f"dnakit/generated/{module}") == _SOURCE.encode()
        assert (
            json.loads(archive.read("dnakit/generated/manifest.json"))["tools"][_NAME]["sha256"]
            == digest
        )
        descriptions = json.loads(archive.read("dnakit/tools/_data/descriptions.json"))
        assert descriptions["native"][_FULL_NAME]["active"] == "packaging-description-v1"
        assert (
            hashlib.sha256(archive.read(f"dnakit/skills/_data/{_SKILL}/SKILL.md")).hexdigest()
            == skill_record["sha256"]
        )
        bundled = json.loads(archive.read("dnakit/skills/_data/catalog.json"))["skills"][0]
        assert bundled["files"] == prepared_entry["files"]
        assert bundled["license"] == entry["license"]

    print(
        "[wheel validation] Pip install into an isolated target and exercise installed APIs",
        flush=True,
    )
    target = tmp_path / "installed"
    _run(
        [
            builder,
            "-I",
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--no-index",
            "--target",
            str(target),
            str(wheel),
        ],
        cwd=tmp_path,
        log=tmp_path / "install.log",
    )
    output = _run(
        [sys.executable, "-I", "-c", _VERIFY, str(target), *arguments],
        cwd=tmp_path,
        log=tmp_path / "installed-validation.log",
    )
    result = json.loads(output.strip().splitlines()[-1])
    assert result["registry_result"] == result["mcp_result"] == 6
    assert result["skill_sha256"] == skill_record["sha256"]
    report = {"wheel": str(wheel), "wheel_bytes": wheel.stat().st_size, **result}
    (tmp_path / "packaging-validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, sort_keys=True), flush=True)
