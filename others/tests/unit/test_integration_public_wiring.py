"""Public integration regressions: manifest exports, aliases and MCP execution.

Version note: covers all 46 packs without installing or running optional models.
"""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import json
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.tools import default_tool_registry

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "paper/write_paper/figure/fig2/task_coverage"
VERIFIER = BASE / "integration/verify_public_integration.py"
if not VERIFIER.is_file():
    # paper/write_paper/ is a local, git-ignored folder, so a fresh checkout has no manifests.
    pytest.skip(
        "needs the local paper/write_paper task-coverage manifests", allow_module_level=True
    )
_spec = importlib.util.spec_from_file_location("public_integration_verifier", VERIFIER)
verifier = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(verifier)


@pytest.mark.parametrize("number", range(46))
def test_every_manifest_export_reaches_public_python_and_agent_catalog(number):
    manifest = json.loads(
        (BASE / f"coverage_gap/manifests/integration_p{number:02d}.json").read_text()
    )
    registry = default_tool_registry(verbose=False)
    for module, name in verifier.declarations(manifest):
        function = getattr(importlib.import_module(module), name)
        area = module.rsplit(".", 1)[0]
        public = verifier.ALIASES.get(f"{module}.{name}", name)
        exported = importlib.import_module(area)
        assert public in exported.__all__
        assert inspect.unwrap(getattr(exported, public)) is inspect.unwrap(function)
        tool = registry.resolve(f"{area}.{public}")
        assert inspect.unwrap(tool.function) is inspect.unwrap(function)
        assert tool.agent_compatible, tool.incompatibilities


def test_new_alias_preserves_existing_trim_api():
    from dnakit.sequence import biopython_trim_sequence, trim_sequence

    assert trim_sequence("AACGTT", left=1, right=2, verbose=False) == "ACG"
    pytest.importorskip("Bio")
    assert biopython_trim_sequence("AACGTT", pattern="AT", verbose=False)["trimmed"] == "CG"


@pytest.mark.parametrize(
    ("name", "effect"),
    [
        ("dnakit.alignment.mmseqs_taxonomy", "external_process"),
        ("dnakit.variants.deepvariant_run", "external_process"),
        ("dnakit.epigenome.modkit_pileup", "external_process"),
        ("dnakit.sequence.write_sequence_file", "filesystem"),
        ("dnakit.predictions.forward_evo2_sequences", "model"),
    ],
)
def test_new_side_effecting_facades_keep_agent_authorization(name, effect):
    registry = default_tool_registry(verbose=False)
    assert registry.describe(name)["effect"] == effect
    with pytest.raises(ConfigurationError) as caught:
        registry.execute(name, {})
    assert caught.value.code == "AGENT_TOOL_CONFIRMATION_REQUIRED"


def test_compact_mcp_executes_new_sequence_tool(monkeypatch):
    pytest.importorskip("Bio")
    from dnakit.tools import server as server_module

    class Server:
        def __init__(self, **kwargs):
            self.functions = {}

        def tool(self, **kwargs):
            def register(function):
                self.functions[kwargs["name"]] = function
                return function

            return register

    monkeypatch.setattr(server_module, "_fastmcp_factory", lambda: Server)
    server = server_module.create_server(verbose=False)
    description = server.functions["describe_dnakit_tool"](
        "dnakit.sequence.biopython_trim_sequence"
    )
    assert description["ok"]
    response = server.functions["call_dnakit_tool"](
        "dnakit.sequence.biopython_trim_sequence",
        {"sequence": "AACGTT", "pattern": "AT", "verbose": False},
    )
    assert response["ok"]
    assert response["result"]["trimmed"] == "CG"


def test_crispresso_missing_dependency_has_working_install_hint(monkeypatch):
    from dnakit.design_crispr import crispresso_p13_cli_functions as module
    from dnakit.exceptions import BackendUnavailableError

    def missing(*args):
        raise BackendUnavailableError("missing backend")

    monkeypatch.setattr(module._common, "require", missing)
    with pytest.raises(BackendUnavailableError) as caught:
        module._require_crispresso("CRISPResso2.CRISPRessoShared")
    assert "conda install" in caught.value.hint
    assert "dnakit[crispresso2]" not in caught.value.hint
