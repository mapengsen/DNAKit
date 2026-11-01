"""The graphical interface: local and public-demo settings, handlers, and the CLI entry point."""

from __future__ import annotations

import json
import re
import sys
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from dnakit.cli.app import app
from dnakit.exceptions import BackendUnavailableError
from dnakit.gui import GuiSettings

runner = CliRunner()


def test_local_and_public_settings_differ_only_in_limits(tmp_path: Path) -> None:
    local = GuiSettings.local(tmp_path / "local")
    public = GuiSettings.public_demo(tmp_path / "public")

    assert local.mode == "local" and not local.is_public
    assert local.allow_side_effects and not local.confine_paths
    assert local.call_timeout_s is None and local.max_upload_mb is None
    assert public.is_public and not public.allow_side_effects and public.confine_paths
    assert public.call_timeout_s == 60 and public.max_upload_mb == 20
    assert local.workdir == (tmp_path / "local").resolve()


def test_cli_gui_help_does_not_need_gui_dependencies() -> None:
    result = runner.invoke(app, ["gui", "--help"])

    assert result.exit_code == 0, result.output
    # Typer forces colored help on GitHub Actions; compare the plain text.
    plain = re.sub(r"\x1b\[[0-9;]*m", "", result.output)
    for option in ("--browser", "--workdir", "--no-open", "--public-demo"):
        assert option in plain


def test_missing_gradio_reports_the_gui_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    from dnakit.gui import build_app

    monkeypatch.setitem(sys.modules, "gradio", None)
    with pytest.raises(BackendUnavailableError) as info:
        build_app(GuiSettings.local())
    assert info.value.code == "GUI_BACKEND_UNAVAILABLE"
    assert info.value.hint is not None and "dnakit[gui]" in info.value.hint


@pytest.fixture
def make_backend(tmp_path: Path) -> Iterator[Callable[[str], Any]]:
    pytest.importorskip("gradio")
    from dnakit.gui._backend import GuiBackend

    backends: list[Any] = []

    def make(mode: str) -> Any:
        folder = tmp_path / mode
        settings = (
            GuiSettings.public_demo(folder) if mode == "public" else GuiSettings.local(folder)
        )
        backend = GuiBackend(settings)
        backends.append(backend)
        return backend

    yield make
    for backend in backends:
        backend.worker.stop()


def test_public_demo_confines_paths_and_blocks_credentials(
    make_backend: Callable[[str], Any],
) -> None:
    import gradio as gr

    backend = make_backend("public")
    _, folder = backend.session(None)

    backend.check_arguments({"path": "input.fasta", "value": "ACGT"}, folder)
    for arguments in ({"path": "/etc/passwd"}, {"path": "../other/x.fa"}, {"url": "http://x"}):
        with pytest.raises(gr.Error):
            backend.check_arguments(arguments, folder)
    assert not backend.runnable("filesystem") and backend.runnable("network")


def test_local_mode_allows_paths_and_confirmed_side_effects(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("local")
    _, folder = backend.session(None)

    assert folder == backend.settings.workdir
    backend.check_arguments({"path": "/etc/hosts", "url": "https://example.org"}, folder)
    assert backend.runnable("filesystem") and backend.runnable("model")


def test_sequence_tab_rejects_sequences_beyond_the_command_line_limit(
    make_backend: Callable[[str], Any],
) -> None:
    import gradio as gr

    backend = make_backend("local")

    assert backend.clean_sequence(">id\nACGT\nacgt\n") == "ACGTacgt"
    with pytest.raises(gr.Error):
        backend.clean_sequence("A" * (backend.settings.max_sequence_length + 1))


def test_sequence_tab_reverse_complements_iupac_text_and_saves_fasta(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("public")

    state, head, reverse, path, code = backend.seq_reverse_complement(">id\nacgRN\nTTA", None)

    assert reverse == "TAANYCGT"
    assert "8 nt" in head
    assert Path(path).read_text(encoding="utf-8") == ">reverse_complement\nTAANYCGT\n"
    assert "reverse_complement(sequence)" in code and "Command line" not in code
    _, head, reverse, path, _ = backend.seq_reverse_complement("../etc/hosts", state)
    assert head.startswith("⚠️") and reverse is None and path is None


def test_dataset_deduplication_saves_results_in_the_local_folder(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("local")

    _, head, groups, files, code = backend.ds_deduplicate(None, True, "reverse_complement", None)

    assert "Removed **2** duplicates" in head
    assert {row[0] for row in groups} == {"seq03", "seq07"}
    assert Path(files[0]).parent == backend.settings.workdir and Path(files[0]).is_file()
    assert "dnakit deduplicate example.fasta" in code


def test_function_explorer_runs_a_tool_and_shows_equivalent_calls(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("public")

    _, status, result, _, _, code = backend.fx_run(
        "dnakit.ops.reverse_complement", json.dumps({"value": "AACG"}), False, None
    )

    assert status.startswith("✅"), status
    assert "CGTT" in json.dumps(result)
    assert "from dnakit.ops import reverse_complement" in code
    assert '"tool_name": "dnakit.ops.reverse_complement"' in code


def test_side_effect_tools_need_confirmation_locally_and_never_run_publicly(
    make_backend: Callable[[str], Any],
) -> None:
    import gradio as gr

    local, public = make_backend("local"), make_backend("public")
    name = next(t.name for t in local.registry.tools if t.effect == "filesystem")

    with pytest.raises(gr.Error, match="tick the box"):
        local.fx_run(name, "{}", False, None)
    with pytest.raises(gr.Error, match="not executed"):
        public.fx_run(name, "{}", True, None)


def test_example_arguments_follow_the_schema() -> None:
    pytest.importorskip("gradio")
    from dnakit.gui._backend import EXAMPLE_SEQUENCE, example_value

    assert example_value("value", {"type": "string"}) == EXAMPLE_SEQUENCE
    assert example_value("values", {"type": "array"}) == []
    assert example_value("task", {"enum": ["classification", "regression"]}) == "classification"
    assert example_value("k", {"type": "integer"}) == 3
    assert example_value("window", {"type": "integer", "default": 200}) == 200


def test_both_modes_build_the_same_page(tmp_path: Path) -> None:
    gr = pytest.importorskip("gradio")
    from dnakit.gui import build_app

    for settings in (GuiSettings.local(tmp_path / "a"), GuiSettings.public_demo(tmp_path / "b")):
        demo = build_app(settings)
        assert isinstance(demo, gr.Blocks)
        demo.close()


def test_page_header_reports_the_shipped_catalogue_counts(tmp_path: Path) -> None:
    pytest.importorskip("gradio")
    from dnakit.gui._app import _hero, catalog_counts
    from dnakit.gui._backend import GuiBackend

    counts = catalog_counts()
    assert counts is not None
    assert counts["operations_total"] == (
        counts["own_functions"] + counts["mcp_bridge"] + counts["agent_skills"]
    )
    header = _hero(GuiBackend(GuiSettings.public_demo(tmp_path)))
    assert f"{counts['operations_total']:,} DNA operations" in header
    assert f"{counts['mcp_tools']:,} tools" in header


def test_function_explorer_tab_uses_the_catalogue_function_count(tmp_path: Path) -> None:
    pytest.importorskip("gradio")
    from dnakit.gui import build_app
    from dnakit.gui._app import catalog_counts

    counts = catalog_counts()
    assert counts is not None
    demo = build_app(GuiSettings.public_demo(tmp_path))
    labels = [getattr(block, "label", None) for block in demo.blocks.values()]
    assert f"Function explorer ({counts['own_functions']:,})" in labels
    demo.close()


def test_mcp_tab_reports_the_catalogue_tool_count_and_searches_bridged_tools(
    tmp_path: Path,
) -> None:
    pytest.importorskip("gradio")
    from dnakit.gui import build_app
    from dnakit.gui._app import catalog_counts
    from dnakit.gui._backend import GuiBackend

    counts = catalog_counts()
    assert counts is not None
    demo = build_app(GuiSettings.public_demo(tmp_path))
    labels = [getattr(block, "label", None) for block in demo.blocks.values()]
    assert f"MCP tools ({counts['mcp_tools']:,})" in labels
    demo.close()
    backend = GuiBackend(GuiSettings.public_demo(tmp_path))
    head, rows = backend.mcp_search("gwas")
    assert rows and "matching third-party tools" in head
    assert backend.mcp_servers()


def test_function_explorer_links_functions_to_their_documentation(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("public")

    _, rows, _ = backend.fx_search("reverse complement", "")
    by_name = {row[0]: row[1] for row in rows}
    link = backend.docs.functions["dnakit.ops.reverse_complement"]
    assert link.startswith("https://mapengsen.github.io/DNAKit/en/api/features/")
    assert by_name["dnakit.ops.reverse_complement"] == f"[docs ↗]({link})"
    doc = backend.fx_describe("dnakit.ops.reverse_complement")[0]
    assert f"[documentation ↗]({link})" in doc


def test_shipped_documentation_links_match_the_docs_folder() -> None:
    import importlib.util

    root = Path(__file__).resolve().parents[3]
    script = root / "others/docs/build_gui_doc_links.py"
    if not (root / "docs").is_dir():
        pytest.skip("documentation sources are not available")
    spec = importlib.util.spec_from_file_location("build_gui_doc_links", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    shipped = json.loads(module.OUTPUT.read_text(encoding="utf-8"))
    rebuilt = module.build(module.registry_tools(), module.catalog_targets())
    assert shipped == rebuilt, "run python others/docs/build_gui_doc_links.py"
    pages = module.pages(skip=False)
    anchors: dict[str, set[str]] = {}
    for link in module.all_links(shipped).values():
        page, _, anchor = link.partition("#")
        assert page in pages, link
        if page not in anchors:
            anchors[page] = module.page_anchors(pages[page])
        assert not anchor or anchor in anchors[page], link


def test_mcp_tab_describes_a_clicked_third_party_tool(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("public")

    detail = backend.mcp_describe("gwas-mcp", "query_gwas_catalog")

    assert detail.startswith("### `query_gwas_catalog`")
    assert "https://github.com/" in detail and "/tree/" in detail
    assert "call_mcp_tool('gwas-mcp', 'query_gwas_catalog'" in detail
    assert backend.mcp_describe("gwas-mcp", "no_such_tool").startswith("⚠️")
    link = backend.docs.mcp_tools["gwas-mcp/query_gwas_catalog"]
    assert link.endswith("llm_tools/mcp_tools/gwas-mcp/#query_gwas_catalog")
    assert f"[documentation ↗]({link})" in detail
    _, rows = backend.mcp_search("query_gwas_catalog")
    row = next(r for r in rows if r[0] == "query_gwas_catalog" and r[2] == "gwas-mcp")
    assert row[1] == f"[docs ↗]({link})"


def test_skills_tab_links_skills_to_the_skills_catalogue(
    make_backend: Callable[[str], Any],
) -> None:
    backend = make_backend("public")

    _, rows, _ = backend.sk_search("variant annotation")
    skill, docs = rows[0][0], rows[0][1]
    link = backend.docs.skills[skill]
    assert "/llm_tools/skills_catalog/" in link and link.endswith(f"/#{skill}")
    assert docs == f"[docs ↗]({link})"
    assert f"[documentation ↗]({link})" in backend.sk_read(skill)
