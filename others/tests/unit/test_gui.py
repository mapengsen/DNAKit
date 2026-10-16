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
