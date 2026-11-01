"""Tests for dnakit.molecular_lab.dnacauldron_gap_functions (simulation PDF report, leveled graph).

Offline: outputs go to tmp_path. The graph test draws with dnacauldron's own layout and checks the
node positions against the layout rule. The PDF writer needs dnacauldron's reports extra
(pdf_reports, which needs WeasyPrint), which is not installed here: its guard is tested, and the
happy path is skipped with that reason.
"""

from __future__ import annotations

import importlib.util
import inspect
import socket
import sys
from pathlib import Path

import pytest

pytest.importorskip("dnacauldron")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab import dnacauldron_gap_functions as gap

@pytest.fixture(autouse=True)
def _close_figures_left_by_other_tests():
    """Other tests in the same process may leave figures open; these tests check that the function closes its own."""
    plt = pytest.importorskip("matplotlib.pyplot")
    plt.close("all")
    yield


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
SIMULATION = [
    {"Outcome": "Valid", "Number of assemblies": 12},
    {"Outcome": "Cancelled", "Number of assemblies": 0},
    {"Outcome": "Errored", "Number of assemblies": 1},
]
PDF_REPORTS_MISSING = importlib.util.find_spec("pdf_reports") is None


class _NoConnectSocket(socket.socket):
    def connect(self, *args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    def connect_ex(self, *args: object, **kwargs: object) -> int:
        raise AssertionError("network access is not allowed in this test")


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(socket, "socket", _NoConnectSocket)


def test_public_functions_take_keyword_only_verbose() -> None:
    assert gap.__all__
    for name in gap.__all__:
        parameter = inspect.signature(getattr(gap, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


# Leveled assembly graph.


def test_leveled_graph_places_nodes_by_the_layout_rule(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt

    target = tmp_path / "graphs" / "plan.png"
    result = gap.dnacauldron_plot_leveled_graph(
        str(target), [["A", "B"], ["C"]], [["A", "C"], ["B", "C"]], verbose=False
    )
    assert result["path"] == str(target.resolve())
    assert target.read_bytes()[:8] == PNG_SIGNATURE
    assert plt.get_fignums() == []
    from matplotlib.image import imread

    pixels = imread(str(target))[..., :3]
    assert (pixels.min(axis=-1) < 0.2).mean() > 0.001  # black box outlines, edges and text
    data = result["data"]
    assert (data["level_count"], data["node_count"], data["edge_count"]) == (2, 3, 2)
    # Level 0 sits at x = 1/3 and level 1 at x = 2/3; level-0 nodes split the height at 1/3, 2/3.
    positions = data["positions"]
    assert positions["A"] == pytest.approx([1 / 3, 1 / 3])
    assert positions["B"] == pytest.approx([1 / 3, 2 / 3])
    assert positions["C"] == pytest.approx([2 / 3, 0.5])


def test_leveled_graph_guards(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"more than one place"):
        gap.dnacauldron_plot_leveled_graph(
            str(tmp_path / "a.png"), [["A"], ["A"]], [], verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"not a node"):
        gap.dnacauldron_plot_leveled_graph(
            str(tmp_path / "b.png"), [["A"]], [["A", "Z"]], verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"pair"):
        gap.dnacauldron_plot_leveled_graph(str(tmp_path / "c.png"), [["A"]], [["A"]], verbose=False)
    with pytest.raises(ConfigurationError, match=r"\.png"):
        gap.dnacauldron_plot_leveled_graph(str(tmp_path / "d.svg"), [["A"]], [], verbose=False)
    existing = tmp_path / "e.png"
    existing.write_bytes(PNG_SIGNATURE)
    with pytest.raises(ConfigurationError, match=r"already exists"):
        gap.dnacauldron_plot_leveled_graph(str(existing), [["A"]], [], verbose=False)


# Simulation PDF report.


def test_pdf_report_validates_the_simulation_rows(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match=r"exactly the keys"):
        gap.dnacauldron_write_simulation_pdf_report(
            str(tmp_path / "r.pdf"), [{"Outcome": "Valid"}], verbose=False
        )
    with pytest.raises(ConfigurationError, match=r"non-negative"):
        gap.dnacauldron_write_simulation_pdf_report(
            str(tmp_path / "r.pdf"),
            [{"Outcome": "Valid", "Number of assemblies": -1}],
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match=r"\.pdf"):
        gap.dnacauldron_write_simulation_pdf_report(
            str(tmp_path / "r.png"), SIMULATION, verbose=False
        )


def test_pdf_report_guard_names_the_reports_extra(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setitem(sys.modules, "pdf_reports", None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        gap.dnacauldron_write_simulation_pdf_report(
            str(tmp_path / "r.pdf"), SIMULATION, verbose=False
        )
    assert "dnacauldron[reports]" in excinfo.value.hint
    assert not (tmp_path / "r.pdf").exists()


@pytest.mark.skipif(
    PDF_REPORTS_MISSING,
    reason="pdf_reports (dnacauldron's reports extra, needs WeasyPrint) is not installed here",
)
def test_pdf_report_happy_path(
    tmp_path: Path,
) -> None:  # pragma: no cover - needs pdf_reports and WeasyPrint
    target = tmp_path / "report.pdf"
    result = gap.dnacauldron_write_simulation_pdf_report(str(target), SIMULATION, verbose=False)
    assert target.read_bytes()[:4] == b"%PDF"
    assert result["data"]["total_assemblies"] == 13
