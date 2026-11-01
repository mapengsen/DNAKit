"""Tests for dnakit.sequence.gimmemotifs_pending_functions (Fig. 2 coverage gap, batch G).

Offline by design: the sample, background, motif and report inputs are small files written into
tmp_path. De novo tools and the prediction tool are replaced by fakes through monkeypatch, so the
tests pass without MEME, Homer or any other program installed. The one test that uses the real tool
lookup is skipped when the tool is configured and installed on this machine. The tests import the
module directly, so they run before the area package is wired.
"""

from __future__ import annotations

import importlib
import inspect
import random
from pathlib import Path

import pytest

pytest.importorskip("gimmemotifs")
pytest.importorskip("pandas")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pending = importlib.import_module("dnakit.sequence.gimmemotifs_pending_functions")

PFM_TEXT = """>MA0001.1
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
>MA0002.1
0.1 0.1 0.7 0.1
0.1 0.7 0.1 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
0.4 0.2 0.2 0.2
0.2 0.4 0.2 0.2
"""

MAELSTROM_TABLE = (
    "\tHypergeom\tMean\tHypergeom % with motif\tcorr Hypergeom\n"
    "MA0001.1\t5.2\t3.1\t80.0\t0.9\n"
    "MA0002.1\t-1.2\t0.2\t20.0\t0.1\n"
)


def _seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records), encoding="utf-8")
    return path


@pytest.fixture()
def gimme_inputs(tmp_path: Path) -> dict[str, Path]:
    rng = random.Random(1)
    sample = []
    for index in range(40):
        seq = _seq(rng, 200)
        if index % 2 == 0:
            seq = seq[:50] + "TGACGTCA" + seq[58:]
        sample.append((f"s{index}", seq))
    background = [(f"b{index}", _seq(rng, 200)) for index in range(80)]
    pfm = tmp_path / "motifs.pfm"
    pfm.write_text(PFM_TEXT, encoding="utf-8")
    return {
        "sample": _write_fasta(tmp_path / "sample.fa", sample),
        "background": _write_fasta(tmp_path / "bg.fa", background),
        "pfm": pfm,
    }


class _FakeProgram:
    """Stands in for a configured and installed gimmemotifs tool program."""

    def __init__(self, pfm: Path | None = None, installed: bool = True) -> None:
        self.pfm = pfm
        self.installed = installed
        self.calls: list[tuple[str, dict]] = []

    def is_installed(self) -> bool:
        return self.installed

    def run(self, fastafile, params=None):
        self.calls.append((fastafile, dict(params or {})))
        from gimmemotifs.motif import read_motifs

        motifs = read_motifs(str(self.pfm), fmt="pfm") if self.pfm else []
        return motifs, "predicted two motifs\n", "a warning\n"


def _patch_tools(monkeypatch: pytest.MonkeyPatch, factory) -> None:
    tools = importlib.import_module("gimmemotifs.tools")
    monkeypatch.setattr(tools, "get_tool", factory)


# ---------------------------------------------------------------- gimme motifs


def test_gimme_motifs_known_motifs_with_explicit_background(tmp_path, gimme_inputs):
    out = tmp_path / "motifs_run"
    result = pending.gimmemotifs_gimme_motifs(
        str(gimme_inputs["sample"]),
        str(gimme_inputs["pfm"]),
        str(out),
        background_fasta=str(gimme_inputs["background"]),
        verbose=False,
    )
    assert result["status"] == "complete"
    assert result["operation"] == "gimme_motifs"
    assert result["backend"] == "gimmemotifs"
    assert result["denovo_tools"] == []
    assert result["motifs_scored"] == 2
    assert Path(result["outputs"]["roc_report"]).is_file()
    assert Path(result["outputs"]["gimme.motifs.html"]).is_file()
    assert (out / "motif_scan_results").is_dir()
    assert "combined.pfm" not in result["outputs"]


def test_gimme_motifs_needs_known_motifs_or_denovo_tools(tmp_path, gimme_inputs):
    with pytest.raises(ConfigurationError, match="known=True"):
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(tmp_path / "out"),
            known=False,
            verbose=False,
        )
    assert not (tmp_path / "out").exists()


def test_gimme_motifs_missing_denovo_tool_raises_before_work(tmp_path, gimme_inputs, monkeypatch):
    def unconfigured(name):
        raise ValueError(f"No configuration found for {name}")

    _patch_tools(monkeypatch, unconfigured)
    out = tmp_path / "out"
    with pytest.raises(BackendUnavailableError) as info:
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(out),
            denovo_tools=["MEME"],
            verbose=False,
        )
    assert "MEME" in str(info.value)
    assert info.value.hint and "gimmemotifs.cfg" in info.value.hint
    assert not out.exists()


def test_gimme_motifs_not_installed_tool_raises(tmp_path, gimme_inputs, monkeypatch):
    _patch_tools(monkeypatch, lambda name: _FakeProgram(installed=False))
    with pytest.raises(BackendUnavailableError, match="not installed"):
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(tmp_path / "out"),
            denovo_tools=["Homer"],
            verbose=False,
        )


def test_gimme_motifs_rejects_unknown_tool_and_analysis(tmp_path, gimme_inputs):
    with pytest.raises(ConfigurationError, match="Unknown gimmemotifs tool"):
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(tmp_path / "a"),
            denovo_tools=["NotATool"],
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="analysis must be one of"):
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(tmp_path / "b"),
            analysis="huge",
            verbose=False,
        )


def test_gimme_motifs_refuses_existing_output_directory(tmp_path, gimme_inputs):
    existing = tmp_path / "exists"
    existing.mkdir()
    with pytest.raises(ConfigurationError, match="already exists"):
        pending.gimmemotifs_gimme_motifs(
            str(gimme_inputs["sample"]),
            str(gimme_inputs["pfm"]),
            str(existing),
            verbose=False,
        )


# ------------------------------------------------------------ gimme prediction


def test_predict_motifs_writes_pfm_from_tool_output(tmp_path, gimme_inputs, monkeypatch):
    fake = _FakeProgram(pfm=gimme_inputs["pfm"])
    seen: list[str] = []

    def factory(name):
        seen.append(name)
        return fake

    _patch_tools(monkeypatch, factory)
    out = tmp_path / "predicted.pfm"
    result = pending.gimmemotifs_predict_motifs(
        "MEME",
        str(gimme_inputs["sample"]),
        str(out),
        params={"width": 8, "nmotifs": 2},
        verbose=False,
    )
    assert seen == ["MEME"]
    assert fake.calls[0][1] == {"width": 8, "nmotifs": 2}
    assert fake.calls[0][0].endswith("sample.fa")
    text = out.read_text(encoding="utf-8")
    assert text.count(">") == 2 and "MA0001.1" in text
    assert result["motifs"] == 2
    assert result["tool"] == "MEME"
    assert "predicted two motifs" in result["stdout_tail"]
    assert "a warning" in result["stderr_tail"]
    assert Path(result["outputs"]["pfm"]) == out.resolve()


def test_predict_motifs_without_params_passes_empty_mapping(tmp_path, gimme_inputs, monkeypatch):
    fake = _FakeProgram(pfm=gimme_inputs["pfm"])
    _patch_tools(monkeypatch, lambda name: fake)
    pending.gimmemotifs_predict_motifs(
        "Homer", str(gimme_inputs["sample"]), str(tmp_path / "p.pfm"), verbose=False
    )
    assert fake.calls[0][1] == {}


def test_predict_motifs_missing_program_raises_with_hint(tmp_path, gimme_inputs, monkeypatch):
    def unconfigured(name):
        raise ValueError(f"No configuration found for {name}")

    _patch_tools(monkeypatch, unconfigured)
    out = tmp_path / "never.pfm"
    with pytest.raises(BackendUnavailableError) as info:
        pending.gimmemotifs_predict_motifs(
            "Homer", str(gimme_inputs["sample"]), str(out), verbose=False
        )
    assert info.value.hint and "Homer" in info.value.hint
    assert not out.exists()


def test_predict_motifs_rejects_unknown_tool_and_existing_output(
    tmp_path, gimme_inputs, monkeypatch
):
    with pytest.raises(ConfigurationError, match="Unknown gimmemotifs tool"):
        pending.gimmemotifs_predict_motifs(
            "NotATool", str(gimme_inputs["sample"]), str(tmp_path / "x.pfm"), verbose=False
        )
    _patch_tools(monkeypatch, lambda name: _FakeProgram(pfm=gimme_inputs["pfm"]))
    existing = tmp_path / "taken.pfm"
    existing.write_text("old", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        pending.gimmemotifs_predict_motifs(
            "MEME", str(gimme_inputs["sample"]), str(existing), verbose=False
        )


def test_predict_motifs_real_tool_lookup_when_not_configured(tmp_path, gimme_inputs):
    tools = importlib.import_module("gimmemotifs.tools")
    try:
        available = tools.get_tool("MEME").is_installed()
    except ValueError:
        available = False
    if available:
        pytest.skip("MEME is configured and installed on this machine")
    with pytest.raises(BackendUnavailableError, match="MEME"):
        pending.gimmemotifs_predict_motifs(
            "MEME", str(gimme_inputs["sample"]), str(tmp_path / "m.pfm"), verbose=False
        )


# ------------------------------------------------------ maelstrom HTML report


def _maelstrom_run(tmp_path: Path, gimme_inputs: dict[str, Path]) -> Path:
    run = tmp_path / "maelstrom_run"
    run.mkdir()
    (run / "final.out.txt").write_text(MAELSTROM_TABLE, encoding="utf-8")
    return run


def test_maelstrom_report_writes_html_and_counts_motifs(tmp_path, gimme_inputs):
    run = _maelstrom_run(tmp_path, gimme_inputs)
    result = pending.gimmemotifs_maelstrom_html_report(
        str(run), str(gimme_inputs["pfm"]), verbose=False
    )
    report = run / "gimme.maelstrom.report.html"
    assert result["status"] == "complete"
    assert result["operation"] == "maelstrom_html_report"
    assert result["motifs_in_report"] == 1
    assert result["methods"] == 2
    assert Path(result["outputs"]["report"]) == report
    html = report.read_text(encoding="utf-8")
    assert "MA0001.1" in html and "MA0002.1" not in html
    assert Path(result["outputs"]["logos"]).is_dir()


def test_maelstrom_report_without_logo_column(tmp_path, gimme_inputs):
    run = _maelstrom_run(tmp_path, gimme_inputs)
    result = pending.gimmemotifs_maelstrom_html_report(
        str(run), str(gimme_inputs["pfm"]), plot_no_motifs=True, verbose=False
    )
    assert (run / "gimme.maelstrom.report.html").is_file()
    assert "logos" not in result["outputs"]


def test_maelstrom_report_threshold_filters_motifs(tmp_path, gimme_inputs):
    run = _maelstrom_run(tmp_path, gimme_inputs)
    result = pending.gimmemotifs_maelstrom_html_report(
        str(run), str(gimme_inputs["pfm"]), threshold=0.1, verbose=False
    )
    assert result["motifs_in_report"] == 2


def test_maelstrom_report_requires_final_table_and_valid_threshold(tmp_path, gimme_inputs):
    empty = tmp_path / "empty_run"
    empty.mkdir()
    with pytest.raises(ConfigurationError, match=r"final\.out\.txt"):
        pending.gimmemotifs_maelstrom_html_report(
            str(empty), str(gimme_inputs["pfm"]), verbose=False
        )
    run = _maelstrom_run(tmp_path, gimme_inputs)
    with pytest.raises(ConfigurationError, match="threshold"):
        pending.gimmemotifs_maelstrom_html_report(
            str(run), str(gimme_inputs["pfm"]), threshold=-1, verbose=False
        )


# ------------------------------------------------------------- module contract


def test_public_functions_take_keyword_only_verbose_default_true():
    names = [
        "gimmemotifs_gimme_motifs",
        "gimmemotifs_predict_motifs",
        "gimmemotifs_maelstrom_html_report",
    ]
    assert sorted(pending.__all__) == sorted(names)
    for name in names:
        parameter = inspect.signature(getattr(pending, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True
