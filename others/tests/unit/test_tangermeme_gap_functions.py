"""Tests for dnakit.regulation.tangermeme_gap_functions (tangermeme plots and attribution figures).

Every test runs offline: socket.socket is replaced by a function that raises, so an accidental
network call fails the test. Happy paths need tangermeme, torch and matplotlib; guard paths
remove the module from sys.modules and must raise BackendUnavailableError.
"""

from __future__ import annotations

import importlib.util
import socket
import sys
from functools import partial
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.regulation import tangermeme_gap_functions as gap

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
MPLD3_MISSING = importlib.util.find_spec("mpld3") is None


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def _blocked(*args: object, **kwargs: object) -> None:
        raise RuntimeError("network access is not allowed in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _attribution_rows(length: int = 12) -> list[list[float]]:
    rows = []
    for position in range(length):
        value = 0.1 * (position + 1)
        rows.append([value, -0.05 * position, 0.0, 0.02 * (position % 3)])
    return rows


def _tiny_model(torch: object, length: int) -> object:
    """Conv + ReLU + Linear: one output, a DeepLIFT-compatible architecture."""

    nn = torch.nn  # type: ignore[attr-defined]
    torch.manual_seed(0)  # type: ignore[attr-defined]
    model = nn.Sequential(
        nn.Conv1d(4, 2, kernel_size=3),
        nn.ReLU(),
        nn.Flatten(),
        nn.Linear(2 * (length - 2), 1),
    )
    model.eval()
    return model


# ---------------------------------------------------------------- attribution logo


def test_attribution_logo_writes_png_and_echoes_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("tangermeme")
    mpl_axes = pytest.importorskip("matplotlib.axes")
    limits: list[tuple[float, float]] = []
    original = mpl_axes.Axes.set_ylim

    def _record(self: object, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        low, high = self.get_ylim()  # type: ignore[attr-defined]
        limits.append((float(low), float(high)))
        return result

    monkeypatch.setattr(mpl_axes.Axes, "set_ylim", _record)
    rows = _attribution_rows(12)
    annotations = [{"name": "motif", "start": 2, "end": 6, "score": 1.5}]
    output = tmp_path / "logos" / "attr.png"
    result = gap.plot_attribution_logo(rows, str(output), ylim=1.0, annotations=annotations)
    assert limits[-1] == (-1.0, 1.0)
    assert result["path"] == str(output.resolve())
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    data = result["data"]
    assert data["position_count"] == 12
    assert data["attributions"] == rows
    assert data["ylim"] == [-1.0, 1.0]
    assert data["annotations"] == annotations
    assert data["alphabet"] == ["A", "C", "G", "T"]


def test_attribution_logo_window_rejects_bad_inputs(tmp_path: Path) -> None:
    rows = _attribution_rows(8)
    with pytest.raises(ConfigurationError, match="four numbers"):
        gap.plot_attribution_logo([[0.1, 0.2, 0.3]], str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="smaller than end"):
        gap.plot_attribution_logo(rows, str(tmp_path / "a.png"), start=5, end=5)
    with pytest.raises(ConfigurationError, match="positive bound"):
        gap.plot_attribution_logo(rows, str(tmp_path / "a.png"), ylim=-1.0)
    with pytest.raises(ConfigurationError, match="score"):
        gap.plot_attribution_logo(
            rows, str(tmp_path / "a.png"), annotations=[{"name": "m", "start": 1, "end": 3}]
        )
    with pytest.raises(ConfigurationError, match="start smaller than end"):
        gap.plot_attribution_logo(
            rows,
            str(tmp_path / "a.png"),
            annotations=[{"name": "m", "start": 4, "end": 4, "score": 1.0}],
        )
    with pytest.raises(ConfigurationError, match=r"\.png"):
        gap.plot_attribution_logo(rows, str(tmp_path / "a.svg"))


def test_attribution_logo_guard_without_tangermeme(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("matplotlib")
    monkeypatch.setitem(sys.modules, "tangermeme.plot", None)
    with pytest.raises(BackendUnavailableError) as info:
        gap.plot_attribution_logo(_attribution_rows(6), str(tmp_path / "logo.png"))
    assert "bio-tangermeme" in info.value.hint
    assert not (tmp_path / "logo.png").exists()


# ---------------------------------------------------------------- PWM logo


def test_pwm_logo_writes_png_and_keeps_rows(tmp_path: Path) -> None:
    pytest.importorskip("tangermeme")
    pwm = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.25, 0.25, 0.25, 0.25],
        [0.5, 0.0, 0.5, 0.0],
    ]
    output = tmp_path / "pwm.png"
    result = gap.plot_pwm_logo(pwm, str(output), name="motif")
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    assert result["data"]["pwm"] == pwm
    assert result["data"]["name"] == "motif"
    assert result["data"]["position_count"] == 4


def test_pwm_logo_rejects_bad_rows(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="sum to 1"):
        gap.plot_pwm_logo([[0.5, 0.5, 0.5, 0.0]], str(tmp_path / "p.png"))
    with pytest.raises(ConfigurationError, match="non-negative"):
        gap.plot_pwm_logo([[1.5, -0.5, 0.0, 0.0]], str(tmp_path / "p.png"))
    with pytest.raises(ConfigurationError, match="eps"):
        gap.plot_pwm_logo([[1.0, 0.0, 0.0, 0.0]], str(tmp_path / "p.png"), eps=0.0)


# ---------------------------------------------------------------- model attributions


def test_sequence_attributions_match_direct_deep_lift_shap(tmp_path: Path) -> None:
    torch = pytest.importorskip("torch")
    pytest.importorskip("tangermeme")
    deep_lift = pytest.importorskip("tangermeme.deep_lift_shap")
    utils = pytest.importorskip("tangermeme.utils")
    length = 12
    model = _tiny_model(torch, length)
    sequences = ["ACGTACGTTTGC", "GGGGAAAACCCC"]
    output = tmp_path / "attr.png"
    result = gap.plot_sequence_attributions(
        model, sequences, str(output), target=0, n_shuffles=5, random_state=3
    )
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    data = result["data"]
    assert data["sequence_count"] == 2
    assert data["length"] == length
    assert data["sequences"] == sequences
    values = torch.tensor(data["attributions"], dtype=torch.float32)
    assert tuple(values.shape) == (2, length, 4)

    batch = torch.stack([utils.one_hot_encode(text) for text in sequences]).float()
    direct = deep_lift.deep_lift_shap(
        model,
        batch,
        target=0,
        n_shuffles=5,
        references=partial(deep_lift.dinucleotide_shuffle, random_state=3),
    )
    assert torch.allclose(values, direct.permute(0, 2, 1), atol=1e-5)


def test_sequence_attributions_validation(tmp_path: Path) -> None:
    torch = pytest.importorskip("torch")
    pytest.importorskip("tangermeme")
    model = _tiny_model(torch, 12)
    with pytest.raises(ConfigurationError, match="same length"):
        gap.plot_sequence_attributions(model, ["ACGTACGTACGT", "ACGT"], str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="1 to 32"):
        gap.plot_sequence_attributions(model, [], str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="non-negative"):
        gap.plot_sequence_attributions(model, ["ACGTACGTACGT"], str(tmp_path / "a.png"), target=-1)
    training = _tiny_model(torch, 12)
    training.train()
    with pytest.raises(ConfigurationError, match="eval"):
        gap.plot_sequence_attributions(training, ["ACGTACGTACGT"], str(tmp_path / "a.png"))
    assert not (tmp_path / "a.png").exists()


# ---------------------------------------------------------------- categorical scatter


def test_categorical_scatter_reports_dominant_base(tmp_path: Path) -> None:
    pytest.importorskip("tangermeme")
    rows = [
        [2.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -3.0],
        [0.0, 0.0, 0.5, 0.25],
    ]
    output = tmp_path / "scatter.png"
    result = gap.plot_categorical_attribution_scatter(rows, str(output), size=30)
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    data = result["data"]
    assert data["dominant_base_index_by_position"] == [0, 3, 2]
    assert data["summed_attribution_by_position"] == pytest.approx([2.0, -3.0, 0.75])
    assert data["alphabet"] == ["A", "C", "G", "T"]


def test_categorical_scatter_rejects_bad_size(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="positive"):
        gap.plot_categorical_attribution_scatter(
            _attribution_rows(4), str(tmp_path / "s.png"), size=0
        )


# ---------------------------------------------------------------- interactive logo


def test_interactive_logo_validates_limits_before_mpld3(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="low < high"):
        gap.render_interactive_logo_png(
            _attribution_rows(6), str(tmp_path / "i.png"), ylim=[1.0, -1.0]
        )
    assert not (tmp_path / "i.png").exists()


def test_interactive_logo_guard_without_mpld3(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("tangermeme")
    monkeypatch.setitem(sys.modules, "mpld3", None)
    output = tmp_path / "interactive.png"
    with pytest.raises(BackendUnavailableError) as info:
        gap.render_interactive_logo_png(_attribution_rows(6), str(output))
    assert "mpld3" in info.value.hint
    assert not output.exists()


@pytest.mark.skipif(
    MPLD3_MISSING,
    reason="mpld3 is not installed; it is an optional dependency of tangermeme.interactive_logo",
)
def test_interactive_logo_writes_static_png(tmp_path: Path) -> None:
    output = tmp_path / "interactive.png"
    result = gap.render_interactive_logo_png(_attribution_rows(10), str(output), start=1, end=9)
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    assert result["data"]["start"] == 1
    assert result["data"]["end"] == 9
