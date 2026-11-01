"""Tests for dnakit.predictions.selene_sdk_gap_functions (selene-sdk plots, pipelines, predictions).

Every test runs offline: socket.socket is replaced by a function that raises, so an accidental
network call fails the test. Happy paths need selene-sdk, matplotlib and torch; guard paths
remove the selene module from sys.modules and must raise BackendUnavailableError.
"""

from __future__ import annotations

import socket
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import selene_sdk_gap_functions as gap

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def _blocked(*args: object, **kwargs: object) -> None:
        raise RuntimeError("network access is not allowed in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _score_matrix() -> list[list[float]]:
    return [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.5, 0.5],
        [0.2, 0.1, 0.3, -0.4],
    ]


def _curve_tables() -> tuple[list[list[float]], list[list[int]]]:
    """Feature 0 is ranked perfectly (AUC 1.0, 4 positives); feature 1 has 1; feature 2 has 0."""

    scores = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
    labels_f0 = [1, 1, 1, 1, 0, 0, 0, 0]
    labels_f1 = [1, 0, 0, 0, 0, 0, 0, 0]
    labels_f2 = [0] * 8
    prediction = [[scores[i], 0.5, 0.1] for i in range(8)]
    target = [[labels_f0[i], labels_f1[i], labels_f2[i]] for i in range(8)]
    return prediction, target


def _count_model(path: Path, length: int = 4) -> None:
    """Whole torch model: feature 0 sums the A channel, feature 1 sums the T channel."""

    torch = pytest.importorskip("torch")
    linear = torch.nn.Linear(4 * length, 2)
    with torch.no_grad():
        linear.weight.zero_()
        linear.bias.zero_()
        for position in range(length):
            linear.weight[0, 0 * length + position] = 1.0  # channel A
            linear.weight[1, 3 * length + position] = 1.0  # channel T
    torch.save(torch.nn.Sequential(torch.nn.Flatten(), linear), path)


# ---------------------------------------------------------------- heatmap


def test_score_heatmap_writes_png_with_mask(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    matrix = _score_matrix()
    mask = [[False, False, False, False]] * 3 + [[False, True, False, False]]
    output = tmp_path / "heat" / "heatmap.png"
    result = gap.selene_score_heatmap(matrix, str(output), mask=mask, cmap="viridis")
    assert output.read_bytes()[:8] == PNG_SIGNATURE
    assert result["data"]["score_matrix"] == matrix
    assert result["data"]["mask"] == mask
    assert result["data"]["row_labels_top_to_bottom"] == ["T", "G", "C", "A"]
    assert result["data"]["cmap"] == "viridis"


def test_score_heatmap_rejects_bad_mask(tmp_path: Path) -> None:
    matrix = _score_matrix()
    with pytest.raises(ConfigurationError, match="same number of rows"):
        gap.selene_score_heatmap(matrix, str(tmp_path / "h.png"), mask=[[False] * 4])
    with pytest.raises(ConfigurationError, match="only True or False"):
        gap.selene_score_heatmap(matrix, str(tmp_path / "h.png"), mask=[[0, 0, 0, 0]] * 4)


# ---------------------------------------------------------------- ROC and precision-recall


def test_roc_curves_known_case(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    matplotlib = pytest.importorskip("matplotlib")
    pyplot = pytest.importorskip("matplotlib.pyplot")
    prediction, target = _curve_tables()
    output_dir = tmp_path / "roc"
    result = gap.selene_roc_curves(
        prediction, target, str(output_dir), report_gt_feature_n_positives=2
    )
    written = output_dir / "roc_curves.svg"
    assert result["path"] == str(written.resolve())
    assert b"<svg" in written.read_bytes()
    data = result["data"]
    assert data["positives_per_feature"] == [4, 1, 0]
    assert data["plotted_features"] == [0]
    assert data["metric"] == "roc_auc"
    assert data["metric_per_plotted_feature"] == pytest.approx([1.0])
    assert data["report_gt_feature_n_positives"] == 2
    # selene switches to the SVG backend and leaves figures open; the wrapper must undo both.
    assert matplotlib.get_backend().lower() == "agg"
    assert pyplot.get_fignums() == []


def test_precision_recall_curves_known_case(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    pytest.importorskip("matplotlib")
    prediction, target = _curve_tables()
    output_dir = tmp_path / "pr"
    result = gap.selene_precision_recall_curves(
        prediction, target, str(output_dir), report_gt_feature_n_positives=2
    )
    assert result["path"] == str((output_dir / "precision_recall_curves.svg").resolve())
    assert b"<svg" in (output_dir / "precision_recall_curves.svg").read_bytes()
    assert result["data"]["plotted_features"] == [0]
    assert result["data"]["metric"] == "average_precision"
    assert result["data"]["metric_per_plotted_feature"] == pytest.approx([1.0])


def test_curve_inputs_are_validated(tmp_path: Path) -> None:
    prediction, target = _curve_tables()
    bad_target = [row[:] for row in target]
    bad_target[0][0] = 2
    with pytest.raises(ConfigurationError, match="0 and 1"):
        gap.selene_roc_curves(prediction, bad_target, str(tmp_path / "a"))
    with pytest.raises(ConfigurationError, match="same number of rows"):
        gap.selene_roc_curves(prediction[:4], target, str(tmp_path / "b"))
    with pytest.raises(ConfigurationError, match="at least two"):
        gap.selene_roc_curves([[0.5]], [[1]], str(tmp_path / "c"))
    with pytest.raises(ConfigurationError, match="report_gt_feature_n_positives"):
        gap.selene_roc_curves(
            prediction, target, str(tmp_path / "d"), report_gt_feature_n_positives=-1
        )


def test_curve_refuses_to_overwrite_existing_svg(tmp_path: Path) -> None:
    prediction, target = _curve_tables()
    output_dir = tmp_path / "exists"
    output_dir.mkdir()
    (output_dir / "roc_curves.svg").write_text("old")
    with pytest.raises(ConfigurationError, match="already exists"):
        gap.selene_roc_curves(prediction, target, str(output_dir))
    assert (output_dir / "roc_curves.svg").read_text() == "old"


def test_curve_guard_without_selene(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("matplotlib")
    prediction, target = _curve_tables()
    monkeypatch.setitem(sys.modules, "selene_sdk.utils.performance_metrics", None)
    with pytest.raises(BackendUnavailableError):
        gap.selene_roc_curves(prediction, target, str(tmp_path / "guard"))


# ---------------------------------------------------------------- pipelines (YAML / parsed dict)


def test_run_selene_configs_writes_run_directory(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    import yaml

    model = tmp_path / "model.pt"
    model.write_bytes(b"not loaded: ops are empty")
    configs = {"ops": [], "model": {"path": str(model)}}
    output = tmp_path / "runs"
    result = gap.run_selene_configs(configs, str(output))
    run_directory = Path(result["run_directory"])
    assert run_directory.parent == output.resolve()
    assert "configuration.yaml" in result["files"]
    assert "model.pt" in result["files"]
    saved = yaml.safe_load((run_directory / "configuration.yaml").read_text())
    assert saved["output_dir"] == str(run_directory)
    assert result["operations"] == []
    # The caller's dict is not changed.
    assert configs == {"ops": [], "model": {"path": str(model)}}


def test_run_selene_config_file_reads_yaml(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    import yaml

    model = tmp_path / "model.pt"
    model.write_bytes(b"placeholder")
    config = tmp_path / "pipeline.yml"
    config.write_text(
        yaml.safe_dump({"ops": [], "model": {"path": str(model)}, "output_dir": "ignored"})
    )
    output = tmp_path / "out"
    result = gap.run_selene_config_file(str(config), str(output), create_subdirectory=False)
    assert result["run_directory"] == str(output.resolve())
    assert result["source"] == str(config.resolve())
    saved = yaml.safe_load((output / "configuration.yaml").read_text())
    assert saved["output_dir"] == str(output.resolve())
    assert "selene_sdk_version" in saved


def test_pipeline_validation_errors(tmp_path: Path) -> None:
    model = tmp_path / "model.pt"
    model.write_bytes(b"x")
    with pytest.raises(ConfigurationError, match="does not exist"):
        gap.run_selene_config_file(str(tmp_path / "missing.yml"), str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="operations"):
        gap.run_selene_configs({"model": {"path": str(model)}}, str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="operations"):
        gap.run_selene_configs({"ops": ["fit"], "model": {"path": str(model)}}, str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="learning rate"):
        gap.run_selene_configs(
            {"ops": ["train"], "model": {"path": str(model)}}, str(tmp_path / "o")
        )
    with pytest.raises(ConfigurationError, match="model"):
        gap.run_selene_configs({"ops": []}, str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="does not exist"):
        gap.run_selene_configs(
            {"ops": [], "model": {"path": str(tmp_path / "nope.pt")}}, str(tmp_path / "o")
        )
    with pytest.raises(ConfigurationError, match="string keys"):
        gap.run_selene_configs(["not", "a", "dict"], str(tmp_path / "o"))  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="lr must be positive"):
        gap.run_selene_configs(
            {"ops": [], "model": {"path": str(model)}}, str(tmp_path / "o"), lr=0.0
        )


def test_run_selene_operations_empty_list_is_a_no_op(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    model = tmp_path / "model.pt"
    model.write_bytes(b"x")
    output = tmp_path / "ops-out"
    result = gap.run_selene_operations([], {"model": {"path": str(model)}}, str(output))
    assert result["operations"] == []
    assert output.is_dir()
    assert result["output_dir"] == str(output.resolve())


def test_run_selene_operations_validation(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="lr"):
        gap.run_selene_operations(["train"], {}, str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="operations"):
        gap.run_selene_operations(["predict"], {}, str(tmp_path / "o"))
    with pytest.raises(ConfigurationError, match="string keys"):
        gap.run_selene_operations([], "not a mapping", str(tmp_path / "o"))  # type: ignore[arg-type]


def test_pipeline_guard_without_selene(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "selene_sdk.utils.config_utils", None)
    with pytest.raises(BackendUnavailableError) as info:
        gap.run_selene_operations([], {}, str(tmp_path / "guard"))
    assert "bio-selene" in info.value.hint


# ---------------------------------------------------------------- whole-model prediction


def test_predict_counts_bases_with_known_answer(tmp_path: Path) -> None:
    pytest.importorskip("selene_sdk")
    model_file = tmp_path / "counts.pt"
    _count_model(model_file)
    result = gap.predict_selene_sequences(str(model_file), ["AAGT", "TTTT", "ACGN"])
    assert result["shape"] == [3, 2]
    assert result["sequence_count"] == 3
    assert result["length"] == 4
    predictions = result["predictions"]
    assert predictions[0] == pytest.approx([2.0, 1.0])
    assert predictions[1] == pytest.approx([0.0, 4.0])
    # N is a uniform 0.25 row in selene's encoding.
    assert predictions[2] == pytest.approx([1.25, 0.25])


def test_predict_rejects_bad_inputs(tmp_path: Path) -> None:
    torch = pytest.importorskip("torch")
    model_file = tmp_path / "counts.pt"
    _count_model(model_file)
    with pytest.raises(ConfigurationError, match="same length"):
        gap.predict_selene_sequences(str(model_file), ["AAAA", "AAA"])
    with pytest.raises(ConfigurationError, match="does not exist"):
        gap.predict_selene_sequences(str(tmp_path / "missing.pt"), ["AAAA"])
    state_file = tmp_path / "state.pt"
    torch.save({"weight": torch.zeros(2)}, state_file)
    with pytest.raises(ConfigurationError, match=r"whole torch\.nn\.Module"):
        gap.predict_selene_sequences(str(state_file), ["AAAA"])
    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        gap.predict_selene_sequences(str(model_file), ["AAXA"])


def test_predict_guard_without_selene(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    model_file = tmp_path / "counts.pt"
    _count_model(model_file)
    monkeypatch.setitem(sys.modules, "selene_sdk.predict.model_predict", None)
    with pytest.raises(BackendUnavailableError):
        gap.predict_selene_sequences(str(model_file), ["AAAA"])
