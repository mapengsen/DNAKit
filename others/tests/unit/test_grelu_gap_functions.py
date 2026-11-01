"""Tests for dnakit.predictions.grelu_gap_functions (gReLU wrappers).

Network-backed gReLU functions are replaced by fakes here, so no test reaches Hugging Face, Weights & Biases or a
genome server. Plots and the tomtom comparison run locally (matplotlib and memelite are part of the environment).
"""

from __future__ import annotations

import inspect
import os
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")
pytest.importorskip("grelu")

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.predictions import grelu_gap_functions as facade  # noqa: E402

_PNG_HEADER = b"\x89PNG\r\n\x1a\n"

_MEME = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF T1 target_one
letter-probability matrix: alength= 4 w= 3 nsites= 20 E= 0
  0.9 0.05 0.03 0.02
  0.05 0.9 0.03 0.02
  0.02 0.03 0.9 0.05

MOTIF T2 target_two
letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 0
  0.02 0.03 0.05 0.9
  0.9 0.05 0.03 0.02
"""


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(facade.__all__) == len(set(facade.__all__))
    for name in facade.__all__:
        verbose = inspect.signature(getattr(facade, name)).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True


def test_missing_package_raises_with_hint() -> None:
    with pytest.raises(BackendUnavailableError) as info:
        facade._upstream("dnakit_missing_package_for_test", "hint text")
    assert info.value.hint == "hint text"


def test_model_zoo_wrappers_return_upstream_results(monkeypatch: pytest.MonkeyPatch) -> None:
    import grelu.resources as resources

    monkeypatch.setattr(resources, "list_models", lambda: ["org/model-a"])
    assert facade.grelu_list_models() == ["org/model-a"]

    seen: dict[str, object] = {}

    def fake_download(repo_id, filename="model.ckpt", **kwargs):
        seen.update(repo_id=repo_id, filename=filename)
        return "/cache/model.ckpt"

    monkeypatch.setattr(resources, "download_model", fake_download)
    assert facade.grelu_download_model("org/model-a") == "/cache/model.ckpt"
    assert seen == {"repo_id": "org/model-a", "filename": "model.ckpt"}

    monkeypatch.setattr(resources, "get_base_models", lambda repo_id: ["base/x"])
    assert facade.grelu_get_base_models("org/model-a") == ["base/x"]


def test_wandb_wrappers_pass_project_and_host(monkeypatch: pytest.MonkeyPatch) -> None:
    import grelu.resources.wandb as zoo

    seen: dict[str, object] = {}

    def fake_artifacts(project, host, type_is=None, type_contains=None):
        seen.update(project=project, host=host, type_is=type_is)
        return ["artifact:v1"]

    monkeypatch.setattr(zoo, "artifacts", fake_artifacts)
    assert facade.grelu_wandb_artifacts("zoo", type_is="model") == ["artifact:v1"]
    assert seen == {"project": "zoo", "host": "https://api.wandb.ai", "type_is": "model"}


def test_genome_wrappers_forward_to_grelu(monkeypatch: pytest.MonkeyPatch) -> None:
    import grelu.io.genome as genome

    monkeypatch.setattr(genome, "read_sizes", lambda genome="hg38": {"chr1": 248956422})
    assert facade.grelu_read_sizes("hg38") == {"chr1": 248956422}

    def failing(genome, **kwargs):
        raise RuntimeError("no network")

    monkeypatch.setattr(genome, "get_genome", failing)
    with pytest.raises(BackendExecutionError, match="no network"):
        facade.grelu_get_genome("hg38")


def test_blank_identifiers_are_rejected() -> None:
    with pytest.raises(ConfigurationError):
        facade.grelu_download_model("  ")
    with pytest.raises(ConfigurationError):
        facade.grelu_wandb_artifacts("")


def test_plot_output_must_be_an_image_path(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.grelu_plot_distribution([1.0, 2.0], str(tmp_path / "plot.txt"))
    with pytest.raises(ConfigurationError):
        facade.grelu_plot_motif([0.1, 0.2, 0.3, 0.4], str(tmp_path / "motif.png"))  # must be a 2-D position matrix


def test_distribution_plot_saves_image_and_data(tmp_path: Path) -> None:
    target = tmp_path / "distribution.png"
    result = facade.grelu_plot_distribution([0.1, 0.4, 0.4, 0.9], str(target), title="test metric")
    assert result["status"] == "completed"
    assert result["function"] == "grelu.visualize.plot_distribution"
    assert Path(result["image_path"]).read_bytes()[:8] == _PNG_HEADER
    assert result["data"]["values"] == [0.1, 0.4, 0.4, 0.9]


def test_prediction_scatter_plot_saves_image(tmp_path: Path) -> None:
    # gReLU expects (sequences, tasks, positions)
    preds = [[[0.1], [0.8]], [[0.4], [0.6]], [[0.7], [0.2]], [[0.9], [0.1]]]
    labels = [[[0.0], [1.0]], [[0.0], [1.0]], [[1.0], [0.0]], [[1.0], [0.0]]]
    target = tmp_path / "scatter.png"
    result = facade.grelu_plot_pred_scatter(preds, labels, str(target), figsize=(3, 3))
    assert Path(result["image_path"]).read_bytes()[:8] == _PNG_HEADER
    assert result["data"]["preds"]["shape"] == [4, 2, 1]


def test_attribution_and_motif_plots(tmp_path: Path) -> None:
    attrs = [[0.1 * (i + j) for j in range(12)] for i in range(4)]
    attr_result = facade.grelu_plot_attributions(attrs, str(tmp_path / "attr.png"), figsize=(6, 2))
    assert Path(attr_result["image_path"]).read_bytes()[:8] == _PNG_HEADER
    motif = [[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1]]
    motif_result = facade.grelu_plot_motif(motif, str(tmp_path / "motif.png"), name="demo")
    assert Path(motif_result["image_path"]).read_bytes()[:8] == _PNG_HEADER


def test_plots_do_not_leave_figures_open(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt

    before = set(plt.get_fignums())
    facade.grelu_plot_distribution([1.0, 2.0, 3.0], str(tmp_path / "left.png"))
    assert set(plt.get_fignums()) == before


def test_tomtom_compares_query_motif_to_local_meme_file(tmp_path: Path) -> None:
    pytest.importorskip("memelite")
    meme = tmp_path / "targets.meme"
    meme.write_text(_MEME, encoding="utf-8")
    query = {"Q1": [[0.9, 0.05, 0.03, 0.02], [0.05, 0.9, 0.03, 0.02], [0.02, 0.03, 0.9, 0.05]]}
    table = facade.grelu_run_tomtom(query, str(meme))
    assert len(table) == 2  # one query against two target motifs
    assert set(table["Query_ID"]) == {"Q1"}
    assert {"Target_ID", "p-value", "q-value", "Orientation"} <= set(table.columns)


def test_tomtom_rejects_bad_motifs(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.grelu_run_tomtom({}, "consensus")
    with pytest.raises(ConfigurationError):
        facade.grelu_run_tomtom({"bad": [[0.5, 0.5]]}, "consensus")
