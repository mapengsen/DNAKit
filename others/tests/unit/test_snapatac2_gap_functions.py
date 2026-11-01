"""Tests for snapatac2 registry downloads and plots (dnakit.epigenome.snapatac2_gap_functions).

Downloads are tested with fakes that replace the upstream loaders, so no network is used. Plot tests build
small synthetic inputs in tmp_path and run the real snapatac2 plotting functions.
"""

from __future__ import annotations

import importlib
import socket
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("snapatac2")

from dnakit.epigenome import snapatac2_gap_functions as module  # noqa: E402
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, DownloadError  # noqa: E402


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


def _motif(identifier: str, name: str, family: str | None = None):
    values = {"id": identifier, "name": name}
    if family is not None:
        values["family"] = family
    return SimpleNamespace(**values)


def test_exports_are_public_functions_with_verbose_last():
    import inspect

    assert len(module.__all__) == 24
    for name in module.__all__:
        function = getattr(module, name)
        assert callable(function)
        parameters = list(inspect.signature(function).parameters.values())
        assert parameters[-1].name == "verbose"
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY
        assert parameters[-1].default is True


def test_load_dataset_forwards_name_and_file(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    calls = []
    target = tmp_path / "fake.tsv.gz"

    def fake_load(dataset, **kwargs):
        calls.append((dataset, kwargs))
        return target

    monkeypatch.setattr(datasets, "load", fake_load)
    result = module.load_snapatac2_dataset("atac_pbmc_500", file="atac_pbmc_500.tsv.gz")
    assert calls == [("atac_pbmc_500", {"file": "atac_pbmc_500.tsv.gz"})]
    assert result["status"] == "completed"
    assert result["dataset"] == "atac_pbmc_500"
    assert result["paths"] == [str(target)]
    assert result["retrieved_utc"].endswith("+00:00")


def test_load_dataset_reports_non_path_objects(monkeypatch):
    import snapatac2.datasets as datasets

    monkeypatch.setattr(datasets, "load", lambda dataset, **kwargs: [_motif("a", "A")] * 3)
    result = module.load_snapatac2_dataset("cisBP_human")
    assert result["paths"] == []
    assert result["data"] == {"object_type": "list", "object_count": 3}


def test_load_dataset_rejects_blank_name():
    with pytest.raises(ConfigurationError):
        module.load_snapatac2_dataset("   ")


def test_download_failure_is_reported_as_download_error(monkeypatch):
    import snapatac2.datasets as datasets

    def fake_load(dataset, **kwargs):
        raise OSError("connection refused")

    monkeypatch.setattr(datasets, "load", fake_load)
    with pytest.raises(DownloadError) as caught:
        module.load_snapatac2_dataset("HEA_cCRE")
    assert "SNAP_DATA_DIR" in (caught.value.hint or "")


def test_unknown_dataset_is_configuration_error(monkeypatch):
    import snapatac2.datasets as datasets

    def fake_load(dataset, **kwargs):
        raise KeyError(dataset)

    monkeypatch.setattr(datasets, "load", fake_load)
    with pytest.raises(ConfigurationError):
        module.load_snapatac2_dataset("no_such_dataset")


def test_pbmc500_forwards_type_and_downsample(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    calls = []
    target = tmp_path / "atac_pbmc_500_downsample.tsv.gz"
    monkeypatch.setattr(datasets, "pbmc500", lambda **kwargs: calls.append(kwargs) or target)
    result = module.fetch_snapatac2_pbmc500(data_type="fragment", downsample=True)
    assert calls == [{"type": "fragment", "downsample": True}]
    assert result["paths"] == [str(target)]
    assert result["parameters"] == {"data_type": "fragment", "downsample": True}


def test_pbmc500_fastq_returns_all_extracted_paths(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    files = [tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"]
    monkeypatch.setattr(datasets, "pbmc500", lambda **kwargs: files)
    result = module.fetch_snapatac2_pbmc500(data_type="fastq")
    assert result["paths"] == [str(path) for path in files]


def test_pbmc500_rejects_downsample_without_fragments():
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_pbmc500(data_type="bam", downsample=True)
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_pbmc500(data_type="vcf")


def test_pbmc5k_and_multiome_validate_and_forward(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    calls = []
    target = tmp_path / "atac_pbmc_5k_annotated.h5ad"
    monkeypatch.setattr(datasets, "pbmc5k", lambda **kwargs: calls.append(("5k", kwargs)) or target)
    monkeypatch.setattr(
        datasets, "pbmc10k_multiome", lambda **kwargs: calls.append(("10k", kwargs)) or target
    )
    assert module.fetch_snapatac2_pbmc5k(data_type="annotated_h5ad")["paths"] == [str(target)]
    assert module.fetch_snapatac2_pbmc10k_multiome(modality="ATAC", data_type="fragment")["paths"] == [str(target)]
    assert calls == [
        ("5k", {"type": "annotated_h5ad"}),
        ("10k", {"modality": "ATAC", "type": "fragment"}),
    ]
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_pbmc5k(data_type="bam")
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_pbmc10k_multiome(modality="protein")


def test_colon_and_hea_return_sample_rows(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    first = tmp_path / "Colon1_fragments.tsv.gz"
    second = tmp_path / "Colon2_fragments.tsv.gz"
    monkeypatch.setattr(datasets, "colon", lambda: [("Colon1", first), ("Colon2", second)])
    monkeypatch.setattr(datasets, "cre_HEA", lambda: tmp_path / "HEA_cCRE.bed.gz")

    colon = module.fetch_snapatac2_colon_transverse()
    assert colon["data"]["samples"] == [
        {"sample": "Colon1", "fragment_file": str(first)},
        {"sample": "Colon2", "fragment_file": str(second)},
    ]
    assert colon["paths"] == [str(first), str(second)]

    hea = module.fetch_snapatac2_hea_ccre()
    assert hea["paths"] == [str(tmp_path / "HEA_cCRE.bed.gz")]


def test_cis_bp_and_meuleman_return_motif_records(monkeypatch):
    import snapatac2.datasets as datasets

    seen = []
    monkeypatch.setattr(
        datasets, "cis_bp", lambda unique: seen.append(unique) or [_motif("MA0001.1+X", "MA0001.1")]
    )
    monkeypatch.setattr(
        datasets,
        "Meuleman_2020",
        lambda: [_motif("ZNF1_fam", "ZNF1", family="fam")],
    )
    cis = module.fetch_snapatac2_cis_bp_motifs(unique=False)
    assert seen == [False]
    assert cis["data"]["motif_count"] == 1
    assert cis["data"]["motifs"] == [{"id": "MA0001.1+X", "name": "MA0001.1"}]
    assert cis["paths"] == []

    meuleman = module.fetch_snapatac2_meuleman_2020_motifs()
    assert meuleman["data"]["motifs"] == [{"id": "ZNF1_fam", "name": "ZNF1", "family": "fam"}]


def test_genome_files_call_genome_loader_with_preset_names(monkeypatch, tmp_path):
    import snapatac2.genome as genome

    calls = []

    def fake_load(dataset, **kwargs):
        calls.append((dataset, kwargs.get("file"), type(kwargs.get("processor")).__name__))
        return tmp_path / kwargs["file"].replace(".gz", "")

    monkeypatch.setattr(genome, "load", fake_load)
    result = module.fetch_snapatac2_genome_files("hg38")
    assert calls == [
        ("gencode_GRCh38", "gencode_v41_GRCh38.fa.gz", "Decompress"),
        ("gencode_GRCh38", "gencode_v41_GRCh38.gff3.gz", "NoneType"),
    ]
    assert result["parameters"]["genome"] == "GRCh38"
    assert result["parameters"]["requested_name"] == "hg38"
    assert result["data"]["fasta"] == str(tmp_path / "gencode_v41_GRCh38.fa")
    assert result["paths"] == [result["data"]["fasta"], result["data"]["annotation"]]


def test_genome_files_single_kind_and_bad_input(monkeypatch, tmp_path):
    import snapatac2.genome as genome

    calls = []
    monkeypatch.setattr(
        genome, "load", lambda dataset, **kwargs: calls.append(kwargs["file"]) or tmp_path / "x.gff3"
    )
    result = module.fetch_snapatac2_genome_files("mm10", kind="annotation")
    assert calls == ["gencode_vM25_GRCm38.gff3.gz"]
    assert list(result["data"]) == ["annotation"]
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_genome_files("T2T")
    with pytest.raises(ConfigurationError):
        module.fetch_snapatac2_genome_files("GRCh38", kind="chromsizes")


def test_missing_snapatac2_raises_backend_unavailable_with_hint(monkeypatch):
    real_import = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name.startswith("snapatac2"):
            raise ImportError("No module named 'snapatac2'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    with pytest.raises(BackendUnavailableError) as caught:
        module.fetch_snapatac2_hea_ccre()
    assert "bio-snapatac2" in (caught.value.hint or "")


def test_verbose_keyword_is_accepted(monkeypatch, tmp_path):
    import snapatac2.datasets as datasets

    monkeypatch.setattr(datasets, "cre_HEA", lambda: tmp_path / "HEA_cCRE.bed.gz")
    assert module.fetch_snapatac2_hea_ccre(verbose=False)["status"] == "completed"


# ---------------------------------------------------------------------------
# Plot functions: real snapatac2 plotting on small synthetic inputs; output PNGs are written to tmp_path.
# ---------------------------------------------------------------------------

_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _png_size(path: Path) -> tuple[int, int]:
    raw = path.read_bytes()
    assert raw[:8] == _PNG_SIGNATURE
    return int.from_bytes(raw[16:20], "big"), int.from_bytes(raw[20:24], "big")


def _trace(result, kind: str) -> list[dict]:
    return [trace for trace in result["data"]["traces"] if trace["type"] == kind]


def _assert_png(result, path: Path) -> None:
    assert result["status"] == "completed"
    assert result["image_path"] == str(path)
    assert path.is_file() and path.stat().st_size > 0
    assert result["image_bytes"] == path.stat().st_size
    assert _png_size(path)[0] > 0


def test_render_figure_png_uses_figure_size(tmp_path):
    target = tmp_path / "figure.png"
    figure = {"data": [{"type": "scatter", "x": [0, 1, 2], "y": [1, 3, 2], "mode": "lines+markers"}],
              "layout": {"title": {"text": "Check"}}}
    result = module.render_snapatac2_figure_png(figure, str(target), width=640, height=320)
    _assert_png(result, target)
    assert _png_size(target) == (640, 320)
    assert result["data"]["title"] == "Check"
    assert _trace(result, "scatter")[0]["y"] == [1, 3, 2]


def test_plot_outputs_must_be_new_png_files(tmp_path):
    existing = tmp_path / "exists.png"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        module.render_snapatac2_figure_png({"data": [{"type": "scatter", "x": [0], "y": [0]}]}, str(existing))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_valid_cells_png([3, 2, 1], str(tmp_path / "figure.pdf"))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_valid_cells_png([3, 2, 1], str(tmp_path / "missing_folder" / "figure.png"))


def test_render_figure_rejects_non_figures(tmp_path):
    with pytest.raises(ConfigurationError):
        module.render_snapatac2_figure_png([1, 2, 3], str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError):
        module.render_snapatac2_figure_png({"layout": {}}, str(tmp_path / "b.png"))


def test_heatmap_png_clusters_and_labels_rows(tmp_path):
    target = tmp_path / "heatmap.png"
    matrix = [[0.0, 1.0, 2.0], [0.1, 1.1, 2.1], [5.0, 0.0, 0.5]]
    result = module.plot_snapatac2_heatmap_png(
        matrix, str(target), row_names=["r1", "r2", "r3"], column_names=["c1", "c2", "c3"],
        cluster_rows=True, cluster_columns=True,
    )
    _assert_png(result, target)
    heat = _trace(result, "heatmap")[0]
    assert len(heat["z"]) == 3 and len(heat["z"][0]) == 3
    assert sorted(result["data"]["yaxis_ticktext"]) == ["r1", "r2", "r3"]
    assert sorted(result["data"]["xaxis_ticktext"]) == ["c1", "c2", "c3"]


def test_heatmap_png_without_clustering_keeps_input_order(tmp_path):
    target = tmp_path / "plain.png"
    result = module.plot_snapatac2_heatmap_png(
        [[1.0, 2.0], [3.0, 4.0]], str(target), row_names=["a", "b"], column_names=["x", "y"],
        cluster_rows=False, cluster_columns=False,
    )
    _assert_png(result, target)
    assert result["data"]["yaxis_ticktext"] == ["a", "b"]
    assert result["data"]["xaxis_ticktext"] == ["x", "y"]
    assert _trace(result, "heatmap")[0]["z"] == [[1.0, 2.0], [3.0, 4.0]]


def test_heatmap_validation(tmp_path):
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_heatmap_png([[1.0, 2.0]], str(tmp_path / "a.png"), cluster_rows=True)
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_heatmap_png([[1.0, 2.0], [3.0]], str(tmp_path / "b.png"))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_heatmap_png([[1.0, 2.0], [3.0, 4.0]], str(tmp_path / "c.png"), linkage="bogus")
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_heatmap_png([[1.0, 2.0], [3.0, 4.0]], str(tmp_path / "d.png"), row_names=["only"])


def test_kde2d_png_draws_density_contour(tmp_path):
    rng = np.random.default_rng(0)
    target = tmp_path / "kde.png"
    result = module.plot_snapatac2_kde2d_png(rng.gamma(2, 1000, 300), rng.normal(5, 1, 300), str(target), log_x=True)
    _assert_png(result, target)
    assert _trace(result, "contour")
    assert result["data"]["xaxis_type"] == "log"
    assert result["parameters"]["points"] == 300


def test_kde2d_rejects_non_positive_values_on_log_axis(tmp_path):
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_kde2d_png([1.0, 0.0, 2.0], [1.0, 2.0, 3.0], str(tmp_path / "k.png"), log_x=True)
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_kde2d_png([1.0, 2.0], [1.0, 2.0, 3.0], str(tmp_path / "k2.png"))


def test_scatter_png_groups_points_by_colour(tmp_path):
    target = tmp_path / "scatter.png"
    result = module.plot_snapatac2_scatter_png(
        [0.0, 1.0, 2.0, 3.0], [0.2, 0.8, 1.4, 0.1], str(target), color=["A", "B", "A", "B"],
        x_label="UMAP-1", y_label="UMAP-2",
    )
    _assert_png(result, target)
    groups = {trace["name"]: trace for trace in _trace(result, "scatter")}
    assert set(groups) == {"A", "B"}
    assert groups["A"]["x"] == [0.0, 2.0] and groups["B"]["y"] == [0.8, 0.1]
    assert result["data"]["xaxis_title"] == "UMAP-1"


def test_scatter_png_numeric_colour_draws_colour_scale(tmp_path):
    target = tmp_path / "numeric.png"
    result = module.plot_snapatac2_scatter_png([0.0, 1.0, 2.0], [0.2, 0.8, 1.4], str(target), color=[1.0, 2.0, 3.0])
    _assert_png(result, target)
    assert _trace(result, "scatter")[0]["y"] == [0.2, 0.8, 1.4]


def test_scatter_png_rejects_mismatched_lengths(tmp_path):
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_scatter_png([0.0, 1.0], [0.2], str(tmp_path / "s.png"))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_scatter_png([0.0, 1.0], [0.2, 0.3], str(tmp_path / "s2.png"), color=["A"])


def test_scatter3d_png_draws_three_dimensional_points(tmp_path):
    target = tmp_path / "scatter3d.png"
    result = module.plot_snapatac2_scatter3d_png(
        [0.0, 1.0, 2.0], [0.2, 0.8, 1.4], [1.0, 0.5, 0.0], str(target), color=["A", "B", "A"]
    )
    _assert_png(result, target)
    groups = _trace(result, "scatter3d")
    assert sorted(trace["name"] for trace in groups) == ["A", "B"]
    assert groups[0]["z"] in ([1.0, 0.0], [0.5])


def _network(edges_extra: list[dict] | None = None):
    edges = [
        {"source": 0, "target": 1, "distance": 1500, "cor_score": 0.7, "regr_score": 0.2},
        {"source": 0, "target": 1, "distance": 50000, "cor_score": 0.3, "regr_score": 0.9},
        {"source": 1, "target": 0, "distance": 800, "cor_score": 0.5, "regr_score": 0.4},
    ]
    return ["peak", "gene"], edges + (edges_extra or [])


def test_network_scores_png_averages_scores_by_distance_bin(tmp_path):
    types, edges = _network()
    target = tmp_path / "scores.png"
    result = module.plot_snapatac2_network_scores_png(types, edges, str(target), score_name="regr_score")
    _assert_png(result, target)
    bars = _trace(result, "bar")[0]
    values = dict(zip(bars["x"], bars["y"]))
    assert values["500 - 2K"] == pytest.approx(0.3)
    assert values["50K - 100K"] == pytest.approx(0.9)


def test_network_edge_stat_png_draws_violins_per_node_type_pair(tmp_path):
    types, edges = _network()
    target = tmp_path / "edges.png"
    result = module.plot_snapatac2_network_edge_stat_png(types, edges, str(target))
    _assert_png(result, target)
    names = sorted(trace["name"] for trace in _trace(result, "violin"))
    assert names == ["gene -> peak", "peak -> gene"]


def test_network_inputs_are_validated(tmp_path):
    types, edges = _network()
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_network_scores_png(types, edges, str(tmp_path / "a.png"), score_name="weight")
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_network_scores_png(
            types, [{"source": 0, "target": 5, "distance": 1.0}], str(tmp_path / "b.png"), score_name="cor_score"
        )
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_network_scores_png(
            types, [{"source": 0, "target": 1, "regr_score": 0.5}], str(tmp_path / "c.png"), score_name="regr_score"
        )
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_network_edge_stat_png(types, [{"source": 0, "target": 1}], str(tmp_path / "d.png"))


def test_valid_cells_png_plots_ranked_counts(tmp_path):
    counts = [10000, 8500, 6200, 1200, 800, 120, 60, 20]
    target = tmp_path / "valid.png"
    result = module.plot_snapatac2_valid_cells_png(counts, str(target), width=500, height=400)
    _assert_png(result, target)
    assert _png_size(target) == (500, 400)
    trace = _trace(result, "scatter")[0]
    assert len(trace["x"]) == len(trace["y"]) == len(set(counts))
    assert result["data"]["xaxis_type"] == "log" and result["data"]["yaxis_type"] == "log"
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_valid_cells_png([5, 0, 1], str(tmp_path / "zero.png"))


def test_tsse_png_counts_cells_above_fragment_threshold(tmp_path):
    rng = np.random.default_rng(1)
    fragments = rng.integers(100, 20000, 300)
    scores = rng.normal(6, 1, 300)
    target = tmp_path / "tsse.png"
    result = module.plot_snapatac2_tsse_png(fragments, scores, str(target), min_fragment=500)
    _assert_png(result, target)
    assert _trace(result, "contour")
    assert result["data"]["cells_plotted"] == int((fragments >= 500).sum())
    assert result["data"]["xaxis_type"] == "log"


def test_frag_size_distr_png_skips_size_zero(tmp_path):
    target = tmp_path / "fragments.png"
    distribution = [0, 5, 12, 18, 9, 3]
    result = module.plot_snapatac2_frag_size_distr_png(distribution, str(target), max_recorded_size=10)
    _assert_png(result, target)
    trace = _trace(result, "scatter")[0]
    assert trace["x"] == [1, 2, 3, 4, 5]
    assert trace["y"] == [5.0, 12.0, 18.0, 9.0, 3.0]
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_frag_size_distr_png([1, -1, 2], str(tmp_path / "neg.png"))


def test_spectral_png_reports_elbow_and_marks_it(tmp_path):
    target = tmp_path / "spectral.png"
    eigenvalues = [10.0, 8.0, 6.0, 1.0, 0.9, 0.85, 0.8]
    result = module.plot_snapatac2_spectral_eigenvalues_png(eigenvalues, str(target))
    _assert_png(result, target)
    assert result["data"]["num_eigen"] == 5
    assert _trace(result, "scatter")[0]["y"] == eigenvalues


def test_spectral_png_without_elbow_reports_none(tmp_path):
    target = tmp_path / "no_elbow.png"
    result = module.plot_snapatac2_spectral_eigenvalues_png([4.0, 2.5, 1.2, 0.4], str(target))
    _assert_png(result, target)
    assert result["data"]["num_eigen"] is None


def test_regions_png_groups_selected_peaks(tmp_path):
    rng = np.random.default_rng(0)
    counts = rng.poisson(2, size=(6, 4)).astype(float)
    peaks = ["chr1:0-5000", "chr1:5000-10000", "chr1:10000-15000", "chr2:0-5000"]
    target = tmp_path / "regions.png"
    result = module.plot_snapatac2_regions_png(
        counts, peaks, ["A", "B"] * 3, {"first": peaks[:2], "second": [peaks[3]]}, str(target), groupby="cell_type",
    )
    _assert_png(result, target)
    heat = _trace(result, "heatmap")[0]
    assert len(heat["x"]) == 2 and len(heat["y"]) == 3
    assert set(heat["x"]) == {"A", "B"}
    assert result["data"]["selected_peaks"] == 3
    assert result["data"]["xaxis_title"] == "cell_type"


def test_regions_png_validates_names_and_selection(tmp_path):
    counts = np.ones((4, 2))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_regions_png(counts, ["chr1-0-5", "chr1:0-5"], list("ABAB"), {"g": ["chr1:0-5"]},
                                          str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_regions_png(counts, ["chr1:0-5", "chr1:5-10"], list("ABAB"), {"g": ["chr9:1-2"]},
                                          str(tmp_path / "b.png"))


def test_umap_png_groups_cells_and_handles_three_dimensions(tmp_path):
    rng = np.random.default_rng(2)
    embedding = rng.normal(size=(20, 2))
    target = tmp_path / "umap.png"
    result = module.plot_snapatac2_umap_png(embedding, str(target), color=["A", "B"] * 10)
    _assert_png(result, target)
    assert sorted(trace["name"] for trace in _trace(result, "scatter")) == ["A", "B"]

    target3 = tmp_path / "umap3.png"
    result3 = module.plot_snapatac2_umap_png(rng.normal(size=(20, 3)), str(target3))
    _assert_png(result3, target3)
    assert _trace(result3, "scatter3d")[0]["name"] == "all cells"
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_umap_png(rng.normal(size=(5, 4)), str(tmp_path / "bad.png"))


def test_motif_enrichment_png_clusters_groups(tmp_path):
    ids = ["MA0001.1", "MA0002.1", "MA0003.1", "MA0004.1"]

    def table(log_fc, adjusted, raw):
        return {"id": ids, "log2(fold change)": log_fc, "adjusted p-value": adjusted, "p-value": raw}

    enrichment = {
        "c1": table([1.5, -0.8, 0.2, 2.0], [0.001, 0.2, 0.5, 0.01], [1e-5, 0.05, 0.3, 1e-3]),
        "c2": table([-1.0, 2.0, 0.1, -1.5], [0.01, 0.001, 0.9, 0.02], [1e-3, 1e-6, 0.8, 1e-4]),
    }
    target = tmp_path / "motif.png"
    result = module.plot_snapatac2_motif_enrichment_png(enrichment, str(target))
    _assert_png(result, target)
    assert result["parameters"] == {"groups": 2, "motifs": 4, "min_log_fc": 1.0, "max_fdr": 0.01}
    assert _trace(result, "heatmap")

    bad_order = dict(enrichment)
    bad_order["c2"] = table([0.0] * 4, [0.5] * 4, [0.5] * 4)
    bad_order["c2"]["id"] = list(reversed(ids))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_motif_enrichment_png(bad_order, str(tmp_path / "bad.png"))
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_motif_enrichment_png({"c1": {"id": ids}, "c2": {"id": ids}}, str(tmp_path / "bad2.png"))


@pytest.fixture()
def fragment_h5ad(tmp_path):
    """A small snapatac2 file with fragments on chr1 and two cell groups (built with snapatac2 itself)."""
    import gzip

    snap = pytest.importorskip("snapatac2")
    rng = np.random.default_rng(3)
    lines = []
    for index in range(400):
        barcode = f"cell{index % 8}"
        start = int(rng.integers(1000, 90000))
        lines.append(f"chr1\t{start}\t{start + 150}\t{barcode}\t1\n")
        lines.append(f"chr1\t{start + 400}\t{start + 600}\t{barcode}\t1\n")
    fragments = tmp_path / "fragments.tsv.gz"
    with gzip.open(fragments, "wt") as handle:
        handle.writelines(lines)
    path = tmp_path / "fragments.h5ad"
    data = snap.pp.import_fragments(str(fragments), chrom_sizes={"chr1": 100000}, file=str(path),
                                    min_num_fragments=1, sorted_by_barcode=False)
    data.obs["group"] = ["A", "B"] * 4
    data.close()
    return path


def test_coverage_png_draws_one_track_per_group(tmp_path, fragment_h5ad):
    target = tmp_path / "coverage.png"
    result = module.plot_snapatac2_coverage_png(str(fragment_h5ad), str(target), region="chr1:1-100000",
                                                groupby="group")
    _assert_png(result, target)
    assert [track["group"] for track in result["data"]["tracks"]] == ["A", "B"]
    assert result["data"]["positions"] == 99999
    assert all(track["max"] > 0 for track in result["data"]["tracks"])


def test_coverage_png_validates_region_and_groupby(tmp_path, fragment_h5ad):
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_coverage_png(str(fragment_h5ad), str(tmp_path / "a.png"), region="chr1:50-10",
                                           groupby="group")
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_coverage_png(str(fragment_h5ad), str(tmp_path / "b.png"), region="chr1:1-100",
                                           groupby="missing")
    with pytest.raises(ConfigurationError):
        module.plot_snapatac2_coverage_png(str(tmp_path / "nope.h5ad"), str(tmp_path / "c.png"),
                                           region="chr1:1-100", groupby="group")


def test_plot_verbose_keyword_is_accepted(tmp_path):
    target = tmp_path / "quiet.png"
    result = module.plot_snapatac2_valid_cells_png([3, 2, 1], str(target), verbose=False)
    _assert_png(result, target)


def test_upstream_type_error_is_reported_as_configuration_error(monkeypatch):
    import snapatac2.datasets as datasets

    def fake_load(dataset, **kwargs):
        raise TypeError("dataset 'x': a `file=` argument is required")

    monkeypatch.setattr(datasets, "load", fake_load)
    with pytest.raises(ConfigurationError) as caught:
        module.load_snapatac2_dataset("atac_pbmc_5k")
    assert "file=" in str(caught.value)
