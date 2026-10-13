"""Contract tests for AnnData, scanpy, muon and scvi-tools transcriptomics backends."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.transcriptomics.backends import (
    aggregate_pseudobulk,
    calculate_cell_qc_metrics,
    calculate_gene_qc_metrics,
    filter_genes_by_min_cells,
    rank_marker_genes,
    read_anndata_summary,
    select_highly_variable_genes,
    summarize_mudata,
    train_scvi_model,
)

# Four cells, four genes. Gene names starting with MT- are mitochondrial.
COUNTS = [[3, 0, 1, 0], [0, 2, 0, 0], [1, 1, 1, 0], [0, 0, 5, 2]]
GENES = ["MT-A", "GENE1", "GENE2", "RPL3"]
CELLS = ["c0", "c1", "c2", "c3"]


def _write_tiny(tmp_path: Path, name: str = "tiny.h5ad") -> Path:
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    obs = pd.DataFrame(
        {"group": ["x", "x", "y", "y"], "sample": ["s1", "s2", "s1", "s2"]}, index=CELLS
    )
    var = pd.DataFrame({"chromosome": ["chrM", "chr1", "chr2", "chr3"]}, index=GENES)
    adata = anndata.AnnData(X=np.array(COUNTS, dtype=np.float32), obs=obs, var=var)
    path = tmp_path / name
    adata.write_h5ad(path)
    return path


def _write_markers(tmp_path: Path) -> Path:
    """Six cells in two groups; gene g0 is high in group A and low in group B."""

    anndata = pytest.importorskip("anndata")
    import pandas as pd

    counts = np.array(
        [[9, 1, 2, 0], [8, 2, 1, 0], [7, 3, 2, 1], [1, 2, 3, 1], [2, 2, 1, 1], [0, 4, 2, 1]],
        dtype=np.float32,
    )
    obs = pd.DataFrame({"ct": ["A", "A", "A", "B", "B", "B"]}, index=[f"c{i}" for i in range(6)])
    path = tmp_path / "markers.h5ad"
    anndata.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=["g0", "g1", "g2", "g3"])).write_h5ad(
        path
    )
    return path


def _write_random(tmp_path: Path, name: str = "random.h5ad") -> Path:
    """Forty cells, thirty genes, seeded; the first 20 cells over-express five genes."""

    anndata = pytest.importorskip("anndata")
    import pandas as pd

    rng = np.random.default_rng(7)
    counts = rng.poisson(3, size=(40, 30)).astype(np.float32)
    counts[:20, :5] += 6
    obs = pd.DataFrame({"ct": ["A"] * 20 + ["B"] * 20}, index=[f"c{i}" for i in range(40)])
    path = tmp_path / name
    anndata.AnnData(
        X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(30)])
    ).write_h5ad(path)
    return path


def test_read_anndata_summary_reports_dimensions_and_gene_metadata(tmp_path: Path) -> None:
    anndata = pytest.importorskip("anndata")
    path = _write_tiny(tmp_path)

    result = read_anndata_summary(path, verbose=False)

    assert result["status"] == "completed"
    assert (result["n_obs"], result["n_vars"]) == anndata.read_h5ad(path).shape
    assert result["obs_columns"] == ["group", "sample"]
    assert result["var_columns"] == ["chromosome"]
    assert result["var_preview"][0] == {"gene": "MT-A", "chromosome": "chrM"}
    assert result["x_dtype"] == "float32"
    assert result["layers"] == []
    assert result["provenance"]["package"].startswith("anndata ")


def test_gene_qc_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    path = _write_tiny(tmp_path)

    result = calculate_gene_qc_metrics(path, top_n=3, verbose=False)

    # Totals: MT-A 4, GENE1 3, GENE2 7, RPL3 2. Cells with a count: 2, 2, 3, 1.
    assert result["n_genes"] == 4
    assert result["genes_detected"] == 4
    assert result["median_cells_per_gene"] == pytest.approx(2.0)
    assert [gene["gene"] for gene in result["top_genes"]] == ["GENE2", "MT-A", "GENE1"]
    gene2, mt_a = result["top_genes"][0], result["top_genes"][1]
    assert gene2["total_counts"] == pytest.approx(7.0)
    assert gene2["n_cells_by_counts"] == 3
    assert gene2["mean_counts"] == pytest.approx(1.75)
    assert gene2["pct_dropout_by_counts"] == pytest.approx(25.0)
    assert mt_a["mean_counts"] == pytest.approx(1.0)
    assert mt_a["pct_dropout_by_counts"] == pytest.approx(50.0)


def test_cell_qc_reports_mitochondrial_fraction(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    path = _write_tiny(tmp_path)

    result = calculate_cell_qc_metrics(path, verbose=False)

    # Cell mito fractions: 3/4, 0/2, 1/3, 0/7 -> 75, 0, 33.33 and 0 percent.
    assert result["mito_genes"] == 1
    assert result["median_genes_per_cell"] == pytest.approx(2.0)  # genes 2, 1, 3, 2
    assert result["median_total_counts"] == pytest.approx(3.5)  # totals 4, 2, 3, 7
    assert result["median_pct_counts_mt"] == pytest.approx(100 / 6)  # median of 0, 0, 33.3, 75
    assert result["max_pct_counts_mt"] == pytest.approx(75.0)
    assert result["cells_above_mito_threshold"] == 2  # 75 and 33.3 exceed 20 percent

    no_prefix = calculate_cell_qc_metrics(path, mito_prefix="mt-", verbose=False)
    assert no_prefix["mito_genes"] == 0
    assert no_prefix["max_pct_counts_mt"] == pytest.approx(0.0)


def test_filter_genes_by_min_cells_writes_subset(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")
    path = _write_tiny(tmp_path)
    output = tmp_path / "filtered.h5ad"

    result = filter_genes_by_min_cells(path, min_cells=2, output_path=output, verbose=False)

    # RPL3 is detected in one cell only, so it is removed.
    assert result["n_genes_before"] == 4
    assert result["n_genes_after"] == 3
    assert result["n_genes_removed"] == 1
    assert anndata.read_h5ad(output).var_names.tolist() == ["MT-A", "GENE1", "GENE2"]
    with pytest.raises(ConfigurationError, match="already exists"):
        filter_genes_by_min_cells(path, min_cells=2, output_path=output, verbose=False)


def test_highly_variable_genes_match_direct_scanpy(tmp_path: Path) -> None:
    sc = pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")
    path = _write_random(tmp_path)

    for flavor in ("seurat", "cell_ranger"):
        result = select_highly_variable_genes(path, n_top_genes=5, flavor=flavor, verbose=False)
        direct = anndata.read_h5ad(path)
        sc.pp.normalize_total(direct, target_sum=1e4)
        sc.pp.log1p(direct)
        table = sc.pp.highly_variable_genes(direct, n_top_genes=5, flavor=flavor, inplace=False)
        selected = table[table["highly_variable"]]
        expected = sorted(
            selected.index.astype(str),
            key=lambda gene: (-selected.loc[gene, "dispersions_norm"], gene),
        )
        assert result["highly_variable_genes"] == expected
        # scanpy flags every gene at the cutoff, so tied dispersions can add genes beyond n_top.
        assert result["n_highly_variable"] == len(expected) >= 5

    with pytest.raises(ConfigurationError, match="exceeds"):
        select_highly_variable_genes(path, n_top_genes=31, verbose=False)
    with pytest.raises(ConfigurationError, match="flavor must be one of"):
        select_highly_variable_genes(path, flavor="seurat_v3", verbose=False)


def test_rank_marker_genes_wilcoxon_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    path = _write_markers(tmp_path)

    result = rank_marker_genes(path, groupby="ct", n_genes=4, log_normalize=False, verbose=False)

    # Gene g0 in group A: all three A values exceed all three B values, so U = 9 of n1*n2 = 9.
    # Normal approximation: z = (U - n1*n2/2) / sqrt(n1*n2*(n1+n2+1)/12).
    z = (9 - 9 / 2) / math.sqrt(3 * 3 * 7 / 12)
    p_value = math.erfc(z / math.sqrt(2))
    top_a = result["groups"]["A"][0]
    assert top_a["gene"] == "g0"
    assert top_a["score"] == pytest.approx(z, rel=1e-6)
    # Benjamini-Hochberg with four genes: the top-ranked p-value is multiplied by 4.
    assert top_a["pval_adj"] == pytest.approx(4 * p_value, rel=1e-5)
    # log2 fold change of the expm1 means: mean A = 8, mean B = (1 + 2 + 0) / 3 = 1.
    assert top_a["logfoldchange"] == pytest.approx(math.log2(math.expm1(8) / math.expm1(1)))
    # Group B mirrors group A: scores flip sign for the same gene.
    by_name_a = {row["gene"]: row["score"] for row in result["groups"]["A"]}
    by_name_b = {row["gene"]: row["score"] for row in result["groups"]["B"]}
    assert by_name_b["g3"] == pytest.approx(-by_name_a["g3"])
    assert list(result["groups"]) == ["A", "B"]


def test_rank_marker_genes_matches_direct_scanpy_after_normalization(tmp_path: Path) -> None:
    sc = pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")
    path = _write_random(tmp_path)

    result = rank_marker_genes(path, groupby="ct", n_genes=3, verbose=False)

    direct = anndata.read_h5ad(path)
    sc.pp.normalize_total(direct, target_sum=1e4)
    sc.pp.log1p(direct)
    sc.tl.rank_genes_groups(
        direct,
        groupby="ct",
        method="wilcoxon",
        n_genes=3,
        use_raw=False,
        corr_method="benjamini-hochberg",
        key_added="rank_genes_groups",
    )
    table = sc.get.rank_genes_groups_df(direct, group="A")
    assert [row["gene"] for row in result["groups"]["A"]] == table["names"].tolist()
    assert [row["score"] for row in result["groups"]["A"]] == pytest.approx(
        table["scores"].tolist()
    )


def test_rank_marker_genes_rejects_bad_groupby(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    path = _write_markers(tmp_path)
    with pytest.raises(ConfigurationError, match="not an obs column"):
        rank_marker_genes(path, groupby="missing", log_normalize=False, verbose=False)
    with pytest.raises(ConfigurationError, match="method must be one of"):
        rank_marker_genes(path, groupby="ct", method="logreg", verbose=False)


def test_aggregate_pseudobulk_sums_match_hand_totals(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    counts = np.array(
        [[1, 0, 2], [3, 1, 0], [0, 4, 1], [2, 2, 2], [5, 0, 0], [1, 1, 1]], dtype=np.float32
    )
    obs = pd.DataFrame(
        {
            "group": ["A", "A", "A", "B", "B", "B"],
            "sample": ["s1", "s1", "s2", "s1", "s2", "s2"],
        },
        index=[f"c{i}" for i in range(6)],
    )
    path = tmp_path / "pseudo.h5ad"
    anndata.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=["g0", "g1", "g2"])).write_h5ad(path)
    output = tmp_path / "pseudobulk.h5ad"

    result = aggregate_pseudobulk(
        path,
        group_key="group",
        sample_key="sample",
        min_cells_per_sample=2,
        output_path=output,
        verbose=False,
    )

    # A/s1 = [4, 1, 2] (2 cells), A/s2 = [0, 4, 1] (1 cell),
    # B/s1 = [2, 2, 2] (1 cell), B/s2 = [6, 1, 1] (2 cells).
    kept = {(row["group"], row["sample"]): row for row in result["samples"]}
    assert set(kept) == {("A", "s1"), ("B", "s2")}
    assert kept[("A", "s1")]["n_cells"] == 2
    assert kept[("A", "s1")]["total_counts"] == pytest.approx(7.0)
    assert kept[("B", "s2")]["total_counts"] == pytest.approx(8.0)
    assert result["dropped_samples"] == 2
    written = anndata.read_h5ad(output)
    assert written.layers["sum"].tolist() == [[4, 1, 2], [6, 1, 1]]
    assert written.obs["n_cells"].tolist() == [2, 2]
    with pytest.raises(ConfigurationError, match="must be different"):
        aggregate_pseudobulk(path, group_key="group", sample_key="group", verbose=False)


def test_summarize_mudata_counts_union_and_shared_cells(tmp_path: Path) -> None:
    anndata = pytest.importorskip("anndata")
    muon = pytest.importorskip("muon")
    import pandas as pd

    rng = np.random.default_rng(0)
    rna = anndata.AnnData(
        X=rng.poisson(1, size=(5, 3)).astype(np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    atac = anndata.AnnData(
        X=rng.poisson(1, size=(4, 6)).astype(np.float32),
        obs=pd.DataFrame(index=["c0", "c1", "c2", "c9"]),
        var=pd.DataFrame(index=[f"p{i}" for i in range(6)]),
    )
    path = tmp_path / "multi.h5mu"
    muon.write_h5mu(path, muon.MuData({"rna": rna, "atac": atac}))

    result = summarize_mudata(path, verbose=False)

    # Cells: rna has c0-c4, atac has c0, c1, c2, c9. Union 6, shared by both 3.
    assert result["n_obs"] == 6
    assert result["n_obs_shared"] == 3
    assert result["modalities"] == {
        "atac": {"n_obs": 4, "n_vars": 6},
        "rna": {"n_obs": 5, "n_vars": 3},
    }


def test_train_scvi_model_one_epoch_is_reproducible(tmp_path: Path) -> None:
    pytest.importorskip("scvi")
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    rng = np.random.default_rng(0)
    counts = rng.poisson(1.5, size=(120, 60)).astype(np.float32)
    obs = pd.DataFrame({"batch": ["b1", "b2"] * 60}, index=[f"c{i}" for i in range(120)])
    path = tmp_path / "scvi.h5ad"
    anndata.AnnData(
        X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(60)])
    ).write_h5ad(path)
    options: dict[str, Any] = {"batch_key": "batch", "n_latent": 5, "n_hidden": 16, "seed": 0}

    first = train_scvi_model(path, verbose=False, **options)
    second = train_scvi_model(path, verbose=False, **options)

    assert first["latent_shape"] == [120, 5]
    assert first["epochs_trained"] == 1
    assert first["accelerator"] == "cpu"
    assert first["latent_finite"] is True
    assert math.isfinite(first["final_train_elbo"])
    assert first["latent_dimension_means"] == second["latent_dimension_means"]
    assert first["final_train_elbo"] == second["final_train_elbo"]


def test_train_scvi_model_requires_counts_and_valid_options(tmp_path: Path) -> None:
    pytest.importorskip("scvi")
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    path = tmp_path / "normalized.h5ad"
    anndata.AnnData(
        X=np.full((4, 3), 0.5, dtype=np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["a", "b", "c"]),
    ).write_h5ad(path)
    with pytest.raises(ConfigurationError, match="integer counts"):
        train_scvi_model(path, verbose=False)
    with pytest.raises(ConfigurationError, match="accelerator must be one of"):
        train_scvi_model(path, accelerator="tpu", verbose=False)
    with pytest.raises(ConfigurationError, match="max_epochs"):
        train_scvi_model(path, max_epochs=0, verbose=False)


def test_bounds_and_input_validation_reject_bad_arguments(tmp_path: Path) -> None:
    pytest.importorskip("scanpy")
    path = _write_tiny(tmp_path)
    with pytest.raises(ConfigurationError, match="top_n must be an integer"):
        calculate_gene_qc_metrics(path, top_n=True, verbose=False)
    with pytest.raises(ConfigurationError, match="max_matrix_elements"):
        calculate_gene_qc_metrics(path, max_matrix_elements=15, verbose=False)
    with pytest.raises(ConfigurationError, match="max_file_bytes"):
        read_anndata_summary(path, max_file_bytes=10, verbose=False)
    with pytest.raises(ConfigurationError, match="NUL"):
        read_anndata_summary("bad\x00name.h5ad", verbose=False)
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("cell,gene\nc0,1\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match=r"must be a \.h5ad file"):
        read_anndata_summary(csv_path, verbose=False)
    with pytest.raises(ConfigurationError, match="min_cells"):
        filter_genes_by_min_cells(path, min_cells=0, verbose=False)
    with pytest.raises(ConfigurationError, match="must be a boolean"):
        select_highly_variable_genes(path, log_normalize="yes", verbose=False)


def test_verbose_flag_controls_dnakit_messages(tmp_path: Path, capsys: Any) -> None:
    pytest.importorskip("anndata")
    path = _write_tiny(tmp_path)
    read_anndata_summary(path, verbose=False)
    assert capsys.readouterr().err == ""
    read_anndata_summary(path)
    assert "Starting" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("module", "extra", "call"),
    [
        ("scanpy", "bio-scanpy", lambda path: calculate_gene_qc_metrics(path, verbose=False)),
        ("anndata", "bio-scanpy", lambda path: read_anndata_summary(path, verbose=False)),
        ("scvi", "bio-scvi", lambda path: train_scvi_model(path, verbose=False)),
        ("muon", "bio-muon", lambda path: summarize_mudata(path, verbose=False)),
    ],
)
def test_missing_optional_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, module: str, extra: str, call: Any
) -> None:
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.raises(BackendUnavailableError, match=extra):
        call(tmp_path / "missing.h5ad")
