from __future__ import annotations

import importlib.util

import pytest

pytest.importorskip("gimmemotifs")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
pytest.importorskip("sklearn")
pytest.importorskip("seaborn")

from gimmemotifs.maelstrom.moap import moap  # noqa: E402
from gimmemotifs.maelstrom.rank import rankagg  # noqa: E402

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.regulation import gimmemotifs_p28_maelstrom_functions as module  # noqa: E402


@pytest.fixture
def tables(tmp_path):
    rng = np.random.RandomState(0)
    regions = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(60)]
    labels = np.array(["up"] * 20 + ["down"] * 20 + ["flat"] * 20)
    scores = pd.DataFrame(
        rng.normal(size=(60, 5)), index=regions, columns=[f"M{i}" for i in range(5)]
    )
    scores.loc[labels == "up", "M0"] += 2
    scores.loc[labels == "down", "M1"] += 2
    scores.to_csv(tmp_path / "scores.tsv", sep="\t")
    pd.DataFrame({"cluster": labels}, index=regions).to_csv(tmp_path / "labels.tsv", sep="\t")
    values = pd.DataFrame(
        {
            "a": (labels == "up") + rng.normal(size=60) * 0.1,
            "b": (labels == "down") + rng.normal(size=60) * 0.1,
        },
        index=regions,
    )
    values.to_csv(tmp_path / "values.tsv", sep="\t")
    return tmp_path


@pytest.mark.parametrize("method", ["mwu", "rf"])
def test_classification_methods_match_upstream(tables, method):
    result = module.gimmemotifs_predict_motif_activity(
        tables / "labels.tsv", tables / "scores.tsv", method=method, seed=1, verbose=False
    )
    expected = moap(
        str(tables / "labels.tsv"),
        method=method,
        scoring="score",
        motiffile=str(tables / "scores.tsv"),
        ncpus=1,
        zscore=False,
        gc=False,
        random_state=np.random.RandomState(1),
    )
    assert result["classes"] == [str(c) for c in expected.columns]
    assert np.allclose(np.array(result["values"]), expected.to_numpy())
    best = {name: row for name, row in zip(result["motifs"], result["values"], strict=True)}
    up = result["classes"].index("up")
    assert best["M0"][up] == max(row[up] for row in best.values())


@pytest.mark.parametrize("method", ["bayesianridge", "multitasklasso", "svr"])
def test_regression_methods_match_upstream(tables, method):
    result = module.gimmemotifs_predict_motif_activity(
        tables / "values.tsv", tables / "scores.tsv", method=method, seed=1, verbose=False
    )
    expected = moap(
        str(tables / "values.tsv"),
        method=method,
        scoring="score",
        motiffile=str(tables / "scores.tsv"),
        ncpus=1,
        zscore=False,
        gc=False,
        random_state=np.random.RandomState(1),
    )
    assert result["classes"] == ["a", "b"]
    assert np.allclose(np.array(result["values"]), expected.to_numpy())
    assert abs(result["values"][0][0]) > abs(result["values"][4][0])  # M0 drives value a


@pytest.mark.parametrize("method", ["multitasklasso", "svr"])
def test_predict_option(tables, method):
    result = module.gimmemotifs_predict_motif_activity(
        tables / "values.tsv",
        tables / "scores.tsv",
        method=method,
        predict=True,
        seed=1,
        verbose=False,
    )
    predicted = np.array(result["predicted"])
    assert predicted.shape == (60, 2) and result["predicted_columns"] == ["a", "b"]
    truth = pd.read_table(tables / "values.tsv", index_col=0).to_numpy()
    assert np.corrcoef(predicted[:, 0], truth[:, 0])[0, 1] > 0.5
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_predict_motif_activity(
            tables / "values.tsv", tables / "scores.tsv", method="mwu", predict=True, verbose=False
        )


def test_xgboost_needs_package(tables):
    if importlib.util.find_spec("xgboost"):
        pytest.skip("xgboost is installed")
    with pytest.raises(BackendUnavailableError):
        module.gimmemotifs_predict_motif_activity(
            tables / "values.tsv", tables / "scores.tsv", method="xgboost", verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_predict_motif_activity(
            tables / "values.tsv", tables / "scores.tsv", method="hypergeom", verbose=False
        )


@pytest.mark.parametrize("method", ["int_stouffer", "stuart"])
def test_aggregate_ranks_match_upstream(tmp_path, method):
    rng = np.random.RandomState(1)
    table = pd.DataFrame(
        rng.normal(size=(12, 3)), index=[f"m{i}" for i in range(12)], columns=["x", "y", "z"]
    )
    table.iloc[0] += 4  # m0 is best in every column
    table.to_csv(tmp_path / "t.tsv", sep="\t")
    result = module.gimmemotifs_aggregate_motif_ranks(
        tmp_path / "t.tsv", method=method, seed=7, verbose=False
    )
    np.random.seed(7)
    import random

    random.seed(7)
    expected = rankagg(table, method=method).iloc[:, 0]
    got = {row["motif"]: row["value"] for row in result["ranking"]}
    assert all(got[name] == pytest.approx(float(value)) for name, value in expected.items())
    assert result["ranking"][0]["motif"] == "m0"
    from_mapping = module.gimmemotifs_aggregate_motif_ranks(
        {name: row.to_dict() for name, row in table.iterrows()},
        method=method,
        seed=7,
        verbose=False,
    )
    assert from_mapping["ranking"] == result["ranking"]
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_aggregate_motif_ranks(tmp_path / "t.tsv", method="mean", verbose=False)


def _maelstrom_dir(path):
    rng = np.random.RandomState(0)
    motifs = [f"M{i}" for i in range(6)]
    pd.DataFrame(
        rng.normal(scale=3, size=(6, 3)),
        index=pd.Index(motifs, name="motif"),
        columns=["up", "down", "flat"],
    ).to_csv(path / "final.out.txt", sep="\t")
    regions = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(30)]
    pd.DataFrame(rng.normal(size=(30, 6)), index=regions, columns=motifs).to_csv(
        path / "motif.score.txt.gz", sep="\t"
    )
    pd.DataFrame(rng.randint(0, 3, size=(30, 6)), index=regions, columns=motifs).to_csv(
        path / "motif.count.txt.gz", sep="\t"
    )
    (path / "nonredundant.pfm").write_text(
        "".join(f">{m}\n10 0 0 0\n0 10 0 0\n0 0 10 0\n" for m in motifs)
    )
    return path


def test_plot_maelstrom_heatmap(tmp_path):
    pytest.importorskip("matplotlib")
    directory = _maelstrom_dir(tmp_path)
    result = module.gimmemotifs_plot_maelstrom_heatmap(
        directory, tmp_path / "h.png", name=False, verbose=False
    )
    final = pd.read_table(directory / "final.out.txt", index_col=0)
    kept = final.loc[(final.abs() >= 2).any(axis=1)]
    assert sorted(result["rows"]) == sorted(kept.index)
    assert result["columns"] == ["up", "down", "flat"]
    shown = pd.DataFrame(result["values"], index=result["rows"], columns=result["columns"])
    assert np.allclose(shown.loc[kept.index].to_numpy(), kept.to_numpy())
    assert (tmp_path / "h.png").read_bytes()[:4] == b"\x89PNG"
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_plot_maelstrom_heatmap(directory, tmp_path / "h.png", verbose=False)
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_plot_maelstrom_heatmap(
            directory, tmp_path / "h3.png", threshold=100, verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.gimmemotifs_plot_maelstrom_heatmap(
            directory, tmp_path / "h4.png", kind="rf.score", verbose=False
        )
