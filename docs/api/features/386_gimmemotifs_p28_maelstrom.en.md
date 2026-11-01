# GimmeMotifs maelstrom motif-activity predictors, rank aggregation and heatmap

These functions in `dnakit.regulation` wrap the motif-activity part of the maelstrom tool of [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4: the individual motif-activity predictors of `gimmemotifs.maelstrom.moap`, the rank aggregation of `gimmemotifs.maelstrom.rank` and the heatmap of a finished maelstrom result. They work on precomputed tables (a region-by-motif score table with a label or value table) or on the output directory of an earlier run; no genome is scanned and nothing is downloaded. The end-to-end `run_maelstrom` pipeline and the score boxplot of `MaelstromResult.plot_scores` are not wrapped, because both fail in the current upstream release (see the disposition file).

GimmeMotifs, scikit-learn, pandas and seaborn are optional dependencies (`pip install "dnakit[bio-gimmemotifs]"`); the `xgboost` method also needs the `xgboost` package. The examples import from `dnakit.regulation.gimmemotifs_p28_maelstrom_functions`; after the area package is wired the same names are available from `dnakit.regulation`.

They cover `BayesianRidgeMoap.fit`, `XgboostRegressionMoap.fit`, `MWUMoap.fit`, `RFMoap.fit`, `MultiTaskLassoMoap.fit` and `predict`, `SVRMoap.fit` and `predict` (all as methods of one function), `rank.rankagg` and `MaelstromResult.plot_heatmap`.

## Example data

The examples build their inputs in a temporary directory: 60 regions, five motifs (M0 enriched in the `up` regions, M1 in the `down` regions) and a three-class label table, with numpy's legacy `RandomState(0)`.

## 1) `dnakit.regulation.gimmemotifs_predict_motif_activity()`

- **Purpose:** Predict motif activity from a region-by-motif table with one maelstrom method.
- **Method:** Call `gimmemotifs.maelstrom.moap.moap` on a precomputed score table (no genome scan, no z-score or GC correction). `method` is `bayesianridge` (Bayesian ridge regression; coefficients), `xgboost` (XGBoost regression; feature importances; needs the `xgboost` package), `mwu` (Mann-Whitney U test per label class; -log10 Benjamini-Hochberg p-value), `rf` (random forest classifier; signed feature importances), `multitasklasso` (multi-task lasso; coefficients) or `svr` (linear support vector regression; weights). `mwu` and `rf` need a label table with one column; the regression methods need numeric columns. With `predict=True`, available for `multitasklasso` and `svr`, the model is fitted again and the values it predicts for the regions are returned. `seed` makes the stochastic methods reproducible.
- **API:** `dnakit.regulation.gimmemotifs_predict_motif_activity(labels_file, motif_table, *, method="bayesianridge", predict=False, seed=None)`.
- **Input:** `labels_file`: tab-separated table with a header, indexed by region. `motif_table`: tab-separated region-by-motif table with the same index. The method, `predict` and `seed`.
- **Output:** `motifs`, `classes` and `values` (one row per motif, one column per class or value column), and with `predict=True` also `predicted_regions`, `predicted_columns` and `predicted`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_maelstrom_functions import gimmemotifs_predict_motif_activity

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.RandomState(0)
folder = Path(tempfile.mkdtemp())
regions = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(60)]
labels = np.array(["up"] * 20 + ["down"] * 20 + ["flat"] * 20)
scores = pd.DataFrame(rng.normal(size=(60, 5)), index=regions, columns=[f"M{i}" for i in range(5)])
scores.loc[labels == "up", "M0"] += 2
scores.loc[labels == "down", "M1"] += 2
scores.to_csv(folder / "scores.tsv", sep="\t")
pd.DataFrame({"cluster": labels}, index=regions).to_csv(folder / "labels.tsv", sep="\t")
values = pd.DataFrame({"a": (labels == "up") + rng.normal(size=60) * 0.1,
                       "b": (labels == "down") + rng.normal(size=60) * 0.1}, index=regions)
values.to_csv(folder / "values.tsv", sep="\t")

for method in ("mwu", "rf"):
    result = gimmemotifs_predict_motif_activity(folder / "labels.tsv", folder / "scores.tsv", method=method, seed=1)
    column = {name: [row[i] for row in result["values"]] for i, name in enumerate(result["classes"])}
    top = {name: result["motifs"][int(np.argmax(values))] for name, values in sorted(column.items())}
    print(method, top)
for method in ("bayesianridge", "multitasklasso", "svr"):
    result = gimmemotifs_predict_motif_activity(folder / "values.tsv", folder / "scores.tsv", method=method, seed=1)
    column = {name: [row[i] for row in result["values"]] for i, name in enumerate(result["classes"])}
    print(method, {name: result["motifs"][int(np.argmax(np.abs(values)))] for name, values in column.items()})
fitted = gimmemotifs_predict_motif_activity(folder / "values.tsv", folder / "scores.tsv", method="svr", predict=True, seed=1)
predicted = np.array(fitted["predicted"])
print(predicted.shape, round(float(np.corrcoef(predicted[:, 0], values["a"])[0, 1]), 1) > 0.5)
```

- **Example output:**

```text
mwu {'down': 'M1', 'flat': 'M2', 'up': 'M0'}
rf {'down': 'M1', 'flat': 'M2', 'up': 'M0'}
bayesianridge {'a': 'M0', 'b': 'M1'}
multitasklasso {'a': 'M0', 'b': 'M1'}
svr {'a': 'M0', 'b': 'M1'}
(60, 2) True
```

## 2) `dnakit.regulation.gimmemotifs_aggregate_motif_ranks()`

- **Purpose:** Aggregate several motif-activity columns into one score per motif by rank aggregation.
- **Method:** Call `gimmemotifs.maelstrom.rank.rankagg`. `int_stouffer` converts each column to ranks, transforms the ranks with the inverse normal distribution and combines them over the columns with Stouffer's method into a z-score per motif. `stuart` is the robust rank aggregation of Kolde et al. (2012) based on Stuart et al. (2003); with `include_reverse` the aggregation of the lowest ranks is subtracted. Ties are broken at random, so `seed` is needed for reproducible values. The motifs are returned best first.
- **API:** `dnakit.regulation.gimmemotifs_aggregate_motif_ranks(scores_table, *, method="int_stouffer", include_reverse=True, seed=None)`.
- **Input:** A tab-separated table with a header (rows are motifs, columns are the activities of different methods) or a mapping `{motif: {column: value}}`; the method, `include_reverse` and `seed`.
- **Output:** `ranking` (list of `motif` and `value`, best first), `value_name`, `n_motifs` and `n_columns`.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_maelstrom_functions import gimmemotifs_aggregate_motif_ranks

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.RandomState(0)
folder = Path(tempfile.mkdtemp())
regions = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(60)]
labels = np.array(["up"] * 20 + ["down"] * 20 + ["flat"] * 20)
scores = pd.DataFrame(rng.normal(size=(60, 5)), index=regions, columns=[f"M{i}" for i in range(5)])
scores.loc[labels == "up", "M0"] += 2
scores.loc[labels == "down", "M1"] += 2
scores.to_csv(folder / "scores.tsv", sep="\t")
pd.DataFrame({"cluster": labels}, index=regions).to_csv(folder / "labels.tsv", sep="\t")
values = pd.DataFrame({"a": (labels == "up") + rng.normal(size=60) * 0.1,
                       "b": (labels == "down") + rng.normal(size=60) * 0.1}, index=regions)
values.to_csv(folder / "values.tsv", sep="\t")

table = pd.DataFrame(rng.normal(size=(12, 3)), index=[f"m{i}" for i in range(12)], columns=["x", "y", "z"])
table.iloc[0] += 4
table.to_csv(folder / "table.tsv", sep="\t")
for method in ("int_stouffer", "stuart"):
    result = gimmemotifs_aggregate_motif_ranks(folder / "table.tsv", method=method, seed=7)
    print(method, result["value_name"], result["ranking"][0]["motif"], result["n_motifs"], result["n_columns"])
mapping = {name: row.to_dict() for name, row in table.iterrows()}
print(gimmemotifs_aggregate_motif_ranks(mapping, seed=7)["ranking"][0]["motif"])
```

- **Example output:**

```text
int_stouffer z-score m0 12 3
stuart score m0 12 3
m0
```

## 3) `dnakit.regulation.gimmemotifs_plot_maelstrom_heatmap()`

- **Purpose:** Plot the clustered heatmap of motif activity of a finished maelstrom run.
- **Method:** Load the directory with `MaelstromResult` and call `MaelstromResult.plot_heatmap`. The directory holds `final.out.txt` (the rank-aggregated activity), `motif.score.txt.gz`, `motif.count.txt.gz` and a PFM file, and optionally `motif.freq.txt` and `activity.*.txt`. `kind` is `final`, `freq` (motif frequencies per cluster) or the name of one method's activity such as `rf.score`. Motifs need an absolute activity of at least `threshold` and a frequency of at least `min_freq`. With `name` the rows carry factor names instead of motif IDs (at most `max_number_factors`; `indirect` adds predicted factors). The rows are ordered by Ward clustering with correlation distance. The PNG is written to a new path and the plotted rows and values are returned.
- **API:** `dnakit.regulation.gimmemotifs_plot_maelstrom_heatmap(maelstrom_dir, output_image, *, kind="final", min_freq=0.01, threshold=2.0, name=True, indirect=True, max_number_factors=5, dpi=150)`.
- **Input:** A maelstrom output directory, a new image path and the options above.
- **Output:** `outputs` with the `image`, and `rows`, `columns` and `values` of the plotted matrix.
- **Example code:**

```python
from dnakit.regulation.gimmemotifs_p28_maelstrom_functions import gimmemotifs_plot_maelstrom_heatmap

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.RandomState(0)
folder = Path(tempfile.mkdtemp())
regions = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(60)]
labels = np.array(["up"] * 20 + ["down"] * 20 + ["flat"] * 20)
scores = pd.DataFrame(rng.normal(size=(60, 5)), index=regions, columns=[f"M{i}" for i in range(5)])
scores.loc[labels == "up", "M0"] += 2
scores.loc[labels == "down", "M1"] += 2
scores.to_csv(folder / "scores.tsv", sep="\t")
pd.DataFrame({"cluster": labels}, index=regions).to_csv(folder / "labels.tsv", sep="\t")
values = pd.DataFrame({"a": (labels == "up") + rng.normal(size=60) * 0.1,
                       "b": (labels == "down") + rng.normal(size=60) * 0.1}, index=regions)
values.to_csv(folder / "values.tsv", sep="\t")

run = folder / "run"
run.mkdir()
motifs = [f"M{i}" for i in range(6)]
final = pd.DataFrame(rng.normal(scale=3, size=(6, 3)), index=pd.Index(motifs, name="motif"), columns=["up", "down", "flat"])
final.to_csv(run / "final.out.txt", sep="\t")
big = [f"chr1:{i * 100}-{i * 100 + 100}" for i in range(30)]
pd.DataFrame(rng.normal(size=(30, 6)), index=big, columns=motifs).to_csv(run / "motif.score.txt.gz", sep="\t")
pd.DataFrame(rng.randint(0, 3, size=(30, 6)), index=big, columns=motifs).to_csv(run / "motif.count.txt.gz", sep="\t")
(run / "nonredundant.pfm").write_text("".join(f">{m}\n10 0 0 0\n0 10 0 0\n0 0 10 0\n" for m in motifs))

result = gimmemotifs_plot_maelstrom_heatmap(run, folder / "heatmap.png", name=False)
print(sorted(result["rows"]), result["columns"])
print(sorted(final.index[(final.abs() >= 2).any(axis=1)]) == sorted(result["rows"]))
print((folder / "heatmap.png").read_bytes()[:4])
```

- **Example output:**

```text
['M0', 'M2', 'M3', 'M4', 'M5'] ['up', 'down', 'flat']
True
b'\x89PNG'
```
