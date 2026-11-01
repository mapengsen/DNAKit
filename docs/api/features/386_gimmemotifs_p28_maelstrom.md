# GimmeMotifs maelstrom 基序活性预测、排名聚合与热图

这一组函数位于 `dnakit.regulation`，封装 [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs) 0.18.4 中 maelstrom 工具的基序活性部分：`gimmemotifs.maelstrom.moap` 的各个基序活性预测器、`gimmemotifs.maelstrom.rank` 的排名聚合，以及已完成 maelstrom 结果的热图。它们处理预先计算好的表（区域×基序的分数表加标签或数值表）或早先运行的输出目录；不扫描基因组，也不下载任何数据。完整的 `run_maelstrom` 流程和 `MaelstromResult.plot_scores` 的分数箱线图未封装，因为二者在当前上游版本中都会失败（见处置文件）。

GimmeMotifs、scikit-learn、pandas 与 seaborn 是可选依赖（`pip install "dnakit[bio-gimmemotifs]"`）；`xgboost` 方法还需要 `xgboost` 包。示例从 `dnakit.regulation.gimmemotifs_p28_maelstrom_functions` 导入；区域包接入后，同名函数也可从 `dnakit.regulation` 导入。

覆盖 `BayesianRidgeMoap.fit`、`XgboostRegressionMoap.fit`、`MWUMoap.fit`、`RFMoap.fit`、`MultiTaskLassoMoap.fit` 与 `predict`、`SVRMoap.fit` 与 `predict`（均作为同一个函数的方法）、`rank.rankagg` 和 `MaelstromResult.plot_heatmap`。

## 示例数据

示例在临时目录中构造输入：60 个区域、5 个基序（M0 在 `up` 区域富集，M1 在 `down` 区域富集）和一个三类标签表，随机数来自 numpy 旧式 `RandomState(0)`。

## 1) `dnakit.regulation.gimmemotifs_predict_motif_activity()`

- **作用：** 用一种 maelstrom 方法，从区域×基序表预测基序活性。
- **处理方法：** 对预先计算的分数表调用 `gimmemotifs.maelstrom.moap.moap`（不扫描基因组，不做 z 分数或 GC 校正）。`method` 取 `bayesianridge`（贝叶斯岭回归，系数）、`xgboost`（XGBoost 回归，特征重要性，需要 `xgboost` 包）、`mwu`（每个标签类别的 Mann-Whitney U 检验，-log10 Benjamini-Hochberg p 值）、`rf`（随机森林分类，带符号的特征重要性）、`multitasklasso`（多任务 lasso，系数）或 `svr`（线性支持向量回归，权重）。`mwu` 与 `rf` 需要只有一列的标签表；回归方法需要数值列。`predict=True`（仅 `multitasklasso` 与 `svr`）会再次拟合模型并返回其对各区域的预测值。`seed` 使随机方法可重复。
- **API：** `dnakit.regulation.gimmemotifs_predict_motif_activity(labels_file, motif_table, *, method="bayesianridge", predict=False, seed=None)`。
- **输入：** `labels_file`：带表头、以区域为索引的制表符分隔表。`motif_table`：索引相同的制表符分隔的区域×基序表。方法、`predict` 与 `seed`。
- **输出：** `motifs`、`classes` 和 `values`（每个基序一行，每个类别或数值列一列）；`predict=True` 时还有 `predicted_regions`、`predicted_columns` 与 `predicted`。
- **示例代码：**

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

- **示例输出：**

```text
mwu {'down': 'M1', 'flat': 'M2', 'up': 'M0'}
rf {'down': 'M1', 'flat': 'M2', 'up': 'M0'}
bayesianridge {'a': 'M0', 'b': 'M1'}
multitasklasso {'a': 'M0', 'b': 'M1'}
svr {'a': 'M0', 'b': 'M1'}
(60, 2) True
```

## 2) `dnakit.regulation.gimmemotifs_aggregate_motif_ranks()`

- **作用：** 用排名聚合把多列基序活性合并成每个基序一个分数。
- **处理方法：** 调用 `gimmemotifs.maelstrom.rank.rankagg`。`int_stouffer` 把每一列转换成名次，用反正态分布变换名次，再用 Stouffer 方法在各列间合并成每个基序一个 z 分数。`stuart` 是基于 Stuart 等（2003）的 Kolde 等（2012）稳健排名聚合；`include_reverse` 时减去最低名次的聚合结果。并列名次随机打破，因此要得到可重复的值需要 `seed`。基序按分数从高到低返回。
- **API：** `dnakit.regulation.gimmemotifs_aggregate_motif_ranks(scores_table, *, method="int_stouffer", include_reverse=True, seed=None)`。
- **输入：** 带表头的制表符分隔表（行是基序，列是不同方法的活性），或映射 `{motif: {column: value}}`；以及方法、`include_reverse` 与 `seed`。
- **输出：** `ranking`（`motif` 与 `value` 的列表，最好的在前）、`value_name`、`n_motifs` 与 `n_columns`。
- **示例代码：**

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

- **示例输出：**

```text
int_stouffer z-score m0 12 3
stuart score m0 12 3
m0
```

## 3) `dnakit.regulation.gimmemotifs_plot_maelstrom_heatmap()`

- **作用：** 绘制已完成 maelstrom 运行的基序活性聚类热图。
- **处理方法：** 用 `MaelstromResult` 加载目录并调用 `MaelstromResult.plot_heatmap`。目录需含 `final.out.txt`（排名聚合后的活性）、`motif.score.txt.gz`、`motif.count.txt.gz` 和一个 PFM 文件，可选 `motif.freq.txt` 与 `activity.*.txt`。`kind` 取 `final`、`freq`（各簇的基序频率）或某个方法活性的名称（如 `rf.score`）。基序需要活性绝对值至少为 `threshold`、频率至少为 `min_freq`。`name` 为真时行显示因子名而非基序 ID（至多 `max_number_factors` 个；`indirect` 加入预测因子）。行按相关距离的 Ward 聚类排序。PNG 写入新路径，并返回所绘行与数值。
- **API：** `dnakit.regulation.gimmemotifs_plot_maelstrom_heatmap(maelstrom_dir, output_image, *, kind="final", min_freq=0.01, threshold=2.0, name=True, indirect=True, max_number_factors=5, dpi=150)`。
- **输入：** maelstrom 输出目录、新的图像路径和上述选项。
- **输出：** `outputs` 中的 `image`，以及所绘矩阵的 `rows`、`columns` 与 `values`。
- **示例代码：**

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

- **示例输出：**

```text
['M0', 'M2', 'M3', 'M4', 'M5'] ['up', 'down', 'flat']
True
b'\x89PNG'
```
