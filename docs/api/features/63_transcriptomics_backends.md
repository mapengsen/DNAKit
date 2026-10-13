# 转录组单细胞后端：AnnData、scanpy、muon 与 scVI

这一组函数是 `dnakit.transcriptomics` 中基于可选依赖的基因水平单细胞适配器。[anndata](https://github.com/scverse/anndata) 读取和写出 `.h5ad`，[scanpy](https://github.com/scverse/scanpy) 计算细胞和基因的 QC、高变基因、标记基因与伪批量求和，[muon](https://github.com/scverse/muon) 汇总 `.h5mu` 多模态文件，[scvi-tools](https://github.com/scverse/scvi-tools) 设置并训练一个 scVI 模型。这些函数不做聚类、降维、细胞类型注释和轨迹分析，也不下载公共数据，只处理显式给出的文件。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的示例按顺序运行，共用同一个工作目录。第一个示例把 `counts.csv` 写成 `tiny.h5ad`，之后的示例复用它；其余示例会自己生成需要的 `.h5ad` 或 `.h5mu` 文件。数据只用于演示接口，统计值不代表真实细胞。

**示例文件：** `counts.csv`

```text
cell,group,sample,MT-A,GENE1,GENE2,RPL3
c0,x,s1,3,0,1,0
c1,x,s2,0,2,0,0
c2,y,s1,1,1,1,0
c3,y,s2,0,0,5,2
```

## 1) `dnakit.transcriptomics.read_anndata_summary()`

- **作用：** 读取 `.h5ad` 的维度、细胞和基因元数据列名、前几个基因的元数据，以及附加数据的键，不加载计数矩阵。
- **处理方法：** 以 backed 模式打开文件，只读取 `obs`、`var`、`layers`、`obsm`、`varm` 和 `uns` 的结构，矩阵只记录 dtype。文件大小超过 `max_file_bytes` 时直接报错。
- **API：** `dnakit.transcriptomics.read_anndata_summary(source, *, max_file_bytes=1073741824)`。
- **输入：** `.h5ad` 文件路径；`max_file_bytes` 是允许读取的文件大小上限（字节）。
- **输出：** 含 `n_obs`、`n_vars`、`x_dtype`、`obs_columns`、`var_columns`、`var_preview`、`layers`、`obsm`、`varm`、`uns` 和 `provenance` 的字典。`var_preview` 最多列出 5 个基因及其元数据。
- **示例代码：**

```python
import anndata as ad
import pandas as pd
from dnakit.transcriptomics import read_anndata_summary

table = pd.read_csv("counts.csv", index_col="cell")
var = pd.DataFrame({"chromosome": ["chrM", "chr1", "chr2", "chr3"]}, index=table.columns[2:])
adata = ad.AnnData(X=table.iloc[:, 2:].to_numpy(dtype="float32"), obs=table[["group", "sample"]], var=var)
adata.write_h5ad("tiny.h5ad")

result = read_anndata_summary("tiny.h5ad")
print(result["n_obs"], result["n_vars"], result["obs_columns"])
print(result["var_preview"][0])
```

- **示例输出：**

```text
4 4 ['group', 'sample']
{'gene': 'MT-A', 'chromosome': 'chrM'}
```

## 2) `dnakit.transcriptomics.calculate_gene_qc_metrics()`

- **作用：** 计算每个基因的检测细胞数、总计数、平均计数和丢失率，并列出总计数最高的基因。
- **处理方法：** 调用 `scanpy.pp.calculate_qc_metrics` 得到基因表。`n_cells_by_counts` 是计数大于 0 的细胞数，`pct_dropout_by_counts` 是未检测到该基因的细胞百分比。总计数相同的基因按名称排序。
- **API：** `dnakit.transcriptomics.calculate_gene_qc_metrics(source, *, top_n=10, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`top_n` 是列出的基因数（1 到 1000）；`max_matrix_elements` 限制细胞数乘基因数，超过时报错。
- **输出：** `n_genes`、`genes_detected`、`median_cells_per_gene`、`top_genes`（每个基因的总计数、检测细胞数、平均计数和丢失率）和 `provenance`。
- **示例代码：**

```python
from dnakit.transcriptomics import calculate_gene_qc_metrics

result = calculate_gene_qc_metrics("tiny.h5ad", top_n=2)
print(result["genes_detected"], result["median_cells_per_gene"])
for gene in result["top_genes"]:
    print(gene["gene"], gene["total_counts"], gene["n_cells_by_counts"])
```

- **示例输出：**

```text
4 2.0
GENE2 7.0 3
MT-A 4.0 2
```

## 3) `dnakit.transcriptomics.calculate_cell_qc_metrics()`

- **作用：** 汇总每个细胞的基因数、总计数和线粒体计数比例，并统计超过阈值的细胞数。
- **处理方法：** 名称以 `mito_prefix` 开头的基因视为线粒体基因（人类常用 `"MT-"`，小鼠常用 `"mt-"`）。调用 `scanpy.pp.calculate_qc_metrics(qc_vars=["mt"])` 得到每个细胞的 `pct_counts_mt`，即线粒体计数占该细胞总计数的百分比，再输出中位数、最大值和超过阈值的细胞数。结果只包含汇总值，不返回逐细胞表。
- **API：** `dnakit.transcriptomics.calculate_cell_qc_metrics(source, *, mito_prefix="MT-", mito_threshold_percent=20.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`mito_threshold_percent` 取 (0, 100] 中的数值，单位为百分比。
- **输出：** `n_cells`、`n_genes`、`mito_genes`、`median_genes_per_cell`、`median_total_counts`、`median_pct_counts_mt`、`max_pct_counts_mt`、`cells_above_mito_threshold` 和 `provenance`。
- **示例代码：**

```python
from dnakit.transcriptomics import calculate_cell_qc_metrics

result = calculate_cell_qc_metrics("tiny.h5ad", mito_threshold_percent=20)
print(result["mito_genes"], result["median_total_counts"], result["cells_above_mito_threshold"])
print(result["max_pct_counts_mt"])
```

- **示例输出：**

```text
1 3.5 2
75.0
```

## 4) `dnakit.transcriptomics.filter_genes_by_min_cells()`

- **作用：** 删除在少于 `min_cells` 个细胞中检测到的基因，并可写出新的 `.h5ad`。
- **处理方法：** 调用 `scanpy.pp.filter_genes(min_cells=...)`。基因在一个细胞中“被检测到”指计数大于 0。`output_path` 指向的文件如果已存在，会报错而不覆盖。
- **API：** `dnakit.transcriptomics.filter_genes_by_min_cells(source, *, min_cells=3, output_path=None, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`min_cells` 为正整数；`output_path` 可选，必须以 `.h5ad` 结尾。
- **输出：** `n_cells`、`n_genes_before`、`n_genes_after`、`n_genes_removed`、`output_path` 和 `provenance`。
- **示例代码：**

```python
from pathlib import Path
import anndata as ad
from dnakit.transcriptomics import filter_genes_by_min_cells

Path("filtered.h5ad").unlink(missing_ok=True)  # the function never overwrites an existing file
result = filter_genes_by_min_cells("tiny.h5ad", min_cells=2, output_path="filtered.h5ad")
print(result["n_genes_before"], result["n_genes_after"], result["n_genes_removed"])
print(ad.read_h5ad("filtered.h5ad").var_names.tolist())
```

- **示例输出：**

```text
4 3 1
['MT-A', 'GENE1', 'GENE2']
```

## 5) `dnakit.transcriptomics.select_highly_variable_genes()`

- **作用：** 选出高变基因。
- **处理方法：** `log_normalize=True` 时，先在副本上做 `scanpy.pp.normalize_total(target_sum)` 和 `scanpy.pp.log1p`，原文件不变。然后调用 `scanpy.pp.highly_variable_genes`，`flavor` 取 `"seurat"` 或 `"cell_ranger"`。结果按归一化离散度降序排列，离散度相同时按基因名排序。离散度相同的基因会一起被标记，因此实际数量可能多于 `n_top_genes`。`seurat_v3` 需要原始计数和 scikit-misc，这里没有包装。
- **API：** `dnakit.transcriptomics.select_highly_variable_genes(source, *, n_top_genes=2000, flavor="seurat", log_normalize=True, target_sum=10000.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`n_top_genes` 不能超过基因数；`flavor` 取 `"seurat"` 或 `"cell_ranger"`；`log_normalize=False` 时要求输入已经是对数归一化数据。
- **输出：** `flavor`、`n_top_genes`、`n_genes_input`、`n_highly_variable`、`highly_variable_genes`（按离散度排序的基因名）、`log_normalized`、`target_sum` 和 `provenance`。
- **示例代码：**

```python
import anndata as ad
import numpy as np
import pandas as pd
from dnakit.transcriptomics import select_highly_variable_genes

rng = np.random.default_rng(0)
counts = rng.poisson(3, size=(30, 20)).astype("float32")
counts[:15, :4] += 8
obs = pd.DataFrame({"group": ["A"] * 15 + ["B"] * 15}, index=[f"c{i}" for i in range(30)])
ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(20)])).write_h5ad("cells.h5ad")

result = select_highly_variable_genes("cells.h5ad", n_top_genes=4)
print(result["n_highly_variable"], result["highly_variable_genes"])
```

- **示例输出：**

```text
7 ['g5', 'g19', 'g10', 'g16', 'g18', 'g2', 'g3']
```

## 6) `dnakit.transcriptomics.rank_marker_genes()`

- **作用：** 对每个组与其余细胞比较，输出排名靠前的标记基因。
- **处理方法：** 调用 `scanpy.tl.rank_genes_groups`，`use_raw=False`，p 值用 Benjamini–Hochberg 校正。`method="wilcoxon"` 使用秩和检验，`method="t-test"` 使用 t 检验。`log_normalize=True` 时先在副本上归一化并取对数。`logfoldchange` 由 scanpy 计算，为 X 的 expm1 均值之比的 log2。
- **API：** `dnakit.transcriptomics.rank_marker_genes(source, *, groupby, method="wilcoxon", n_genes=10, log_normalize=True, target_sum=10000.0, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`groupby` 是 `obs` 中的列名，该列不能有缺失值且至少两组；`n_genes` 是每组输出的基因数。
- **输出：** `groupby`、`method`、`n_groups` 和 `groups`（组名到 `{gene, score, logfoldchange, pval_adj}` 列表的映射，组名排序）以及 `provenance`。
- **示例代码：**

```python
from dnakit.transcriptomics import rank_marker_genes

result = rank_marker_genes("cells.h5ad", groupby="group", n_genes=2)
for group, rows in result["groups"].items():
    print(group, [row["gene"] for row in rows])
```

- **示例输出：**

```text
A ['g0', 'g1']
B ['g13', 'g15']
```

## 7) `dnakit.transcriptomics.aggregate_pseudobulk()`

- **作用：** 按（分组, 样本）对原始计数求和，得到伪批量（pseudobulk）计数。
- **处理方法：** 求和由 `scanpy.get.aggregate(func="sum")` 完成，不做归一化。每个组合的细胞数由 pandas 对 `obs` 分组计数。细胞数少于 `min_cells_per_sample` 的组合被丢弃，并计入 `dropped_samples`。输出按（分组, 样本）排序。
- **API：** `dnakit.transcriptomics.aggregate_pseudobulk(source, *, group_key, sample_key, min_cells_per_sample=1, output_path=None, max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`group_key` 与 `sample_key` 是两个不同的 `obs` 列；`output_path` 可选，写出的 `.h5ad` 在 `layers["sum"]` 中保存求和结果，并在 `obs["n_cells"]` 中记录细胞数。
- **输出：** `n_pseudobulk_samples`、`dropped_samples`、`n_genes`、`samples`（每个保留组合的 `group`、`sample`、`n_cells` 和 `total_counts`）、`output_path` 和 `provenance`。
- **示例代码：**

```python
from dnakit.transcriptomics import aggregate_pseudobulk

result = aggregate_pseudobulk("tiny.h5ad", group_key="group", sample_key="sample")
print(result["n_pseudobulk_samples"], result["dropped_samples"])
for sample in result["samples"]:
    print(sample["group"], sample["sample"], sample["n_cells"], sample["total_counts"])
```

- **示例输出：**

```text
4 0
x s1 1 4.0
x s2 1 2.0
y s1 1 3.0
y s2 1 7.0
```

## 8) `dnakit.transcriptomics.summarize_mudata()`

- **作用：** 汇总 `.h5mu` 多模态文件中的模态、每个模态的维度，以及细胞数量。
- **处理方法：** 用 `muon.read_h5mu` 完整读取文件（受 `max_file_bytes` 限制）。`n_obs` 是各模态细胞的并集，`n_obs_shared` 是在所有模态中都有数据的细胞数。
- **API：** `dnakit.transcriptomics.summarize_mudata(source, *, max_file_bytes=1073741824)`。
- **输入：** `.h5mu` 文件路径。
- **输出：** `n_obs`、`n_obs_shared`、`modalities`（每个模态的 `n_obs` 和 `n_vars`）、`obs_columns` 和 `provenance`。
- **示例代码：**

```python
import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
from dnakit.transcriptomics import summarize_mudata

rng = np.random.default_rng(0)
rna = ad.AnnData(
    X=rng.poisson(1, size=(5, 3)).astype("float32"),
    obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]),
    var=pd.DataFrame(index=["G1", "G2", "G3"]),
)
atac = ad.AnnData(
    X=rng.poisson(1, size=(4, 6)).astype("float32"),
    obs=pd.DataFrame(index=["c0", "c1", "c2", "c9"]),
    var=pd.DataFrame(index=[f"p{i}" for i in range(6)]),
)
mu.write_h5mu("multi.h5mu", mu.MuData({"rna": rna, "atac": atac}))

result = summarize_mudata("multi.h5mu")
print(result["n_obs"], result["n_obs_shared"])
print(result["modalities"])
```

- **示例输出：**

```text
6 3
{'atac': {'n_obs': 4, 'n_vars': 6}, 'rna': {'n_obs': 5, 'n_vars': 3}}
```

## 9) `dnakit.transcriptomics.train_scvi_model()`

- **作用：** 在原始整数计数上设置并训练一个 scVI 模型，返回潜变量的形状和训练摘要。
- **处理方法：** 先检查矩阵是有限的非负整数，否则报错。然后依次调用 `SCVI.setup_anndata`、`SCVI(...)` 和 `train(...)`，最后用 `get_latent_representation()` 取后验均值。`seed` 会设置 `scvi.settings.seed`，这是进程级的全局设置。`accelerator` 默认 `"cpu"`，取 `"gpu"` 或 `"auto"` 时由 Lightning 决定设备，没有 GPU 时会报错。
- **API：** `dnakit.transcriptomics.train_scvi_model(source, *, layer=None, batch_key=None, n_latent=10, n_hidden=128, max_epochs=1, seed=0, accelerator="cpu", max_file_bytes=1073741824, max_matrix_elements=50000000)`。
- **输入：** `.h5ad` 文件路径；`layer` 指定计数所在的 layer；`batch_key` 是 `obs` 中的批次列；`n_latent` 取 1 到 256，`n_hidden` 取 1 到 4096，`max_epochs` 取 1 到 1000。
- **输出：** `n_cells`、`n_genes`、`latent_shape`、`latent_finite`、`latent_dimension_means`、`epochs_trained`、`final_train_elbo`、`seed`、`accelerator` 和 `provenance`。
- **示例代码：**

```python
import anndata as ad
import numpy as np
import pandas as pd
from dnakit.transcriptomics import train_scvi_model

rng = np.random.default_rng(0)
counts = rng.poisson(1.5, size=(60, 40)).astype("float32")
obs = pd.DataFrame({"batch": ["b1", "b2"] * 30}, index=[f"c{i}" for i in range(60)])
ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(40)])).write_h5ad(
    "scvi.h5ad"
)

result = train_scvi_model(
    "scvi.h5ad", batch_key="batch", n_latent=3, n_hidden=16, verbose=False
)
print(result["latent_shape"], result["epochs_trained"], result["latent_finite"])
```

- **示例输出：**

```text
[60, 3] 1 True
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-scanpy]"   # anndata 与 scanpy：摘要、QC、高变基因、标记基因、伪批量、filter
python -m pip install "dnakit[bio-muon]"     # muon：summarize_mudata
python -m pip install "dnakit[bio-scvi]"     # scvi-tools 与 anndata：train_scvi_model，需要 PyTorch
```

项目页面：[anndata](https://github.com/scverse/anndata)、[scanpy](https://github.com/scverse/scanpy)、[muon](https://github.com/scverse/muon)、[scvi-tools](https://github.com/scverse/scvi-tools)。

## 验证方式

测试文件是 `others/tests/unit/test_transcriptomics_backends.py`，共 18 个测试。在安装了 scanpy、muon 和 scvi-tools 的环境中全部运行；缺少某个包时，相关测试会跳过。测试检查以下内容：

- 基因 QC 和细胞 QC 的数值与手工计算一致，包括基因检测比例、中位数和线粒体比例。
- `filter_genes_by_min_cells` 的保留结果、写出文件的内容，以及拒绝覆盖已有文件。
- 伪批量求和与手工相加的结果一致，并按最小细胞数丢弃组合。
- 标记基因的 Wilcoxon z 值、Benjamini–Hochberg 校正后的 p 值和 log2 倍数变化与手工公式一致；归一化后的结果与直接调用 scanpy 一致。
- 高变基因（`seurat` 与 `cell_ranger`）的结果与直接调用 scanpy 一致。
- `summarize_mudata` 的并集细胞数与共享细胞数与手工计数一致。
- scVI 在 CPU 上训练一个 epoch，潜变量形状正确、结果有限，相同种子的两次运行结果一致；非整数计数和非法选项被拒绝。
- 缺少 `scanpy`、`anndata`、`scvi`、`muon` 时抛出 `BackendUnavailableError`，并给出对应的 extra 名称。

本批次不需要网络，因此没有 live 测试。

## 当前范围

本批次包装了 AnnData 摘要与 `.h5ad` 写出、scanpy 的基因和细胞 QC、高变基因（`seurat`、`cell_ranger`）、标记基因（`wilcoxon`、`t-test`）、伪批量求和、muon 的 `.h5mu` 摘要，以及 scVI 的设置与训练。

未包装的内容包括：scanpy 的聚类、降维和邻接图（如 `pca`、`neighbors`、`umap`、`leiden`）；回归校正和细胞周期评分；`seurat_v3` 高变基因；scVI 之外的模型（如 totalVI、scANVI）及其差异表达；muon 的多模态预处理；10x、loom 等格式的读取。这些内容需要单独的批次，或者由上游函数直接调用。
