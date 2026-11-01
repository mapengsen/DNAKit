# scikit-bio 表增强与嵌入排序

这一组函数位于 `dnakit.microbiome.skbio_p21_functions`，用 [scikit-bio](https://scikit.bio/) 对样本×特征表（如微生物组或表达计数）做数据增强（mixup、Aitchison mixup、组成 cutmix、PhyloMix），并对嵌入向量做 SVD 排序。依赖为可选的 `dnakit[bio-skibio]`；所有随机操作都接受整数 `seed`。

## 1) `dnakit.microbiome.skbio_p21_functions.mixup_table()`

- **作用：** 由随机样本对按 Beta 系数线性组合生成合成样本。
- **处理方法：** Each synthetic sample is a linear combination of a random sample pair with a Beta(alpha, alpha) coefficient (`skbio.table.mixup`). `intra_class` restricts mixing to one class and `append` also returns the original samples.
- **API：** `labels=[0, 1, 0, 1], `。
- **输入：** 样本×特征的表（至少 2 个样本）、合成样本数 `n`，可选类别索引 `labels`。
- **输出：** `matrix`（合成样本）、`labels`（独热编码，提供标签时）、`n_synthetic_rows`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import mixup_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = mixup_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **示例输出：**

```text
2 [[2.127, 1.859, 2.718, 3.577, 4.436], [3.431, 2.954, 2.477, 2.0, 1.523]]
```

## 2) `dnakit.microbiome.skbio_p21_functions.aitchison_mixup_table()`

- **作用：** 在 Aitchison 几何中做 mixup 生成合成组成。
- **处理方法：** Mixup in Aitchison geometry (log-ratio space) (`skbio.table.aitchison_mixup`); `normalize` closes the rows first when needed.
- **API：** ``。
- **输入：** 样本×特征的表（至少 2 个样本）、合成样本数 `n`，可选类别索引 `labels`。
- **输出：** `matrix`（合成样本）、`labels`（独热编码，提供标签时）、`n_synthetic_rows`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import aitchison_mixup_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = aitchison_mixup_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **示例输出：**

```text
2 [[0.105, 0.139, 0.197, 0.253, 0.306], [0.264, 0.237, 0.207, 0.17, 0.122]]
```

## 3) `dnakit.microbiome.skbio_p21_functions.compositional_cutmix_table()`

- **作用：** 在同类样本之间做组成 cutmix 生成合成样本。
- **处理方法：** For a pair of samples of one class each feature is taken from one of them at random and the result is closed to sum 1 (`skbio.table.compos_cutmix`).
- **API：** `labels=[0, 1, 0, 1], `。
- **输入：** 样本×特征的表（至少 2 个样本）、合成样本数 `n`，可选类别索引 `labels`。
- **输出：** `matrix`（合成样本）、`labels`（独热编码，提供标签时）、`n_synthetic_rows`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import compositional_cutmix_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = compositional_cutmix_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **示例输出：**

```text
2 [[0.067, 0.133, 0.2, 0.267, 0.333], [0.2, 0.2, 0.2, 0.2, 0.2]]
```

## 4) `dnakit.microbiome.skbio_p21_functions.phylomix_table()`

- **作用：** 利用系统发育树混合相关类群来生成合成样本（PhyloMix）。
- **处理方法：** 随机选取树上的节点，在两个样本间混合该节点下所有末端类群的丰度（`skbio.table.phylomix`）。
- **API：** `phylomix_table(table, n, tree_newick, taxa, *, labels=None, intra_class=False, alpha=2.0, append=False, seed=None)`。
- **输入：** 样本×特征的表、合成样本数 `n`、Newick 树和与列对应的末端名称 `taxa`。
- **输出：** `matrix`、`labels`、`n_synthetic_rows`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import phylomix_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = phylomix_table(table, 2, "(((a:1,b:1):1,c:2):1,(d:1,e:1):2);", list("abcde"), seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **示例输出：**

```text
2 [[9.0, 1.0, 1.0, 1.0, 1.0], [2.667, 4.0, 3.0, 2.667, 2.667]]
```

## 5) `dnakit.microbiome.skbio_p21_functions.ordinate_embedding_vectors()`

- **作用：** 对嵌入向量（例如序列模型的嵌入）做奇异值分解排序。
- **处理方法：** 调用 `skbio.embedding.embed_vec_to_ordination`：对向量矩阵做 SVD，样本坐标为 U·S，解释比例为奇异值平方占比。
- **API：** `ordinate_embedding_vectors(vectors)`。
- **输入：** 对象 ID 到等长嵌入列表的映射（至少 3 个）。
- **输出：** `ids`、`coordinates`、`eigvals`、`proportion_explained`。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p21_functions import ordinate_embedding_vectors
vectors = {"a": [1, 2, 3], "b": [2, 1, 0], "c": [0, 5, 1], "d": [3, 3, 3]}
result = ordinate_embedding_vectors(vectors, verbose=False)
print(result["ids"], [round(v, 4) for v in result["proportion_explained"]])
```

- **示例输出：**

```text
['a', 'b', 'c', 'd'] [0.8086, 0.1395, 0.052]
```
