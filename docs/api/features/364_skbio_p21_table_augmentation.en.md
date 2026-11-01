# scikit-bio table augmentation and embedding ordination

These functions in `dnakit.microbiome.skbio_p21_functions` use [scikit-bio](https://scikit.bio/) to augment sample x feature tables (for example microbiome or expression counts) with mixup, Aitchison mixup, compositional cutmix and PhyloMix, and to ordinate embedding vectors by SVD. The dependency is the optional `dnakit[bio-skibio]`; every random operation takes an integer `seed`.

## 1) `dnakit.microbiome.skbio_p21_functions.mixup_table()`

- **Purpose:** 由随机样本对按 Beta(alpha, alpha) 系数线性组合生成合成样本（`skbio.table.mixup`）。`intra_class` 限制在同类内混合，`append` 把原样本一并返回。
- **Method:** mixup_table(table, n, *, labels=None, intra_class=False, alpha=2.0, append=False, seed=None)
- **API:** `labels=[0, 1, 0, 1], `.
- **Input:** A samples x features table (at least 2 samples), the number `n` of synthetic samples, and optional class indices `labels`.
- **Output:** `matrix` (synthetic samples), `labels` (one-hot, when labels are given) and `n_synthetic_rows`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import mixup_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = mixup_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **Example output:**

```text
2 [[2.127, 1.859, 2.718, 3.577, 4.436], [3.431, 2.954, 2.477, 2.0, 1.523]]
```

## 2) `dnakit.microbiome.skbio_p21_functions.aitchison_mixup_table()`

- **Purpose:** 在 Aitchison 几何（对数比空间）中做 mixup（`skbio.table.aitchison_mixup`）；`normalize` 在需要时先做闭合。
- **Method:** aitchison_mixup_table(table, n, *, labels=None, intra_class=False, alpha=2.0, normalize=True, append=False, seed=None)
- **API:** ``.
- **Input:** A samples x features table (at least 2 samples), the number `n` of synthetic samples, and optional class indices `labels`.
- **Output:** `matrix` (synthetic samples), `labels` (one-hot, when labels are given) and `n_synthetic_rows`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import aitchison_mixup_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = aitchison_mixup_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **Example output:**

```text
2 [[0.105, 0.139, 0.197, 0.253, 0.306], [0.264, 0.237, 0.207, 0.17, 0.122]]
```

## 3) `dnakit.microbiome.skbio_p21_functions.compositional_cutmix_table()`

- **Purpose:** 在同类样本对之间逐特征随机选取来源，再闭合为和为 1 的组成（`skbio.table.compos_cutmix`）。
- **Method:** compositional_cutmix_table(table, n, *, labels=None, normalize=True, append=False, seed=None)
- **API:** `labels=[0, 1, 0, 1], `.
- **Input:** A samples x features table (at least 2 samples), the number `n` of synthetic samples, and optional class indices `labels`.
- **Output:** `matrix` (synthetic samples), `labels` (one-hot, when labels are given) and `n_synthetic_rows`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import compositional_cutmix_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = compositional_cutmix_table(table, 2, seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **Example output:**

```text
2 [[0.067, 0.133, 0.2, 0.267, 0.333], [0.2, 0.2, 0.2, 0.2, 0.2]]
```

## 4) `dnakit.microbiome.skbio_p21_functions.phylomix_table()`

- **Purpose:** Generate synthetic samples by PhyloMix, mixing the taxa below a random tree node.
- **Method:** A random tree node is chosen and the abundances of all tips below it are mixed between two samples (`skbio.table.phylomix`).
- **API:** `phylomix_table(table, n, tree_newick, taxa, *, labels=None, intra_class=False, alpha=2.0, append=False, seed=None)`.
- **Input:** A samples x features table, the number `n`, a Newick tree and the tip names `taxa` for the columns.
- **Output:** `matrix`, `labels` and `n_synthetic_rows`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import phylomix_table
table = [[1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [2, 2, 2, 2, 2], [9, 1, 1, 1, 1]]
result = phylomix_table(table, 2, "(((a:1,b:1):1,c:2):1,(d:1,e:1):2);", list("abcde"), seed=1, verbose=False)
print(result["n_synthetic_rows"], [[round(v, 3) for v in row] for row in result["matrix"]])
```

- **Example output:**

```text
2 [[9.0, 1.0, 1.0, 1.0, 1.0], [2.667, 4.0, 3.0, 2.667, 2.667]]
```

## 5) `dnakit.microbiome.skbio_p21_functions.ordinate_embedding_vectors()`

- **Purpose:** Ordinate embedding vectors (for example sequence-model embeddings) by singular value decomposition.
- **Method:** Calls `skbio.embedding.embed_vec_to_ordination`: an SVD of the vector matrix, with sample coordinates U·S and proportion explained from the squared singular values.
- **API:** `ordinate_embedding_vectors(vectors)`.
- **Input:** A mapping of object id to equal-length embedding lists (at least 3).
- **Output:** `ids`, `coordinates`, `eigvals` and `proportion_explained`.
- **Example code:**

```python
from dnakit.microbiome.skbio_p21_functions import ordinate_embedding_vectors
vectors = {"a": [1, 2, 3], "b": [2, 1, 0], "c": [0, 5, 1], "d": [3, 3, 3]}
result = ordinate_embedding_vectors(vectors, verbose=False)
print(result["ids"], [round(v, 4) for v in result["proportion_explained"]])
```

- **Example output:**

```text
['a', 'b', 'c', 'd'] [0.8086, 0.1395, 0.052]
```
