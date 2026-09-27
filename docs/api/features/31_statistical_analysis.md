# 统计与差异分析

## 1) 用距离矩阵检验组间差异 {#bio-35}

- **作用：** 用 PERMANOVA、ANOSIM 或 PERMDISP 比较样本组。
- **计算方法：** 先验证距离矩阵，再置换分组标签计算 p 值。
- **API：** `dnakit.bioinformatics.statistics.distance_matrix(values)`; `dnakit.bioinformatics.statistics.distance_test(distances, groups, method='permanova', permutations=999, seed=0)`
- **输入：** 对称非负距离矩阵、对应组标签；PERMDISP 要求欧氏距离；需要 `bio-analysis`。
- **输出：** `distance_matrix` 返回距离矩阵；`distance_test` 返回统计量、p 值和检验方法。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import distance_matrix, distance_test

distances = distance_matrix([[0, 1, 4, 5], [1, 0, 3, 4], [4, 3, 0, 1], [5, 4, 1, 0]])
result = distance_test(distances, ["a", "a", "b", "b"], permutations=99, seed=0)
print(result.method, result.statistic, result.pvalue)
```

- **示例结果：**

```text
permanova 32.0 0.36
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 2) 比较距离矩阵和寻找关联环境变量 {#bio-36}

- **作用：** 用 Mantel 检验两个距离矩阵的关联，用 BIOENV 筛选环境变量子集。
- **计算方法：** Mantel 对矩阵上三角相关系数做置换；BIOENV 枚举小规模变量子集。
- **API：** `dnakit.bioinformatics.statistics.mantel(left, right, method='pearson', permutations=999, seed=0)`; `dnakit.bioinformatics.statistics.bioenv(distances, environmental, max_columns=12)`
- **输入：** 相同样本顺序的距离矩阵，或距离矩阵与样本×环境变量矩阵；需要 `bio-analysis`。
- **输出：** `mantel` 返回相关统计量和 p 值；变量选择返回关联变量及分数。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import mantel, bioenv

distances = [[0, 1, 3], [1, 0, 2], [3, 2, 0]]
test = mantel(distances, distances, permutations=99, seed=0)
selection = bioenv(distances, [[0, 0], [1, 2], [3, 1]])
print(round(test.statistic, 6))
print(selection[0])
```

- **示例结果：**

```text
1.0
{'columns': (0,), 'correlation': 1.0}
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 3) 运行 PCoA、CA、CCA 或 RDA 排序 {#bio-37}

- **作用：** 把样本间结构转换为低维坐标。
- **计算方法：** PCoA 对距离矩阵特征分解；CA/CCA/RDA 对数据矩阵进行无约束或约束排序。
- **API：** `dnakit.bioinformatics.statistics.pcoa(distances, dimensions=None)`; `dnakit.bioinformatics.statistics.ordination(values, method='ca', constraints=None)`
- **输入：** PCoA 输入距离矩阵；CA/CCA 输入非负矩阵；CCA/RDA 另需约束矩阵；需要 `bio-analysis`。
- **输出：** 返回 `OrdinationResult`，包含样本坐标、解释方差和轴信息。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import pcoa, ordination

coordinates = pcoa([[0, 1, 3], [1, 0, 2], [3, 2, 0]], dimensions=1)
ca = ordination([[1, 2, 3], [3, 2, 1], [2, 4, 1]], method="ca")
print(coordinates.samples.shape, ca.samples.shape)
```

- **示例结果：**

```text
(3, 1) (3, 2)
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 4) 进行 CLR、ILR 或 ALR 组成变换 {#bio-38}

- **作用：** 将正组成数据转换到对数比例空间，并可逆变换。
- **计算方法：** CLR 减去行内平均对数；ILR 使用 Helmert 正交基；ALR 相对指定列取比值。
- **API：** `dnakit.bioinformatics.statistics.composition_transform(values, method='clr', reference=0, inverse=False)`
- **输入：** 严格正值 samples×features 矩阵；零值需自行处理；逆变换返回和为 1 的组成；需要 `bio-analysis`。
- **输出：** 返回变换后的数值矩阵；逆变换返回恢复后的矩阵。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import composition_transform

encoded = composition_transform([[1, 2, 1]], method="clr")
recovered = composition_transform(encoded, method="clr", inverse=True)
print(encoded.round(6).tolist())
print(recovered.round(6).tolist())
```

- **示例结果：**

```text
[[-0.231049, 0.462098, -0.231049]]
[[0.25, 0.5, 0.25]]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 5) 使用 ANCOM 比较组成数据 {#bio-39}

- **作用：** 统计每个 taxon 与其他 taxon 的 log-ratio 差异数量。
- **计算方法：** 检验 pairwise log-ratio，按 taxon 做 Holm 校正并返回 W。
- **API：** `dnakit.bioinformatics.statistics.ancom(values, groups, alpha=0.05, pseudocount=0.0)`
- **输入：** 非负丰度矩阵和组标签；含零时显式设置 pseudocount；需要 `bio-analysis`。W 的判定阈值需自行确定。
- **输出：** 返回字典，包含 ANCOM 的 `W` 统计量、p 值和显著性结果。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import ancom

result = ancom([[1, 2, 3], [2, 3, 4], [8, 2, 3], [9, 3, 4]],
               ["a", "a", "b", "b"])
print(result["W"].tolist())
print(result["pvalues"].shape)
```

- **示例结果：**

```text
[0, 0, 0]
(3, 3)
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 6) 计算 alpha 和 beta 多样性 {#bio-40}

- **作用：** 计算样本内部多样性或样本之间的组成距离。
- **计算方法：** alpha 支持 observed/Shannon/Simpson/inverse Simpson；beta 支持 Bray–Curtis/Jaccard/Euclidean 等。
- **API：** `dnakit.bioinformatics.statistics.alpha_diversity(counts, metric='shannon')`; `dnakit.bioinformatics.statistics.beta_diversity(counts, metric='braycurtis')`
- **输入：** 非负 samples×taxa 计数矩阵；需要 `bio-analysis`。
- **输出：** `alpha_diversity` 返回每个样本的多样性向量；`beta_diversity` 返回样本间距离矩阵。
- **示例代码：**

```python
from dnakit.bioinformatics.statistics import alpha_diversity, beta_diversity

counts = [[1, 1], [2, 0]]
print(alpha_diversity(counts).round(6).tolist())
print(beta_diversity(counts).round(6).tolist())
```

- **示例结果：**

```text
[0.693147, -0.0]
[[0.0, 0.5], [0.5, 0.0]]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 7) 计算系统发育 UniFrac 距离 {#bio-41}

- **作用：** 利用系统发育树比较群落组成。
- **计算方法：** 按枝上出现情况计算 unweighted UniFrac，或按丰度计算 weighted UniFrac。
- **API：** `dnakit.bioinformatics.phylogeny.unifrac(counts, taxa, tree, weighted=False, normalized=True)`
- **输入：** 非负丰度矩阵、与树叶完全一致的 taxa 列、定根且枝长非负的树；需要 `bio-analysis`。
- **输出：** 返回样本间 UniFrac 距离矩阵。
- **示例代码：**

```python
from dnakit.bioinformatics.phylogeny import TreeNode, unifrac

tree = TreeNode(children=(TreeNode("a", 1), TreeNode("b", 1)))
print(unifrac([[1, 0], [0, 1]], ["a", "b"], tree).tolist())
```

- **示例结果：**

```text
[[0.0, 1.0], [1.0, 0.0]]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 8) 使用 ANCOM-BC 估计丰度差异 {#bio-42}

- **作用：** 估计设计变量对应的丰度效应、标准误和 p/q 值。
- **计算方法：** 使用 OLS、稳健 sandwich 方差和三高斯分量 EM 校正采样偏差。
- **API：** `dnakit.bioinformatics.counts.ancom_bc(counts, design, tolerance=1e-5, max_iterations=1000)`
- **输入：** 正值 samples×taxa 矩阵（至少 10 taxa）、含截距的满秩 design；需要 `bio-analysis`。示例为合成数据；[推断范围](../../faq/bioinformatics.md#statistics)。
- **输出：** 返回 `AbundanceResult`，包含系数、p 值和校正后的 q 值。
- **示例代码：**

```python
import numpy as np
from dnakit.bioinformatics.counts import ancom_bc

rng = np.random.default_rng(13)
design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
effects = np.r_[np.full(6, 1.2), np.full(6, -1.0), np.zeros(48)]
abundance = np.exp(rng.normal(4, 0.7, (30, 60)) + design[:, 1, None] * effects)
result = ancom_bc(abundance, design)
print(result.coefficients.shape, result.qvalues.shape)
```

- **示例结果：**

```text
(60, 2) (60, 2)
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 9) 使用 NB2 GLM 分析基因表达差异 {#bio-43}

- **作用：** 估计指定 contrast 的 log2 fold change、标准误和 p/q 值。
- **计算方法：** 拟合未收缩 NB2 GLM，执行 Wald 检验并进行 BH 校正。
- **API：** `dnakit.bioinformatics.counts.differential_expression(counts, design, contrast, size_factors=None, dispersions=None)`
- **输入：** 整数 samples×genes counts、数值 design、contrast，至少两个残差自由度；真实分析需生物学重复。需要 `bio-analysis`；[与 DESeq2/edgeR 的区别](../../faq/bioinformatics.md#statistics)。
- **输出：** 返回 `DifferentialExpressionResult`，包含 log2 倍数变化、p 值和收敛状态。
- **示例代码：**

```python
import numpy as np
from dnakit.bioinformatics.counts import differential_expression

rng = np.random.default_rng(18)
design = np.column_stack([np.ones(30), np.repeat([0, 1], 15)])
mean = np.exp(design @ np.array([[3, 4, 5], [0.8, -1.0, 0.0]]))
counts = rng.negative_binomial(10, 10 / (mean + 10))
result = differential_expression(counts, design, [0, 1],
                                 size_factors=np.ones(30), dispersions=0.1)
print(result.log2_fold_change.round(3).tolist())
print(result.converged.tolist())
```

- **示例结果：**

```text
[1.319, -1.501, -0.011]
[True, True, True]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。

## 10) 校正多重检验 p 值 {#bio-44}

- **作用：** 将一组 p 值转换为 Benjamini–Hochberg q 值。
- **计算方法：** 按 p 值排序后调整，并恢复原有顺序。
- **API：** `dnakit.bioinformatics.counts.adjust_pvalues(values)`
- **输入：** 一维 p 值数组；需要 `bio-analysis`。
- **输出：** 返回校正后的 p 值数组。
- **示例代码：**

```python
from dnakit.bioinformatics.counts import adjust_pvalues

print(adjust_pvalues([0.01, 0.04, 0.03]).tolist())
```

- **示例结果：**

```text
[0.03, 0.04, 0.04]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#statistics)。
