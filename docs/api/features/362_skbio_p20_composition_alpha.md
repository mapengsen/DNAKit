# scikit-bio 组成数据与 Alpha 多样性

本页介绍 P20 组成数据分析接口，位于 `dnakit.microbiome.skbio_p20_composition_functions` 与 `dnakit.microbiome.skbio_p20_stats_functions`。依赖可选项 `bio-skibio`（scikit-bio 0.7.4）。组成矩阵按行表示样本、按列表示特征；输入保持原始数值，不会自动将序列字母大小写转换等同处理。

## 1) 组成数据变换

- **作用：** 对组成矩阵执行闭合、零值替代、扰动、幂变换、对数比变换及结构零/差异丰度分析。
- **处理方法：** 包装 scikit-bio composition 算法，并对矩阵维度、有限值与参数范围进行校验。
- **API：** `dnakit.microbiome.skbio_p20_composition_functions.closure_compositions(matrix, *, verbose=True)`；零值替代接口为 `replace_zeros_multiplicative(matrix, *, delta=None, verbose=True)`。
- **输入：** 非负有限值组成矩阵；部分模型需要样本元数据或分组标签。
- **输出：** JSON 友好的变换矩阵、统计量和来源信息。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p20_composition_functions import (
    closure_compositions,
    replace_zeros_multiplicative,
)

print(closure_compositions([[2, 3, 5], [1, 1, 2]], verbose=False)["compositions"])
print(replace_zeros_multiplicative([[0, 2, 3]], verbose=False)["compositions"])
```

- **示例输出：**

```text
[[0.2, 0.3, 0.5], [0.25, 0.25, 0.5]]
[0.1111111111111111, 0.35555555555555557, 0.5333333333333333]
```

## 2) Alpha 多样性

- **作用：** 计算丰富度、均匀度和系统发育多样性等 alpha 指标。
- **处理方法：** `alpha_diversity_metric` 根据 `metric` 选择 scikit-bio 指标，并可传入特征树计算 Faith's PD。
- **API：** `dnakit.microbiome.skbio_p20_stats_functions.alpha_diversity_metric(counts, metric, *, sample_ids=None, feature_ids=None, tree_newick=None, verbose=True)`。
- **输入：** 每行一个样本的非负整数计数，及指标名；Faith's PD 还需 Newick 树和特征 ID。
- **输出：** 每个样本对应的多样性值与来源信息。
- **示例代码：**

```python
from dnakit.microbiome.skbio_p20_stats_functions import alpha_diversity_metric

result = alpha_diversity_metric(
    [[1, 2, 0], [0, 1, 3]], "dominance",
    sample_ids=["s1", "s2"], verbose=False,
)
print(result["values"])
```

- **示例输出：**

```text
[0.5555555555555556, 0.625]
```
