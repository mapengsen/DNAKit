# 二级结构性质

解析 DNA 二级结构并计算碱基配对、茎环、配对概率、可接近性、集合缺陷和目标结构概率。

原生功能支持 dot-bracket 解析和配对概率派生指标；NUPACK 分析只在用户已获许可并独立安装 NUPACK 4 后执行，项目不会自动下载、安装或静默调用 NUPACK。

## 1) Dot-bracket 结构解析

- **作用：** 解析 dot-bracket 字符串，恢复碱基配对关系并统计 stem、hairpin、未配对碱基等结构元素，用于检查和汇总已有二级结构预测。
- **API：** `dnakit.secondary_structure.analyze_dot_bracket(strands[必须], dot_bracket[必须], three_prime_window[可选])`
- **输入：** 必填 DNA 链集合和长度一致的 dot-bracket；多链用 `+` 分隔，支持 `() [] {} <>`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import analyze_dot_bracket

result = analyze_dot_bracket(
    (DNASequence("ATCCTAGTTATAGGAT"),),
    "((((((....))))))",
)
print(result.structure_type, result.base_pair_count)
print(result.stem_lengths, result.hairpin_loop_lengths)
```

- **示例结果：**

```text
hairpin 6
(6,) (4,)
```

- **限制：** 该 API 只解析已提供的结构注释，不预测 MFE 或自由能。扩展括号只保留配对关系，不等同于完整 pseudoknot loop decomposition。

## 2) 配对概率与窗口可接近性

- **作用：** 校验配对概率矩阵的维度、对称性和概率范围，计算每个位置的配对与未配对概率，为 ensemble 指标提供可靠输入。
- **API：** `dnakit.secondary_structure.pair_probability_metrics(strands[必须], probability_matrix[必须], accessibility_window_size[可选])`
- **输入：** 必填 DNA 链集合和对称方阵；对角线是未配对概率，每行含对角线后的和必须为 1；可选可接近性窗口长度。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import pair_probability_metrics

matrix = (
    (0.7, 0.0, 0.0, 0.3),
    (0.0, 0.8, 0.2, 0.0),
    (0.0, 0.2, 0.8, 0.0),
    (0.3, 0.0, 0.0, 0.7),
)
result = pair_probability_metrics(
    (DNASequence("ACGT"),),
    matrix,
    accessibility_window_size=2,
)
print(tuple(round(value, 1) for value in result.pairing_probabilities_by_base))
print(tuple(round(value, 1) for value in result.unpaired_probabilities_by_base))
print(result.most_accessible_window_start)
```

- **示例结果：**

```text
(0.3, 0.2, 0.2, 0.3)
(0.7, 0.8, 0.8, 0.7)
1
```

- **限制：** 窗口可接近性是边缘未配对概率的算术平均，不是窗口所有碱基同时未配对的联合概率。

## 3) 目标结构的集合缺陷

- **作用：** 根据目标配对关系和 ensemble 配对概率计算 ensemble defect，估算未按目标结构正确配对的碱基期望数量及比例。
- **API：** `dnakit.secondary_structure.ensemble_defect_from_probabilities(target[必须], probabilities[必须])`
- **输入：** 必填 `SecondaryStructureSummary` 目标和描述相同 DNA 链的 `PairProbabilityResult`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import (
    analyze_dot_bracket,
    ensemble_defect_from_probabilities,
    pair_probability_metrics,
)

matrix = (
    (0.7, 0.0, 0.0, 0.3),
    (0.0, 0.8, 0.2, 0.0),
    (0.0, 0.2, 0.8, 0.0),
    (0.3, 0.0, 0.0, 0.7),
)
strands = (DNASequence("ACGT"),)
target = analyze_dot_bracket(strands, "(())")
probabilities = pair_probability_metrics(strands, matrix, accessibility_window_size=2)
print(ensemble_defect_from_probabilities(target, probabilities))
```

- **示例结果：**

```text
0.75
```

- **限制：** 目标和概率结果必须对应完全相同的链与顺序；这是由输入概率派生的指标，不会自行计算概率矩阵。

## 4) 目标结构热力学概率

- **作用：** 根据目标结构自由能与配分函数计算该结构在热力学 ensemble 中的理论概率，用于比较候选结构的相对占比。
- **API：** `dnakit.secondary_structure.target_structure_probability(target_free_energy_kcal_per_mol[必须], ensemble_free_energy_kcal_per_mol[必须], temperature_celsius[可选])`
- **输入：** 必填目标结构和集合的 kcal/mol 自由能；可选 0–100 °C 温度。
- **示例代码：**

```python
from dnakit.secondary_structure import target_structure_probability

same_energy = target_structure_probability(-2.0, -2.0)
higher_target_energy = target_structure_probability(-1.0, -2.0)
print(same_energy, round(higher_target_energy, 6))
```

- **示例结果：**

```text
1.0 0.197404
```

- **限制：** 两个自由能必须使用相同温度、标准态和模型得到；该 API 不会评估自由能来源是否一致。

## 5) NUPACK 被动可用性检查

- **作用：** 检查当前环境能否导入并调用兼容的 NUPACK 后端，返回版本和能力信息，使结构计算在运行前明确报告可用性。
- **API：** `dnakit.secondary_structure.probe_nupack()`
- **输入：** 无。
- **示例代码：**

```python
from dnakit.secondary_structure import probe_nupack

status = probe_nupack()
print(status.available, status.version, status.metadata["import_executed"])
```

- **当前本地示例结果：**

```text
False None False
```

- **限制：** `available=False` 只说明当前环境没有可执行 NUPACK，不是二级结构计算结果。探测过程不验证数值准确性。

## 6) NUPACK 单复合物集合分析

- **作用：** 调用 NUPACK 对一个或多个 DNA 链计算 MFE 结构、自由能、配对概率或结构抽样，并转换为 DNAKit 的统一结果对象。
- **API：** `dnakit.secondary_structure.NupackAdapter.analyze_complex(strands[必须], conditions[可选], ensemble[可选], target_structure[可选], suboptimal_energy_gap_kcal_per_mol[可选], num_samples[可选], accessibility_window_size[可选])`
- **输入：** 必填 DNA 链集合；可选 NUPACK 条件、ensemble、目标 dot-bracket、0–20 kcal/mol 次优差值和 0–100000 个抽样。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import NupackAdapter, probe_nupack

status = probe_nupack()
if status.available:
    result = NupackAdapter().analyze_complex(
        (DNASequence("ATCCTAGTTATAGGAT"),),
        target_structure="((((((....))))))",
        num_samples=100,
    )
    print(result.mfe_structures[0].summary.dot_bracket)
    print(result.ensemble_free_energy_kcal_per_mol)
else:
    print("NUPACK unavailable")
```

- **当前本地示例结果：**

```text
NUPACK unavailable
```

- **限制：** 必须由用户另行获得许可并安装 NUPACK 4。当前环境没有真实 NUPACK 数值结果；受控替身测试只验证 adapter 字段、边界和错误处理。

## 7) NUPACK tube 多复合物平衡

- **作用：** 调用 NUPACK 在给定链浓度和温度下分析 tube 内多个复合物的平衡组成，返回各复合物浓度及质量守恒信息。
- **API：** `dnakit.secondary_structure.NupackAdapter.analyze_tube(strands[必须], concentrations_molar[必须], target_strand_names[必须], conditions[可选], max_complex_size[可选])`
- **输入：** 必填命名 DNA 链映射、键集合相同的摩尔浓度映射和目标复合物链名；最多 20 种链，`max_complex_size` 为 1–4。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import NupackAdapter, probe_nupack

status = probe_nupack()
if status.available:
    result = NupackAdapter().analyze_tube(
        {"a": DNASequence("CCC"), "b": DNASequence("GGG")},
        {"a": 1e-6, "b": 1e-6},
        target_strand_names=("a", "b"),
        max_complex_size=2,
    )
    print(result.target_complex_concentration_molar)
    print(result.target_complex_fraction, result.non_target_complex_fraction)
else:
    print("NUPACK unavailable")
```

- **当前本地示例结果：**

```text
NUPACK unavailable
```

- **限制：** 目标比例的分母是本次枚举结果中所有复合物平衡浓度之和，并保存在 `complex_fraction_denominator_molar`；它不是按碱基数或投料链总浓度归一化的产率。

NUPACK 的外部功能边界可查阅其[官方分析文档](https://docs.nupack.org/analysis/)、[实用函数文档](https://docs.nupack.org/utilities/)和[模型文档](https://docs.nupack.org/model/)；安装前须遵守其单独的[下载与许可要求](https://www.nupack.org/download/overview)。

!!! warning
    不得把 `probe_nupack().available=False`、dot-bracket 注释或 Primer3 结果写成 NUPACK 计算结果。
