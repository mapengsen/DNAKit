# DNA 综合评价体系

从合法性、唯一性、多样性、新颖性、Fréchet 表征分布距离、片段分布、最近邻相似度、模糊度和冗余度等方面综合评价 DNA 序列或数据集。

`EVAL-011` 和 `EVAL-012` 在数据评价下单列于[相似度计算](09_similarity_alignment.md)。所有“参考相关”结果都必须绑定调用方提供的本地 `ReferenceLibrary`；novelty 不是实验结论或任务模型预测。

以下 API 括号内列出全部直接调用参数；`[必须]` 和 `[可选]` 是参数说明，不是 Python 语法。

## 1) `EVAL-001` Validity

- **作用：** 按字母表、长度、Gap、metadata 等规则检查 DNA 序列和记录，返回有效比例及逐条问题，用于在评价或建模前确认输入是否合法。
- **API**：`dnakit.evaluation.evaluate_validity(value[必须], config[可选], limits[可选])`；`config` 使用 `dnakit.ValidationConfig`，`limits` 使用 `dnakit.evaluation.EvaluationLimits`。
- **输入**：`DNASequence`、`DNARecord`、`DNASet` 或相应迭代器；可选 `ValidationConfig` 和资源上限。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_validity

report = evaluate_validity(DNASequence("ACGT"))
print(report.metrics["valid_fraction"])  # 1.0
```

- **示例结果：**

```text
1.0
```

- **限制**：这是对象与规则合法性检查，不证明生物学功能正确。

## 2) `EVAL-005` Uniqueness

- **作用：** 按精确、反向互补、IUPAC 或近似等价规则对序列分组，返回唯一序列数量、比例及重复组，用于衡量集合中独立记录的占比。
- **API**：`dnakit.evaluation.evaluate_uniqueness(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.UniquenessEvaluationConfig`。
- **输入**：`DNASet`；可选等价规则、方法和阈值。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_uniqueness

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAA"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
])
report = evaluate_uniqueness(records)
print(report.metrics["uniqueness_score"])  # 0.666666...
print(report.metrics["duplicate_groups"])  # (("a", "b"),)
```

- **示例结果：**

```text
0.6666666666666666
(('a', 'b'),)
```

- **限制**：近似模式是有界 pairwise 计算，大型数据库应使用已验证的外部索引。

## 3) `EVAL-006` Diversity

- **作用：** 计算序列集合的两两距离、最近邻距离和阈值 cluster 数量，返回整体及局部差异指标，用于判断样本覆盖是否广泛。
- **API**：`dnakit.evaluation.evaluate_diversity(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.DiversityEvaluationConfig`。
- **输入**：`DNASet`；可选相似度方法、k 和聚类阈值。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_diversity

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
])
report = evaluate_diversity(records)
print(report.metrics["mean_pair_distance"])
print(report.metrics["cluster_count"])
```

- **示例结果：**

```text
1.0
3
```

- **限制**：结果完全依赖所选相似度和阈值，不是绝对多样性。

## 4) `EVAL-008` Novelty

**参考库示例准备：**

本项使用版本化本地参考库。最小构造方式如下：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import create_reference_library

reference_records = DNASet([
    DNARecord(DNASequence("AAAA"), "ref-a"),
    DNARecord(DNASequence("CCCC"), "ref-c"),
])
reference = create_reference_library(
    reference_records,
    name="training",
    version="1",
    source="local:example",
)
queries = DNASet([
    DNARecord(DNASequence("AAAA"), "copy"),
    DNARecord(DNASequence("GGGG"), "query-new"),
])
```

- **作用：** 把每条查询序列与版本化参考库比较，按相似度和覆盖度阈值判断是否新颖，返回逐条最近命中及集合的新颖序列比例。
- **API**：`dnakit.evaluation.evaluate_novelty(queries[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.ReferenceSearchConfig`。
- **输入**：query 集合和 `ReferenceLibrary`；可选方法、阈值、k 和覆盖率。
- **示例代码**：接在“参考库示例准备”之后运行。

```python
from dnakit.evaluation import ReferenceSearchConfig, evaluate_novelty

report = evaluate_novelty(
    queries,
    reference,
    config=ReferenceSearchConfig(method="identity", copy_threshold=0.9),
)
print(report.metrics["novel_fraction"])  # 0.5
```

- **示例结果：**

```text
0.5
```

- **限制**：novelty 永远相对于给定参考库，不能脱离库版本声称“全局新颖”。

## 5) `EVAL-016` Fréchet DNA distance

- **作用：** 使用同一个 DNA 基础模型分别表示两个序列集合，将两组向量近似为多元高斯分布，再计算均值和协方差的 Fréchet 距离。数值越小表示两个集合在该表征空间中的分布越接近，完全相同的表征分布趋近于 `0`。
- **API**：`dnakit.evaluation.evaluate_frechet_distance(left[必须], right[必须], config[可选], backend[可选])`；`config` 使用 `dnakit.evaluation.FrechetDistanceConfig`，其中的表征配置使用 `dnakit.representations.RepresentationConfig`。
- **输入**：两个各含至少 2 条 `DNARecord` 的集合；默认复用 `DATA-027` 的 `lucaone`、mean pooling 和 L2 归一化。可选择其他已注册模型、checkpoint、pooling、设备、dtype 和批大小。
- **公式**：`||μ_left - μ_right||² + Tr(Σ_left + Σ_right - 2(Σ_left^(1/2) Σ_right Σ_left^(1/2))^(1/2))`。实现使用数学等价的样本空间 cross-Gram 核范数，不显式创建高维稠密协方差矩阵。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import FrechetDistanceConfig, evaluate_frechet_distance
from dnakit.representations import RepresentationConfig

reference = DNASet([
    DNARecord(DNASequence("ACGTACGT"), "ref-1"),
    DNARecord(DNASequence("AACCGGTT"), "ref-2"),
])
generated = DNASet([
    DNARecord(DNASequence("ACGTTCGT"), "gen-1"),
    DNARecord(DNASequence("AACCAGTT"), "gen-2"),
])
report = evaluate_frechet_distance(
    generated,
    reference,
    config=FrechetDistanceConfig(
        representation=RepresentationConfig(allow_remote_code=True),
    ),
)
print(report.metrics["frechet_distance"])
```

- **示例结果：**

```text
非负浮点数；具体值取决于模型、checkpoint、pooling、归一化和输入集合
```

- **进度：** checkpoint 下载和逐序列表征提取默认显示进度条；可用 `RepresentationConfig(show_progress=False)` 关闭。
- **限制**：这是借鉴 FCD/FID 数学形式的 **Fréchet DNA 表征距离**，不是使用 ChemNet 的分子 FCD，二者数值不可比较。有限样本的均值和协方差估计存在偏差；不同模型、checkpoint、pooling、归一化或精度生成的值也不可横向比较。LucaOne 标准后端包含 checkpoint 自带代码，仍须审查后显式设置 `allow_remote_code=True`。该距离不是实验功能、生成质量或生物安全结论。
- **依据**：[Preuer 等人的 FCD 原始论文](https://doi.org/10.1021/acs.jcim.8b00234)；[LucaOne 原始论文](https://doi.org/10.1038/s42256-025-01044-4)。

## 6) `EVAL-017` Frag

- **作用：** 参考分子生成评价中的 Frag，把生成集合和参考集合分别转换为片段出现次数向量，再计算余弦相似度。数值范围为 `[0,1]`，越高表示两个集合的局部片段分布越接近。
- **API**：`dnakit.evaluation.evaluate_fragment_similarity(generated[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.FragmentSimilarityConfig`。
- **DNA 适配：** 分子 Frag 使用 BRICS 片段；DNA 没有对应的化学断键规则，因此本实现明确使用 overlapping fixed-length k-mer。默认 `k=3`、canonical k-mer、忽略含 IUPAC 模糊碱基的窗口。
- **公式**：`Σ_f c_generated(f)c_reference(f) / sqrt(Σ_f c_generated(f)² × Σ_f c_reference(f)²)`。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import FragmentSimilarityConfig, evaluate_fragment_similarity

generated = DNASet([
    DNARecord(DNASequence("AAAA"), "gen-a"),
    DNARecord(DNASequence("CCCC"), "gen-c"),
])
reference = DNASet([
    DNARecord(DNASequence("AAAA"), "ref-a"),
    DNARecord(DNASequence("GGGG"), "ref-g"),
])
report = evaluate_fragment_similarity(
    generated,
    reference,
    config=FragmentSimilarityConfig(k=2, canonical=False, show_progress=False),
)
print(report.metrics["frag"])  # 0.5
```

- **示例结果：**

```text
0.5
```

- **进度：** 默认显示生成集合和参考集合的 k-mer 统计进度；可用 `FragmentSimilarityConfig(show_progress=False)` 关闭。
- **限制**：这是 MOSES Frag 的 DNA k-mer 适配，不是 BRICS Frag，数值不能与分子指标直接比较。结果依赖 `k`、canonical、模糊碱基和 Gap 策略；两个集合即使没有相同完整序列，也可能得到 `1`。
- **依据**：[MOSES 对 Frag 的定义](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full)。

## 7) `EVAL-018` SNN

- **作用：** 对每条生成 DNA，在参考集合中查找二进制指纹 Tanimoto 相似度最高的序列，再对这些最近邻相似度取算术平均。数值范围为 `[0,1]`，越高表示生成样本越接近参考集合覆盖的表征空间。
- **API**：`dnakit.evaluation.evaluate_snn(generated[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.SNNConfig`。
- **DNA 适配：** 分子 SNN 使用 1024 位、radius 2 Morgan 指纹；本实现默认使用 canonical 7-mer 经 SHA-256 映射得到的 1024 位 DNA 二进制指纹，并计算等价于二进制 Jaccard 的 Tanimoto 相似度。
- **公式**：`mean_g max_r Tanimoto(fp(g), fp(r))`；该指标以生成集合为查询，因此通常不对称。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import SNNConfig, evaluate_snn

generated = DNASet([
    DNARecord(DNASequence("AAAAAAA"), "copy"),
    DNARecord(DNASequence("ATATATA"), "far"),
])
reference = DNASet([
    DNARecord(DNASequence("AAAAAAA"), "ref-a"),
    DNARecord(DNASequence("CCCCCCC"), "ref-c"),
])
report = evaluate_snn(
    generated,
    reference,
    config=SNNConfig(show_progress=False),
)
print(report.metrics["snn"])  # 0.5
print(report.entries[0].metrics["nearest_reference_id"])  # ref-a
```

- **示例结果：**

```text
0.5
ref-a
```

- **进度：** 默认显示两组指纹构建和最近邻扫描进度；可用 `SNNConfig(show_progress=False)` 关闭。
- **限制**：这是 MOSES SNN 的 DNA k-mer 指纹适配，不是 Morgan-fingerprint SNN。哈希碰撞、`k`、位数和反向互补折叠会影响结果；两条都没有可用 k-mer 的短序列按两个空指纹相同处理。SNN 高也可能表示训练集记忆，应结合 Novelty、Uniqueness 和 FCD 一起解释。
- **依据**：[MOSES 对 SNN 的定义](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full)。

## 8) `EVAL-002` Ambiguity

- **作用：** 统计每条序列中 IUPAC 模糊碱基的数量、位置和比例，并按配置判断是否超限，用于量化数据的不确定碱基负担。
- **API**：`dnakit.evaluation.evaluate_ambiguity(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.AmbiguityEvaluationConfig`。
- **输入**：一条或一组 DNA；可配置最大比例、符号权重和 Gap 分母策略。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_ambiguity

report = evaluate_ambiguity(DNASequence("ACNT", alphabet="iupac"))
entry = report.entries[0]
print(entry.metrics["ambiguity_count"])     # 1
print(entry.metrics["ambiguity_fraction"])  # 0.25
```

- **示例结果：**

```text
1
0.25
```

- **限制**：未知长度 Gap 参与分母时，比例可能明确返回 `None`。

## 9) `EVAL-007` Redundancy

- **作用：** 统计序列集合中完全重复和近似重复序列所占的比例，并报告近似序列对和 cluster 压缩比例，用于判断数据是否存在过度重复。
- **API**：`dnakit.evaluation.evaluate_redundancy(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.DiversityEvaluationConfig`。
- **输入**：`DNASet`；可选相似度方法和阈值。
- **示例代码**：

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.evaluation import evaluate_redundancy

records = DNASet([
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAA"), "b"),
    DNARecord(DNASequence("AAAT"), "c"),
])
report = evaluate_redundancy(records)
print(report.metrics["score"])
print(report.metrics["exact_duplicate_fraction"])
```

- **示例结果：**

```text
0.3333333333333333
0.3333333333333333
```

- **限制**：综合分数是三个透明分项的算术平均，不是训练数据质量标签。
