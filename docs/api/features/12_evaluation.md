# 常用评价指标

从合法性、唯一性、多样性、新颖性、Fréchet 表征分布距离、片段分布、最近邻相似度、模糊度和冗余度等方面综合评价 DNA 序列或数据集。

参考相似度和分布相似度在数据评价下单列于[相似度计算](09_similarity_alignment.md)。参考库相似度、新颖性等指标使用调用方提供的 `ReferenceLibrary`；参考基因组一致性使用本地 FASTA 或 PAF。novelty 不是实验结论或任务模型预测。

以下 API 括号内列出全部直接调用参数；`[必须]` 和 `[可选]` 是参数说明，不是 Python 语法。

对 DNA 序列评价，普通批量输入统一使用 `DNA(["ACGT", "ACGA", ...])`；需要保留自定义 ID 或 metadata 时，可用 `DNA` 的记录字典输入。直接传原始 `list[str]` 不是评价函数的统一接口。逐条评价的结果在 `report.entries`，整体统计在 `report.metrics`。唯一性、多样性、片段/分布距离等是集合级指标，其整体值不等于逐条分数的平均。与参考数据比较时，仍需额外提供参考库或第二组 `DNA`；数值矩阵、FASTQ 路径和 scorecard 则继续使用各自的专用输入类型。

| 输入形式 | 适用的评价函数 | 结果形态 |
| --- | --- | --- |
| 一组 `DNA` | `evaluate_validity`、`evaluate_ambiguity`、`evaluate_complexity`、`evaluate_quality`、`evaluate_synthesis_risk` | 逐条 `entries` 与汇总 `metrics` |
| 一组 `DNA` | `evaluate_uniqueness`、`evaluate_diversity`、`evaluate_kmer_diversity`、`evaluate_redundancy` | 集合级 `metrics` |
| `DNA` + `ReferenceLibrary` | `nearest_reference`、`evaluate_novelty`、`evaluate_memorization`、`evaluate_reference_similarity` | 逐条最近参考结果与汇总指标 |
| 两组 `DNA` | `evaluate_fragment_similarity`、`evaluate_snn`、`evaluate_distribution_similarity`、`evaluate_frechet_distance` | SNN 有逐条结果；其余主要是集合级结果 |
| `DNA` + 参考 FASTA/PAF | `evaluate_reference_consistency` | 逐条一致率及全体平均值 |

`evaluate_scorecard` 接收的是指标名到数值或报告的映射，不直接接收 DNA 序列。
其他只处理单条序列的操作可将 `DNA([...])` 交给[批量计算 `run_batch`](14_engineering.md#3-eng-006)，由回调逐条处理；它们不会自动改成返回集合平均值。

## 1) Validity {#1-eval-001-validity}

Validity 包含语法合法性与参考基因组一致性两个子模块。两者分别返回报告：序列通过语法检查，不意味着它与目标参考基因组一致。

### 1.1 语法合法性 {#1-1-syntax-validity}

- **作用：** 按字母表、长度、Gap、metadata 等规则检查 DNA 序列和记录，返回有效比例及逐条问题，用于在评价或建模前确认输入是否合法。
- **计算方法：** 逐条调用验证规则，使用有效记录数除以总记录数得到有效比例。
- **API：** `dnakit.evaluation.evaluate_validity(value[必须], sequence_length[可选], config[可选], limits[可选])`；`config` 使用 `dnakit.ValidationConfig`，`limits` 使用 `dnakit.evaluation.EvaluationLimits`。
- **长度要求：** `sequence_length` 为非负整数，要求每条合法序列的碱基符号数精确等于该值（不包含 `Gap` 跨度）。非空输入的所有序列都满足长度及其他验证规则时，`valid_fraction` 才为 `1.0`；不匹配会产生 `STD_SEQUENCE_LENGTH_MISMATCH`。默认 `None` 沿用配置；显式传入时覆盖 `config.sequence_length`，保留其他规则且不修改原配置。
- **输入：** `DNA`、`DNASequence`、`DNARecord`、`DNASet` 或记录迭代器；可选 `ValidationConfig` 和资源上限。
- **输出：** 返回 `EvaluationReport`；`metrics["valid_fraction"]` 是有效比例，`entries` 保存逐条合法性和问题。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_validity

sequences = DNA(["ACGT", "ACGA", "ACG"])
report = evaluate_validity(sequences, sequence_length=4)
print(report.metrics["valid_fraction"])  # 2 / 3
print([(entry.subject_id, entry.metrics["is_valid"]) for entry in report.entries])
```

- **示例结果：**

```text
0.6666666666666666
[('sequence_1', True), ('sequence_2', True), ('sequence_3', False)]
```

### 1.2 参考基因组一致性 {#1-2-reference-consistency}

- **作用：** 衡量一组 DNA 序列与指定参考基因组的匹配程度，返回每条序列的一致率及所有输入序列的平均一致率。
- **计算方法：** 为每条序列选取匹配碱基比例最高的主比对。单条一致率 = `PAF 匹配碱基数 / 该序列长度`；没有主比对时记为 `0`。总体 `consistency_rate` 是所有输入序列的一致率平均值，空输入时为 `None`。未参与比对的碱基仍包含在序列长度中。
- **API：** `dnakit.evaluation.evaluate_reference_consistency(queries[必须], paf_path[可选], *, reference_fasta[可选], reference_dir[可选], reference_name[可选], reference_version[可选], reference_sha256[可选], config[可选], verbose[可选])`；`config` 使用 `dnakit.evaluation.ReferenceConsistencyConfig`。
- **输入：** `DNA([...])` 或 DNA 记录集合。默认下载并复用 hg38 参考 FASTA；也可传入本地 `reference_fasta`，这两种方式需要 minimap2。传入已有 `paf_path` 时，需同时指定 `reference_name` 和 `reference_version`，且 PAF 的序列 ID 和长度须与输入一致。
- **输出：** 返回 `EvaluationReport`。`report.metrics["consistency_rate"]` 是全部输入序列的平均一致率；`report.entries` 中每条记录的 `metrics["consistency_rate"]` 是该序列的一致率。
- **示例代码：** 以下 PAF 是演示数据，用来展示平均值的计算。

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNA
from dnakit.evaluation import evaluate_reference_consistency

sequences = DNA(["ACGTACGT", "TGCATGCA", "CCCCCCCC"])
with TemporaryDirectory() as folder:
    paf = Path(folder) / "example.paf"
    paf.write_text(
        "sequence_1\t8\t0\t8\t+\tchr1\t100\t20\t28\t8\t8\t60\ttp:A:P\tcg:Z:8M\n"
        "sequence_2\t8\t0\t8\t-\tchr2\t100\t30\t38\t7\t8\t30\ttp:A:P\tcg:Z:8M\n",
        encoding="utf-8",
    )
    report = evaluate_reference_consistency(
        sequences, paf, reference_name="example", reference_version="demo-1",
        verbose=False,
    )

print(report.metrics["consistency_rate"])
print([(entry.subject_id, entry.metrics["consistency_rate"]) for entry in report.entries])
```

- **示例结果：**

```text
0.625
[('sequence_1', 1.0), ('sequence_2', 0.875), ('sequence_3', 0.0)]
```

## 2) Uniqueness {#2-eval-005-uniqueness}

- **作用：** 按精确、反向互补、IUPAC 或近似等价规则对序列分组，返回唯一序列数量、比例及重复组，用于衡量集合中独立记录的占比。
- **计算方法：** 按等价规则去重；唯一记录数除以总记录数得到唯一性分数。
- **API：** `dnakit.evaluation.evaluate_uniqueness(value[必须], config[可选], plot[可选])`；`config` 使用 `dnakit.evaluation.UniquenessEvaluationConfig`。
- **输入：** `DNA([...])`；也支持 `DNASet` 等记录集合。可选等价规则、方法和阈值。
- **输出：** 返回 `EvaluationReport`；`metrics` 包含唯一性分数和重复组。
- **可选绘图：** `plot=True` 在当前工作目录保存 `uniqueness_ratio.png`（默认 `False`，同名文件会覆盖）。默认等价规则为 `exact`，图题为“Exact uniqueness”；其他等价规则的图题会标明实际规则。绘图需要 `dnakit[bio-plots]`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_uniqueness

records = DNA(["AAAA", "AAAA", "CCCC"])
report = evaluate_uniqueness(records)
print(report.metrics["uniqueness_score"])  # 0.666666...
print(report.metrics["duplicate_groups"])  # (("sequence_1", "sequence_2"),)
```

- **示例结果：**

```text
0.6666666666666666
(('sequence_1', 'sequence_2'),)
```

## 3) Diversity {#3-eval-006-diversity}

- **作用：** 按论文公式计算全部序列对的平均原始 Levenshtein 距离。
- **计算方法：** `Σ(i≠j) Levenshtein(xᵢ,xⱼ) / [n(n−1)]`，等价于所有无序序列对的原始编辑距离均值；少于两条序列时返回 `None`。
- **API：** `dnakit.evaluation.evaluate_diversity(value[必须], config[可选], plot[可选])`；`config` 使用 `dnakit.evaluation.DiversityEvaluationConfig`。
- **输入：** `DNA([...])` 或其他 DNA 记录集合；可设置 `show_progress=True` 显示两两比较进度。默认最多使用 16 个进程并行计算；可通过 `DiversityEvaluationConfig(workers=8)` 调整，或设为 `workers=1` 使用单进程。结果不受进程数影响。
- **输出：** 返回 `EvaluationReport`；`metrics` 包含平均成对 Levenshtein 距离。
- **可选绘图：** `plot=True` 在当前工作目录保存 `diversity_mean_bar.png` 和 `diversity_pairwise_distance_histogram.png`（默认 `False`，同名文件会覆盖）。直方图统计所有无序序列对；绘图需要 `dnakit[bio-plots]`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_diversity

records = DNA(["AAAA", "AAAT", "CCCC"])
report = evaluate_diversity(records)
print(report.metrics["mean_pairwise_levenshtein_distance"])
```

默认即为多进程。调整进程数时，从 `dnakit.evaluation` 导入 `DiversityEvaluationConfig`，例如 `evaluate_diversity(records, config=DiversityEvaluationConfig(workers=8))`。纯 Python 编辑距离受 GIL 限制，因此这里使用多进程。

- **示例结果：**

```text
3
```

[FAQ：补充说明](../../faq/tutorial_basics.md#levenshtein-3)

需要规范化 k-mer 的 Jaccard 距离时，使用 `evaluate_kmer_diversity(records, k=6)`。返回报告中的 `mean_pair_distance` 是所有无序序列对的平均 `1 − Jaccard`，`mean_nearest_neighbor_distance` 是每条序列到集合内其他序列的最近距离均值；少于两条序列时两者均为 `None`。默认使用正反互补规范化的 k-mer，并遵守 `EvaluationLimits` 的成对比较上限。它与上述原始 Levenshtein 距离是不同指标。

```python
from dnakit.evaluation import evaluate_kmer_diversity

report = evaluate_kmer_diversity(records, k=6)
print(report.metrics["mean_nearest_neighbor_distance"])
```

## 4) Novelty {#4-eval-008-novelty}

- **作用：** 并列报告连续的新颖性值与完整序列精确一致比例。默认每条序列的 `novelty = 1 - 最近参考相似度`，汇总 `score` 为其平均值；`exact_match_rate` 是与参考库任意一条完整记录逐字相同的查询比例。第二种方法按论文公式报告每条查询到参考库的最小原始 Levenshtein 距离及其平均值，同时报告精确一致比例。
- **计算方法：** 精确一致独立检查所有参考记录，不依赖最近相似序列的排序；不按 `0.9` 等阈值生成 `is_novel` 或 `novel_fraction`。`copy_threshold` 只用于记忆/近复制评价，对 Novelty 的两个结果无效。
- **API：** `dnakit.evaluation.evaluate_novelty(queries[必须], reference[必须], config[可选], plot[可选])`；`config` 使用 `dnakit.evaluation.ReferenceSearchConfig`。
- **输入：** `DNA([...])` 查询集合和非空 `ReferenceLibrary`；默认 `novelty_calculation="similarity"`，论文方法使用 `novelty_calculation="levenshtein"`；可设置 `show_progress=True`。
- **输出：** 返回 `EvaluationReport`；`metrics` 包含平均新颖性分数和精确匹配率。
- **可选绘图：** `plot=True` 在当前工作目录保存 `novelty_batch_summary.png`（平均新颖性与精确匹配率）、`novelty_distribution.png`（逐条新颖性分布）和 `novelty_nearest_levenshtein.png`（最近参考原始编辑距离分布）；默认 `False`，同名文件会覆盖。查询记录的 `metadata["batch"]` 若为非空文本，汇总柱状图会按批次分组；否则显示一个查询集。相似度模式绘制第三张图时会额外计算 Levenshtein 距离。绘图需要 `dnakit[bio-plots]`。Levenshtein 模式下，新颖性本身就是最近参考原始编辑距离，因此后两张直方图展示同一组距离。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import ReferenceSearchConfig, create_reference_library, evaluate_novelty

reference_records = DNA(["AAAA", "CCCC"])
reference = create_reference_library(
    reference_records,
    name="training",
    version="1",
    source="local:example",
)
queries = DNA(["AAAA", "GGGG"])

report = evaluate_novelty(
    queries,
    reference,
    config=ReferenceSearchConfig(method="identity"),
)
print(report.metrics["score"])             # 平均 novelty：0.5
print(report.metrics["exact_match_rate"])  # 完整序列精确一致比例：0.5
```

- **示例结果：**

```text
0.5
0.5
```

[FAQ：补充说明](../../faq/tutorial_basics.md#levenshtein-4)

## 5) Fréchet DNA distance {#5-eval-016-frechet-dna-distance}

- **作用：** 使用同一个 DNA 基础模型分别表示两个序列集合，将两组向量近似为多元高斯分布，再计算均值和协方差的 Fréchet 距离。数值越小表示两个集合在该表征空间中的分布越接近，完全相同的表征分布趋近于 `0`。
- **计算方法：** 计算两组模型向量的均值差平方与协方差 Fréchet 项；采用等价的样本空间 cross-Gram 核范数实现。
- **API：** `dnakit.evaluation.evaluate_frechet_distance(left[必须], right[必须], config[可选], backend[可选], plot_distance[可选], plot_pca[可选], plot_output_dir[可选], plot_overwrite[可选])`；两个绘图开关默认均为 `False`。`config` 使用 `dnakit.evaluation.FrechetDistanceConfig`，其中的表征配置使用 `dnakit.representations.RepresentationConfig`。
- **输入：** 两组各含至少 2 条序列的 `DNA([...])`，也支持其他记录集合；此指标默认使用 `lucaone`、mean pooling 和 L2 归一化。可选择其他已注册模型、checkpoint、pooling、设备、dtype 和批大小。
- **输出：** 返回 `EvaluationReport`；`metrics["frechet_distance"]` 是 Fréchet 距离。启用绘图后，在 `plot_output_dir`（默认当前目录）保存 `frechet_distance.png` 和/或 `frechet_pca.png`，绝对路径记录在 `parameters["plot_files"]`。已有文件默认不覆盖，可设置 `plot_overwrite=True`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import FrechetDistanceConfig, evaluate_frechet_distance
from dnakit.representations import RepresentationConfig

reference = DNA(["ACGTACGT", "AACCGGTT"])
generated = DNA(["ACGTTCGT", "AACCAGTT"])
report = evaluate_frechet_distance(
    generated,
    reference,
    config=FrechetDistanceConfig(
        representation=RepresentationConfig(allow_remote_code=True),
    ),
    plot_distance=True,
    plot_pca=True,
    plot_output_dir="frechet_plots",
)
print(report.metrics["frechet_distance"])
print(report.parameters["plot_files"])
```

- **示例结果：**

```text
非负浮点数；具体值取决于模型、checkpoint、pooling、归一化和输入集合。示例另保存两张 PNG 图。
```



[FAQ：补充说明](../../faq/tutorial_basics.md#evaluation-background-5)

## 6) Frag {#6-eval-017-frag}

- **作用：** 参考分子生成评价中的 Frag，把生成集合和参考集合分别转换为片段出现次数向量，再计算余弦相似度。数值范围为 `[0,1]`，越高表示两个集合的局部片段分布越接近。
- **计算方法：** 分别统计两组 DNA 的重叠 k-mer 次数，再计算计数向量的余弦相似度。
- **API：** `dnakit.evaluation.evaluate_fragment_similarity(generated[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.FragmentSimilarityConfig`。
- **输入：** 生成集合和参考集合均可用 `DNA([...])`；默认 k=3、canonical k-mer，跳过包含 IUPAC 模糊碱基的窗口。
- **输出：** 返回 `EvaluationReport`；`metrics["frag"]` 是片段分布分数。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import FragmentSimilarityConfig, evaluate_fragment_similarity

generated = DNA(["AAAA", "CCCC"])
reference = DNA(["AAAA", "GGGG"])
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



[FAQ：补充说明](../../faq/tutorial_basics.md#evaluation-background-6)

## 7) SNN {#7-eval-018-snn}

- **作用：** 对每条生成 DNA，在参考集合中查找二进制指纹 Tanimoto 相似度最高的序列，再对这些最近邻相似度取算术平均。数值范围为 `[0,1]`，越高表示生成样本越接近参考集合覆盖的表征空间。
- **计算方法：** 为每条查询取参考集中的最高 Tanimoto 指纹相似度，再对查询取平均；该指标通常不对称。
- **API：** `dnakit.evaluation.evaluate_snn(generated[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.SNNConfig`。
- **输入：** 生成集合和参考集合均可用 `DNA([...])`；默认 canonical 7-mer 经哈希映射为 1024 位指纹。
- **输出：** 返回 `EvaluationReport`；`metrics["snn"]` 是最近邻相似度指标。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import SNNConfig, evaluate_snn

generated = DNA(["AAAAAAA", "ATATATA"])
reference = DNA(["AAAAAAA", "CCCCCCC"])
report = evaluate_snn(
    generated,
    reference,
    config=SNNConfig(show_progress=False),
)
print(report.metrics["snn"])  # 0.5
print(report.entries[0].metrics["nearest_reference_id"])  # sequence_1
```

- **示例结果：**

```text
0.5
sequence_1
```



[FAQ：补充说明](../../faq/tutorial_basics.md#evaluation-background-7)

## 8) Ambiguity {#8-eval-002-ambiguity}

- **作用：** 统计每条序列中 IUPAC 模糊碱基的数量、位置和比例，并按配置判断是否超限，用于量化数据的不确定碱基负担。
- **计算方法：** 计数非 A/C/G/T 的 IUPAC 符号，按指定 Gap 分母策略计算比例并检查阈值。
- **API：** `dnakit.evaluation.evaluate_ambiguity(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.AmbiguityEvaluationConfig`。
- **输入：** 一条或一组 DNA，批量可用 `DNA([...], alphabet="iupac")`；可配置最大比例、符号权重和 Gap 分母策略。
- **输出：** 返回 `EvaluationReport`；逐条 `entries` 包含模糊碱基数量和比例，`metrics` 包含汇总比例。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_ambiguity

report = evaluate_ambiguity(DNA(["ACNT", "ACGT"], alphabet="iupac"))
entry = report.entries[0]
print(entry.metrics["ambiguity_count"])     # 1
print(entry.metrics["ambiguity_fraction"])  # 0.25
print(report.metrics["ambiguity_fraction_summary"]["mean"])  # 0.125
```

- **示例结果：**

```text
1
0.25
0.125
```

## 9) Redundancy {#9-eval-007-redundancy}

- **作用：** 统计序列集合中完全重复和近似重复序列所占的比例，并报告近似序列对和 cluster 压缩比例，用于判断数据是否存在过度重复。
- **计算方法：** 结合精确重复分组与阈值近似分组，汇总重复比例及簇压缩比例。
- **API：** `dnakit.evaluation.evaluate_redundancy(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.DiversityEvaluationConfig`。
- **输入：** `DNA([...])` 或其他 DNA 记录集合；可选相似度方法和阈值。
- **输出：** 返回 `EvaluationReport`；`metrics` 包含冗余分数和精确重复比例。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.evaluation import evaluate_redundancy

records = DNA(["AAAA", "AAAA", "AAAT"])
report = evaluate_redundancy(records)
print(report.metrics["score"])
print(report.metrics["exact_duplicate_fraction"])
```

- **示例结果：**

```text
0.3333333333333333
0.3333333333333333
```

<span id="eval-008-novelty"></span>

## 10) 生成 FastQC 质量报告 {#bio-19}

- **作用：** 查看 reads 的质量分布和测序质控报告。
- **计算方法：** 对提供的 FASTQ 文件运行 FastQC。
- **API：** `dnakit.bioinformatics.ngs.quality_control(reads, output_directory, config=None)`
- **输入：** FASTQ 路径列表；FastQC 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import quality_control

result = quality_control(["R1.fq.gz", "R2.fq.gz"], "new_qc")
print(sorted(result.artifacts))
```

- **示例结果：**

返回 `ToolResult`，`artifacts` 包含各输入文件对应的 HTML/ZIP 报告。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

## 11) 批量计算序列性质 {#bio-26}

- **作用：** 用自定义预测函数生成后续评价需要的性质矩阵。
- **计算方法：** 对每条序列依次调用命名的标量 predictor，保持 ID 和列名顺序。
- **API：** `dnakit.bioinformatics.metrics.sequence_properties(sequences, predictors)`
- **输入：** `DNA([...])` 或 ID→序列字典，以及名称→Python 标量函数字典；需要 `bio-analysis`。`DNA` 中的记录 ID 和顺序会保留在结果中。
- **输出：** 返回字典，包含序列 ID、性质名称和逐条性质值矩阵。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.bioinformatics.metrics import sequence_properties

result = sequence_properties(DNA(["ACGT", "GGGG"]),
                             {"gc": lambda s: sum(x in "GC" for x in s) / len(s)})
print(result["ids"], result["names"], result["values"].tolist())
```

- **示例结果：**

```text
('sequence_1', 'sequence_2') ('gc',) [[0.5], [1.0]]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 12) 计算 MMD² 或 KID 分布距离 {#bio-27}

- **作用：** 比较生成样本和参考样本的 embedding 分布。
- **计算方法：** MMD 使用 Gaussian kernel；KID 使用三次多项式 kernel，默认无偏估计。
- **API：** `dnakit.bioinformatics.metrics.kernel_distance(generated, reference, metric='mmd', sigma=1.0, unbiased=True)`
- **输入：** 两组列数相同的 embedding；无偏估计可为负；需要 `bio-analysis`。
- **输出：** 返回 MMD² 或 KID 的浮点数值。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import kernel_distance

values = [[0, 0], [1, 1], [2, 0]]
print(kernel_distance(values, values, unbiased=False))
```

- **示例结果：**

```text
0.0
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 13) 计算 PRDC 精确率、召回率、密度与覆盖率 {#bio-28}

- **作用：** 评价生成样本与参考分布的局部覆盖关系。
- **计算方法：** 用 k 近邻确定半径，以严格开球统计两组点的覆盖。
- **API：** `dnakit.bioinformatics.metrics.manifold_metrics(generated, reference, neighbors=5)`
- **输入：** 同维 embedding；每组样本数大于 `neighbors`；固定 embedding 模型与缩放；需要 `bio-analysis`。
- **输出：** 返回字典，包含 `precision`、`recall`、`density` 和 `coverage`。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import manifold_metrics

values = [[0], [1], [3], [6]]
print(manifold_metrics(values, values, neighbors=1))
```

- **示例结果：**

```text
{'precision': 1.0, 'recall': 1.0, 'density': 1.0, 'coverage': 1.0}
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 14) 计算 clipped density 与 coverage {#bio-29}

- **作用：** 计算带截断和有限样本校准的覆盖指标。
- **计算方法：** 使用 median-clipped 半径、leave-one-out 密度和 beta-binomial coverage 校准。
- **API：** `dnakit.bioinformatics.metrics.clipped_manifold_metrics(generated, reference, neighbors=5)`
- **输入：** 两组样本数相同、维度相同且每组大于 `neighbors` 的 embedding；需要 `bio-analysis`。
- **输出：** 返回字典，包含 `clipped_density` 和 `clipped_coverage`。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import clipped_manifold_metrics

values = [[0], [1], [3], [6]]
print(clipped_manifold_metrics(values, values, neighbors=1))
```

- **示例结果：**

```text
{'clipped_density': 1.0, 'clipped_coverage': 1.0}
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 15) 检查生成样本是否过于接近训练集 {#bio-30}

- **作用：** 估计生成样本中未过度贴近训练样本的比例。
- **计算方法：** 比较生成点到最近训练点的距离与该训练点到其他训练点的最近距离。
- **API：** `dnakit.bioinformatics.metrics.authenticity(generated, training)`
- **输入：** 生成 embedding 和生成器实际训练集 embedding；至少两个训练样本；需要 `bio-analysis`。
- **输出：** 返回真实性分数，数值类型为 `float`。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import authenticity

print(authenticity([[0], [20]], [[0], [1], [2]]))
```

- **示例结果：**

```text
0.5
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 16) 计算 Pareto 层级和超体积 {#bio-31}

- **作用：** 比较同时优化多个目标的候选结果。
- **计算方法：** 按各列目标方向逐层非支配排序，并计算相对较差参考点的精确超体积。
- **API：** `dnakit.bioinformatics.metrics.pareto_ranks(values, maximize=...)`; `dnakit.bioinformatics.metrics.hypervolume(values, reference_point, maximize=...)`
- **输入：** samples×objectives 数组、各列最大化/最小化方向；超体积最多 5 维、1000 点；需要 `bio-analysis`。
- **输出：** `pareto_ranks` 返回各样本层级；`hypervolume` 返回超体积数值。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import pareto_ranks, hypervolume

values = [[1, 1], [2, 2], [3, 1]]
print(pareto_ranks(values, maximize=[True, True]).tolist())
print(hypervolume(values, [0, 0], maximize=[True, True]))
```

- **示例结果：**

```text
[1, 0, 0]
5.0
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 17) 统计目标性质命中率 {#bio-32}

- **作用：** 查看生成样本有多少同时满足各性质范围。
- **计算方法：** 按各列闭区间判断，返回联合命中率、逐列命中率和 mask。
- **API：** `dnakit.bioinformatics.metrics.property_hit_rate(values, lower=..., upper=...)`
- **输入：** 性质矩阵和每列上下界；需要 `bio-analysis`。
- **输出：** 返回字典，包含联合命中率和各性质命中率。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import property_hit_rate

result = property_hit_rate([[0.5, 10], [0.7, 15]], lower=[0.4, 5], upper=[0.6, 20])
print(result["joint"])
print(sorted(result))
```

- **示例结果：**

```text
0.5
['joint', 'mask', 'per_property']
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 18) 计算性质分布 KL 散度 {#bio-33}

- **作用：** 比较生成与参考数据的各列性质分布。
- **计算方法：** 共享直方图边界并加 pseudocount，计算自然对数 KL(generated∥reference)。
- **API：** `dnakit.bioinformatics.metrics.property_kl(generated, reference, bins=20, pseudocount=1e-8)`
- **输入：** 两组同维性质矩阵；分箱和 pseudocount 在比较前固定；需要 `bio-analysis`。
- **输出：** 返回各性质的 KL 散度数值。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import property_kl

values = [[0], [1], [2], [3]]
print(property_kl(values, values, bins=2).tolist())
```

- **示例结果：**

```text
[0.0]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。

## 19) 计算分布一致性分数 {#bio-34}

- **作用：** 使用参考分布密度评价生成样本的一致性。
- **计算方法：** 固定带宽 Gaussian KDE，交叉拟合参考集，再按验证密度经验分布评分。
- **API：** `dnakit.bioinformatics.metrics.conformity_score(generated, reference, folds=5, bandwidth=1.0, seed=0)`
- **输入：** 生成和独立参考 embedding；样本数须足够划分 folds；需要 `bio-analysis`。此分数不是类别概率。
- **输出：** 返回字典，包含一致性分数、标准误和交叉验证分数。
- **示例代码：**

```python
from dnakit.bioinformatics.metrics import conformity_score

result = conformity_score([[0], [1], [2]], [[0], [1], [2], [3], [4], [5]], folds=3)
print(sorted(result))
print(round(result["value"], 6))
```

- **示例结果：**

```text
['fold_scores', 'standard_error', 'value']
0.62963
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#metrics)。
