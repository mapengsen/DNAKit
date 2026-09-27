# DNA 描述符

计算 DNA 序列的长度、碱基组成、GC/AT、偏斜、CpG、k-mer、熵、复杂度、重复和密码子等基础描述符。 [FAQ：详细说明](../../faq/tutorial_integrations.md#descriptors)。

## 1) 模糊碱基统计 {#1-std-005}

- **作用：** 统计 `N` 及其他 IUPAC 模糊碱基的数量、比例和具体位置，用于评估序列确定性并按阈值筛选数据。
- **处理方法：** 按 IUPAC 字母表识别模糊符号，记录数量、比例和位置。
- **API：** `dnakit.normalize(raw[必须], keep_ambiguous[可选], keep_u[可选], keep_other[可选], config[可选])`、`dnakit.validate(sequence[必须], config[可选])`；结果字段为 `ambiguity`（`AmbiguityReport`）。
- **输入：** 必填 IUPAC DNA；可选比例分母是否包含 Gap，以及字母表和模糊碱基策略。
- **输出：** 返回 `NormalizationResult`；`ambiguity` 包含模糊碱基数量、比例和各符号统计。
- **示例代码：**

```python
from dnakit import normalize

result = normalize("ANRY")
print(result.ambiguity.total_count)
print(result.ambiguity.fraction)
print([(item.symbol, item.count) for item in result.ambiguity.by_symbol])
```

- **示例结果：**

```text
3
0.75
[('N', 1), ('R', 1), ('Y', 1)]
```

## 2) 长度特征 {#2-desc-001}

- **作用：** 计算 DNA 序列的碱基长度，并区分符号数量与 Gap 对坐标跨度的影响，用于长度筛选、分箱和窗口参数设置。
- **处理方法：** 分别统计序列符号数和包含 Gap 的坐标跨度。
- **API：** `dnakit.descriptors.length_features(value[必须])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`。
- **输出：** 返回 `LengthResult`，示例展示 `canonical_base_count`。
- **示例代码：**

```python
import dnakit
from dnakit.descriptors import length_features

seq = dnakit.normalize("ACGTACGT").sequence
assert seq is not None
result = length_features(seq)
print(len(seq))                     # 8
print(result.canonical_base_count)  # 8
```

- **示例结果：**

```text
8
8
```

## 3) 碱基组成 {#3-desc-002}

- **作用：** 统计 DNA 序列中 A、C、G、T 的数量和比例，形成基础组成向量，供 GC、偏斜、数据分布比较和建模使用。
- **处理方法：** 逐位统计 A/C/G/T，并按所选模糊碱基策略确定比例分母。
- **API：** `dnakit.descriptors.base_composition(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选 `ambiguity_policy="error"|"ignore"`。
- **输出：** 返回 `CompositionResult`，示例展示 `counts`、`fractions`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_composition

result = base_composition(DNASequence("AACGT"))
print(dict(result.counts), result.fractions["A"])
```

- **示例结果：**

```text
{'A': 2, 'C': 1, 'G': 1, 'T': 1} 0.4
```

## 4) GC/AT特征 {#4-desc-003-gcat}

- **作用：** 计算 DNA 序列的 GC、AT 数量及比例，用于比较序列组成、评估扩增或测序偏好，并作为其他描述符的基础输入。
- **处理方法：** 分别求 G+C、A+T 的计数，再除以有效碱基数。
- **API：** `dnakit.descriptors.gc_at_content(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选模糊碱基处理策略。
- **输出：** 返回 `ContentResult`，示例展示 `gc_fraction`、`at_fraction`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import gc_at_content

result = gc_at_content(DNASequence("AACG"))
print(result.gc_fraction, result.at_fraction)
```

- **示例结果：**

```text
0.5 0.5
```

## 5) 碱基偏斜 {#5-desc-004}

- **作用：** 分别计算 `(G-C)/(G+C)` 和 `(A-T)/(A+T)`，量化互补碱基在当前链上的不对称程度，用于观察局部或整体组成偏向。
- **处理方法：** 计算 `(G-C)/(G+C)` 与 `(A-T)/(A+T)`。
- **API：** `dnakit.descriptors.base_skew(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选模糊碱基处理策略。
- **输出：** 返回 `SkewResult`，示例展示 `gc_skew`、`at_skew`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import base_skew

result = base_skew(DNASequence("AACG"))
print(result.gc_skew, result.at_skew)
```

- **示例结果：**

```text
0.0 1.0
```

## 6) CpG特征 {#6-desc-005-cpg}

- **作用：** 统计 CpG 位点数量、密度及观测值/期望值，量化 CpG 富集或缺失程度，供 CpG 区域筛查和集合比较使用。
- **处理方法：** 统计相邻 CpG，结合有效长度及 C/G 计数计算密度与观测/期望比。
- **API：** `dnakit.descriptors.cpg_features(value[必须], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选 `ambiguity_policy`、`cross_gaps`。
- **输出：** 返回 `CpGResult`，示例展示 `cpg_count`、`density`、`observed_expected`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import cpg_features

result = cpg_features(DNASequence("ACGCGT"))
print(result.cpg_count, result.density, result.observed_expected)
```

- **示例结果：**

```text
2 0.4 3.0
```

## 7) k-mer统计 {#7-desc-006-k-mer}

- **作用：** 枚举指定长度的 k-mer，返回数量、频率和存在性，可选择合并反向互补 k-mer，用于组成分析、相似度计算和特征建模。
- **处理方法：** 扫描长度为 k 的窗口；按配置允许重叠或合并反向互补词。
- **API：** `dnakit.descriptors.kmer_statistics(value[必须], k[必须], overlapping[可选], canonical[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord` 和 `k`；可选重叠、canonical、模糊碱基及 Gap 策略。
- **输出：** 返回 `KmerResult`，示例展示 `counts`、`denominator`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import kmer_statistics

result = kmer_statistics(DNASequence("ACGT"), 2)
print(dict(result.counts), result.denominator)
```

- **示例结果：**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 3
```

## 8) 序列熵 {#8-desc-007}

- **作用：** 根据碱基或 k-mer 的概率分布计算 Shannon entropy，量化分布均匀程度，用于识别信息量较低或组成单一的序列。
- **处理方法：** 按碱基或 k-mer 频率计算 `−Σ p log(p)`。
- **API：** `dnakit.descriptors.shannon_entropy(value[必须], unit[可选], k[可选], log_base[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 `unit`、`k`、`log_base`、模糊碱基及 Gap 策略。
- **输出：** 返回 `EntropyResult`，示例展示 `entropy`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import shannon_entropy

result = shannon_entropy(DNASequence("ACGT"))
print(result.entropy)
```

- **示例结果：**

```text
2.0
```

## 9) 序列复杂度 {#9-desc-008}

- **作用：** 比较不同长度短序列的实际种类数与理论可出现种类数，得到 linguistic complexity，用于发现重复较多或模式种类不足的序列。
- **处理方法：** 对各词长计算观察到的词种数与可能词种数之比，再取乘积。
- **API：** `dnakit.descriptors.linguistic_complexity(value[必须], max_word_size[可选], ambiguity_policy[可选], cross_gaps[可选], max_observations[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选最大词长、模糊碱基、Gap 和工作量上限。
- **输出：** 返回 `ComplexityResult`；`score` 是序列复杂度分数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import linguistic_complexity

simple = linguistic_complexity(DNASequence("AAAAAAAA"), max_word_size=3)
diverse = linguistic_complexity(DNASequence("ACGTAGCT"), max_word_size=3)
print(simple.score < diverse.score)
```

- **示例结果：**

```text
True
```

## 10) Homopolymer {#10-desc-009-homopolymer}

- **作用：** 查找连续重复的同一种碱基，返回各区段的位置、碱基和长度，并汇总最长区段，用于筛查测序、合成和扩增风险。
- **处理方法：** 连续扫描相同碱基，保留达到最短长度的区段。
- **API：** `dnakit.descriptors.homopolymer_runs(value[必须], min_run_length[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选最短连续长度、模糊碱基及 Gap 策略。
- **输出：** 返回 `HomopolymerResult`，示例展示 `longest_length`、`runs`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import homopolymer_runs

result = homopolymer_runs(DNASequence("AAACCGTTTT"), min_run_length=2)
print(result.longest_length, result.runs[-1].base)
```

- **示例结果：**

```text
4 T
```

## 11) 重复比例 {#11-desc-010}

- **作用：** 查找相邻重复的序列单元，返回重复区段和覆盖比例，用于量化高重复区域以及筛选可能影响比对或合成的序列。
- **处理方法：** 识别精确串联重复，将覆盖区间取并集后除以有效长度。
- **API：** `dnakit.descriptors.exact_repeat_fraction(value[必须], min_unit_length[可选], max_unit_length[可选], min_repeats[可选], ambiguity_policy[可选], cross_gaps[可选], max_comparisons[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选重复单元长度、最少次数、Gap 策略和比较上限。
- **输出：** 返回 `ExactRepeatResult`，示例展示 `repeat_fraction`、`runs`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import exact_repeat_fraction

result = exact_repeat_fraction(DNASequence("ATATATGC"), min_unit_length=2)
print(result.repeat_fraction, result.runs[0].unit)
```

- **示例结果：**

```text
0.75 AT
```

## 12) 窗口描述符 {#12-desc-011}

- **作用：** 按窗口长度和步长切分 DNA，逐窗计算 GC、entropy 等描述符并保留位置，用于观察序列特征沿坐标的局部变化。
- **处理方法：** 按窗口长度与步长切片，在每个窗口计算所选描述符。
- **API：** `dnakit.descriptors.window_descriptors(value[必须], descriptors[必须], window_size[必须], step[可选], include_partial[可选], entropy_log_base[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填序列、描述符列表和 `window_size`；可选步长、`include_partial`、熵底数、模糊碱基和 Gap 策略。
- **输出：** 返回 `WindowDescriptorResult`，示例展示 `windows`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import window_descriptors

result = window_descriptors(
    DNASequence("ACGT"), ["gc", "entropy", "cpg"], window_size=2, step=2
)
print(result.windows[0].symbol_start, result.windows[0].values["gc_fraction"])
```

- **示例结果：**

```text
0 0.5
```

## 13) 密码子统计 {#13-desc-012}

- **作用：** 按指定阅读框统计密码子及起始、终止密码子，返回计数和位置，用于分析编码组成、阅读框和密码子使用情况。
- **处理方法：** 从指定阅读框按三联体计数，并用遗传密码表标记起始和终止密码子。
- **API：** `dnakit.descriptors.codon_statistics(value[必须], frame[可选], genetic_code[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 `frame=0|1|2`、`genetic_code=1`、模糊碱基和 Gap 策略。
- **输出：** 返回 `CodonResult`，示例展示 `codon_count`、`start_count`、`stop_count`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import codon_statistics

result = codon_statistics(DNASequence("ATGAAATAA"))
print(result.codon_count, result.start_count, result.stop_count)
```

- **示例结果：**

```text
3 1 1
```

---


## 14) 全部描述符计算（240 项） {#all-descriptors}

- **作用：** 一次得到固定顺序的组成、统计和理化特征向量。
- **计算方法：** 按 `descriptor_schema_v1` 计算；不适用的字段返回 `None` 并记录原因。
- **API：** `dnakit.descriptors.all_descriptors(value, *, ambiguity_policy="ignore", conditions=None, dinucleotide_property_table=None)`
- **输入：** `DNASequence` 或 `DNARecord`；后 60 项二核苷酸属性需要用户提供参数表。各项适用范围与默认条件见[FAQ：详细说明](../../faq/tutorial_integrations.md#descriptors)。
- **输出：** 返回 `AllDescriptorsResult`，示例展示 `values`、`unavailable_reasons`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

- **示例结果：**

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

<!-- BEGIN GENERATED API: sequence-descriptors -->

<span id="sequence-descriptors"></span>

## 15) 统计指定碱基或氨基酸的数量与比例 {#dnakit.databases.sequence_count_residues}

- **作用：** 统计指定碱基或氨基酸的数量与比例。
- **处理方法：** 逐位统计目标残基，输出数量、比例及从 1 开始的位置。
- **API：** `dnakit.databases.sequence_count_residues(*, operation=None, sequence=None, uniprot_id=None, residue, verbose=True)`
- **输入：** 必须提供 `residue`。 另需提供 `sequence` 或 `uniprot_id`；使用 UniProt 编号时需要联网。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sequence_count_residues(
    operation='count_residues',
    uniprot_id='P00533',
    residue='C',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "count": "<integer>",
    "sequence_length": "<integer>",
    "fraction": "<number>",
    "percent": "<number>",
    "positions_1based": [
      "<value>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_count_residues-parameters)

## 16) 计算核酸序列的 GC 含量 {#dnakit.databases.sequence_gc_content}

- **作用：** 计算核酸序列的 GC 含量。
- **处理方法：** 统计 G、C 数量，除以序列长度减 N 数量得到 GC 比例。
- **API：** `dnakit.databases.sequence_gc_content(*, operation=None, sequence, verbose=True)`
- **输入：** 必须提供 `sequence`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sequence_gc_content(
    operation='gc_content',
    sequence='ATGCGATCGATCGATCGATCG',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "gc_fraction": "<number>",
    "gc_percent": "<number>",
    "length": "<integer>",
    "composition": "<object>"
  }
}
```

[完整参数与使用说明](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_gc_content-parameters)

## 17) 统计序列长度、组成与类型特征 {#dnakit.databases.sequence_stats}

- **作用：** 统计序列长度、组成与类型特征。
- **处理方法：** 统计序列长度和字符组成，再按核酸或蛋白类型计算对应性质。
- **API：** `dnakit.databases.sequence_stats(*, operation=None, sequence=None, uniprot_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 另需提供 `sequence` 或 `uniprot_id`；使用 UniProt 编号时需要联网。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sequence_stats(
    operation='stats',
    sequence='MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAK',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "sequence_type": "<string>",
    "length": "<integer>",
    "composition": "<object>"
  }
}
```

[完整参数与使用说明](../../faq/api/05_all_descriptors.md#dnakit.databases.sequence_stats-parameters)

<!-- END GENERATED API: sequence-descriptors -->
