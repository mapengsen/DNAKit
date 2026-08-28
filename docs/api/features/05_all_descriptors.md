# DNA 描述符

计算 DNA 序列的长度、碱基组成、GC/AT、偏斜、CpG、k-mer、熵、复杂度、重复和密码子等基础描述符。

结果对象可通过 `to_dict()` 转为可序列化字典。凡接口提供 `cross_gaps`，`True` 也只允许跨越可跨 Gap；显式设置 `Gap(crossable=False)` 的位置始终是硬边界。接受环状序列的描述符仍以当前 origin 为扫描边界，不自动补算末端到开头的相邻对或窗口；需要跨原点语义时应先旋转或显式线性化序列。

以下先列出 `STD-005` 和 `DESC-001`–`DESC-012` 共 13 个可单独调用的基础功能；需要一次输出更完整的固定特征向量时，请继续查看本页后面的[“全部描述符计算”](#all-descriptors)。

## 1) STD-005 模糊碱基统计

- **作用：** 统计 `N` 及其他 IUPAC 模糊碱基的数量、比例和具体位置，用于评估序列确定性并按阈值筛选数据。
- **API：** `dnakit.normalize(raw[必须], keep_ambiguous[可选], keep_u[可选], keep_other[可选], config[可选])`、`dnakit.validate(sequence[必须], config[可选])`；结果字段为 `ambiguity`（`AmbiguityReport`）。
- **输入：** 必填 IUPAC DNA；可选比例分母是否包含 Gap，以及字母表和模糊碱基策略。
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

## 2) DESC-001 长度特征

- **作用：** 计算 DNA 序列的碱基长度，并区分符号数量与 Gap 对坐标跨度的影响，用于长度筛选、分箱和窗口参数设置。
- **API：** `dnakit.descriptors.length_features(value[必须])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`。
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

## 3) DESC-002 碱基组成

- **作用：** 统计 DNA 序列中 A、C、G、T 的数量和比例，形成基础组成向量，供 GC、偏斜、数据分布比较和建模使用。
- **API：** `dnakit.descriptors.base_composition(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选 `ambiguity_policy="error"|"ignore"`。
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

## 4) DESC-003 GC/AT特征

- **作用：** 计算 DNA 序列的 GC、AT 数量及比例，用于比较序列组成、评估扩增或测序偏好，并作为其他描述符的基础输入。
- **API：** `dnakit.descriptors.gc_at_content(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选模糊碱基处理策略。
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

## 5) DESC-004 碱基偏斜

- **作用：** 分别计算 `(G-C)/(G+C)` 和 `(A-T)/(A+T)`，量化互补碱基在当前链上的不对称程度，用于观察局部或整体组成偏向。
- **API：** `dnakit.descriptors.base_skew(value[必须], ambiguity_policy[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选模糊碱基处理策略。
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

## 6) DESC-005 CpG特征

- **作用：** 统计 CpG 位点数量、密度及观测值/期望值，量化 CpG 富集或缺失程度，供 CpG 区域筛查和集合比较使用。
- **API：** `dnakit.descriptors.cpg_features(value[必须], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence` 或 `DNARecord`；可选 `ambiguity_policy`、`cross_gaps`。
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

## 7) DESC-006 k-mer统计

- **作用：** 枚举指定长度的 k-mer，返回数量、频率和存在性，可选择合并反向互补 k-mer，用于组成分析、相似度计算和特征建模。
- **API：** `dnakit.descriptors.kmer_statistics(value[必须], k[必须], overlapping[可选], canonical[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord` 和 `k`；可选重叠、canonical、模糊碱基及 Gap 策略。
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

## 8) DESC-007 序列熵

- **作用：** 根据碱基或 k-mer 的概率分布计算 Shannon entropy，量化分布均匀程度，用于识别信息量较低或组成单一的序列。
- **API：** `dnakit.descriptors.shannon_entropy(value[必须], unit[可选], k[可选], log_base[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 `unit`、`k`、`log_base`、模糊碱基及 Gap 策略。
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

## 9) DESC-008 序列复杂度

- **作用：** 比较不同长度短序列的实际种类数与理论可出现种类数，得到 linguistic complexity，用于发现重复较多或模式种类不足的序列。
- **API：** `dnakit.descriptors.linguistic_complexity(value[必须], max_word_size[可选], ambiguity_policy[可选], cross_gaps[可选], max_observations[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选最大词长、模糊碱基、Gap 和工作量上限。
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

## 10) DESC-009 Homopolymer

- **作用：** 查找连续重复的同一种碱基，返回各区段的位置、碱基和长度，并汇总最长区段，用于筛查测序、合成和扩增风险。
- **API：** `dnakit.descriptors.homopolymer_runs(value[必须], min_run_length[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选最短连续长度、模糊碱基及 Gap 策略。
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

## 11) DESC-010 重复比例

- **作用：** 查找相邻重复的序列单元，返回重复区段和覆盖比例，用于量化高重复区域以及筛选可能影响比对或合成的序列。
- **API：** `dnakit.descriptors.exact_repeat_fraction(value[必须], min_unit_length[可选], max_unit_length[可选], min_repeats[可选], ambiguity_policy[可选], cross_gaps[可选], max_comparisons[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选重复单元长度、最少次数、Gap 策略和比较上限。
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

## 12) DESC-011 窗口描述符

- **作用：** 按窗口长度和步长切分 DNA，逐窗计算 GC、entropy 等描述符并保留位置，用于观察序列特征沿坐标的局部变化。
- **API：** `dnakit.descriptors.window_descriptors(value[必须], descriptors[必须], window_size[必须], step[可选], include_partial[可选], entropy_log_base[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填序列、描述符列表和 `window_size`；可选步长、`include_partial`、熵底数、模糊碱基和 Gap 策略。
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

## 13) DESC-012 密码子统计

- **作用：** 按指定阅读框统计密码子及起始、终止密码子，返回计数和位置，用于分析编码组成、阅读框和密码子使用情况。
- **API：** `dnakit.descriptors.codon_statistics(value[必须], frame[可选], genetic_code[可选], ambiguity_policy[可选], cross_gaps[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 `frame=0|1|2`、`genetic_code=1`、模糊碱基和 Gap 策略。
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



## 14) 全部描述符计算（240项）

**完整 240 项描述符**

一次性计算 DNA 序列的组成、统计、理化和二核苷酸等描述符，从而得到固定顺序的 240 项特征结果。

`all_descriptors()` 使用固定版本 `descriptor_schema_v1`，并同时返回不可计算原因、计算条件和来源信息。字段顺序不会随输入改变；前 180 项按输入及各自适用域计算，DNAKit 不内置二核苷酸数值表，因此后 60 项默认均为 `None`。

<span id="1"></span>**最短用法**

```python
from dnakit import DNASequence
from dnakit.descriptors import all_descriptors

result = all_descriptors(DNASequence("ACGT"))
print(len(result.values))
print(result.values["epsilon260_ss_m_inverse_cm_inverse"])
print(sum(name.startswith("diprodb_") and value is None for name, value in result.values.items()))
print(result.unavailable_reasons["diprodb_twist_mean"])
```

```text
240
40300.0
60
requires an explicit user-supplied DinucleotidePropertyTable; DNAKit bundles no DiProDB numerical values
```

CLI 默认输出同一套完整字段；如需旧版四组精简结果，可显式使用 `--compact`：

```bash
dnakit describe ACGT
dnakit describe ACGT --compact
```

<span id="2-api"></span>**API**

- `dnakit.descriptors.all_descriptors(value[必须], ambiguity_policy[可选], conditions[可选], dinucleotide_property_table[可选])`：返回全部字段；只有显式提供表时才计算后 60 项。
- `dnakit.descriptors.load_dinucleotide_property_table(path[必须])`：有界读取并严格验证用户 JSON 表，记录文件 SHA-256。
- `dnakit.descriptors.descriptor_schema_v1()`：返回 240 个不可变 `DescriptorField`，每项包含 `index`、`name`、`category`、`unit`、`formula` 和 `source`。
- `dnakit.descriptors.DESCRIPTOR_SCHEMA_V1`：同一固定 schema 的常量形式。
- `dnakit.descriptors.DESCRIPTOR_NAMES_V1`：仅包含有序字段名。

结果对象结构：

| 属性                    | 含义                                                              |
| ----------------------- | ----------------------------------------------------------------- |
| `schema_version`      | 固定为`descriptor_schema_v1`                                    |
| `sequence_id`         | 输入为`DNARecord` 时保留记录 ID，否则为 `None`                |
| `values`              | 严格按 schema 排列的 240 个值                                     |
| `unavailable_reasons` | 仅包含值为`None` 的字段，且每个字段恰有一个原因                 |
| `conditions`          | IUPAC、Gap、k-mer、ORF、重复、热力学和参数表条件                  |
| `provenance`          | DNAKit 版本、实现类型，以及用户表声明的名称、版本、来源和 SHA-256 |

<span id="3"></span>**统一计算口径**

- 仅 A/C/G/T 计入 canonical 分母；默认 `ambiguity_policy="ignore"`，也可设为 `"error"`。
- k-mer、CpG/GpC、用户二核苷酸表和其他相邻词统计均允许重叠，但不会跨越 IUPAC 模糊符号或显式 Gap。
- 所有比例在分母为 0 时返回 `None`，不会用伪造的 `0` 代替；原因写入 `unavailable_reasons`。
- 用户表统计的 `sd` 使用总体标准差；表必须提供规定的 15 组属性和每组 16 个 DNA 二核苷酸有限数值。
- frame 0 使用 NCBI standard genetic code 1；完整 ORF 使用 `ATG` 起始和 `TAA/TAG/TGA` 终止，并扫描正反链共六个阅读框。
- 热力学默认条件为 37 °C、Na⁺ 0.05 M、链浓度 250 nM；Wallace 仅限 2–13 nt，SantaLucia NN 仅限 2–60 nt。
- 分子量和 ε260 只用于线性、无 Gap、无修饰、A/C/G/T DNA；ε260 是理论系数，不是实验 A260。
- 当前 LZ76 对线性、无 Gap、无模糊符号且不超过 10,000 nt 的序列计算；更长序列返回 `None`，避免二次时间复杂度失控。

<span id="4"></span>**类别与字段范围**

| 序号     | 类别                                     | 字段数 |
| -------- | ---------------------------------------- | -----: |
| 1–12    | 基础长度与数据质量                       |     12 |
| 13–28   | 碱基化学分组组成                         |     16 |
| 29–112  | 1/2/3-mer 频率                           |     84 |
| 113–128 | skew、CpG/GpC 与 Chargaff 指标           |     16 |
| 129–148 | 熵、复杂度、homopolymer 与重复           |     20 |
| 149–164 | 密码子与六阅读框 ORF                     |     16 |
| 165–180 | 分子量、ε260、Tm 与 NN 热力学           |     16 |
| 181–240 | 15 组用户二核苷酸参数 × mean/sd/min/max |     60 |

<span id="5-240"></span>**240 项完整字段表**

|   # | 字段                                         | 类别                      | 单位            | 公式                                                                                         | 来源                                                               |
| --: | -------------------------------------------- | ------------------------- | --------------- | -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
|   1 | `symbol_length`                            | `basic`                 | `nt`          | `number of nucleotide symbols; explicit gaps excluded`                                     | DNAKit descriptor_schema_v1                                        |
|   2 | `coordinate_span`                          | `basic`                 | `nt`          | `symbol_length + sum(known gap lengths); undefined with an unknown gap`                    | DNAKit descriptor_schema_v1                                        |
|   3 | `canonical_base_count`                     | `basic`                 | `count`       | `count(A,C,G,T)`                                                                           | DNAKit descriptor_schema_v1                                        |
|   4 | `ambiguity_symbol_count`                   | `basic`                 | `count`       | `symbol_length - canonical_base_count`                                                     | DNAKit descriptor_schema_v1                                        |
|   5 | `gap_object_count`                         | `basic`                 | `count`       | `number of explicit Gap objects`                                                           | DNAKit descriptor_schema_v1                                        |
|   6 | `known_gap_nt`                             | `basic`                 | `nt`          | `sum(length of known explicit gaps)`                                                       | DNAKit descriptor_schema_v1                                        |
|   7 | `unknown_gap_count`                        | `basic`                 | `count`       | `number of explicit gaps with unknown length`                                              | DNAKit descriptor_schema_v1                                        |
|   8 | `canonical_symbol_fraction`                | `basic`                 | `fraction`    | `canonical_base_count / symbol_length`                                                     | DNAKit descriptor_schema_v1                                        |
|   9 | `ambiguity_symbol_fraction`                | `basic`                 | `fraction`    | `ambiguity_symbol_count / symbol_length`                                                   | DNAKit descriptor_schema_v1                                        |
|  10 | `known_gap_fraction`                       | `basic`                 | `fraction`    | `known_gap_nt / coordinate_span`                                                           | DNAKit descriptor_schema_v1                                        |
|  11 | `canonical_run_count`                      | `basic`                 | `count`       | `number of uninterrupted A/C/G/T runs split by ambiguity or gaps`                          | DNAKit descriptor_schema_v1                                        |
|  12 | `longest_canonical_run_nt`                 | `basic`                 | `nt`          | `maximum uninterrupted A/C/G/T run length`                                                 | DNAKit descriptor_schema_v1                                        |
|  13 | `purine_count`                             | `composition`           | `count`       | `count(A)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  14 | `purine_fraction`                          | `composition`           | `fraction`    | `purine_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  15 | `pyrimidine_count`                         | `composition`           | `count`       | `count(C)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  16 | `pyrimidine_fraction`                      | `composition`           | `fraction`    | `pyrimidine_count / canonical_base_count`                                                  | DNAKit descriptor_schema_v1                                        |
|  17 | `amino_count`                              | `composition`           | `count`       | `count(A)+count(C)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  18 | `amino_fraction`                           | `composition`           | `fraction`    | `amino_count / canonical_base_count`                                                       | DNAKit descriptor_schema_v1                                        |
|  19 | `keto_count`                               | `composition`           | `count`       | `count(G)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  20 | `keto_fraction`                            | `composition`           | `fraction`    | `keto_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  21 | `weak_count`                               | `composition`           | `count`       | `count(A)+count(T)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  22 | `weak_fraction`                            | `composition`           | `fraction`    | `weak_count / canonical_base_count`                                                        | DNAKit descriptor_schema_v1                                        |
|  23 | `strong_count`                             | `composition`           | `count`       | `count(C)+count(G)`                                                                        | DNAKit descriptor_schema_v1                                        |
|  24 | `strong_fraction`                          | `composition`           | `fraction`    | `strong_count / canonical_base_count`                                                      | DNAKit descriptor_schema_v1                                        |
|  25 | `purine_pyrimidine_skew`                   | `composition`           | `ratio`       | `(purine_count-pyrimidine_count)/(purine_count+pyrimidine_count)`                          | DNAKit descriptor_schema_v1                                        |
|  26 | `amino_keto_skew`                          | `composition`           | `ratio`       | `(amino_count-keto_count)/(amino_count+keto_count)`                                        | DNAKit descriptor_schema_v1                                        |
|  27 | `weak_strong_skew`                         | `composition`           | `ratio`       | `(weak_count-strong_count)/(weak_count+strong_count)`                                      | DNAKit descriptor_schema_v1                                        |
|  28 | `gc_at_ratio`                              | `composition`           | `ratio`       | `strong_count / weak_count`                                                                | DNAKit descriptor_schema_v1                                        |
|  29 | `k1_A_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(A) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  30 | `k1_C_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(C) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  31 | `k1_G_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(G) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  32 | `k1_T_frequency`                           | `kmer`                  | `fraction`    | `overlapping count(T) / valid canonical 1-mer positions`                                   | DNAKit exact overlapping k-mer definition                          |
|  33 | `k2_AA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  34 | `k2_AC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  35 | `k2_AG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  36 | `k2_AT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(AT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  37 | `k2_CA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  38 | `k2_CC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  39 | `k2_CG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  40 | `k2_CT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(CT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  41 | `k2_GA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  42 | `k2_GC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  43 | `k2_GG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  44 | `k2_GT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(GT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  45 | `k2_TA_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TA) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  46 | `k2_TC_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TC) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  47 | `k2_TG_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TG) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  48 | `k2_TT_frequency`                          | `kmer`                  | `fraction`    | `overlapping count(TT) / valid canonical 2-mer positions`                                  | DNAKit exact overlapping k-mer definition                          |
|  49 | `k3_AAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  50 | `k3_AAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  51 | `k3_AAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  52 | `k3_AAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  53 | `k3_ACA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  54 | `k3_ACC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  55 | `k3_ACG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  56 | `k3_ACT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ACT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  57 | `k3_AGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  58 | `k3_AGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  59 | `k3_AGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  60 | `k3_AGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(AGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  61 | `k3_ATA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  62 | `k3_ATC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  63 | `k3_ATG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  64 | `k3_ATT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(ATT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  65 | `k3_CAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  66 | `k3_CAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  67 | `k3_CAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  68 | `k3_CAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  69 | `k3_CCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  70 | `k3_CCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  71 | `k3_CCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  72 | `k3_CCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  73 | `k3_CGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  74 | `k3_CGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  75 | `k3_CGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  76 | `k3_CGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  77 | `k3_CTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  78 | `k3_CTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  79 | `k3_CTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  80 | `k3_CTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(CTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  81 | `k3_GAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  82 | `k3_GAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  83 | `k3_GAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  84 | `k3_GAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  85 | `k3_GCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  86 | `k3_GCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  87 | `k3_GCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  88 | `k3_GCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  89 | `k3_GGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  90 | `k3_GGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  91 | `k3_GGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  92 | `k3_GGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  93 | `k3_GTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  94 | `k3_GTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  95 | `k3_GTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  96 | `k3_GTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(GTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  97 | `k3_TAA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  98 | `k3_TAC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
|  99 | `k3_TAG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 100 | `k3_TAT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TAT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 101 | `k3_TCA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 102 | `k3_TCC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 103 | `k3_TCG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 104 | `k3_TCT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TCT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 105 | `k3_TGA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 106 | `k3_TGC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 107 | `k3_TGG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 108 | `k3_TGT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TGT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 109 | `k3_TTA_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTA) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 110 | `k3_TTC_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTC) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 111 | `k3_TTG_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTG) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 112 | `k3_TTT_frequency`                         | `kmer`                  | `fraction`    | `overlapping count(TTT) / valid canonical 3-mer positions`                                 | DNAKit exact overlapping k-mer definition                          |
| 113 | `gc_skew`                                  | `skew_cpg`              | `ratio`       | `(count(G)-count(C))/(count(G)+count(C))`                                                  | DNAKit descriptor_schema_v1                                        |
| 114 | `at_skew`                                  | `skew_cpg`              | `ratio`       | `(count(A)-count(T))/(count(A)+count(T))`                                                  | DNAKit descriptor_schema_v1                                        |
| 115 | `cpg_count`                                | `skew_cpg`              | `count`       | `overlapping count(CG)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 116 | `cpg_density`                              | `skew_cpg`              | `fraction`    | `count(CG) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 117 | `cpg_observed_expected`                    | `skew_cpg`              | `ratio`       | `count(CG)*canonical_base_count/(count(C)*count(G))`                                       | DNAKit descriptor_schema_v1                                        |
| 118 | `gpc_count`                                | `skew_cpg`              | `count`       | `overlapping count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 119 | `gpc_density`                              | `skew_cpg`              | `fraction`    | `count(GC) / valid canonical dinucleotide positions`                                       | DNAKit descriptor_schema_v1                                        |
| 120 | `cpg_gpc_ratio`                            | `skew_cpg`              | `ratio`       | `count(CG) / count(GC)`                                                                    | DNAKit descriptor_schema_v1                                        |
| 121 | `cumulative_gc_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 122 | `cumulative_gc_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where G=+1,C=-1,A/T=0`                                        | DNAKit descriptor_schema_v1                                        |
| 123 | `cumulative_gc_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_gc_skew_max - cumulative_gc_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 124 | `cumulative_at_skew_max`                   | `skew_cpg`              | `count`       | `max prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 125 | `cumulative_at_skew_min`                   | `skew_cpg`              | `count`       | `min prefix cumulative score where A=+1,T=-1,C/G=0`                                        | DNAKit descriptor_schema_v1                                        |
| 126 | `cumulative_at_skew_range`                 | `skew_cpg`              | `count`       | `cumulative_at_skew_max - cumulative_at_skew_min`                                          | DNAKit descriptor_schema_v1                                        |
| 127 | `dinucleotide_rc_total_variation`          | `skew_cpg`              | `fraction`    | `0.5*sum_xy(abs(f_xy-f_reverse_complement(xy)))`                                           | DNAKit descriptor_schema_v1                                        |
| 128 | `mono_chargaff_l1_distance`                | `skew_cpg`              | `fraction`    | `abs(f_A-f_T)+abs(f_C-f_G)`                                                                | DNAKit descriptor_schema_v1                                        |
| 129 | `shannon_entropy_k1_bits`                  | `complexity`            | `bits`        | `-sum(p(1-mer)*log2(p(1-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 130 | `shannon_entropy_k2_bits`                  | `complexity`            | `bits`        | `-sum(p(2-mer)*log2(p(2-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 131 | `shannon_entropy_k3_bits`                  | `complexity`            | `bits`        | `-sum(p(3-mer)*log2(p(3-mer)))`                                                            | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 132 | `normalized_entropy_k1`                    | `complexity`            | `fraction`    | `shannon_entropy_k1_bits / log2(4**1)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 133 | `normalized_entropy_k2`                    | `complexity`            | `fraction`    | `shannon_entropy_k2_bits / log2(4**2)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 134 | `normalized_entropy_k3`                    | `complexity`            | `fraction`    | `shannon_entropy_k3_bits / log2(4**3)`                                                     | Shannon 1948; DOI 10.1002/j.1538-7305.1948.tb01338.x               |
| 135 | `linguistic_complexity_k2`                 | `complexity`            | `fraction`    | `unique 2-mers / min(4**2, valid 2-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 136 | `linguistic_complexity_k3`                 | `complexity`            | `fraction`    | `unique 3-mers / min(4**3, valid 3-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 137 | `linguistic_complexity_k4`                 | `complexity`            | `fraction`    | `unique 4-mers / min(4**4, valid 4-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 138 | `linguistic_complexity_k5`                 | `complexity`            | `fraction`    | `unique 5-mers / min(4**5, valid 5-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 139 | `linguistic_complexity_k6`                 | `complexity`            | `fraction`    | `unique 6-mers / min(4**6, valid 6-mer positions)`                                         | Observed/possible k-word linguistic complexity                     |
| 140 | `linguistic_complexity_product_k1_k6`      | `complexity`            | `fraction`    | `product of defined linguistic_complexity_k values for k=1..6`                             | Observed/possible k-word linguistic complexity                     |
| 141 | `lz76_complexity`                          | `complexity`            | `count`       | `number of phrases in exhaustive Lempel-Ziv 1976 parsing`                                  | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 142 | `normalized_lz76_complexity`               | `complexity`            | `ratio`       | `lz76_complexity*log_base4(canonical_base_count)/canonical_base_count`                     | Lempel and Ziv 1976; DOI 10.1109/TIT.1976.1055501                  |
| 143 | `longest_homopolymer_nt`                   | `complexity`            | `nt`          | `max canonical homopolymer run`                                                            | DNAKit descriptor_schema_v1                                        |
| 144 | `longest_homopolymer_a_nt`                 | `complexity`            | `nt`          | `max homopolymer run of A`                                                                 | DNAKit descriptor_schema_v1                                        |
| 145 | `longest_homopolymer_c_nt`                 | `complexity`            | `nt`          | `max homopolymer run of C`                                                                 | DNAKit descriptor_schema_v1                                        |
| 146 | `longest_homopolymer_g_nt`                 | `complexity`            | `nt`          | `max homopolymer run of G`                                                                 | DNAKit descriptor_schema_v1                                        |
| 147 | `longest_homopolymer_t_nt`                 | `complexity`            | `nt`          | `max homopolymer run of T`                                                                 | DNAKit descriptor_schema_v1                                        |
| 148 | `exact_tandem_repeat_coverage_fraction`    | `complexity`            | `fraction`    | `union bases covered by exact tandem repeats / canonical_base_count`                       | DNAKit exact tandem repeat scanner; units 1..20; minimum repeats 2 |
| 149 | `frame0_codon_count`                       | `coding`                | `count`       | `number of valid forward frame-0 codons`                                                   | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 150 | `frame0_unique_codon_count`                | `coding`                | `count`       | `number of distinct forward frame-0 codons`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 151 | `frame0_start_codon_count`                 | `coding`                | `count`       | `count(ATG) in forward frame 0`                                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 152 | `frame0_stop_codon_count`                  | `coding`                | `count`       | `count(TAA,TAG,TGA) in forward frame 0`                                                    | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 153 | `frame0_start_codon_fraction`              | `coding`                | `fraction`    | `frame0_start_codon_count / frame0_codon_count`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 154 | `frame0_stop_codon_fraction`               | `coding`                | `fraction`    | `frame0_stop_codon_count / frame0_codon_count`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 155 | `frame0_codon_entropy_bits`                | `coding`                | `bits`        | `-sum(frame0 codon frequency*log2(frequency))`                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 156 | `frame0_effective_number_of_codons`        | `coding`                | `count`       | `2**frame0_codon_entropy_bits`                                                             | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 157 | `frame0_gc1_fraction`                      | `coding`                | `fraction`    | `GC bases at position 1 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 158 | `frame0_gc2_fraction`                      | `coding`                | `fraction`    | `GC bases at position 2 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 159 | `frame0_gc3_fraction`                      | `coding`                | `fraction`    | `GC bases at position 3 / valid frame-0 codons`                                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 160 | `six_frame_complete_orf_count`             | `coding`                | `count`       | `complete start-to-next-stop ORFs across three frames on both strands`                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 161 | `six_frame_forward_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three forward frames`                                                | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 162 | `six_frame_reverse_complete_orf_count`     | `coding`                | `count`       | `complete ORFs across three reverse-complement frames`                                     | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 163 | `six_frame_longest_complete_orf_nt`        | `coding`                | `nt`          | `maximum complete six-frame ORF length including terminal stop`                            | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 164 | `six_frame_complete_orf_coverage_fraction` | `coding`                | `fraction`    | `union of symbol positions covered by complete six-frame ORFs / canonical_base_count`      | NCBI standard genetic code table 1; DNAKit ORF rules               |
| 165 | `mw_ss_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime OH`                                       | DNAKit anhydrous DNA residue mass table v1                         |
| 166 | `mw_ss_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of one ssDNA strand with 5-prime phosphate`                                | DNAKit anhydrous DNA residue mass table v1                         |
| 167 | `mw_ds_oh_da`                              | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement with 5-prime OH`              | DNAKit anhydrous DNA residue mass table v1                         |
| 168 | `mw_ds_5p_phosphate_da`                    | `physicochemical`       | `Da`          | `anhydrous mass of sequence plus complete reverse complement; both 5-prime phosphorylated` | DNAKit anhydrous DNA residue mass table v1                         |
| 169 | `epsilon260_ss_m_inverse_cm_inverse`       | `physicochemical`       | `M^-1 cm^-1`  | `nearest-neighbor epsilon260 pair sum minus internal-base sum`                             | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 170 | `nmol_per_a260_1ml_1cm`                    | `physicochemical`       | `nmol`        | `1e6 / epsilon260 for A260=1, volume=1 mL, path=1 cm`                                      | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 171 | `ug_per_a260_1ml_1cm`                      | `physicochemical`       | `ug`          | `1000*mw_ss_oh_da/epsilon260 for A260=1, volume=1 mL, path=1 cm`                           | Warshaw-Tinoco 1966 and Cantor-Warshaw-Shapiro 1970                |
| 172 | `tm_wallace_c`                             | `physicochemical`       | `degree C`    | `2*(A+T)+4*(G+C)`                                                                          | Wallace short-oligo 2AT+4GC rule                                   |
| 173 | `stacking_delta_h_kcal_per_mol`            | `physicochemical`       | `kcal/mol`    | `sum SantaLucia nearest-neighbor stacking delta H`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 174 | `stacking_delta_s_cal_per_k_mol_k`         | `physicochemical`       | `cal/(K mol)` | `sum SantaLucia nearest-neighbor stacking delta S`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 175 | `stacking_delta_g37_kcal_per_mol`          | `physicochemical`       | `kcal/mol`    | `stacking delta H - 310.15*stacking delta S/1000`                                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 176 | `nn_delta_h_kcal_per_mol`                  | `physicochemical`       | `kcal/mol`    | `SantaLucia complete-duplex delta H with initiation and symmetry`                          | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 177 | `nn_delta_s_cal_per_mol_k`                 | `physicochemical`       | `cal/(K mol)` | `SantaLucia complete-duplex delta S with initiation, symmetry, and salt`                   | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 178 | `nn_delta_g37_kcal_per_mol`                | `physicochemical`       | `kcal/mol`    | `complete-duplex delta H - 310.15*delta S/1000`                                            | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 179 | `nn_tm_c`                                  | `physicochemical`       | `degree C`    | `SantaLucia concentration- and sodium-adjusted Tm`                                         | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 180 | `self_complementary`                       | `physicochemical`       | `boolean`     | `sequence == reverse_complement(sequence)`                                                 | SantaLucia 1998; DOI 10.1073/pnas.95.4.1460                        |
| 181 | `diprodb_twist_mean`                       | `dinucleotide_property` | `degree`      | `population mean of Twist values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 182 | `diprodb_twist_sd`                         | `dinucleotide_property` | `degree`      | `population sd of Twist values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 183 | `diprodb_twist_min`                        | `dinucleotide_property` | `degree`      | `population min of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 184 | `diprodb_twist_max`                        | `dinucleotide_property` | `degree`      | `population max of Twist values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 185 | `diprodb_tilt_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Tilt values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 186 | `diprodb_tilt_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Tilt values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 187 | `diprodb_tilt_min`                         | `dinucleotide_property` | `degree`      | `population min of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 188 | `diprodb_tilt_max`                         | `dinucleotide_property` | `degree`      | `population max of Tilt values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 189 | `diprodb_roll_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Roll values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 190 | `diprodb_roll_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Roll values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 191 | `diprodb_roll_min`                         | `dinucleotide_property` | `degree`      | `population min of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 192 | `diprodb_roll_max`                         | `dinucleotide_property` | `degree`      | `population max of Roll values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 193 | `diprodb_shift_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Shift values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 194 | `diprodb_shift_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Shift values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 195 | `diprodb_shift_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 196 | `diprodb_shift_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Shift values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 197 | `diprodb_slide_mean`                       | `dinucleotide_property` | `angstrom`    | `population mean of Slide values over valid overlapping dinucleotides`                     | Caller-supplied table; DNAKit bundles no numerical values          |
| 198 | `diprodb_slide_sd`                         | `dinucleotide_property` | `angstrom`    | `population sd of Slide values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 199 | `diprodb_slide_min`                        | `dinucleotide_property` | `angstrom`    | `population min of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 200 | `diprodb_slide_max`                        | `dinucleotide_property` | `angstrom`    | `population max of Slide values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 201 | `diprodb_rise_mean`                        | `dinucleotide_property` | `angstrom`    | `population mean of Rise values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 202 | `diprodb_rise_sd`                          | `dinucleotide_property` | `angstrom`    | `population sd of Rise values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 203 | `diprodb_rise_min`                         | `dinucleotide_property` | `angstrom`    | `population min of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 204 | `diprodb_rise_max`                         | `dinucleotide_property` | `angstrom`    | `population max of Rise values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 205 | `diprodb_bend_mean`                        | `dinucleotide_property` | `degree`      | `population mean of Bend values over valid overlapping dinucleotides`                      | Caller-supplied table; DNAKit bundles no numerical values          |
| 206 | `diprodb_bend_sd`                          | `dinucleotide_property` | `degree`      | `population sd of Bend values over valid overlapping dinucleotides`                        | Caller-supplied table; DNAKit bundles no numerical values          |
| 207 | `diprodb_bend_min`                         | `dinucleotide_property` | `degree`      | `population min of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 208 | `diprodb_bend_max`                         | `dinucleotide_property` | `degree`      | `population max of Bend values over valid overlapping dinucleotides`                       | Caller-supplied table; DNAKit bundles no numerical values          |
| 209 | `diprodb_inclination_mean`                 | `dinucleotide_property` | `degree`      | `population mean of Inclination values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 210 | `diprodb_inclination_sd`                   | `dinucleotide_property` | `degree`      | `population sd of Inclination values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 211 | `diprodb_inclination_min`                  | `dinucleotide_property` | `degree`      | `population min of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 212 | `diprodb_inclination_max`                  | `dinucleotide_property` | `degree`      | `population max of Inclination values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 213 | `diprodb_direction_mean`                   | `dinucleotide_property` | `degree`      | `population mean of Direction values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 214 | `diprodb_direction_sd`                     | `dinucleotide_property` | `degree`      | `population sd of Direction values over valid overlapping dinucleotides`                   | Caller-supplied table; DNAKit bundles no numerical values          |
| 215 | `diprodb_direction_min`                    | `dinucleotide_property` | `degree`      | `population min of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 216 | `diprodb_direction_max`                    | `dinucleotide_property` | `degree`      | `population max of Direction values over valid overlapping dinucleotides`                  | Caller-supplied table; DNAKit bundles no numerical values          |
| 217 | `diprodb_propeller_twist_mean`             | `dinucleotide_property` | `degree`      | `population mean of Propeller twist values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 218 | `diprodb_propeller_twist_sd`               | `dinucleotide_property` | `degree`      | `population sd of Propeller twist values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 219 | `diprodb_propeller_twist_min`              | `dinucleotide_property` | `degree`      | `population min of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 220 | `diprodb_propeller_twist_max`              | `dinucleotide_property` | `degree`      | `population max of Propeller twist values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 221 | `diprodb_major_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Major groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 222 | `diprodb_major_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Major groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 223 | `diprodb_major_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 224 | `diprodb_major_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Major groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 225 | `diprodb_minor_groove_width_mean`          | `dinucleotide_property` | `angstrom`    | `population mean of Minor groove width values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 226 | `diprodb_minor_groove_width_sd`            | `dinucleotide_property` | `angstrom`    | `population sd of Minor groove width values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 227 | `diprodb_minor_groove_width_min`           | `dinucleotide_property` | `angstrom`    | `population min of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 228 | `diprodb_minor_groove_width_max`           | `dinucleotide_property` | `angstrom`    | `population max of Minor groove width values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 229 | `diprodb_persistence_length_mean`          | `dinucleotide_property` | `nanometer`   | `population mean of Persistence length values over valid overlapping dinucleotides`        | Caller-supplied table; DNAKit bundles no numerical values          |
| 230 | `diprodb_persistence_length_sd`            | `dinucleotide_property` | `nanometer`   | `population sd of Persistence length values over valid overlapping dinucleotides`          | Caller-supplied table; DNAKit bundles no numerical values          |
| 231 | `diprodb_persistence_length_min`           | `dinucleotide_property` | `nanometer`   | `population min of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 232 | `diprodb_persistence_length_max`           | `dinucleotide_property` | `nanometer`   | `population max of Persistence length values over valid overlapping dinucleotides`         | Caller-supplied table; DNAKit bundles no numerical values          |
| 233 | `diprodb_stacking_energy_mean`             | `dinucleotide_property` | `kcal/mol`    | `population mean of Stacking energy values over valid overlapping dinucleotides`           | Caller-supplied table; DNAKit bundles no numerical values          |
| 234 | `diprodb_stacking_energy_sd`               | `dinucleotide_property` | `kcal/mol`    | `population sd of Stacking energy values over valid overlapping dinucleotides`             | Caller-supplied table; DNAKit bundles no numerical values          |
| 235 | `diprodb_stacking_energy_min`              | `dinucleotide_property` | `kcal/mol`    | `population min of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 236 | `diprodb_stacking_energy_max`              | `dinucleotide_property` | `kcal/mol`    | `population max of Stacking energy values over valid overlapping dinucleotides`            | Caller-supplied table; DNAKit bundles no numerical values          |
| 237 | `diprodb_free_energy_mean`                 | `dinucleotide_property` | `kcal/mol`    | `population mean of Free energy values over valid overlapping dinucleotides`               | Caller-supplied table; DNAKit bundles no numerical values          |
| 238 | `diprodb_free_energy_sd`                   | `dinucleotide_property` | `kcal/mol`    | `population sd of Free energy values over valid overlapping dinucleotides`                 | Caller-supplied table; DNAKit bundles no numerical values          |
| 239 | `diprodb_free_energy_min`                  | `dinucleotide_property` | `kcal/mol`    | `population min of Free energy values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |
| 240 | `diprodb_free_energy_max`                  | `dinucleotide_property` | `kcal/mol`    | `population max of Free energy values over valid overlapping dinucleotides`                | Caller-supplied table; DNAKit bundles no numerical values          |

方法、论文、数据库与网址已统一移至[致谢与主要来源](../../acknowledgements.md#methods-and-references)；许可和用户责任见[第三方声明](../../acknowledgements.md#third-party-notices)。



<span id="all-descriptors"></span>
