# 无外部依赖的序列基础算法

这些函数只处理内存中的 DNA/RNA 字符串，不解析 FASTA/FASTQ，不执行比对，也不推断测序质量。它们适合在论文方法中记录为确定性的序列预处理或描述性统计步骤。

## 1) 计算序列指标

- **作用：** 计算序列长度、GC 比例、N 碱基数、IUPAC 模糊碱基数、熵和碱基组成。
- **处理方法：** 去除空白并统一为大写，按 IUPAC DNA 字母表校验；模糊碱基计入长度，但不会被计入 GC。
- **API：** `dnakit.sequence.calculate_sequence_metrics(sequence[必须])`。
- **输入：** 一个 DNA 字符串，可包含大小写字母、空白和 IUPAC 模糊符号。
- **输出：** 返回可序列化映射，包含 `length`、`gc_fraction`、`gc_percent`、`n_count`、`ambiguous_count`、`entropy_bits` 和 `composition`。
- **示例代码：**

```python
from dnakit.sequence import calculate_sequence_metrics

print(calculate_sequence_metrics("ACGTN", verbose=False))
```

- **示例结果：**

```json
{
  "length": 5,
  "gc_fraction": 0.4,
  "gc_percent": 40.0,
  "n_count": 1,
  "ambiguous_count": 1,
  "entropy_bits": 2.321928094887362,
  "composition": {"A": 1, "C": 1, "G": 1, "N": 1, "T": 1}
}
```

## 2) 统计重叠 k-mer

- **作用：** 统计 DNA 序列中所有重叠的 k-mer。
- **处理方法：** 按从左到右的窗口枚举序列，保留每个 k-mer 的出现次数；不做规范化、纠错或去重。
- **API：** `dnakit.sequence.count_kmers(sequence[必须], k[可选])`。
- **输入：** DNA 字符串和正整数 `k`。
- **输出：** 返回按 k-mer 字典序排列的计数字典。
- **示例代码：**

```python
from dnakit.sequence import count_kmers

print(count_kmers("ACACA", k=2, verbose=False))
```

- **示例结果：**

```json
{"AC": 2, "CA": 2}
```

## 3) 去除重复序列

- **作用：** 去除重复 DNA 序列并保留首次出现顺序。
- **处理方法：** 先去除空白并统一大小写；`canonical=True` 时使用序列和反向互补序列中字典序较小者作为去重键。
- **API：** `dnakit.sequence.deduplicate_sequences(sequences[必须], canonical[可选])`。
- **输入：** DNA 字符串序列；可选 `canonical` 布尔值。
- **输出：** 返回规范化后的唯一序列列表。
- **示例代码：**

```python
from dnakit.sequence import deduplicate_sequences

print(deduplicate_sequences(["ACG", " ACG ", "CGT"], canonical=True, verbose=False))
```

- **示例结果：**

```json
["ACG"]
```

## 4) 计算 DNA 反向互补

- **作用：** 返回 DNA 的反向互补序列，并支持 IUPAC 模糊符号。
- **处理方法：** 先按 IUPAC DNA 字母表校验，再反转字符顺序并应用互补表；该函数只处理序列字符串。
- **API：** `dnakit.sequence.reverse_complement_dna(sequence[必须])`。
- **输入：** DNA 字符串，可包含 IUPAC 模糊符号。
- **输出：** 返回大写 DNA 字符串。
- **示例代码：**

```python
from dnakit.sequence import reverse_complement_dna

print(reverse_complement_dna("ACGTRN", verbose=False))
```

- **示例结果：**

```text
NYACGT
```

## 5) 汇总序列长度和 GC 比例

- **作用：** 对一组 DNA/RNA 序列一次遍历汇总数量、长度分布和平均 GC 比例。
- **处理方法：** 忽略空白、统一大小写并把 U 转为 T；IUPAC 模糊符号计入长度但不计入 GC；输入迭代器只消费一次。
- **API：** `dnakit.sequence.sequence_length_gc_summary(sequences[必须])`。
- **输入：** 可迭代的 DNA/RNA 字符串。
- **输出：** 返回可序列化映射，包含序列数、每条长度、最小/最大/平均/中位长度和平均 GC 比例。
- **示例代码：**

```python
from dnakit.sequence import sequence_length_gc_summary

print(sequence_length_gc_summary(["ACGT", "GGNN"], verbose=False))
```

- **示例结果：**

```json
{
  "count": 2,
  "lengths": [4, 4],
  "min_length": 4,
  "max_length": 4,
  "mean_length": 4,
  "median_length": 4.0,
  "mean_gc_fraction": 0.5
}
```
