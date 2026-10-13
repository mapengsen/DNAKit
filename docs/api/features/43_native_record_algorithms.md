# 原生 FASTA 与变异记录算法

这些函数只处理调用者明确提供的 FASTA/FASTQ 或 VCF-like 记录，不执行比对、变异检测、参考基因组校正或注释查询。坐标和基因型规则写入方法与 provenance，适合用于可复现的论文预处理或描述性统计。

方法引用可分别使用 [Biopython](https://doi.org/10.1093/bioinformatics/btp163) 的序列处理论文和 [SAMtools/BCFtools/HTSlib](https://doi.org/10.1093/gigascience/giab008) 的变异记录工具论文；这些 DOI 证明相关软件方法有期刊来源，不替代本项目对版本、参数和数据的独立验证。

## 1) 计算 FASTA 记录 GC profile

- **作用：** 为每条 FASTA 记录计算长度、GC 分数和 GC 百分比。
- **处理方法：** 解析文本或内存记录，统一 DNA 字母大小写后逐条计数；空序列的 GC 分数为 0。
- **API：** `dnakit.fasta_fastq.calculate_gc_profile(source[必须])`。
- **输入：** FASTA 文本、路径或 `(id, sequence)` 记录序列。
- **输出：** 返回包含 `id`、`length`、`gc_fraction` 和 `gc_percent` 的列表。
- **示例代码：**

```python
from dnakit.fasta_fastq import calculate_gc_profile

print(calculate_gc_profile(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **示例结果：**

```json
[
  {"id": "a", "length": 4, "gc_fraction": 0.5, "gc_percent": 50.0},
  {"id": "b", "length": 4, "gc_fraction": 0.5, "gc_percent": 50.0}
]
```

## 2) 统计 FASTA/FASTQ 长度分布

- **作用：** 汇总记录数量、最小长度、最大长度、平均长度和长度计数。
- **处理方法：** 只使用每条记录的序列长度，不推断质量或测序平台。
- **API：** `dnakit.fasta_fastq.calculate_length_distribution(source[必须])`。
- **输入：** FASTA/FASTQ 文本、路径或已物化记录。
- **输出：** 返回 `count`、`min`、`max`、`mean` 和 `counts`。
- **示例代码：**

```python
from dnakit.fasta_fastq import calculate_length_distribution

print(calculate_length_distribution(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **示例结果：**

```json
{"count": 2, "min": 4, "max": 4, "mean": 4.0, "counts": {"4": 2}}
```

## 3) 去除重复 FASTA 序列

- **作用：** 按规范化后的序列内容去重并保留首次出现的记录 ID。
- **处理方法：** 逐条解析并统一序列，再以序列字符串作为键；不会把反向互补序列当成重复项。
- **API：** `dnakit.fasta_fastq.remove_duplicate_sequences(source[必须])`。
- **输入：** FASTA 文本、路径或 `(id, sequence)` 记录序列。
- **输出：** 返回首次出现的 `(id, sequence)` 列表。
- **示例代码：**

```python
from dnakit.fasta_fastq import remove_duplicate_sequences

print(remove_duplicate_sequences(">a\nACGT\n>dup\nACGT\n>b\nGGAA\n", verbose=False))
```

- **示例结果：**

```json
[["a", "ACGT"], ["b", "GGAA"]]
```

## 4) 生成记录的反向互补

- **作用：** 为每条 FASTA 记录生成反向互补序列并保留原 ID。
- **处理方法：** 使用 DNA/IUPAC 互补表反转每条序列；该操作不改变记录顺序。
- **API：** `dnakit.fasta_fastq.reverse_complement_records(source[必须])`。
- **输入：** FASTA 文本、路径或 `(id, sequence)` 记录序列。
- **输出：** 返回反向互补后的 `(id, sequence)` 列表。
- **示例代码：**

```python
from dnakit.fasta_fastq import reverse_complement_records

print(reverse_complement_records(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **示例结果：**

```json
[["a", "ACGT"], ["b", "TTCC"]]
```

## 5) 校验 FASTA 序列 ID

- **作用：** 检查 FASTA 记录 ID 是否非空且唯一。
- **处理方法：** 解析记录头部并按输入顺序收集 ID；重复项会列在 `duplicates` 中。
- **API：** `dnakit.fasta_fastq.validate_sequence_ids(source[必须])`。
- **输入：** FASTA 文本、路径或 `(id, sequence)` 记录序列。
- **输出：** 返回 `valid`、`ids` 和 `duplicates`。
- **示例代码：**

```python
from dnakit.fasta_fastq import validate_sequence_ids

print(validate_sequence_ids(">a\nACGT\n>b\nGGAA\n", verbose=False))
```

- **示例结果：**

```json
{"valid": true, "ids": ["a", "b"], "duplicates": []}
```

## 6) 计算 VCF-like 等位基因频率

- **作用：** 根据显式基因型计算每条变异的 ALT 计数、called allele 数和频率。
- **处理方法：** 解析 `0/1`、`1|1` 等基因型，忽略缺失等位基因；不进行基因型调用、填补或参考序列校正。
- **API：** `dnakit.variants.calculate_allele_frequency(input, *, options=None, verbose=True)`。
- **输入：** `{"variants": [{"chrom", "pos", "ref", "alt", "genotypes"}]}`，`pos` 使用 VCF-like 1-based 坐标。
- **输出：** 完成信封中的 `result.variants`、总 ALT/called allele 数和聚合频率。
- **示例代码：**

```python
from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency(
    {"variants": [{
        "chrom": "chr1", "pos": 1, "ref": "A", "alt": "G",
        "genotypes": ["0/1", "1/1", "0/0"],
    }]},
    verbose=False,
)
print(result["result"])
```

- **示例结果：**

```json
{
  "variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": ["G"], "alt_count": 3, "called_alleles": 6, "frequency": 0.5}],
  "total_alt_alleles": 3,
  "total_called_alleles": 6,
  "aggregate_frequency": 0.5
}
```

## 7) 统计 transition/transversion

- **作用：** 按显式 REF/ALT 碱基把 SNP 分类为 transition 或 transversion，并计算 Ti/Tv 比。
- **处理方法：** 只接受单碱基 REF/ALT；非 SNP 记录计入 `excluded`，不进行左对齐或规范化。
- **API：** `dnakit.variants.calculate_transition_transversion(input, *, options=None, verbose=True)`。
- **输入：** `input.variants` 中的记录包含 `chrom`、`pos`、`ref` 和 `alt`。
- **输出：** 完成信封中的分类记录、`transitions`、`transversions`、`excluded` 和 `ti_tv_ratio`。
- **示例代码：**

```python
from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion(
    {"variants": [
        {"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G"},
        {"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C"},
    ]},
    verbose=False,
)
print(result["result"])
```

- **示例结果：**

```json
{
  "variants": [
    {"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "classification": "transition"},
    {"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C", "classification": "transversion"}
  ],
  "transitions": 1,
  "transversions": 1,
  "excluded": 0,
  "ti_tv_ratio": 1.0
}
```

## 8) 按测序深度筛选变异

- **作用：** 按显式 `DP`、`depth` 或 `read_depth` 保留 VCF-like 记录。
- **处理方法：** 对每条记录计算深度并应用 `min_depth`、`max_depth` 和 `keep_missing`。
- **API：** `dnakit.variants.filter_variants_by_depth(input, *, options=None)`。
- **输入：** `input.variants` 中包含 `chrom`、`pos`、`ref`、`alt` 和深度字段的记录。
- **输出：** 完成信封中的保留记录、逐记录指标、输入/输出计数和筛选条件。
- **示例代码：** `dnakit.variants.filter_variants_by_depth({"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "DP": 30}]}, options={"min_depth": 20})`。
- **示例结果：** `{"output_record_count": 1, "removed_record_count": 0}`。

## 9) 按基因型质量筛选变异

- **作用：** 按显式 GQ 值筛选变异记录。
- **处理方法：** 从记录或 sample 映射读取 GQ，使用该记录的最低观测 GQ 与阈值比较。
- **API：** `dnakit.variants.filter_variants_by_genotype_quality(input, *, options=None)`。
- **输入：** 包含 `GQ`、`genotype_quality`、`genotype_qualities` 或 sample GQ 的记录。
- **输出：** 完成信封中的保留记录、GQ 指标、输入/输出计数和筛选条件。
- **示例代码：** `dnakit.variants.filter_variants_by_genotype_quality({"variants": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "G", "GQ": 45}]}, options={"min_genotype_quality": 20})`。
- **示例结果：** `{"output_record_count": 1, "removed_record_count": 0}`。
