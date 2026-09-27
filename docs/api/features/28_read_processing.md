# 测序数据整理

## 1) 匹配双端 reads {#bio-10}

- **作用：** 按名称匹配 R1/R2，并识别缺少配对的记录。
- **计算方法：** 去掉 ID 末尾的 `/1`、`/2` 后匹配，缺失的一端返回 `None`。
- **API：** `dnakit.bioinformatics.reads.pair_reads(first, second, max_pending=100000, include_orphans=True)`
- **输入：** 两个 `Read` 迭代器；每端规范化后的 ID 不能重复。
- **输出：** 返回 `ReadPair` 迭代器，每项包含配对 ID 和两端 reads。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Read, pair_reads

pairs = pair_reads([Read("r1/1", "ACGT")], [Read("r1/2", "TGCA")])
print([(pair.id, pair.first.sequence, pair.second.sequence) for pair in pairs])
```

- **示例结果：**

```text
[('r1', 'ACGT', 'TGCA')]
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 2) 恢复损坏 FASTQ 中的合法记录 {#bio-11}

- **作用：** 从损坏文件中保留可识别的四行 FASTQ 记录。
- **计算方法：** 扫描记录头、序列、分隔符与等长质量行，并记录丢弃行数。
- **API：** `dnakit.bioinformatics.reads.repair_fastq(source, destination, phred_offset=33)`
- **输入：** 损坏 FASTQ 和新输出路径；多行 FASTQ 应用 `read_fastx` 严格读取。
- **输出：** 返回 `RepairResult`，包含修复后 FASTQ 路径、恢复的 reads 数和丢弃的行数。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import repair_fastq

with TemporaryDirectory() as folder:
    source = Path(folder) / "broken.fq"
    source.write_text("invalid\n@r1\nACGT\n+\nIIII\n")
    result = repair_fastq(source, Path(folder) / "repaired.fq")
    print(result.recovered_reads, result.discarded_lines)
```

- **示例结果：**

```text
1 1
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 3) 排序、取交集、重命名和校验 reads {#bio-12}

- **作用：** 整理记录顺序、提取共有序列并计算序列摘要。
- **计算方法：** 排序或固定种子打乱；按 ID/序列取交集；正则重命名；按记录边界计算 SHA256。
- **API：** `dnakit.bioinformatics.reads.arrange_reads(reads, operation='sort')`; `dnakit.bioinformatics.reads.common_reads(collections)`; `dnakit.bioinformatics.reads.rename_reads(reads, pattern, replacement)`; `dnakit.bioinformatics.reads.sequence_digest(reads)`
- **输入：** `Read` 集合；摘要与记录顺序有关，与 ID 无关。
- **输出：** 返回排序、交集或重命名后的 `Read` 记录；校验函数返回比较结果。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Read, arrange_reads, common_reads, rename_reads, sequence_digest

reads = [Read("b", "TGCA"), Read("a", "ACGT")]
ordered = arrange_reads(reads)
shared = common_reads([ordered, [Read("a", "ACGT")]])
renamed = list(rename_reads(shared, "^a$", "sample_a"))
print([read.id for read in ordered])
print([read.id for read in renamed])
print(sequence_digest(shared) == sequence_digest(renamed))
```

- **示例结果：**

```text
['a', 'b']
['sample_a']
True
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 4) 整理、标记重复并索引 BAM {#bio-16}

- **作用：** 把 SAM/BAM 整理为可用于变异检测的 BAM。
- **计算方法：** samtools 依次执行 name-sort、fixmate、坐标排序、markdup、index 和 stats。
- **API：** `dnakit.bioinformatics.ngs.prepare_bam(source, output_directory, mark_duplicates=True, config=None)`
- **输入：** 已有 `aligned.sam` 或 BAM；samtools 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import prepare_bam

result = prepare_bam("aligned.sam", "new_bam")
print(result.artifacts["aligned.bam"])
```

- **示例结果：**

返回 `ToolResult`，含 `aligned.bam`、索引和统计文件。外部程序示例未执行。

[依赖与实现说明](../../faq/bioinformatics.md#reads)。
