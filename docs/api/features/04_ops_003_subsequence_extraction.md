# 子序列提取

从 DNA 序列中提取指定片段。

## 1) 提取子序列

- **作用：** 按起止坐标截取 DNA 片段并保留来源范围，也支持环状序列跨原点截取，用于提取基因、功能区域或指定窗口。
- **计算方法：** 使用 0-based 半开区间 [start,end) 截取；环状跨原点片段按尾段与首段顺序拼接。
- **API：** `dnakit.ops.subsequence(sequence[必须], start[必须], end[必须], allow_gaps[可选])`、`dnakit.ops.circular_subsequence(sequence[必须], start[必须], end[必须], allow_gaps[可选])`。
- **输入：** 必填 `DNASequence`、`start`、`end`；可选 `allow_gaps=True` 保留或截短已知 Gap。
- **输出：** 返回处理后的 `DNASequence`，序列内容可从 `symbols` 读取。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.ops import subsequence

seq = DNASequence("AACCGG")
selected = subsequence(seq, 1, 5)
print(selected.symbols)
```

- **示例结果：**

```text
ACCG
```

## 2) 读取区间并随机提取参考序列 {#bio-13}

- **作用：** 根据 BED/GTF 区间读取参考 FASTA，支持负链反向互补。
- **计算方法：** 将区间统一为 0-based 半开坐标，用 pysam 的 FASTA 索引随机读取。
- **API：** `dnakit.bioinformatics.reads.read_regions(path, format='bed', feature='exon')`; `dnakit.bioinformatics.reads.index_reference(path)`; `dnakit.bioinformatics.reads.IndexedReference(path).fetch(region)`
- **输入：** 参考 FASTA、BED/GTF 或 `Region`；需要 `bio-reads`。`index_reference` 在参考文件旁创建 `.fai`。
- **输出：** `read_regions` 返回 `Region` 迭代器；`index_reference` 生成 `.fai` 索引；`IndexedReference.fetch` 返回含目标区间序列的 `Read`。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import read_regions, index_reference, IndexedReference

with TemporaryDirectory() as folder:
    fasta, bed = Path(folder) / "ref.fa", Path(folder) / "regions.bed"
    fasta.write_text(">chr1\nAACCGGTT\n")
    bed.write_text("chr1\t0\t4\tregion1\t0\t-\n")
    index_reference(fasta)
    region = next(read_regions(bed))
    with IndexedReference(fasta) as reference:
        print(region.start, region.end, reference.fetch(region).sequence)
```

- **示例结果：**

```text
0 4 GGTT
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 3) 求区间并集、交集或差集 {#bio-60}

- **作用：** 组合两个区间集合，提取重叠或独有部分。
- **计算方法：** 按染色体和链分组，扫描区间端点并合并相邻片段。
- **API：** `dnakit.bioinformatics.workflows.query_regions(left, right, operation='intersection')`
- **输入：** 两个 `Region` 集合；所有坐标使用 0-based 半开区间。
- **输出：** 返回区间元组，每个区间包含参考名、起点和终点。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Region
from dnakit.bioinformatics.workflows import query_regions

result = query_regions([Region("chr1", 0, 10)], [Region("chr1", 5, 15)])
print([(region.reference, region.start, region.end) for region in result])
```

- **示例结果：**

```text
[('chr1', 5, 10)]
```

[依赖与实现说明](../../faq/bioinformatics.md#workflows)。
