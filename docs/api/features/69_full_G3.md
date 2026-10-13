# 第 69 号：G3 区间、注释与 bigWig 包装

本文件记录 G3 组的 24 个 DNAKit 包装函数，分别来自 bioframe（6）、pybedtools（4）、gffutils（4）、gtfparse（3）、pyBigWig（3）和 pyranges1（4）。这些函数都从 `dnakit.annotation` 导出，并且只处理显式给出的本地文件，不访问网络。

例如 `from dnakit.annotation import overlap_bioframe_intervals`；API 行给出完整签名。pybedtools 相关函数需要 `bedtools` 可执行文件（本环境为 `/usr/bin/bedtools`）。可选依赖只在调用时导入，缺少时会给出安装提示。示例输出由实际运行生成。

## 示例数据

下面的示例共用这些文件，每个文件的用途如下。pyBigWig 的示例在运行时把 bedGraph 写成临时 bigWig，因此无需额外文件。

BED6 intervals for bioframe, pybedtools, and pyranges1 examples.

**示例文件：** `bf_a.bed`

```text
chr1	10	20	a	5	+
chr1	15	30	b	7	+
chr1	40	50	c	1	-
chr2	5	9	d	2	+
```

BED4 intervals used as the second input of overlap-type examples.

**示例文件：** `bf_b.bed`

```text
chr1	15	18	x
chr1	45	60	y
chr2	100	120	z
```

One BEDPE record (10 columns) for the tagging example.

**示例文件：** `pbt_pairs.bedpe`

```text
chr1	10	20	chr1	100	110	p1	0	+	-
```

BED4 query intervals named q1 and q2.

**示例文件：** `pbt_q.bed`

```text
chr1	12	18	q1
chr1	101	105	q2
```

Small GFF3 file: one gene with one mRNA and two exons, and one more gene.

**示例文件：** `gff_genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	exon	300	500	.	+	.	ID=exon2;Parent=tx1
chr2	src	gene	50	90	.	-	.	ID=gene2
```

Small GTF file: gene, transcript, and one exon.

**示例文件：** `gtp_genes.gtf`

```text
chr1	src	gene	100	500	.	+	.	gene_id "g1"; gene_name "GENE1";
chr1	src	transcript	100	500	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	200	.	+	.	gene_id "g1"; transcript_id "t1"; exon_number 1;
```

bedGraph rows (chrom, start, end, value) for the bigWig examples.

**示例文件：** `bw_cov.bedgraph`

```text
chr1	0	10	1.5
chr1	10	20	3
chr2	5	9	2
```

## 1) `dnakit.annotation.overlap_bioframe_intervals()`

- **作用：** 找出两个 BED 文件中相互重叠的区间对。
- **处理方法：** 调用 `bioframe.overlap`。`how` 为 `inner` 时只保留有配对的行；`left`、`right`、`outer` 同时保留未配对的行，其配对列为空。
- **API：** `dnakit.annotation.overlap_bioframe_intervals(source: 'str | Path', other: 'str | Path', *, how: 'str' = 'inner', max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 两个 BED3 到 BED6 文件；`how` 取 inner、left、right 或 outer；`max_rows` 限制每个文件的行数。
- **输出：** `records` 列表，每项是一对区间。左侧列无后缀，右侧列以 `_` 结尾；另含 `record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import overlap_bioframe_intervals

result = overlap_bioframe_intervals("bf_a.bed", "bf_b.bed")
print(result["record_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
3
[('chr1', 10, 20), ('chr1', 15, 30), ('chr1', 40, 50)]
```

## 2) `dnakit.annotation.subtract_bioframe_intervals()`

- **作用：** 从每个区间中扣除另一文件覆盖的碱基，返回剩余片段。
- **处理方法：** 调用 `bioframe.subtract`。被扣除的部分不再出现；一个区间可能被拆成多段。
- **API：** `dnakit.annotation.subtract_bioframe_intervals(source: 'str | Path', other: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 两个 BED3 到 BED6 文件；`max_rows` 限制每个文件的行数。
- **输出：** `intervals` 列表（chrom、start、end 等列），`interval_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import subtract_bioframe_intervals

result = subtract_bioframe_intervals("bf_a.bed", "bf_b.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 10, 15), ('chr1', 18, 20), ('chr1', 18, 30), ('chr1', 40, 45), ('chr2', 5, 9)]
```

## 3) `dnakit.annotation.expand_bioframe_intervals()`

- **作用：** 把每个区间向两侧各扩展 `pad` 个碱基。
- **处理方法：** 调用 `bioframe.expand`。起点减去 `pad`，终点加上 `pad`；结果不会截断到 0。
- **API：** `dnakit.annotation.expand_bioframe_intervals(source: 'str | Path', *, pad: 'int', max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 BED3 到 BED6 文件；`pad` 为非负整数；`max_rows` 限制输入行数。
- **输出：** `intervals` 列表与 `interval_count`，`pad` 原样返回，并附 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import expand_bioframe_intervals

result = expand_bioframe_intervals("bf_a.bed", pad=5)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 5, 25), ('chr1', 10, 35), ('chr1', 35, 55), ('chr2', 0, 14)]
```

## 4) `dnakit.annotation.cluster_bioframe_intervals()`

- **作用：** 给每个区间分配簇号；重叠或间距不超过 `min_dist` 的区间属于同一簇。
- **处理方法：** 调用 `bioframe.cluster`。同一染色体上的区间按顺序合并成簇，并返回 `cluster`、`cluster_start`、`cluster_end` 列。
- **API：** `dnakit.annotation.cluster_bioframe_intervals(source: 'str | Path', *, min_dist: 'int' = 0, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 BED3 到 BED6 文件；`min_dist` 为非负整数；`max_rows` 限制输入行数。
- **输出：** `records` 列表（每个输入区间一行，含簇号），`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import cluster_bioframe_intervals

result = cluster_bioframe_intervals("bf_a.bed")
print([(row["start"], row["end"], row["cluster"]) for row in result["records"]])
```

- **示例输出：**

```text
[(10, 20, 0), (15, 30, 0), (40, 50, 1), (5, 9, 2)]
```

## 5) `dnakit.annotation.binnify_bioframe_chromosomes()`

- **作用：** 把每条染色体按固定宽度切成 bin，用于基因组分箱。
- **处理方法：** 调用 `bioframe.binnify`。每条染色体从 0 开始切分；最后一个 bin 可能比 `binsize` 短。
- **API：** `dnakit.annotation.binnify_bioframe_chromosomes(chromsizes: 'Mapping[str, int]', *, binsize: 'int', max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 染色体长度映射（名称到长度）；`binsize` 为正整数；`max_rows` 限制 bin 数量。
- **输出：** `bins` 列表（chrom、start、end），`bin_count`、`binsize` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import binnify_bioframe_chromosomes

result = binnify_bioframe_chromosomes({"chr1": 25, "chr2": 10}, binsize=10)
print(result["bin_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["bins"]])
```

- **示例输出：**

```text
4
[('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 25), ('chr2', 0, 10)]
```

## 6) `dnakit.annotation.parse_bioframe_region()`

- **作用：** 解析 UCSC 样式的区域字符串，例如 `chr1:101-200`。
- **处理方法：** 调用 `bioframe.parse_region_string`。千位分隔逗号会被去掉；数字原样返回，不做 1 基到 0 基的转换。
- **API：** `dnakit.annotation.parse_bioframe_region(region: 'str', *, verbose: bool = True)`
- **输入：** 一个字符串（长度不超过 200）。
- **输出：** `chrom`、`start`、`end` 字段，以及 `provenance`。无法解析时抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.annotation import parse_bioframe_region

result = parse_bioframe_region("chr1:1,101-2,000")
print(result["chrom"], result["start"], result["end"])
```

- **示例输出：**

```text
chr1 1101 2000
```

## 7) `dnakit.annotation.five_prime_pybedtools_coordinates()`

- **作用：** 返回每个特征 5' 端两侧的窗口。
- **处理方法：** 调用 `pybedtools.featurefuncs.five_prime`，通过 `BedTool.each` 应用到每个区间。+ 链的 5' 端是起点，- 链的 5' 端是终点。
- **API：** `dnakit.annotation.five_prime_pybedtools_coordinates(source: 'str | Path', *, upstream: 'int' = 500, downstream: 'int' = 500, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 BED3 到 BED6 文件；`upstream` 与 `downstream` 为非负整数；`max_rows` 限制输入行数。需要 BEDTools 可执行文件。
- **输出：** `records` 列表（chrom、start、end、name、strand），`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import five_prime_pybedtools_coordinates

result = five_prime_pybedtools_coordinates("bf_a.bed", upstream=5, downstream=5)
print([(row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
[(5, 15), (10, 20), (45, 55), (0, 10)]
```

## 8) `dnakit.annotation.midpoint_pybedtools_coordinates()`

- **作用：** 返回每个 BED 特征的 1 bp 中点。
- **处理方法：** 调用 `pybedtools.featurefuncs.midpoint`，中点为 `floor((start + end) / 2)`。
- **API：** `dnakit.annotation.midpoint_pybedtools_coordinates(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 BED3 到 BED6 文件；`max_rows` 限制输入行数。需要 BEDTools 可执行文件。
- **输出：** `records` 列表（每项是一个 1 bp 区间），`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import midpoint_pybedtools_coordinates

result = midpoint_pybedtools_coordinates("bf_a.bed")
print([(row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
[(15, 16), (22, 23), (45, 46), (7, 8)]
```

## 9) `dnakit.annotation.tag_pybedtools_bedpe()`

- **作用：** 标记 BEDPE 每一端与哪些命名查询区间重叠。
- **处理方法：** 调用 `pybedtools.contrib.long_range_interaction.tag_bedpe`。两端分别与查询 BED4 求交，输出每个查询名和重叠碱基数，无重叠的行被丢弃。
- **API：** `dnakit.annotation.tag_pybedtools_bedpe(bedpe: 'str | Path', queries: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 10 列的 BEDPE 文件和一个 BED4 查询文件（第 4 列为名称）。需要 BEDTools 可执行文件。
- **输出：** `pairs` 列表，每项含 `name`、`left`、`right`（每端的 `query` 与 `overlap_bp`），`pair_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import tag_pybedtools_bedpe

result = tag_pybedtools_bedpe("pbt_pairs.bedpe", "pbt_q.bed")
print(result["pairs"])
```

- **示例输出：**

```text
[{'name': 'p1', 'left': [{'query': 'q1', 'overlap_bp': 6}], 'right': [{'query': 'q2', 'overlap_bp': 4}]}]
```

## 10) `dnakit.annotation.write_pybedtools_chromsizes()`

- **作用：** 把染色体长度映射写成两列的 chromsizes 文件。
- **处理方法：** 调用 `pybedtools.chromsizes_to_file`；行按染色体名排序，保证输出确定。输出路径必须是新文件。
- **API：** `dnakit.annotation.write_pybedtools_chromsizes(chromsizes: 'dict[str, int]', output: 'str | Path', *, verbose: bool = True)`
- **输入：** 染色体长度映射（名称到长度）；`output` 为新文件路径，所在目录必须存在。
- **输出：** `path`、`chromosomes` 列表（chrom、length）、`chromosome_count` 与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_pybedtools_chromsizes

out = Path(tempfile.mkdtemp()) / "sizes.txt"
result = write_pybedtools_chromsizes({"chr2": 50, "chr1": 100}, out)
print(result["chromosomes"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'length': 100}, {'chrom': 'chr2', 'length': 50}]
```

## 11) `dnakit.annotation.to_bed12_gffutils_transcript()`

- **作用：** 把一个转录本及其子特征转换为一行 BED12。
- **处理方法：** 先从数据库取出转录本，再调用 `gffutils.convert.to_bed12`。BED 起点为 0 基，块来自 `child_type` 指定的子特征（默认 exon）。
- **API：** `dnakit.annotation.to_bed12_gffutils_transcript(database: 'str | Path', transcript_id: 'str', *, child_type: 'str' = 'exon', name_field: 'str' = 'ID', verbose: bool = True)`
- **输入：** gffutils 数据库路径（需已存在）；`transcript_id` 为转录本 ID；`child_type` 与 `name_field` 为字符串。
- **输出：** `bed12` 字符串（制表符分隔，不含换行）、`transcript_id` 与 `provenance`。未找到 ID 时抛出 `ConfigurationError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import gffutils

db = Path(tempfile.mkdtemp()) / "genes.db"
gffutils.create_db("gff_genes.gff3", str(db), force=True, verbose=False, merge_strategy="create_unique")
from dnakit.annotation import to_bed12_gffutils_transcript

print(to_bed12_gffutils_transcript(db, "tx1")["bed12"])
```

- **示例输出：**

```text
chr1	99	500	tx1	.	+	100	500	0,0,0	2	101,201	0,200
```

## 12) `dnakit.annotation.tss_gffutils_intervals()`

- **作用：** 列出数据库中每个转录本的 1 bp 转录起始位点（TSS）。
- **处理方法：** 调用 `gffutils.pybedtools_integration.tsses`，并以 BED6 输出。起点为 0 基。
- **API：** `dnakit.annotation.tss_gffutils_intervals(database: 'str | Path', *, verbose: bool = True)`
- **输入：** gffutils 数据库路径（需已存在）。需要 pybedtools 与 BEDTools 可执行文件。
- **输出：** `records` 列表（chrom、start、end、name、strand），`record_count` 与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import gffutils

db = Path(tempfile.mkdtemp()) / "genes.db"
gffutils.create_db("gff_genes.gff3", str(db), force=True, verbose=False, merge_strategy="create_unique")
from dnakit.annotation import tss_gffutils_intervals

print(tss_gffutils_intervals(db)["records"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 99, 'end': 100, 'name': 'tx1', 'strand': '+'}]
```

## 13) `dnakit.annotation.parse_gffutils_feature_line()`

- **作用：** 把一行 GFF3 记录解析为特征字段和属性。
- **处理方法：** 调用 `gffutils.feature.feature_from_line`（严格模式）。坐标按 GFF 原样保留（1 基、闭区间）。
- **API：** `dnakit.annotation.parse_gffutils_feature_line(line: 'str', *, verbose: bool = True)`
- **输入：** 一行制表符分隔的 9 列 GFF3 文本（不超过 65536 字符）。
- **输出：** `seqid`、`source`、`featuretype`、`start`、`end`、`score`、`strand`、`frame`、`attributes`（键排序，值为列表）与 `provenance`。格式错误时抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.annotation import parse_gffutils_feature_line

line = "chr1\tsrc\texon\t300\t500\t.\t+\t.\tID=exon2;Parent=tx1"
result = parse_gffutils_feature_line(line)
print(result["featuretype"], result["start"], result["end"], result["attributes"])
```

- **示例输出：**

```text
exon 300 500 {'ID': ['exon2'], 'Parent': ['tx1']}
```

## 14) `dnakit.annotation.inspect_gffutils_source()`

- **作用：** 统计 GFF 或 GTF 文件中的特征类型、染色体、属性键和特征数。
- **处理方法：** 调用 `gffutils.inspect.inspect`（不输出日志）。计数字典按键排序。
- **API：** `dnakit.annotation.inspect_gffutils_source(path: 'str | Path', *, verbose: bool = True)`
- **输入：** 一个 GFF 或 GTF 文件路径。
- **输出：** `featuretype`、`chrom`、`attribute_keys` 计数字典，`feature_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import inspect_gffutils_source

result = inspect_gffutils_source("gff_genes.gff3")
print(result["featuretype"], result["feature_count"])
```

- **示例输出：**

```text
{'exon': 2, 'gene': 2, 'mRNA': 1} 5
```

## 15) `dnakit.annotation.parse_gtfparse_fixed_columns()`

- **作用：** 读取 GTF 的 9 个固定列，属性列保留原始文本。
- **处理方法：** 调用 `gtfparse.parse_gtf`（`split_attributes=False`）。坐标为 1 基、闭区间。
- **API：** `dnakit.annotation.parse_gtfparse_fixed_columns(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 GTF 文件；`max_rows` 限制行数。
- **输出：** `records` 列表（seqname、source、feature、start、end、score、strand、frame、attribute），`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import parse_gtfparse_fixed_columns

result = parse_gtfparse_fixed_columns("gtp_genes.gtf")
print(result["record_count"], result["records"][2]["feature"], result["records"][2]["start"])
```

- **示例输出：**

```text
3 exon 100
```

## 16) `dnakit.annotation.expand_gtfparse_attributes()`

- **作用：** 读取 GTF 并把属性列展开为每个键一列。
- **处理方法：** 调用 `gtfparse.parse_gtf_and_expand_attributes`。缺失的属性键记为 None。
- **API：** `dnakit.annotation.expand_gtfparse_attributes(source: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 一个 GTF 文件；`max_rows` 限制行数。
- **输出：** `columns`（排序后的列名）、`records`、`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import expand_gtfparse_attributes

result = expand_gtfparse_attributes("gtp_genes.gtf")
print(result["columns"])
print(result["records"][2]["transcript_id"], result["records"][2]["exon_number"])
```

- **示例输出：**

```text
['end', 'exon_number', 'feature', 'frame', 'gene_id', 'gene_name', 'score', 'seqname', 'source', 'start', 'strand', 'transcript_id']
t1 1
```

## 17) `dnakit.annotation.write_gtfparse_table()`

- **作用：** 解析 GTF 后写出为一个新的 GTF 文件。
- **处理方法：** 调用 `gtfparse.parse_gtf` 与 `gtfparse.write_gtf`。`frame` 按解析结果写出，`.` 写回后为 0。输出路径必须是新文件。
- **API：** `dnakit.annotation.write_gtfparse_table(source: 'str | Path', output: 'str | Path', *, verbose: bool = True)`
- **输入：** 源 GTF 文件；`output` 为新文件路径，所在目录必须存在。
- **输出：** `path`、`line_count` 与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_gtfparse_table

out = Path(tempfile.mkdtemp()) / "round.gtf"
result = write_gtfparse_table("gtp_genes.gtf", out)
print(result["line_count"])
print(out.read_text(encoding="utf-8").splitlines()[2])
```

- **示例输出：**

```text
3
chr1	src	exon	100	200	.	+	0	gene_id "g1"; transcript_id "t1"; exon_number 1;
```

## 18) `dnakit.annotation.read_bigwig_header()`

- **作用：** 读取 bigWig 文件的染色体长度与头部摘要。
- **处理方法：** 调用 `pyBigWig.open` 打开只读句柄，读取 `chroms()` 与 `header()`。非 bigWig 文件会被拒绝。
- **API：** `dnakit.annotation.read_bigwig_header(path: 'str | Path', *, verbose: bool = True)`
- **输入：** 一个 bigWig 文件路径。
- **输出：** `chromosomes` 列表（chrom、length，按名称排序）、`chromosome_count`、`header` 字典与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import read_bigwig_header, write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
result = read_bigwig_header(out)
print(result["chromosomes"])
print(result["header"]["nBasesCovered"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'length': 100}, {'chrom': 'chr2', 'length': 50}]
24
```

## 19) `dnakit.annotation.write_bigwig_intervals()`

- **作用：** 由 bedGraph 行（chrom、start、end、value）写出一个新的 bigWig 文件。
- **处理方法：** 调用 `pyBigWig.open(path, 'w')`，用 `addHeader` 写入染色体长度，再用 `addEntries` 写入区间。区间必须互不重叠且位于染色体范围内。
- **API：** `dnakit.annotation.write_bigwig_intervals(bedgraph: 'str | Path', output: 'str | Path', *, chromsizes: 'Mapping[str, int]', verbose: bool = True)`
- **输入：** bedGraph 文件；`chromsizes` 为长度映射，须包含所有出现的染色体；`output` 为新文件路径。
- **输出：** `path`、`interval_count`、`chromosomes` 与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
result = write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
print(result["interval_count"], result["chromosomes"])
```

- **示例输出：**

```text
3 ['chr1', 'chr2']
```

## 20) `dnakit.annotation.summarize_bigwig_bins()`

- **作用：** 把 bigWig 的一段区域按等宽分箱并计算统计量。
- **处理方法：** 调用 `pyBigWig` 句柄的 `stats(chrom, start, end, type=statistic, nBins=bins)`。空 bin 返回 None。
- **API：** `dnakit.annotation.summarize_bigwig_bins(path: 'str | Path', *, chrom: 'str', start: 'int', end: 'int', bins: 'int', statistic: 'str' = 'mean', verbose: bool = True)`
- **输入：** bigWig 路径；`chrom`、`start`、`end`（须在染色体范围内）；`bins` 为正整数；`statistic` 取 coverage、max、mean、min、std 或 sum。
- **输出：** `values` 列表（每个 bin 一项）、`chrom`、`start`、`end`、`bins`、`statistic` 与 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import summarize_bigwig_bins, write_bigwig_intervals

out = Path(tempfile.mkdtemp()) / "cov.bw"
write_bigwig_intervals("bw_cov.bedgraph", out, chromsizes={"chr1": 100, "chr2": 50})
result = summarize_bigwig_bins(out, chrom="chr1", start=0, end=20, bins=2)
print(result["values"])
```

- **示例输出：**

```text
[1.5, 3.0]
```

## 21) `dnakit.annotation.tile_pyranges1_genome()`

- **作用：** 把染色体长度映射切成固定宽度的 tile。
- **处理方法：** 调用 `pyranges1.tile_genome`。每条染色体从 0 开始切分；`full_last_tile` 原样传给 pyranges1。
- **API：** `dnakit.annotation.tile_pyranges1_genome(chromsizes: 'Mapping[str, int]', *, tile_size: 'int', full_last_tile: 'bool' = False, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 染色体长度映射；`tile_size` 为正整数；`full_last_tile` 为布尔值；`max_rows` 限制 tile 数量。
- **输出：** `tiles` 列表（chrom、start、end）、`tile_count`、`tile_size` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import tile_pyranges1_genome

result = tile_pyranges1_genome({"chr1": 25, "chr2": 10}, tile_size=10)
print(result["tile_count"])
print([(row["chrom"], row["start"], row["end"]) for row in result["tiles"]])
```

- **示例输出：**

```text
4
[('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 25), ('chr2', 0, 10)]
```

## 22) `dnakit.annotation.count_pyranges1_overlaps()`

- **作用：** 统计另一文件中与每个源区间重叠的区间数。
- **处理方法：** 调用 `pyranges1.count_overlaps`，以源区间为 features，另一文件为计数对象。
- **API：** `dnakit.annotation.count_pyranges1_overlaps(source: 'str | Path', other: 'str | Path', *, max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 两个 BED 文件；`max_rows` 限制每个文件的行数。
- **输出：** `records` 列表（chrom、start、end、strand 与 `overlaps` 计数），`record_count` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import count_pyranges1_overlaps

result = count_pyranges1_overlaps("bf_a.bed", "bf_b.bed")
print([(row["chrom"], row["start"], row["overlaps"]) for row in result["records"]])
```

- **示例输出：**

```text
[('chr1', 10, 1), ('chr1', 15, 1), ('chr1', 40, 1), ('chr2', 5, 0)]
```

## 23) `dnakit.annotation.parse_pyranges1_string()`

- **作用：** 解析带表头的多行表格文本为区间。
- **处理方法：** 调用 `pyranges1.from_string`。第一行为表头，必须含 Chromosome、Start、End；Strand 可选。
- **API：** `dnakit.annotation.parse_pyranges1_string(text: 'str', *, verbose: bool = True)`
- **输入：** 一段文本（不超过 65536 字符）。
- **输出：** `records` 列表（chrom、start、end，有 Strand 列时含 strand）、`record_count` 与 `provenance`。无法解析时抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.annotation import parse_pyranges1_string

text = "Chromosome Start End Strand\nchr1 10 20 +\nchr1 40 50 -\n"
result = parse_pyranges1_string(text)
print(result["records"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 10, 'end': 20, 'strand': '+'}, {'chrom': 'chr1', 'start': 40, 'end': 50, 'strand': '-'}]
```

## 24) `dnakit.annotation.jaccard_pyranges1_intervals()`

- **作用：** 计算两个 BED 文件按碱基计的 Jaccard 指数。
- **处理方法：** 调用 `pyranges1.stats.jaccard`。分子为两者共有的碱基数，分母为任一文件覆盖的碱基数。
- **API：** `dnakit.annotation.jaccard_pyranges1_intervals(source: 'str | Path', other: 'str | Path', *, strand_behavior: 'str' = 'auto', max_rows: 'int' = 100000, verbose: bool = True)`
- **输入：** 两个 BED 文件；`strand_behavior` 取 auto、same、opposite 或 ignore；`max_rows` 限制每个文件的行数。
- **输出：** `jaccard`（浮点数）、`strand_behavior` 与 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import jaccard_pyranges1_intervals

result = jaccard_pyranges1_intervals("bf_a.bed", "bf_b.bed")
print(round(result["jaccard"], 4))
```

- **示例输出：**

```text
0.125
```
