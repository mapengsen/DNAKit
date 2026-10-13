# 可选区间、注释与区域集后端

这一组函数是 `dnakit.annotation` 中基于可选依赖的适配器，覆盖 [bioframe](https://github.com/open2c/bioframe)、[pybedtools](https://github.com/daler/pybedtools)、[gtfparse](https://github.com/openvax/gtfparse)、[gffutils](https://github.com/daler/gffutils)、[pyBigWig](https://github.com/deeptools/pyBigWig)、[pyranges1](https://github.com/pyranges/pyranges1)、[polars-bio](https://github.com/biodatageeks/polars-bio)、[gtars](https://github.com/databio/gtars)、[ncls](https://github.com/pyranges/ncls) 和 [genomicranges](https://github.com/BiocPy/genomicranges)。它们计算 BED 区间的合并、最近邻、覆盖度、交集、补集和序列提取，GTF/GFF 的解析与数据库查询，bigWig 数值读取，区域集的 Jaccard 相似度、缺口与分词，以及区间索引查询。所有函数只处理显式给出的本地文件，不访问网络，也不做序列比对、峰检测或注释数据库下载。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令；pybedtools 相关函数还需要单独安装 `bedtools` 可执行文件。

## 示例数据

下面的示例共用这些文件。`a.bed` 和 `b.bed` 是两组小区间，`seqs.bed` 与 `ref.fa` 用于序列提取，`genes.gff3` 与 `genes.gtf` 用于注释，`universe.bed` 是分词用的 BED3 区域集。bigWig 文件和 gffutils 数据库由示例代码在运行时生成。

**示例文件：** `a.bed`

```text
chr1	10	20	a	5	+
chr1	15	30	b	7	+
chr1	40	50	c	1	-
chr2	5	9	d	2	+
```

**示例文件：** `b.bed`

```text
chr1	12	18	x
chr1	45	60	y
chr2	100	120	z
```

**示例文件：** `seqs.bed`

```text
chr1	10	20	p1	0	+
chr1	0	5	p2	0	-
chr2	5	9	p3	0	+
```

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTAC
GGGGTTTTCC
>chr2
ATATATATAT
```

**示例文件：** `universe.bed`

```text
chr1	12	18
chr1	45	60
chr2	100	120
```

**示例文件：** `genes.gff3`

```text
##gff-version 3
chr1	src	gene	100	500	.	+	.	ID=gene1;Name=GENE1
chr1	src	mRNA	100	500	.	+	.	ID=tx1;Parent=gene1;Name=GENE1-201
chr1	src	exon	100	200	.	+	.	ID=exon1;Parent=tx1
chr1	src	exon	300	500	.	+	.	ID=exon2;Parent=tx1
chr1	src	gene	800	900	.	-	.	ID=gene2;Name=GENE2
chr1	src	mRNA	800	900	.	-	.	ID=tx2;Parent=gene2
chr1	src	exon	800	900	.	-	.	ID=exon3;Parent=tx2
```

**示例文件：** `genes.gtf`

```text
chr1	ens	gene	100	500	.	+	.	gene_id "ENSG1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1	ens	transcript	100	500	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; transcript_version "2";
chr1	ens	exon	100	200	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1	ens	exon	300	500	.	+	.	gene_id "ENSG1"; transcript_id "ENST1"; exon_number "2"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1	ens	gene	800	900	.	-	.	gene_id "ENSG2"; gene_name "GENE2"; gene_biotype "lncRNA";
chr1	ens	transcript	800	900	.	-	.	gene_id "ENSG2"; transcript_id "ENST2"; gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";
chr1	ens	exon	800	900	.	-	.	gene_id "ENSG2"; transcript_id "ENST2"; exon_number "1"; gene_name "GENE2"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";
```

## 1) `dnakit.annotation.merge_bioframe_intervals()`

- **作用：** 合并重叠或相邻的 BED 区间。
- **处理方法：** 调用 `bioframe.merge`。`min_distance` 为 0 时，相接的区间也会合并；大于 0 时，间隔不超过该值的区间会合并。
- **API：** `dnakit.annotation.merge_bioframe_intervals(source, *, min_distance=0, max_rows=100000)`。
- **输入：** BED3 到 BED6 文件；`min_distance` 为非负整数；`max_rows` 限制输入行数。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `n_intervals`，并附 `provenance`。
- **示例代码：**

```python
from dnakit.annotation import merge_bioframe_intervals

result = merge_bioframe_intervals("a.bed")
print([(row["chrom"], row["start"], row["end"], row["n_intervals"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 10, 30, 2), ('chr1', 40, 50, 1), ('chr2', 5, 9, 1)]
```

## 2) `dnakit.annotation.closest_bioframe_intervals()`

- **作用：** 为每个源区间找到另一文件中最近的 k 个区间。
- **处理方法：** 调用 `bioframe.closest`。距离按半开区间计算，即 `start_b - end_a`；重叠时为 0。`ignore_overlaps=True` 时跳过重叠区间。
- **API：** `dnakit.annotation.closest_bioframe_intervals(source, other, *, k=1, ignore_overlaps=False, max_rows=100000)`。
- **输入：** 两个 BED 文件；`k` 取 1 到 100；`ignore_overlaps` 为布尔值。
- **输出：** `records` 列表，每项包含源区间和最近区间的列（带 `_` 后缀的为 `other` 的列）以及 `distance`。
- **示例代码：**

```python
from dnakit.annotation import closest_bioframe_intervals

result = closest_bioframe_intervals("a.bed", "b.bed", k=1)
print([(row["name"], row["name_"], row["distance"]) for row in result["records"]])
```

- **示例输出：**

```text
[('a', 'x', 0), ('b', 'x', 0), ('c', 'y', 0), ('d', 'z', 91)]
```

## 3) `dnakit.annotation.coverage_bioframe_intervals()`

- **作用：** 统计另一文件覆盖了源区间中的多少碱基。
- **处理方法：** 调用 `bioframe.coverage`，结果列为 `coverage`，单位是碱基数。
- **API：** `dnakit.annotation.coverage_bioframe_intervals(source, other, *, max_rows=100000)`。
- **输入：** 两个 BED 文件。
- **输出：** `records` 列表，每个源区间附带 `coverage`。
- **示例代码：**

```python
from dnakit.annotation import coverage_bioframe_intervals

result = coverage_bioframe_intervals("a.bed", "b.bed")
print([(row["name"], row["coverage"]) for row in result["records"]])
```

- **示例输出：**

```text
[('a', 6), ('b', 3), ('c', 5), ('d', 0)]
```

## 4) `dnakit.annotation.merge_pybedtools_intervals()`

- **作用：** 用 BEDTools 合并 BED 区间。
- **处理方法：** 先由 pybedtools 排序，再调用 `bedtools merge -d`。距离为 `distance` 时，间隔不超过它的区间会合并。输出只保留 chrom、start 和 end 三列。
- **API：** `dnakit.annotation.merge_pybedtools_intervals(source, *, distance=0, max_rows=100000)`。
- **输入：** BED 文件；`distance` 为非负整数。需要 `bedtools` 可执行文件。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end`。
- **示例代码：**

```python
from dnakit.annotation import merge_pybedtools_intervals

result = merge_pybedtools_intervals("a.bed", distance=10)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 10, 50), ('chr2', 5, 9)]
```

## 5) `dnakit.annotation.closest_pybedtools_intervals()`

- **作用：** 为每个源区间找到另一文件中最近的 k 个区间，并返回 BEDTools 的距离。
- **处理方法：** 调用 `bedtools closest -d -t first`。与 `bioframe` 不同，BEDTools 对不重叠区间的距离加 1，即 `start_b - end_a + 1`。
- **API：** `dnakit.annotation.closest_pybedtools_intervals(source, other, *, k=1, max_rows=100000)`。
- **输入：** 两个 BED 文件；`k` 取 1 到 100。需要 `bedtools` 可执行文件。
- **输出：** `records` 列表，每项含 `a`（源区间字段）、`b`（最近区间字段，缺失时为 `None`）和 `distance`。
- **示例代码：**

```python
from dnakit.annotation import closest_pybedtools_intervals

result = closest_pybedtools_intervals("a.bed", "b.bed", k=1)
print([(row["a"][3], row["b"][3], row["distance"]) for row in result["records"]])
```

- **示例输出：**

```text
[('a', 'x', 0), ('b', 'x', 0), ('c', 'y', 0), ('d', 'z', 92)]
```

## 6) `dnakit.annotation.coverage_pybedtools_intervals()`

- **作用：** 统计源区间被另一文件覆盖的碱基数和比例。
- **处理方法：** 调用 `bedtools coverage`。`fraction` 为覆盖碱基数除以区间长度。
- **API：** `dnakit.annotation.coverage_pybedtools_intervals(source, other, *, max_rows=100000)`。
- **输入：** 两个 BED 文件。需要 `bedtools` 可执行文件。
- **输出：** `records` 列表，每项含 `a`、`overlaps`、`bases_covered`、`length` 和 `fraction`。
- **示例代码：**

```python
from dnakit.annotation import coverage_pybedtools_intervals

result = coverage_pybedtools_intervals("a.bed", "b.bed")
print([(row["a"][3], row["bases_covered"], row["fraction"]) for row in result["records"]])
```

- **示例输出：**

```text
[('a', 6, 0.6), ('b', 3, 0.2), ('c', 5, 0.5), ('d', 0, 0.0)]
```

## 7) `dnakit.annotation.extract_pybedtools_sequences()`

- **作用：** 按 BED 区间从 FASTA 中取出序列。
- **处理方法：** 调用 `bedtools getfasta -tab`。`use_strand=True` 时，负链区间取反向互补序列。区间超出序列长度时，bedtools 会跳过它，本函数将其报错而不是静默丢弃。如果 FASTA 没有 `.fai` 索引，bedtools 会在同目录生成 `<fasta>.fai`。
- **API：** `dnakit.annotation.extract_pybedtools_sequences(source, fasta, *, use_strand=False, max_intervals=1000, max_sequence_length=100000)`。
- **输入：** BED 文件和 FASTA 文件；每个区间的长度不超过 `max_sequence_length`。需要 `bedtools` 可执行文件。
- **输出：** `sequences` 列表，每项含 `name`（`chrom:start-end`）、`sequence` 和 `length`。
- **示例代码：**

```python
from dnakit.annotation import extract_pybedtools_sequences

result = extract_pybedtools_sequences("seqs.bed", "ref.fa", use_strand=True)
print([(row["name"], row["sequence"]) for row in result["sequences"]])
```

- **示例输出：**

```text
[('chr1:10-20(+)', 'GGGGTTTTCC'), ('chr1:0-5(-)', 'TACGT'), ('chr2:5-9(+)', 'TATA')]
```

## 8) `dnakit.annotation.read_gtfparse_table()`

- **作用：** 把 GTF 解析为记录列表，每个属性键展开为一列。
- **处理方法：** 调用 `gtfparse.read_gtf`。`features` 只保留指定的特征类型。`normalize_gencode_biotypes=True` 时，GENCODE 的 `gene_type` 和 `transcript_type` 会改名为 `gene_biotype` 和 `transcript_biotype`。gtfparse 中缺失的属性以 `None` 表示。
- **API：** `dnakit.annotation.read_gtfparse_table(source, *, features=(), normalize_gencode_biotypes=False, max_rows=100000)`。
- **输入：** GTF 文件（可为 `.gz`）；`max_rows` 限制文件的数据行数。
- **输出：** `columns`、`records`、`record_count` 和 `provenance`。坐标为 GTF 的 1-based 闭区间。
- **示例代码：**

```python
from dnakit.annotation import read_gtfparse_table

result = read_gtfparse_table("genes.gtf", features=["exon"])
print(result["record_count"], [(row["gene_name"], row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
3 [('GENE1', 100, 200), ('GENE1', 300, 500), ('GENE2', 800, 900)]
```

## 9) `dnakit.annotation.select_gtfparse_records()`

- **作用：** 选出基因 ID、基因名或转录本 ID 在给定列表中的 GTF 记录。
- **处理方法：** 读取 GTF 后，只要记录的 `gene_id`、`gene_name` 或 `transcript_id` 之一出现在对应列表中，就保留该记录。至少要给出一个列表。
- **API：** `dnakit.annotation.select_gtfparse_records(source, *, gene_ids=(), gene_names=(), transcript_ids=(), features=(), max_rows=100000)`。
- **输入：** GTF 文件；三个 ID 或名称列表；可选的 `features` 过滤。
- **输出：** `records`、`record_count`、`gene_ids` 和 `transcript_ids`。
- **示例代码：**

```python
from dnakit.annotation import select_gtfparse_records

result = select_gtfparse_records("genes.gtf", gene_names=["GENE2"])
print(result["gene_ids"], [(row["feature"], row["transcript_id"]) for row in result["records"]])
```

- **示例输出：**

```text
['ENSG2'] [('gene', None), ('transcript', 'ENST2'), ('exon', 'ENST2')]
```

## 10) `dnakit.annotation.summarize_gtfparse_annotation()`

- **作用：** 统计 GTF 的特征类型、基因数、转录本数、基因生物类型和染色体分布。
- **处理方法：** 基因数与转录本数按不同的 `gene_id`、`transcript_id` 计数。基因生物类型按每个基因计数一次。
- **API：** `dnakit.annotation.summarize_gtfparse_annotation(source, *, max_rows=100000)`。
- **输入：** GTF 文件。
- **输出：** `feature_counts`、`gene_count`、`transcript_count`、`gene_biotype_counts` 和 `chromosome_counts`。
- **示例代码：**

```python
from dnakit.annotation import summarize_gtfparse_annotation

result = summarize_gtfparse_annotation("genes.gtf")
print(result["feature_counts"], result["gene_biotype_counts"])
```

- **示例输出：**

```text
{'gene': 2, 'transcript': 2, 'exon': 3} {'protein_coding': 1, 'lncRNA': 1}
```

## 11) `dnakit.annotation.create_gffutils_database()`

- **作用：** 从 GFF3 或 GTF 文件创建 gffutils 的 SQLite 数据库。
- **处理方法：** 调用 `gffutils.create_db`，格式由 gffutils 自动识别。数据库路径必须是新文件，已存在时报错，不会覆盖旧数据库。输入文件大小受 `max_input_bytes` 限制。
- **API：** `dnakit.annotation.create_gffutils_database(source, database, *, max_input_bytes=100000000)`。
- **输入：** GFF3 或 GTF 文件；`database` 为新数据库路径，其父目录必须存在。
- **输出：** `database`、`input_format`、`featuretype_counts` 和 `feature_count`。坐标为 GFF 的 1-based 闭区间。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import create_gffutils_database

target = Path(tempfile.mkdtemp()) / "genes.db"
result = create_gffutils_database("genes.gff3", target)
print(result["input_format"], result["featuretype_counts"], result["feature_count"])
```

- **示例输出：**

```text
gff3 {'exon': 3, 'gene': 2, 'mRNA': 2} 7
```

## 12) `dnakit.annotation.list_gffutils_children()`

- **作用：** 列出已有 gffutils 数据库中某个特征的子特征。
- **处理方法：** 调用 `FeatureDB.children`。`level=1` 只返回直接子特征，`level=2` 返回孙特征。GTF 生成的数据库会把转录本记为自身的子特征，本函数会去掉这一条自引用。
- **API：** `dnakit.annotation.list_gffutils_children(database, parent_id, *, featuretype=None, level=1, max_features=1000)`。
- **输入：** gffutils 数据库路径；父特征 ID；可选的 `featuretype`；`level` 取 1 到 10。
- **输出：** `records` 列表，每项含 `id`、`featuretype`、`seqid`、`start`、`end` 和 `strand`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import create_gffutils_database, list_gffutils_children

database = Path(tempfile.mkdtemp()) / "genes.db"
create_gffutils_database("genes.gff3", database)
result = list_gffutils_children(database, "gene1", level=2)
print([(row["id"], row["featuretype"], row["start"], row["end"]) for row in result["records"]])
```

- **示例输出：**

```text
[('exon1', 'exon', 100, 200), ('exon2', 'exon', 300, 500)]
```

## 13) `dnakit.annotation.count_gffutils_featuretypes()`

- **作用：** 统计已有 gffutils 数据库中每种特征类型的数量。
- **处理方法：** 读取数据库的 `featuretypes()`，再对每种类型调用 `count_features_of_type`。
- **API：** `dnakit.annotation.count_gffutils_featuretypes(database)`。
- **输入：** gffutils 数据库路径。
- **输出：** `featuretype_counts` 字典和 `feature_count` 总数。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import count_gffutils_featuretypes, create_gffutils_database

database = Path(tempfile.mkdtemp()) / "genes.db"
create_gffutils_database("genes.gff3", database)
result = count_gffutils_featuretypes(database)
print(result["featuretype_counts"], result["feature_count"])
```

- **示例输出：**

```text
{'exon': 3, 'gene': 2, 'mRNA': 2} 7
```

## 14) `dnakit.annotation.read_pybigwig_values()`

- **作用：** 读取 bigWig 在指定区域内每个碱基的数值。
- **处理方法：** 调用 `pyBigWig` 的 `values`。没有数据的碱基返回 `None`。区域长度受 `max_bases` 限制。
- **API：** `dnakit.annotation.read_pybigwig_values(source, chromosome, start, end, *, max_bases=100000)`。
- **输入：** bigWig 文件；染色体名；0-based 半开区间 `[start, end)`，且 `end` 不超过染色体长度。
- **输出：** `values` 列表（长度等于区域长度）、`covered_bases` 和 `length`。
- **示例代码：**

```python
import pyBigWig

from dnakit.annotation import read_pybigwig_values

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = read_pybigwig_values("signal.bw", "chr1", 8, 30)
print(result["covered_bases"], result["values"][:8])
```

- **示例输出：**

```text
10 [None, None, 1.0, 1.0, 1.0, 1.0, 1.0, None]
```

## 15) `dnakit.annotation.read_pybigwig_intervals()`

- **作用：** 返回与指定区域重叠的 bigWig 存储区间。
- **处理方法：** 调用 `pyBigWig` 的 `intervals`。返回的区间不会被裁剪到查询区域。
- **API：** `dnakit.annotation.read_pybigwig_intervals(source, chromosome, start, end, *, max_bases=1000000, max_intervals=10000)`。
- **输入：** bigWig 文件；染色体名；0-based 半开区间；`max_bases` 和 `max_intervals` 为显式上限。
- **输出：** `records` 列表，每项含 `start`、`end` 和 `value`。
- **示例代码：**

```python
import pyBigWig

from dnakit.annotation import read_pybigwig_intervals

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = read_pybigwig_intervals("signal.bw", "chr1", 8, 30)
print([(row["start"], row["end"], row["value"]) for row in result["records"]])
```

- **示例输出：**

```text
[(10, 15, 1.0), (20, 25, 3.0)]
```

## 16) `dnakit.annotation.summarize_pybigwig_bed_intervals()`

- **作用：** 对 BED 文件中的每个区间计算一个 bigWig 统计量。
- **处理方法：** 对每个区间调用 `pyBigWig` 的 `stats`。`statistic` 取 `mean`、`min`、`max`、`coverage`、`std` 或 `sum`。区间内没有数据时，`value` 为 `None`。
- **API：** `dnakit.annotation.summarize_pybigwig_bed_intervals(source, bed, *, statistic="mean", max_intervals=10000)`。
- **输入：** bigWig 文件和 BED 文件；每个区间必须位于 bigWig 的染色体范围内。
- **输出：** `records` 列表，每项含 `chrom`、`start`、`end`、`name` 和 `value`。
- **示例代码：**

```python
import pyBigWig

from dnakit.annotation import summarize_pybigwig_bed_intervals

handle = pyBigWig.open("signal.bw", "w")
handle.addHeader([("chr1", 100), ("chr2", 50)])
handle.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 25], values=[1.0, 3.0])
handle.addEntries(["chr2"], [0], ends=[10], values=[2.0])
handle.close()
result = summarize_pybigwig_bed_intervals("signal.bw", "seqs.bed", statistic="sum")
print([(row["name"], row["value"]) for row in result["records"]])
```

- **示例输出：**

```text
[('p1', 5.0), ('p2', None), ('p3', 8.0)]
```

## 17) `dnakit.annotation.intersect_pyranges1_intervals()`

- **作用：** 返回源区间中与另一文件区间重叠的部分。
- **处理方法：** 调用 `pyranges1` 的 `intersect_overlaps`。输出的区间裁剪到重叠范围内，列名沿用 pyranges 的 `Chromosome`、`Start`、`End`。
- **API：** `dnakit.annotation.intersect_pyranges1_intervals(source, other, *, max_rows=100000)`。
- **输入：** 两个 BED 文件。
- **输出：** `records` 列表。
- **示例代码：**

```python
from dnakit.annotation import intersect_pyranges1_intervals

result = intersect_pyranges1_intervals("a.bed", "b.bed")
print([(row["Chromosome"], row["Start"], row["End"], row["Name"]) for row in result["records"]])
```

- **示例输出：**

```text
[('chr1', 12, 18, 'a'), ('chr1', 15, 18, 'b'), ('chr1', 45, 50, 'c')]
```

## 18) `dnakit.annotation.nearest_pyranges1_intervals()`

- **作用：** 为每个源区间找到另一文件中最近的 k 个区间。
- **处理方法：** 调用 `pyranges1` 的 `nearest_ranges`。距离列为 `Distance`，采用与 BEDTools 相同的约定，即不重叠时加 1。
- **API：** `dnakit.annotation.nearest_pyranges1_intervals(source, other, *, k=1, max_rows=100000)`。
- **输入：** 两个 BED 文件；`k` 取 1 到 100。
- **输出：** `records` 列表，每项含源区间列、`other` 的 `_b` 后缀列以及 `Distance`。
- **示例代码：**

```python
from dnakit.annotation import nearest_pyranges1_intervals

result = nearest_pyranges1_intervals("a.bed", "b.bed", k=1)
print([row["Distance"] for row in result["records"]])
```

- **示例输出：**

```text
[0, 0, 0, 92]
```

## 19) `dnakit.annotation.complement_pyranges1_intervals()`

- **作用：** 返回每条染色体上未被源区间覆盖的部分。
- **处理方法：** 调用 `pyranges1` 的 `complement_ranges`，忽略链方向，并包含从 0 开始的第一段。染色体长度由 `chrom_sizes` 给出。
- **API：** `dnakit.annotation.complement_pyranges1_intervals(source, *, chrom_sizes, max_rows=100000)`。
- **输入：** BED 文件；`chrom_sizes` 是 `{染色体名: 长度}` 字典。
- **输出：** `records` 列表，每项含 `Chromosome`、`Start` 和 `End`。
- **示例代码：**

```python
from dnakit.annotation import complement_pyranges1_intervals

result = complement_pyranges1_intervals("a.bed", chrom_sizes={"chr1": 100, "chr2": 50})
print([(row["Chromosome"], row["Start"], row["End"]) for row in result["records"]])
```

- **示例输出：**

```text
[('chr1', 0, 10), ('chr1', 30, 40), ('chr1', 50, 100), ('chr2', 0, 5), ('chr2', 9, 50)]
```

## 20) `dnakit.annotation.merge_polars_bio_intervals()`

- **作用：** 用 polars-bio 合并 BED 区间。
- **处理方法：** 以 `use_zero_based=True` 读取 BED，保持与 BED 相同的坐标。polars-bio 只合并间隔严格小于 `min_dist` 的区间，因此本函数传入 `min_distance + 1`，使其与其他合并函数一致：间隔不超过 `min_distance` 时合并。
- **API：** `dnakit.annotation.merge_polars_bio_intervals(source, *, min_distance=0, max_rows=100000)`。
- **输入：** BED 文件；`min_distance` 为非负整数。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `n_intervals`。
- **示例代码：**

```python
from dnakit.annotation import merge_polars_bio_intervals

result = merge_polars_bio_intervals("a.bed", min_distance=10)
print([(row["chrom"], row["start"], row["end"], row["n_intervals"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 10, 50, 3), ('chr2', 5, 9, 1)]
```

## 21) `dnakit.annotation.nearest_polars_bio_intervals()`

- **作用：** 为每个源区间找到另一文件中最近的 k 个区间。
- **处理方法：** 调用 `polars_bio.nearest`，`overlap=True` 时重叠区间的距离为 0。距离采用 0-based 半开坐标，即 `start_b - end_a`。
- **API：** `dnakit.annotation.nearest_polars_bio_intervals(source, other, *, k=1, max_rows=100000)`。
- **输入：** 两个 BED 文件；`k` 取 1 到 100。
- **输出：** `records` 列表，每项含两侧的 `_1`、`_2` 后缀列以及 `distance`。
- **示例代码：**

```python
from dnakit.annotation import nearest_polars_bio_intervals

result = nearest_polars_bio_intervals("a.bed", "b.bed", k=1)
print([row["distance"] for row in result["records"]])
```

- **示例输出：**

```text
[0, 0, 0, 91]
```

## 22) `dnakit.annotation.coverage_polars_bio_intervals()`

- **作用：** 统计另一文件覆盖源区间的碱基数。
- **处理方法：** 调用 `polars_bio.coverage`，两个输入都以 0-based 坐标读取。
- **API：** `dnakit.annotation.coverage_polars_bio_intervals(source, other, *, max_rows=100000)`。
- **输入：** 两个 BED 文件。
- **输出：** `records` 列表，每个源区间附带 `coverage`。
- **示例代码：**

```python
from dnakit.annotation import coverage_polars_bio_intervals

result = coverage_polars_bio_intervals("a.bed", "b.bed")
print([row["coverage"] for row in result["records"]])
```

- **示例输出：**

```text
[6, 3, 5, 0]
```

## 23) `dnakit.annotation.jaccard_gtars_regions()`

- **作用：** 计算两个 BED 区域集的碱基级 Jaccard 相似度。
- **处理方法：** 调用 `gtars` 的 `RegionSet.jaccard`。该值为两个集合共同覆盖的碱基数除以任一集合覆盖的碱基数。
- **API：** `dnakit.annotation.jaccard_gtars_regions(source, other, *, max_rows=100000)`。
- **输入：** 两个 BED 文件。
- **输出：** `jaccard`、`source_regions` 和 `other_regions`。
- **示例代码：**

```python
from dnakit.annotation import jaccard_gtars_regions

result = jaccard_gtars_regions("a.bed", "b.bed")
print(result["jaccard"])
```

- **示例输出：**

```text
0.171875
```

## 24) `dnakit.annotation.gaps_gtars_regions()`

- **作用：** 返回区域集中各染色体的缺口。
- **处理方法：** 调用 `gtars` 的 `RegionSet.gaps`。缺口包括首段、区域之间的间隔和末段；没有任何区域的染色体会给出整条染色体的缺口。
- **API：** `dnakit.annotation.gaps_gtars_regions(source, chrom_sizes, *, max_rows=100000)`。
- **输入：** BED 文件；`chrom_sizes` 是 `{染色体名: 长度}` 字典。
- **输出：** `gaps` 列表，每项含 `chrom`、`start` 和 `end`（0-based 半开）。
- **示例代码：**

```python
from dnakit.annotation import gaps_gtars_regions

result = gaps_gtars_regions("a.bed", {"chr1": 100, "chr2": 50})
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

- **示例输出：**

```text
[('chr1', 0, 10), ('chr1', 30, 40), ('chr1', 50, 100), ('chr2', 0, 5), ('chr2', 9, 50)]
```

## 25) `dnakit.annotation.tokenize_gtars_regions()`

- **作用：** 把 BED 区域映射到 BED3 通用区域集上，返回每个重叠对应的词元。
- **处理方法：** 用 `gtars.tokenizers.Tokenizer.from_bed` 读取通用区域集，再调用 `tokenize`。每个与通用区域重叠的查询区域都会产生一个词元；没有重叠的查询区域被丢弃，重复的词元保留。
- **API：** `dnakit.annotation.tokenize_gtars_regions(source, universe, *, max_rows=100000)`。
- **输入：** 查询 BED 文件；`universe` 必须是 BED3 文件。
- **输出：** `tokens`（形如 `chr1:12-18`）、`token_ids` 和 `token_count`。
- **示例代码：**

```python
from dnakit.annotation import tokenize_gtars_regions

result = tokenize_gtars_regions("a.bed", "universe.bed")
print(result["tokens"])
```

- **示例输出：**

```text
['chr1:12-18', 'chr1:12-18', 'chr1:45-60']
```

## 26) `dnakit.annotation.query_ncls_overlaps()`

- **作用：** 返回源区间与查询区间之间所有重叠的配对。
- **处理方法：** 按染色体建立 NCLS 索引，源区间的 id 为其行号，然后调用 `all_overlaps_both`。坐标为 0-based 半开，因此首尾相接的区间不算重叠。
- **API：** `dnakit.annotation.query_ncls_overlaps(source, queries, *, max_rows=100000, max_pairs=100000)`。
- **输入：** 源 BED 文件和查询 BED 文件；`max_pairs` 限制输出配对数。
- **输出：** `records` 列表，每项含 `chrom`、`query_index`、`target_index` 以及两个区间的坐标。索引为各文件中的行号，从 0 开始。
- **示例代码：**

```python
from dnakit.annotation import query_ncls_overlaps

result = query_ncls_overlaps("a.bed", "b.bed")
print([(row["query_index"], row["target_index"]) for row in result["records"]])
```

- **示例输出：**

```text
[(0, 0), (0, 1), (1, 2)]
```

## 27) `dnakit.annotation.count_ncls_overlaps()`

- **作用：** 统计源区间中与每个查询区间重叠的数量。
- **处理方法：** 与 `query_ncls_overlaps` 使用同样的索引，并按查询行统计。
- **API：** `dnakit.annotation.count_ncls_overlaps(source, queries, *, max_rows=100000)`。
- **输入：** 源 BED 文件和查询 BED 文件。
- **输出：** `records` 列表，每个查询区间附带 `overlap_count`，顺序与输入一致。
- **示例代码：**

```python
from dnakit.annotation import count_ncls_overlaps

result = count_ncls_overlaps("a.bed", "b.bed")
print([(row["name"], row["overlap_count"]) for row in result["records"]])
```

- **示例输出：**

```text
[('x', 2), ('y', 1), ('z', 0)]
```

## 28) `dnakit.annotation.find_ncls_containing_intervals()`

- **作用：** 返回完整包含每个查询区间的源区间。
- **处理方法：** 与 `query_ncls_overlaps` 使用同样的索引，但调用 `all_containments_both`。只有当源区间的起点不晚于查询起点、终点不早于查询终点时才匹配。
- **API：** `dnakit.annotation.find_ncls_containing_intervals(source, queries, *, max_rows=100000, max_pairs=100000)`。
- **输入：** 源 BED 文件和查询 BED 文件。
- **输出：** `records` 列表，每项含 `query_index` 和 `target_index`。
- **示例代码：**

```python
from dnakit.annotation import find_ncls_containing_intervals

result = find_ncls_containing_intervals("a.bed", "b.bed")
print([(row["query_index"], row["target_index"]) for row in result["records"]])
```

- **示例输出：**

```text
[(0, 0)]
```

## 29) `dnakit.annotation.reduce_genomicranges_intervals()`

- **作用：** 合并重叠或相接的区间。
- **处理方法：** 将 BED 转为 `genomicranges.GenomicRanges`（1-based 闭区间），调用 `reduce`，再转回 0-based 半开坐标。`min_gap_width=1` 时，重叠和相接的区间会合并；`0` 时只合并重叠的区间。`use_strand=False` 时忽略链方向。
- **API：** `dnakit.annotation.reduce_genomicranges_intervals(source, *, min_gap_width=1, use_strand=False, max_rows=100000)`。
- **输入：** BED 文件。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `strand`。元数据列不会保留。
- **示例代码：**

```python
from dnakit.annotation import reduce_genomicranges_intervals

result = reduce_genomicranges_intervals("a.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 10, 30), ('chr1', 40, 50), ('chr2', 5, 9)]
```

## 30) `dnakit.annotation.intersect_genomicranges_intervals()`

- **作用：** 返回两组区间共同覆盖的位置。
- **处理方法：** 调用 `GenomicRanges.intersect`，结果是两组覆盖范围的集合交集，而不是逐对求交。
- **API：** `dnakit.annotation.intersect_genomicranges_intervals(source, other, *, use_strand=False, max_rows=100000)`。
- **输入：** 两个 BED 文件。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `strand`。
- **示例代码：**

```python
from dnakit.annotation import intersect_genomicranges_intervals

result = intersect_genomicranges_intervals("a.bed", "b.bed")
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 12, 18), ('chr1', 45, 50)]
```

## 31) `dnakit.annotation.flank_genomicranges_intervals()`

- **作用：** 为每个区间生成指定宽度的上游或下游侧翼区间。
- **处理方法：** 调用 `GenomicRanges.flank`。`side="start"` 取左侧（或在 `use_strand=True` 时取 5' 端），`side="end"` 取右侧（或 3' 端）。侧翼不会被裁剪到染色体长度；起点小于 0 时报错。
- **API：** `dnakit.annotation.flank_genomicranges_intervals(source, width, *, side="start", use_strand=False, max_rows=100000)`。
- **输入：** BED 文件；`width` 为正整数；`side` 取 `"start"` 或 `"end"`。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `strand`。
- **示例代码：**

```python
from dnakit.annotation import flank_genomicranges_intervals

result = flank_genomicranges_intervals("a.bed", 5)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

- **示例输出：**

```text
[('chr1', 5, 10), ('chr1', 10, 15), ('chr1', 35, 40), ('chr2', 0, 5)]
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-bioframe]"        # bioframe
python -m pip install "dnakit[bio-pybedtools]"      # pybedtools，另需 bedtools 可执行文件
python -m pip install "dnakit[bio-gtfparse]"        # gtfparse
python -m pip install "dnakit[bio-gffutils]"        # gffutils
python -m pip install "dnakit[bio-pybigwig]"        # pyBigWig
python -m pip install "dnakit[bio-pyranges1]"       # pyranges1
python -m pip install "dnakit[bio-polarsbio]"       # polars-bio
python -m pip install "dnakit[bio-gtars]"           # gtars
python -m pip install "dnakit[bio-ncls]"            # ncls 与 numpy
python -m pip install "dnakit[bio-genomicranges]"   # genomicranges 与 pandas
```

## 验证方式

测试文件是 `others/tests/unit/test_annotation_backends.py`，共 57 项，主要检查以下内容：

- 每个公开函数至少调用一次真实的上游包；结果与手工计算的值或直接调用上游 API 的结果一致。例如 bioframe 与 pybedtools 的覆盖碱基数、gtars 的 Jaccard 值 11/64、ncls 的重叠与包含结果（随机区间，与暴力枚举对比）。
- BEDTools 的距离比 bioframe 多 1，polars-bio 合并时间隔要求与其他函数统一；这些约定在测试中分别检查。
- gffutils 的 GTF 自引用、gtfparse 的缺失属性、pyBigWig 的无数据碱基都有对应测试。
- 每个模块都有一项缺少可选依赖时的测试，要求抛出 `BackendUnavailableError` 并包含安装的 extra 名称。
- 非法的布尔值、越界的 `k` 与 `max_rows`、反向区间、缺少区间的 BED、超出 FASTA 范围的区间和非 BED3 的分词区域集都会被拒绝。

本批没有网络服务，因此没有 `DNAKIT_LIVE_TESTS` 相关的测试。

## 当前范围

本批覆盖 10 个包的 31 个函数，都是本地文件上的 DNA 或基因区间操作。以下操作没有包装：

- bioframe 的 `complement`、`subtract` 和 `cluster`，pybedtools 的 `intersect`，以及 gffutils 的 `FeatureDB` 查询：分别与已有的 `intersect_pybedtools_bed`、`query_gffutils_database` 功能重复，或超出本批的 DNA 区间范围。
- pyranges1 的 `join_overlaps` 与 `merge_overlaps`，polars-bio 的 `overlap`、`subtract` 和 `cluster`：与 `overlap_polars_bio_bed` 或 pyranges 0.x 的 `merge_pyranges_bed` 重复，或在本批之外。
- pyBigWig 的 `addEntries` 等写入接口，gtars 的 `lola`、`refget` 和 `genomic_distributions`：写入和远程功能不在本批的只读范围内，或不属于 DNA 区间操作。
- gtfparse 的 `write_gtf`，gffutils 的 `GFFWriter`：这些写出功能超出当前的只读范围。
- ncls 的 `k_overlaps_both` 与 `first_overlap_both`：与本批已有的重叠查询功能重复。
