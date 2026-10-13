# 93 G18g：pybedtools、gffutils 与 gtfparse 的区间与注释函数

本页记录 G18g 组封装的 22 个函数。它们来自 pybedtools、gffutils 与 gtfparse 的 DNA 级公开函数，全部只读取本地输入，离线运行；网络下载、绘图和 UCSC 二进制工具没有封装。示例输出由实际运行生成。封装函数分别位于 `dnakit.annotation` 的 `pybedtools_review_functions`、`gffutils_review_functions` 与 `gtfparse_review_functions` 模块中，示例直接从这些模块导入。这些模块尚未接入 `dnakit.annotation` 的 `__init__`，因此 API 行给出的是模块内的完整名称。

## 覆盖统计

| 包 | 状态 | 数量 | 处理前 other 与 cap |
| --- | --- | ---: | ---: |
| pybedtools | internal | 231 | 62 |
| pybedtools | wrapped | 13 |  |
| pybedtools | needs-binary | 7 |  |
| pybedtools | plot | 7 |  |
| pybedtools | needs-network | 5 |  |
| pybedtools | duplicate-of | 4 |  |
| pybedtools | family-variant | 4 |  |
| pybedtools | deprecated | 1 |  |
| pybedtools | not-dna | 1 |  |
| gffutils | internal | 103 | 35 |
| gffutils | wrapped | 13 |  |
| gffutils | family-variant | 10 |  |
| gffutils | duplicate-of | 5 |  |
| gffutils | already-wrapped | 1 |  |
| gffutils | cli | 1 |  |
| gtfparse | duplicate-of | 7 | 5 |
| gtfparse | wrapped | 5 |  |
| gtfparse | internal | 3 |  |
| gtfparse | family-variant | 2 |  |
| gtfparse | already-wrapped | 1 |  |

状态说明：`wrapped` 为本文档中的函数；`duplicate-of` 指向同一上游函数的已封装名称；`family-variant` 由封装函数的参数选择；`internal` 为状态、全局配置、临时文件或构造器类工具；`needs-binary` 需要 samtools，本机未安装；`plot` 为绘图或显示用途；`not-dna` 为网络可达性检查。处理前每包的 other 与 cap 行数之和为 102（pybedtools 62、gffutils 35、gtfparse 5），处理后已全部改为新状态。

## 已封装与别名的说明

- `pybedtools.cbedtools.overlap` 与 `pybedtools.overlap` 为同一上游函数，别名行指向 `overlap_pybedtools_bp`。
- `gffutils.DataIterator` 与 `gffutils.iterators.DataIterator` 为同一函数，别名行指向 `read_gffutils_records`。
- `gffutils` 的 9 个 `merge_criteria` 与 `overlap_any_threshold` 通过 `merge_gffutils_features` 的 `criteria` 参数选择，状态为 `family-variant`。
- `gtfparse` 的 `attribute_parsing`、`create_missing_features` 与 `read_gtf` 中的别名指向对应的已封装名称。

## 1) `overlap_pybedtools_bp`

- **作用：** 计算两个 0 基半开区间的有符号重叠长度，判断是否重叠。
- **处理方法：** 调用 pybedtools.overlap。该上游函数的参数顺序是两个起点在前、两个终点在后，封装已按此转换。负值是两区间之间的间隔长度，首尾相接时为 0。
- **API：** `dnakit.annotation.overlap_pybedtools_bp(start1, end1, start2, end2, *, verbose=True)`
- **输入：** start1、end1、start2、end2 为非负整数，要求 start < end。
- **输出：** overlap_bp（有符号整数）和 overlaps（布尔值）。

- **示例代码：**

```python
from dnakit.annotation import overlap_pybedtools_bp

for start1, end1, start2, end2 in [(0, 10, 5, 15), (0, 10, 10, 20), (0, 10, 20, 30)]:
    result = overlap_pybedtools_bp(start1, end1, start2, end2, verbose=False)
    print(result["overlap_bp"], result["overlaps"])
```

- **示例输出：**

```text
5 True
0 False
-10 False
```

## 2) `normalize_pybedtools_scores_to_length`

- **作用：** 把 BED 指定列的数值除以区间长度并乘以系数，例如 reads 计数转为每 kb 速率。
- **处理方法：** 对每个区间调用 featurefuncs.normalized_to_length：值 × scalar ÷ (end − start)。默认 index=4、scalar=0.001 即每 kb。
- **API：** `dnakit.annotation.normalize_pybedtools_scores_to_length(source, *, index=4, scalar=0.001, max_rows=100000, verbose=True)`
- **输入：** BED3 至 BED6 文件；index 列（0 基）必须是有限数值。
- **输出：** records 列表（chrom、start、end、name、score、strand），按坐标排序。

- **示例文件：** `bg.bedgraph`

```text
chr1	0	100	10
chr1	100	200	20
```

- **示例代码：**

```python
from dnakit.annotation import normalize_pybedtools_scores_to_length

result = normalize_pybedtools_scores_to_length("bg.bedgraph", index=3, scalar=1000, verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **示例输出：**

```text
chr1 0 100 100.0
chr1 100 200 200.0
```

## 3) `scale_pybedtools_bedgraph_values`

- **作用：** 把 bedGraph 的数值列乘以系数，用于信号轨道的归一化或单位换算。
- **处理方法：** 对每个区间调用 featurefuncs.bedgraph_scale，把第 4 列乘以 scalar 并写回字符串形式的浮点数。
- **API：** `dnakit.annotation.scale_pybedtools_bedgraph_values(source, *, scalar, max_rows=100000, verbose=True)`
- **输入：** 恰好 4 列的 bedGraph 文件，第 4 列为有限数值。
- **输出：** records 列表；数值以 name 字段返回，按坐标排序。

- **示例文件：** `bg.bedgraph`

```text
chr1	0	100	10
chr1	100	200	20
```

- **示例代码：**

```python
from dnakit.annotation import scale_pybedtools_bedgraph_values

result = scale_pybedtools_bedgraph_values("bg.bedgraph", scalar=0.5, verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **示例输出：**

```text
chr1 0 100 5.0
chr1 100 200 10.0
```

## 4) `center_pybedtools_coordinates`

- **作用：** 取每个区间中点两侧固定宽度的窗口，常用于峰中心的标准化。
- **处理方法：** 调用 featurefuncs.center。窗口长度为 2 × (width // 2) bp；短于 width 的区间保持不变；width=1 时取中点一个碱基。上游把窗口起点下限设为 1 而不是 0，封装保留该行为。
- **API：** `dnakit.annotation.center_pybedtools_coordinates(source, *, width=100, max_rows=100000, verbose=True)`
- **输入：** BED3 至 BED6 文件；width 为 1 至 1e9 的整数。
- **输出：** records 列表和 width；按坐标排序。

- **示例文件：** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **示例代码：**

```python
from dnakit.annotation import center_pybedtools_coordinates

result = center_pybedtools_coordinates("peaks.bed", width=40, verbose=False)
for record in result["records"]:
    print(record["name"], record["start"], record["end"])
```

- **示例输出：**

```text
peakA 130 170
peakB 300 340
peakC 50 60
```

## 5) `three_prime_pybedtools_coordinates`

- **作用：** 取每个区间 3' 端周围的窗口，例如 polyA 位点或 3' UTR 附近的区域。
- **处理方法：** 调用 featurefuncs.three_prime。+ 链的 3' 端为 stop，- 链为 start；upstream 计入 3' 端靠向基因内部一侧的碱基，downstream 计入另一侧。
- **API：** `dnakit.annotation.three_prime_pybedtools_coordinates(source, *, upstream=500, downstream=500, max_rows=100000, verbose=True)`
- **输入：** BED3 至 BED6 文件；upstream 与 downstream 为非负整数。
- **输出：** records 列表、upstream 和 downstream；按坐标排序。

- **示例文件：** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **示例代码：**

```python
from dnakit.annotation import three_prime_pybedtools_coordinates

result = three_prime_pybedtools_coordinates("peaks.bed", upstream=10, downstream=20, verbose=False)
for record in result["records"]:
    print(record["name"], record["start"], record["end"])
```

- **示例输出：**

```text
peakA 190 220
peakB 280 310
peakC 50 80
```

## 6) `filter_pybedtools_intervals_by_length`

- **作用：** 按区间长度筛选 BED 行，可设下限、上限或两者同时使用。
- **处理方法：** 使用 BedTool.filter 依次调用 featurefuncs.greater_than（长度 > 下限）和 featurefuncs.less_than（长度 < 上限）。两个阈值均为严格不等式。
- **API：** `dnakit.annotation.filter_pybedtools_intervals_by_length(source, *, greater_than=None, less_than=None, max_rows=100000, verbose=True)`
- **输入：** BED3 至 BED6 文件；greater_than 与 less_than 至少给出一个。
- **输出：** records 列表、greater_than、less_than；按坐标排序。

- **示例文件：** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **示例代码：**

```python
from dnakit.annotation import filter_pybedtools_intervals_by_length

result = filter_pybedtools_intervals_by_length("peaks.bed", greater_than=30, less_than=100, verbose=False)
for record in result["records"]:
    print(record["name"], record["end"] - record["start"])
```

- **示例输出：**

```text
peakB 40
```

## 7) `gff_to_bed_pybedtools`

- **作用：** 把 GFF 特征转换为 BED6，便于与其他区间工具组合使用。
- **处理方法：** 对每行调用 featurefuncs.gff2bed。name_field 可指定属性名或第几列作为名称；缺省时使用 ID，属性不存在时名称为 '.'。
- **API：** `dnakit.annotation.gff_to_bed_pybedtools(source, *, name_field=None, max_rows=100000, verbose=True)`
- **输入：** 九列 GFF 文件，start 与 end 为 1 基整数且 start ≤ end。
- **输出：** records 列表（BED6 字段），起点为 0 基。

- **示例文件：** `g.gff`

```text
chr1	src	gene	101	200	.	+	.	ID=g1;Name=GENE1
chr1	src	mRNA	101	180	5	+	.	ID=t1;Parent=g1
```

- **示例代码：**

```python
from dnakit.annotation import gff_to_bed_pybedtools

result = gff_to_bed_pybedtools("g.gff", name_field="Name", verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"], record["score"], record["strand"])
```

- **示例输出：**

```text
chr1 100 180 . 5 +
chr1 100 200 GENE1 . +
```

## 8) `bed_to_gff_pybedtools`

- **作用：** 把 BED 区间写成 GFF 记录，用于需要 GFF 的下游工具。
- **处理方法：** 对每行调用 featurefuncs.bed2gff：起点加 1 变为 1 基，名称写入 Name 属性，source 与 feature 为 '.'。
- **API：** `dnakit.annotation.bed_to_gff_pybedtools(source, *, max_rows=100000, verbose=True)`
- **输入：** BED3 至 BED6 文件。
- **输出：** records 列表（九列 GFF 字段），按坐标排序。

- **示例文件：** `peaks.bed`

```text
chr1	100	200	peakA	5	+
chr1	300	340	peakB	9	-
chr2	50	60	peakC	2	+
```

- **示例代码：**

```python
from dnakit.annotation import bed_to_gff_pybedtools

result = bed_to_gff_pybedtools("peaks.bed", verbose=False)
for record in result["records"]:
    print(record["seqid"], record["start"], record["end"], record["strand"], record["attributes"])
```

- **示例输出：**

```text
chr1 101 200 + Name="peakA";
chr1 301 340 - Name="peakB";
chr2 51 60 + Name="peakC";
```

## 9) `cis_trans_pybedtools_interactions`

- **作用：** 在染色质互作（BEDPE）中，找出每一端命中的命名区间，以及同一端（cis）和另一端（distal）的伙伴。
- **处理方法：** 调用 tag_bedpe 将 BEDPE 两端与命名的查询集合相交，再调用 cis_trans_interactions 汇总。至少需要两个查询集合：只给一个 -b 文件时，bedtools 不输出名称列。
- **API：** `dnakit.annotation.cis_trans_pybedtools_interactions(bedpe, queries, *, max_rows=100000, verbose=True)`
- **输入：** 十列 BEDPE 文件；queries 为两到八个标签到 BED4 文件的映射，标签只能含字母、数字、下划线、点和短横线。
- **输出：** records 列表（target、cis、distal 与 label 字段，缺失为 '.'），按字段排序。

- **示例文件：** `pairs.bedpe`

```text
chr1	1	10	chr1	50	90	pair1	5	+	-
chr1	2	15	chr1	200	210	pair2	1	+	+
```

- **示例文件：** `tss.bed`

```text
chr1	5	6	gene1
chr1	60	61	gene2
chr1	88	89	gene3
```

- **示例文件：** `pk.bed`

```text
chr1	3	4	peak1	50	.
```

- **示例代码：**

```python
from dnakit.annotation import cis_trans_pybedtools_interactions

result = cis_trans_pybedtools_interactions(
    "pairs.bedpe", {"tss": "tss.bed", "pk": "pk.bed"}, verbose=False
)
print(result["record_count"])
for row in result["records"]:
    if row["label"] == "pair1":
        print(row["target_label"], row["target_name"], row["cis_name"], row["distal_name"])
```

- **示例输出：**

```text
10
tss gene1 peak1 gene2
tss gene1 peak1 gene3
pk peak1 gene1 gene2
pk peak1 gene1 gene3
tss gene3 gene2 peak1
tss gene3 gene2 gene1
tss gene2 gene3 peak1
tss gene2 gene3 gene1
```

## 10) `random_jaccard_pybedtools_intervals`

- **作用：** 比较两个 BED 文件的 Jaccard 指数，并给出一次随机打乱的参照值。
- **处理方法：** 对 a 在基因组上用 bedtools shuffle 打乱一次（-seed 固定），与 b 计算 Jaccard；同时给出未打乱的观测值。这只是零模型的一次抽样，不计算 p 值。
- **API：** `dnakit.annotation.random_jaccard_pybedtools_intervals(a, b, genome, *, seed=0, max_rows=100000, verbose=True)`
- **输入：** 两个 BED3 至 BED6 文件和两列的染色体长度文件；seed 为非负整数。
- **输出：** observed 与 shuffled（intersection_bp、union_bp、jaccard、n_intersections）。

- **示例文件：** `a.bed`

```text
chr1	100	300	peakA
```

- **示例文件：** `b.bed`

```text
chr1	150	350	b1
chr1	600	700	b2
```

- **示例文件：** `genome.tsv`

```text
chr1	1000
chr2	500
```

- **示例代码：**

```python
from dnakit.annotation import random_jaccard_pybedtools_intervals

result = random_jaccard_pybedtools_intervals("a.bed", "b.bed", "genome.tsv", seed=7, verbose=False)
print(result["observed"])
print(result["shuffled"])
```

- **示例输出：**

```text
{'intersection_bp': 150, 'union_bp': 350, 'jaccard': 0.428571, 'n_intersections': 1}
{'intersection_bp': 12, 'union_bp': 488, 'jaccard': 0.0245902, 'n_intersections': 1}
```

## 11) `read_gffutils_records`

- **作用：** 逐条读取 GFF3 或 GTF 特征，保留坐标、类型和属性。
- **处理方法：** 调用 gffutils.DataIterator 按文件顺序读取；属性值排序为列表；max_rows 限制条目数量。
- **API：** `dnakit.annotation.read_gffutils_records(source, *, max_rows=100000, verbose=True)`
- **输入：** GFF3 或 GTF 文件。
- **输出：** records 列表（seqid、source、featuretype、start、end、score、strand、frame、attributes）。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
from dnakit.annotation import read_gffutils_records

result = read_gffutils_records("small.gff3", verbose=False)
print(result["record_count"])
first = result["records"][0]
print(first["featuretype"], first["start"], first["end"], first["attributes"]["Name"])
```

- **示例输出：**

```text
11
gene 101 400 ['GENE1']
```

## 12) `bin_gffutils_interval`

- **作用：** 计算 UCSC 基因组分箱编号，用于按区间快速查找特征。
- **处理方法：** 调用 gffutils.bins.bins。one=True 返回能完整包含该区间的最小 bin；one=False 返回所有与该区间重叠的 bin。
- **API：** `dnakit.annotation.bin_gffutils_interval(start, stop, *, fmt='gff', one=True, verbose=True)`
- **输入：** 起点与终点为整数；fmt='gff' 时起点为 1 基且至少为 1，fmt='bed' 时起点为 0 基。
- **输出：** bins 列表（升序）。

- **示例代码：**

```python
from dnakit.annotation import bin_gffutils_interval

print(bin_gffutils_interval(100, 200, verbose=False)["bins"])
print(bin_gffutils_interval(1, 70000, one=False, verbose=False)["bins"])
```

- **示例输出：**

```text
[4681]
[1, 9, 73, 585, 4681]
```

## 13) `asinterval_gffutils_feature`

- **作用：** 把数据库中的一个特征转换为 pybedtools 区间对象的坐标与名称。
- **处理方法：** 按 ID 取出特征（GFF 或 GTF 文件先建立内存数据库），调用 gffutils.helpers.asinterval。
- **API：** `dnakit.annotation.asinterval_gffutils_feature(database, feature_id, *, verbose=True)`
- **输入：** GFF3 或 GTF 文件、gffutils 数据库（.db），以及特征 ID。
- **输出：** chrom、start（0 基）、end、name、score、strand、featuretype 与 fields。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
from dnakit.annotation import asinterval_gffutils_feature

result = asinterval_gffutils_feature("small.gff3", "g1", verbose=False)
print(result["chrom"], result["start"], result["end"], result["name"], result["strand"])
```

- **示例输出：**

```text
chr1 100 400 g1 +
```

## 14) `features_gffutils_to_bedtool`

- **作用：** 把数据库中的特征（可限定类型）批量转换为区间，便于做区间运算。
- **处理方法：** 调用 gffutils.pybedtools_integration.to_bedtool，按 seqid 与起点排序读取特征。
- **API：** `dnakit.annotation.features_gffutils_to_bedtool(database, *, featuretype=None, max_rows=100000, verbose=True)`
- **输入：** GFF3 或 GTF 文件或 .db 数据库；featuretype 可选。
- **输出：** records 列表（chrom、start、end、name、score、strand、featuretype），按坐标排序。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
from dnakit.annotation import features_gffutils_to_bedtool

result = features_gffutils_to_bedtool("small.gff3", featuretype="gene", verbose=False)
for record in result["records"]:
    print(record["chrom"], record["start"], record["end"], record["name"])
```

- **示例输出：**

```text
chr1 100 400 g1
chr1 500 700 g2
chr2 10 40 g3
```

## 15) `infer_gffutils_attribute_dialect`

- **作用：** 根据九列中的属性文本判断 GFF3 与 GTF 的方言（分隔符、引号与键顺序）。
- **处理方法：** 调用 gffutils.helpers.infer_dialect。
- **API：** `dnakit.annotation.infer_gffutils_attribute_dialect(attributes, *, verbose=True)`
- **输入：** 属性列文本，非空，最长 65536 字符。
- **输出：** fmt（gff3 或 gtf）与 dialect 字典；集合转为排序列表。

- **示例代码：**

```python
from dnakit.annotation import infer_gffutils_attribute_dialect

print(infer_gffutils_attribute_dialect("ID=g1;Name=x", verbose=False)["fmt"])
print(infer_gffutils_attribute_dialect('gene_id "A"; transcript_id "B";', verbose=False)["dialect"]["keyval separator"])
```

- **示例输出：**

```text
gff3
```

## 16) `canonical_gffutils_transcript_sequences`

- **作用：** 为每个基因选出一个典型转录本，并取出其外显子序列。
- **处理方法：** 调用 gffutils.helpers.canonical_transcripts：选 CDS 最长的转录本，没有 CDS 时选最长的转录本；用 pyfaidx 读取序列，- 链做反向互补。上游会打印选择结果，封装屏蔽了该输出。pyfaidx 会在 FASTA 旁写入 .fai 索引。
- **API：** `dnakit.annotation.canonical_gffutils_transcript_sequences(database, fasta, *, max_rows=100000, verbose=True)`
- **输入：** GFF3 或 GTF 文件（或 .db），以及覆盖所有外显子的本地 FASTA；每个基因都必须有转录本。
- **输出：** records 列表（transcript_id、gene_ids、seqid、start、end、strand、length、sequence），按转录本 ID 排序。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
from pathlib import Path
from dnakit.annotation import canonical_gffutils_transcript_sequences

pattern = ("ACGT" * 300)[:1000]
Path("ref.fa").write_text(f">chr1\n{pattern}\n>chr2\n{'G' * 100}\n", encoding="utf-8")
result = canonical_gffutils_transcript_sequences("small.gff3", "ref.fa", verbose=False)
for record in result["records"]:
    print(record["transcript_id"], record["gene_ids"], record["length"], record["sequence"][:12])
```

- **示例输出：**

```text
t1 ['g1'] 50 GTACGTACGTAC
t2 ['g2'] 200 ACGTACGTACGT
t3 ['g3'] 30 GGGGGGGGGGGG
```

## 17) `gffutils_line_to_seqfeature`

- **作用：** 把一行 GFF 转换为 Biopython 的 SeqFeature，便于在 Biopython 生态中使用。
- **处理方法：** 调用 gffutils.biopython_integration.to_seqfeature，返回类型、标识、位置与限定符。
- **API：** `dnakit.annotation.gffutils_line_to_seqfeature(line, *, verbose=True)`
- **输入：** 九列、以制表符分隔的 GFF 行。
- **输出：** type、id、location（起点 0 基、终点开区间、链）与 qualifiers。

- **示例代码：**

```python
from dnakit.annotation import gffutils_line_to_seqfeature

line = "chr1\tens\tgene\t101\t400\t.\t+\t.\tID=g1;Name=GENE1"
result = gffutils_line_to_seqfeature(line, verbose=False)
print(result["type"], result["location"], result["qualifiers"]["Name"])
```

- **示例输出：**

```text
gene {'start': 100, 'end': 400, 'strand': 1} ['GENE1']
```

## 18) `seqfeature_to_gffutils_line`

- **作用：** 把区间描述写成一行 GFF，相当于 gffutils 与 Biopython 的反向转换。
- **处理方法：** 构建 Bio.SeqFeature，调用 gffutils.biopython_integration.from_seqfeature 得到 Feature，再输出 GFF 行。
- **API：** `dnakit.annotation.seqfeature_to_gffutils_line(seqid, feature_type, start, end, *, strand='+', feature_id=None, source='.', score='.', frame='.', attributes=None, verbose=True)`
- **输入：** seqid、feature_type、0 基起点、终点（包含）与链；可选 feature_id、source、score、frame 与属性字典。属性名不能是 GFF 列名。
- **输出：** line（一行 GFF）、seqid、featuretype、start（1 基）与 end。

- **示例代码：**

```python
from dnakit.annotation import seqfeature_to_gffutils_line

result = seqfeature_to_gffutils_line(
    "chr1", "gene", 99, 400, feature_id="g9", attributes={"Name": "Z"}, verbose=False
)
print(result["line"])
```

- **示例输出：**

```text
chr1	.	gene	100	400	.	+	.	ID=g9;Name=Z
```

## 19) `sanitize_gffutils_file`

- **作用：** 生成一个经过清理的 GFF 副本：修正坐标、给每条记录加 gid，并丢弃不在基因内的特征。
- **处理方法：** 把输入复制到新路径，再调用 gffutils.helpers.sanitize_gff_file（in_place 只改副本）。输出文件头含时间戳；源文件不被修改。
- **API：** `dnakit.annotation.sanitize_gffutils_file(source, output, *, verbose=True)`
- **输入：** GFF3 文件和一个不存在的输出路径（输出目录必须存在）。
- **输出：** path、feature_count、gene_count。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
import os
from dnakit.annotation import sanitize_gffutils_file

if os.path.exists("sanitised.gff3"):
    os.remove("sanitised.gff3")
result = sanitize_gffutils_file("small.gff3", "sanitised.gff3", verbose=False)
print(result["feature_count"], result["gene_count"])
data = [line for line in open("sanitised.gff3", encoding="utf-8") if not line.startswith("#")]
print(data[0].rstrip("\n"))
```

- **示例输出：**

```text
11 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding;gid=g1
```

## 20) `merge_gffutils_features`

- **作用：** 合并同一类型且满足条件的特征，例如把同一转录本的相邻外显子合并为一段。
- **处理方法：** 按 seqid 与起点排序后调用 FeatureDB.merge。criteria 从上游 merge_criteria 中选择，所有条件都满足才合并；阈值类条件使用 threshold。
- **API：** `dnakit.annotation.merge_gffutils_features(source, *, featuretype, criteria=('seqid', 'overlap_end_inclusive', 'strand', 'feature_type'), threshold=0, max_rows=100000, verbose=True)`
- **输入：** GFF3 或 GTF 文件（或 .db）、featuretype，以及可选的 criteria 与 threshold。
- **输出：** merged_count 与 merged 列表（id、seqid、start、end、strand、member_ids）。

- **示例文件：** `small.gff3`

```text
##gff-version 3
chr1	ens	gene	101	400	.	+	.	ID=g1;Name=GENE1;gene_biotype=protein_coding
chr1	ens	mRNA	101	400	.	+	.	ID=t1;Parent=g1;Name=T1
chr1	ens	exon	101	200	.	+	.	ID=e1;Parent=t1
chr1	ens	exon	301	400	.	+	.	ID=e2;Parent=t1
chr1	ens	CDS	151	200	.	+	0	ID=c1;Parent=t1
chr1	ens	gene	501	700	.	-	.	ID=g2;Name=GENE2
chr1	ens	mRNA	501	700	.	-	.	ID=t2;Parent=g2
chr1	ens	exon	501	700	.	-	.	ID=e3;Parent=t2
chr2	ens	gene	11	40	.	+	.	ID=g3;Name=GENE3
chr2	ens	mRNA	11	40	.	+	.	ID=t3;Parent=g3
chr2	ens	exon	11	40	.	+	.	ID=e4;Parent=t3
```

- **示例代码：**

```python
from dnakit.annotation import merge_gffutils_features

result = merge_gffutils_features(
    "small.gff3", featuretype="exon", criteria=["seqid", "strand", "feature_type"], verbose=False
)
print(result["merged_count"])
for record in result["merged"]:
    print(record["seqid"], record["start"], record["end"], record["member_ids"])
```

- **示例输出：**

```text
3
chr1 101 400 ['e1', 'e2']
chr1 501 700 []
chr2 11 40 []
```

## 21) `expand_gtfparse_attribute_strings`

- **作用：** 把 GTF 第九列的属性文本展开为每个键一列，便于做表格式分析。
- **处理方法：** 调用 gtfparse.attribute_parsing.expand_attribute_strings。每个输入行对应列表中的一项；同一行的重复键用逗号连接；缺失的键为空字符串；None 行表示没有属性。
- **API：** `dnakit.annotation.expand_gtfparse_attribute_strings(attribute_strings, *, quote_char="'", usecols=None, max_rows=100000, verbose=True)`
- **输入：** 属性文本列表（字符串或 None）、可选的引号字符与需要保留的键。
- **输出：** row_count、column_count 与 columns 字典（键排序）。

- **示例代码：**

```python
from dnakit.annotation import expand_gtfparse_attribute_strings

rows = ['gene_id "G1"; transcript_id "T1"; gene_name "A";', None, 'gene_id "G2"; tag "x"; tag "y";']
result = expand_gtfparse_attribute_strings(rows, verbose=False)
for key, values in result["columns"].items():
    print(key, values)
```

- **示例输出：**

```text
gene_id ['G1', '', 'G2']
gene_name ['A', '', '']
tag ['', '', 'x,y']
transcript_id ['T1', '', '']
```

## 22) `create_missing_gtfparse_features`

- **作用：** 当 GTF 只有外显子时，生成缺失的转录本或基因行。
- **处理方法：** 先调用 gtfparse.parse_gtf_and_expand_attributes，再调用 create_missing_features 按键（如 transcript_id）分组；起点取成员最小值，终点取最大值。若该类型已存在，则不生成任何行。
- **API：** `dnakit.annotation.create_missing_gtfparse_features(source, *, feature='transcript', key='transcript_id', extra_columns=(), max_rows=100000, verbose=True)`
- **输入：** GTF 文件、feature（transcript 或 gene）、key 列名与可选的额外列。
- **输出：** feature、key、created_count 与 records 列表。

- **示例文件：** `exons.gtf`

```text
chr1	HAVANA	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	HAVANA	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	HAVANA	exon	501	600	.	-	.	gene_id "G2"; transcript_id "T2"; gene_name "GENEB";
```

- **示例代码：**

```python
from dnakit.annotation import create_missing_gtfparse_features

result = create_missing_gtfparse_features(
    "exons.gtf", feature="transcript", key="transcript_id", extra_columns=["gene_name"], verbose=False
)
for record in result["records"]:
    print(record["feature"], record["transcript_id"], record["start"], record["end"], record["strand"])
```

- **示例输出：**

```text
transcript T1 101 400 +
transcript T2 501 600 -
```
