# pybedtools：格式转换、格式检测与区间变换

这组函数位于 `dnakit.annotation`，封装了 [pybedtools](https://daler.github.io/pybedtools/) 中 DNAKit 以前没有的 15 个操作：BED12 转 BED6、BAM 转 BED、BED 或 BEDPE 转 BAM、区间表转记录、tabix 染色体名、BAM/CRAM/BGZF 检测、单个区间与文件的重叠检测、slop、shift、固定窗口、拼接并合并、按染色体边界截断、删除无效记录、由基因与外显子特征求内含子，以及用区间文件给 BAM 加标签。

pybedtools 调用 [BEDTools](https://bedtools.readthedocs.io/) 程序，因此必须安装 BEDTools（已用 BEDTools 2.27.1 与 pybedtools 0.12.1 测试；安装命令 `pip install "dnakit[bio-pybedtools]"` 与 `conda install -c bioconda bedtools`）。坐标与 BED 一致，为 0 起点、半开区间。所有函数接受文件路径，返回可序列化的字典，包含前 `max_rows` 条记录；`output_path`（必须是尚不存在的路径）接收完整结果。示例从模块 `dnakit.annotation.pybedtools_p33_functions` 导入。

pybedtools 清单中的其他行已由 DNAKit 现有函数覆盖：窗口搜索 `bedtools_window`、打乱 `bedtools_shuffle`、多文件相交 `bedtools_multiinter`、BEDPE 对 BEDPE `bedtools_pairtopair`、bedGraph 并集 `bedtools_unionbedg`、多 BAM 覆盖 `bedtools_multicov`，以及分数按长度归一化 `normalize_pybedtools_scores_to_length`。

## 示例数据

所有示例使用下列文件，它们在运行代码前写入工作目录。

**示例文件：** `a.bed`

```text
chr1	100	200	a1	5	+
chr1	500	600	a2	6	-
chr2	10	90	a3	7	+
```

**示例文件：** `b.bed`

```text
chr1	150	250	b1	1	+
chr1	900	950	b2	2	-
chr2	100	150	b3	3	+
```

**示例文件：** `genome.sizes`

```text
chr1	2000
chr2	1000
```

**示例文件：** `bed12.bed`

```text
chr1	100	400	g1	0	+	100	400	0	2	50,60	0,240
```

**示例文件：** `pairs.bedpe`

```text
chr1	100	150	chr1	300	350	p1	5	+	-
chr2	10	60	chr2	200	250	p2	6	+	+
```

**示例文件：** `genes.gff`

```text
chr1	t	gene	101	600	.	+	.	ID=g1
chr1	t	exon	101	200	.	+	.	Parent=g1
chr1	t	exon	301	400	.	+	.	Parent=g1
chr1	t	exon	501	600	.	+	.	Parent=g1
```

**示例文件：** `wide.bed`

```text
chr1	1900	2500	w1
chr3	1	5	w2
chr2	5	50	w3
```

**示例文件：** `bad.bed`

```text
chr1	10	100
chr1	100	10
chr1	5	6
```

**示例文件：** `make_bam.py`

```python
import pysam

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 1000}]}
with pysam.AlignmentFile("reads.bam", "wb", header=header) as out:
    for i, (ref, pos) in enumerate([(0, 120), (0, 520), (0, 1500), (1, 20)]):
        read = pysam.AlignedSegment()
        read.query_name = f"read{i}"
        read.query_sequence = "A" * 40
        read.flag = 0
        read.reference_id = ref
        read.reference_start = pos
        read.mapping_quality = 60
        read.cigar = [(0, 40)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 40)
        out.write(read)
```

## 1) `dnakit.annotation.convert_bed12_to_bed6_pybedtools()`

- **作用：** 把每条 BED12 记录按区块拆成多条 BED6 记录。
- **处理方法：** 调用 `BedTool.bed6()`，即运行 `bedtools bed12tobed6`。`number_names=True` 时把从 1 开始的区块编号写入得分列（`-n`）。
- **API：** `dnakit.annotation.convert_bed12_to_bed6_pybedtools(bed12_path: 'str | Path', *, number_names: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED12 文件路径，可选的 `number_names`、`output_path` 和 `max_rows`。
- **输出：** `n_records` 与 `records`，每条记录是六个 BED 字段的字符串列表。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import convert_bed12_to_bed6_pybedtools

result = convert_bed12_to_bed6_pybedtools("bed12.bed", verbose=False)
print(result["n_records"])
for record in result["records"]:
    print(record)
numbered = convert_bed12_to_bed6_pybedtools("bed12.bed", number_names=True, verbose=False)
print([r[4] for r in numbered["records"]])
```

- **示例输出：**

```text
2
['chr1', '100', '150', 'g1', '0', '+']
['chr1', '340', '400', 'g1', '0', '+']
['1', '2']
```

## 2) `dnakit.annotation.bam_to_bed_pybedtools()`

- **作用：** 把 BAM 比对转换为 BED、BED12 或拆分后的 BED。
- **处理方法：** 调用 `BedTool.bam_to_bed()`（`bedtools bamtobed`）。`bed12`、`split`、`split_deletions`、`edit_distance`、`tag`、`color`、`cigar` 对应 `-bed12`、`-split`、`-splitD`、`-ed`、`-tag`、`-color`、`-cigar`。
- **API：** `dnakit.annotation.bam_to_bed_pybedtools(bam_path: 'str | Path', *, bed12: 'bool' = False, split: 'bool' = False, split_deletions: 'bool' = False, edit_distance: 'bool' = False, tag: 'str | None' = None, color: 'str | None' = None, cigar: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BAM 文件路径及上述选项。
- **输出：** `n_records` 与 `records`（染色体、起点、终点、读段名、比对质量、链）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import bam_to_bed_pybedtools

exec(open("make_bam.py").read())
result = bam_to_bed_pybedtools("reads.bam", verbose=False)
print(result["n_records"])
for record in result["records"]:
    print(record)
```

- **示例输出：**

```text
4
['chr1', '120', '160', 'read0', '60', '+']
['chr1', '520', '560', 'read1', '60', '+']
['chr1', '1500', '1540', 'read2', '60', '+']
['chr2', '20', '60', 'read3', '60', '+']
```

## 3) `dnakit.annotation.convert_intervals_to_bam_pybedtools()`

- **作用：** 把 BED 或 BEDPE 区间写成未比对的 BAM 文件。
- **处理方法：** `input_format="bed"` 调用 `BedTool.to_bam()`（`bedtools bedtobam`）；`input_format="bedpe"` 调用 `BedTool.bedpe_to_bam()`（`bedpeToBam`）。`mapq` 设置比对质量，`bed12` 把 BED12 区块写成 CIGAR 区块（仅 BED），`uncompressed` 写出未压缩 BAM。
- **API：** `dnakit.annotation.convert_intervals_to_bam_pybedtools(intervals_path: 'str | Path', genome_path: 'str | Path', output_path: 'str | Path', *, input_format: 'str' = 'bed', mapq: 'int' = 255, bed12: 'bool' = False, uncompressed: 'bool' = False)`。
- **输入：** 区间文件、染色体长度文件、新的输出 BAM 路径及上述选项。
- **输出：** `output_path`、`n_alignments`（安装 pysam 时计数）和 `input_format`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import convert_intervals_to_bam_pybedtools

result = convert_intervals_to_bam_pybedtools("a.bed", "genome.sizes", "a.bam", verbose=False)
print(result["input_format"], result["n_alignments"])
pe = convert_intervals_to_bam_pybedtools("pairs.bedpe", "genome.sizes", "pairs.bam", input_format="bedpe", verbose=False)
print(pe["input_format"], pe["n_alignments"])
```

- **示例输出：**

```text
bed 3
bedpe 4
```

## 4) `dnakit.annotation.intervals_to_records_pybedtools()`

- **作用：** 把 BED、GFF、VCF 或 SAM 记录解析为带列名、带类型的值。
- **处理方法：** 调用 `BedTool.to_dataframe()`。列名遵循 pybedtools（BED 为 `chrom start end name score strand ...`）；`disable_auto_names=True` 时改用数字列名。
- **API：** `dnakit.annotation.intervals_to_records_pybedtools(path: 'str | Path', *, disable_auto_names: 'bool' = False, max_rows: 'int' = 1000)`。
- **输入：** 区间文件路径，可选的 `disable_auto_names` 和 `max_rows`。
- **输出：** `columns`、`dtypes`、`n_records` 与 `records`（字典）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import intervals_to_records_pybedtools

result = intervals_to_records_pybedtools("a.bed", verbose=False)
print(result["columns"])
print(result["dtypes"]["start"], result["n_records"])
print(result["records"][0])
```

- **示例输出：**

```text
['chrom', 'start', 'end', 'name', 'score', 'strand']
int64 3
{'chrom': 'chr1', 'start': 100, 'end': 200, 'name': 'a1', 'score': 5, 'strand': '+'}
```

## 5) `dnakit.annotation.list_tabix_contigs_pybedtools()`

- **作用：** 列出 bgzip 压缩并带 tabix 索引的文件中的染色体（contig）名称。
- **处理方法：** 用 `pysam.TabixFile` 打开文件（`.tbi` 或 `.csi` 索引须与文件相邻），与 `BedTool.tabix_contigs()` 相同。
- **API：** `dnakit.annotation.list_tabix_contigs_pybedtools(path: 'str | Path')`。
- **输入：** 带索引的 bgzip 文件路径。需要 pysam。
- **输出：** `contigs` 与 `n_contigs`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import list_tabix_contigs_pybedtools

import pysam

compressed = pysam.tabix_index("a.bed", preset="bed", force=True, keep_original=True)
result = list_tabix_contigs_pybedtools(compressed, verbose=False)
print(result["contigs"], result["n_contigs"])
```

- **示例输出：**

```text
['chr1', 'chr2'] 2
```

## 6) `dnakit.annotation.detect_alignment_file_format_pybedtools()`

- **作用：** 根据文件头魔数检测 BAM、CRAM、BGZF 与 gzip 文件。
- **处理方法：** 调用 `pybedtools.helpers.isBAM`、`isCRAM`、`isBGZIP` 和 `isGZIP`。BAM 本身也是 BGZF 压缩，所以 BAM 的 `is_bgzip` 为真；`format` 给出最具体的类型。
- **API：** `dnakit.annotation.detect_alignment_file_format_pybedtools(path: 'str | Path')`。
- **输入：** 文件路径。
- **输出：** `is_bam`、`is_cram`、`is_bgzip`、`is_gzip` 与 `format`（`bam`、`cram`、`bgzf`、`gzip` 或 `other`）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import detect_alignment_file_format_pybedtools

exec(open("make_bam.py").read())
for path in ("reads.bam", "a.bed"):
    result = detect_alignment_file_format_pybedtools(path, verbose=False)
    print(path, result["format"], result["is_bam"], result["is_bgzip"])
```

- **示例输出：**

```text
reads.bam bam True True
a.bed other False False
```

## 7) `dnakit.annotation.interval_has_hits_pybedtools()`

- **作用：** 检测单个查询区间是否与文件中的任一区间重叠。
- **处理方法：** 用 `chrom`、`start`、`end`、`strand` 构造 pybedtools 区间，再调用 `BedTool.any_hits()`。`same_strand` 要求链相同，`min_overlap_fraction` 要求查询区间被重叠的比例。
- **API：** `dnakit.annotation.interval_has_hits_pybedtools(intervals_path: 'str | Path', chrom: 'str', start: 'int', end: 'int', *, strand: 'str' = '.', same_strand: 'bool' = False, min_overlap_fraction: 'float' = 0.0)`。
- **输入：** 区间文件路径、查询坐标（BED 风格）及上述选项。
- **输出：** `has_hits` 与 `query`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import interval_has_hits_pybedtools

print(interval_has_hits_pybedtools("a.bed", "chr1", 150, 160, verbose=False)["has_hits"])
print(interval_has_hits_pybedtools("a.bed", "chr1", 300, 400, verbose=False)["has_hits"])
print(interval_has_hits_pybedtools("a.bed", "chr1", 150, 160, strand="-", same_strand=True, verbose=False)["has_hits"])
```

- **示例输出：**

```text
True
False
False
```

## 8) `dnakit.annotation.slop_pybedtools_intervals()`

- **作用：** 按两侧侧翼碱基数扩展区间，并在染色体末端截断。
- **处理方法：** 调用 `BedTool.slop()`（`bedtools slop`）。`both` 同时扩展两侧；`left` 与 `right`（须同时给出）分别扩展两侧，`strand_aware` 时按上游/下游。`as_fraction` 把数值当作区间长度的比例。
- **API：** `dnakit.annotation.slop_pybedtools_intervals(bed_path: 'str | Path', genome_path: 'str | Path', *, both: 'float | None' = None, left: 'float | None' = None, right: 'float | None' = None, strand_aware: 'bool' = False, as_fraction: 'bool' = False, header: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件、染色体长度文件及扩展选项。
- **输出：** `n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import slop_pybedtools_intervals

result = slop_pybedtools_intervals("a.bed", "genome.sizes", both=20, verbose=False)
for record in result["records"]:
    print(record[:3])
stranded = slop_pybedtools_intervals("a.bed", "genome.sizes", left=10, right=0, strand_aware=True, verbose=False)
print(stranded["records"][1][:3])
```

- **示例输出：**

```text
['chr1', '80', '220']
['chr1', '480', '620']
['chr2', '0', '110']
['chr1', '500', '610']
```

## 9) `dnakit.annotation.shift_pybedtools_intervals()`

- **作用：** 沿染色体平移区间，并在染色体末端截断。
- **处理方法：** 调用 `BedTool.shift()`（`bedtools shift`）。`shift` 平移所有区间；`plus_shift` 与 `minus_shift` 分别平移正链和负链区间。`as_fraction` 把平移量当作区间长度的比例。
- **API：** `dnakit.annotation.shift_pybedtools_intervals(bed_path: 'str | Path', genome_path: 'str | Path', *, shift: 'float | None' = None, plus_shift: 'float | None' = None, minus_shift: 'float | None' = None, as_fraction: 'bool' = False, header: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件、染色体长度文件及平移选项。
- **输出：** `n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import shift_pybedtools_intervals

result = shift_pybedtools_intervals("a.bed", "genome.sizes", shift=30, verbose=False)
for record in result["records"]:
    print(record[:3])
split = shift_pybedtools_intervals("a.bed", "genome.sizes", plus_shift=5, minus_shift=-5, verbose=False)
print([r[1] for r in split["records"]])
```

- **示例输出：**

```text
['chr1', '130', '230']
['chr1', '530', '630']
['chr2', '40', '120']
['105', '495', '15']
```

## 10) `dnakit.annotation.make_windows_pybedtools()`

- **作用：** 把基因组或 BED 区间切成固定大小或固定数量的窗口。
- **处理方法：** 调用 `BedTool.window_maker()`（`bedtools makewindows`）。`source_type` 为 `genome`（染色体长度文件）或 `bed`。给出 `window_size`（可加 `step_size` 做滑动窗口）或 `number_of_windows`；`id_naming` 为 `winnum` 或 `srcwinnum`。
- **API：** `dnakit.annotation.make_windows_pybedtools(source_path: 'str | Path', *, source_type: 'str' = 'genome', window_size: 'int | None' = None, step_size: 'int | None' = None, number_of_windows: 'int | None' = None, reverse_numbering: 'bool' = False, id_naming: 'str | None' = None, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** 染色体长度文件或 BED 文件及窗口选项。
- **输出：** `n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import make_windows_pybedtools

result = make_windows_pybedtools("genome.sizes", window_size=700, verbose=False)
print(result["n_records"], result["records"][:2])
bed = make_windows_pybedtools("a.bed", source_type="bed", number_of_windows=2, id_naming="winnum", verbose=False)
print(bed["records"][:2])
```

- **示例输出：**

```text
5 [['chr1', '0', '700'], ['chr1', '700', '1400']]
[['chr1', '100', '150', '1'], ['chr1', '150', '200', '2']]
```

## 11) `dnakit.annotation.cat_pybedtools_intervals()`

- **作用：** 拼接多个区间文件，默认再合并结果。
- **处理方法：** 调用 `BedTool.cat()`。`postmerge` 为真时拼接、排序并合并（`merge_distance` 为最大间隙，`same_strand` 按链分别合并），只保留前三列；为假时仅原样拼接。
- **API：** `dnakit.annotation.cat_pybedtools_intervals(paths: 'Sequence[str | Path]', *, postmerge: 'bool' = True, force_truncate: 'bool' = False, merge_distance: 'int' = 0, same_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** 两个或更多区间文件及上述选项。
- **输出：** `n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import cat_pybedtools_intervals

merged = cat_pybedtools_intervals(["a.bed", "b.bed"], verbose=False)
print(merged["records"])
plain = cat_pybedtools_intervals(["a.bed", "b.bed"], postmerge=False, verbose=False)
print(plain["n_records"])
```

- **示例输出：**

```text
[['chr1', '100', '250'], ['chr1', '500', '600'], ['chr1', '900', '950'], ['chr2', '10', '90'], ['chr2', '100', '150']]
6
```

## 12) `dnakit.annotation.truncate_pybedtools_to_chrom()`

- **作用：** 按大小文件中的染色体边界截断区间。
- **处理方法：** 以文件中的大小作为 `(0, 长度)` 边界调用 `BedTool.truncate_to_chrom()`。超出染色体末端的区间被截断，位于未知染色体上的区间被丢弃。
- **API：** `dnakit.annotation.truncate_pybedtools_to_chrom(bed_path: 'str | Path', genome_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件和染色体长度文件。
- **输出：** `n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import truncate_pybedtools_to_chrom

result = truncate_pybedtools_to_chrom("wide.bed", "genome.sizes", verbose=False)
for record in result["records"]:
    print(record)
```

- **示例输出：**

```text
['chr1', '1900', '2000', 'w1']
['chr2', '5', '50', 'w3']
```

## 13) `dnakit.annotation.remove_invalid_pybedtools_intervals()`

- **作用：** 删除坐标无效的记录。
- **处理方法：** 调用 `BedTool.remove_invalid()`，删除起点大于终点的记录、负坐标和格式错误的行。`n_removed` 统计被删除的数据行数。
- **API：** `dnakit.annotation.remove_invalid_pybedtools_intervals(bed_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件路径。
- **输出：** `n_records`、`records` 与 `n_removed`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import remove_invalid_pybedtools_intervals

result = remove_invalid_pybedtools_intervals("bad.bed", verbose=False)
print(result["n_records"], result["n_removed"])
print(result["records"])
```

- **示例输出：**

```text
2 1
[['chr1', '10', '100'], ['chr1', '5', '6']]
```

## 14) `dnakit.annotation.introns_pybedtools()`

- **作用：** 由 GFF 的 gene 与 exon 特征（或 BED12）推导内含子和 UTR 区间。
- **处理方法：** 调用 `BedTool.introns()`：每个基因与其内部同链的外显子配对；外显子之间的间隙为 `intron`，首个外显子之前和末个外显子之后的部分为 `utr5` 或 `utr3`。基因重叠时可能产生起点大于终点的记录，这些记录被丢弃并计入 `n_invalid_dropped`。
- **API：** `dnakit.annotation.introns_pybedtools(annotation_path: 'str | Path', *, gene: 'str' = 'gene', exon: 'str' = 'exon', output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** GFF 文件（或 BED12）及特征类型 `gene` 与 `exon`。
- **输出：** `n_records`、`records`（BED6，名称为基因号与类型）和 `n_invalid_dropped`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import introns_pybedtools

result = introns_pybedtools("genes.gff", verbose=False)
for record in result["records"]:
    print(record)
print(result["n_invalid_dropped"])
```

- **示例输出：**

```text
['chr1', '200', '300', 'g1', 'intron', '+']
['chr1', '400', '500', 'g1', 'intron', '+']
0
```

## 15) `dnakit.annotation.tag_bam_pybedtools()`

- **作用：** 给 BAM 比对添加标签，标签值为其重叠的区间文件（或区间名称）。
- **处理方法：** 调用 `BedTool.tag_bam()`（`bedtools tag`）。`labels` 为每个区间文件给出一个标签；`report` 选择以标签、区间 `names`、`scores` 或 `intervals`（chrom:start-end）作为标签值。带标签的比对写入新的 BAM。
- **API：** `dnakit.annotation.tag_bam_pybedtools(bam_path: 'str | Path', interval_files: 'Sequence[str | Path]', output_path: 'str | Path', *, labels: 'Sequence[str] | None' = None, tag: 'str' = 'YB', same_strand: 'bool' = False, opposite_strand: 'bool' = False, min_overlap_fraction: 'float | None' = None, report: 'str' = 'labels')`。
- **输入：** BAM 文件、区间文件列表、新的输出 BAM 路径及标签选项。
- **输出：** `output_path`、`n_alignments`、`tag` 和 `labels`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import tag_bam_pybedtools

import pysam

exec(open("make_bam.py").read())
result = tag_bam_pybedtools("reads.bam", ["a.bed"], "tagged.bam", labels=["A"], tag="XA", verbose=False)
print(result["n_alignments"], result["tag"])
with pysam.AlignmentFile("tagged.bam", "rb") as handle:
    print([(read.query_name, read.get_tag("XA") if read.has_tag("XA") else None) for read in handle])
```

- **示例输出：**

```text
4 XA
[('read0', 'A'), ('read1', 'A'), ('read2', None), ('read3', 'A')]
```
