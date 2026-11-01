# pybedtools：区间统计、随机化检验与报告

这组函数位于 `dnakit.annotation`，在 [BEDTools](https://bedtools.readthedocs.io/) 之上封装了 [pybedtools](https://daler.github.io/pybedtools/) 的另外 14 个操作：`map`、`annotate`、全基因组覆盖、总覆盖碱基数、BEDPE 对 BED 的相交、区间间距、Fisher 精确检验、核苷酸组成、随机化重叠分布（个数或碱基数）及其显著性摘要和两两矩阵、特征轨迹图、IGV 批处理脚本和 UCSC 链接页。

要求与上一篇文档相同（`PATH` 中有 BEDTools，安装 `dnakit[bio-pybedtools]`；已用 BEDTools 2.27.1 与 pybedtools 0.12.1 测试）。`map`、`spacing` 与 `fisher` 的输入会先排序，因为 BEDTools 需要排序输入。随机化函数带随机种子：第 `i` 次打乱使用 `bedtools shuffle -seed seed+i`，`seed=0` 时与 pybedtools 的 `debug=True` 一致，结果可重复。绘图函数在 matplotlib `Figure` 上绘制（`dnakit[bio-plots]`），返回图片路径与所绘区间。以上函数均不访问网络；UCSC 链接页与 IGV 脚本仅是文本。

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

**示例文件：** `pairs.bedpe`

```text
chr1	100	150	chr1	300	350	p1	5	+	-
chr2	10	60	chr2	200	250	p2	6	+	+
```

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
```

## 1) `dnakit.annotation.map_pybedtools_values()`

- **作用：** 对每个 A 区间，汇总与之重叠的 B 区间的某一列。
- **处理方法：** 先对 A、B 排序，再调用 `BedTool.map()`（`bedtools map`）。`operation` 为 sum、min、max、absmin、absmax、mean、median、mode、antimode、stdev、sstdev、collapse、distinct、distinct_sort_num、distinct_sort_num_desc、count、count_distinct、first 或 last 之一；无重叠的 A 区间写入 `null`。
- **API：** `dnakit.annotation.map_pybedtools_values(a_path: 'str | Path', b_path: 'str | Path', *, column: 'int' = 5, operation: 'str' = 'sum', null: 'str' = '.', delimiter: 'str' = ',', min_overlap_a: 'float | None' = None, min_overlap_b: 'float | None' = None, reciprocal: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, split: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** 两个区间文件、B 中从 1 开始的 `column`、`operation` 及重叠选项（`min_overlap_a`、`min_overlap_b`、`reciprocal`、`same_strand`、`opposite_strand`、`split`）。
- **输出：** `n_records` 与 `records`：A 的各字段后接汇总值。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import map_pybedtools_values

result = map_pybedtools_values("a.bed", "b.bed", column=5, operation="sum", verbose=False)
for record in result["records"]:
    print(record)
```

- **示例输出：**

```text
['chr1', '100', '200', 'a1', '5', '+', '1']
['chr1', '500', '600', 'a2', '6', '-', '.']
['chr2', '10', '90', 'a3', '7', '+', '.']
```

## 2) `dnakit.annotation.annotate_pybedtools_intervals()`

- **作用：** 报告每个 A 区间被多个文件中特征覆盖的程度。
- **处理方法：** 调用 `BedTool.annotate()`（`bedtools annotate`）。每个特征文件对应一组列：覆盖比例（默认）、重叠个数（`counts`）或两者（`both`）。`names` 为文件命名。
- **API：** `dnakit.annotation.annotate_pybedtools_intervals(a_path: 'str | Path', feature_files: 'Sequence[str | Path]', *, names: 'Sequence[str] | None' = None, counts: 'bool' = False, both: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** A 文件、特征文件列表以及链与报告选项。
- **输出：** `n_records` 与 `records`：A 的各字段后接每个特征文件一个值。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import annotate_pybedtools_intervals

result = annotate_pybedtools_intervals("a.bed", ["b.bed"], counts=True, verbose=False)
for record in result["records"]:
    print(record)
```

- **示例输出：**

```text
['chr1', '100', '200', 'a1', '5', '+', '1']
['chr1', '500', '600', 'a2', '6', '-', '0']
['chr2', '10', '90', 'a3', '7', '+', '0']
```

## 3) `dnakit.annotation.genome_coverage_pybedtools()`

- **作用：** 计算 BED 区间或 BAM 比对对基因组的覆盖。
- **处理方法：** 调用 `BedTool.genome_coverage()`（`bedtools genomecov`）。`report` 为 `histogram`（默认）、`bedgraph`、`bedgraph_all`、`per_base` 或 `per_base_nonzero`。区间输入会先排序；BAM 输入必须按坐标排序，且不需要长度文件。
- **API：** `dnakit.annotation.genome_coverage_pybedtools(path: 'str | Path', genome_path: 'str | Path | None' = None, *, report: 'str' = 'histogram', scale: 'float' = 1.0, strand: 'str | None' = None, five_prime: 'bool' = False, three_prime: 'bool' = False, split: 'bool' = False, max_depth: 'int | None' = None, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 或 BAM 文件、BED 输入所需的染色体长度文件，以及选项（`scale`、`strand`、`five_prime`、`three_prime`、`split`、`max_depth`）。
- **输出：** `n_records`、`records` 与 `report`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import genome_coverage_pybedtools

graph = genome_coverage_pybedtools("a.bed", "genome.sizes", report="bedgraph", verbose=False)
print(graph["records"])
hist = genome_coverage_pybedtools("a.bed", "genome.sizes", verbose=False)
print(hist["records"][:3])
```

- **示例输出：**

```text
[['chr1', '100', '200', '1'], ['chr1', '500', '600', '1'], ['chr2', '10', '90', '1']]
[['chr1', '0', '1800', '2000', '0.9'], ['chr1', '1', '200', '2000', '0.1'], ['chr2', '0', '920', '1000', '0.92']]
```

## 4) `dnakit.annotation.total_coverage_pybedtools()`

- **作用：** 统计至少被一个区间覆盖的碱基数。
- **处理方法：** 调用 `BedTool.total_coverage()`，先合并区间再累加合并后的长度。
- **API：** `dnakit.annotation.total_coverage_pybedtools(path: 'str | Path')`。
- **输入：** 区间文件路径。
- **输出：** `total_bases` 与 `n_intervals`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import total_coverage_pybedtools

result = total_coverage_pybedtools("a.bed", verbose=False)
print(result["total_bases"], result["n_intervals"])
```

- **示例输出：**

```text
280 3
```

## 5) `dnakit.annotation.pair_to_bed_pybedtools()`

- **作用：** 保留与 BED 区间重叠的 BEDPE 配对。
- **处理方法：** 调用 `BedTool.pair_to_bed()`（`bedtools pairtobed`）。`pair_type` 为 either、neither、xor、both、notboth、ispan、ospan、notispan 或 notospan；`min_overlap_fraction`、`same_strand`、`opposite_strand` 进一步限定重叠。
- **API：** `dnakit.annotation.pair_to_bed_pybedtools(bedpe_path: 'str | Path', bed_path: 'str | Path', *, pair_type: 'str' = 'either', min_overlap_fraction: 'float | None' = None, same_strand: 'bool' = False, opposite_strand: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BEDPE 文件、BED 文件及上述选项。
- **输出：** `n_records` 与 `records`（BEDPE 行）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import pair_to_bed_pybedtools

hit = pair_to_bed_pybedtools("pairs.bedpe", "a.bed", verbose=False)
print(hit["n_records"], hit["records"][0][6])
none = pair_to_bed_pybedtools("pairs.bedpe", "b.bed", pair_type="neither", verbose=False)
print(none["n_records"])
```

- **示例输出：**

```text
2 p1
2
```

## 6) `dnakit.annotation.spacing_pybedtools_intervals()`

- **作用：** 报告每个区间与前一个区间之间的间隙。
- **处理方法：** 先对输入排序，再调用 `BedTool.spacing()`（`bedtools spacing`）。每条染色体的第一个区间为 `.`，重叠或相邻区间按 BEDTools 的写法给出负值或零。
- **API：** `dnakit.annotation.spacing_pybedtools_intervals(bed_path: 'str | Path', *, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件路径。
- **输出：** `n_records` 与 `records`：BED 各字段后接间隙。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import spacing_pybedtools_intervals

result = spacing_pybedtools_intervals("a.bed", verbose=False)
print([record[-1] for record in result["records"]])
```

- **示例输出：**

```text
['.', '300', '.']
```

## 7) `dnakit.annotation.fisher_pybedtools_overlap()`

- **作用：** 用 Fisher 精确检验判断两组区间的重叠是否多于随机预期。
- **处理方法：** 对 A、B 排序，调用 `BedTool.fisher()`（`bedtools fisher`）并解析其报告。可能的区间总数由染色体长度文件估计。
- **API：** `dnakit.annotation.fisher_pybedtools_overlap(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, min_overlap_a: 'float | None' = None, reciprocal: 'bool' = False, same_strand: 'bool' = False, opposite_strand: 'bool' = False, split: 'bool' = False)`。
- **输入：** 两个区间文件、染色体长度文件及重叠选项。
- **输出：** `table`（2x2 计数）、`left_tail`、`right_tail`、`two_tail`、`ratio` 与原始 `report`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import fisher_pybedtools_overlap

result = fisher_pybedtools_overlap("a.bed", "b.bed", "genome.sizes", verbose=False)
print(result["table"]["in -a"]["in -b"], result["table"]["not in -a"]["not in -b"])
print(result["two_tail"] < 0.001, result["right_tail"] < 0.001)
```

- **示例输出：**

```text
1 13
False False
```

## 8) `dnakit.annotation.nucleotide_content_pybedtools()`

- **作用：** 统计 FASTA 中各区间的核苷酸组成。
- **处理方法：** 调用 `BedTool.nucleotide_content()`（`bedtools nuc`）。每个区间得到 AT、GC 比例，碱基计数（A、C、G、T、N、其他）和序列长度；`pattern` 增加模式计数，`include_sequence` 增加序列，`strand_aware` 先对负链区间取反向互补。
- **API：** `dnakit.annotation.nucleotide_content_pybedtools(bed_path: 'str | Path', fasta_path: 'str | Path', *, strand_aware: 'bool' = False, pattern: 'str | None' = None, case_insensitive: 'bool' = False, include_sequence: 'bool' = False, output_path: 'str | Path | None' = None, max_rows: 'int' = 1000)`。
- **输入：** BED 文件、FASTA 文件及上述选项。
- **输出：** `header`（BEDTools 写出的列名）、`n_records` 与 `records`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import nucleotide_content_pybedtools

result = nucleotide_content_pybedtools("a.bed", "ref.fa", verbose=False)
print(result["header"][6:])
print(result["records"][0][6:])
```

- **示例输出：**

```text
['7_pct_at', '8_pct_gc', '9_num_A', '10_num_C', '11_num_G', '12_num_T', '13_num_N', '14_num_oth', '15_seq_len']
['0.500000', '0.500000', '25', '25', '25', '25', '0', '0', '100']
```

## 9) `dnakit.annotation.random_overlap_distribution_pybedtools()`

- **作用：** 计算随机打乱的 A 区间与固定 B 区间之间重叠的分布。
- **处理方法：** 在基因组内以 `bedtools shuffle -seed seed+i` 打乱 A 共 `iterations` 次，并与 B 求交。`metric="count"` 统计重叠的 A 区间数（`BedTool.randomintersection`）；`metric="bp"` 累加重叠碱基数（`BedTool.randomintersection_bp`）。同时返回未打乱 A 的重叠值作比较。
- **API：** `dnakit.annotation.random_overlap_distribution_pybedtools(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, iterations: 'int' = 100, metric: 'str' = 'count', seed: 'int' = 0)`。
- **输入：** 两个区间文件、染色体长度文件、`iterations`（最多 10000）、`metric` 与 `seed`。
- **输出：** `actual`、`distribution`（每次迭代一个值）、`metric`、`iterations` 与 `seed`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_distribution_pybedtools

count = random_overlap_distribution_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=5, verbose=False)
print(count["actual"], count["distribution"])
bp = random_overlap_distribution_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=5, metric="bp", verbose=False)
print(bp["actual"], bp["distribution"])
```

- **示例输出：**

```text
1 [0, 1, 2, 0, 1]
50 [0, 50, 25, 0, 84]
```

## 10) `dnakit.annotation.random_overlap_summary_pybedtools()`

- **作用：** 通过打乱 A 来评估 A 与 B 重叠的显著性。
- **处理方法：** 按上文计算随机化重叠个数，并像 `BedTool.randomstats()` 那样汇总：实际个数、随机个数的中位数、随机个数高于和低于实际值的比例、实际值除以中位数、实际值的百分位以及 2.5% 与 97.5% 分位数。
- **API：** `dnakit.annotation.random_overlap_summary_pybedtools(a_path: 'str | Path', b_path: 'str | Path', genome_path: 'str | Path', *, iterations: 'int' = 100, seed: 'int' = 0, include_distribution: 'bool' = False)`。
- **输入：** 两个区间文件、染色体长度文件、`iterations`、`seed` 与 `include_distribution`。
- **输出：** `actual`、`self`、`other`、`frac_randomized_above_actual`、`frac_randomized_below_actual`、`median_randomized`、`normalized`、`percentile`、`lower_2.5th`、`upper_97.5th`（及 `distribution`）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_summary_pybedtools

result = random_overlap_summary_pybedtools("a.bed", "b.bed", "genome.sizes", iterations=20, verbose=False)
print(result["actual"], result["self"], result["other"])
print(result["median_randomized"], result["percentile"])
print(round(result["normalized"], 3))
```

- **示例输出：**

```text
1 3 3
1.0 72.5
1.0
```

## 11) `dnakit.annotation.random_overlap_matrix_pybedtools()`

- **作用：** 对多个区间文件的每个有序文件对计算随机化重叠统计。
- **处理方法：** 对每一对 (A, B) 打乱 A `iterations` 次，并按 `random_overlap_summary_pybedtools` 汇总，等同于 `IntersectionMatrix.create_matrix()`。去掉扩展名的文件名作为矩阵标签，必须唯一。
- **API：** `dnakit.annotation.random_overlap_matrix_pybedtools(bed_paths: 'Sequence[str | Path]', genome_path: 'str | Path', *, iterations: 'int' = 100, seed: 'int' = 0)`。
- **输入：** 2 至 12 个区间文件、染色体长度文件、`iterations` 与 `seed`。
- **输出：** `names` 与 `matrix`（`name_a -> name_b -> 摘要` 的映射）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import random_overlap_matrix_pybedtools

result = random_overlap_matrix_pybedtools(["a.bed", "b.bed"], "genome.sizes", iterations=5, verbose=False)
print(result["names"])
for first in result["names"]:
    print(first, [result["matrix"][first][second]["actual"] for second in result["names"]])
```

- **示例输出：**

```text
['a', 'b']
a [3, 1]
b [1, 3]
```

## 12) `dnakit.annotation.plot_interval_tracks_pybedtools()`

- **作用：** 把多个 BED 文件在同一区域上画成堆叠的特征轨迹，并保存为 PNG。
- **处理方法：** 读取每个文件中落在 `chrom:start-end` 内的区间，并像 `TrackCollection.plot()` 一样在 matplotlib `Figure` 上画成条形。`label_intervals` 在条形上方写出 BED 名称。
- **API：** `dnakit.annotation.plot_interval_tracks_pybedtools(tracks: 'Mapping[str, str | Path]', chrom: 'str', start: 'int', end: 'int', output_path: 'str | Path', *, label_intervals: 'bool' = False)`。
- **输入：** 轨迹标签到 BED 文件的映射、区域、新的 PNG 路径与 `label_intervals`。
- **输出：** `image_path`、`region` 与 `tracks`（每条轨迹所绘区间）。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import plot_interval_tracks_pybedtools

result = plot_interval_tracks_pybedtools({"a": "a.bed", "b": "b.bed"}, "chr1", 0, 1000, "tracks.png", verbose=False)
print(open(result["image_path"], "rb").read(4))
print({label: len(items) for label, items in result["tracks"].items()})
```

- **示例输出：**

```text
b'\x89PNG'
{'a': 2, 'b': 2}
```

## 13) `dnakit.annotation.build_igv_script_pybedtools()`

- **作用：** 生成对每个区间截图的 IGV 批处理脚本。
- **处理方法：** 调用 `BedTool.igv()`（`bedtools igv`）并返回脚本文本。`slop` 为每个位点加边，`sort` 与 `collapse` 设置显示方式，`name_in_file` 用 BED 名称命名截图。不会启动 IGV。
- **API：** `dnakit.annotation.build_igv_script_pybedtools(bed_path: 'str | Path', *, snapshot_dir: 'str' = '/tmp/snapshots', image_format: 'str' = 'png', slop: 'int' = 0, sort: 'str | None' = None, collapse: 'bool' = False, name_in_file: 'bool' = False, output_path: 'str | Path | None' = None)`。
- **输入：** BED 文件和截图选项（`snapshot_dir`、`image_format`、`slop` 等）。
- **输出：** `script`、`n_commands` 与 `output_path`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import build_igv_script_pybedtools

result = build_igv_script_pybedtools("a.bed", snapshot_dir="/snaps", slop=5, verbose=False)
print(result["n_commands"])
print(result["script"].splitlines()[:3])
```

- **示例输出：**

```text
3
['snapshotDirectory /snaps', 'goto chr1:95-205', 'snapshot chr1_100_200_slop5.png']
```

## 14) `dnakit.annotation.build_ucsc_links_page_pybedtools()`

- **作用：** 生成 HTML 表格，把每个区间链接到 UCSC 基因组浏览器的对应位置。
- **处理方法：** 运行 `bedtools links`（即 `BedTool.links()`）。页面在本地生成；链接指向 `base_url`/cgi-bin/hgTracks，使用 `organism` 与 `database`。
- **API：** `dnakit.annotation.build_ucsc_links_page_pybedtools(bed_path: 'str | Path', *, base_url: 'str' = 'http://genome.ucsc.edu', organism: 'str' = 'human', database: 'str' = 'hg19', output_path: 'str | Path | None' = None)`。
- **输入：** BED 文件、`base_url`（http 或 https）、`organism` 与 `database`。
- **输出：** `html`、`n_links` 与 `output_path`。
- **示例代码：**

```python
from dnakit.annotation.pybedtools_p33_functions import build_ucsc_links_page_pybedtools

result = build_ucsc_links_page_pybedtools("a.bed", database="hg38", verbose=False)
print(result["n_links"])
print([line.strip() for line in result["html"].splitlines() if "<a href=" in line][0])
```

- **示例输出：**

```text
3
<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=chr1:100-200>chr1:101-200</a>
```
