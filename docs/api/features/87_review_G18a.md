# G18a 功能：bioframe 区间、区域、序列、文件与组装工具

本页记录 `annotation` 区域第 G18a 组的包装函数。每个函数位于 `dnakit.annotation`，可直接导入，例如 `from dnakit.annotation import validate_bioframe_bed`。
- 本组 36 个包装函数全部只读本地输入（BED、FASTA、SAM、tabix、bigBed、chrom.sizes 与内存映射），不需要网络。UCSC 下载、远程着丝粒表与 UCSC mRNA 抓取不在本组范围内。
- 示例输入以 `示例文件` 标签列在下文「示例输入文件」一节。`tab.bed.gz` 与 `.tbi` 由示例代码用 pysam 生成；bigBed 示例使用 pyBigWig 测试包自带的 `test.bigBed`。
- 坐标默认为 0 基半开区间。区域字符串按原样读取，不做坐标平移，与已有的 `parse_bioframe_region` 一致。
- `digest_bioframe_fasta` 对 bioframe.digest 的切割位置做了修正，详见该节。

## 已有包装与未包装的函数

- 本组处理了台账中 153 个 `other` 与 `cap` 行：新包装 36 个；别名（同一上游对象）62 行，标记为 `duplicate-of`；内部辅助与构造函数 44 行；非 DNA 颜色格式化 2 行；需要 pairix 绑定 3 行；同族变体 6 行（`select_indices`、`select_labels`、`select_mask`，与 `select_bioframe_region` 同一重叠判断）。
- 没有 `other` 或 `cap` 行残留。台账为 `/root/scratch_dnakit/ledger/bioframe.json`。

## 示例输入文件

- **示例文件：** `a.bed`

```text
chr1	0	100	a	5	+
chr1	50	150	b	3	+
chr1	300	400	c	1	-
chr2	10	60	d	2	+
```

BED6，四个区间，两条染色体

- **示例文件：** `view.bed`

```text
chr1	0	200	v1
chr2	0	100	v2
```

BED4 视图，名称 v1 与 v2

- **示例文件：** `labelled.bed`

```text
chr1	0	10	v1
chr2	0	5	v2
```

BED4，名称与视图名一致

- **示例文件：** `inside.bed`

```text
chr1	10	20
chr2	0	50
```

BED3，位于视图内

- **示例文件：** `tiles.bed`

```text
chr1	0	100
chr1	100	200
```

BED3，两个相接的区间

- **示例文件：** `chr1_view.bed`

```text
chr1	0	200	v
```

BED4 单个视图区域 chr1:0-200

- **示例文件：** `s.bed`

```text
chr10	5	9
chr2	1	3
chr1	7	8
chr2	0	1
```

BED3，未排序

- **示例文件：** `bad.bed`

```text
chr1	20	10	n
chr1	0	5	m
```

BED4，第一行 start 大于 end

- **示例文件：** `b.bed`

```text
chr1	80	120	x
chr2	0	60	y
```

BED4，两个区间

- **示例文件：** `runs.bed`

```text
chr1	0	100	x
chr1	100	200	x
chr1	200	300	y
```

BED4，相邻区间名称相同则成为一段

- **示例文件：** `g.fa`

```text
>chr1
ACGTNNGGCCGATTACAGAATTCGC
>chr2
GGGGCCCCAAAATTTTGCGCGCAT
```

两条染色体的 FASTA，已由示例生成 .fai 索引

- **示例文件：** `small.bed`

```text
chr1	0	10	r1
chr1	10	20	r2
chr2	0	12	r3
```

BED4，三个小区间

- **示例文件：** `x.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	11	60	8M	*	0	0	ACGTACGT	IIIIIIII	NM:i:0
r2	16	chr1	101	30	4M2D4M	*	0	0	ACGTACGT	IIIIIIII
```

SAM，两条比对记录，POS 为 SAM 的 1 基

- **示例文件：** `mrna.tsv`

```text
matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName	qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts
100	0	0	0	0	0	0	0	+	NM1	500	0	100	chr1	1000	20	120	1	100,	0,	20,
```

带表头的本地 mRNA 表，tName chr1，tStart 20，tEnd 120

- **示例文件：** `tab.bed`

```text
chr1	10	20	n
chr1	120	130	m
chr2	90	99	k
```

BED4，示例代码会压缩并建立 tabix 索引

- **示例文件：** `sizes.txt`

```text
chr1	1000
chr2	500
chrUn_x	10
```

chrom.sizes 文件

## 1) `validate_bioframe_bed`

- **作用：** 检查 BED3 至 BED6 文件是否满足 bioframe 的 bedframe 规则。
- **处理方法：** 用 `read_table` 按列数选择 schema 读取，再调用 `is_bedframe(raise_errors=False)`。不合法时 `passed` 为 false，不抛异常，也不改写文件。
- **API：** `validate_bioframe_bed(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件；行数不超过 `max_rows`。
- **输出：** `passed`（布尔值）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import validate_bioframe_bed

result = validate_bioframe_bed("a.bed", verbose=False)
print(result["passed"], result["record_count"])
```

**示例输出：**

```text
True 4
```

## 2) `check_bioframe_cataloged`

- **作用：** 检查 BED4 以上文件的名称列是否都是视图中的区域名。
- **处理方法：** 调用 `is_cataloged(df_view_col="name")`，只比较名称，不比较坐标。
- **API：** `check_bioframe_cataloged(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** BED4 以上的源文件，以及名称为区域名的视图文件（BED4 以上）。
- **输出：** `passed`、`check`（固定为 cataloged）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_cataloged

labelled = check_bioframe_cataloged("labelled.bed", "view.bed", verbose=False)
named = check_bioframe_cataloged("a.bed", "view.bed", verbose=False)
print(labelled["passed"], named["passed"])
```

**示例输出：**

```text
True False
```

## 3) `check_bioframe_contained`

- **作用：** 检查每个区间是否完全落在某个视图区域之内。
- **处理方法：** 调用 `is_contained`。区间与视图都按半开区间比较，区间端点等于视图端点仍算包含。
- **API：** `check_bioframe_contained(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `passed`、`check`（固定为 contained）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_contained

inside = check_bioframe_contained("inside.bed", "view.bed", verbose=False)
outside = check_bioframe_contained("a.bed", "view.bed", verbose=False)
print(inside["passed"], outside["passed"])
```

**示例输出：**

```text
True False
```

## 4) `check_bioframe_covering`

- **作用：** 检查源文件的区间是否无缝覆盖视图中的每个区域。
- **处理方法：** 调用 `is_covering`。只有并集没有空隙地铺满视图区域时才通过。
- **API：** `check_bioframe_covering(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `passed`、`check`（固定为 covering）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_covering

result = check_bioframe_covering("tiles.bed", "chr1_view.bed", verbose=False)
print(result["passed"], result["record_count"])
```

**示例输出：**

```text
True 2
```

## 5) `check_bioframe_overlapping`

- **作用：** 报告文件中是否存在共享碱基的两个区间。
- **处理方法：** 调用 `is_overlapping`。只共享一个端点（半开区间相接）不算重叠。
- **API：** `check_bioframe_overlapping(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件。
- **输出：** `passed`（true 表示存在重叠）、`check`、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_overlapping

print(check_bioframe_overlapping("a.bed", verbose=False)["passed"])
```

**示例输出：**

```text
True
```

## 6) `check_bioframe_sorted`

- **作用：** 判断 BED 文件是否已经是 bioframe 的排序顺序。
- **处理方法：** 调用 `is_sorted`。若排序后文件不变则 `passed` 为 true。
- **API：** `check_bioframe_sorted(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件。
- **输出：** `passed`、`check`（固定为 sorted）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_sorted

print(check_bioframe_sorted("s.bed", verbose=False)["passed"])
```

**示例输出：**

```text
False
```

## 7) `check_bioframe_tiling`

- **作用：** 检查源文件的区间是否首尾相接地铺满视图中每个区域。
- **处理方法：** 调用 `is_tiling(df_view_col=None)`。区间之间既不能重叠也不能有空隙。
- **API：** `check_bioframe_tiling(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `passed`、`check`（固定为 tiling）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_tiling

print(check_bioframe_tiling("tiles.bed", "chr1_view.bed", verbose=False)["passed"])
```

**示例输出：**

```text
True
```

## 8) `check_bioframe_viewframe`

- **作用：** 检查文件是否是合法的视图：有唯一名称且区间有效。
- **处理方法：** 调用 `is_viewframe(raise_errors=False)`。没有名称列（BED3）的文件不通过。
- **API：** `check_bioframe_viewframe(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 以上文件，视图应为 BED4 至 BED6。
- **输出：** `passed`、`check`（固定为 viewframe）、`record_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import check_bioframe_viewframe

print(check_bioframe_viewframe("view.bed", verbose=False)["passed"])
```

**示例输出：**

```text
True
```

## 9) `resolve_bioframe_region`

- **作用：** 把区域字符串解析为 chrom、start、end；只写染色体名时用 chromsizes 补全终点。
- **处理方法：** 调用 `parse_region`，数字按原样读取并去掉逗号，不做坐标平移。没有终点且未提供 chromsizes 时报错。
- **API：** `resolve_bioframe_region(region[必须], chromsizes[可选], check_bounds[可选], verbose[可选])`
- **输入：** 区域字符串，例如 `chr1:100-200` 或 `chr1`；可选的 chromsizes 映射。
- **输出：** `region`、`chrom`、`start`、`end`、`length`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import resolve_bioframe_region

span = resolve_bioframe_region("chr1:100-200", verbose=False)
whole = resolve_bioframe_region("chr1", chromsizes={"chr1": 1000}, verbose=False)
print(span["chrom"], span["start"], span["end"], whole["start"], whole["end"])
```

**示例输出：**

```text
chr1 100 200 0 1000
```

## 10) `format_bioframe_region`

- **作用：** 把 chrom、start、end 写成 UCSC 区域字符串，是 `resolve_bioframe_region` 的反向操作。
- **处理方法：** 调用 `to_ucsc_string`，坐标按原样写出，不做平移。要求 start 小于 end，且 start 不为负。
- **API：** `format_bioframe_region(chrom[必须], start[必须], end[必须], verbose[可选])`
- **输入：** 染色体名、起点和终点（整数）。
- **输出：** `region`、`chrom`、`start`、`end`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import format_bioframe_region

print(format_bioframe_region("chr1", 0, 100, verbose=False)["region"])
```

**示例输出：**

```text
chr1:0-100
```

## 11) `sanitize_bioframe_intervals`

- **作用：** 清理 BED 区间：处理 start 大于 end 的行，并可删除缺失值行。
- **处理方法：** 调用 `sanitize_bedframe`。`start_exceed_end_action` 可为 None（遇到错误即报错）、`flip`（交换端点）或 `drop`（删除该行）。
- **API：** `sanitize_bioframe_intervals(source[必须], start_exceed_end_action[可选], drop_null[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件。
- **输出：** `intervals`、`interval_count`、`start_exceed_end_action`、`drop_null`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import sanitize_bioframe_intervals

result = sanitize_bioframe_intervals("bad.bed", start_exceed_end_action="flip", verbose=False)
print([(row["start"], row["end"]) for row in result["intervals"]])
```

**示例输出：**

```text
[(10, 20), (0, 5)]
```

## 12) `sort_bioframe_intervals`

- **作用：** 按染色体名（文本顺序）与起点排序 BED 区间。
- **处理方法：** 调用 `sort_bedframe`。染色体名按文本比较，因此 chr10 排在 chr2 之前。
- **API：** `sort_bioframe_intervals(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件。
- **输出：** `intervals`、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import sort_bioframe_intervals

result = sort_bioframe_intervals("s.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('chr1', 7, 8), ('chr10', 5, 9), ('chr2', 0, 1), ('chr2', 1, 3)]
```

## 13) `setdiff_bioframe_intervals`

- **作用：** 删除与另一文件任一区间有重叠的整条区间。
- **处理方法：** 调用 `setdiff`。被删除的是整条区间，部分重叠的区间不会被裁剪（裁剪请用已有的 `subtract_bioframe_intervals`）。
- **API：** `setdiff_bioframe_intervals(source[必须], other[必须], max_rows[可选], verbose[可选])`
- **输入：** 两个本地 BED3 至 BED6 文件。
- **输出：** `intervals`、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import setdiff_bioframe_intervals

result = setdiff_bioframe_intervals("a.bed", "b.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('chr1', 300, 400)]
```

## 14) `trim_bioframe_intervals`

- **作用：** 把每个区间裁剪到视图区域之内；完全落在视图外的区间被删除。
- **处理方法：** 调用 `trim`。被删除的行数写入 `dropped_count`；bioframe 对视图外的行返回缺失值，包装会去掉这些行。
- **API：** `trim_bioframe_intervals(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `intervals`、`interval_count`、`dropped_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import trim_bioframe_intervals

result = trim_bioframe_intervals("a.bed", "view.bed", verbose=False)
print(result["interval_count"], result["dropped_count"], [(row["start"], row["end"]) for row in result["intervals"]])
```

**示例输出：**

```text
3 1 [(0, 100), (50, 150), (10, 60)]
```

## 15) `complement_bioframe_intervals`

- **作用：** 返回视图区域中没有被任何区间覆盖的部分。
- **处理方法：** 调用 `complement`。必须提供视图；没有视图时 bioframe 会使用极大的占位终点，因此包装不接受这种输入。
- **API：** `complement_bioframe_intervals(source[必须], view[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `gaps`、`gap_count`、`provenance`；每个缺口的 `view_region` 为所属视图名。

**示例代码：**

```python
from dnakit.annotation import complement_bioframe_intervals

result = complement_bioframe_intervals("b.bed", "view.bed", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["gaps"]])
```

**示例输出：**

```text
[('chr1', 0, 80), ('chr1', 120, 200), ('chr2', 60, 100)]
```

## 16) `count_bioframe_overlaps`

- **作用：** 统计另一文件中与每个区间重叠的区间数。
- **处理方法：** 调用 `count_overlaps`，结果保持源文件的行顺序，并新增 `count` 列。
- **API：** `count_bioframe_overlaps(source[必须], other[必须], max_rows[可选], verbose[可选])`
- **输入：** 两个本地 BED3 至 BED6 文件。
- **输出：** `intervals`（含 `count`）、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import count_bioframe_overlaps

result = count_bioframe_overlaps("a.bed", "b.bed", verbose=False)
print([(row["chrom"], row["start"], row["count"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('chr1', 0, 1), ('chr1', 50, 1), ('chr1', 300, 0), ('chr2', 10, 1)]
```

## 17) `select_bioframe_region`

- **作用：** 返回与一个区域字符串有重叠的区间。
- **处理方法：** 调用 `select`。只写染色体名时返回该染色体上的全部区间。
- **API：** `select_bioframe_region(source[必须], region[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件与区域字符串。
- **输出：** `region`、`intervals`、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import select_bioframe_region

result = select_bioframe_region("a.bed", "chr1:90-120", verbose=False)
print([row["name"] for row in result["intervals"]])
```

**示例输出：**

```text
['a', 'b']
```

## 18) `assign_bioframe_views`

- **作用：** 给每个区间标注包含它的视图区域名。
- **处理方法：** 调用 `assign_view`；与任何视图都不重叠的区间得到缺失的 `view_region`，`drop_unassigned` 为真时删除它们。
- **API：** `assign_bioframe_views(source[必须], view[必须], drop_unassigned[可选], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及视图文件（BED4 以上）。
- **输出：** `intervals`（含 `view_region`）、`interval_count`、`drop_unassigned`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import assign_bioframe_views

result = assign_bioframe_views("a.bed", "view.bed", verbose=False)
print([(row["name"], row["view_region"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('a', 'v1'), ('b', 'v1'), ('c', None), ('d', 'v2')]
```

## 19) `mark_bioframe_runs`

- **作用：** 给相邻且列值相同的区间编号，形成一段一段的 run。
- **处理方法：** 调用 `mark_runs`。区间不得重叠（除非 `allow_overlaps` 为真）；列只能取 name、score 或 strand，且文件必须具有该列。
- **API：** `mark_bioframe_runs(source[必须], column[必须], allow_overlaps[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED 文件，以及列名（name、score 或 strand）。
- **输出：** `column`、`intervals`（含 `run`）、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import mark_bioframe_runs

result = mark_bioframe_runs("runs.bed", "name", verbose=False)
print([(row["name"], row["run"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('x', 0), ('x', 0), ('y', 1)]
```

## 20) `merge_bioframe_runs`

- **作用：** 把相邻且列值相同的区间合并为一个 run。
- **处理方法：** 调用 `merge_runs`。规则与 `mark_bioframe_runs` 相同。
- **API：** `merge_bioframe_runs(source[必须], column[必须], allow_overlaps[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED 文件，以及列名（name、score 或 strand）。
- **输出：** `column`、`intervals`、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import merge_bioframe_runs

result = merge_bioframe_runs("runs.bed", "name", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('chr1', 0, 200), ('chr1', 200, 300)]
```

## 21) `pair_bioframe_intervals_by_distance`

- **作用：** 列出间距落在给定窗口内的所有不重复区间对。
- **处理方法：** 调用 `pair_by_distance`，间距默认按中点计算，满足 `min_sep <= 距离 < max_sep`。要求 min_sep 小于 max_sep。
- **API：** `pair_bioframe_intervals_by_distance(source[必须], min_sep[可选], max_sep[可选], min_intervening[可选], max_intervening[可选], relative_to[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 BED3 至 BED6 文件；可选的区间之间的中间距离限制。
- **输出：** `relative_to`、`pairs`（列名带 `_1` 与 `_2` 后缀）、`pair_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import pair_bioframe_intervals_by_distance

result = pair_bioframe_intervals_by_distance("a.bed", min_sep=0, max_sep=280, verbose=False)
print(result["pair_count"], [(row["name_1"], row["name_2"]) for row in result["pairs"]])
```

**示例输出：**

```text
2 [('a', 'b'), ('b', 'c')]
```

## 22) `sequence_gc_fraction`

- **作用：** 计算一条 DNA 字符串中 G 与 C 的比例。
- **处理方法：** 用 `sequence_text` 检查字符（允许 IUPAC 碱基代码），再调用 `seq_gc`。`mapped_only` 为真时 N 不计入分母；没有可计数碱基时返回 null。
- **API：** `sequence_gc_fraction(sequence[必须], mapped_only[可选], verbose[可选])`
- **输入：** DNA 字符串（大小写不敏感，长度不超过一千万）。
- **输出：** `length`、`gc_fraction`、`mapped_only`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import sequence_gc_fraction

result = sequence_gc_fraction("GGCCAATT", verbose=False)
print(result["gc_fraction"], result["length"])
```

**示例输出：**

```text
0.5 8
```

## 23) `digest_bioframe_fasta`

- **作用：** 用限制性内切酶把 FASTA 中每条序列切成片段。
- **处理方法：** 调用 `bioframe.digest`（依赖 Biopython）。bioframe 把 Biopython 的切割位置多加了 1，使内部切割点右移 2 个碱基；本包装减去这 2，使 EcoRI 在 G^AATTC 的 G 之后切开。测试中用 Biopython 的 `EcoRI.search` 直接核对。
- **API：** `digest_bioframe_fasta(source[必须], enzyme[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 FASTA 与 Biopython 的酶名，例如 EcoRI。
- **输出：** `enzyme`、`fragments`、`fragment_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import digest_bioframe_fasta

result = digest_bioframe_fasta("g.fa", "EcoRI", verbose=False)
print([(row["chrom"], row["start"], row["end"]) for row in result["fragments"]])
```

**示例输出：**

```text
[('chr1', 0, 18), ('chr1', 18, 25), ('chr2', 0, 24)]
```

## 24) `frac_gc_bioframe_intervals`

- **作用：** 给每个区间计算其 FASTA 序列的 GC 比例。
- **处理方法：** 调用 `frac_gc`，并把上游列 `GC` 改名为 `frac_gc`。`mapped_only` 为真时 N 不计入分母。区间与序列不在同一染色体时报错。
- **API：** `frac_gc_bioframe_intervals(source[必须], fasta[必须], mapped_only[可选], max_rows[可选], verbose[可选])`
- **输入：** BED3 至 BED6 文件与 FASTA 文件。
- **输出：** `mapped_only`、`intervals`（含 `frac_gc`）、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import frac_gc_bioframe_intervals

result = frac_gc_bioframe_intervals("small.bed", "g.fa", verbose=False)
print([round(row["frac_gc"], 3) for row in result["intervals"]])
```

**示例输出：**

```text
[0.75, 0.3, 0.667]
```

## 25) `frac_mapped_bioframe_intervals`

- **作用：** 给每个区间计算序列中非 N 碱基的比例。
- **处理方法：** 调用 `frac_mapped`。比例的分母是区间在序列内的碱基数；落在序列之外的区间得 0。
- **API：** `frac_mapped_bioframe_intervals(source[必须], fasta[必须], max_rows[可选], verbose[可选])`
- **输入：** BED3 至 BED6 文件与 FASTA 文件。
- **输出：** `intervals`（含 `frac_mapped`）、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import frac_mapped_bioframe_intervals

result = frac_mapped_bioframe_intervals("small.bed", "g.fa", verbose=False)
print([round(row["frac_mapped"], 3) for row in result["intervals"]])
```

**示例输出：**

```text
[0.8, 1.0, 1.0]
```

## 26) `frac_gene_coverage_bioframe`

- **作用：** 统计每个区间被 mRNA 比对覆盖的碱基数与重叠的 mRNA 数。
- **处理方法：** 调用 `frac_gene_coverage`，mRNA 表必须是本地带表头的制表符文件，且含 tName、tStart、tEnd 列。上游只返回碱基数（coverage）与重叠数（count），不返回比例。
- **API：** `frac_gene_coverage_bioframe(source[必须], mrna[必须], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件，以及本地 mRNA 表（UCSC all_mrna 格式，带表头）。
- **输出：** `intervals`（含 `coverage` 与 `count`）、`interval_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import frac_gene_coverage_bioframe

result = frac_gene_coverage_bioframe("a.bed", "mrna.tsv", verbose=False)
print([(row["name"], row["coverage"], row["count"]) for row in result["intervals"]])
```

**示例输出：**

```text
[('a', 80, 1), ('b', 70, 1), ('c', 0, 0), ('d', 0, 0)]
```

## 27) `make_bioframe_chromarms`

- **作用：** 把每条染色体在着丝粒中点处切成 p 臂与 q 臂。
- **处理方法：** 调用 `make_chromarms`。每条染色体都必须有中点，且中点位于其长度之内；缺少中点时上游会把整条染色体当作一个臂，因此包装拒绝这种输入。
- **API：** `make_bioframe_chromarms(chromsizes[必须], midpoints[必须], max_rows[可选], verbose[可选])`
- **输入：** 染色体长度映射与着丝粒中点映射（两者的染色体集合必须相同）。
- **输出：** `arms`（每项含 `name`，如 chr1_p）、`arm_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import make_bioframe_chromarms

result = make_bioframe_chromarms({"chr1": 1000, "chr2": 500}, {"chr1": 400, "chr2": 200}, verbose=False)
print([(row["name"], row["start"], row["end"]) for row in result["arms"]])
```

**示例输出：**

```text
[('chr1_p', 0, 400), ('chr1_q', 400, 1000), ('chr2_p', 0, 200), ('chr2_q', 200, 500)]
```

## 28) `list_bioframe_assemblies`

- **作用：** 列出 bioframe 随包提供的基因组组装元数据。
- **处理方法：** 调用 `assemblies_available`，读取包内 YAML 文件，不访问网络。结果按提供方与组装名排序。
- **API：** `list_bioframe_assemblies(verbose[可选])`
- **输入：** 无输入参数。
- **输出：** `assemblies`（含 provider、provider_build、organism、默认 roles 与 units、url）、`assembly_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import list_bioframe_assemblies

result = list_bioframe_assemblies(verbose=False)
hg38 = [row for row in result["assemblies"] if row["provider_build"] == "hg38"]
print(result["assembly_count"] > 10, hg38[0]["organism"], hg38[0]["provider"])
```

**示例输出：**

```text
True homo sapiens ucsc
```

## 29) `describe_bioframe_assembly`

- **作用：** 返回一个随包组装的序列名与长度，例如 ucsc.hg38。
- **处理方法：** 调用 `assembly_info`。未指定 roles 与 units 时使用组装的默认值；名称不存在或不唯一时报错。
- **API：** `describe_bioframe_assembly(name[必须], roles[可选], units[可选], max_rows[可选], verbose[可选])`
- **输入：** 组装名（可写成 provider.build），以及可选的 roles 与 units（None、`all` 或字符串列表）。
- **输出：** `name`、`organism`、`provider`、`provider_build`、`release_year`、`url`、`sequences`、`sequence_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import describe_bioframe_assembly

result = describe_bioframe_assembly("ucsc.hg38", verbose=False)
first = result["sequences"][0]
print(result["sequence_count"], first["name"], first["length"])
```

**示例输出：**

```text
25 chr1 248956422
```

## 30) `summarize_bioframe_fasta`

- **作用：** 列出 FASTA 文件中每条序列的名称与长度，不把序列读入内存。
- **处理方法：** 调用 `load_fasta`，长度来自 pysam 的索引。缺少 .fai 索引时，pysam 会在原文件旁建立该索引。
- **API：** `summarize_bioframe_fasta(source[必须], max_rows[可选], verbose[可选])`
- **输入：** 本地 FASTA 文件。
- **输出：** `sequences`（含 name 与 length）、`sequence_count`、`total_length`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import summarize_bioframe_fasta

result = summarize_bioframe_fasta("g.fa", verbose=False)
print(result["sequences"], result["total_length"])
```

**示例输出：**

```text
[{'name': 'chr1', 'length': 25}, {'name': 'chr2', 'length': 24}] 49
```

## 31) `read_bioframe_chromsizes`

- **作用：** 读取本地 chrom.sizes 或 chromInfo 文件。
- **处理方法：** 调用 `read_chromsizes`，只保留染色体名与长度。`filter_chroms` 为真时只保留 chr1 至 chr22、chrX、chrY 与 chrM；`natsort` 为真时按自然顺序排列。
- **API：** `read_bioframe_chromsizes(source[必须], filter_chroms[可选], natsort[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 chrom.sizes 文件；不接受 UCSC 数据库名（那需要网络）。
- **输出：** `chromsizes`（每项含 chrom 与 length）、`chromosome_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import read_bioframe_chromsizes

result = read_bioframe_chromsizes("sizes.txt", verbose=False)
print(result["chromsizes"])
```

**示例输出：**

```text
[{'chrom': 'chr1', 'length': 1000}, {'chrom': 'chr2', 'length': 500}]
```

## 32) `read_bioframe_tabix`

- **作用：** 读取经 bgzip 压缩并建立 tabix 索引的 BED 文件，可整体读取或按区域读取。
- **处理方法：** 调用 `read_tabix`（pysam）。没有 .tbi 或 .csi 索引时报错，并提示用 `pysam.tabix_index` 建立。按区域读取时必须同时给出 chrom。
- **API：** `read_bioframe_tabix(source[必须], chrom[可选], start[可选], end[可选], max_rows[可选], verbose[可选])`
- **输入：** bgzip 压缩的 BED 文件与其 tabix 索引。示例代码会从 tab.bed 生成 tab.bed.gz 与 .tbi。
- **输出：** `intervals`、`interval_count`、`provenance`。

**示例代码：**

```python
import pysam
from dnakit.annotation import read_bioframe_tabix

pysam.tabix_compress("tab.bed", "tab.bed.gz", force=True)
pysam.tabix_index("tab.bed.gz", preset="bed", force=True)
whole = read_bioframe_tabix("tab.bed.gz", verbose=False)
region = read_bioframe_tabix("tab.bed.gz", chrom="chr1", start=0, end=100, verbose=False)
print(whole["interval_count"], [row["name"] for row in region["intervals"]])
```

**示例输出：**

```text
3 ['n']
```

## 33) `read_bioframe_alignments`

- **作用：** 读取 SAM 或 BAM 比对记录。
- **处理方法：** 调用 `read_alignments`（pysam）。返回的 POS 为 0 基，即 SAM 列值减 1；按区域读取需要带索引的 BAM，SAM 文件会报错。
- **API：** `read_bioframe_alignments(source[必须], chrom[可选], start[可选], end[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 SAM 或 BAM 文件；可选 chrom、start、end。
- **输出：** `alignments`（含 QNAME、FLAG、RNAME、POS、MAPQ、CIGAR 等）、`alignment_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import read_bioframe_alignments

result = read_bioframe_alignments("x.sam", verbose=False)
print(result["alignment_count"], [(row["QNAME"], row["POS"]) for row in result["alignments"]])
```

**示例输出：**

```text
2 [('r1', 10), ('r2', 100)]
```

## 34) `read_bioframe_bigbed`

- **作用：** 读取本地 bigBed 文件中一条染色体的区间。
- **处理方法：** 调用 `read_bigbed`，使用 pyBigWig 引擎。染色体不存在时报错。示例使用 pyBigWig 测试包自带的 `test.bigBed`，因为二进制文件无法放入文档。
- **API：** `read_bioframe_bigbed(source[必须], chrom[必须], start[可选], end[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地 bigBed 文件、染色体名，以及可选的 start 与 end。
- **输出：** `chrom`、`intervals`（含 chrom、start、end、rest）、`interval_count`、`provenance`。

**示例代码：**

```python
import os
import pyBigWigTest
from dnakit.annotation import read_bioframe_bigbed

path = os.path.join(os.path.dirname(pyBigWigTest.__file__), "test.bigBed")
result = read_bioframe_bigbed(path, "chr1", verbose=False)
print(result["chrom"], result["interval_count"] > 0, result["intervals"][0]["start"])
```

**示例输出：**

```text
chr1 True 10009333
```

## 35) `read_bioframe_table`

- **作用：** 读取制表符分隔的 BED3 至 BED12 表格，并检查 schema 与列数一致。
- **处理方法：** 先逐行检查列数（所有行必须列数相同），再调用 `read_table`。以 # 开头的行被跳过；track 与 browser 行不会被跳过。
- **API：** `read_bioframe_table(source[必须], schema[可选], max_rows[可选], verbose[可选])`
- **输入：** 本地表格文件；schema 可省略（按列数推断 bed3 至 bed6、bed9、bed12），或显式给出。
- **输出：** `schema`、`rows`、`row_count`、`provenance`。

**示例代码：**

```python
from dnakit.annotation import read_bioframe_table

result = read_bioframe_table("a.bed", verbose=False)
print(result["schema"], result["row_count"], result["rows"][0]["name"])
```

**示例输出：**

```text
bed6 4 a
```

## 36) `write_bioframe_bed`

- **作用：** 把 BED 文件经 `to_bed` 校验后写到一个新路径。
- **处理方法：** 调用 `to_bed`。输出路径必须不存在，父目录必须已存在；`require_sorted` 为真时拒绝未排序的输入。写入失败时会删除不完整的输出。
- **API：** `write_bioframe_bed(source[必须], output[必须], schema[可选], require_sorted[可选], max_rows[可选], verbose[可选])`
- **输入：** 源 BED3 至 BED6 文件、新的输出路径、schema（infer 或 bed3 至 bed12）与 require_sorted。
- **输出：** `output`、`schema`、`record_count`、`provenance`。输出文件无表头，以制表符分隔。

**示例代码：**

```python
import pathlib
from dnakit.annotation import write_bioframe_bed

out = pathlib.Path("sorted_out.bed")
out.unlink(missing_ok=True)
result = write_bioframe_bed("a.bed", out, verbose=False)
print(result["record_count"], out.read_text(encoding="utf-8").splitlines()[0].split("\t")[:3])
```

**示例输出：**

```text
4 ['chr1', '0', '100']
```
