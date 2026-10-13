# 可选表观基因组后端：接触矩阵、配对统计、覆盖度与甲基化

这一组函数是 `dnakit.epigenome` 中基于可选依赖的适配器，分别调用 [cooler](https://github.com/open2c/cooler)、[cooltools](https://github.com/open2c/cooltools)、[pairtools](https://github.com/mirnylab/pairtools)、[MACS3](https://github.com/macs3-project/MACS)、[methylpy](https://github.com/yupenghe/methylpy) 和 [deeptools](https://github.com/deeptools/deepTools) 的 Python 接口。它们只处理显式给出的文件或数值，不比对原始测序读段，不调用参考数据库，也不运行上游的完整分析流程（例如 `cooler cload` 的 BAM 输入、`bismark` 式比对或 MACS3 的 `callpeak` 统计模型）。依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的示例使用这几个小文件。`tiny.pairs` 是 7 个接触、6 个 bin 对的 pairs 文件，位置为 1-based；`typed.pairs` 带 `pair_type` 列；`dup.pairs` 含一组 PCR 重复；`tiny.allc` 是 6 个胞嘧啶位点。这些文件只用于演示接口，数值不代表真实数据。需要二进制文件（cooler、BAM、bigWig）的示例会在代码里自行生成。

**示例文件：** `tiny.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
r1	chr1	10	chr1	150	+	-
r2	chr1	20	chr1	160	+	-
r3	chr1	15	chr1	250	+	-
r4	chr1	510	chr1	520	+	-
r5	chr1	120	chr2	30	+	-
r6	chr2	40	chr2	90	+	-
r7	chr2	200	chr2	420	-	+
```

**示例文件：** `typed.pairs`

```text
## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a	chr1	10	chr1	150	+	-	UU
b	chr1	20	chr1	300	+	+	UU
c	chr1	500	chr1	510	-	+	UR
d	chr1	600	chr2	40	+	-	UU
e	!	0	chr1	700	+	-	NU
```

**示例文件：** `dup.pairs`

```text
## pairs format v1.0
#sorted: chr1-pos1-chr2-pos2
#chromsize: chr1 1000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
a	chr1	10	chr1	150	+	-	UU
b	chr1	11	chr1	151	+	-	UU
c	chr1	14	chr1	152	+	-	UU
d	chr1	300	chr1	400	+	-	UU
e	chr1	300	chr1	400	-	+	UU
f	chr1	500	chr1	600	+	-	UU
g	!	0	chr1	700	+	-	NU
```

**示例文件：** `tiny.allc`

```text
chrL	1	+	CGN	0	10	0
chrL	2	+	CGG	1	10	0
chrL	3	-	CHH	0	4	0
chr1	5	+	CGN	2	2	1
chr1	9	+	CHG	1	2	1
chr1	12	+	CGN	3	10	1
```

## 1) `dnakit.epigenome.build_contact_matrix()`

- **作用：** 把 `.pairs` 接触列表分箱，写入一个新的 `.cool` 矩阵。
- **处理方法：** 读取第 2 到 5 列（`chrom1`、`pos1`、`chrom2`、`pos2`），位置按 1-based 解释。用 cooler 的 pairs 管线按 `binsize` 分箱并把同一对 bin 的计数求和，下三角反射到上三角。文件按块流式读取。`cool_path` 已存在时拒绝覆盖。
- **API：** `dnakit.epigenome.build_contact_matrix(pairs_path, cool_path, *, chromsizes, binsize, max_bins=100000, chunksize=1000000)`。
- **输入：** pairs 文件路径；新的 `cool_path`；`chromsizes` 为 `{染色体名: 长度}`；`binsize` 为 bp；`max_bins` 限制分箱总数。
- **输出：** 包含 `cool_path`、`binsize`、`chromosomes`、`bins`、`nonzero_pixels`、`total_contacts` 和 `provenance` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix

work = Path(tempfile.mkdtemp())
result = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)
print(result["nonzero_pixels"], result["total_contacts"], result["bins"])
```

- **示例输出：**

```text
6 7 16
```

## 2) `dnakit.epigenome.fetch_contact_matrix()`

- **作用：** 读取一个区域的稠密接触矩阵。
- **处理方法：** 调用 cooler 的 `Cooler.matrix(balance=...).fetch(region)`。`balance=True` 时使用已保存的权重，没有权重列则报错。区域跨越的 bin 数超过 `max_region_bins` 时报错，以限制输出大小。
- **API：** `dnakit.epigenome.fetch_contact_matrix(cool_path, region, *, balance=False, max_region_bins=200)`。
- **输入：** `.cool` 文件路径；UCSC 风格区域，如 `chr1` 或 `chr1:0-300`；是否使用平衡权重。
- **输出：** `matrix`（列表的列表，缺失值为 `null`）、`bins`（每行的染色体、起点、终点）、`binsize` 和 `balanced`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, fetch_contact_matrix

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
block = fetch_contact_matrix(cool, "chr1:0-300")
print(block["matrix"])
```

- **示例输出：**

```text
[[0.0, 2.0, 1.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
```

## 3) `dnakit.epigenome.balance_contact_matrix()`

- **作用：** 计算迭代校正权重，并把它们保存为 `.cool` 文件的 `weight` 列。
- **处理方法：** 调用 cooler 的 `balance_cooler`。先计算权重并检查是否有有限值；如果所有 bin 都被过滤掉，则报错且不修改文件。已有的 `weight` 列会被替换。
- **API：** `dnakit.epigenome.balance_contact_matrix(cool_path, *, cis_only=False, ignore_diags=2, min_nnz=10, min_count=0, mad_max=5, tol=1e-5, max_iters=200)`。
- **输入：** `.cool` 文件路径，以及 cooler 的过滤和收敛参数。`min_nnz` 与 `min_count` 决定哪些 bin 参与平衡。
- **输出：** `finite_weights`（有效权重数）、`bins`、`scale`、`converged`、`variance`、`parameters` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, balance_contact_matrix

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
result = balance_contact_matrix(cool, min_nnz=1)
print(result["finite_weights"], result["converged"], result["scale"])
```

- **示例输出：**

```text
6 True 1.0
```

## 4) `dnakit.epigenome.calculate_coverage()`

- **作用：** 计算每个 bin 的顺式接触和与全部接触和。
- **处理方法：** 调用 cooltools 的 `coverage`。对称矩阵的对角线像素在两侧各计一次。`cis` 只包括同一染色体的伙伴，`total` 包括所有伙伴。`balanced=True` 时使用 `weight` 列。
- **API：** `dnakit.epigenome.calculate_coverage(cool_path, *, balanced=False, ignore_diags=0, max_bins=100000)`。
- **输入：** `.cool` 文件路径；`balanced`；`ignore_diags` 为忽略的对角线数；`max_bins` 限制 bin 数。
- **输出：** `bins` 列表，每项包含 `chrom`、`start`、`end`、`cis`、`total`，以及 `balanced`、`ignore_diags` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, calculate_coverage

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
rows = calculate_coverage(cool, ignore_diags=0)["bins"]
print([(row["start"], row["cis"], row["total"]) for row in rows[:3]])
```

- **示例输出：**

```text
[(0, 3.0, 3.0), (100, 2.0, 3.0), (200, 1.0, 1.0)]
```

## 5) `dnakit.epigenome.calculate_expected_contacts()`

- **作用：** 计算同一染色体内，每个基因组距离上的平均接触频率（P(s)）。
- **处理方法：** 调用 cooltools 的 `expected_cis`，关闭平滑。每条对角线的 `sum` 是接触和，`mean` 是 `sum / valid_pixels`，`valid_pixels` 是该对角线上有效像素数。`balanced=True` 时使用平衡值。
- **API：** `dnakit.epigenome.calculate_expected_contacts(cool_path, *, balanced=False, ignore_diags=2, max_bins=100000)`。
- **输入：** `.cool` 文件路径；`balanced`；`ignore_diags`；`max_bins`。
- **输出：** `rows` 列表，每项包含 `chrom`、`distance_bins`、`distance_bp`、`valid_pixels`、`sum` 和 `mean`，以及 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, calculate_expected_contacts

work = Path(tempfile.mkdtemp())
cool = build_contact_matrix("tiny.pairs", work / "tiny.cool", chromsizes={"chr1": 1000, "chr2": 600}, binsize=100)["cool_path"]
rows = calculate_expected_contacts(cool, ignore_diags=0)["rows"]
print(rows[1])
```

- **示例输出：**

```text
{'chrom': 'chr1', 'distance_bins': 1, 'distance_bp': 100, 'valid_pixels': 9, 'sum': 2.0, 'mean': 0.2222222222222222}
```

## 6) `dnakit.epigenome.calculate_insulation_score()`

- **作用：** 计算菱形插入分数（diamond insulation score），并标出插入边界。
- **处理方法：** 调用 cooltools 的 `insulation`。输入必须是已平衡的 `.cool`（含 `weight` 列）。边界阈值使用 Li 方法（scikit-image 的 `threshold_li`）。`ignore_diags=None` 时沿用平衡时记录的对角线数。
- **API：** `dnakit.epigenome.calculate_insulation_score(cool_path, window_bp, *, ignore_diags=None, max_bins=100000)`。
- **输入：** 已平衡的 `.cool` 文件路径；`window_bp` 为菱形窗口的 bp 长度；`max_bins` 限制 bin 数。
- **输出：** `window_bp`、`boundary_count` 和 `bins` 列表。每个 bin 包含 `chrom`、`start`、`end`、`is_bad_bin`、`log2_insulation_score`、`valid_pixels`、`boundary_strength` 和 `is_boundary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import build_contact_matrix, balance_contact_matrix, calculate_insulation_score

work = Path(tempfile.mkdtemp())
lines = ["## pairs format v1.0", "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2"]
count = 0
for a in range(12):
    for b in range(a, 12):
        if b - a <= 2 and (a < 6) == (b < 6):
            for _ in range(3):
                count += 1
                lines.append(f"r{count}\tchr1\t{a * 100 + 10}\tchr1\t{b * 100 + 10}\t+\t-")
(work / "dense.pairs").write_text("\n".join(lines) + "\n", encoding="utf-8")
cool = build_contact_matrix(work / "dense.pairs", work / "dense.cool", chromsizes={"chr1": 1200}, binsize=100)["cool_path"]
balance_contact_matrix(cool, min_nnz=1)
result = calculate_insulation_score(cool, 300)
print(result["boundary_count"], [round(row["log2_insulation_score"], 3) for row in result["bins"][4:8]])
```

- **示例输出：**

```text
0 [-0.585, -1.585, -1.585, -0.585]
```

## 7) `dnakit.epigenome.read_pairs_header()`

- **作用：** 读取 `.pairs` 头部的列名、染色体长度和排序信息。
- **处理方法：** 用 pairtools 的 `get_header` 分离以 `#` 开头的头部行，再用 `extract_column_names` 和 `extract_chromsizes` 解析。没有 `#chromsize:` 行时 `chromsizes` 为空字典。
- **API：** `dnakit.epigenome.read_pairs_header(pairs_path)`。
- **输入：** `.pairs` 文件路径。
- **输出：** `header_lines`、`columns`、`sorted`（`#sorted:` 的值，没有则为 `null`）、`chromsizes` 和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import read_pairs_header

result = read_pairs_header("dup.pairs")
print(result["columns"], result["chromsizes"], result["sorted"])
```

- **示例输出：**

```text
['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type'] {'chr1': 1000} chr1-pos1-chr2-pos2
```

## 8) `dnakit.epigenome.summarize_pairs()`

- **作用：** 统计映射对、单侧映射对、顺式和反式接触，以及各 `pair_type` 的数量。
- **处理方法：** 用 pairtools 的 `PairCounter`，按 `chunksize` 行分块读取并调用 `add_pairs_from_dataframe`。文件必须有 `pair_type` 列，缺少时报错。
- **API：** `dnakit.epigenome.summarize_pairs(pairs_path, *, chunksize=100000)`。
- **输入：** 带 `pair_type` 列的 `.pairs` 文件路径；`chunksize` 为每块行数。
- **输出：** `total`、`total_mapped`、`total_single_sided_mapped`、`total_unmapped`、`cis`、`trans`、`pair_types` 和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import summarize_pairs

result = summarize_pairs("typed.pairs")
print(result["total"], result["cis"], result["trans"], result["pair_types"])
```

- **示例输出：**

```text
5 3 1 {'UU': 3, 'UR': 1, 'NU': 1}
```

## 9) `dnakit.epigenome.deduplicate_pairs()`

- **作用：** 识别 PCR 重复对并统计，可选地写出去重后的 pairs 文件。
- **处理方法：** 调用 pairtools 的 `streaming_dedup`（scipy 后端）。两端位置都在 `max_mismatch` bp 以内的对组成一簇，每簇保留一个。链方向不同的对不算重复。单侧映射对只计入统计，只有给出 `unmapped_path` 时才写出。`duplicate_fraction` 为 `total_dups / total_mapped`。
- **API：** `dnakit.epigenome.deduplicate_pairs(pairs_path, *, output_path=None, unmapped_path=None, max_mismatch=3, method="max", chunksize=10000, carryover=100)`。
- **输入：** 带链方向列的 `.pairs` 文件路径；可选的 `output_path` 与 `unmapped_path`（必须是新文件）；`max_mismatch`；`method` 取 `max` 或 `sum`。
- **输出：** `total`、`total_mapped`、`total_single_sided_mapped`、`total_nodups`、`total_dups`、`duplicate_fraction`、`output`、`unmapped_output` 和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import deduplicate_pairs

result = deduplicate_pairs("dup.pairs", max_mismatch=3)
print(result["total"], result["total_dups"], result["total_nodups"], result["duplicate_fraction"])
```

- **示例输出：**

```text
7 2 4 0.3333333333333333
```

## 10) `dnakit.epigenome.build_fragment_pileup()`

- **作用：** 把一条染色体上的片段累加为分段常数的覆盖度。
- **处理方法：** 调用 MACS3 的 `Signal.Pileup.quick_pileup`。输出的连续段从位置 0 开始，到最大片段终点结束。`baseline_value` 是每段的最小值，`scale_factor` 乘以覆盖度。起点和终点会分别排序，要求每个片段的终点大于起点。
- **API：** `dnakit.epigenome.build_fragment_pileup(starts, ends, *, scale_factor=1.0, baseline_value=0.0, max_fragments=1000000)`。
- **输入：** 等长的片段起点列表和终点列表（0-based 半开坐标）；缩放因子；基线值；最多片段数。
- **输出：** `fragments`、`segments`（每段包含 `start`、`end`、`value`）和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import build_fragment_pileup

result = build_fragment_pileup([10, 50, 100], [40, 90, 140])
print([(row["start"], row["end"], row["value"]) for row in result["segments"]])
```

- **示例输出：**

```text
[(0, 10, 0.0), (10, 40, 1.0), (40, 50, 0.0), (50, 90, 1.0), (90, 100, 0.0), (100, 140, 1.0)]
```

## 11) `dnakit.epigenome.build_extended_tag_pileup()`

- **作用：** 把单端 tag 位置向两侧各延伸固定长度后，计算覆盖度。
- **处理方法：** 调用 MACS3 的 `Signal.Pileup.naive_quick_pileup`。每个 tag 覆盖 `[position - extension, position + extension)`，小于 0 的部分截断到 0。位置会排序后再计算。
- **API：** `dnakit.epigenome.build_extended_tag_pileup(positions, extension, *, max_tags=1000000)`。
- **输入：** tag 位置列表（0-based）；延伸长度；最多 tag 数。
- **输出：** `tags`、`extension`、`segments`（分段覆盖度）和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import build_extended_tag_pileup

result = build_extended_tag_pileup([100, 200, 205], 50)
print([(row["start"], row["end"], row["value"]) for row in result["segments"]])
```

- **示例输出：**

```text
[(0, 50, 0.0), (50, 155, 1.0), (155, 250, 2.0), (250, 255, 1.0)]
```

## 12) `dnakit.epigenome.call_fragment_peaks()`

- **作用：** 在片段覆盖度上找出超过阈值的峰顶（summit）。
- **处理方法：** 先用 `quick_pileup` 计算覆盖度，再调用 MACS3 的 `naive_call_peaks`：覆盖度不低于 `min_value` 的相邻段，若间隔不超过 `max_gap` 则合并；合并后长度不少于 `min_length` 的区域保留，峰顶取最高值的位置。这是朴素规则，不是 MACS3 的统计峰检出模型，也不计算 p 值。
- **API：** `dnakit.epigenome.call_fragment_peaks(starts, ends, *, min_value, max_gap=50, min_length=200, max_fragments=1000000)`。
- **输入：** 片段起点和终点列表；`min_value`；`max_gap`；`min_length`；最多片段数。
- **输出：** `fragments`、`min_value`、`max_gap`、`min_length`、`peaks`（每项包含 `summit` 和 `height`）和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import call_fragment_peaks

starts = [100] * 5 + [400] * 2 + [1000] * 4
ends = [300] * 5 + [500] * 2 + [1200] * 4
print(call_fragment_peaks(starts, ends, min_value=1.0)["peaks"])
```

- **示例输出：**

```text
[{'summit': 200, 'height': 5.0}, {'summit': 1100, 'height': 4.0}]
```

## 13) `dnakit.epigenome.filter_allc_sites()`

- **作用：** 按染色体、胞嘧啶上下文和覆盖度筛选 allc 记录。
- **处理方法：** 调用 methylpy 的 `utilities.filter_allc_files`。`mc_type` 如 `CGN`、`CHH` 会展开为具体的三核苷酸上下文。覆盖度是 allc 第 6 列的总读数。
- **API：** `dnakit.epigenome.filter_allc_sites(allc_path, output_path, *, chromosomes=None, mc_type=None, min_coverage=0, max_coverage=None)`。
- **输入：** allc 文件路径（可以是纯文本）；新的输出路径；可选的染色体列表；`mc_type`；覆盖度范围。
- **输出：** `input_sites`、`kept_sites`、`output`、`chromosomes`、`mc_type`、`min_coverage`、`max_coverage` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import filter_allc_sites

work = Path(tempfile.mkdtemp())
result = filter_allc_sites("tiny.allc", work / "filtered.allc", chromosomes=["chrL", "chr1"], mc_type="CGN", min_coverage=2)
print(result["input_sites"], result["kept_sites"])
```

- **示例输出：**

```text
6 4
```

## 14) `dnakit.epigenome.calculate_non_conversion_rate()`

- **作用：** 用非甲基化对照区域（例如 lambda DNA）估计亚硫酸氢盐非转化率。
- **处理方法：** 调用 methylpy 的 `calculate_non_conversion_rate`：对照区域内甲基化读数之和除以全部胞嘧啶读数之和，所有上下文都计入。methylpy 会在输入文件旁写索引，因此这里先把文件链接到临时目录再计算，输入目录不会出现 `.idx` 文件。染色体名按 allc 文件中的写法精确匹配。
- **API：** `dnakit.epigenome.calculate_non_conversion_rate(allc_path, control_region)`。
- **输入：** allc 文件路径；对照区域，形如 `chrom` 或 `chrom:start-end`（1-based 闭区间）。
- **输出：** `control_region`、`non_conversion_rate` 和 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import calculate_non_conversion_rate

print(calculate_non_conversion_rate("tiny.allc", "chrL")["non_conversion_rate"])
```

- **示例输出：**

```text
0.041666666666666664
```

## 15) `dnakit.epigenome.calculate_site_binomial_pvalues()`

- **作用：** 对每个胞嘧啶做单侧二项检验，判断甲基化是否高于非转化率。
- **处理方法：** 调用 methylpy 的 `allc_run_binom_tests`。p 值为 `P(X >= 甲基化读数)`，其中 `X ~ Binomial(覆盖度, 非转化率)`。没有甲基化读数或覆盖度低于 `min_coverage` 的位点记为 p = 2.0，不参与检验。输出文件按 p 值排序（使用 unix `sort`）。
- **API：** `dnakit.epigenome.calculate_site_binomial_pvalues(allc_path, output_path, *, non_conversion, min_coverage=1, sort_buffer="64M", max_report=20)`。
- **输入：** allc 文件路径；新的输出路径（只允许字母、数字和 `._/@+=,:-`）；`non_conversion`（0 到 1）；`min_coverage`；`sort_buffer`（如 `64M`）；`max_report` 为返回的最小 p 值位点数。
- **输出：** `sites_written`、`tested_sites`、`mc_class_counts`、`lowest_p_values`（每项包含 `chrom`、`pos`、`strand`、`mc_class`、`methylated`、`coverage` 和 `p_value`）、`non_conversion_rate` 和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
from dnakit.epigenome import calculate_site_binomial_pvalues

work = Path(tempfile.mkdtemp())
result = calculate_site_binomial_pvalues("tiny.allc", work / "pvalues.tsv", non_conversion=0.01, min_coverage=2, max_report=2)
print(result["tested_sites"], result["lowest_p_values"][0])
```

- **示例输出：**

```text
4 {'chrom': 'chr1', 'pos': 5, 'strand': '+', 'mc_class': 'CGN', 'methylated': 2, 'coverage': 2, 'p_value': 0.0001}
```

## 16) `dnakit.epigenome.summarize_bigwig_bins()`

- **作用：** 计算一个或多个 bigWig 文件在等宽基因组窗口内的平均信号。
- **处理方法：** 调用 deeptools 的 `getScorePerBigWigBin.getScorePerBin`。只使用所有文件中名称和长度都相同的染色体。窗口值是该窗口内有数据碱基的平均值，没有数据的窗口为 `null`。最后一个窗口可能较短。
- **API：** `dnakit.epigenome.summarize_bigwig_bins(bigwig_files, *, bin_size=1000, max_bins=100000, number_of_processors=1)`。
- **输入：** 一个或多个 bigWig 路径；`bin_size`（bp）；`max_bins` 限制窗口总数；`number_of_processors`。
- **输出：** `bin_size`、`files`、`chromosomes`、`bins`（每项包含 `chrom`、`start`、`end` 和每个文件的 `values`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
import pyBigWig
from dnakit.epigenome import summarize_bigwig_bins

work = Path(tempfile.mkdtemp())
bw = pyBigWig.open(str(work / "a.bw"), "w")
bw.addHeader([("chr1", 1000)])
bw.addEntries(["chr1", "chr1"], [0, 500], ends=[500, 1000], values=[2.0, 4.0])
bw.close()
result = summarize_bigwig_bins([work / "a.bw"], bin_size=200)
print([row["values"] for row in result["bins"]])
```

- **示例输出：**

```text
[[2.0], [2.0], [3.0], [4.0], [4.0]]
```

## 17) `dnakit.epigenome.count_bam_reads_per_bin()`

- **作用：** 统计与每个等宽窗口重叠的比对记录数。
- **处理方法：** 调用 deeptools 的 `countReadsPerBin.CountReadsPerBin`，不延伸读段。每条比对按自身的位置单独计数，一条记录可以计入它跨越的多个窗口。`min_mapping_quality` 之下的记录被丢弃。输入 BAM 必须有 `.bai` 索引。
- **API：** `dnakit.epigenome.count_bam_reads_per_bin(bam_files, *, bin_size=1000, max_bins=100000, min_mapping_quality=0, number_of_processors=1)`。
- **输入：** 一个或多个已索引的 BAM 路径；`bin_size`；`max_bins`；`min_mapping_quality`（0 到 255）；`number_of_processors`。
- **输出：** `bin_size`、`files`、`chromosomes`、`bins`（每项包含 `chrom`、`start`、`end` 和每个文件的 `counts`）和 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
import pysam
from dnakit.epigenome import count_bam_reads_per_bin

work = Path(tempfile.mkdtemp())
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
bam = work / "reads.bam"
with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
    for index, start in enumerate([100, 110, 120, 300, 310, 320, 400, 410]):
        read = pysam.AlignedSegment()
        read.query_name = f"q{index}"
        read.query_sequence = "A" * 50
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, 50)]
        out.write(read)
pysam.index(str(bam))
result = count_bam_reads_per_bin([bam], bin_size=100)
print([row["counts"][0] for row in result["bins"][:5]])
```

- **示例输出：**

```text
[0, 3, 0, 3, 2]
```

## 18) `dnakit.epigenome.estimate_fragment_lengths()`

- **作用：** 从已索引的 BAM 中抽样，估计片段长度和读长的分布。
- **处理方法：** 调用 deeptools 的 `get_read_and_fragment_length`。片段长度取 proper pair 中 read1 的模板长度；如果文件中没有 proper pair，则改用全部映射记录。估计基于抽样的基因组区域，不是对每条记录的统计。
- **API：** `dnakit.epigenome.estimate_fragment_lengths(bam_file, *, number_of_processors=1)`。
- **输入：** 已索引的 BAM 路径；`number_of_processors`。
- **输出：** `bam`、`fragment_length` 和 `read_length`（每个都是包含 `sample_size`、`min`、`median`、`mean`、`std` 和若干分位数的字典，文件无数据时为 `null`）以及 `provenance`。
- **示例代码：**

```python
import tempfile
from pathlib import Path
import pysam
from dnakit.epigenome import estimate_fragment_lengths

work = Path(tempfile.mkdtemp())
def make_read(name, flag, start, tlen=0):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = [(0, 50)]
    read.next_reference_id = 0
    read.template_length = tlen
    return read

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
reads = []
for index, start in enumerate([100, 110, 120, 400, 410]):
    reads.append(make_read(f"q{index}", 99, start, 250))
    reads.append(make_read(f"q{index}", 147, start + 200, -250))
reads.sort(key=lambda read: read.reference_start)
bam = work / "pairs.bam"
with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
    for read in reads:
        out.write(read)
pysam.index(str(bam))
result = estimate_fragment_lengths(bam)
print(result["fragment_length"]["sample_size"], result["fragment_length"]["median"], result["read_length"]["median"])
```

- **示例输出：**

```text
5 250.0 50.0
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-cooler]"          # cooler: build_contact_matrix / fetch_contact_matrix / balance_contact_matrix
python -m pip install "dnakit[bio-cooltools]"       # cooltools: calculate_* on .cool files
python -m pip install "dnakit[bio-pairtools]"       # pairtools: read_pairs_header / summarize_pairs / deduplicate_pairs
python -m pip install "dnakit[bio-macs3]"           # MACS3: build_* pileups and call_fragment_peaks
python -m pip install "dnakit[bio-methylpy]"        # methylpy: allc filtering and binomial tests
python -m pip install "dnakit[bio-deeptools]"       # deeptools; also needs pyBigWig and pysam (bio-pybigwig, bio-reads)
```

各包的项目地址：[cooler](https://github.com/open2c/cooler)、[cooltools](https://github.com/open2c/cooltools)、[pairtools](https://github.com/mirnylab/pairtools)、[MACS3](https://github.com/macs3-project/MACS)、[methylpy](https://github.com/yupenghe/methylpy)、[deeptools](https://github.com/deeptools/deepTools)、[pyBigWig](https://github.com/deeptools/pyBigWig)、[pysam](https://github.com/pysam-developers/pysam)。

methylpy 1.4.7 在导入时需要 `pkg_resources`，它随 setuptools 80 及更早版本提供；如果环境中缺少该模块，请安装 `setuptools<81`。

## 验证方式

测试文件是 `others/tests/unit/test_epigenome_backends.py`，共 34 项，主要检查：

- 接触矩阵：`build_contact_matrix()` 的计数与手工结果一致，包括像素 (0,1)=2、(0,2)=1、(1,10)=1、(5,5)=1、(10,10)=1、(11,14)=1；`fetch_contact_matrix()` 与直接调用 `Cooler.matrix().fetch()` 相同；`balance_contact_matrix()` 写入的权重与直接调用 `balance_cooler` 相同。
- cooltools：`calculate_coverage()` 的 cis/total 与手工求和及直接调用 `cooltools.coverage` 一致；`calculate_expected_contacts()` 的对角线均值与手工结果及直接调用 `expected_cis` 一致；`calculate_insulation_score()` 与直接调用 `insulation` 一致，并在合成数据中检出约 2 kb 处的边界。
- pairtools：`summarize_pairs()` 的计数；`deduplicate_pairs()` 的链式重复簇（a、b、c 为一簇）、链方向不同的对不算重复，以及单侧映射对的输出，均按手工结果核对。
- MACS3：分段覆盖度、延伸 tag 覆盖度和峰顶与手工结果一致。
- methylpy：非转化率为手工比值 1/24；二项 p 值与手工公式（0.01²、1 − 0.99²）及 scipy 一致；allc 过滤结果手工核对。
- deeptools：bigWig 窗口均值、缺失碱基的处理、BAM 窗口计数（含 MAPQ 过滤）和片段长度与手工结果一致。
- 缺少可选依赖时，每个模块都给出对应的 `bio-*` 安装提示；缺失的输入、已存在的输出、超过 `max_bins` 的 bin 数、不安全的 shell 输出路径等都会被拒绝。
- 测试不访问网络，也没有 `DNAKIT_LIVE_TESTS` 测试，因为这一批没有网络服务。

## 当前范围

这一批包装了六个包中每个包的三个操作，共 18 个函数。未包装的内容如下。

- cooler：`cooler cload` 的 BAM 和 `.hic` 输入、`.mcool` 多分辨率文件、`zoomify`、`merge` 和 `coarsen`。
- cooltools：compartment 特征向量、dot 检出、saddle、pileup/snipping 和 directionality。
- pairtools：`parse`（SAM 转 pairs）、`sort`、`select`、`split`、`restrict` 和 `flip`；统计的 YAML 输出与过滤表达式。输入只测试了未压缩的纯文本 pairs。
- MACS3：`callpeak` 的统计模型。这里只有朴素的峰检出。BAM 输入的 `callpeak` 已由 `dnakit.ngs.call_chip_peaks()` 通过 CLI 提供，未重复包装。
- methylpy：比对流程、基于 BAM 的 `call_methylated_sites`、DMR 检出（DMRfind）以及区域甲基化水平（`get_methylation_levels_DMRfind`）。
- deeptools：`bamCoverage` 写 bigWig、`computeMatrix` 与热图、相关性和 PCA，以及 GC 偏好校正。
- 与已有适配器的关系：`dnakit.integrations.summarize_pybigwig_region()` 是 pyBigWig 的单区域统计，与这里的 deeptools 窗口矩阵不同。`dnakit.epigenome.extract_methylation_methyldackel()` 调用 MethylDackel 的 CLI，与 methylpy 的函数无关。
