# G10：cooler、cooltools 与 pairtools 的文件与矩阵函数

这一组函数包装了三个上游包中与 DNA 或基因组坐标相关、可离线运行的公开函数，共 42 个。`cooler_functions.py` 包含 11 个函数（染色体长度、基因组区间、FASTA 长度、接触矩阵文件的截取、合并、压缩与多分辨率生成）；`cooltools_functions.py` 包含 20 个函数（期望值、区室、saddle、virtual 4C、pileup、视窗、稠密矩阵的校正与特征分解）；`pairtools_functions.py` 包含 11 个函数（限制性片段、scaling 曲线、基因组区域、读段错配、SAM 比对解析、.pairs 头部与单倍型）。

这些函数位于 `dnakit.epigenome` 的三个模块中（`dnakit.epigenome` 等），调用时才导入上游包；没有安装对应可选依赖时抛出 `BackendUnavailableError`。输入先检查，上游的 `ValueError`、`KeyError`、`OSError` 转为 `ConfigurationError`。结果为普通字典，NaN 返回为 `None`，分组结果显式排序。`fetch_chrom_sizes` 需要网络，默认测试使用模拟客户端，文档示例只展示参数检查。

示例中的接触矩阵由已有的 `build_contact_matrix()` 从下面的 `contacts.pairs` 构建，并由 `balance_contact_matrix()` 平衡；每段示例在当前目录下的 `g10_out/` 中写入输出。

## 示例文件

以下文件在运行示例前写入示例目录。

**示例文件：** `chrom.sizes`

```text
chr1	1000
chr2	800
```

**示例文件：** `contacts.pairs`

```text
## pairs format v1.0
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r0	chr2	649	chr1	86	+	-	UU
r1	chr1	182	chr1	447	+	-	UU
r2	chr1	801	chr1	843	+	-	UU
r3	chr2	266	chr2	319	+	-	UU
r4	chr1	433	chr1	446	+	-	UU
r5	chr1	691	chr1	699	+	-	UU
r6	chr2	587	chr2	642	+	-	UU
r7	chr1	887	chr1	919	+	-	UU
r8	chr1	517	chr1	672	+	-	UU
r9	chr2	139	chr2	384	+	-	UU
r10	chr1	738	chr1	847	+	-	UU
r11	chr1	320	chr2	519	+	-	UU
r12	chr2	695	chr2	696	+	-	UU
r13	chr1	293	chr1	433	+	-	UU
r14	chr1	943	chr1	999	+	-	UU
r15	chr2	239	chr2	318	+	-	UU
r16	chr1	663	chr1	734	+	-	UU
r17	chr1	585	chr1	656	+	-	UU
r18	chr2	380	chr2	639	+	-	UU
r19	chr1	31	chr1	201	+	-	UU
r20	chr1	254	chr1	429	+	-	UU
r21	chr2	73	chr2	441	+	-	UU
r22	chr1	928	chr2	166	+	-	UU
r23	chr1	248	chr1	396	+	-	UU
r24	chr2	239	chr2	495	+	-	UU
r25	chr1	654	chr1	817	+	-	UU
r26	chr1	219	chr1	358	+	-	UU
r27	chr2	4	chr2	217	+	-	UU
r28	chr1	683	chr1	730	+	-	UU
r29	chr1	957	chr1	999	+	-	UU
r30	chr2	607	chr2	631	+	-	UU
r31	chr1	593	chr1	679	+	-	UU
r32	chr1	849	chr1	912	+	-	UU
r33	chr2	309	chr1	147	+	-	UU
r34	chr1	799	chr1	975	+	-	UU
r35	chr1	292	chr1	366	+	-	UU
r36	chr2	419	chr2	558	+	-	UU
r37	chr1	562	chr1	999	+	-	UU
r38	chr1	931	chr1	954	+	-	UU
r39	chr2	158	chr2	375	+	-	UU
```

**示例文件：** `expected.tsv`

```text
region1	region2	dist	n_valid	count.avg	balanced.avg
chr1	chr1	0	10	1.2	0.9
chr1	chr1	1	9	0.8	0.7
chr1	chr1	2	8	0.5	0.6
chr2	chr2	0	6	1.0	0.8
```

**示例文件：** `features.bed`

```text
chr1	300	400
chr1	500	600
chr2	200	300
```

**示例文件：** `genome.fa`

```text
>chr1
ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA
>chr2
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT
```

**示例文件：** `rearrange_view.tsv`

```text
chrom	start	end	name	new_chrom	strand
chr1	0	1000	chr1	chrX	+
chr2	0	800	chr2	chrX	-
```

**示例文件：** `regions.bed`

```text
chr1	0	150	x
chr1	200	300	y
chr2	10	20	z
```

**示例文件：** `sample.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:800
r1	65	chr1	101	60	10M	=	101	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r1	129	chr2	301	60	10M	=	301	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
```

**示例文件：** `track.bed`

```text
chr1	0	100	-1.000000
chr1	100	200	-0.777778
chr1	200	300	-0.555556
chr1	300	400	-0.333333
chr1	400	500	-0.111111
chr1	500	600	0.111111
chr1	600	700	0.333333
chr1	700	800	0.555556
chr1	800	900	0.777778
chr1	900	1000	1.000000
chr2	0	100	1.000000
chr2	100	200	0.714286
chr2	200	300	0.428571
chr2	300	400	0.142857
chr2	400	500	-0.142857
chr2	500	600	-0.428571
chr2	600	700	-0.714286
chr2	700	800	-1.000000
```

**示例文件：** `view.bed`

```text
chr1	0	1000	chr1
chr2	0	800	chr2
```


## cooler 函数（11 个）

## 1) `dnakit.epigenome.read_chrom_sizes()`

- **作用：** 读取 UCSC 风格的 chrom.sizes 或 chromInfo 文件，返回染色体长度。
- **处理方法：** 调用 cooler.util.read_chromsizes。默认只保留 chr1-22、chrX、chrY、chrM 形式的名字；`all_names=True` 保留全部；`name_patterns` 可自定义正则。结果保持文件中的顺序。
- **API：** `dnakit.epigenome.read_chrom_sizes(path, *, all_names=False, name_patterns=None, verbose=True)`。
- **输入：** `path`：两列（名字、长度）或 UCSC chromInfo 文件；`all_names`、`name_patterns`（正则列表）。
- **输出：** 含 `chromosome_count` 和 `chromosomes`（每项为 `name`、`length`）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import read_chrom_sizes

result = read_chrom_sizes("chrom.sizes", all_names=True, verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **示例输出：**

```text
chr1 1000
chr2 800
```

## 2) `dnakit.epigenome.fetch_chrom_sizes()`

- **作用：** 下载 UCSC 基因组版本的染色体长度表（需要网络）。
- **处理方法：** 先检查 `database` 是否为 2 到 64 个字母、数字或下划线组成的 UCSC 名字，然后用 cooler.util.fetch_chromsizes 请求 `chromInfo.txt.gz`。网络错误转为 `BackendUnavailableError`。测试用模拟的 HTTP 客户端，默认测试不访问网络。
- **API：** `dnakit.epigenome.fetch_chrom_sizes(database, *, all_names=False, name_patterns=None, verbose=True)`。
- **输入：** `database`，例如 `hg38`；`all_names`、`name_patterns` 与 `read_chrom_sizes` 相同。
- **输出：** 与 `read_chrom_sizes` 相同的字典，`source` 字段记录下载地址。
- **示例代码：**

```python
from dnakit.epigenome import fetch_chrom_sizes
from dnakit.exceptions import ConfigurationError

try:
    fetch_chrom_sizes("hg 38", verbose=False)
except ConfigurationError as exc:
    print("rejected:", str(exc).split("] ", 1)[-1])
```

- **示例输出：**

```text
rejected: database must be a UCSC assembly name such as hg38.
```

## 3) `dnakit.epigenome.parse_genomic_region()`

- **作用：** 把 `chr1:100-200` 这样的区域字符串解析为半开区间（起点从 0 开始，终点从 1 开始）。
- **处理方法：** 调用 cooler.util.parse_region。只写染色体名时，若给出 `chromsizes`，区间为整条染色体。超出染色体长度时报错。
- **API：** `dnakit.epigenome.parse_genomic_region(region, *, chromsizes=None, verbose=True)`。
- **输入：** `region`：字符串；`chromsizes`：可选的染色体长度映射。
- **输出：** 含 `chrom`、`start`、`end`、`length` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import parse_genomic_region

region = parse_genomic_region("chr1:100-200", verbose=False)
print(region["chrom"], region["start"], region["end"], region["length"])
whole = parse_genomic_region("chr2", chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)
print(whole["start"], whole["end"])
```

- **示例输出：**

```text
chr1 100 200 100
0 800
```

## 4) `dnakit.epigenome.query_bed_intervals()`

- **作用：** 从按起点排序、互不重叠的 BED 文件中取出某个区域内的区间。
- **处理方法：** 检查区间的起点和终点、同一染色体内的排序与重叠，然后调用 cooler.util.bedslice。返回与区域有交叠的区间，按起点和终点排序。
- **API：** `dnakit.epigenome.query_bed_intervals(bed_path, region, *, chromsizes, verbose=True)`。
- **输入：** `bed_path`：至少三列（chrom、start、end）的 BED 文件；`region`；`chromsizes`。
- **输出：** 含 `interval_count` 和 `intervals`（每项为 `chrom`、`start`、`end`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import query_bed_intervals

result = query_bed_intervals("regions.bed", "chr1:100-250", chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)
for row in result["intervals"]:
    print(row["chrom"], row["start"], row["end"])
```

- **示例输出：**

```text
chr1 0 150
chr1 200 300
```

## 5) `dnakit.epigenome.rename_cooler_chromosomes()`

- **作用：** 复制一个 cooler 文件并在副本中改写染色体名，原文件不变。
- **处理方法：** 检查新名字唯一且旧名字存在，复制文件，然后调用 cooler.rename_chroms 修改副本。
- **API：** `dnakit.epigenome.rename_cooler_chromosomes(cool_path, output_path, rename_map, *, verbose=True)`。
- **输入：** `cool_path`：输入的 .cool 文件；`output_path`：新文件（不能已存在）；`rename_map`：旧名字到新名字的映射。
- **输出：** 含 `chromosomes_before`、`chromosomes_after` 和 `renamed` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import rename_cooler_chromosomes

result = rename_cooler_chromosomes(out / "in.cool", out / "renamed.cool", {"chr1": "1"}, verbose=False)
print(result["chromosomes_before"], result["chromosomes_after"])
```

- **示例输出：**

```text
['chr1', 'chr2'] ['1', 'chr2']
```

## 6) `dnakit.epigenome.coarsen_contact_matrix()`

- **作用：** 把接触矩阵按整数倍数降低分辨率（k×k 相邻像素求和）。
- **处理方法：** 调用 cooler.coarsen_cooler，逐条染色体分块处理。输出的分箱大小为输入的 `factor` 倍，接触总数保持不变。
- **API：** `dnakit.epigenome.coarsen_contact_matrix(cool_path, output_path, *, factor, chunksize=10000000, verbose=True)`。
- **输入：** `cool_path`；`output_path`（新文件）；`factor`（2 到 1000 的整数）；`chunksize`。
- **输出：** 含 `input_binsize`、`output_binsize`、`output_bins` 和接触数统计的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import coarsen_contact_matrix

result = coarsen_contact_matrix(out / "in.cool", out / "coarse.cool", factor=2, verbose=False)
print(result["input_binsize"], result["output_binsize"], result["output_total_contacts"])
```

- **示例输出：**

```text
100 200 40
```

## 7) `dnakit.epigenome.merge_contact_matrices()`

- **作用：** 把两个或更多分箱完全相同的 cooler 文件相加为一个文件。
- **处理方法：** 检查输入不重复且数量为 2 到 64，然后调用 cooler.merge_coolers。输出接触数为各输入之和。
- **API：** `dnakit.epigenome.merge_contact_matrices(cool_paths, output_path, *, mergebuf=10000000, verbose=True)`。
- **输入：** `cool_paths`：文件路径列表；`output_path`（新文件）；`mergebuf`。
- **输出：** 含 `input_total_contacts`、`output_total_contacts` 和 `output_nonzero_pixels` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome.backends import build_contact_matrix
from dnakit.epigenome import merge_contact_matrices

build_contact_matrix("contacts.pairs", out / "second.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
result = merge_contact_matrices([out / "in.cool", out / "second.cool"], out / "merged.cool", verbose=False)
print(result["input_total_contacts"], result["output_total_contacts"])
```

- **示例输出：**

```text
80 80
```

## 8) `dnakit.epigenome.zoomify_contact_matrix()`

- **作用：** 从一个基础 cooler 生成多分辨率（.mcool）文件。
- **处理方法：** 调用 cooler.zoomify_cooler，逐级递归合并，每个分辨率存放在 `/resolutions/<分箱大小>` 下。返回每一级的分箱数和接触数。
- **API：** `dnakit.epigenome.zoomify_contact_matrix(cool_paths, output_path, *, resolutions, chunksize=10000000, verbose=True)`。
- **输入：** `cool_paths`（一个或多个基础文件）；`output_path`（新文件）；`resolutions`（升序、互不相同的分箱大小列表）。
- **输出：** 含 `resolutions` 列表（每项为 `resolution`、`binsize`、`bins`、`total_contacts`）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import zoomify_contact_matrix

result = zoomify_contact_matrix(out / "in.cool", out / "multi.mcool", resolutions=[100, 200, 400], verbose=False)
for level in result["resolutions"]:
    print(level["resolution"], level["bins"], level["total_contacts"])
```

- **示例输出：**

```text
100 18 40
200 9 40
400 5 40
```

## 9) `dnakit.epigenome.load_fasta_lengths()`

- **作用：** 读取 FASTA 中指定序列的长度，不把碱基读入内存。
- **处理方法：** 调用 cooler.util.load_fasta（基于 pyfaidx 的懒加载索引）。返回顺序与 `names` 相同。名字不存在时报错。
- **API：** `dnakit.epigenome.load_fasta_lengths(fasta_path, names, *, verbose=True)`。
- **输入：** `fasta_path`：FASTA 文件；`names`：序列名字列表。
- **输出：** 含 `sequence_count` 和 `sequences`（每项为 `name`、`length`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import load_fasta_lengths

result = load_fasta_lengths("genome.fa", ["chr2", "chr1"], verbose=False)
for row in result["sequences"]:
    print(row["name"], row["length"])
```

- **示例输出：**

```text
chr2 800
chr1 1000
```

## 10) `dnakit.epigenome.read_cooler_pixels()`

- **作用：** 读取 cooler 文件中前若干个非零像素（上三角接触），可选地连接分箱坐标。
- **处理方法：** 调用 cooler.Cooler.pixels（api.pixels），按存储顺序（bin1_id、bin2_id）切片。`join=True` 时加入两端的染色体、起点和终点。
- **API：** `dnakit.epigenome.read_cooler_pixels(cool_path, *, limit=100, join=False, verbose=True)`。
- **输入：** `cool_path`；`limit`（1 到 10000）；`join`（布尔值）。
- **输出：** 含 `nonzero_pixels`、`returned`、`columns` 和 `pixels` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import read_cooler_pixels

result = read_cooler_pixels(out / "in.cool", limit=3, join=True, verbose=False)
print(result["nonzero_pixels"], result["returned"])
for row in result["pixels"]:
    print(row["chrom1"], row["start1"], row["chrom2"], row["start2"], row["count"])
```

- **示例输出：**

```text
28 3
chr1 0 chr1 200 1
chr1 0 chr2 600 1
chr1 100 chr1 400 1
```

## 11) `dnakit.epigenome.read_cooler_chroms()`

- **作用：** 列出 cooler 文件中存储的染色体（按接触矩阵的顺序）。
- **处理方法：** 调用 cooler.Cooler.chroms（api.chroms）读取染色体表。
- **API：** `dnakit.epigenome.read_cooler_chroms(cool_path, *, verbose=True)`。
- **输入：** `cool_path`：.cool 文件。
- **输出：** 含 `chromosome_count` 和 `chromosomes`（每项为 `name`、`length`）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "in.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
from dnakit.epigenome import read_cooler_chroms

result = read_cooler_chroms(out / "in.cool", verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **示例输出：**

```text
chr1 1000
chr2 800
```

## cooltools 函数（20 个）

## 12) `dnakit.epigenome.calculate_expected_trans()`

- **作用：** 计算每一对染色体之间的平均接触频率（跨染色体期望值）。
- **处理方法：** 调用 cooltools.expected_trans，使用平衡权重列（默认 `weight`）。结果按 region1、region2 排序。
- **API：** `dnakit.epigenome.calculate_expected_trans(cool_path, *, weight='weight', verbose=True)`。
- **输入：** `cool_path`：已平衡的 .cool 文件；`weight`：权重列名。
- **输出：** 含 `pair_count` 和 `pairs`（每项含 `region1`、`region2`、`count.sum`、`balanced.sum` 等）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_expected_trans

result = calculate_expected_trans(out / "bal.cool", verbose=False)
print(result["pair_count"])
for row in result["pairs"]:
    print(row["region1"], row["region2"], row["count.sum"], round(row["balanced.sum"], 4))
```

- **示例输出：**

```text
1
chr1 chr2 4.0 1.8912
```

## 13) `dnakit.epigenome.calculate_compartments()`

- **作用：** 计算顺式 A/B 区室的特征值和特征向量（每条染色体）。
- **处理方法：** 调用 cooltools.eigs_cis，视图为整条染色体。在 ARPACK 起始向量和 C 随机种子固定的条件下运行，使结果可重复。每个特征向量的绝对值最大的元素为正。
- **API：** `dnakit.epigenome.calculate_compartments(cool_path, *, n_eigs=3, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`n_eigs`（1 到 10）；`weight`。
- **输出：** 含 `eigenvalues`（每条染色体一行）和 `eigenvectors`（每个分箱一行，列 E1、E2 ...）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_compartments

result = calculate_compartments(out / "bal.cool", n_eigs=2, verbose=False)
for row in result["eigenvalues"]:
    print(row["name"], round(row["eigval1"], 4), round(row["eigval2"], 4))
print(len(result["eigenvectors"]))
```

- **示例输出：**

```text
chr1 8.0892 -6.6156
chr2 4.3574 -4.3272
18
```

## 14) `dnakit.epigenome.calculate_trans_compartments()`

- **作用：** 从跨染色体接触计算全基因组的反式区室特征向量。
- **处理方法：** 调用 cooltools.eigs_trans，分区为全部染色体。固定随机输入（同上）。特征向量符号按绝对值最大元素为正统一。
- **API：** `dnakit.epigenome.calculate_trans_compartments(cool_path, *, n_eigs=3, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`n_eigs`；`weight`。
- **输出：** 含 `eigenvalues`、`eigenvectors` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_trans_compartments

result = calculate_trans_compartments(out / "bal.cool", n_eigs=1, verbose=False)
print(round(result["eigenvalues"][0]["eigval1"], 4), len(result["eigenvectors"]))
```

- **示例输出：**

```text
8.0001 18
```

## 15) `dnakit.epigenome.calculate_saddle()`

- **作用：** 计算 saddle 图：按基因组轨道的分位数分箱后，各箱之间的平均平衡接触。
- **处理方法：** 用 cooltools.expected_cis 或 expected_trans 计算期望值，再调用 cooltools.saddle（`qrange` 给出分位范围）。输出矩阵比分箱数多两行两列，外侧一圈表示低于或高于分位范围的值。
- **API：** `dnakit.epigenome.calculate_saddle(cool_path, track_path, *, contact_type='cis', n_bins=5, qrange=(0.02, 0.98), min_diag=3, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`track_path`：含 chrom、start、end、value 四列的轨道文件；`n_bins`；`qrange`（低、高分位）；`min_diag`。
- **输出：** 含 `matrix_size`、`interaction_sum` 和 `interaction_count` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_saddle

result = calculate_saddle(out / "bal.cool", "track.bed", n_bins=3, qrange=(0.0, 1.0), min_diag=1, verbose=False)
print(result["matrix_size"])
print([round(value, 4) if value is not None else None for value in result["interaction_sum"][1]])
```

- **示例输出：**

```text
5
[0.0, 6.902, 26.0064, 0.0, 0.0]
```

## 16) `dnakit.epigenome.sample_contacts()`

- **作用：** 随机抽取一部分接触，写入新的 cooler 文件。
- **处理方法：** 给出 `count` 时用 cooltools.sample 的精确抽样；给出 `fraction` 时用二项抽样。两者只能给一个。随机种子只在调用期间设定，之后恢复 numpy 的全局随机状态。
- **API：** `dnakit.epigenome.sample_contacts(cool_path, output_path, *, count=None, fraction=None, seed=0, verbose=True)`。
- **输入：** `cool_path`；`output_path`（新文件）；`count` 或 `fraction`；`seed`。
- **输出：** 含 `method`、`seed`、`total_contacts_before` 和 `total_contacts_after` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import sample_contacts

result = sample_contacts(out / "bal.cool", out / "sampled.cool", count=20, seed=1, verbose=False)
print(result["method"], result["total_contacts_before"], result["total_contacts_after"])
```

- **示例输出：**

```text
exact 40 20
```

## 17) `dnakit.epigenome.calculate_virtual_4c()`

- **作用：** 计算某个视点（viewpoint）与全基因组各分箱的接触轮廓（virtual 4C）。
- **处理方法：** 解析 `chr1:0-100` 形式的视点，调用 cooltools.virtual4c（上游标注为实验性接口）。结果中缺失值返回为 None。
- **API：** `dnakit.epigenome.calculate_virtual_4c(cool_path, viewpoint, *, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`viewpoint`：区域字符串；`weight`。
- **输出：** 含 `viewpoint`、`bin_count` 和 `profile`（每个分箱一行）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_virtual_4c

result = calculate_virtual_4c(out / "bal.cool", "chr1:0-100", verbose=False)
print(result["bin_count"], result["viewpoint"]["chrom"], result["viewpoint"]["end"])
print([row["balanced"] is None or round(row["balanced"], 4) for row in result["profile"][:3]])
```

- **示例输出：**

```text
18 chr1 100
[True, True, 0.5738]
```

## 18) `dnakit.epigenome.calculate_directionality_index()`

- **作用：** 计算每个分箱的方向性指数（DI），用于找出 TAD 边界。
- **处理方法：** 调用 cooltools.directionality，窗口为 `window_bp` 碱基（必须不小于分箱大小）。缺失值返回为 None。
- **API：** `dnakit.epigenome.calculate_directionality_index(cool_path, *, window_bp=100000, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`window_bp`；`weight`。
- **输出：** 含 `bin_count` 和 `bins`（每项含 `chrom`、`start`、`end` 与 DI 列）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_directionality_index

result = calculate_directionality_index(out / "bal.cool", window_bp=300, verbose=False)
print(result["bin_count"])
print([row["directionality_index_300"] is None or round(row["directionality_index_300"], 4) for row in result["bins"][:4]])
```

- **示例输出：**

```text
18
[True, True, 0.0069, True]
```

## 19) `dnakit.epigenome.pileup_features()`

- **作用：** 围绕 BED 文件中的每个特征截取平衡接触小片段，并求平均（pileup）。
- **处理方法：** 用 cooltools.pileup 取得每个特征的片段，然后对有限值求平均（忽略 NaN）。不做期望值归一化。
- **API：** `dnakit.epigenome.pileup_features(cool_path, features_path, *, flank_bp=100000, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`features_path`：至少含 chrom、start、end 的 BED 文件；`flank_bp`；`weight`。
- **输出：** 含 `feature_count`、`features_with_data`、`snippet_size` 和 `mean_snippet` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import pileup_features

result = pileup_features(out / "bal.cool", "features.bed", flank_bp=200, verbose=False)
print(result["feature_count"], result["features_with_data"], result["snippet_size"])
```

- **示例输出：**

```text
3 3 5
```

## 20) `dnakit.epigenome.calculate_logbin_expected()`

- **作用：** 把期望接触频率随距离的变化按对数间隔分箱，并给出各段的斜率。
- **处理方法：** 先用 cooltools.expected_cis 计算期望值，再调用 cooltools.api.expected.logbin_expected。返回分箱表、斜率表和分箱边界。
- **API：** `dnakit.epigenome.calculate_logbin_expected(cool_path, *, bins_per_order_magnitude=10, min_nvalid=200, min_count=50, weight='weight', verbose=True)`。
- **输入：** `cool_path`；`bins_per_order_magnitude`（每个数量级的分箱数）；`min_nvalid`；`min_count`；`weight`。
- **输出：** 含 `binned_rows`、`binned`、`slopes` 和 `bin_edges` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import calculate_logbin_expected

result = calculate_logbin_expected(out / "bal.cool", bins_per_order_magnitude=4, min_nvalid=1, min_count=1, verbose=False)
print(result["binned_rows"], result["bin_edges"])
```

- **示例输出：**

```text
4 [1, 2, 3, 6]
```

## 21) `dnakit.epigenome.rearrange_contact_matrix()`

- **作用：** 按视窗（viewframe）给出的新顺序重排 cooler 的分箱和像素，写入新文件。
- **处理方法：** 读取含 chrom、start、end、name 以及新染色体名和链方向列的制表符表格，调用 cooltools.api.rearrange.rearrange_cooler。`-` 链的区域会被反转。
- **API：** `dnakit.epigenome.rearrange_contact_matrix(cool_path, view_path, output_path, *, new_chrom_column='new_chrom', orientation_column='strand', verbose=True)`。
- **输入：** `cool_path`；`view_path`：视窗表；`output_path`（新文件）；`new_chrom_column`；`orientation_column`。
- **输出：** 含 `chromosomes_before`、`chromosomes_after`、`output_bins` 和 `output_total_contacts` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import rearrange_contact_matrix

result = rearrange_contact_matrix(out / "bal.cool", "rearrange_view.tsv", out / "rearranged.cool", verbose=False)
print(result["chromosomes_before"], result["chromosomes_after"], result["output_total_contacts"])
```

- **示例输出：**

```text
['chr1', 'chr2'] ['chrX'] 40
```

## 22) `dnakit.epigenome.align_windows_to_bins()`

- **作用：** 把基因组位点转换为与分箱边界对齐的窗口。
- **处理方法：** 以位点为中心，左右各取 `flank_bp`，再把窗口两端移到分箱边界上。调用 cooltools.api.snipping.make_bin_aligned_windows。`lo` 和 `hi` 是相对于 `region_start_bp` 的分箱序号。
- **API：** `dnakit.epigenome.align_windows_to_bins(binsize, chroms, centers_bp, *, flank_bp=0, region_start_bp=0, verbose=True)`。
- **输入：** `binsize`；`chroms`；`centers_bp`（与 chroms 等长）；`flank_bp`；`region_start_bp`。
- **输出：** 含 `window_count` 和 `windows`（每项为 `chrom`、`start`、`end`、`lo`、`hi`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import align_windows_to_bins

result = align_windows_to_bins(100, ["chr1"], [450], flank_bp=200, verbose=False)
window = result["windows"][0]
print(window["chrom"], window["start"], window["end"], window["lo"], window["hi"])
```

- **示例输出：**

```text
chr1 200 700 2 7
```

## 23) `dnakit.epigenome.assign_features_to_view()`

- **作用：** 为 BED 文件中的每个特征标注它所在的视窗区域名。
- **处理方法：** 用 cooltools.lib.common.assign_view_auto 做区间重叠判断。不重叠的特征返回 `region` 为 None，并计入 `unassigned_count`。特征顺序保持不变。
- **API：** `dnakit.epigenome.assign_features_to_view(features_path, view_path, *, verbose=True)`。
- **输入：** `features_path`：BED 文件；`view_path`：BED 文件，第四列为区域名。
- **输出：** 含 `feature_count`、`unassigned_count` 和 `features`（每项含 `chrom`、`start`、`end`、`region`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import assign_features_to_view

result = assign_features_to_view("features.bed", "view.bed", verbose=False)
for row in result["features"]:
    print(row["chrom"], row["start"], row["region"])
```

- **示例输出：**

```text
chr1 300 chr1
chr1 500 chr1
chr2 200 chr2
```

## 24) `dnakit.epigenome.make_chromosome_view()`

- **作用：** 由 cooler 的染色体长度生成全染色体视窗（每条染色体一个区域）。
- **处理方法：** 调用 cooltools.lib.common.make_cooler_view。`ucsc_names=True` 时使用 UCSC 风格的区域名。
- **API：** `dnakit.epigenome.make_chromosome_view(cool_path, *, ucsc_names=False, verbose=True)`。
- **输入：** `cool_path`；`ucsc_names`。
- **输出：** 含 `region_count` 和 `regions`（每项含 `chrom`、`start`、`end`、`name`）的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import make_chromosome_view

result = make_chromosome_view(out / "bal.cool", verbose=False)
for row in result["regions"]:
    print(row["name"], row["start"], row["end"])
```

- **示例输出：**

```text
chr1 0 1000
chr2 0 800
```

## 25) `dnakit.epigenome.read_viewframe_bed()`

- **作用：** 读取并检查视窗 BED 文件（区域名唯一、互不重叠）。
- **处理方法：** 调用 cooltools.lib.io.read_viewframe_from_file。给出 `cool_path` 时，还检查区域是否落在染色体长度之内。
- **API：** `dnakit.epigenome.read_viewframe_bed(view_path, *, cool_path=None, verbose=True)`。
- **输入：** `view_path`：BED 文件；`cool_path`（可选）。
- **输出：** 含 `region_count`、`regions` 和 `verified_against` 的字典。
- **示例代码：**

```python
import shutil
from pathlib import Path
from dnakit.epigenome.backends import balance_contact_matrix, build_contact_matrix

out = Path("g10_out")
shutil.rmtree(out, ignore_errors=True)
out.mkdir()
build_contact_matrix("contacts.pairs", out / "bal.cool", chromsizes={"chr1": 1000, "chr2": 800}, binsize=100, verbose=False)
balance_contact_matrix(out / "bal.cool", min_nnz=1, verbose=False)
from dnakit.epigenome import read_viewframe_bed

result = read_viewframe_bed("view.bed", cool_path=out / "bal.cool", verbose=False)
print(result["region_count"], [row["name"] for row in result["regions"]])
```

- **示例输出：**

```text
2 ['chr1', 'chr2']
```

## 26) `dnakit.epigenome.read_expected_table()`

- **作用：** 读取 cooltools v1.0 格式的期望接触表，并检查列的类型。
- **处理方法：** 调用 cooltools.lib.io.read_expected_from_file（顺式表的列为 region1、region2、dist、n_valid 与数值列）。返回的行数不超过 `limit`。
- **API：** `dnakit.epigenome.read_expected_table(expected_path, *, contact_type='cis', limit=100, verbose=True)`。
- **输入：** `expected_path`：制表符分隔文件；`contact_type`（`cis` 或 `trans`）；`limit`。
- **输出：** 含 `row_count`、`columns` 和 `rows` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import read_expected_table

result = read_expected_table("expected.tsv", limit=2, verbose=False)
print(result["row_count"], result["columns"])
print(result["rows"][1]["dist"], result["rows"][1]["balanced.avg"])
```

- **示例输出：**

```text
4 ['region1', 'region2', 'dist', 'n_valid', 'count.avg', 'balanced.avg']
1 0.7
```

## 27) `dnakit.epigenome.observed_over_expected_matrix()`

- **作用：** 把稠密的对称接触矩阵除以各距离的平均值，得到观测/期望（O/E）矩阵。
- **处理方法：** 调用 cooltools.lib.numutils.observed_over_expected。距离按指数增长的区间分组，每条对角线除以其所在区间的平均值。输入必须对称。
- **API：** `dnakit.epigenome.observed_over_expected_matrix(matrix, *, dist_bin_edge_ratio=1.03, verbose=True)`。
- **输入：** `matrix`：方阵（列表的列表）；`dist_bin_edge_ratio`（大于 1，且不超过 2）。
- **输出：** 含 `observed_over_expected`、`distance_bin_edges`、`sum_pixels` 和 `n_pixels` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import observed_over_expected_matrix

matrix = [[1.0, 2.0, 3.0], [2.0, 1.0, 2.0], [3.0, 2.0, 1.0]]
result = observed_over_expected_matrix(matrix, verbose=False)
print([[round(value, 4) for value in row] for row in result["observed_over_expected"]])
```

- **示例输出：**

```text
[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
```

## 28) `dnakit.epigenome.iterative_correction_dense()`

- **作用：** 对稠密接触矩阵做迭代校正（矩阵平衡），对称或非对称均可。
- **处理方法：** 对称矩阵调用 cooltools.lib.numutils.iterative_correction_symmetric（可设 `ignore_diags`）；非对称矩阵调用 `iterative_correction_asymmetric`（此时 `ignore_diags` 必须为 0）。
- **API：** `dnakit.epigenome.iterative_correction_dense(matrix, *, symmetric=True, max_iter=1000, ignore_diags=0, tol=1e-05, verbose=True)`。
- **输入：** `matrix`；`symmetric`；`max_iter`；`ignore_diags`；`tol`。
- **输出：** 含 `converged`、`iterations`、`corrected`、`bias` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import iterative_correction_dense

matrix = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
result = iterative_correction_dense(matrix, max_iter=500, verbose=False)
print(result["converged"], [round(value, 4) for value in result["bias"]])
```

- **示例输出：**

```text
True [0.5618, 1.0536, 1.3846]
```

## 29) `dnakit.epigenome.dense_cis_eigenvectors()`

- **作用：** 对一个稠密顺式接触矩阵做特征分解，得到区室特征值和特征向量。
- **处理方法：** 调用 cooltools.api.eigdecomp.cis_eig，去掉靠近对角线的 `ignore_diags` 条对角线。ARPACK 起始向量固定，特征向量的符号按绝对值最大元素为正统一。
- **API：** `dnakit.epigenome.dense_cis_eigenvectors(matrix, *, n_eigs=3, ignore_diags=2, clip_percentile=0.0, verbose=True)`。
- **输入：** `matrix`：对称方阵；`n_eigs`；`ignore_diags`；`clip_percentile`。
- **输出：** 含 `eigenvalues` 和 `eigenvectors`（每个特征向量一行）的字典。
- **示例代码：**

```python
from dnakit.epigenome import dense_cis_eigenvectors

matrix = [[4.0, 3.0, 1.0, 0.5], [3.0, 4.0, 1.0, 0.5], [1.0, 1.0, 4.0, 3.0], [0.5, 0.5, 3.0, 4.0]]
result = dense_cis_eigenvectors(matrix, n_eigs=1, ignore_diags=0, verbose=False)
print([round(value, 4) for value in result["eigenvalues"]])
```

- **示例输出：**

```text
[0.8095]
```

## 30) `dnakit.epigenome.expected_from_dense_matrix()`

- **作用：** 统计稠密矩阵每条对角线上的有效像素数与接触总和，即顺式期望值的原始汇总。
- **处理方法：** 调用 cooltools.api.expected.diagsum_from_array，按距离排序输出。
- **API：** `dnakit.epigenome.expected_from_dense_matrix(matrix, *, ignore_diags=2, region_name=None, verbose=True)`。
- **输入：** `matrix`：方阵；`ignore_diags`；`region_name`（可选）。
- **输出：** 含 `diagonals`（每项含 `dist`、`n_valid`、`balanced.sum`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import expected_from_dense_matrix

matrix = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
result = expected_from_dense_matrix(matrix, ignore_diags=0, verbose=False)
for row in result["diagonals"]:
    print(row["dist"], row["n_valid"], row["balanced.sum"])
```

- **示例输出：**

```text
0 3 11.0
1 2 7.0
2 1 3.0
```

## 31) `dnakit.epigenome.digitize_track()`

- **作用：** 把基因组信号轨道按分位数分为整数类别，供 saddle 分析使用。
- **处理方法：** 调用 cooltools.api.saddle.digitize。类别 1 到 n 对应分位范围内的值；0 为低于范围，n+1 为高于范围，-1 为缺失值。
- **API：** `dnakit.epigenome.digitize_track(track_path, *, n_bins=5, qrange=(0.02, 0.98), verbose=True)`。
- **输入：** `track_path`：含 chrom、start、end、value 四列的文件；`n_bins`；`qrange`（低、高分位）。
- **输出：** 含 `intervals`（每项含 `chrom`、`start`、`end`、`class`）与 `bin_edges` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import digitize_track

result = digitize_track("track.bed", n_bins=3, qrange=(0.0, 1.0), verbose=False)
print([row["class"] for row in result["intervals"][:8]])
```

- **示例输出：**

```text
[1, 1, 1, 2, 2, 2, 2, 3]
```

## pairtools 函数（11 个）

## 32) `dnakit.epigenome.find_restriction_fragment()`

- **作用：** 找出包含某个位置的限制性酶切片段。
- **处理方法：** 在排序的酶切位点中做二分查找，调用 pairtools.lib.restrict.find_rfrag。片段是相邻两个位点之间的半开区间。位置超出位点范围时报错。
- **API：** `dnakit.epigenome.find_restriction_fragment(sites, chrom, position, *, verbose=True)`。
- **输入：** `sites`：染色体名到升序位点列表的映射；`chrom`；`position`。
- **输出：** 含 `fragment_index`、`start`、`end` 和 `length` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import find_restriction_fragment

sites = {"chr1": [0, 300, 700, 1000], "chr2": [0, 800]}
result = find_restriction_fragment(sites, "chr1", 450, verbose=False)
print(result["fragment_index"], result["start"], result["end"], result["length"])
```

- **示例输出：**

```text
1 300 700 400
```

## 33) `dnakit.epigenome.compute_contact_scaling()`

- **作用：** 从 .pairs 文件计算接触频率随基因组距离的变化（scaling 曲线）。
- **处理方法：** 调用 pairtools.lib.scaling.compute_scaling，再把各区域和链组合的结果按距离区间相加。每个区间的接触频率等于接触数除以可能的碱基对数。
- **API：** `dnakit.epigenome.compute_contact_scaling(pairs_path, *, chromsizes, min_distance=1, max_distance=1000000000, bins_per_decade=8, ignore_trans=False, verbose=True)`。
- **输入：** `pairs_path`；`chromsizes`；`min_distance`；`max_distance`；`bins_per_decade`；`ignore_trans`。
- **输出：** 含 `cis_pairs`、`trans_pairs` 和 `bins`（每项含 `min_dist`、`max_dist`、`n_pairs`、`contact_frequency`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import compute_contact_scaling

result = compute_contact_scaling("contacts.pairs", chromsizes={"chr1": 1000, "chr2": 800}, max_distance=1000, bins_per_decade=2, verbose=False)
print(result["cis_pairs"], result["trans_pairs"])
print(result["bins"][1]["min_dist"], result["bins"][1]["max_dist"], result["bins"][1]["n_pairs"])
```

- **示例输出：**

```text
36 4
1 3 1
```

## 34) `dnakit.epigenome.normalize_scaling_curve()`

- **作用：** 把接触频率曲线除以某个距离窗口内的平均值，使窗口内的平均值为 1。
- **处理方法：** 调用 pairtools.lib.scaling.norm_scaling_factor 和 norm_scaling。距离必须唯一且升序。
- **API：** `dnakit.epigenome.normalize_scaling_curve(distances, frequencies, *, norm_window=(10000, 1000000000), verbose=True)`。
- **输入：** `distances`：整数列表；`frequencies`：同长度的非负数；`norm_window`：低、高两个距离。
- **输出：** 含 `factor` 和 `points`（每项含 `distance`、`frequency`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import normalize_scaling_curve

result = normalize_scaling_curve([1, 10, 100], [0.5, 0.25, 0.1], norm_window=(1, 100), verbose=False)
print(round(result["factor"], 4))
print([round(point["frequency"], 4) for point in result["points"]])
```

- **示例输出：**

```text
0.2833
[1.7647, 0.8824, 0.3529]
```

## 35) `dnakit.epigenome.assign_genomic_regions()`

- **作用：** 判断每个位置落在哪个区域内（例如某条染色体的一段）。
- **处理方法：** 调用 pairtools.lib.scaling.assign_regs。区域索引是同一染色体内按起点排序后的序号，未落入任何区域时为 None。
- **API：** `dnakit.epigenome.assign_genomic_regions(positions, regions, *, verbose=True)`。
- **输入：** `positions`：含 `chrom` 与 `pos` 的映射列表；`regions`：含 `chrom`、`start`、`end` 的映射列表。
- **输出：** 含 `assigned_count` 和 `positions`（每项含 `region_index`、`region_start`、`region_end`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import assign_genomic_regions

positions = [{"chrom": "chr1", "pos": 50}, {"chrom": "chr1", "pos": 1500}, {"chrom": "chr2", "pos": 10}, {"chrom": "chr3", "pos": 5}]
regions = [{"chrom": "chr1", "start": 0, "end": 1000}, {"chrom": "chr1", "start": 1000, "end": 2000}, {"chrom": "chr2", "start": 0, "end": 800}]
result = assign_genomic_regions(positions, regions, verbose=False)
print([row["region_index"] for row in result["positions"]], result["assigned_count"])
```

- **示例输出：**

```text
[0, 1, 0, None] 3
```

## 36) `dnakit.epigenome.find_read_mismatches()`

- **作用：** 找出读段中与参考序列不一致的碱基位置。
- **处理方法：** 调用 pairtools.lib.parse_pysam.get_mismatches_c，比较读段碱基与比对的参考碱基（忽略大小写）。结果按读段位置、参考位置排序。
- **API：** `dnakit.epigenome.find_read_mismatches(sequence, qualities, aligned_pairs, *, verbose=True)`。
- **输入：** `sequence`：读段序列（只含 A、C、G、T、N）；`qualities`：与碱基数相同的 PHRED 值；`aligned_pairs`：(读段位置, 参考位置, 参考碱基) 三元组。
- **输出：** 含 `mismatch_count` 和 `mismatches`（每项含 `ref_base`、`read_base`、`quality`、`ref_pos`、`read_pos`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import find_read_mismatches

pairs = [(0, 100, "A"), (1, 101, "C"), (2, 102, "G"), (3, 103, "T")]
result = find_read_mismatches("ACTT", [30, 31, 32, 33], pairs, verbose=False)
print(result["mismatch_count"], result["mismatches"][0])
```

- **示例输出：**

```text
1 {'ref_base': 'G', 'read_base': 'T', 'quality': 32, 'ref_pos': 102, 'read_pos': 2}
```

## 37) `dnakit.epigenome.parse_sam_alignments()`

- **作用：** 解析 SAM 或 BAM 文件前若干条记录，得到链方向、5' 端与 3' 端坐标。
- **处理方法：** 用 pairtools 的 AlignmentFilePairtoolized 读取，再对每条记录调用 pairtools.lib.parse.parse_pysam_entry。输入必须带 @SQ 头部。
- **API：** `dnakit.epigenome.parse_sam_alignments(sam_path, *, min_mapq=30, max_records=100, verbose=True)`。
- **输入：** `sam_path`；`min_mapq`；`max_records`（1 到 100000）。
- **输出：** 含 `record_count` 和 `alignments`（每项含 `read_id`、`chrom`、`pos5`、`pos3`、`strand`、`mapq`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import parse_sam_alignments

result = parse_sam_alignments("sample.sam", min_mapq=30, max_records=2, verbose=False)
for row in result["alignments"]:
    print(row["read_id"], row["chrom"], row["pos5"], row["pos3"], row["strand"])
```

- **示例输出：**

```text
r1 chr1 101 110 +
r1 chr2 301 310 +
```

## 38) `dnakit.epigenome.read_pairs_table()`

- **作用：** 读取 .pairs 文件：头部、染色体长度、行数和前若干行。
- **处理方法：** 调用 pairtools.lib.pairsio.read_pairs（单进程）。行数为全部数据行；`rows` 只返回前 `limit` 行。
- **API：** `dnakit.epigenome.read_pairs_table(pairs_path, *, limit=100, verbose=True)`。
- **输入：** `pairs_path`：.pairs 文件；`limit`（1 到 1000）。
- **输出：** 含 `row_count`、`columns`、`chromsizes`、`header` 与 `rows` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import read_pairs_table

result = read_pairs_table("contacts.pairs", limit=2, verbose=False)
print(result["row_count"], result["chromsizes"])
print(result["rows"][0]["readID"], result["rows"][0]["chrom1"], result["rows"][0]["pos1"])
```

- **示例输出：**

```text
40 {'chr1': 1000, 'chr2': 800}
r0 chr2 649
```

## 39) `dnakit.epigenome.make_pairs_header()`

- **作用：** 为某个基因组生成标准的 .pairs 头部行。
- **处理方法：** 调用 pairtools.lib.headerops.make_standard_pairsheader，写入版本、形状、装配名、染色体长度和列名。
- **API：** `dnakit.epigenome.make_pairs_header(*, chromsizes, assembly='unknown', columns=None, shape='upper triangle', verbose=True)`。
- **输入：** `chromsizes`（必填）；`assembly`；`columns`（可选）；`shape`（`upper triangle` 或 `whole matrix`）。
- **输出：** 含 `chromosome_count` 与 `lines` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import make_pairs_header

result = make_pairs_header(chromsizes={"chr1": 1000, "chr2": 800}, assembly="hg38", verbose=False)
for line in result["lines"][:4]:
    print(line)
```

- **示例输出：**

```text
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: hg38
#chromsize: chr1 1000
```

## 40) `dnakit.epigenome.subset_pairs_header()`

- **作用：** 只保留 .pairs 头部中指定染色体的 chromsize 与 chromosomes 行。
- **处理方法：** 调用 pairtools.lib.headerops.subset_chroms_in_pairsheader，其余行保持不变。
- **API：** `dnakit.epigenome.subset_pairs_header(header_lines, chroms, *, verbose=True)`。
- **输入：** `header_lines`：头部行列表；`chroms`：要保留的染色体名。
- **输出：** 含 `input_line_count`、`output_line_count` 和 `lines` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import make_pairs_header, subset_pairs_header

header = make_pairs_header(chromsizes={"chr1": 1000, "chr2": 800}, verbose=False)["lines"]
result = subset_pairs_header(header, ["chr2"], verbose=False)
print(result["input_line_count"], result["output_line_count"])
print([line for line in result["lines"] if line.startswith("#chromsize")])
```

- **示例输出：**

```text
6 5
['#chromsize: chr2 800']
```

## 41) `dnakit.epigenome.read_chrom_order()`

- **作用：** 按染色体文件中的顺序为每条染色体编号（从 1 开始）。
- **处理方法：** 调用 pairtools.lib.headerops.get_chrom_order。若给出 `sam_chroms`，只编号其中的染色体，其余的按名字排序后追加在末尾。
- **API：** `dnakit.epigenome.read_chrom_order(chroms_path, *, sam_chroms=None, verbose=True)`。
- **输入：** `chroms_path`：两列（名字、长度）的文件；`sam_chroms`（可选）。
- **输出：** 含 `chromosome_count` 与 `chromosomes`（每项含 `chrom`、`order`）的字典。
- **示例代码：**

```python
from dnakit.epigenome import read_chrom_order

result = read_chrom_order("chrom.sizes", sam_chroms=["chr2", "chr1", "chrC"], verbose=False)
print([(row["chrom"], row["order"]) for row in result["chromosomes"]])
```

- **示例输出：**

```text
[('chr1', 1), ('chr2', 2), ('chrC', 3)]
```

## 42) `dnakit.epigenome.assign_haplotype_phase()`

- **作用：** 根据 XA 或 XB 备选比对标签，判断一条读段属于哪一个单倍型。
- **处理方法：** 调用 pairtools.lib.phase.phase_side_XB（`variant='XB'`）或 phase_side_XA（`variant='XA'`，需给出 `mismatches`）。染色体名的后缀（如 `_1`、`_2`）决定单倍型编号。
- **API：** `dnakit.epigenome.assign_haplotype_phase(chrom, *, tag_value, primary_score, secondary_score, variant='XB', mismatches=None, phase_suffixes=('_1', '_2'), verbose=True)`。
- **输入：** `chrom`；`tag_value`：标签文本；`primary_score` 与 `secondary_score`：主、次比对分数；`variant`；`mismatches`；`phase_suffixes`。
- **输出：** 含 `phase`（`0`、`1`、`.` 或 `!`）、`chromosome_base` 与 `scores` 的字典。
- **示例代码：**

```python
from dnakit.epigenome import assign_haplotype_phase

result = assign_haplotype_phase("chr1_1", tag_value="", primary_score=50, secondary_score=40, variant="XB", verbose=False)
print(result["phase"], result["chromosome_base"], result["scores"])
```

- **示例输出：**

```text
0 chr1 [50, 40, -1]
```
