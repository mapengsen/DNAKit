# 基因组区间与注释：pybedtools、bioframe、pyranges、genomicranges、gtars 和 pyGenomeViz

这一组函数是基于可选包和外部程序的适配器：[pybedtools](https://github.com/daler/pybedtools)（配合 bedtools、samtools 和 UCSC 的 bigBed/bigWig 工具）、[bioframe](https://github.com/open2c/bioframe)、[pyranges1](https://github.com/pyranges/pyranges1)（配合 [pysam](https://github.com/pysam-developers/pysam)）、[genomicranges](https://github.com/BiocPy/genomicranges)、[gtars](https://github.com/databio/gtars) 和 [pyGenomeViz](https://github.com/moshi4/pyGenomeViz)（配合 MMseqs2、MUMmer 和 progressiveMauve）。它们可以获取染色体大小和着丝粒，在 BED、bedGraph、wiggle、bigBed 和 bigWig 之间转换，统计比对上的读段数，构建并绘制区间集合（Venn 集合、二值热图、区间轨道、itemRgb 颜色），把 BAM 比对读成区间，读取 UCSC 基因模型 GTF，把 HGVS 变异映射为 GA4GH VRS 标识符，以及绘制基因组特征图和基因组比较图。26 个区间与注释函数通过 `from dnakit.annotation import ...` 导入；6 个 pyGenomeViz 函数通过 `from dnakit.visualization import ...` 导入。

已安装的 pybedtools、bioframe、pyranges1、genomicranges、pyGenomeViz 和 pysam 的元数据声明为 MIT 许可；已安装的 gtars 包在元数据中没有声明许可。DNAKit 只调用这些包，不复制其任何代码。UCSC 工具、bedtools、samtools、MMseqs2、MUMmer 和 progressiveMauve 是需要自行安装的独立程序，各有各的许可；DNAKit 只用固定的参数列表启动它们，从不通过 shell。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时，函数会抛出带安装命令的 `BackendUnavailableError`；缺少程序时，抛出带安装提示的 `BackendUnavailableError`。输出文件必须是新文件：这些函数不会覆盖已有文件。下面有若干示例展示缺少程序时的这个错误，因为生成输出的环境里没有安装 UCSC 工具、samtools、pypairix、MMseqs2、MUMmer 和 progressiveMauve；对这些函数，另有一个标注为“安装了依赖后的用法”的代码块，给出未实际运行的真实调用。

## 示例数据

示例使用下面这些很小的文件。每个示例把输出写到临时目录。

**示例文件：** `peaks_a.bed`

```text
chr1	100	200	a1	10	+
chr1	300	400	a2	50	-
chr2	50	150	a3	90	+
```

**示例文件：** `peaks_b.bed`

```text
chr1	150	250	b1	20	+
chr1	700	800	b2	60	+
chr2	100	200	b3	70	-
```

**示例文件：** `chrom.sizes`

```text
chr1	1000
chr2	800
```

**示例文件：** `signal.bedgraph`

```text
chr1	0	100	1.5
chr1	100	200	3
```

**示例文件：** `signal.wig`

```text
variableStep chrom=chr1 span=10
101	2.0
111	3.0
```

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	101	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	16	chr1	201	30	8M	*	0	0	ACGTACGT	IIIIIIII
r3	1024	chr1	301	60	8M	*	0	0	ACGTACGT	IIIIIIII
r4	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
```

**示例文件：** `pairs.txt`

```text
r1	chr1	100	chr1	500	+	-
r2	chr1	120	chr2	300	+	+
```

**示例文件：** `genes.gtf`

```text
chr1	refGene	exon	101	200	.	+	.	gene_id "GENE1"; transcript_id "TX1";
chr1	refGene	exon	301	400	.	+	.	gene_id "GENE1"; transcript_id "TX1";
chr2	refGene	exon	51	150	.	-	.	gene_id "GENE2"; transcript_id "TX2";
```

**示例文件：** `tiny.fa`

```text
>chr1
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

**示例文件：** `genomeA.gbk`

```text
LOCUS       genomeA                   84 bp    DNA              UNK 01-JAN-1980
DEFINITION  demo genome.
ACCESSION   genomeA
VERSION     genomeA
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             11..37
                     /product="protein 10"
                     /protein_id="genomeA_10"
                     /translation="MAKGLTEV"
     CDS             complement(46..72)
                     /product="protein 45"
                     /protein_id="genomeA_45"
                     /translation="MKRLGTLD"
ORIGIN
        1 actgactgac atggctaaag gtctgaccga agtttaaact gactgttaat ccagggtacc
       61 cagacgtttc atactgactg actg
//
```

**示例文件：** `genomeB.gbk`

```text
LOCUS       genomeB                   86 bp    DNA              UNK 01-JAN-1980
DEFINITION  demo genome.
ACCESSION   genomeB
VERSION     genomeB
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             15..41
                     /product="protein 14"
                     /protein_id="genomeB_14"
                     /translation="MAKGLTEV"
     CDS             complement(48..74)
                     /product="protein 47"
                     /protein_id="genomeB_47"
                     /translation="MKRLGTLD"
ORIGIN
        1 tgcatgcatg catgatggct aaaggtctga ccgaagttta atgcatgtta atccagggta
       61 cccagacgtt tcattgcatg catgca
//
```

## 1) `dnakit.annotation.add_color()`

- **作用：** 写出一个 BED 文件，其 itemRgb 颜色按每个特征的分数经 matplotlib 颜色映射得到。
- **处理方法：** 从第 5 列读取分数，用 `matplotlib.colors.Normalize` 在最小值和最大值（或 `score_min` 与 `score_max`）之间线性缩放，再对每个特征以所选颜色映射调用 `pybedtools.featurefuncs.add_color`，用 `saveas()` 保存。返回的颜色是从新文件第 9 列读回的。
- **API：** `add_color(bed: 'str', output_bed: 'str', *, colormap: 'str' = 'viridis', score_min: 'float | None' = None, score_max: 'float | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`）。
- **输入：** `bed` 是已存在的 BED 文件，每个特征的第 5 列都必须是数值分数；`output_bed` 是已存在文件夹中的新文件；`colormap` 是 matplotlib 颜色映射名（未知名称抛出 `ConfigurationError`）；`score_min` 和 `score_max` 是可选的有限数，且 `score_min <= score_max`（默认取文件中的最小值和最大值）。最多 100,000 个特征。
- **输出：** 包含 `status`、`backend`、`output_path` 和 `data` 的字典；`data` 中有 `feature_count`、`colormap`、`score_min`、`score_max` 和 `rgb`（按文件顺序排列的 itemRgb 字符串列表）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import add_color

out = Path(tempfile.mkdtemp()) / "colored.bed"
result = add_color("peaks_a.bed", str(out), colormap="viridis", verbose=False)
data = result["data"]
print(data["feature_count"], data["score_min"], data["score_max"])
print(data["rgb"])
```

- **示例输出：**

```text
3 10.0 90.0
['68,1,84', '32,144,140', '253,231,36']
```

## 2) `dnakit.annotation.bam_to_bigwig()`

- **作用：** 从 BAM 文件写出读段覆盖度的 bigWig，可选按每百万比对读段归一化。
- **处理方法：** 运行 `BedTool(bam).genome_coverage(bg=True, split=True)`（即 pybedtools 所用的 bedtools genomecov）生成临时 bedGraph，再用 UCSC `bedGraphToBigWig` 转换。`scale=True` 时，缩放因子 `1e6 / 比对读段数` 传给 genomecov，读段数来自 `mapped_read_count()`（samtools）。如果 bedGraphToBigWig 报告顺序问题，则用 bedtools 对 bedGraph 排序后重试。临时 bedGraph 文件之后会被删除。`chromsizes` 文件只传给 bedGraphToBigWig；对 BAM 输入，bedtools 会忽略 `-g`。
- **API：** `bam_to_bigwig(bam: 'str', chromsizes: 'str', output_bigwig: 'str', *, scale: 'bool' = False, samtools_path: 'str | None' = None, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）；需要外部程序 bedtools 和 UCSC bedGraphToBigWig 在 PATH 中（或使用 `bedgraphtobigwig_path`）；`scale=True` 时还需要 samtools 在 PATH 中。
- **输入：** `bam` 是已存在且按坐标排序的 BAM；`chromsizes` 是两列的染色体大小文件；`output_bigwig` 是新路径；`scale` 是布尔值。BAM 中没有比对读段且 `scale=True` 时抛出 `ConfigurationError`。
- **输出：** 包含 `status`（`completed` 或 `failed`）、`backend`、`executable`、`command`、`return_code`、`outputs`、`log_tail`、`mapped_reads`（不缩放时为 None）、`scale_factor` 的字典；若使用了重试，还有 `sorted_before_conversion`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import bam_to_bigwig
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
try:
    bam_to_bigwig(
        str(work / "reads.bam"),
        "chrom.sizes",
        str(work / "coverage.bw"),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

BAM 由 `reads.sam` 经 pysam 生成。本环境没有安装 bedGraphToBigWig，所以示例把 `bedgraphtobigwig_path` 指向一个不存在的程序并展示错误。（对 bedtools 的检查在先，因此在没有 bedtools 的机器上，错误会提到 bedtools。）

安装了依赖后的用法：

```python
from dnakit.annotation import bam_to_bigwig

result = bam_to_bigwig("reads.bam", "chrom.sizes", "coverage.bw", scale=True)
print(result["status"], result["mapped_reads"], result["outputs"])
```

## 3) `dnakit.annotation.bedgraph_to_bigwig()`

- **作用：** 把排好序的 bedGraph 文件转换为带索引的 bigWig 文件。
- **处理方法：** 通过 DNAKit 的有界命令执行器（不经过 shell）运行 `bedGraphToBigWig <bedgraph> <chromsizes> <output>`。返回码为 0 且生成了输出文件才算成功。
- **API：** `bedgraph_to_bigwig(bedgraph: 'str', chromsizes: 'str', output_bigwig: 'str', *, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bedGraphToBigWig 在 PATH 中（或使用 `bedgraphtobigwig_path`）。不需要 Python 包。
- **输入：** `bedgraph` 是已存在、按染色体和起点排序的 bedGraph；`chromsizes` 是染色体大小文件；`output_bigwig` 是已存在文件夹中的新路径；可选 `bedgraphtobigwig_path`。
- **输出：** 包含 `status`（`completed` 或 `failed`）、`backend`、`executable`、`command`、`return_code`、`outputs`（生成时为 bigWig 路径）和 `log_tail`（程序输出的最后 4000 个字符）的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bedgraph_to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    bedgraph_to_bigwig(
        "signal.bedgraph",
        "chrom.sizes",
        str(out),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.annotation import bedgraph_to_bigwig

result = bedgraph_to_bigwig("signal.bedgraph", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 4) `dnakit.annotation.bigbed()`

- **作用：** 把 BED 文件转换为带索引的 bigBed 文件。
- **处理方法：** 运行 UCSC `bedToBigBed`，带 `-blockSize`、`-itemsPerSlot` 和 `-type`（按需加 `-unc`、`-tab` 和 `-as=`）。未给 `bed_type` 时，类型为 `bed<N>`，N 是第一条记录以制表符分隔的列数。
- **API：** `bigbed(bed: 'str', chromsizes: 'str', output_bigbed: 'str', *, bed_type: 'str | None' = None, block_size: 'int' = 256, items_per_slot: 'int' = 512, autosql: 'str | None' = None, unc: 'bool' = False, tab: 'bool' = False, bedtobigbed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bedToBigBed 在 PATH 中（或使用 `bedtobigbed_path`）。不需要 Python 包。
- **输入：** `bed` 是已存在、按染色体和起点排序的 BED 文件；`chromsizes` 是染色体大小文件；`output_bigbed` 是新路径；`bed_type` 如 `bed3`、`bed6` 或 `bed6+4`；`block_size` 和 `items_per_slot` 是正整数（不超过 1,000,000）；`autosql` 是可选的已存在 `.as` 文件；`unc` 和 `tab` 是布尔值。
- **输出：** 与 `bedgraph_to_bigwig()` 相同的字典：`status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigbed
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "peaks.bb"
try:
    bigbed(
        "peaks_a.bed",
        "chrom.sizes",
        str(out),
        bed_type="bed6",
        bedtobigbed_path="/nonexistent/bedToBigBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedToBigBed path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.annotation import bigbed

result = bigbed("peaks_a.bed", "chrom.sizes", "peaks.bb", bed_type="bed6")
print(result["status"], result["command"])
```

## 5) `dnakit.annotation.bigbed_to_bed()`

- **作用：** 把带索引的 bigBed 文件中的区间提取为 BED 文件。
- **处理方法：** 运行 UCSC `bigBedToBed`，可选 `-chrom`、`-start`、`-end` 和 `-maxItems`，后面接输入和输出路径。
- **API：** `bigbed_to_bed(bigbed: 'str', output_bed: 'str', *, chrom: 'str | None' = None, start: 'int | None' = None, end: 'int | None' = None, max_items: 'int | None' = None, bigbedtobed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bigBedToBed 在 PATH 中（或使用 `bigbedtobed_path`）。不需要 Python 包。
- **输入：** `bigbed` 是已存在的 bigBed 文件；`output_bed` 是新路径；`chrom`、`start` 和 `end` 用于限定区域（给了 `start` 或 `end` 却没有 `chrom` 时抛出 `ConfigurationError`）；`max_items` 是正整数。
- **输出：** 包含 `status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigbed_to_bed
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
(work / "peaks.bb").write_bytes(b"")
try:
    bigbed_to_bed(
        str(work / "peaks.bb"),
        str(work / "peaks.bed"),
        chrom="chr1",
        bigbedtobed_path="/nonexistent/bigBedToBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bigBedToBed path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。示例用一个空的占位文件代替 bigBed，因为程序是在打开文件之前查找的。

安装了依赖后的用法：

```python
from dnakit.annotation import bigbed_to_bed

result = bigbed_to_bed("peaks.bb", "peaks_chr1.bed", chrom="chr1", start=0, end=500)
print(result["status"], result["outputs"])
```

## 6) `dnakit.annotation.bigwig_to_bedgraph()`

- **作用：** 把带索引的 bigWig 文件中的数值提取为 bedGraph 文件。
- **处理方法：** 运行 UCSC `bigWigToBedGraph`，可选 `-chrom`、`-start`、`-end` 和 `-udcDir`，后面接输入和输出路径。
- **API：** `bigwig_to_bedgraph(bigwig: 'str', output_bedgraph: 'str', *, chrom: 'str | None' = None, start: 'int | None' = None, end: 'int | None' = None, udc_dir: 'str | None' = None, bigwigtobedgraph_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bigWigToBedGraph 在 PATH 中（或使用 `bigwigtobedgraph_path`）。不需要 Python 包。
- **输入：** `bigwig` 是已存在的 bigWig 文件；`output_bedgraph` 是新路径；`chrom`、`start` 和 `end` 用于限定区域（给了 `start` 或 `end` 却没有 `chrom` 时抛出 `ConfigurationError`）；`udc_dir` 是传给 `-udcDir` 的缓存文件夹名。
- **输出：** 包含 `status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import bigwig_to_bedgraph
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
(work / "signal.bw").write_bytes(b"")
try:
    bigwig_to_bedgraph(
        str(work / "signal.bw"),
        str(work / "signal.bedgraph"),
        bigwigtobedgraph_path="/nonexistent/bigWigToBedGraph",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bigWigToBedGraph path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。示例用一个空的占位文件代替 bigWig，因为程序是在打开文件之前查找的。

安装了依赖后的用法：

```python
from dnakit.annotation import bigwig_to_bedgraph

result = bigwig_to_bedgraph("signal.bw", "signal_chr1.bedgraph", chrom="chr1", start=0, end=300)
print(result["status"], result["outputs"])
```

## 7) `dnakit.annotation.binary_heatmap()`

- **作用：** 绘制二值热图，显示每个区间出现在哪些 BED 文件中，并返回各类别的计数。
- **处理方法：** 运行 `BedTool().multi_intersect(i=beds, names=names, cluster=True)`（bedtools multiinter），统计每一行的类别标签（如 `a` 或 `a,b`），由各文件所在的列构成 0/1 存在矩阵，用与 `sort_binary_matrix()` 相同的列加权分数给行排序，再用 matplotlib `imshow`（颜色映射 `binary`）以 120 dpi 画成 PNG。
- **API：** `binary_heatmap(beds: 'Sequence[str]', names: 'Sequence[str]', output_png: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`），以及外部程序 bedtools 在 PATH 中。
- **输入：** `beds` 是 2 到 30 个已存在的 BED 文件；`names` 为每个文件一个唯一名称（字母、数字、`_`、`.`、`-` 和空格，最多 40 个字符）；`output_png` 是新路径。如果这些文件没有共享区间，抛出 `ConfigurationError`。
- **输出：** 包含 `status`、`backend`、`image_path` 和 `data` 的字典；`data` 中有 `names`、`interval_count`、`class_counts`（先按计数、再按名称排序）、`matrix` 和 `row_order`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import binary_heatmap

out = Path(tempfile.mkdtemp()) / "heatmap.png"
result = binary_heatmap(["peaks_a.bed", "peaks_b.bed"], ["a", "b"], str(out), verbose=False)
data = result["data"]
print(data["interval_count"], data["class_counts"])
print(data["matrix"], data["row_order"])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **示例输出：**

```text
4 {'a,b': 2, 'a': 1, 'b': 1}
[[1, 1], [1, 0], [0, 1], [1, 1]] [3, 0, 1, 2]
True
```

## 8) `dnakit.annotation.chromsizes()`

- **作用：** 返回某个基因组组装的染色体大小。
- **处理方法：** 调用 `pybedtools.chromsizes(genome)`，它先查 pybedtools 内置的组装注册表，否则依次使用 genomepy 和 UCSC MySQL 服务器。`service` 字段说明走了哪条路径：注册表是离线的，其他名称需要 genomepy 或网络。大小以“名称到长度”的映射返回；给了 `output_path` 时还会写成染色体大小文件。
- **API：** `chromsizes(genome: 'str', *, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）。内置注册表中的组装（如 `hg19`）可离线使用；其他组装需要 genomepy 或访问 UCSC 的网络。
- **输入：** `genome` 是由字母、数字、`_`、`.` 和 `-` 组成的组装名（最多 64 个字符）；`output_path` 是可选的新文件路径。查询失败时抛出 `BackendExecutionError`。
- **输出：** 包含 `status`、`service`、`request`、`retrieved_at`（UTC 时间）、`pybedtools_version`、`row_count`、`data`（染色体名到长度）和 `output_path` 的字典。
- **示例代码：**

```python
from dnakit.annotation import chromsizes

result = chromsizes("hg19", verbose=False)
print(result["status"], result["service"])
print(result["row_count"], result["data"]["chr1"], result["data"]["chrM"])
```

- **示例输出：**

```text
completed pybedtools built-in genome registry (offline)
93 249250621 16571
```

## 9) `dnakit.annotation.cleaned_intersect()`

- **作用：** 为 2 到 4 个 BED 文件构建区间字符串集合，其中相互重叠的区间被统一为相同的区间。
- **处理方法：** 调用 `pybedtools.contrib.venn_maker.cleaned_intersect()`，它对每一对重叠的文件做交集，使不同文件中相互重叠的区间变成同一个区间。随后把每个区间转成字符串（制表符换成 `|`），每个集合排序并去重。
- **API：** `cleaned_intersect(beds: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）和外部程序 bedtools 在 PATH 中。
- **输入：** `beds` 是 2 到 4 个已存在的 BED 文件的列表。
- **输出：** 包含 `status`、`backend`、`pybedtools_version`、`files`、`sets`（每个文件一个排好序的字符串列表）和 `set_sizes` 的字典。
- **示例代码：**

```python
from dnakit.annotation import cleaned_intersect

result = cleaned_intersect(["peaks_a.bed", "peaks_b.bed"], verbose=False)
print(result["set_sizes"])
for interval_set in result["sets"]:
    print(interval_set)
```

- **示例输出：**

```text
[3, 3]
['chr1|100|200', 'chr1|300|400', 'chr2|50|150']
['chr1|100|200', 'chr1|700|800', 'chr2|50|150']
```

## 10) `dnakit.annotation.get_chromsizes_from_genomepy()`

- **作用：** 读取已用 genomepy 安装的基因组的染色体大小。
- **处理方法：** 调用 `pybedtools.helpers.get_chromsizes_from_genomepy(genome)`。当 genomepy 没有该基因组的大小文件时，函数不会报错，而是返回 `status` 为 `not_found` 的结果，带一条消息，且 `data` 为 None。
- **API：** `get_chromsizes_from_genomepy(genome: 'str', *, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 genomepy（额外依赖 `dnakit[bio-genomepy]`）和 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）。基因组必须已用 genomepy 安装；该函数自己不下载任何东西。
- **输入：** `genome` 是组装名；`output_path` 是可选的新文件路径（只在找到大小时才写入）。genomepy 的其他失败会抛出 `BackendExecutionError`。
- **输出：** 包含 `status`（`completed` 或 `not_found`）、`service`、`request`、`retrieved_at`、`pybedtools_version`、`row_count`、`data`（染色体名到长度，或 None）、`output_path` 的字典；`not_found` 时还有 `message`。
- **示例代码：**

```python
from dnakit.annotation import get_chromsizes_from_genomepy

result = get_chromsizes_from_genomepy("no_such_genome_demo", verbose=False)
print(result["status"], result["row_count"], result["data"])
print(result["message"])
```

- **示例输出：**

```text
not_found 0 None
genomepy has no local sizes file for this genome; install it with genomepy install.
```

本环境没有安装这个名称的基因组，所以结果是 `not_found`。

安装了依赖后的用法：

```python
from dnakit.annotation import get_chromsizes_from_genomepy

# after: genomepy install hg38
result = get_chromsizes_from_genomepy("hg38", output_path="hg38.chrom.sizes")
print(result["status"], result["row_count"])
```

## 11) `dnakit.annotation.get_chromsizes_from_ucsc()`

- **作用：** 从 UCSC MySQL 服务器下载 UCSC 组装的染色体大小。
- **处理方法：** 调用 `pybedtools.get_chromsizes_from_ucsc(genome, saveas=..., timeout=..., host_url=...)`，它针对 UCSC 服务器运行 `mysql` 和 `fetchChromSizes` 程序。`OSError` 和 `ValueError` 会被重新抛出为 `BackendExecutionError`。
- **API：** `get_chromsizes_from_ucsc(genome: 'str', *, host_url: 'str' = 'genome-mysql.cse.ucsc.edu', timeout: 'float | None' = None, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）、能访问 UCSC MySQL 服务器的网络（不需要密钥），以及 PATH 中的 `mysql` 和 `fetchChromSizes` 程序（由 pybedtools 在查询时启动）。
- **输入：** `genome` 是组装名，如 `hg38`（字母、数字、`_`、`.` 和 `-`）；`host_url` 是主机名；`timeout` 是 (0, 3600] 内的秒数或 None；`output_path` 是可选的新文件路径。
- **输出：** 包含 `status`、`service`、`request`（含 `genome` 和 `host_url`）、`retrieved_at`、`pybedtools_version`、`row_count`、`data` 和 `output_path` 的字典。
- **示例代码：**

```python
from dnakit.annotation import get_chromsizes_from_ucsc
from dnakit.exceptions import ConfigurationError

try:
    get_chromsizes_from_ucsc("hg38; drop", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "genome must be an assembly name" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

示例不访问网络，只展示无效的组装名如何在发出任何请求之前被拒绝。

安装了依赖后的用法：

```python
from dnakit.annotation import get_chromsizes_from_ucsc

result = get_chromsizes_from_ucsc("hg38", timeout=60, output_path="hg38.chrom.sizes")
print(result["row_count"], result["data"]["chr1"])
```

## 12) `dnakit.annotation.mapped_read_count()`

- **作用：** 统计 BAM 文件中比对上的读段数。
- **处理方法：** 运行 `samtools view -c -F 0x4 <bam>`，把输出的最后一行读成整数。pybedtools 会把计数缓存在 BAM 旁边的 `<bam>.scale` 文件里，而这个函数不读也不写该缓存，所以输入文件夹保持不变。
- **API：** `mapped_read_count(bam: 'str', *, samtools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 samtools 在 PATH 中（或使用 `samtools_path`）。不需要 Python 包。
- **输入：** `bam` 是已存在的 BAM 文件；可选 `samtools_path`。
- **输出：** 包含 `status`（`completed` 或 `failed`）、`backend`、`executable`、`command`、`return_code`、`mapped_reads`（整数；无法读出计数时为 None）、`outputs`（为空）和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import mapped_read_count
from dnakit.exceptions import BackendUnavailableError

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
try:
    mapped_read_count(str(work / "reads.bam"), samtools_path="/nonexistent/samtools", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "samtools path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装 samtools，所以示例展示程序路径不存在时的错误。BAM 由 `reads.sam` 经 pysam 生成。

安装了依赖后的用法：

```python
from dnakit.annotation import mapped_read_count

result = mapped_read_count("reads.bam")
print(result["status"], result["mapped_reads"])
```

## 13) `dnakit.annotation.sort_binary_matrix()`

- **作用：** 返回 0/1 存在矩阵按列加权分数从高到低排序后的行顺序。
- **处理方法：** 设矩阵有 C 列，第 j 列（从左起，从 1 计数）的权重为 `2 ** (C - j + 1)`，所以最左列权重最大。把每一行的加权值求和，用 `numpy.argsort(...)[::-1]` 返回行下标，与 pybedtools 的 `sort_binary_matrix` 相同。空矩阵返回空列表。
- **API：** `sort_binary_matrix(matrix: 'Sequence[Sequence[int]]', *, verbose: bool = True) -> 'list[int]'`.
- **运行要求：** 需要 numpy（随 `dnakit[bio-pybedtools]` 环境安装；缺少 numpy 时函数会报告该额外依赖）。
- **输入：** `matrix` 是等长、非空行组成的列表，只能含整数 0 和 1（布尔值会被拒绝）；最多 100,000 行。其他情况抛出 `ConfigurationError`。
- **输出：** 行下标列表（Python 整数）。
- **示例代码：**

```python
from dnakit.annotation import sort_binary_matrix

matrix = [
    [0, 1, 1],
    [1, 0, 0],
    [1, 1, 0],
    [0, 0, 1],
]
print(sort_binary_matrix(matrix, verbose=False))
```

- **示例输出：**

```text
[2, 1, 0, 3]
```

## 14) `dnakit.annotation.venn_maker()`

- **作用：** 写出一个用 VennDiagram 包绘制 2 到 4 个 BED 文件 Venn 图的 R 脚本，并在有 Rscript 时渲染。
- **处理方法：** 用 `cleaned_intersect()` 构建区间集合，调用 `pybedtools.contrib.venn_maker.venn_maker(..., run=False)` 写出 R 脚本，并统计每个精确集合组合（如 `a` 或 `a&b`）中的条目数。如果 `render` 为真、给了 `figure_path` 且 PATH 中有 `Rscript`，就运行脚本；否则只写出脚本，并由 `render_note` 说明原因。
- **API：** `venn_maker(beds: 'Sequence[str]', output_script: 'str', *, names: 'Sequence[str] | None' = None, figure_path: 'str | None' = None, additional_args: 'Sequence[str] | None' = None, render: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pybedtools（额外依赖 `dnakit[bio-pybedtools]`）和外部程序 bedtools 在 PATH 中。渲染图片还需要装有 VennDiagram 包的 R，并且 `Rscript` 在 PATH 中；只写脚本和计数则不需要。
- **输入：** `beds` 是 2 到 4 个已存在的 BED 文件；`output_script` 是 R 脚本的新路径；`names` 是唯一的标签（默认 `a`、`b`、`c`、`d`）；`figure_path` 是可选的新路径，不能含引号或反斜杠；`additional_args` 是 R 的命名参数，如 `euler.d=TRUE`（不能含分号或换行）；`render` 是布尔值。
- **输出：** 包含 `status`、`backend`、`script_path`、`figure_path`、`rendered`、`render_note`、`render`（Rscript 的返回码和日志尾部，或 None）以及 `data`（含 `names`、`set_sizes` 和 `regions`）的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import venn_maker

script = Path(tempfile.mkdtemp()) / "venn.R"
result = venn_maker(
    ["peaks_a.bed", "peaks_b.bed"],
    str(script),
    names=["A", "B"],
    render=False,
    verbose=False,
)
print(result["data"])
print(result["rendered"], result["render_note"])
print(script.is_file(), "venn.diagram" in script.read_text())
```

- **示例输出：**

```text
{'names': ['A', 'B'], 'set_sizes': {'A': 3, 'B': 3}, 'regions': {'A': 1, 'A&B': 2, 'B': 1}}
False Rendering was not requested; the R script was written only.
True True
```

示例使用 `render=False`，所以有没有 R 输出都相同。

安装了依赖后的用法：

```python
from dnakit.annotation import venn_maker

result = venn_maker(["peaks_a.bed", "peaks_b.bed"], "venn.R", names=["A", "B"], figure_path="venn.png")
print(result["rendered"], result["render_note"])
```

## 15) `dnakit.annotation.wig_to_bigwig()`

- **作用：** 把 wiggle 文件转换为带索引的 bigWig 文件。
- **处理方法：** 通过 DNAKit 的有界命令执行器（不经过 shell）运行 UCSC `wigToBigWig <wig> <chromsizes> <output>`。
- **API：** `wig_to_bigwig(wig: 'str', chromsizes: 'str', output_bigwig: 'str', *, wigtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC wigToBigWig 在 PATH 中（或使用 `wigtobigwig_path`）。不需要 Python 包。
- **输入：** `wig` 是已存在的 wiggle 文件；`chromsizes` 是染色体大小文件；`output_bigwig` 是新路径；可选 `wigtobigwig_path`。
- **输出：** 包含 `status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import wig_to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    wig_to_bigwig(
        "signal.wig",
        "chrom.sizes",
        str(out),
        wigtobigwig_path="/nonexistent/wigToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "wigToBigWig path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.annotation import wig_to_bigwig

result = wig_to_bigwig("signal.wig", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 16) `dnakit.annotation.fetch_centromeres()`

- **作用：** 获取某个组装每条染色体的着丝粒位置。
- **处理方法：** 调用 `bioframe.fetch_centromeres(db, provider=provider)`。`provider="ucsc"` 时表来自 UCSC 基因组浏览器；`provider="local"` 时使用 bioframe 自带的表，不需要网络。行以记录列表返回；给了 `output_path` 时还会把表写成带表头、以制表符分隔的文件。
- **API：** `fetch_centromeres(db: 'str', *, provider: 'str' = 'ucsc', output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 bioframe（额外依赖 `dnakit[bio-bioframe]`）。`provider="local"` 对 bioframe 自带的组装可离线使用；`provider="ucsc"`（默认）需要访问网络。
- **输入：** `db` 是组装名，如 `hg38`；`provider` 是 `ucsc` 或 `local`；`output_path` 是可选的新文件路径。查询失败（bioframe 抛出 `OSError`、`ValueError` 或 `KeyError`）时抛出 `BackendExecutionError`。
- **输出：** 包含 `status`、`service`、`request`、`retrieved_at`、`bioframe_version`、`row_count`、`data`（bioframe 返回的各列组成的记录列表，如 `chrom`、`start`、`end`、`mid`）和 `output_path` 的字典。
- **示例代码：**

```python
from dnakit.annotation import fetch_centromeres

result = fetch_centromeres("hg38", provider="local", verbose=False)
print(result["service"])
print(result["row_count"], result["data"][0])
```

- **示例输出：**

```text
bioframe bundled assembly tables (no network access)
24 {'chrom': 'chr1', 'start': 121700000, 'end': 125100000, 'mid': 123400000}
```

## 17) `dnakit.annotation.fetch_chromsizes()`

- **作用：** 以“名称到长度”的映射或 BED3 表的形式获取某个组装的染色体长度。
- **处理方法：** 调用 `bioframe.fetch_chromsizes(db, provider=..., as_bed=..., filter_chroms=...)`。`filter_chroms=True` 只保留主要染色体。给了 `output_path` 时，表不带表头写出（两列大小文件，`as_bed` 为真时是 BED3）。
- **API：** `fetch_chromsizes(db: 'str', *, provider: 'str' = 'ucsc', as_bed: 'bool' = False, filter_chroms: 'bool' = True, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 bioframe（额外依赖 `dnakit[bio-bioframe]`）。`provider="local"` 对 bioframe 自带的组装可离线使用；`provider="ucsc"`（默认）需要访问网络。
- **输入：** `db` 是组装名；`provider` 是 `ucsc` 或 `local`；`as_bed` 和 `filter_chroms` 是布尔值；`output_path` 是可选的新文件路径。查询失败时抛出 `BackendExecutionError`。
- **输出：** 包含 `status`、`service`、`request`、`retrieved_at`、`bioframe_version`、`row_count`、`data`（染色体到长度的映射，`as_bed` 为真时是 `chrom`/`start`/`end` 记录列表）和 `output_path` 的字典。
- **示例代码：**

```python
from dnakit.annotation import fetch_chromsizes

sizes = fetch_chromsizes("hg38", provider="local", verbose=False)
print(sizes["row_count"], sizes["data"]["chr1"], sizes["data"]["chrX"])
bed = fetch_chromsizes("hg38", provider="local", as_bed=True, verbose=False)
print(bed["data"][0])
```

- **示例输出：**

```text
25 248956422 156040895
{'chrom': 'chr1', 'start': 0, 'end': 248956422}
```

## 18) `dnakit.annotation.plot_intervals()`

- **作用：** 把类 BED 的区间文件画成 PNG，每条染色体一个堆叠的轨道。
- **处理方法：** 读取没有表头、以制表符分隔的文件（第 0 到 2 列是染色体、起点和终点），按染色体分组，用 `bioframe.core.arrops.stack_intervals` 给每个区间分配一行，使重叠的区间不会互相遮挡，再在 matplotlib 的 Figure（而不是 pyplot）上把每个区间画成矩形，可选用某一列作文字标签。
- **API：** `plot_intervals(intervals: 'str', output_png: 'str', *, label_column: 'int | None' = None, color: 'str' = 'skyblue', show_coords: 'bool' = False, figsize: 'Sequence[float]' = (10.0, 2.0), verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 bioframe（额外依赖 `dnakit[bio-bioframe]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`），以及 bioframe 一并安装的 pandas 和 numpy。
- **输入：** `intervals` 是已存在、以制表符分隔、至少三列的文件，起点和终点为整数，满足 `0 <= start <= end`，共 1 到 100,000 行；`output_png` 是新路径；`label_column` 是列下标或 None；`color` 是 matplotlib 颜色名；`show_coords` 是布尔值；`figsize` 是两个正数（宽度、每个轨道的高度）。
- **输出：** 包含 `status`、`backend`、`image_path` 和 `data` 的字典；`data` 中有 `source`、`interval_count`、`chromosomes`、`tracks`（染色体、区间数和最高堆叠层级）和 `intervals`（每个区间的染色体、起点、终点、层级和标签）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import plot_intervals

out = Path(tempfile.mkdtemp()) / "intervals.png"
result = plot_intervals("peaks_a.bed", str(out), label_column=3, verbose=False)
data = result["data"]
print(data["interval_count"], data["chromosomes"], data["tracks"])
print(data["intervals"][0])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **示例输出：**

```text
3 ['chr1', 'chr2'] [{'chrom': 'chr1', 'interval_count': 2, 'max_level': 0}, {'chrom': 'chr2', 'interval_count': 1, 'max_level': 0}]
{'chrom': 'chr1', 'start': 100, 'end': 200, 'level': 0, 'label': 'a1'}
True
```

## 19) `dnakit.annotation.plot_intervals_arr()`

- **作用：** 把用起点和终点数组给出的一组区间画成 PNG。
- **处理方法：** 使用你给出的堆叠层级，或用 `bioframe.core.arrops.stack_intervals` 计算，然后在 matplotlib 的 Figure 上为每个区间画一个矩形，可选在中央加标签。除非给了 `xlim`，x 范围是 `-0.5` 到 `max(ends) + 0.5`。
- **API：** `plot_intervals_arr(starts: 'Sequence[float]', ends: 'Sequence[float]', output_png: 'str', *, levels: 'Sequence[float] | None' = None, labels: 'Sequence[str | None] | None' = None, colors: 'str | Sequence[str] | None' = None, xlim: 'Sequence[float] | None' = None, show_coords: 'bool' = False, figsize: 'Sequence[float]' = (10.0, 2.0), verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 bioframe（额外依赖 `dnakit[bio-bioframe]`）和 matplotlib（额外依赖 `dnakit[bio-plots]`），以及 bioframe 一并安装的 numpy。
- **输入：** `starts` 和 `ends` 是等长、非空的有限数列表，且 `start <= end`；`output_png` 是新路径；可选 `levels`（每个区间一个）、`labels`（每个区间一个字符串或 None）、`colors`（一个颜色名或每个区间一个）、`xlim`（两个数，小的在前）、`show_coords` 和 `figsize`。
- **输出：** 包含 `status`、`backend`、`image_path` 和 `data` 的字典；`data` 中有 `interval_count`、`starts`、`ends`、`levels`、`labels`、`colors`、`xlim` 和 `ylim`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import plot_intervals_arr

out = Path(tempfile.mkdtemp()) / "arr.png"
result = plot_intervals_arr(
    [0, 10, 25],
    [20, 30, 40],
    str(out),
    labels=["x", "y", None],
    colors="salmon",
    verbose=False,
)
data = result["data"]
print(data["levels"], data["labels"], data["xlim"], data["ylim"])
print(out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **示例输出：**

```text
[0, 1, 0] ['x', 'y', None] [-0.5, 40.5] [-0.5, 1.5]
True
```

## 20) `dnakit.annotation.read_pairix()`

- **作用：** 读取 pairix 索引文件中落在某个区域内（或两个区域之间）的接触对。
- **处理方法：** 调用 `bioframe.read_pairix(path, region1, region2=..., chromsizes=..., usecols=...)`，它用 pypairix 查询索引，并把表以记录列表返回。
- **API：** `read_pairix(path: 'str', region1: 'str', *, region2: 'str | None' = None, chromsizes: 'Mapping[str, int] | None' = None, usecols: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pypairix 和 cytoolz（额外依赖 `dnakit[bio-pairix]`）以及 bioframe（额外依赖 `dnakit[bio-bioframe]`）。文件必须是 bgzip 压缩并已用 pairix 建立索引的。
- **输入：** `path` 是已存在的带索引接触对文件；`region1` 是区域字符串，如 `chr1:1-1000`；`region2` 是可选的第二个区域；`chromsizes` 是可选的“染色体名到正整数长度”的映射；`usecols` 是可选的非空列名列表。pypairix 和读取失败会抛出 `BackendExecutionError`。
- **输出：** 包含 `status`、`backend`、`bioframe_version`、`request`、`row_count`、`columns` 和 `data`（记录列表）的字典。
- **示例代码：**

```python
from dnakit.annotation import read_pairix
from dnakit.exceptions import BackendUnavailableError

try:
    read_pairix("pairs.txt", "chr1:1-1000", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "pypairix" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装 pypairix，所以示例展示 `BackendUnavailableError`。普通文本文件代替带索引的文件，因为对 pypairix 的检查发生在读取文件之前。

安装了依赖后的用法：

```python
from dnakit.annotation import read_pairix

# pairs.sorted.gz: bgzip-compressed and indexed with `pairix -p pairs`
result = read_pairix("pairs.sorted.gz", "chr1:1-1000", region2="chr2:1-1000")
print(result["row_count"], result["columns"])
```

## 21) `dnakit.annotation.to_bigbed()`

- **作用：** 用显式给出的染色体大小文件，把排好序的 BED 文件写成带索引的 bigBed 文件。
- **处理方法：** 通过 DNAKit 的有界命令执行器（不经过 shell）运行 UCSC `bedToBigBed [-type=<bed_format>] <bed> <chromsizes> <output>`。与 `bigbed()` 不同，它不传块大小或 autoSql 选项，也不猜测类型。
- **API：** `to_bigbed(bed: 'str', chromsizes: 'str', output_bigbed: 'str', *, bed_format: 'str | None' = None, bedtobigbed_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bedToBigBed 在 PATH 中（或使用 `bedtobigbed_path`）。不需要 Python 包。
- **输入：** `bed` 是已存在、按染色体和起点排序的 BED 文件；`chromsizes` 是染色体大小文件；`output_bigbed` 是新路径；`bed_format` 是可选类型，如 `bed3`、`bed6` 或 `bed6+4`。
- **输出：** 包含 `status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import to_bigbed
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "peaks.bb"
try:
    to_bigbed(
        "peaks_a.bed",
        "chrom.sizes",
        str(out),
        bed_format="bed6",
        bedtobigbed_path="/nonexistent/bedToBigBed",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedToBigBed path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.annotation import to_bigbed

result = to_bigbed("peaks_a.bed", "chrom.sizes", "peaks.bb", bed_format="bed6")
print(result["status"], result["command"])
```

## 22) `dnakit.annotation.to_bigwig()`

- **作用：** 用显式给出的染色体大小文件，把排好序的 bedGraph 文件写成带索引的 bigWig 文件。
- **处理方法：** 通过 DNAKit 的有界命令执行器（不经过 shell）运行 UCSC `bedGraphToBigWig <bedgraph> <chromsizes> <output>`。它与 `bedgraph_to_bigwig()` 做的事情相同，是后者在 bioframe 命名下的对应函数。
- **API：** `to_bigwig(bedgraph: 'str', chromsizes: 'str', output_bigwig: 'str', *, bedgraphtobigwig_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要外部程序 UCSC bedGraphToBigWig 在 PATH 中（或使用 `bedgraphtobigwig_path`）。不需要 Python 包。
- **输入：** `bedgraph` 是已存在、按染色体和起点排序的 bedGraph；`chromsizes` 是染色体大小文件；`output_bigwig` 是新路径；可选 `bedgraphtobigwig_path`。
- **输出：** 包含 `status`、`backend`、`executable`、`command`、`return_code`、`outputs` 和 `log_tail` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import to_bigwig
from dnakit.exceptions import BackendUnavailableError

out = Path(tempfile.mkdtemp()) / "signal.bw"
try:
    to_bigwig(
        "signal.bedgraph",
        "chrom.sizes",
        str(out),
        bedgraphtobigwig_path="/nonexistent/bedGraphToBigWig",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "bedGraphToBigWig path does not exist" in str(error))
```

- **示例输出：**

```text
BackendUnavailableError True
```

本环境没有安装该程序，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.annotation import to_bigwig

result = to_bigwig("signal.bedgraph", "chrom.sizes", "signal.bw")
print(result["status"], result["outputs"])
```

## 23) `dnakit.annotation.read_pyranges1_bam()`

- **作用：** 把 BAM 比对读成 pyranges1 区间，列为 `Chromosome`、`Start`、`End`、`Strand` 和 `Flag`。
- **处理方法：** 用 pysam 读取 BAM（`pyranges1.read_bam` 需要另外的 bamread 包，所以不调用它），应用过滤条件，`Start` 取 0 起始的参考起点，`End` 取 CIGAR 的参考跨度（终点开区间），标志 16 时 `Strand` 为 `-`，并用该表构建 `pyranges1.PyRanges`。返回的记录是前 `max_rows` 个区间。
- **API：** `read_pyranges1_bam(path: 'str | Path', *, mapq: 'int' = 0, required_flag: 'int' = 0, filter_flag: 'int' = 1540, sparse: 'bool' = True, max_rows: 'int' = 1000, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pysam（额外依赖 `dnakit[bio-reads]`）、pyranges1 和 pandas（额外依赖 `dnakit[bio-pyranges1]`）。
- **输入：** `path` 是已存在的 BAM 文件；`mapq`（0 到 255）保留比对质量不低于该值的读段；`required_flag` 保留带有所有这些 SAM 标志位的读段；`filter_flag` 丢弃带有其中任一标志位的读段（默认 1540：未比对、质控失败、重复）；`sparse` 是布尔值（`False` 时增加 `QueryStart`、`QueryEnd`、`QuerySequence`、`Name`、`Cigar` 和 `Quality`）；`max_rows` 为 1 到 1,000,000。整个文件都会被读取；没有参考跨度的读段被跳过。
- **输出：** 包含 `status`、`operation`、`provenance`（包、后端和坐标约定）、`columns`、`record_count`（保留的全部区间数）、`truncated`（保留的区间多于 `max_rows` 时为真）和 `records` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.annotation import read_pyranges1_bam

work = Path(tempfile.mkdtemp())
pysam.view("-b", "-o", str(work / "reads.bam"), "reads.sam", catch_stdout=False)
result = read_pyranges1_bam(str(work / "reads.bam"), verbose=False)
print(result["columns"], result["record_count"], result["truncated"])
for record in result["records"]:
    print(record)
strict = read_pyranges1_bam(str(work / "reads.bam"), mapq=50, verbose=False)
print(strict["record_count"])
```

- **示例输出：**

```text
['Chromosome', 'Start', 'End', 'Strand', 'Flag'] 2 False
{'Chromosome': 'chr1', 'Start': 100, 'End': 110, 'Strand': '+', 'Flag': 0}
{'Chromosome': 'chr1', 'Start': 200, 'End': 208, 'Strand': '-', 'Flag': 16}
1
```

BAM 由 `reads.sam` 经 pysam 生成。读段 `r3` 是重复，`r4` 未比对，所以默认过滤会丢弃这两条。

## 24) `dnakit.annotation.read_ucsc_gtf()`

- **作用：** 把 UCSC 基因模型 GTF 读成表，并写成 TSV 文件。
- **处理方法：** 用 `genomicranges.io.gtf.parse_gtf` 解析 GTF，用 `to_csv(sep="\t", index=False)` 写出表并汇总。不给 `local_gtf` 时，从 `ucsc_gtf_url()` 给出的地址下载；给了 `local_gtf`（普通文件或 `.gz`）时读取本地文件，不下载任何东西，`genome` 和 `annotation` 只用来标注结果。`source_date` 对本地文件是文件的修改日期，对下载是当天日期。
- **API：** `read_ucsc_gtf(genome: 'str', output_table: 'str', *, annotation: 'str' = 'refGene', local_gtf: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 genomicranges 和 pandas（额外依赖 `dnakit[bio-genomicranges]`）。不给 `local_gtf` 时还需要访问 UCSC 的网络。
- **输入：** `genome` 是 UCSC 基因组代码，如 `hg38`（不含 `/`，不以 `-` 开头）；`output_table` 是已存在文件夹中的新文件；`annotation` 是 `refGene`、`ensGene`、`knownGene` 和 `ncbiRefSeq` 之一；`local_gtf` 是可选的已存在文件。
- **输出：** 包含 `status`、`backend`、`genome`、`annotation`、`source`、`source_date`、`row_count`、`columns`、`rows_per_chromosome` 和 `output_table` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation import read_ucsc_gtf

out = Path(tempfile.mkdtemp()) / "genes.tsv"
result = read_ucsc_gtf("hg38", str(out), local_gtf="genes.gtf", verbose=False)
print(result["row_count"], result["rows_per_chromosome"])
print(result["columns"])
print(len(out.read_text().splitlines()))
```

- **示例输出：**

```text
3 {'chr1': 2, 'chr2': 1}
['seqnames', 'source', 'feature', 'starts', 'ends', 'score', 'strand', 'frame', 'gene_id', 'transcript_id']
4
```

示例读取本地的 `genes.gtf`，所以不会下载任何东西。

## 25) `dnakit.annotation.ucsc_gtf_url()`

- **作用：** 返回某个基因组的基因模型 GTF 在 UCSC 的下载地址。
- **处理方法：** 拼出字符串 `http://hgdownload.cse.ucsc.edu/goldenPath/<genome>/bigZips/genes/<genome>.<annotation>.gtf.gz`。genomicranges 自己构造的地址在文件名前有重复的斜杠，而这个函数只写一个斜杠。不发出任何请求。
- **API：** `ucsc_gtf_url(genome: 'str', annotation: 'str' = 'refGene', *, verbose: bool = True) -> 'str'`.
- **运行要求：** 直接 pip 安装后即可使用（它只拼接字符串，不导入任何可选包）。
- **输入：** `genome` 是 UCSC 基因组代码（非空，不含 `/`，不以 `-` 开头）；`annotation` 是 `refGene`、`ensGene`、`knownGene` 和 `ncbiRefSeq` 之一。其他值抛出 `ConfigurationError`。
- **输出：** 字符串形式的地址。
- **示例代码：**

```python
from dnakit.annotation import ucsc_gtf_url

print(ucsc_gtf_url("hg38", verbose=False))
print(ucsc_gtf_url("mm10", "ncbiRefSeq", verbose=False))
```

- **示例输出：**

```text
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
```

## 26) `dnakit.annotation.hgvs_to_vrs_allele_id()`

- **作用：** 把 HGVS 变异映射为它的 GA4GH VRS 等位基因标识符。
- **处理方法：** 打开转录本存储（`gtars.reftx.ReadonlyTxStore.open`）并包装成 `ReftxProvider`，用 `gtars.refget.RefgetStore.open_local` 打开 refget 存储目录，再调用 `gtars.vrs.hgvs.hgvs_to_vrs_id`。参考不匹配（HGVS 中的参考碱基与序列不同）或其他映射失败，都会抛出 `ConfigurationError`；不会下载任何东西。
- **API：** `hgvs_to_vrs_allele_id(hgvs: 'str', transcript_store: 'str', refget_store: 'str', collection_digest: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 gtars（额外依赖 `dnakit[bio-gtars]`），以及你事先一次性构建好的两个本地输入：转录本存储文件和 refget 存储目录。
- **输入：** `hgvs` 是 HGVS 字符串，如 `NM_TEST.1:c.4G>A`；`transcript_store` 是已存在的存储文件；`refget_store` 是已存在的存储目录；`collection_digest` 是该存储中参考集合的摘要。
- **输出：** 包含 `status`、`backend`、`version`、`hgvs`、`vrs_allele_id`（`ga4gh:VA.<digest>`）、`transcript_store`、`refget_store` 和 `collection_digest` 的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from gtars import refget, reftx

from dnakit.annotation import hgvs_to_vrs_allele_id
from dnakit.exceptions import ConfigurationError

work = Path(tempfile.mkdtemp())
store = refget.RefgetStore.on_disk(str(work / "refget"))
store.set_quiet(True)
meta, _ = store.add_sequence_collection_from_fasta("tiny.fa")
sequence = store.list_sequences()[0]
accession = "SQ." + sequence.sha512t24u

builder = reftx.TxStoreBuilder()
builder.add_chrom_mapping("chr1", accession)
builder.add_transcript(
    {
        "accession": "NM_TEST.1",
        "gene": "TEST",
        "chrom": accession,
        "strand": "+",
        "cds_start": 0,
        "cds_end": 39,
        "exons": [{"start": 0, "end": 39}],
        "mane": None,
    }
)
builder.build(str(work / "tx.reftx"))

args = (str(work / "tx.reftx"), str(work / "refget"), meta.digest)
result = hgvs_to_vrs_allele_id("NM_TEST.1:c.4G>A", *args, verbose=False)
print(result["status"], result["vrs_allele_id"])
try:
    hgvs_to_vrs_allele_id("NM_TEST.1:c.4C>A", *args, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "REF mismatch" in str(error))
```

- **示例输出：**

```text
completed ga4gh:VA.931BbxClFtmksLkac20ozAYn01o1uv_f
ConfigurationError True
```

示例用 gtars 从 `tiny.fa` 构建只含一条序列的 refget 存储和只含一个转录本的存储，然后映射 `c.4G>A`（编码序列的第 4 个碱基，是 G）和一个故意写错参考碱基的变异。

## 27) `dnakit.visualization.compare_genomes_mmseqs()`

- **作用：** 用 MMseqs2 的互为最佳的蛋白命中把相邻基因组连起来，并把连线画在 `result.png` 中。
- **处理方法：** 用 pyGenomeViz 读取 GenBank 文件，把每个基因组翻译后的 CDS 蛋白写入 `<name>.faa`，对每一对相邻的基因组运行 `mmseqs easy-rbh <a.faa> <b.faa> <result.tsv> <tmp> --threads N -e EVALUE -v 0`。命中被转成带百分比一致性的比对坐标；全部写入 `align_coords.tsv`，一致性不低于 `min_identity` 的才作为 CDS 轨道之间的连线绘出，并配一致性颜色条。MMseqs2 的临时文件夹会被删除。
- **API：** `compare_genomes_mmseqs(gbk_resources: 'Sequence[str]', output_dir: 'str', *, evalue: 'float' = 0.001, min_identity: 'float' = 0.0, threads: 'int' = 1, fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`）和外部程序 MMseqs2（`mmseqs`）在 PATH 中或使用 `mmseqs_path`。
- **输入：** `gbk_resources` 是至少两个 GenBank 文件（`file:start-end` 语法可选取范围；基因组名取自记录，只能含字母、数字、`.`、`_` 和 `-`，且必须唯一）；`output_dir` 是新的或空的文件夹；`evalue` 和 `fig_width` 是正数；`min_identity` 是 0 到 100 的百分比；`threads` 和 `dpi` 是正整数；`timeout_seconds` 是每次程序调用的时限。基因组中没有任何带翻译的 CDS 时抛出 `ConfigurationError`。
- **输出：** 包含 `status`（`completed`，程序调用失败时为 `failed`）、`backend`、`commands`、`return_codes`、`outputs`（命中表、`align_coords.tsv` 和 `result.png`）、`image_path`、`log_tail` 和 `data`（`genomes`、`alignment_count` 和 `alignments`）的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_mmseqs

out = Path(tempfile.mkdtemp()) / "mmseqs_out"
try:
    compare_genomes_mmseqs(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        mmseqs_path="/nonexistent/mmseqs",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "mmseqs was not found" in str(error), out.exists())
```

- **示例输出：**

```text
BackendUnavailableError True False
```

本环境没有安装 MMseqs2，所以示例展示程序路径不存在时的错误。GenBank 文件会先被读取，因此必须有效。

安装了依赖后的用法：

```python
from dnakit.visualization import compare_genomes_mmseqs

result = compare_genomes_mmseqs(["genomeA.gbk", "genomeB.gbk"], "mmseqs_out", min_identity=30, dpi=150)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 28) `dnakit.visualization.compare_genomes_mummer()`

- **作用：** 用 MUMmer 比对相邻基因组，并把比对连线画在 `result.png` 中。
- **处理方法：** 把每个基因组写成 FASTA 文件，对每一对相邻基因组运行 `promer`（`seqtype="protein"`）或 `nucmer`（`"nucleotide"`）并带 `--mum`，再运行 `delta-filter -m`（`one_to_one` 时用 `-1`），再运行 `show-coords -H -T -r -k`。坐标由 pyGenomeViz 解析，平移到基因组坐标，写入 `align_coords.tsv`，长度不小于 `min_length` 且一致性不低于 `min_identity`（百分比）的比对被画成连线。临时的 delta 和坐标文件会被删除。
- **API：** `compare_genomes_mummer(gbk_resources: 'Sequence[str]', output_dir: 'str', *, seqtype: 'str' = 'protein', min_length: 'int' = 0, min_identity: 'float' = 0.0, one_to_one: 'bool' = False, fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, mummer_dir: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`）以及 MUMmer 4 程序 `promer` 或 `nucmer`、`delta-filter` 和 `show-coords` 在 PATH 中，或位于文件夹 `mummer_dir` 中。
- **输入：** `gbk_resources` 是至少两个 GenBank 文件；`output_dir` 是新的或空的文件夹；`seqtype` 是 `protein` 或 `nucleotide`；`min_length` 是非负整数；`min_identity` 是 0 到 100 的百分比；`one_to_one` 是布尔值；`fig_width`、`dpi`、`timeout_seconds`；`mummer_dir` 是存放这些程序的可选文件夹。
- **输出：** 包含 `status`（`completed` 或 `failed`）、`backend`、`commands`、`return_codes`、`outputs`（`align_coords.tsv` 和 `result.png`）、`image_path`、`log_tail` 和 `data`（`genomes`、`alignment_count` 和 `alignments`）的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_mummer

out = Path(tempfile.mkdtemp()) / "mummer_out"
try:
    compare_genomes_mummer(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        seqtype="nucleotide",
        mummer_dir="/nonexistent/mummer",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "nucmer was not found" in str(error), out.exists())
```

- **示例输出：**

```text
BackendUnavailableError True False
```

本环境没有安装 MUMmer，所以示例把 `mummer_dir` 指向一个不存在的文件夹并展示错误。缺少程序时，输出文件夹不会被创建。

安装了依赖后的用法：

```python
from dnakit.visualization import compare_genomes_mummer

result = compare_genomes_mummer(["genomeA.gbk", "genomeB.gbk"], "mummer_out", seqtype="nucleotide", min_identity=70)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 29) `dnakit.visualization.compare_genomes_progressive_mauve()`

- **作用：** 用 progressiveMauve 比对序列文件，并把保守块及其连线画在 `result.png` 中。
- **处理方法：** 把输入文件复制到 `seqfiles` 文件夹，运行 `progressiveMauve --output=mauve.xmfa --backbone-output=mauve_bbone.tsv <files>`，把骨架文件解析为块坐标（短于 20 bp 的块和缺失于至少一个基因组的块会被丢弃，块按参考基因组 `refid` 排序），写入 `align_coords.tsv`，并为每个基因组画一条轨道，每个块一个彩色方框，相邻基因组之间有连线。
- **API：** `compare_genomes_progressive_mauve(seq_files: 'Sequence[str]', output_dir: 'str', *, refid: 'int' = 0, cmap: 'str' = 'hsv', fig_width: 'float' = 15.0, dpi: 'int' = 300, timeout_seconds: 'float' = 3600.0, progressive_mauve_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`）和外部程序 progressiveMauve 在 PATH 中或使用 `progressive_mauve_path`。
- **输入：** `seq_files` 是至少两个 FASTA 或 GenBank 文件（基因组以去掉后缀的文件名命名；名称必须唯一，且只能含字母、数字、`.`、`_` 和 `-`）；`output_dir` 是新的或空的文件夹；`refid` 是参考文件的下标；`cmap` 是 matplotlib 颜色映射名；`fig_width`、`dpi`、`timeout_seconds`。
- **输出：** 包含 `status`（`completed` 或 `failed`）、`backend`、`commands`、`return_codes`、`outputs`（`mauve.xmfa`、`mauve_bbone.tsv`、`align_coords.tsv` 和 `result.png`）、`image_path`、`log_tail` 和 `data`（带各自序列文件的 `genomes`、`alignment_count` 和 `alignments`）的字典。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.visualization import compare_genomes_progressive_mauve

out = Path(tempfile.mkdtemp()) / "mauve_out"
try:
    compare_genomes_progressive_mauve(
        ["genomeA.gbk", "genomeB.gbk"],
        str(out),
        progressive_mauve_path="/nonexistent/progressiveMauve",
        verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "progressiveMauve was not found" in str(error), out.exists())
```

- **示例输出：**

```text
BackendUnavailableError True False
```

本环境没有安装 progressiveMauve，所以示例展示程序路径不存在时的错误。

安装了依赖后的用法：

```python
from dnakit.visualization import compare_genomes_progressive_mauve

result = compare_genomes_progressive_mauve(["genomeA.gbk", "genomeB.gbk"], "mauve_out", refid=0)
print(result["status"], result["data"]["alignment_count"], result["image_path"])
```

## 30) `dnakit.visualization.load_example_dataset()`

- **作用：** 从 GitHub 下载 pyGenomeViz 示例 GenBank 数据集，并列出其中的文件和连线。
- **处理方法：** 先对照 `pygenomeviz.config.DATASETS` 检查名称，再调用 `pygenomeviz.utils.load_example_dataset(name, cache_dir=..., overwrite_cache=...)`，它从 pygenomeviz-data 仓库把 GenBank 文件下载到缓存文件夹（或复用缓存），并返回文件路径和连线定义。
- **API：** `load_example_dataset(name: 'str', *, cache_dir: 'str | None' = None, overwrite_cache: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`）和访问 GitHub 的网络（不需要密钥）；已缓存的数据集可离线复用。
- **输入：** `name` 是 pyGenomeViz 中的数据集名之一（未知名称在下载之前抛出 `ConfigurationError`）；`cache_dir` 是可选文件夹；`overwrite_cache` 是布尔值。
- **输出：** 包含 `status`、`service`、`source_url`、`dataset`、`request`、`retrieved_at_utc`、`gbk_files`（路径）、`link_count` 和 `links` 的字典。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import load_example_dataset

try:
    load_example_dataset("no_such_dataset", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Unknown pyGenomeViz example dataset" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

示例不访问网络，只展示未知的数据集名会先被拒绝。

安装了依赖后的用法：

```python
from dnakit.visualization import load_example_dataset

result = load_example_dataset("escherichia_coli", cache_dir="pgv_cache")
print(len(result["gbk_files"]), result["link_count"])
```

## 31) `dnakit.visualization.load_example_gff()`

- **作用：** 从 GitHub 把 pyGenomeViz 示例 GFF 文件下载到缓存文件夹，并返回其路径。
- **处理方法：** 先对照 `pygenomeviz.config.EXAMPLE_GFF_FILES` 检查文件名，再调用 `pygenomeviz.utils.load_example_gff(filename, cache_dir=..., overwrite_cache=...)`。
- **API：** `load_example_gff(filename: 'str', *, cache_dir: 'str | None' = None, overwrite_cache: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`）和访问 GitHub 的网络（不需要密钥）；已缓存的文件可离线复用。
- **输入：** `filename` 是 pyGenomeViz 中的示例 GFF 名之一（未知名称在下载之前抛出 `ConfigurationError`）；`cache_dir` 是可选文件夹；`overwrite_cache` 是布尔值。
- **输出：** 包含 `status`、`service`、`source_url`、`filename`、`request`、`retrieved_at_utc` 和 `gff_path`（缓存文件的绝对路径）的字典。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import load_example_gff

try:
    load_example_gff("no_such_file.gff", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "Unknown pyGenomeViz example GFF" in str(error))
```

- **示例输出：**

```text
ConfigurationError True
```

示例不访问网络，只展示未知的文件名会先被拒绝。

安装了依赖后的用法：

```python
from dnakit.visualization import load_example_gff

result = load_example_gff("enterobacteria_phage.gff", cache_dir="pgv_cache")
print(result["gff_path"])
```

## 32) `dnakit.visualization.plot_genome_features()`

- **作用：** 把 GenBank 基因组的 CDS 特征画成 pyGenomeViz 的堆叠轨道，保存在一个 PNG 中。
- **处理方法：** 用 pyGenomeViz 读取 GenBank 文件，每个基因组添加一条特征轨道，把每个带 `translation` 限定符且位于基因组范围内的 CDS 画成沿其链方向的箭头（`bigarrow` 或 `arrow`）。链方向读自 `feature.location.strand`，因为上游代码读取的 `feature.strand` 在 Biopython 1.88 中已不存在。图以 `bbox_inches="tight"` 保存。
- **API：** `plot_genome_features(gbk_resources: 'Sequence[str]', image_path: 'str', *, fig_width: 'float' = 15.0, fig_track_height: 'float' = 1.0, feature_color: 'str' = 'orange', feature_plotstyle: 'str' = 'bigarrow', align_type: 'str' = 'left', tick_style: 'str | None' = None, dpi: 'int' = 300, verbose: bool = True) -> 'dict[str, Any]'`.
- **运行要求：** 需要可选包 pygenomeviz（额外依赖 `dnakit[bio-pygenomeviz]`），它会带来 matplotlib 和 Biopython。不需要外部程序。
- **输入：** `gbk_resources` 是至少一个 GenBank 文件（`file:start-end` 可选取范围）；`image_path` 是已存在文件夹中以 `.png` 结尾的新路径；`fig_width` 和 `fig_track_height` 是正数；`feature_color` 是 matplotlib 颜色；`feature_plotstyle` 是 `bigarrow` 或 `arrow`；`align_type` 是 `left`、`center` 或 `right`；`tick_style` 是 None、`bar` 或 `axis`；`dpi` 是正整数。
- **输出：** 包含 `status`、`image_path` 和 `data` 的字典；`data["genomes"]` 中每个基因组有 `name`、`min_range`、`max_range`、`range_size`、`reverse`、`full_genome_length` 和 `cds_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.visualization import plot_genome_features

image = Path(tempfile.mkdtemp()) / "features.png"
result = plot_genome_features(
    ["genomeA.gbk", "genomeB.gbk"], str(image), dpi=80, verbose=False
)
for genome in result["data"]["genomes"]:
    print(genome)
print(image.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
```

- **示例输出：**

```text
{'name': 'genomeA', 'min_range': 0, 'max_range': 84, 'range_size': 84, 'reverse': False, 'full_genome_length': 84, 'cds_count': 2}
{'name': 'genomeB', 'min_range': 0, 'max_range': 86, 'range_size': 86, 'reverse': False, 'full_genome_length': 86, 'cds_count': 2}
True
```
