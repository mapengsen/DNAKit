# deepTools：覆盖度、比较、汇总与 reads 过滤

`dnakit.epigenome` 中的这组函数是对 [deepTools](https://github.com/deeptools/deepTools)（4 版）的带类型适配器。其中 21 个各自运行一个 deepTools 命令行程序（共 20 个不同的程序，`computeMatrix` 被两个函数使用：`bamCoverage`、`bamCompare`、`bigwigCompare`、`bigwigAverage`、`computeMatrix`、`multiBamSummary`、`multiBigwigSummary`、`plotHeatmap`、`plotProfile`、`plotCorrelation`、`plotPCA`、`plotFingerprint`、`plotCoverage`、`plotEnrichment`、`estimateReadFiltering`、`estimateScaleFactor`、`bamPEFragmentSize`、`alignmentSieve`、`computeGCBias`、`correctGCBias`），使用固定的、带类型的选项列表，不经 shell；`bam_mapping_statistics()` 则直接导入 deepTools，在进程内调用 `openBam()`。它们用于生成覆盖度与比较轨道、汇总并绘制样本图、计算并绘制区域矩阵、估计缩放因子和 reads 过滤比例，以及计算并校正 GC 偏好。它们没有暴露这些工具的全部选项，也不做 peak 调用或结果解读。

deepTools 以 MIT 许可证发布（已在已安装的 4.0.0 包的 `LICENSE.txt` 中核对）。DNAKit 把这些程序当作独立进程运行，不复制 deepTools 的任何代码。

使用 `from dnakit.epigenome import ...` 导入。deepTools 是可选依赖（`pip install "dnakit[bio-deeptools]"`），仅在调用函数时才使用。命令行适配器从 PATH 或 `deeptools_path`（某个可执行文件的完整路径）查找程序；找不到时抛出带安装提示的 `BackendUnavailableError`。每个适配器都会检查输入文件存在、输出路径尚不存在（否则抛出 `ConfigurationError`）。工具以错误退出时不会抛出异常：返回的 `status` 为 `failed`，错误信息在 `log_tail` 中。除非另有说明，每个适配器都返回「输出」条目所描述的字典。

## 示例数据

两个示例文件都很小。`reads.sam` 是按坐标排序的 SAM 文件，含 `chr1`（3000 碱基）和 `chr2`（2000 碱基）上的 24 对 reads（48 条比对，不含序列）；示例用 pysam 把它转换为带索引的 BAM 文件。`genes.bed` 含两个区域。需要 bigWig 文件的示例用 pyBigWig 自己写出，两个 GC 偏好示例则自己构建随机 2bit 基因组和 reads。每个示例先把当前 Python 所在目录加入 PATH，以便找到同一环境中安装的 deepTools 命令；如果这些命令已在 PATH 中，则不需要这一行。

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:3000
@SQ	SN:chr2	LN:2000
p1	99	chr1	100	60	50M	=	250	200	*	*
p1	147	chr1	250	60	50M	=	100	-200	*	*
p2	99	chr1	262	60	50M	=	412	200	*	*
p2	147	chr1	412	60	50M	=	262	-200	*	*
p3	99	chr1	425	60	50M	=	575	200	*	*
p3	147	chr1	575	60	50M	=	425	-200	*	*
p4	99	chr1	587	60	50M	=	737	200	*	*
p4	147	chr1	737	60	50M	=	587	-200	*	*
p5	99	chr1	750	60	50M	=	900	200	*	*
p5	147	chr1	900	60	50M	=	750	-200	*	*
p6	99	chr1	912	60	50M	=	1062	200	*	*
p6	147	chr1	1062	60	50M	=	912	-200	*	*
p7	99	chr1	1075	60	50M	=	1225	200	*	*
p7	147	chr1	1225	60	50M	=	1075	-200	*	*
p8	99	chr1	1237	60	50M	=	1387	200	*	*
p8	147	chr1	1387	60	50M	=	1237	-200	*	*
p9	99	chr1	1400	60	50M	=	1550	200	*	*
p9	147	chr1	1550	60	50M	=	1400	-200	*	*
p10	99	chr1	1562	60	50M	=	1712	200	*	*
p10	147	chr1	1712	60	50M	=	1562	-200	*	*
p11	99	chr1	1725	60	50M	=	1875	200	*	*
p11	147	chr1	1875	60	50M	=	1725	-200	*	*
p12	99	chr1	1887	60	50M	=	2037	200	*	*
p12	147	chr1	2037	60	50M	=	1887	-200	*	*
p13	99	chr1	2050	60	50M	=	2200	200	*	*
p13	147	chr1	2200	60	50M	=	2050	-200	*	*
p14	99	chr1	2212	60	50M	=	2362	200	*	*
p14	147	chr1	2362	60	50M	=	2212	-200	*	*
p15	99	chr1	2375	60	50M	=	2525	200	*	*
p15	147	chr1	2525	60	50M	=	2375	-200	*	*
p16	99	chr1	2537	60	50M	=	2687	200	*	*
p16	147	chr1	2687	60	50M	=	2537	-200	*	*
p17	99	chr2	100	60	50M	=	250	200	*	*
p17	147	chr2	250	60	50M	=	100	-200	*	*
p18	99	chr2	300	60	50M	=	450	200	*	*
p18	147	chr2	450	60	50M	=	300	-200	*	*
p19	99	chr2	500	60	50M	=	650	200	*	*
p19	147	chr2	650	60	50M	=	500	-200	*	*
p20	99	chr2	700	60	50M	=	850	200	*	*
p20	147	chr2	850	60	50M	=	700	-200	*	*
p21	99	chr2	900	60	50M	=	1050	200	*	*
p21	147	chr2	1050	60	50M	=	900	-200	*	*
p22	99	chr2	1100	60	50M	=	1250	200	*	*
p22	147	chr2	1250	60	50M	=	1100	-200	*	*
p23	99	chr2	1300	60	50M	=	1450	200	*	*
p23	147	chr2	1450	60	50M	=	1300	-200	*	*
p24	99	chr2	1500	60	50M	=	1650	200	*	*
p24	147	chr2	1650	60	50M	=	1500	-200	*	*
```

**示例文件：** `genes.bed`

```text
chr1	1000	2000	g1	0	+
chr1	2200	2800	g2	0	-
```

## 1) `dnakit.epigenome.alignment_sieve()`

- **作用：** 过滤 BAM 文件中的比对，并可选地移动其端点，输出新的 BAM 或 BEDPE 文件。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `alignmentSieve`（超时 300 秒）。参数为 `--bam`、`-o`、`-p`；给出 `shift` 时加 `--shift`，`bed_output` 为真时加 `--BED`。没有暴露其他过滤选项，因此使用默认值时该工具只是重写比对。
- **API：** `dnakit.epigenome.alignment_sieve(bam: 'str', output: 'str', *, shift: 'Sequence[int] | None' = None, bed_output: 'bool' = False, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `alignmentSieve`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bam` 必须是已存在的文件，`output` 必须是新路径（不能已存在，其目录必须存在）。`shift` 为 `None` 或恰好四个整数，否则抛出 `ConfigurationError`。与 deepTools 一致，前两个数是 read 1 在左侧的片段的左端和右端位移，后两个数是 read 1 在右侧的片段的左端和右端位移（正值表示端点向右移）。适配器的报错信息把这四个值描述为正负链的起点和终点，与 deepTools 的帮助不符。`bed_output=True` 时输出 BEDPE 而不是 BAM。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import alignment_sieve

result = alignment_sieve("a.bam", "shifted.bam", shift=[1, -1, 1, -1], verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
with pysam.AlignmentFile("a.bam") as before, pysam.AlignmentFile("shifted.bam") as after:
    first_before = [read.reference_start for read in before][:3]
    first_after = [read.reference_start for read in after][:3]
print(first_before, first_after)
```

- **示例输出：**

```text
completed alignmentSieve shifted.bam
[99, 249, 261] [100, 249, 262]
```

## 2) `dnakit.epigenome.bam_compare()`

- **作用：** 逐个 bin 比较两个 BAM 文件，结果写成 bigWig 或 bedGraph 轨道。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `bamCompare`（超时 300 秒）。参数为 `-b1`、`-b2`、`-o`、`-of`、`--binSize`、`--operation`、`--pseudocount`、`-p`。由于没有暴露缩放选项，使用 bamCompare 的默认缩放方法（`readCount`）。
- **API：** `dnakit.epigenome.bam_compare(bam1: 'str', bam2: 'str', output: 'str', *, operation: '_OPERATIONS' = 'log2', output_format: "Literal['bigwig', 'bedgraph']" = 'bigwig', bin_size: 'int' = 50, pseudocount: 'float' = 1.0, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bamCompare`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bam1` 和 `bam2` 必须是已存在且已建索引的 BAM 文件，`output` 必须是新路径。`operation` 取 `log2`、`ratio`、`subtract`、`add`、`mean`、`reciprocal_ratio`、`first`、`second` 之一；`output_format` 为 `bigwig` 或 `bedgraph`；`bin_size` 为 bin 宽度（碱基）；`pseudocount` 在计算比值前加到两个样本上。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_compare

with pysam.AlignmentFile("a.bam") as source, pysam.AlignmentFile("b.bam", "wb", template=source) as half:
    for number, read in enumerate(source):
        if number % 2 == 0:
            half.write(read)
pysam.index("b.bam")
result = bam_compare("a.bam", "b.bam", "diff.bedgraph", operation="subtract", output_format="bedgraph", bin_size=500, verbose=False)
print(result["status"], result["return_code"])
print(open("diff.bedgraph").read().splitlines()[:3])
```

- **示例输出：**

```text
completed 0
['chr1\t0\t500\t-0.5', 'chr1\t500\t2500\t0', 'chr1\t2500\t3000\t0.5']
```

## 3) `dnakit.epigenome.bam_coverage()`

- **作用：** 从一个 BAM 文件生成覆盖度轨道（bigWig 或 bedGraph）。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `bamCoverage`（超时 300 秒）。参数为 `-b`、`-o`、`-of`、`--binSize`、`--normalizeUsing`、`-p`；给出时再加 `--effectiveGenomeSize`、`--extendReads`、`--scaleFactor`、`--region`。
- **API：** `dnakit.epigenome.bam_coverage(bam: 'str', output: 'str', *, output_format: "Literal['bigwig', 'bedgraph']" = 'bigwig', bin_size: 'int' = 50, normalize_using: '_NORMALIZATION' = 'None', effective_genome_size: 'int | None' = None, extend_reads: 'int | None' = None, scale_factor: 'float | None' = None, region: 'str | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bamCoverage`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bam` 必须是已存在且已建索引的 BAM 文件，`output` 必须是新路径。`bin_size` 至少为 1。`normalize_using` 取 `RPGC`、`CPM`、`BPM`、`RPKM` 或 `None`（默认字符串 `"None"`）；`RPGC` 需要 `effective_genome_size`（缺少时 deepTools 在日志中报错，`status` 变为 `failed`）。`extend_reads`、`scale_factor` 和 `region`（如 `chr1:1000:2000`）可选。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_coverage

result = bam_coverage("a.bam", "cov.bedgraph", output_format="bedgraph", bin_size=500, verbose=False)
print(result["status"], result["tool"], result["command"][0])
print(open("cov.bedgraph").read().splitlines()[:4])
```

- **示例输出：**

```text
completed bamCoverage bamCoverage
['chr1\t0\t500\t5', 'chr1\t500\t2500\t6', 'chr1\t2500\t3000\t3', 'chr2\t0\t1000\t5']
```

## 4) `dnakit.epigenome.bam_pe_fragment_size()`

- **作用：** 估计双端 BAM 文件的片段长度分布，并画成直方图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `bamPEFragmentSize`（超时 300 秒）。参数为 `-b`、`-hist`、`-p`；给出 `table_output` 时加 `--table`，给出 `max_fragment_length` 时加 `--maxFragmentLength`。该工具不会拒绝单端 BAM 文件：它照常完成，只在日志中提示数据可能不是双端。
- **API：** `dnakit.epigenome.bam_pe_fragment_size(bams: 'Sequence[str]', histogram_output: 'str', *, table_output: 'str | None' = None, max_fragment_length: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bamPEFragmentSize`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 是非空的、已存在的双端 BAM 文件列表；`histogram_output` 是新的图片路径；`table_output` 是可选的新路径，用于保存制表符分隔的汇总表；`max_fragment_length` 可选。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_pe_fragment_size

result = bam_pe_fragment_size(["a.bam"], "hist.png", table_output="sizes.tab", verbose=False)
print(result["status"], [os.path.basename(path) for path in result["outputs"]])
print(open("sizes.tab").read().splitlines()[0].split("\t")[1:5])
print(open("sizes.tab").read().splitlines()[1].split("\t")[1:5])
```

- **示例输出：**

```text
completed ['hist.png', 'sizes.tab']
['Frag. Sampled', 'Frag. Len. Min.', 'Frag. Len. 1st. Qu.', 'Frag. Len. Mean']
['24', '200', '200.0', '200.0']
```

## 5) `dnakit.epigenome.bigwig_compare()`

- **作用：** 逐个 bin 比较两个 bigWig 轨道，输出 bigWig。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `bigwigCompare`（超时 300 秒）。参数为 `-b1`、`-b2`、`-o`、`--binSize`、`--operation`、`--pseudocount`、`-p`。输出始终是 bigWig。
- **API：** `dnakit.epigenome.bigwig_compare(bigwig1: 'str', bigwig2: 'str', output: 'str', *, operation: '_OPERATIONS' = 'log2', bin_size: 'int' = 50, pseudocount: 'float' = 1.0, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bigwigCompare`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bigwig1` 和 `bigwig2` 必须是已存在的 bigWig 文件，`output` 必须是新路径。`operation` 的取值与 `bam_compare` 相同；`bin_size` 为 bin 宽度（碱基）；`pseudocount` 在计算比值前加上。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import bigwig_compare

write_bigwig("one.bw", [4.0, 6.0, 8.0, 2.0, 2.0])
write_bigwig("two.bw", [2.0, 2.0, 2.0, 1.0, 1.0])
result = bigwig_compare("one.bw", "two.bw", "diff.bw", operation="subtract", bin_size=1000, verbose=False)
print(result["status"], result["return_code"])
print(pyBigWig.open("diff.bw").intervals("chr1"))
```

- **示例输出：**

```text
completed 0
((0, 1000, 2.0), (1000, 2000, 4.0), (2000, 3000, 6.0))
```

## 6) `dnakit.epigenome.estimate_scale_factor()`

- **作用：** 根据采样的 reads 计数，估计两个或多个 BAM 文件的缩放因子。
- **处理方法：** 检查输入，构造参数列表（`-b`、`-l`、`-n`、`-nl`、`-p`），在第一个 BAM 文件所在目录下不经 shell 运行 `estimateScaleFactor`。该工具不写文件，它打印的缩放因子包含在返回结果中。
- **API：** `dnakit.epigenome.estimate_scale_factor(bams: 'Sequence[str]', *, sample_window_length: 'int' = 1000, number_of_samples: 'int' = 1000000, normalization_length: 'int' = 1000, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `estimateScaleFactor`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 至少包含两个已存在的 BAM 文件（否则抛出 `ConfigurationError`）。`sample_window_length`（`-l`）、`number_of_samples`（`-n`）和 `normalization_length`（`-nl`）控制采样；`number_of_processors` 对应 `-p`。
- **输出：** `status`（退出码为 0 时为 `completed`，否则为 `failed`）、`backend`、`tool`、`command`、`return_code` 和 `stdout_tail`（工具输出的最后 4000 个字符，包含缩放因子和统计量，例如 `size_factors_based_on_median`）。没有 `outputs` 字段。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import estimate_scale_factor

shutil.copy("a.bam", "b.bam")
shutil.copy("a.bam.bai", "b.bam.bai")
result = estimate_scale_factor(["a.bam", "b.bam"], sample_window_length=500, number_of_samples=100, verbose=False)
print(result["status"], result["return_code"])
for line in result["stdout_tail"].splitlines():
    if line.startswith(("size_factors_based_on_median", "sites_sampled")):
        print(line)
```

- **示例输出：**

```text
completed 0
size_factors_based_on_median: [1. 1.]
sites_sampled: 82
```

## 7) `dnakit.epigenome.plot_coverage()`

- **作用：** 绘制一个或多个 BAM 文件的覆盖度分布图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotCoverage`（超时 300 秒）。参数为 `-b`、`-o`、`-n`；给出 `labels` 时加 `--labels`。该工具还会打印覆盖度统计表，内容包含在 `log_tail` 中。
- **API：** `dnakit.epigenome.plot_coverage(bams: 'Sequence[str]', plot_output: 'str', *, labels: 'Sequence[str] | None' = None, number_of_samples: 'int' = 1000000, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotCoverage`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 是非空的、已存在且已建索引的 BAM 文件列表；`plot_output` 是新的图片路径；`labels` 可选，数量必须与 BAM 文件一致；`number_of_samples` 是采样的基因组位置数。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_coverage

result = plot_coverage(["a.bam"], "coverage.png", labels=["sample"], number_of_samples=100, verbose=False)
print(result["status"], os.path.isfile("coverage.png"))
print(result["log_tail"].splitlines()[0].split("	"))
```

- **示例输出：**

```text
completed True
['sample', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
```

## 8) `dnakit.epigenome.plot_fingerprint()`

- **作用：** 绘制类 ChIP BAM 文件的累积 reads 覆盖度指纹图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotFingerprint`（超时 300 秒）。参数为 `-b`、`-plot`、`-p`；给出时再加 `--extendReads` 和 `--minMappingQuality`。
- **API：** `dnakit.epigenome.plot_fingerprint(bams: 'Sequence[str]', plot_output: 'str', *, extend_reads: 'int | None' = None, min_mapping_quality: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotFingerprint`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 是非空的、已存在且已建索引的 BAM 文件列表；`plot_output` 是新的图片路径；`extend_reads` 和 `min_mapping_quality` 为可选整数。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_fingerprint

result = plot_fingerprint(["a.bam"], "fingerprint.png", min_mapping_quality=10, verbose=False)
print(result["status"], os.path.isfile("fingerprint.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **示例输出：**

```text
completed True
['-b', '-plot', '-p', '--minMappingQuality']
```

## 9) `dnakit.epigenome.bigwig_average()`

- **作用：** 将多个 bigWig 轨道逐 bin 取平均，合并成一个 bigWig。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `bigwigAverage`（超时 300 秒）。参数为 `-b`、`-o`、`--binSize`、`-p`。
- **API：** `dnakit.epigenome.bigwig_average(bigwigs: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 50, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bigwigAverage`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bigwigs` 是非空的、已存在的 bigWig 文件列表；`output` 是新路径；`bin_size` 为 bin 宽度（碱基）。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import bigwig_average

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [3.0, 2.0, 1.0, 0.0, 1.0])
result = bigwig_average(["one.bw", "two.bw"], "mean.bw", bin_size=1000, verbose=False)
print(result["status"], result["return_code"])
print(pyBigWig.open("mean.bw").intervals("chr1"))
```

- **示例输出：**

```text
completed 0
((0, 3000, 2.0),)
```

## 10) `dnakit.epigenome.compute_matrix_reference_point()`

- **作用：** 围绕每个区域的参考点（TSS、TES 或中心）计算信号值矩阵。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `computeMatrix`（超时 300 秒）。子命令为 `reference-point`，参数为 `--referencePoint`、`-S`、`-R`、`-o`、`-b`、`-a`、`--binSize`、`-p`。输出是 deepTools 的 gzip 矩阵文件，可作为 `plot_heatmap` 和 `plot_profile` 的输入。
- **API：** `dnakit.epigenome.compute_matrix_reference_point(score_files: 'Sequence[str]', regions: 'Sequence[str]', output_matrix: 'str', *, reference_point: "Literal['TSS', 'TES', 'center']" = 'TSS', before_region_length: 'int' = 500, after_region_length: 'int' = 500, bin_size: 'int' = 10, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `computeMatrix`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `score_files`（bigWig 文件）和 `regions`（BED 或 GTF 文件）是非空的、已存在的文件列表；`output_matrix` 是新路径；`reference_point` 取 `TSS`、`TES` 或 `center`；`before_region_length` 和 `after_region_length` 是参考点两侧的距离（碱基）；`bin_size` 为 bin 宽度。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import gzip
import json

from dnakit.epigenome import compute_matrix_reference_point

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
result = compute_matrix_reference_point(
    ["one.bw"], ["genes.bed"], "matrix.gz", reference_point="TSS",
    before_region_length=500, after_region_length=500, bin_size=100, verbose=False,
)
print(result["status"], result["return_code"])
with gzip.open("matrix.gz", "rt") as handle:
    meta = json.loads(handle.readline()[1:])
    rows = handle.read().splitlines()
print(meta["ref point"], meta["sample_boundaries"], meta["group_boundaries"], len(rows))
```

- **示例输出：**

```text
completed 0
['TSS'] [0, 10] [0, 2] 2
```

## 11) `dnakit.epigenome.compute_matrix_scale_regions()`

- **作用：** 将区域缩放到统一长度并带上两侧侧翼，计算信号值矩阵。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `computeMatrix`（超时 300 秒）。子命令为 `scale-regions`，参数为 `-S`、`-R`、`-o`、`-m`、`-b`、`-a`、`--binSize`、`-p`。每个区域都被拉伸或压缩到 `region_body_length`。
- **API：** `dnakit.epigenome.compute_matrix_scale_regions(score_files: 'Sequence[str]', regions: 'Sequence[str]', output_matrix: 'str', *, region_body_length: 'int' = 1000, before_region_length: 'int' = 500, after_region_length: 'int' = 500, bin_size: 'int' = 10, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `computeMatrix`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `score_files` 和 `regions` 是非空的、已存在的文件列表；`output_matrix` 是新路径；`region_body_length` 是每个区域缩放后的长度；`before_region_length` 和 `after_region_length` 是侧翼长度；`bin_size` 为 bin 宽度。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import gzip
import json

from dnakit.epigenome import compute_matrix_scale_regions

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
result = compute_matrix_scale_regions(
    ["one.bw"], ["genes.bed"], "scaled.gz", region_body_length=300,
    before_region_length=100, after_region_length=100, bin_size=50, verbose=False,
)
print(result["status"], result["return_code"])
with gzip.open("scaled.gz", "rt") as handle:
    meta = json.loads(handle.readline()[1:])
    rows = handle.read().splitlines()
print(meta["body"], meta["sample_boundaries"], len(rows))
```

- **示例输出：**

```text
completed 0
[300] [0, 10] 2
```

## 12) `dnakit.epigenome.multi_bam_summary()`

- **作用：** 在连续的基因组 bin 中统计多个 BAM 文件的 reads 数，并保存结果表供相关性或 PCA 绘图使用。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `multiBamSummary`（超时 300 秒）。子命令为 `bins`，参数为 `-b`、`-o`、`--binSize`、`-n`（bin 之间的间隔）、`-p`；给出 `labels` 时加 `--labels`。输出是 NumPy 的 `.npz` 文件。
- **API：** `dnakit.epigenome.multi_bam_summary(bams: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 10000, distance_between_bins: 'int' = 0, labels: 'Sequence[str] | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `multiBamSummary`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 是非空的、已存在且已建索引的 BAM 文件列表；`output` 是新路径；`bin_size` 为 bin 宽度；`distance_between_bins` 是被采样 bin 之间的间隔（0 表示连续的 bin）；`labels` 可选，数量必须与 BAM 文件一致。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

import numpy as np

from dnakit.epigenome import multi_bam_summary

shutil.copy("a.bam", "b.bam")
shutil.copy("a.bam.bai", "b.bam.bai")
result = multi_bam_summary(["a.bam", "b.bam"], "summary.npz", bin_size=500, labels=["A", "B"], verbose=False)
print(result["status"], result["return_code"])
with np.load("summary.npz") as data:
    print(sorted(data.files), data["matrix"].shape, data["matrix"].sum(axis=0).tolist())
```

- **示例输出：**

```text
completed 0
['labels', 'matrix'] (10, 2) [50.0, 50.0]
```

## 13) `dnakit.epigenome.multi_bigwig_summary()`

- **作用：** 在基因组 bin 中汇总多个 bigWig 轨道的值，并保存结果表。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `multiBigwigSummary`（超时 300 秒）。子命令为 `bins`，参数为 `-b`、`-o`、`--binSize`、`-p`；给出 `labels` 时加 `--labels`。输出是 NumPy 的 `.npz` 文件。
- **API：** `dnakit.epigenome.multi_bigwig_summary(bigwigs: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 10000, labels: 'Sequence[str] | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `multiBigwigSummary`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bigwigs` 是非空的、已存在的 bigWig 文件列表；`output` 是新路径；`bin_size` 为 bin 宽度；`labels` 可选，数量必须与 bigWig 文件一致。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

import numpy as np

from dnakit.epigenome import multi_bigwig_summary

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
result = multi_bigwig_summary(["one.bw", "two.bw"], "summary.npz", bin_size=1000, labels=["one", "two"], verbose=False)
print(result["status"], result["return_code"])
with np.load("summary.npz") as data:
    print(sorted(data.files), [str(label) for label in data["labels"]], data["matrix"].shape)
```

- **示例输出：**

```text
completed 0
['labels', 'matrix'] ['one', 'two'] (5, 2)
```

## 14) `dnakit.epigenome.plot_correlation()`

- **作用：** 根据 multiBamSummary 或 multiBigwigSummary 的结果表绘制样本间相关性图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotCorrelation`（超时 300 秒）。参数为 `-in`、`-c`（相关系数方法）、`-p`（图类型）、`-o`。
- **API：** `dnakit.epigenome.plot_correlation(summary_matrix: 'str', output_figure: 'str', *, method: "Literal['spearman', 'pearson']" = 'spearman', plot_type: "Literal['heatmap', 'scatterplot']" = 'heatmap', deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotCorrelation`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `summary_matrix` 是 `multi_bam_summary` 或 `multi_bigwig_summary` 生成的已存在 `.npz` 文件；`output_figure` 是新的图片路径；`method` 取 `spearman` 或 `pearson`；`plot_type` 取 `heatmap` 或 `scatterplot`。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import multi_bigwig_summary, plot_correlation

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
write_bigwig("three.bw", [5.0, 4.0, 3.0, 2.0, 1.5])
multi_bigwig_summary(["one.bw", "two.bw", "three.bw"], "summary.npz", bin_size=1000, labels=["one", "two", "three"], verbose=False)
result = plot_correlation("summary.npz", "correlation.png", method="pearson", plot_type="heatmap", verbose=False)
print(result["status"], os.path.isfile("correlation.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **示例输出：**

```text
completed True
['-in', '-c', '-p', '-o']
```

## 15) `dnakit.epigenome.plot_heatmap()`

- **作用：** 根据 computeMatrix 的矩阵文件绘制热图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotHeatmap`（超时 300 秒）。参数为 `-m`、`-out`、`--dpi`；给出 `kmeans` 时加 `--kmeans`。
- **API：** `dnakit.epigenome.plot_heatmap(matrix: 'str', output_figure: 'str', *, kmeans: 'int | None' = None, dpi: 'int' = 200, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotHeatmap`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `matrix` 是由 `compute_matrix_reference_point` 或 `compute_matrix_scale_regions` 生成的已存在文件；`output_figure` 是新的图片路径；`kmeans` 为 `None` 或正整数（否则抛出 `ConfigurationError`）；`dpi` 为图像分辨率。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import compute_matrix_reference_point, plot_heatmap

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
compute_matrix_reference_point(["one.bw"], ["genes.bed"], "matrix.gz", bin_size=100, verbose=False)
result = plot_heatmap("matrix.gz", "heatmap.png", dpi=100, verbose=False)
print(result["status"], os.path.isfile("heatmap.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **示例输出：**

```text
completed True
['-m', '-out', '--dpi']
```

## 16) `dnakit.epigenome.plot_pca()`

- **作用：** 对 multiBamSummary 或 multiBigwigSummary 结果表中的样本做主成分分析并作图。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotPCA`（超时 300 秒）。参数为 `-in` 和 `-o`；没有暴露 plotPCA 的其他选项。
- **API：** `dnakit.epigenome.plot_pca(summary_matrix: 'str', output_figure: 'str', *, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotPCA`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `summary_matrix` 是 `multi_bam_summary` 或 `multi_bigwig_summary` 生成的已存在 `.npz` 文件；`output_figure` 是新的图片路径。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import multi_bigwig_summary, plot_pca

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
write_bigwig("two.bw", [2.0, 1.0, 4.0, 3.0, 6.0])
write_bigwig("three.bw", [5.0, 4.0, 3.0, 2.0, 1.5])
multi_bigwig_summary(["one.bw", "two.bw", "three.bw"], "summary.npz", bin_size=1000, labels=["one", "two", "three"], verbose=False)
result = plot_pca("summary.npz", "pca.png", verbose=False)
print(result["status"], os.path.isfile("pca.png"))
print(result["command"][0], result["return_code"])
```

- **示例输出：**

```text
completed True
plotPCA 0
```

## 17) `dnakit.epigenome.plot_profile()`

- **作用：** 根据 computeMatrix 的矩阵文件绘制平均信号曲线。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotProfile`（超时 300 秒）。参数为 `-m`、`-out`、`--dpi`；给出 `kmeans` 时加 `--kmeans`。
- **API：** `dnakit.epigenome.plot_profile(matrix: 'str', output_figure: 'str', *, kmeans: 'int | None' = None, dpi: 'int' = 200, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotProfile`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `matrix` 是由 `compute_matrix_reference_point` 或 `compute_matrix_scale_regions` 生成的已存在文件；`output_figure` 是新的图片路径；`kmeans` 为 `None` 或正整数（否则抛出 `ConfigurationError`）；`dpi` 为图像分辨率。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("genes.bed", work)
os.chdir(work)
import pyBigWig


def write_bigwig(path, values):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([("chr1", 3000), ("chr2", 2000)])
    bw.addEntries(
        ["chr1"] * 3 + ["chr2"] * 2,
        [0, 1000, 2000, 0, 1000],
        ends=[1000, 2000, 3000, 1000, 2000],
        values=values,
    )
    bw.close()

from dnakit.epigenome import compute_matrix_reference_point, plot_profile

write_bigwig("one.bw", [1.0, 2.0, 3.0, 4.0, 5.0])
compute_matrix_reference_point(["one.bw"], ["genes.bed"], "matrix.gz", bin_size=100, verbose=False)
result = plot_profile("matrix.gz", "profile.png", dpi=100, verbose=False)
print(result["status"], os.path.isfile("profile.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **示例输出：**

```text
completed True
['-m', '-out', '--dpi']
```

## 18) `dnakit.epigenome.estimate_read_filtering()`

- **作用：** 估计每个 BAM 文件中有多少 reads 被 deepTools 的 reads 过滤条件去掉，并按原因分类。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `estimateReadFiltering`（超时 300 秒）。参数为 `-b`、`-o`、`--binSize`、`-p`。没有暴露过滤选项，因此使用工具的默认值。输出是制表符分隔的表，每个 BAM 文件一行。
- **API：** `dnakit.epigenome.estimate_read_filtering(bams: 'Sequence[str]', output: 'str', *, bin_size: 'int' = 1000, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `estimateReadFiltering`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 是非空的、已存在且已建索引的 BAM 文件列表；`output` 是新路径；`bin_size` 是采样 bin 的宽度。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import estimate_read_filtering

result = estimate_read_filtering(["a.bam"], "filtering.tab", bin_size=500, verbose=False)
print(result["status"], result["return_code"])
rows = open("filtering.tab").read().splitlines()
print(rows[0].split("\t")[1:3])
print(rows[1].split("\t")[1:3])
```

- **示例输出：**

```text
completed 0
['Total Reads', 'Mapped Reads']
['48', '48']
```

## 19) `dnakit.epigenome.plot_enrichment()`

- **作用：** 绘制落入各组区域的 reads 所占比例。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `plotEnrichment`（超时 300 秒）。参数为 `-b`、`--BED`、`-o`；给出 `labels` 时加 `--labels`。
- **API：** `dnakit.epigenome.plot_enrichment(bams: 'Sequence[str]', regions: 'Sequence[str]', output_figure: 'str', *, labels: 'Sequence[str] | None' = None, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `plotEnrichment`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bams` 和 `regions`（BED 或 GTF 文件）是非空的、已存在的文件列表；`output_figure` 是新的图片路径；`labels` 可选，数量必须与 BAM 文件一致。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

import pysam

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
shutil.copy("genes.bed", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import plot_enrichment

result = plot_enrichment(["a.bam"], ["genes.bed"], "enrichment.png", labels=["sample"], verbose=False)
print(result["status"], os.path.isfile("enrichment.png"))
print([argument for argument in result["command"] if argument.startswith("-")])
```

- **示例输出：**

```text
completed True
['-b', '--BED', '-o', '--labels']
```

## 20) `dnakit.epigenome.compute_gc_bias()`

- **作用：** 针对 2bit 基因组计算 BAM 文件的 GC 偏好频率表。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `computeGCBias`（超时 300 秒）。参数为 `-b`、`-g`、`--effectiveGenomeSize`、`-freq`、`-p`；给出 `fragment_length` 时加 `-l`。该表是 `correct_gc_bias` 的输入；没有暴露该工具的绘图选项。
- **API：** `dnakit.epigenome.compute_gc_bias(bam: 'str', genome_2bit: 'str', effective_genome_size: 'int', frequencies_output: 'str', *, fragment_length: 'int | None' = None, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `computeGCBias`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bam` 是已存在且已建索引的 BAM 文件，`genome_2bit` 是已存在的 2bit 文件；`effective_genome_size` 必须为正整数；`frequencies_output` 是新路径；`fragment_length` 为 `None` 或正整数（单端数据需要给出）。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import random
import struct

import pysam


def write_twobit(path, name, seq):
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(seq) + 3) // 4)
    for i, base in enumerate(seq):
        packed[i // 4] |= code[base] << (6 - 2 * (i % 4))
    record = struct.pack("<IIII", len(seq), 0, 0, 0) + bytes(packed)
    index = struct.pack("B", len(name)) + name.encode() + struct.pack("<I", 16 + 1 + len(name) + 4)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<IIII", 0x1A412743, 0, 1, 0) + index + record)


rng = random.Random(7)
genome = "".join(rng.choice("ACGT") for _ in range(20000))
write_twobit("genome.2bit", "chr1", genome)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(genome)}]}
starts = sorted(rng.randrange(0, 19800) for _ in range(1500))
with pysam.AlignmentFile("g.bam", "wb", header=header) as out:
    for i, start in enumerate(starts):
        read = pysam.AlignedSegment()
        read.query_name = f"r{i}"
        read.query_sequence = genome[start:start + 50]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = "50M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(read)
pysam.index("g.bam")

from dnakit.epigenome import compute_gc_bias

result = compute_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", fragment_length=200, verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
rows = open("freq.txt").read().splitlines()
print(len(rows), len(rows[0].split()))
```

- **示例输出：**

```text
completed computeGCBias freq.txt
201 3
```

## 21) `dnakit.epigenome.correct_gc_bias()`

- **作用：** 使用 `compute_gc_bias` 生成的频率表写出经 GC 校正的 BAM 文件。
- **处理方法：** 检查输入，构造固定的参数列表，在输出文件所在目录下不经 shell 运行 `correctGCBias`（超时 300 秒）。参数为 `-b`、`-g`、`--effectiveGenomeSize`、`-freq`、`-o`、`-p`。
- **API：** `dnakit.epigenome.correct_gc_bias(bam: 'str', genome_2bit: 'str', effective_genome_size: 'int', frequencies: 'str', output_bam: 'str', *, number_of_processors: 'int' = 1, deeptools_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `correctGCBias`（随 deepTools 安装，extra 为 `dnakit[bio-deeptools]`）在 PATH 中，或在 `deeptools_path` 中给出完整路径。示例还用到 pysam 和 pyBigWig，二者随 deepTools 一起安装。
- **输入：** `bam` 是已存在且已建索引的 BAM 文件，`genome_2bit` 是已存在的 2bit 文件，`frequencies` 是 `compute_gc_bias` 生成的表；`effective_genome_size` 必须为正整数；`output_bam` 是新路径。
- **输出：** 字典，包含 `status`（退出码为 0 且所有输出文件都存在时为 `completed`，否则为 `failed`）、`backend`（`deeptools`）、`tool`、`command`（工具名加参数）、`return_code`、`outputs`（实际存在的输出文件）和 `log_tail`（工具输出的最后 4000 个字符）。
- **示例代码：**

```python
import os
import shutil
import sys
import tempfile

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
work = tempfile.mkdtemp(dir=".")
os.chdir(work)
import random
import struct

import pysam


def write_twobit(path, name, seq):
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytearray((len(seq) + 3) // 4)
    for i, base in enumerate(seq):
        packed[i // 4] |= code[base] << (6 - 2 * (i % 4))
    record = struct.pack("<IIII", len(seq), 0, 0, 0) + bytes(packed)
    index = struct.pack("B", len(name)) + name.encode() + struct.pack("<I", 16 + 1 + len(name) + 4)
    with open(path, "wb") as handle:
        handle.write(struct.pack("<IIII", 0x1A412743, 0, 1, 0) + index + record)


rng = random.Random(7)
genome = "".join(rng.choice("ACGT") for _ in range(20000))
write_twobit("genome.2bit", "chr1", genome)
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(genome)}]}
starts = sorted(rng.randrange(0, 19800) for _ in range(1500))
with pysam.AlignmentFile("g.bam", "wb", header=header) as out:
    for i, start in enumerate(starts):
        read = pysam.AlignedSegment()
        read.query_name = f"r{i}"
        read.query_sequence = genome[start:start + 50]
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = "50M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(read)
pysam.index("g.bam")

from dnakit.epigenome import compute_gc_bias, correct_gc_bias

compute_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", fragment_length=200, verbose=False)
result = correct_gc_bias("g.bam", "genome.2bit", 20000, "freq.txt", "corrected.bam", verbose=False)
print(result["status"], result["tool"], os.path.basename(result["outputs"][0]))
print(os.path.getsize("corrected.bam") > 0)
```

- **示例输出：**

```text
completed correctGCBias corrected.bam
True
```

## 22) `dnakit.epigenome.bam_mapping_statistics()`

- **作用：** 统计已建索引的 BAM 或 CRAM 文件中比对上和未比对上的 reads 数，并给出每条参考序列的计数。
- **处理方法：** 导入 `deeptools.plotEnrichment` 并调用其 `openBam(returnStats=True)`：总数（`mapped`、`unmapped`）来自文件索引，每条参考序列的计数来自 `get_index_statistics()`。函数返回前会关闭比对文件句柄。deepTools 在文件缺失、缺少索引或文件不是比对文件时会调用 `sys.exit`，这些退出被转换为 `ConfigurationError`。
- **API：** `dnakit.epigenome.bam_mapping_statistics(bam_path: 'str | Path', *, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 deepTools（extra 为 `dnakit[bio-deeptools]`）；仅在调用函数时才导入，不运行外部程序。示例还用到 pysam，它随 deepTools 一起安装。
- **输入：** `bam_path` 是已建索引的 BAM 或 CRAM 文件；`threads` 是 1 到 256 的整数，仅对 CRAM 有效（此时由 deepTools 自己数 reads）。
- **输出：** `status`、`operation`、`bam_path`、`format`（`BAM` 或 `CRAM`）、`mapped_reads`、`unmapped_reads`、`references`（按头部顺序，每条参考序列含 `name`、`mapped`、`unmapped`）和 `provenance`（包名与版本，以及方法说明）。
- **示例代码：**

```python
import os
import shutil
import tempfile

import pysam

work = tempfile.mkdtemp(dir=".")
shutil.copy("reads.sam", work)
os.chdir(work)
pysam.view("-b", "-o", "a.bam", "reads.sam", catch_stdout=False)
pysam.index("a.bam")

from dnakit.epigenome import bam_mapping_statistics

result = bam_mapping_statistics("a.bam", verbose=False)
print(result["status"], result["format"], result["mapped_reads"], result["unmapped_reads"])
print(result["references"])
```

- **示例输出：**

```text
completed BAM 48 0
[{'name': 'chr1', 'mapped': 32, 'unmapped': 0}, {'name': 'chr2', 'mapped': 16, 'unmapped': 0}]
```
